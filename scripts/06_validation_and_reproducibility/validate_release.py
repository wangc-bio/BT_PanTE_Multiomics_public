#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Validate source-release hygiene and syntax; this does not rerun biology."""

import argparse
import ast
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

REQUIRED_DOCS = [
    "README.md",
    "docs/analysis_workflow.md",
    "docs/data_management_and_release.md",
    "docs/reproducibility.md",
    "scripts/01_genome_and_sample_qc/README.md",
    "scripts/02_phased_F1_assembly_and_parental_origin/README.md",
    "scripts/03_pan_TE_catalogue/README.md",
    "scripts/04A_TE_insertion_inheritance_and_genomic_divergence/README.md",
    "scripts/04B_TE_transcription_and_epigenomic_profiling/README.md",
    "scripts/05_integrated_statistical_modelling/README.md",
    "scripts/06_validation_and_reproducibility/README.md",
]

TEXT_SUFFIXES = {
    ".md", ".py", ".sh", ".slurm", ".r", ".yml", ".yaml",
    ".tsv", ".csv", ".cff", ".toml", ".txt", ".example", ".json",
}
SOURCE_SUFFIXES = {".py", ".sh", ".slurm", ".r"}
RAW_SUFFIXES = (
    ".fastq", ".fastq.gz", ".fq", ".fq.gz",
    ".bam", ".cram", ".sam", ".sra", ".pod5",
)

PRIVATE_PATH_PATTERNS = [
    re.compile(r"/public/home/[A-Za-z0-9_.-]+/", re.IGNORECASE),
    re.compile(r"/public" + r"/software/", re.IGNORECASE),
    re.compile(r"/home/(?!runner/)[A-Za-z0-9_.-]+/", re.IGNORECASE),
    re.compile(r"/Users/[A-Za-z0-9_.-]+/", re.IGNORECASE),
    re.compile(r"[A-Za-z]:\\Users\\[^\\]+", re.IGNORECASE),
]
SENSITIVE_PATTERNS = [
    re.compile(r"ghp_[A-Za-z0-9]{20,}"),
    re.compile(r"github_pat_[A-Za-z0-9_]{20,}"),
    re.compile(r"AKIA[0-9A-Z]{16}"),
    re.compile(r"xox[baprs]-[A-Za-z0-9-]{10,}"),
    re.compile(r"-----BEGIN (?:RSA|OPENSSH|EC) PRIVATE KEY-----"),
]
SAMPLE_ID_PATTERN = re.compile(r"\bYPQ\d{6}\b", re.IGNORECASE)


def add(rows, name, status, detail):
    rows.append({"Check": name, "Status": status, "Detail": detail})


def iter_files(root):
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if ".git" in path.parts or "__pycache__" in path.parts:
            continue
        if path.suffix == ".pyc":
            continue
        yield path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--output", default="")
    parser.add_argument("--max-file-mb", type=float, default=100.0)
    args = parser.parse_args()

    root = Path(args.repo_root).resolve()
    files = list(iter_files(root))
    rows = []
    failures = 0
    warnings = 0

    missing_docs = [item for item in REQUIRED_DOCS if not (root / item).is_file()]
    if missing_docs:
        add(rows, "required_documentation", "FAIL", "; ".join(missing_docs))
        failures += 1
    else:
        add(rows, "required_documentation", "PASS",
            "All six module READMEs and core documentation are present.")

    private_hits = []
    secret_hits = []
    sample_id_hits = []
    raw_hits = []
    oversized_hits = []
    local_env_hits = []

    python_files = []
    shell_files = []
    r_files = []

    for path in files:
        rel = path.relative_to(root).as_posix()
        lower_name = path.name.lower()
        lower_path = rel.lower()

        if lower_name.endswith(RAW_SUFFIXES):
            raw_hits.append(rel)
        if path.stat().st_size > args.max_file_mb * 1024 * 1024:
            oversized_hits.append("{0} ({1:.1f} MB)".format(
                rel, path.stat().st_size / (1024 * 1024)
            ))
        if lower_name == ".env":
            local_env_hits.append(rel)

        if path.suffix.lower() in TEXT_SUFFIXES:
            text = path.read_text(encoding="utf-8", errors="replace")
            for pattern in PRIVATE_PATH_PATTERNS:
                if pattern.search(text):
                    private_hits.append(rel)
                    break
            for pattern in SENSITIVE_PATTERNS:
                if pattern.search(text):
                    secret_hits.append(rel)
                    break
            if SAMPLE_ID_PATTERN.search(text):
                sample_id_hits.append(rel)

        if path.suffix.lower() == ".py":
            python_files.append(path)
        elif path.suffix.lower() in {".sh", ".slurm"}:
            shell_files.append(path)
        elif path.suffix.lower() == ".r":
            r_files.append(path)

    checks = [
        ("private_absolute_paths", private_hits,
         "No personal or cluster-specific absolute paths detected."),
        ("credential_patterns", secret_hits,
         "No common API-token or private-key patterns detected."),
        ("individual_sample_identifiers", sample_id_hits,
         "No individual sample IDs detected in public text files."),
        ("raw_sequencing_or_alignment_files", raw_hits,
         "No raw reads or alignment files are included."),
        ("oversized_files", oversized_hits,
         "No file exceeds the configured {:.1f} MB limit.".format(args.max_file_mb)),
        ("local_environment_files", local_env_hits,
         "No local .env files are included; .env.example files are allowed."),
    ]
    for name, hits, detail in checks:
        if hits:
            add(rows, name, "FAIL", "; ".join(sorted(set(hits))))
            failures += 1
        else:
            add(rows, name, "PASS", detail)

    python_errors = []
    for path in python_files:
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except Exception as exc:
            python_errors.append("{}: {}".format(path.relative_to(root), exc))
    if python_errors:
        add(rows, "python_ast", "FAIL", " | ".join(python_errors))
        failures += 1
    else:
        add(rows, "python_ast", "PASS",
            "AST parsing passed for {} Python files.".format(len(python_files)))

    bash = shutil.which("bash")
    if not bash:
        add(rows, "shell_syntax", "FAIL", "bash is unavailable.")
        failures += 1
    else:
        shell_errors = []
        for path in shell_files:
            result = subprocess.run(
                [bash, "-n", str(path)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if result.returncode:
                shell_errors.append("{}: {}".format(
                    path.relative_to(root), result.stderr.strip()
                ))
        if shell_errors:
            add(rows, "shell_syntax", "FAIL", " | ".join(shell_errors))
            failures += 1
        else:
            add(rows, "shell_syntax", "PASS",
                "bash -n passed for {} shell/Slurm files.".format(len(shell_files)))

    rscript = shutil.which("Rscript")
    if r_files and not rscript:
        add(rows, "R_syntax", "FAIL", "Rscript is required but unavailable.")
        failures += 1
    elif not r_files:
        add(rows, "R_syntax", "NOT_APPLICABLE", "No R files are included.")
    else:
        r_errors = []
        for path in r_files:
            result = subprocess.run(
                [rscript, "-e", "parse(file=commandArgs(TRUE)[1])", str(path)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if result.returncode:
                r_errors.append("{}: {}".format(
                    path.relative_to(root), result.stderr.strip()
                ))
        if r_errors:
            add(rows, "R_syntax", "FAIL", " | ".join(r_errors))
            failures += 1
        else:
            add(rows, "R_syntax", "PASS",
                "Rscript parse passed for {} R files.".format(len(r_files)))

    if not (root / "LICENSE").exists():
        add(rows, "software_license", "WARN",
            "No license has been assigned; public visibility does not grant reuse rights.")
        warnings += 1
    else:
        add(rows, "software_license", "PASS", "A license file is present.")

    if not (root / "CITATION.cff").exists():
        add(rows, "citation_metadata", "WARN",
            "CITATION.cff is pending confirmation of final bibliographic details and authors.")
        warnings += 1
    else:
        add(rows, "citation_metadata", "PASS", "CITATION.cff is present.")

    if args.output:
        output = Path(args.output)
        if not output.is_absolute():
            output = root / output
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=["Check", "Status", "Detail"], delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(rows)

    print("Public repository validation")
    print("Repository: {}".format(root))
    for row in rows:
        print("{:<34} {:<14} {}".format(
            row["Check"], row["Status"], row["Detail"]
        ))
    print("\nChecks: {}; failures: {}; warnings: {}".format(
        len(rows), failures, warnings
    ))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
