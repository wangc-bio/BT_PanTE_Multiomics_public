# Module 06 — Validation and reproducibility

Module 06 defines the source-release checks used before publishing this repository. It separates repository hygiene and static syntax checks from data-dependent reruns and biological validation.

## Included tools

- `validate_release.py`: checks required documentation, raw/oversized files, local environment files, recognized personal absolute paths, Python syntax, Bash/Slurm syntax and R syntax when Rscript is available.
- `generate_sha256_manifest.py`: generates a SHA-256 manifest for the repository snapshot.
- `software_versions.tsv`: records software versions and indicates where a version was not established.
- Root-level `CITATION.cff`: provides the software citation and confirmed author order; it does not yet encode an unconfirmed manuscript title or DOI.
- `RELEASE_CHECKLIST.md`: records the checks required before a tagged release.

## Run locally

From the repository root:

```bash
python scripts/06_validation_and_reproducibility/validate_release.py \
  --repo-root . \
  --output release_validation_report.tsv

python scripts/06_validation_and_reproducibility/generate_sha256_manifest.py \
  --repo-root . \
  --output release_checksums.tsv
```

The validator requires Python 3. It uses `bash -n` for shell syntax and uses Rscript to parse R files when Rscript is installed. A missing Rscript is reported as SKIP rather than as a passing R check.

## GitHub Actions

The `Module 06 public-release validation` workflow runs on pushes and pull requests to `main`, and can also be started manually. It produces a validation report and SHA-256 manifest as a workflow artifact. Retain the artifact associated with the exact release commit.

## Interpretation

A passing validator establishes repository-level checks only. It does not rerun the biological analyses, validate numerical results against the original study data, confirm that all software environments are reproducible, or establish biological causality. The Module 03 synthetic test exercises the classification workflow using synthetic records only.
