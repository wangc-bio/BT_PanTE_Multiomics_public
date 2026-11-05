#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Generate a repository-wide SHA-256 manifest for a release snapshot."""

import argparse
import csv
import hashlib
from pathlib import Path

EXCLUDED_DIRS = {".git", "__pycache__", ".venv", ".conda"}
EXCLUDED_SUFFIXES = {".pyc"}
EXCLUDED_NAMES = {"module05.env"}


def digest(path, chunk_size=1024 * 1024):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            value.update(block)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--output", default="release_checksums.tsv")
    args = parser.parse_args()

    root = Path(args.repo_root).resolve()
    output = Path(args.output)
    if not output.is_absolute():
        output = root / output
    output = output.resolve()

    rows = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.resolve() == output:
            continue
        if any(part in EXCLUDED_DIRS for part in path.parts):
            continue
        if path.suffix in EXCLUDED_SUFFIXES or path.name in EXCLUDED_NAMES:
            continue
        if path.name.endswith(".env") and not path.name.endswith(".env.example"):
            continue
        rows.append({
            "File": str(path.relative_to(root)),
            "Size_bytes": path.stat().st_size,
            "SHA256": digest(path)
        })

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["File", "Size_bytes", "SHA256"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)
    print("Wrote {} file checksums to {}".format(len(rows), output))


if __name__ == "__main__":
    main()
