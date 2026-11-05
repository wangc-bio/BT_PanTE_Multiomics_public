# Public-release validation record

Validation is performed for each commit by the GitHub Actions workflow **Module 06 public-release validation**.

## Authoritative record

For a formal release, open the repository's [GitHub Actions runs](../../actions/workflows/module06-release-validation.yml) and select the successful run whose **head SHA exactly matches the release commit**. Retain the artifact named `public-release-validation` from that run. It contains:

- `release_validation_report.tsv` — check-by-check results;
- `release_checksums.tsv` — SHA-256 checksums and file sizes for the repository snapshot.

Do not use artifacts from intermediate migration commits. The report and checksum manifest must correspond to the exact commit being cited or tagged.

## Checks performed

The Module 06 workflow:

1. verifies the required documentation;
2. scans public text files for personal/cluster-specific paths, credential patterns and individual sample IDs;
3. checks that raw sequencing/alignment files, local environment files and oversized files are absent;
4. parses every Python source file;
5. runs `bash -n` on shell and Slurm scripts;
6. parses R scripts with Rscript;
7. generates the repository SHA-256 manifest.

The separate Module 03 workflow compiles the classification scripts and runs the synthetic end-to-end classification test.

## Scope and limitations

These are source-release and software-regression checks. They do not constitute a data-dependent rerun of the complete manuscript analyses, validation of all numerical results against the original study data, or experimental validation. The external input datasets and exact execution environments are still required for full reproduction.

The software license and formal citation metadata remain pending confirmation of ownership, final bibliographic details and the author list.
