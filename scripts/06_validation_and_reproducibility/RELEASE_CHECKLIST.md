# Public-release checklist

## Repository contents

- [ ] Include all numbered module directories (01–06; Module 04A and 04B), core method documentation and reusable source scripts.
- [ ] Exclude raw sequencing reads, genome assemblies, BAM/CRAM files and large intermediate matrices.
- [ ] Exclude private sample-level QC records, individual sample identifiers and internal evidence bundles.
- [ ] Confirm that every distributed file is owned by the project or has permission for redistribution.
- [ ] Confirm that no credentials, tokens, private absolute paths or unrelated personal information remain.

## Workflow and methods

- [ ] Confirm that the four-assembly Pan-TE analysis scope is BSB, TC, BT-M and BT-P.
- [ ] Check that module documentation agrees with the manuscript Methods.
- [ ] Check that scripts use configurable input/output paths and runtime manifests.
- [ ] Record exact software versions where supported by execution records.
- [ ] Keep unverified software versions explicitly marked as not recorded or not established.
- [ ] Distinguish source-level review, static validation, synthetic tests and data-dependent reruns.

## Validation

- [ ] Run the Module 03 synthetic classification test.
- [ ] Run Python AST validation for all released Python scripts.
- [ ] Run Bash syntax validation for all shell and Slurm scripts.
- [ ] Parse R scripts with Rscript.
- [ ] Run the Module 06 GitHub Actions workflow on the intended release commit.
- [ ] Review every FAIL, WARN and SKIP item; do not relabel skipped checks as passed.
- [ ] Generate and retain the SHA-256 manifest for the release snapshot.

## Release metadata

- [ ] Confirm the manuscript title, author list and citation details.
- [x] Add CITATION.cff for the software repository using the confirmed manuscript author order.
- [ ] Add the manuscript's preferred-citation metadata after its final title, journal, publication year and DOI are confirmed.
- [ ] Select and add a software license only after confirming ownership and third-party redistribution terms.
- [ ] Add only verified sequence-archive accessions and dataset DOIs.
- [ ] Tag the exact commit whose validation artifact and checksum manifest have been retained.

## Release decision

A public repository may be visible before a formal software release, but it should not be described as fully reproducible until the required input datasets, software environments and data-dependent analyses have been independently verified.
