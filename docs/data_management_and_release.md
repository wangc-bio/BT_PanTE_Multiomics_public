# Data management and release

## Contents distributed through GitHub

This repository is intended for source code, workflow documentation, blank metadata templates, environment examples, software-version records and lightweight validation tools. It does not distribute raw sequencing data, genome assemblies, large BAM/CRAM files, private sample-level QC tables or intermediate result matrices.

## Inputs and external data

For each analysis, record the exact input file, sample/library identifier, reference assembly, annotation release, file checksum and archive accession where applicable. Verify the mapping between BioProject, BioSample, experiment, run and library identifiers before using archive records. Do not infer a run-to-sample mapping from a tissue or genotype label alone.

Large input and output files should be deposited in an appropriate sequence archive or general-purpose data repository, subject to data ownership, consent, institutional policy and journal requirements. Add only verified accession numbers and DOI links to the repository and manuscript.

## Public-release checks

1. Confirm ownership and permission to distribute each file.
2. Remove personal server paths, credentials, tokens and unrelated personal information.
3. Exclude raw reads, alignments, large genome assemblies and private sample-level QC outputs.
4. Check that scripts use configurable paths and external manifests rather than private absolute paths.
5. Verify file integrity and generate a release-specific SHA-256 manifest.
6. Check that software versions and parameters agree with the manuscript Methods.
7. Distinguish source-level review, static syntax checks, synthetic tests and data-dependent reruns.
8. Record only verified data accessions, archive DOIs and software releases.

## Release limitations

A passing GitHub Actions run is a source-code and repository-hygiene check, not proof that the complete study analysis can be reproduced without the original input datasets and software environments. Large data files and archive identifiers should be linked only after their final versions and mappings have been verified.
