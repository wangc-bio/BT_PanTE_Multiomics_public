# Module 01 — Genome and sample quality control

Module 01 contains reusable command-line examples for read-mapping summaries and assembly-quality assessment.

## Included workflows

- Read mapping and mapping-rate summaries from alignment logs.
- STAR log collection and summary-table generation.
- Assembly statistics with QUAST.
- k-mer counting with Meryl and assembly evaluation with Merqury.

## Inputs and configuration

The scripts use environment variables or command-line paths for input FASTQ files, genome assemblies, mapping logs, software executables and output directories. Review each script's configuration section before submission to a cluster. Create log/output directories before running Slurm jobs, and confirm that the selected software versions are available in the execution environment.

## Release scope

Raw reads, BAM/CRAM files, sample-level QC tables, assembly-specific contamination evidence and private evidence bundles are not included. The scripts are supplied as workflow examples; their successful syntax validation does not imply that QC was rerun on the original study data.
