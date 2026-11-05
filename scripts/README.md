# Analysis scripts

This directory contains the executable workflows for Modules 01–06.

## Module map

- `01_genome_and_sample_qc/` — read mapping summaries, assembly statistics and Merqury evaluation.
- `02_phased_F1_assembly_and_parental_origin/` — whole-genome alignment, synteny, phylogenetic validation, parent-specific k-mer and parental-depth analyses.
- `03_pan_TE_catalogue/` — TE discovery/classification, unified-library construction, Pan-TE clustering and genomic occupancy.
- `04A_TE_insertion_inheritance_and_genomic_divergence/` — Kimura-divergence calculation and landscape comparisons.
- `04B_TE_transcription_and_epigenomic_profiling/` — RNA-seq, WGBS and ATAC-seq processing and Pan-TE quantification.
- `05_integrated_statistical_modelling/` — cross-assay integration and Methods 2.11–2.12 statistical analyses.
- `06_validation_and_reproducibility/` — release validation, software records and checksums.

## Execution conventions

1. Read the README in the relevant module and subdirectory before running a script.
2. Supply reference files, sample manifests and output directories explicitly.
3. Check the software versions in `scripts/06_validation_and_reproducibility/software_versions.tsv`.
4. Create log directories before submitting Slurm jobs.
5. Review the script's input validation and output naming before launching an array job.
6. Preserve the distinction among biological individuals, tissues, libraries and sequencing runs.

The scripts are source workflows, not a substitute for the study's external input datasets. Some modules require large reference assemblies, TE libraries, alignment files or assay-specific matrices that are not distributed in this repository.
