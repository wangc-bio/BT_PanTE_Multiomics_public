# Analysis workflow

This document summarizes the six numbered computational modules (01–06), with Module 04 divided into 04A and 04B, and the evidence each module is intended to produce. The repository contains source workflows and parameter documentation; it does not bundle the raw sequencing data or all intermediate analysis outputs.

## Module 01 — Genome and sample quality control

This module covers read-level quality control, mapping summaries and assembly-level quality assessment. The scripts provide examples for collecting STAR mapping statistics, generating assembly statistics and evaluating assembly k-mer completeness with Meryl/Merqury. Input FASTQ files, reference assemblies, alignments and sample manifests are external inputs. The repository does not distribute sample-level QC tables or contamination-screening evidence bundles.

## Module 02 — Phased F1 assembly and parental-origin validation

The module documents whole-genome alignment, synteny analysis, chromosome-wise phylogenetic validation, parent-specific k-mer mapping, parental resequencing-depth validation and homologous-gene identification. The scripts are parameterized through working-directory and input-file variables. Assembly FASTA/GFF files, parental read data, k-mer databases and external phylogenetic inputs must be supplied by the user. The scripts do not imply that every step can be rerun without these data and software dependencies.

## Module 03 — Pan-TE catalogue construction

The module covers de novo repeat discovery and consensus curation, structural and homology evidence, TE classification, unified repeat-library construction, sequence-based Pan-TE clustering and genomic occupancy quantification.

The manuscript-level comparison uses four focal assemblies: BSB, TC, BT-M and BT-P. BT-M and BT-P are the phased F1 haplotypes. Runtime manifests are used to specify the inputs. The BC1 backcross assemblies are outside this four-assembly Pan-TE analysis.

Pan-TE clusters are operational sequence-homology units constructed from consensus sequences. They should not be equated with genomic insertion loci or population insertion frequencies. The Module 03 synthetic test validates the classification scripts with synthetic records; it does not establish the biological accuracy of production classifications or replace a run using the study's original data.

## Module 04A — TE sequence-divergence landscapes

This module documents RepeatMasker alignment generation using the unified repeat library, Kimura-divergence parsing, divergence-bin coverage, genome-size normalization, parental/F1 landscape comparisons, non-zero-bin sensitivity analysis and LTR contrast analysis. Kimura divergence is interpreted as sequence divergence from consensus, not as an absolute insertion-time estimate. The workflow does not perform insertion-site presence/absence genotyping or population-frequency estimation.

## Module 04B — TE transcription and epigenomic profiling

This module contains RNA-seq quality control and mapping, Pan-TE GTF preparation, TEcount quantification, differential transcription and transcriptional-inheritance classification, WGBS processing and weighted methylation quantification, ATAC-seq processing and Pan-TE accessibility quantification. Assay-specific filtering and QC steps are kept separate.

The metadata files used for the original sequencing libraries are not distributed. Users must construct manifests for their own samples and preserve the distinction among individual, tissue, library and sequencing run. Because no external ATAC spike-in normalization was used in the study, ATAC-derived Pan-TE accessibility should be interpreted as a relative measure.

## Module 05 — Integrated statistical modelling

This module contains the analysis scripts for the cross-assay master matrix, parental-state correspondence, genomic-divergence/epigenetic-deviation associations, abundance-adjusted transcriptional comparisons, transgressive epigenetic-deviation analyses, directional concordance and background-centered permutation tests. It corresponds to the integrated analyses described in Methods 2.11 and 2.12.

Large input matrices are external. The environment template and script arguments must be populated with the exact released input files. The source code and static validation alone do not constitute an independent rerun of the complete multi-omics analysis.

## Module 06 — Validation and reproducibility

This module provides a repository validator, SHA-256 manifest generator, software-version inventory, release checklist and validation status record. GitHub Actions checks syntax and repository hygiene and runs the synthetic Module 03 classification test. The workflow artifact associated with the final release commit should be retained as the formal validation record.

## Scope of validation

Static syntax checks establish that scripts can be parsed; synthetic tests establish that tested code paths operate on the supplied fixtures. Neither establishes numerical correctness on the study data, successful execution on a particular HPC system, or biological validity. Data-dependent reproduction requires the exact input files, software versions, parameters and reference assemblies used for the manuscript.
