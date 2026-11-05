# Assembly-level QC scripts

This directory contains path-portable scripts that supplement the assembly-QC methods. The scripts use runtime input variables and do not embed cluster-specific absolute paths, usernames, or private project directories.

## FCS-GX contamination-candidate review

- [Contamination-screening report](../../../01_genome_and_sample_qc/assembly_qc/01B_contamination_screening/README.md) — interpretation, limitations, and released evidence.
- [Validated evidence bundle](../../../01_genome_and_sample_qc/assembly_qc/01B_contamination_screening/evidence_bundle.zip) — source FCS-GX reports, derived candidate evidence, alignments, methods, and checksums.
- [chr07M coordinate verification record](../../../01_genome_and_sample_qc/assembly_qc/01B_contamination_screening/results/chr07M_coordinate_verification.tsv) — nested-coordinate conversion and sequence-identity check (PASS).
- [Public script guide](01B_contamination_screening/README.md) — links to the path-portable chr07M coordinate verifier and parental-reference minimap2 workflow.

## Read mapping and BAM filtering

- [01_read_mapping_workflow_public.sh](01_read_mapping_workflow_public.sh) — path-portable templates for resequencing, HiFi, and Hi-C mapping. The Hi-C mode reproduces the recorded BWA `-5SP -t 24`, samblaster, samtools `-F 3340`, and HapHiC `filter_bam` steps. Supply the `filter_bam` executable through `HAPHIC_FILTER_BAM`.
- [02_summarize_star_logs.sh](02_summarize_star_logs.sh) — extracts per-sample STAR mapping categories and calculates a read-count-weighted aggregate from STAR `Log.final.out` files. The exact STAR alignment command-line parameters were not present in the supplied run records and are not reconstructed.

The current read-mapping summary is provisional. Do not treat its aggregate mapping rates as final until the pending mapping-statistics review and Supplementary Table S4 reconciliation are complete.

## QUAST assembly statistics

- [03_quast_assembly_statistics.sh](03_quast_assembly_statistics.sh) — path-portable QUAST v5.2.0 command template using five threads. The assembly statistics are already summarized in Supplementary Table S2.

## Merqury assembly evaluation

- [Merqury workflow README](merqury/README.md) — documents Meryl v1.4.1 k=21 counting and Merqury v1.4.1 evaluation.
- [01_meryl_k21_count.sh](merqury/01_meryl_k21_count.sh) — creates the k=21 Meryl database from PacBio HiFi reads.
- [02_merqury_assembly_evaluation.sh](merqury/02_merqury_assembly_evaluation.sh) — evaluates the final assembly using the Meryl database.

BUSCO v5.8.0 is documented in the assembly-QC methods note. No original BUSCO execution script was supplied, so no unverified command script has been reconstructed.
