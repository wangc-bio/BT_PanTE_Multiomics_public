# Module 04B Pan-TE RNA source-script inventory

Updated: 2026-11-04.

| Stage | Public source files | Verified method-defining settings |
|---|---|---|
| RNA-seq QC | 01_fastp_QC.sh | fastp v0.23.2; default filtering settings |
| RNA-seq alignment | 02_STAR_mapping.sh | STAR v2.7.10b; --outFilterMultimapNmax 200; --winAnchorMultimapNmax 200 |
| Metadata/QC collection | 03_build_metadata.py; 04_collect_STAR_QC.sh | genotype/tissue metadata and STAR mapping-QC extraction |
| Pan-TE GTF preparation | build_PanTE_TEcount_GTF.py | unified TE ID to Pan-TE cluster mapping; genotype-specific TE GTFs; BT merges BT-M + BT-P loci |
| Pan-TE quantification | build_TEcount_jobs.py; run_TEcount_multi_array.sh; check_and_merge_TEcount.py | TEcount v2.2.4; --stranded no; --mode multi; integer-count and missing-output QC |
| Differential transcription | 01_RNA_PanTE_QC_DESeq2.R; run_01_RNA_PanTE_QC_DESeq2.sh; check_within_group_correlation.py | tissue-wise DESeq2 v1.46.0; design ~Group; count >=10 in >=3 samples; BH FDR |
| Inheritance classification | classify_PanTE_inheritance.py | Comparable_all3; Additive positional definition; parental-expression dominance; above/below parental range; strict transgressive requires FDR <0.05 vs both parents |
| Cross-tissue recurrence | summarize_transgressive_recurrence.py | recurrence count by tissue; direction consistency; recurrent manuscript definition >=2 tissues |
| Genomic-RNA integration | integrate_genomic_abundance_RNA.py; prioritize_transgressive_PanTE.py | genomic abundance per Mb; transgressive genomic context; downstream candidate prioritization |

## Deliberate exclusions from the Pan-TE RNA method

- RSEM expression quantification is not used for Pan-TE counts and is not part of this module's Pan-TE quantification workflow.
- Server paths, Slurm resources and thread counts are execution details and are not treated as scientific method parameters.
- ATAC-seq and WGBS are documented in their own sections of Module 04B.


## WGBS source-script inventory

| Stage | Public source files | Verified method-defining settings |
|---|---|---|
| WGBS read QC | 07_WGBS_processing/01_fastp_array.sh | fastp v0.23.2; paired-end adapter detection; polyG/polyX trimming; -q 20; -u 30; -n 5; -l 30 |
| Bismark index | 07_WGBS_processing/02a_bismark_index.sh | Bismark v0.24.2; Bowtie2 v2.3.5.1 |
| WGBS alignment | 07_WGBS_processing/02b_bismark_align.sh | directional mode; --score_min L,0,-0.6; paired-end duplicate removal |
| M-bias / R2 trimming QC | 07_WGBS_processing/run_Mbias_QC.sh; check_R2_stable_position.py | CpG positional M-bias assessment; final --ignore_r2 13 |
| Methylation extraction | 07_WGBS_processing/03_methylation_extract_ignoreR2_13.sh | --paired-end; --comprehensive; --merge_non_CpG; --bedGraph; --counts; --report; --gzip; --ignore_r2 13 |
| Pan-TE interval preparation | 08_PanTE_mCG_quantification/check_WGBS_PanTE_coordinates.sh; prepare_PanTE_WGBS_intervals.py | coordinate compatibility QC; BT-M/BT-P retained separately; overlapping intervals for the same Pan-TE merged |
| Pan-TE mCG quantification | 08_PanTE_mCG_quantification/quantify_PanTE_weighted_mCG.py; run_PanTE_weighted_mCG_array.sh | ambiguous cross-cluster CpGs excluded; weighted mCG = mC/(mC+uC) |
| Coverage/QC summaries | summarize_PanTE_mCG_coverage.py; summarize_PanTE_mCG_QC.py | sample-level coverage and assignment summaries |
| Replicate filtering / pooling | 09_PanTE_mCG_comparison/build_PanTE_mCG_comparison_matrix.py | main: >=5 CpGs, mean depth >=3, 3/3 replicates; relaxed and stringent sensitivity schemes |
| Parent/haplotype retention | 09_PanTE_mCG_comparison/analyze_parental_mCG_difference_retention.py | BSB-BTM and TC-BTP matched comparisons; parental versus F1 haplotypic mCG contrasts; thresholds 0.05/0.10/0.15 |


## ATAC-seq source-script inventory

| Stage | Public source files | Method-defining details |
|---|---|---|
| Alignment and duplicate removal | 10_ATAC_processing_and_PanTE_accessibility/01_ATAC_alignment_and_QC/02_ATAC_array.sh | BWA-MEM; Picard MarkDuplicates; duplicate removal |
| Metadata and reference audit | User-supplied sample manifest and chromosome-length metadata | Sample-specific metadata mappings and BAM audit outputs are external inputs and are not distributed |
| Mapping QC | run_ATAC_mapping_QC.sh; run_ATAC_mapping_QC_array.sh; collect_ATAC_mapping_QC.sh; collect_ATAC_mapping_QC.py | MAPQ distribution and mapping summary |
| Clean BAM QC | run_ATAC_clean_BAM_QC.sh | clean BAM QC |
| Assignable Pan-TE segments | 02_PanTE_assignable_segments/build_ATAC_PanTE_assignable_segments.py | original Pan-TE genomic intervals; main chromosomes; unique cluster-level assignability; BT-M/BT-P separated |
| Fragment-count preflight | 03_PanTE_fragment_quantification/preflight_ATAC_PanTE_count.sh | coordinate and input checks before counting |
| Fragment quantification | 03_PanTE_fragment_quantification/run_ATAC_PanTE_count_array.sh | proper pairs; primary; mapped; QC-pass; non-duplicate; non-supplementary; one fragment per pair; ambiguous multi-cluster fragments excluded |
| Accessibility normalization | 03_PanTE_fragment_quantification/normalize_ATAC_PanTE_counts.py | FPKM = fragment_count × 10^9 / (assignable_bp × total state fragments); log2(FPKM+1) |
| Pan-TE accessibility QC | 03_PanTE_fragment_quantification/summarize_ATAC_PanTE_QC.py | main-chromosome fraction; unique assignment; assignable-bp distribution; replicate correlations; filter sensitivity |
| Parent/haplotype retention | 04_Accessibility_retention_analysis/analyze_ATAC_parental_state_retention.py | primary >=500 bp and >=3 fragments in each of 3 replicates; four-state common Pan-TE set; matched BSB-BTM and TC-BTP; Spearman contrasts |
| Filter sensitivity | 04_Accessibility_retention_analysis/analyze_ATAC_filter_sensitivity_retention.py | relaxed 100 bp/1 fragment; main 500 bp/3 fragments; stringent 1,000 bp/5 fragments |
| Retention effect-size sensitivity | analyze_ATAC_parental_state_retention.py | parental contrast cutoffs 0.5, 1.0 and 1.5 on log2 accessibility scale; effect-size categories, not significance tests |

## ATAC interval provenance

The three interval inputs named `BSB.PanTE.WGBS_intervals.tsv`, `TC.PanTE.WGBS_intervals.tsv` and `BT.PanTE.WGBS_intervals.tsv` are derived from the original Pan-TE cluster genomic annotations. They are not prefiltered by WGBS CpG counts, methylation depth or replicate availability. ATAC assignable-segment construction is therefore independent of WGBS methylation-coverage filtering.

## ATAC interpretation limitation

No external spike-in normalization was used. The reported Pan-TE accessibility values are relative within each genomic state and do not support claims of genome-wide absolute chromatin opening or closing.
