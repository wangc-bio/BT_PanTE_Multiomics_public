# Module 04B — Pan-TE transcription and epigenomic profiling

## Current scope

This release covers the **Pan-TE RNA-seq analysis** used to quantify Pan-TE transcription, test genotype effects within tissues, classify F1 transcriptional inheritance, summarize cross-tissue transgressive recurrence, and integrate transcription with genomic Pan-TE abundance.

WGBS processing, Pan-TE methylation quantification and ATAC-seq Pan-TE accessibility analyses are included below.

The RNA analysis includes five tissues (brain, gill, kidney, liver and muscle) from blunt snout bream (BSB), topmouth culter (TC) and the F1 hybrid (BT).

## Workflow

### 1. RNA-seq quality control and alignment

Public scripts:
- `01_RNAseq_QC_and_mapping/01_fastp_QC.sh`
- `01_RNAseq_QC_and_mapping/02_STAR_mapping.sh`
- `01_RNAseq_QC_and_mapping/03_build_metadata.py`
- `01_RNAseq_QC_and_mapping/04_collect_STAR_QC.sh`

Paired-end reads were quality-filtered with **fastp v0.23.2** using default filtering settings. Filtered reads were aligned to the corresponding genotype-specific reference genome with **STAR v2.7.10b**.

For downstream Pan-TE quantification, STAR retained multi-mapping reads by setting:

```text
--outFilterMultimapNmax 200
--winAnchorMultimapNmax 200
```

Other execution details such as server paths and thread counts are not method-defining and are omitted from the public description.

### 2. Pan-TE GTF construction

Script:
- `02_PanTE_GTF_preparation/build_PanTE_TEcount_GTF.py`

The refined Pan-TE membership table links each unified TE consensus to a Pan-TE cluster. RepeatMasker GFF records were parsed and TE loci carrying unified v5.2 TE identifiers were assigned to their corresponding Pan-TE clusters.

For TEcount-compatible GTF files:
- `gene_id` = Pan-TE cluster ID;
- `transcript_id` = locus-specific TE identifier;
- `family_id` = representative TE subclass;
- `class_id` = broad TE class.

Separate Pan-TE GTF files were generated for BSB and TC. The BT Pan-TE GTF combines loci from both BT-M and BT-P haplotypes.

The script also writes locus metadata and QC summaries for each genotype.

### 3. Pan-TE transcription quantification

Scripts:
- `03_PanTE_quantification/build_TEcount_jobs.py`
- `03_PanTE_quantification/run_TEcount_multi_array.sh`
- `03_PanTE_quantification/check_and_merge_TEcount.py`

Pan-TE transcription was quantified with **TEcount v2.2.4** in the TEtranscripts framework. The method-defining options were:

```text
--stranded no
--mode multi
```

Each sample was quantified using the corresponding gene GTF and genotype-specific Pan-TE GTF. TEcount outputs were checked for missing samples, unexpected Pan-TE IDs and non-integer count values before constructing the Pan-TE raw-count matrix.

A genomic-presence matrix was independently derived from the Pan-TE GTF annotations. This matrix records whether each Pan-TE cluster is represented in the BSB, TC and BT genomic annotations; it is **not** an RNA-detection matrix.

### 4. Pan-TE QC and differential transcription

Scripts:
- `04_Differential_transcription/01_RNA_PanTE_QC_DESeq2.R`
- `04_Differential_transcription/run_01_RNA_PanTE_QC_DESeq2.sh`
- `04_Differential_transcription/check_within_group_correlation.py`

Differential transcription was analyzed separately within each tissue using **DESeq2 v1.46.0** with genotype as the design factor (`~Group`).

Pan-TE clusters were retained within a tissue when raw counts were at least 10 in at least three biological samples.

Three contrasts were tested within each tissue:
- BSB versus TC;
- BT versus BSB;
- BT versus TC.

DESeq2-normalized counts were used to calculate genotype mean expression values. VST-transformed counts were used separately for PCA and sample-correlation QC.

The general differential-expression output includes a significance annotation requiring adjusted P < 0.05 and |log2FC| >= 1. This effect-size flag is **not** the criterion used for the main inheritance-class assignment described below.

### 5. Transcriptional inheritance classification

Scripts:
- `05_Transcriptional_inheritance/classify_PanTE_inheritance.py`
- `05_Transcriptional_inheritance/summarize_transgressive_recurrence.py`

Inheritance classification is performed only for Pan-TE clusters with genomic annotation in all three genotype groups (`Comparable_all3 = 1`).

Five principal classes are assigned:

- **Additive**: BT mean transcription lies strictly between the BSB and TC means. Statistical significance relative to the parents is not required.
- **BSB expression dominance**: BT is not significantly different from BSB but is significantly different from TC.
- **TC expression dominance**: BT is not significantly different from TC but is significantly different from BSB.
- **Above parental range**: BT mean transcription exceeds both parental means.
- **Below parental range**: BT mean transcription is below both parental means.

Exact numerical ties that cannot be assigned by these rules are labeled unresolved.

Pairwise significance for inheritance classification is based on DESeq2 Benjamini-Hochberg adjusted P values with FDR < 0.05.

#### Strict transgressive transcription

A Pan-TE event is classified as statistically supported transgressive transcription only when:

1. the BT mean lies outside the parental expression range; and
2. BT differs significantly from BSB (adjusted P < 0.05); and
3. BT differs significantly from TC (adjusted P < 0.05).

No additional fold-change threshold is applied to this main transgressive definition.

A separate **strong-transgressive** annotation is provided for events that also show |BT-parent log2FC| >= 1 in both parental comparisons. This stronger annotation is not used to define the main transgressive set.

The classifier also records a 12-pattern subdivision of the five principal inheritance classes.

### 6. Cross-tissue recurrence and directional consistency

Script:
- `05_Transcriptional_inheritance/summarize_transgressive_recurrence.py`

Strict transgressive tissue-Pan-TE events are summarized across tissues. For each Pan-TE cluster, the analysis records the number and identity of tissues showing statistically supported transgression and whether all transgressive events have the same direction.

For manuscript interpretation:
- recurrent = detected as strict transgressive in at least two tissues;
- directionally consistent = all recurrent transgressive events are either Up or all are Down.

The script additionally summarizes exact tissue combinations and recurrence separately for Up and Down events.

### 7. Genomic-abundance/RNA integration

Scripts:
- `06_Genomic_RNA_integration/integrate_genomic_abundance_RNA.py`
- `06_Genomic_RNA_integration/prioritize_transgressive_PanTE.py`

This downstream analysis integrates Pan-TE genomic abundance with RNA inheritance results. It is kept separate from the primary transcriptional-inheritance classification.

Genomic abundance is normalized as assigned TE bp per Mb of the corresponding genome. The analysis evaluates genomic-versus-RNA correlations and identifies transgressive RNA events that occur without a corresponding increase or decrease in genomic Pan-TE abundance.

The prioritization script further summarizes recurrence, direction consistency and availability of epigenomic tissues for downstream ATAC/WGBS integration.

## WGBS processing and Pan-TE DNA methylation

### 8. WGBS preprocessing

Scripts:
- `07_WGBS_processing/01_fastp_array.sh`
- `07_WGBS_processing/02a_bismark_index.sh`
- `07_WGBS_processing/02b_bismark_align.sh`
- `07_WGBS_processing/03_methylation_extract_ignoreR2_13.sh`
- `07_WGBS_processing/run_Mbias_QC.sh`
- `07_WGBS_processing/check_R2_stable_position.py`

Paired-end WGBS reads were quality filtered with **fastp v0.23.2** using paired-end adapter detection, polyG/polyX trimming, `-q 20`, `-u 30`, `-n 5` and `-l 30`. Bismark indices were prepared with **Bismark v0.24.2** using **Bowtie2 v2.3.5.1**. Reads were aligned in directional mode with `--score_min L,0,-0.6`, followed by paired-end duplicate removal with `deduplicate_bismark`.

CpG methylation calls were extracted with `bismark_methylation_extractor` using `--paired-end --comprehensive --merge_non_CpG --bedGraph --counts --report --gzip --ignore_r2 13`. The 13-bp R2 trimming was selected after M-bias inspection of CpG methylation across read positions.

### 9. Pan-TE methylation quantification

Scripts:
- `08_PanTE_mCG_quantification/check_WGBS_PanTE_coordinates.sh`
- `08_PanTE_mCG_quantification/prepare_PanTE_WGBS_intervals.py`
- `08_PanTE_mCG_quantification/quantify_PanTE_weighted_mCG.py`
- `08_PanTE_mCG_quantification/run_PanTE_weighted_mCG_array.sh`
- `08_PanTE_mCG_quantification/summarize_PanTE_mCG_coverage.py`
- `08_PanTE_mCG_quantification/summarize_PanTE_mCG_QC.py`

Reference/Pan-TE coordinate compatibility is checked before quantification. Pan-TE intervals are derived from the genotype-specific Pan-TE annotations; BT-M and BT-P are retained as separate haplotypic states. CpG calls overlapping more than one Pan-TE cluster are discarded, whereas overlaps among multiple intervals belonging to the same Pan-TE cluster remain assignable.

For each sample and Pan-TE state, weighted methylation is calculated as `mC / (mC + uC)`. The main filter requires at least 5 informative CpGs and mean CpG depth >= 3 in each biological replicate, with all three replicates required. Methylated and unmethylated CpG counts are then pooled across the three replicates to calculate the state-level pooled weighted mCG. Relaxed (>=3 CpGs, mean depth >=2) and stringent (>=10 CpGs, mean depth >=5) schemes are retained as sensitivity analyses.

### 10. Parent/haplotype methylation comparisons

Scripts:
- `09_PanTE_mCG_comparison/build_PanTE_mCG_comparison_matrix.py`
- `09_PanTE_mCG_comparison/analyze_parental_mCG_difference_retention.py`

Matched comparisons pair BSB with BT-M and TC with BT-P. Parental methylation contrast is defined as BSB minus TC, and the corresponding F1 haplotypic contrast as BT-M minus BT-P. Their correspondence is evaluated using Spearman correlation. Directional retention is additionally summarized for absolute parental mCG-difference thresholds of 0.05, 0.10 and 0.15.


## ATAC-seq processing and Pan-TE chromatin accessibility

### 11. ATAC-seq alignment and BAM quality control

Scripts:
- `10_ATAC_processing_and_PanTE_accessibility/01_ATAC_alignment_and_QC/02_ATAC_array.sh`
- `10_ATAC_processing_and_PanTE_accessibility/01_ATAC_alignment_and_QC/run_ATAC_clean_BAM_QC.sh`
- `10_ATAC_processing_and_PanTE_accessibility/01_ATAC_alignment_and_QC/run_ATAC_mapping_QC.sh`
- `10_ATAC_processing_and_PanTE_accessibility/01_ATAC_alignment_and_QC/run_ATAC_mapping_QC_array.sh`
- `10_ATAC_processing_and_PanTE_accessibility/01_ATAC_alignment_and_QC/collect_ATAC_mapping_QC.sh`
- `10_ATAC_processing_and_PanTE_accessibility/01_ATAC_alignment_and_QC/collect_ATAC_mapping_QC.py`

Paired-end ATAC-seq reads were aligned to the corresponding reference genomes using BWA-MEM v0.7.17-r1188. PCR duplicates were removed with Picard MarkDuplicates v2.25.5 using `REMOVE_DUPLICATES=true`. Clean alignments were filtered for mapping quality >=30 and proper pairing, with unmapped reads, unmapped mates, secondary alignments, PCR duplicates and QC-failed reads excluded. Supplementary alignments were additionally excluded during the downstream Pan-TE fragment-counting step. The Pan-TE analysis was restricted to the main nuclear chromosomes: 24 chromosomes for each parental genome and the 24 M-derived plus 24 P-derived chromosomes for BT.

The released QC scripts collect mapping statistics and support BAM-level checks. Sample-specific metadata construction and chromosome-length audit inputs are supplied externally and are not distributed with this repository.

### 12. Pan-TE uniquely assignable sequence

Script:
- `10_ATAC_processing_and_PanTE_accessibility/02_PanTE_assignable_segments/build_ATAC_PanTE_assignable_segments.py`

The input interval tables (`BSB.PanTE.WGBS_intervals.tsv`, `TC.PanTE.WGBS_intervals.tsv` and `BT.PanTE.WGBS_intervals.tsv`) are coordinate tables derived from the original Pan-TE cluster genomic annotations. They are **not** filtered by WGBS CpG counts, methylation depth or replicate availability. Thus, ATAC assignable intervals are constructed independently of WGBS methylation-coverage filtering.

Within each main chromosome, overlapping Pan-TE intervals were partitioned into non-overlapping segments. A segment was retained only when it was covered by exactly one distinct Pan-TE cluster. Overlaps between multiple loci belonging to the same Pan-TE cluster remain assignable; segments overlapping different Pan-TE clusters are excluded. BT M- and P-haplotype chromosomes are assigned separately to BTM and BTP.

### 13. Pan-TE fragment counting and accessibility normalization

Scripts:
- `10_ATAC_processing_and_PanTE_accessibility/03_PanTE_fragment_quantification/preflight_ATAC_PanTE_count.sh`
- `10_ATAC_processing_and_PanTE_accessibility/03_PanTE_fragment_quantification/run_ATAC_PanTE_count_array.sh`
- `10_ATAC_processing_and_PanTE_accessibility/03_PanTE_fragment_quantification/normalize_ATAC_PanTE_counts.py`
- `10_ATAC_processing_and_PanTE_accessibility/03_PanTE_fragment_quantification/summarize_ATAC_PanTE_QC.py`

Paired-end fragments, rather than individual read alignments, were counted. Only proper-pair, primary, mapped, QC-passing, non-duplicate and non-supplementary alignments were considered; MAPQ >=30 was retained from the clean BAM filtering. One fragment record was generated per pair using the positive template-length record. A fragment was assigned to a Pan-TE cluster only when its overlap involved exactly one distinct Pan-TE cluster; fragments overlapping multiple distinct clusters were excluded. Fragment counts were normalized by the assignable sequence length and the total retained fragments in the corresponding genomic state:

```text
FPKM = fragment_count × 10^9 / (assignable_bp × total_state_fragments)
accessibility = log2(FPKM + 1)
```

For BT, BTM and BTP were normalized independently using their corresponding state-specific retained-fragment totals.

### 14. Filtering and parental/haplotype accessibility comparisons

Scripts:
- `10_ATAC_processing_and_PanTE_accessibility/04_Accessibility_retention_analysis/analyze_ATAC_parental_state_retention.py`
- `10_ATAC_processing_and_PanTE_accessibility/04_Accessibility_retention_analysis/analyze_ATAC_filter_sensitivity_retention.py`

The primary analysis required >=500 bp of uniquely assignable Pan-TE sequence and >=3 fragments in **each** of the three biological replicates. The mean accessibility state was calculated as the mean of replicate-level `log2(FPKM + 1)` values. Only Pan-TE clusters passing the primary criteria in all four genomic states (BSB, TC, BTM and BTP) were retained for four-state comparisons within each tissue.

Matched parent-to-haplotype comparisons were BSB versus BTM and TC versus BTP. Parental contrast was calculated as BSB minus TC, and the corresponding F1 haplotypic contrast as BTM minus BTP. Spearman correlation was used to assess the correspondence between parental and F1 haplotypic contrasts.

Filter sensitivity was assessed using relaxed (>=100 bp and >=1 fragment per replicate), primary (>=500 bp and >=3 fragments per replicate) and stringent (>=1,000 bp and >=5 fragments per replicate) criteria. The retention classification additionally evaluated parental-contrast effect-size cutoffs of 0.5, 1.0 and 1.5 on the log2 accessibility scale. These are effect-size thresholds, not statistical-significance tests.

### 15. Interpretation limitation

No external spike-in normalization was used. Accordingly, these analyses describe relative Pan-TE accessibility states after within-genomic-state library-size and assignable-length normalization; they do not establish genome-wide absolute chromatin opening or closing.

## Software versions reported for the RNA, WGBS and ATAC analyses

- fastp v0.23.2
- STAR v2.7.10b
- TEcount / TEtranscripts v2.2.4
- DESeq2 v1.46.0
- Bismark v0.24.2
- Bowtie2 v2.3.5.1
- BWA-MEM v0.7.17-r1188
- Picard v2.25.5

## Interpretation notes

- Multi-mapping reads were intentionally retained for Pan-TE quantification.
- Genomic presence is based on reference TE annotation and should not be interpreted as RNA expression.
- The main additive class is a positional expression category: the F1 mean lies between parental means. It is not defined by equality to the mid-parent value.
- The main strict-transgressive definition uses parental-range position plus significance versus both parents and does not use an additional fold-change cutoff.
- The |log2FC| >= 1 threshold is used only for the separate strong-transgressive annotation and for general differential-expression effect-size summaries.
- RSEM gene/transcript quantification is not part of the Pan-TE transcription workflow documented here.
- ATAC assignable intervals originate from the original Pan-TE genomic annotations and are independent of WGBS CpG-coverage filters.
- ATAC accessibility values are relative and should not be interpreted as genome-wide absolute chromatin changes without external spike-in normalization.

## Source-release notes

The released scripts use configurable project roots or executable variables where applicable. Cluster-specific environment activation commands and absolute software paths have been removed; users must make the required tools available in PATH and supply local input/output paths. Sample-level metadata and private QC records are not distributed.

Results already reported in the manuscript or supplementary materials are not duplicated in this code directory.


## Cross-assay and multi-omics statistical analyses

The assay-specific WGBS and ATAC processing/retention scripts above generate the inputs used for manuscript Methods 2.11 and 2.12. Cross-assay Figure 5 summary preparation, pooled OLS/PanTE_ID cluster-bootstrap analysis, multi-omics master-matrix construction, genomic-divergence models, abundance-adjusted logistic regression and transgressive epigenetic analyses are organized in Module 05. These integrated analyses are kept separate from assay-level quantification to avoid duplicating Module 04B workflows.
