# Module 05 — Integrated statistical modelling

## Scope and boundary

This module contains the cross-assay and multi-omics statistical analyses described in manuscript Methods 2.11 and 2.12. It starts from processed Pan-TE matrices produced by Module 04B and genomic-occupancy/RNA tables produced by Modules 03 and 04B. It does not replace assay-specific read processing, quantification or QC workflows in Module 04B.

The statistical scripts are based on the source analyses used for the manuscript, with machine-specific absolute paths replaced by environment variables. The analytical calculations and manuscript-level settings are retained. One release-time change is documented: prepare_Figure5_summary.py now recalculates WGBS correlations from the six source matrices rather than embedding frozen summary values; its outputs were cross-checked against the previously recorded Figure 5 tables. Large production matrices and raw sequencing files remain external inputs.

## Directory structure

~~~text
05_integrated_statistical_modelling/
├── 00_build_multiomics_master/
├── 01_parental_state_correspondence/
├── 02_genomic_divergence_and_epigenetic_deviation/
├── 03_parental_divergence_and_transgressive_transcription/
├── 04_transgressive_epigenetic_deviation/
├── module05.env.example
├── SOURCE_SCRIPT_INVENTORY.md
└── README.md
~~~

## Methods 2.11 — Matched-parent state correspondence and contrast remodeling

### Inputs from Module 04B

- Six main-filter WGBS matrices: main.{gill,kidney,muscle}.{BSB_to_BTM,TC_to_BTP}.PanTE_mCG.tsv.
- Three WGBS per-tissue parental-contrast tables: {gill,kidney,muscle}.parental_mCG_difference_retention.tsv.
- ATAC four-state matrix: ATAC_parental_state_retention.four_state_matrix.tsv.

The WGBS matrices contain pooled weighted mCG values for the matched parental/haplotype pairs. The ATAC matrix contains BSB, TC, BTM and BTP values on the log2(FPKM + 1) scale.

### Scripts

1. prepare_Figure5_summary.py recalculates BSB–BT-M and TC–BT-P Spearman correlations separately for WGBS and ATAC in each tissue, and recalculates parental contrast (BSB − TC) versus F1 haplotypic contrast (BT-M − BT-P). It writes Figure5_absolute_state_retention.tsv, Figure5_parental_contrast_retention.tsv and Figure5_layer_summary.tsv. It uses the six WGBS matrices directly rather than embedding fixed correlation values.
2. plot_Figure5_retention_landscape_v5.py generates the Figure 5 retention landscape and pooled contrast-transmission panels. It fits descriptive unweighted OLS relationships to pooled Pan-TE × tissue observations, using 1,000 bootstrap replicates and percentile 95% confidence intervals. PanTE_ID clusters are resampled as units; row-level resampling is the fallback only when a usable PanTE_ID column with at least 100 unique clusters is unavailable. WGBS seed: 20260908; ATAC seed: 20260909. The existing analysis output records PanTE_ID cluster resampling for both assays.

The WGBS and ATAC per-assay retention calculations remain in Module 04B. They are not duplicated here.

## Methods 2.12 — Multi-omics integration

### 1. Build the multi-omics master matrix

build_Result6_multiomics_master.py joins resolved Pan-TE genomic occupancy, tissue-specific RNA inheritance, WGBS parental-contrast tables and the ATAC four-state matrix. It writes all integrated tissue rows, the RNA-testable subset, the complete RNA/WGBS/ATAC subset and input/QC summaries. Downstream models use the RNA-testable matrix unless otherwise specified.

### 2. Parental genomic divergence and epigenetic deviation

analyze_genomic_divergence_vs_epigenetic_deviation.py defines D_genomic = |BSB_bp_per_Mb − TC_bp_per_Mb| / (BSB_bp_per_Mb + TC_bp_per_Mb), uses matched-parent epigenetic deviation (D_mCG_mean or D_ATAC_mean) as the response, and calculates tissue-wise Spearman correlations for assay-available and parent-informative backgrounds (Parent_max_bp >= 10,000). BH correction is applied to the 12 primary tissue/background/assay tests and the additional corrections retained in the source analysis.

### 3. TE-class and abundance-controlled sensitivity

analyze_TEclass_abundance_controlled_sensitivity.py evaluates DNA, LTR, LINE, SINE and RC classes using partial Spearman correlation controlling for log10(BSB_bp_per_Mb + TC_bp_per_Mb + 1). It requires at least 30 observations for a class-level test; SINE results below this threshold are not interpreted as supported tests. Raw/partial correlations, sample sizes and BH-adjusted P values within tissue and assay are reported.

### 4. Parent-divergence and strict-transgressive transcription

reanalyze_Figure6C_abundance_adjusted.py restricts the main analysis to RNA-testable, parent-informative Pan-TE × tissue observations (Parent_max_bp >= 10,000). It uses binomial GLM logistic regression with strict-transgressive status as the response and parent-divergent status plus total parental genomic abundance as predictors. Tissue-specific models use HC0 robust standard errors; the pooled model includes tissue and PanTE_ID-clustered robust standard errors. The abundance-decile model is retained as a sensitivity analysis, and BH correction is applied to tissue-wise model P values.

### 5. Strict-transgressive epigenetic deviation and direction

- analyze_transgressive_epigenetic_deviation.py: two-sided Mann–Whitney U comparisons of matched-parent mCG/ATAC deviations between strict-transgressive and non-transgressive Pan-TEs, with Cliff's delta effect sizes and BH correction.
- analyze_transgressive_directional_epigenetic_concordance.py: signed matched-parent shifts for BT-M versus BSB and BT-P versus TC; background-centers shifts on non-transgressive Pan-TEs; evaluates concordance using Wilcoxon signed-rank and exact binomial tests. Both-haplotype concordance requires both haplotypes to agree with the assay-specific direction rule.
- analyze_background_centered_directional_permutation.py: 10,000 within-tissue/assay permutations of RNA Up/Down labels, preserving observed epigenetic shifts and Up/Down counts. Seed: 20260909.

For concordance scoring, ATAC shifts are evaluated in the same direction as RNA change, whereas mCG shifts are evaluated in the opposite direction. These analyses quantify association and directional agreement; they do not establish causation.

## Execution

1. Enter the Module 05 directory, copy module05.env.example to module05.env, edit paths and source it:

~~~bash
cd scripts/05_integrated_statistical_modelling
cp module05.env.example module05.env
source module05.env
~~~

2. Build the multi-omics master matrix:

~~~bash
python 00_build_multiomics_master/build_Result6_multiomics_master.py
~~~

3. Prepare and plot Figure 5:

~~~bash
export PANTE_OUTPUT_DIR="$PANTE_FIG5_SUMMARY_DIR"
python 01_parental_state_correspondence/prepare_Figure5_summary.py
export PANTE_OUTPUT_DIR="$PANTE_INTEGRATION_DIR/03_Figure5_retention_landscape_v5"
python 01_parental_state_correspondence/plot_Figure5_retention_landscape_v5.py
~~~

4. Run each Methods 2.12 analysis with a distinct output directory. Reset PANTE_OUTPUT_DIR before each run. The master matrix defaults to $PANTE_INTEGRATION_DIR/01_matrix/PanTE_multiomics_master.RNA_testable.tsv.

~~~bash
export PANTE_OUTPUT_DIR="$PANTE_INTEGRATION_DIR/02_statistics/05_genomic_divergence_vs_epigenetic_deviation"
python 02_genomic_divergence_and_epigenetic_deviation/analyze_genomic_divergence_vs_epigenetic_deviation.py
export PANTE_OUTPUT_DIR="$PANTE_INTEGRATION_DIR/02_statistics/06_TEclass_abundance_controlled_sensitivity"
python 02_genomic_divergence_and_epigenetic_deviation/analyze_TEclass_abundance_controlled_sensitivity.py
export PANTE_OUTPUT_DIR="$PANTE_INTEGRATION_DIR/02_statistics/07_Figure6C_abundance_adjusted_reanalysis"
python 03_parental_divergence_and_transgressive_transcription/reanalyze_Figure6C_abundance_adjusted.py
export PANTE_OUTPUT_DIR="$PANTE_INTEGRATION_DIR/02_statistics/02_transgressive_vs_epigenetic_deviation"
python 04_transgressive_epigenetic_deviation/analyze_transgressive_epigenetic_deviation.py
export PANTE_OUTPUT_DIR="$PANTE_INTEGRATION_DIR/02_statistics/03_transgressive_directional_epigenetic_concordance"
python 04_transgressive_epigenetic_deviation/analyze_transgressive_directional_epigenetic_concordance.py
export PANTE_OUTPUT_DIR="$PANTE_INTEGRATION_DIR/02_statistics/04_background_centered_directional_permutation"
python 04_transgressive_epigenetic_deviation/analyze_background_centered_directional_permutation.py
~~~

## Software dependencies

Module 05 requires Python 3, NumPy, pandas and SciPy; the GLM script also requires statsmodels; Figure 5 plotting requires Matplotlib. Package versions were not pinned in the supplied source record. See Module 06 software_versions.tsv; do not infer exact versions from the current environment.

## Release and interpretation notes

- Large production matrices, raw FASTQ/BAM files and manuscript result tables are not copied into this module.
- The final Methods 2.12 logistic-regression analysis is the abundance-adjusted GLM script above. The earlier unadjusted Fisher/CMH exploratory script is not included in the primary workflow.
- Bootstrap confidence intervals are descriptive uncertainty intervals for pooled OLS slopes, not causal-model estimates.
- Interpret results within the documented data availability and filtering backgrounds.
