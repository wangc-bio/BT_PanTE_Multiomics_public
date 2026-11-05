# Module 05 source-script inventory

| Manuscript section | Public release path | Source basis | Function and retained settings |
|---|---|---|---|
| 2.11 | 01_parental_state_correspondence/prepare_Figure5_summary.py | Reworked from previous summary-preparation script | Calculates WGBS and ATAC matched-parent and parental-contrast Spearman summaries from Module 04B matrices; removes embedded fixed WGBS values. |
| 2.11 | 01_parental_state_correspondence/plot_Figure5_retention_landscape_v5.py | plot_Figure5_retention_landscape_v5.py | Figure 5 plots and pooled descriptive OLS/cluster bootstrap; 1,000 replicates; PanTE_ID cluster resampling; percentile 95% CI; WGBS seed 20260908, ATAC seed 20260909. |
| 2.12 | 00_build_multiomics_master/build_Result6_multiomics_master.py | build_Result6_multiomics_master.py | Integrates genomic occupancy, RNA inheritance, WGBS and ATAC; writes all-row, RNA-testable and complete multi-omics matrices plus QC summaries. |
| 2.12 | 02_genomic_divergence_and_epigenetic_deviation/analyze_genomic_divergence_vs_epigenetic_deviation.py | analyze_genomic_divergence_vs_epigenetic_deviation.py | Tissue-wise Dgenomic versus mCG/ATAC deviation Spearman tests; assay-available and parent-informative backgrounds; Parent_max_bp >= 10,000; BH correction. |
| 2.12 | 02_genomic_divergence_and_epigenetic_deviation/analyze_TEclass_abundance_controlled_sensitivity.py | analyze_TEclass_abundance_controlled_sensitivity.py | TE-class partial Spearman controlling total parental occupancy; DNA/LTR/LINE/SINE/RC; minimum N=30; BH correction within tissue and assay. |
| 2.12 | 03_parental_divergence_and_transgressive_transcription/reanalyze_Figure6C_abundance_adjusted.py | reanalyze_Figure6C_abundance_adjusted.py | Abundance-adjusted strict-transgressive logistic regression; binomial GLM; HC0 tissue SE; pooled PanTE_ID-clustered SE with tissue covariate; abundance-decile sensitivity; BH correction. |
| 2.12 | 04_transgressive_epigenetic_deviation/analyze_transgressive_epigenetic_deviation.py | analyze_transgressive_epigenetic_deviation.py | Mann–Whitney U; Cliff's delta; BH correction. |
| 2.12 | 04_transgressive_epigenetic_deviation/analyze_transgressive_directional_epigenetic_concordance.py | analyze_transgressive_directional_epigenetic_concordance.py | Signed shifts; non-transgressive median background; Wilcoxon signed-rank; exact binomial; both-haplotype concordance. |
| 2.12 | 04_transgressive_epigenetic_deviation/analyze_background_centered_directional_permutation.py | analyze_background_centered_directional_permutation.py | 10,000 within-tissue/assay Up/Down label permutations; seed 20260909; preserves observed Up/Down counts and epigenetic shifts. |

## Source provenance and exclusions

- WGBS and ATAC per-assay retention scripts remain in Module 04B; Module 05 consumes their outputs and does not duplicate assay-level processing.
- Figure 5 summary preparation recalculates correlations from input matrices; its outputs were cross-checked against previously recorded summary tables.
- The earlier unadjusted Fisher/CMH script analyze_parent_divergence_vs_transgressive.py is not included in the primary release because final Methods 2.12 uses the abundance-adjusted binomial GLM.
- Machine-specific paths were replaced with environment variables. The directional-concordance script includes a SciPy binom_test compatibility wrapper returning the scalar P value expected by the original calculation.
