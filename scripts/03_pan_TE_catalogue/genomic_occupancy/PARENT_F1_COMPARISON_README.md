# Pan-TE parental and F1 occupancy comparisons

This subworkflow implements the parent-informative, parent-biased and parental-versus-F1 contrast analyses described in Methods 2.6. It is parameterized by role rather than by fixed sample names.

## Input role manifest

Copy `occupancy_roles.tsv.example` to a local TSV file and map the four analysis roles to the corresponding sample prefixes in the Pan-TE abundance matrix:

- `parent_A`: first parental genome;
- `parent_B`: second parental genome;
- `f1_haplotype_A`: F1 haplotype corresponding to parent A;
- `f1_haplotype_B`: F1 haplotype corresponding to parent B.

The four sample IDs must be distinct and must match the matrix column prefixes.

## Run
```bash
python3 03_parental_f1_occupancy_contrast.py \
    --matrix /path/to/PanTE_resolved_abundance_matrix.tsv \
    --role-manifest /path/to/occupancy_roles.tsv \
    --outdir /path/to/parent_f1_comparison \
    --presence-bp 10000 \
    --parent-fold 2 \
    --pseudocount 1
```

Dependencies: Python 3, NumPy, pandas and SciPy. Exact package versions were not recorded in the source execution record.

## Definitions

1. A Pan-TE cluster is parent-informative when its raw occupied base-pair count reaches the presence threshold in either parent.
2. Parental bias is calculated from assembly-normalized occupancy (bp/Mb), without a pseudocount. A ratio >=2 identifies parent-A bias; a ratio <=0.5 identifies parent-B bias; intermediate ratios are classified as similar.
3. The parental and F1 contrasts are calculated as log2 ratios of bp/Mb values after adding the specified pseudocount (default 1 bp/Mb).
4. Spearman rank correlation and sign concordance are calculated among parent-informative clusters. Sign concordance excludes pairs where either contrast equals zero.
5. Directional retention among parent-biased clusters is assessed by comparing the two F1 haplotype bp/Mb values. Ties are counted separately and excluded from the directional-concordance denominator.

## Outputs

- `parent_f1_contrast_per_PanTE.tsv`: parent-informative status, parental bias, log2 contrasts and F1 directional-retention label for every Pan-TE cluster.
- `parent_f1_contrast_statistics.tsv`: Spearman association and sign-concordance summary for parent-informative clusters.
- `parent_biased_directional_retention.tsv`: concordant, discordant and tied F1 directions for each parental-bias class and combined parent-biased set.
- `run_parameters.tsv`: analysis-role mapping and thresholds.

No population-level insertion frequencies or insertion-site presence/absence calls are inferred from these reference-genome occupancy comparisons.
