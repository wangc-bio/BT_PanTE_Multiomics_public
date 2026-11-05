# Pan-TE genomic occupancy quantification

This subworkflow completes the genomic-occupancy component of Methods 2.6. It uses the refined Pan-TE membership table and RepeatMasker `.out` files generated with the same unified v5.2 repeat library for cross-assembly comparison.

## 1. Resolve RepeatMasker overlaps and quantify Pan-TE occupancy

Script: `01_resolve_pante_abundance.py`

Run once per genome:

```bash
python3 01_resolve_pante_abundance.py \
    /path/to/PanTE_refined_cluster_members.tsv \
    /path/to/RepeatMasker.out \
    genome_01 \
    /path/to/output/genome_01
```

The script maps RepeatMasker consensus IDs to refined Pan-TE IDs and resolves overlapping annotations at base-pair resolution. All RepeatMasker hits participate in the overlap competition, including hits that cannot be assigned to a classified Pan-TE cluster. At each genomic segment, a non-`*` primary hit takes precedence over a `*` secondary hit; ties are resolved by higher Smith–Waterman score and then greater hit length. Each genomic base is assigned to one winning annotation only.

The script writes:

- `*.PanTE_resolved_abundance.tsv`: occupied base pairs and resolved segments per Pan-TE cluster;
- `*.resolution_qc.tsv`: total hits, primary/secondary counts, RepeatMasker union coverage, Pan-TE-resolved coverage and OTHER coverage.

The `Repeat_union_bp` QC value should be compared with the `bases masked` value in the corresponding RepeatMasker `.tbl` file.

## 2. Build the cross-assembly abundance matrix

Script: `02_build_resolved_matrix.py`

Copy `occupancy_manifest.tsv.example` to a local manifest and supply one row per assembly. The manifest contains the assembly identifier, total FASTA length in base pairs, and the abundance table generated in step 1. The script does not hard-code sample names, assembly sizes or the number of assemblies.

```bash
python3 02_build_resolved_matrix.py \
    --members /path/to/PanTE_refined_cluster_members.tsv \
    --manifest /path/to/occupancy_manifest.tsv \
    --out /path/to/PanTE_resolved_abundance_matrix.tsv \
    --presence-bp 10000
```

For each Pan-TE cluster and assembly, the matrix reports occupied base pairs, resolved segments and assembly-normalized occupancy:

`bp_per_Mb = occupied_bp / genome_size_bp × 1,000,000`

Clusters absent from a per-assembly abundance table are assigned zero occupied base pairs. A companion summary table reports total occupied base pairs, detected Pan-TE clusters and the number meeting the supplied occupied-bp threshold.

## 3. Compare parental and corresponding F1 haplotype occupancy

Script: `03_parental_f1_occupancy_contrast.py`

Copy `occupancy_roles.tsv.example` to a local role manifest and map the four roles (`parent_A`, `parent_B`, `f1_haplotype_A`, `f1_haplotype_B`) to the sample prefixes in the abundance matrix.

```bash
python3 03_parental_f1_occupancy_contrast.py \
    --matrix /path/to/PanTE_resolved_abundance_matrix.tsv \
    --role-manifest /path/to/occupancy_roles.tsv \
    --outdir /path/to/parent_f1_comparison \
    --presence-bp 10000 \
    --parent-fold 2 \
    --pseudocount 1
```

Parent-informative status is determined from raw occupied base pairs in either parent. Parent bias is determined from the ratio of assembly-normalized occupancy (bp/Mb), without a pseudocount. The log2 parental and F1 haplotypic contrasts use the specified pseudocount. Spearman correlation and sign concordance are evaluated among parent-informative clusters. Directional retention among parent-biased clusters is assessed by comparing the two F1 haplotype bp/Mb values; ties are reported separately.

See [the parent/F1 comparison details](PARENT_F1_COMPARISON_README.md) for the role manifest and output definitions.

## Scope and interpretation

The abundance values describe reference-genome occupancy of homologous Pan-TE clusters. They are not population-level insertion frequencies and do not establish insertion-site presence/absence in a population. The same unified repeat library must be used for all assemblies being compared so that differences in library composition do not become an additional source of ascertainment bias.

## Inputs not included

Genome FASTA files, RepeatMasker `.out`/`.tbl` files, the unified repeat library, production manifests and large abundance matrices are supplied externally and are not included in this repository.
