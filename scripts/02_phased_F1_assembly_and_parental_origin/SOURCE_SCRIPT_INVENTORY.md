# Module 02 source-script inventory

Updated: 2026-11-04. This inventory records the original source files reviewed and their corresponding public scripts.

| Component | Original source material reviewed | Public scripts |
|---|---|---|
| Whole-genome alignment | FB_M_vs_P.sh; user-confirmed identical MUMmer commands/parameters across BT-M/BT-P versus BSB/TC reference comparisons | `01_whole_genome_alignment/01_MUMmer_alignment.sh` |
| Synteny analysis | 01-03_gff2bed_deduplication_extract.sh; 04_compara.sh; 05_build_simpe.sh; 06_plot.sh | `02_synteny_analysis/01_prepare_bed_cds.sh`; `02_run_mcscan.sh`; `03_filter_anchors.sh`; `04_plot_karyotype.sh` |
| Phylogenetic validation | 01_orthofinder_chr.sh; 02_muscle.sh; 03_Gblock.sh; 04-05_seqkit.sh; 06_prottest.sh; 07_raxml_chr.sh | `03_phylogenetic_validation/01_orthofinder_chr.sh`; `02_muscle_alignment.sh`; `03_gblocks_filter.sh`; `04_seqkit_format.sh`; `05_prottest.sh`; `06_raxml.sh` |
| Parent-specific k-mer | 01_jellyfish.sh; 02_k31_1Mb_map2BT.sh; 03_validate_5Mb.sh | `04_parent_specific_kmer/01_jellyfish_count.sh`; `02_parent_specific_kmer_mapping.sh`; `03_validate_5Mb.sh` |
| Parental depth validation | 01_BSB_mapping.sh; 02_step2_filter_bam_BSB.sh; 02_step2_filter_bam_TC.sh; 03_depth_by_window_1Mb.sh; 03_depth_by_window_5Mb.sh; 04_step4_merge_matrix_1Mb.sh; 04_step4_merge_matrix_5Mb_1.sh; 05_depth_validate_1Mb.sh; 05_depth_validate_5Mb.sh; 05_depth_validate_nd.py; 06_bam.sh | `05_parental_depth_validation/01_map_parental_reads.sh`; `01_window_definitions.sh`; `02_filter_bam.sh`; `03_depth_by_window.sh`; `04_merge_ND_matrix.sh`; `05_depth_validate_nd.py`; `05_validate_depth_vs_kmer.sh` |
| Homologous gene pairs | BT_hom_gene_pairs.sh | `06_homologous_gene_identification/01_MMseqs2_RBH.sh` |

## Confirmed method details

- MUMmer: nucmer `-c 90 -l 40 -t 30`; delta-filter `-1 -i 90 -o 90 -l 5000`.
- JCVI/MCscan: six pairwise comparisons, C-score 0.80, anchor screen `--minspan=30 --simple`.
- Phylogeny: chr01–chr24; OrthoFinder DIAMOND search; MUSCLE v3; Gblocks `-b4=5 -b5=h -t=p`; SeqKit sorting/formatting and paste/sed merge followed by manual taxon-ID calibration and matrix inspection; ProtTest and RAxML settings are documented in the phylogenetic README.
- Parent-specific k-mer: Jellyfish k=31 and k=21 counts; k=31 used for ancestry; copy-number range 1–5; Bowtie exact unique mapping; 1-Mb and 5-Mb calls; minimum 50 hits and parental fraction 0.70. Both window-level and chromosome-level source classifications are included.
- Parental depth: MAPQ >=20, proper-pair filtering, duplicates retained by default; ND is window mean depth divided by length-weighted genome-wide mean depth; Welch t-test, BH-FDR, FDR 0.05 and |log2FC| 0.3 in the study wrapper. The validation Python source is included.
- Homologous pairs: MMseqs2 easy-rbh, sensitivity 7.5, E-value 1e-5, minimum identity 0.5, coverage 0.6, coverage mode 0.

## Release notes

- Scripts use configurable variables rather than embedding the original server-specific absolute paths.
- Results already present in the manuscript or supplementary materials are intentionally not duplicated in GitHub.
- The study-specific JCVI `seqids` and `layout` files are required only to recreate the exact karyotype visualization; they are not needed to describe the underlying synteny comparisons.
- Software versions are included only where verified from the supplied source materials.
- This is a source-level reproducibility review. End-to-end execution on an independent clean environment has not been performed as part of this repository audit.
