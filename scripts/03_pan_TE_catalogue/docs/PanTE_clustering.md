# Pan-TE clustering (Methods 2.6)

This document covers the consensus-sequence clustering component of Methods 2.6. RepeatMasker overlap resolution, genomic occupancy calculation and parental/F1 occupancy comparisons are documented in the [Module 03 genomic-occupancy subworkflow](../genomic_occupancy/README.md).

## Input and class partitioning

The input was the unified curated repeat library generated for cross-assembly comparison. It contained 25,130 exact-sequence-deduplicated entries. The library was constructed using `unified_library/01_build_unified_library.py`, which applies 100% exact-sequence MD5 deduplication and reconciles classification metadata for duplicate sequences; no near-identity clustering was applied. When exactly one known TE fine-level classification conflicted only with Unknown, the known classification was retained; other unresolved broad-class conflicts were assigned Unknown. High-confidence classified consensus sequences were grouped into five broad TE classes (DNA, LTR, LINE, SINE and RC). The five class-specific FASTA files contained 19,132 sequences in total. The input manifest and production FASTA files are not included in the repository.

## Within-class all-versus-all BLASTn

BLASTn v2.15.0+ searches were performed separately within each broad class. The initial search used `-task blastn`, `-perc_identity 70`, `-evalue 1e-10`, `-dust no`, `-soft_masking false` and `-max_hsps 50`.

## High-confidence relationships

Multiple HSPs for each query-subject pair were integrated by merging their covered intervals separately on the query and subject sequences. A pair was retained as a high-confidence relationship when all of the following criteria were met:

- alignment-length-weighted sequence identity >=80%;
- merged query coverage >=80%;
- merged subject coverage >=80%;
- at least 100 bp of merged HSP coverage on both the query and subject sequences;
- best E-value <=1e-20.

The query and subject coverage values in this clustering step were calculated from merged HSP intervals; they were not taken from the BLASTn `qcovs` field used in the separate homology-based classification workflow.

## Initial clustering and refinement

Initial clusters were defined as connected components in the high-confidence sequence-similarity network. Consensus sequences without a qualifying edge were retained as singleton clusters.

Each initial component was then partitioned using a representative-centered greedy strategy. At each iteration, the remaining sequence with the largest number of direct high-confidence neighbours was selected as the representative. Ties were resolved by summed bit score and then sequence identifier. The representative and all directly connected remaining sequences formed one refined cluster; the process was repeated until all sequences in the component had been assigned. Thus, every non-representative member had a direct high-confidence relationship with its representative, but members of a refined cluster were not required to be pairwise connected.

## Final cluster summary

| Broad class | Consensus sequences | Initial components | Refined clusters | Singleton clusters | Multi-member clusters | Maximum cluster size |
|---|---:|---:|---:|---:|---:|---:|
| DNA | 12,350 | 8,070 | 8,164 | 6,169 | 1,995 | 19 |
| LTR | 5,701 | 2,822 | 2,927 | 1,830 | 1,097 | 56 |
| LINE | 802 | 644 | 646 | 550 | 96 | 6 |
| SINE | 61 | 34 | 34 | 21 | 13 | 6 |
| RC | 218 | 175 | 179 | 159 | 20 | 8 |
| **Total** | **19,132** | **11,745** | **11,950** | **8,729** | **3,221** | **56** |

Pan-TE clusters are operational sequence-homology units and are not assumed to correspond directly to biologically defined TE families.

**Four-sample provenance check:** the `19,132` consensus / `11,950` refined-cluster summary currently stored in this repository is traceable to a historical six-sample workflow. The current manuscript scope is BSB, TC, BT-M and BT-P only, with BTB-A and BTB-B excluded. Four-sample-specific class FASTA, BLAST, cluster-membership and QC records have not yet been located. Until those inputs are verified, the stored counts must be treated as historical and not as an independently confirmed four-sample result.

## Generated files

- `PanTE_cluster_members.tsv`: initial connected-component membership.
- `PanTE_cluster_summary.tsv`: initial clustering summary by broad class.
- `PanTE_refined_cluster_members.tsv`: final refined cluster membership.
- `PanTE_refined_cluster_summary.tsv`: final cluster counts and size summary.
- `PanTE_refined_component_audit.tsv`: record of initial components split during refinement.
- `PanTE_refined_cluster_qc.tsv`: representative-to-member identity and reciprocal-coverage QC.
- `PanTE_large_cluster_graph_qc.tsv`: optional graph QC for initial components containing at least 10 members.
