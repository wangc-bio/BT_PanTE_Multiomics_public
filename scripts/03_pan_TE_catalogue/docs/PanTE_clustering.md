# Pan-TE clustering (Methods 2.6)

This document covers the consensus-sequence clustering component of Methods 2.6. RepeatMasker overlap resolution, genomic occupancy calculation and parental/F1 occupancy comparisons are documented in the [Module 03 genomic-occupancy subworkflow](../genomic_occupancy/README.md).

## Input and class partitioning

The input was the Unified v5.2 curated repeat library constructed by combining exactly four genome-specific curated repeat libraries: BSB, TC, BT-M and BT-P. BTB-A and BTB-B were not included. The library-building workflow normalizes sequences to uppercase and removes exact sequence duplicates using MD5 hashing only; no CD-HIT or other near-identity clustering is applied. Classification metadata are reconciled for exact duplicate sequences: when exactly one known TE fine-level classification conflicts only with Unknown, the known classification is retained; other unresolved broad-class conflicts are assigned Unknown.

High-confidence classified consensus sequences from this four-source unified library were grouped into five broad TE classes (DNA, LTR, LINE, SINE and RC) before within-class clustering. The production manifest, unified FASTA and class-specific FASTA files are not included in the public repository. Consequently, the number of input records, classified consensus sequences and final Pan-TE clusters should be taken from the validated four-library production outputs rather than copied from historical summaries.

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

The final number of classified consensus sequences and refined Pan-TE clusters is specific to the four-library Unified v5.2 input and must be reported from the corresponding production outputs. Historical numerical summaries generated from a six-library input set are not applicable to this manuscript analysis and are not reproduced here. The public repository provides the parameterized clustering scripts but does not include the production input FASTAs or membership/QC tables needed to independently verify those counts.

Pan-TE clusters are operational sequence-homology units and are not assumed to correspond directly to biologically defined TE families.

## Generated files

- `PanTE_cluster_members.tsv`: initial connected-component membership.
- `PanTE_cluster_summary.tsv`: initial clustering summary by broad class.
- `PanTE_refined_cluster_members.tsv`: final refined cluster membership.
- `PanTE_refined_cluster_summary.tsv`: final cluster counts and size summary.
- `PanTE_refined_component_audit.tsv`: record of initial components split during refinement.
- `PanTE_refined_cluster_qc.tsv`: representative-to-member identity and reciprocal-coverage QC.
- `PanTE_large_cluster_graph_qc.tsv`: optional graph QC for initial components containing at least 10 members.
