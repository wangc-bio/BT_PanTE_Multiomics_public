# Module 03 — Pan-TE catalogue construction

## Scope

This module documents TE identification, consensus-library curation, genome-wide annotation, sequence-based Pan-TE clustering and genomic-occupancy quantification across the four focal genome assemblies. Each assembly is processed independently using the same parameterized workflow. No assembly-specific script copies are maintained.

## Workflow

| Manuscript section | Analysis | Documentation |
|---|---|---|
| 2.5 | TE identification, classification, consensus curation and genome-wide annotation | [TE identification and annotation](docs/TE_identification_and_annotation.md) |
| 2.6 | Unified-library construction, within-class sequence-similarity search, representative-centered refinement, RepeatMasker overlap resolution and parental/F1 occupancy comparison | [Unified library](unified_library/README.md); [Pan-TE clustering](docs/PanTE_clustering.md); [Genomic occupancy](genomic_occupancy/README.md) |

Genome-specific curated repeat libraries are used for the primary RepeatMasker annotations. A unified curated repeat library is constructed for cross-assembly quantitative comparison and sensitivity analysis; it also provides the consensus sequences for Pan-TE clustering. These are distinct library uses and should not be conflated. Unified-library construction uses 100% exact-sequence MD5 deduplication; near-identity clustering (for example, CD-HIT) was not used. For exact duplicates with one known TE classification and an otherwise Unknown classification, the known classification is retained; other unresolved broad-class conflicts are assigned Unknown.

## Pan-TE catalogue

High-confidence classified consensus sequences were assigned to five broad classes (DNA, LTR, LINE, SINE and RC) and clustered within class. The currently archived summary contains 19,132 classified consensus sequences and 11,950 refined Pan-TE clusters. The four-sample Pan-TE catalogue contains 19,132 classified consensus sequences and 11,950 refined Pan-TE clusters. The unified library contains 25,540 input consensus records and 25,130 exact-sequence-unique records before broad-class filtering. BTB-A and BTB-B are not included in this manuscript analysis.

| Broad class | Consensus sequences | Initial components | Refined clusters | Singleton clusters | Multi-member clusters | Maximum cluster size |
|---|---:|---:|---:|---:|---:|---:|
| DNA | 12,350 | 8,070 | 8,164 | 6,169 | 1,995 | 19 |
| LTR | 5,701 | 2,822 | 2,927 | 1,830 | 1,097 | 56 |
| LINE | 802 | 644 | 646 | 550 | 96 | 6 |
| SINE | 61 | 34 | 34 | 21 | 13 | 6 |
| RC | 218 | 175 | 179 | 159 | 20 | 8 |
| **Total** | **19,132** | **11,745** | **11,950** | **8,729** | **3,221** | **56** |

Pan-TE clusters are operational sequence-homology units. They are not assumed to correspond directly to biologically defined TE families and are not genomic insertion loci.

## Module boundary

Module 03 includes the refined consensus-level Pan-TE catalogue, base-level RepeatMasker overlap resolution, assembly-normalized genomic occupancy and parental/F1 occupancy comparisons. Module 04A documents the separate Kimura-divergence landscape analysis. Reference-genome occupancy values should not be interpreted as population-level insertion frequencies.

## Scripts and provenance

- [TE annotation scripts](te_annotation/README.md) provide one parameterized script per software/workflow stage and a generic external configuration template.
- [Unified-library construction](unified_library/README.md) provides the manifest-driven exact-sequence deduplication and classification-reconciliation script.
- [Evidence integration](evidence_integration/README.md) documents reciprocal structural-evidence filtering, family-level evidence summaries, evidence-master construction, strict and operational classification, and the v5.1/v5.2 Unknown-rescue stages.
- [Pan-TE clustering scripts](pante_clustering/README.md) provide the five scripts used for consensus splitting, BLASTn searches, initial clustering, optional graph QC and representative-centered refinement.
- [Genomic occupancy scripts](genomic_occupancy/README.md) provide the generic sweep-line overlap resolver, matrix builder and parent/F1 occupancy comparison.
- [Source-script inventory](SOURCE_SCRIPT_INVENTORY.md) records the source basis, retained parameters and provenance limitations for each workflow component.
- [Classification synthetic test](tests/test_classification_release.py) exercises strict arbitration, operational arbitration and the two-stage rescue scripts with synthetic inputs. The test passes in GitHub Actions (run 41864679408); it validates the release code paths, not production classifications.
- The final cluster summary is provided in [results/PanTE_refined_cluster_summary.tsv](results/PanTE_refined_cluster_summary.tsv). The minimal public-release policy for result files is described in [results/README.md](results/README.md).

The TE annotation and clustering scripts use generic input variables and do not embed assembly-specific FASTA paths. The clustering release copies retain the original thresholds and clustering logic.

The initial RepBase consensus screen and subsequent RepBase/Dfam rescue searches produce hit-level evidence. The `evidence_integration/` directory includes the Dfam, RepBase and TEsorter family-summary generators, reciprocal structural-evidence processing and the evidence-master builder. The evidence-integration directory now includes the portable, manifest-driven strict arbitration (`04_final_arbitration.py`), operational arbitration (`05_operational_arbitration.py`), v5.1 RepBase rescue (`06_relaxed_repbase_rescue_v5_1.py`) and v5.2 RepBase/Dfam rescue (`07_relaxed_homology_rescue_v5_2.py`) scripts. A synthetic end-to-end test is provided in `tests/test_classification_release.py`; it tests release logic with synthetic inputs and does not reproduce production classifications. For the manuscript workflow, the production sample set is BSB, TC, BTM and BTP; the manifest-based release scripts are configured with these four inputs.

**Four-sample publication scope:** all Pan-TE results reported in the manuscript are based on BSB, TC, BT-M and BT-P. BTB-A and BTB-B are excluded. The reported 19,132 classified consensus sequences and 11,950 refined clusters are the four-sample catalogue used for the manuscript analysis.
