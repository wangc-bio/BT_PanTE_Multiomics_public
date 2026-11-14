# Module 03 — Pan-TE catalogue construction

## Scope

This module documents TE identification, consensus-library curation, genome-wide annotation, sequence-based Pan-TE clustering and genomic-occupancy quantification across the four focal genome assemblies. Each assembly is processed independently using the same parameterized workflow. No assembly-specific script copies are maintained.

## Workflow

| Manuscript section | Analysis | Documentation |
|---|---|---|
| 2.5 | TE identification, classification, consensus curation and genome-wide annotation | [TE identification and annotation](docs/TE_identification_and_annotation.md) |
| 2.5 | Genome-specific TE annotation and construction of the four-source Unified v5.2 library | [TE identification and annotation](docs/TE_identification_and_annotation.md); [Unified library](unified_library/README.md) |
| 2.6 | Within-class sequence-similarity search, representative-centered refinement, RepeatMasker overlap resolution and parental/F1 occupancy comparison | [Pan-TE clustering](docs/PanTE_clustering.md); [Genomic occupancy](genomic_occupancy/README.md) |

Genome-specific curated repeat libraries are used for the primary RepeatMasker annotations. For the manuscript analysis, Unified v5.2 was constructed by combining exactly four genome-specific curated repeat libraries: BSB, TC, BT-M and BT-P. BTB-A and BTB-B were not included. The unified library provides a common consensus reference for Pan-TE clustering and the unified-library RepeatMasker annotations used in cross-assembly genomic-occupancy, Kimura-divergence and Pan-TE transcript-quantification analyses. These uses are distinct from the primary genome-specific annotations. Unified-library construction uses uppercase sequence normalization and 100% exact-sequence MD5 deduplication; near-identity clustering (for example, CD-HIT) was not used. For exact duplicates with one known TE classification and an otherwise Unknown classification, the known classification is retained; other unresolved broad-class conflicts are assigned Unknown.

## Unified v5.2 library and Pan-TE catalogue

For the manuscript analysis, Unified v5.2 was constructed from exactly four genome-specific curated repeat libraries: BSB, TC, BT-M and BT-P. BTB-A and BTB-B were not included as input libraries. The library-building script performs uppercase sequence normalization and 100% exact-sequence deduplication using MD5; it does not perform near-identity clustering. Classification reconciliation is described in the [unified-library workflow](unified_library/README.md).

The resulting unified library provides the common consensus reference for Pan-TE clustering and for the unified-library RepeatMasker annotations used in cross-assembly genomic-occupancy, Kimura-divergence and Pan-TE transcript-quantification analyses. The primary genome-specific annotations described in Methods 2.5 remain separate.

The production input manifest and FASTA files are not distributed in this repository. Therefore, input-record totals, unique-consensus counts, broad-class counts and refined Pan-TE cluster counts are not asserted here. These values should be reported only after checking the final four-library production FASTA, sequence map, class-specific FASTAs and cluster outputs. Historical counts from a six-library workflow must not be reused for the four-library manuscript analysis.

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
- Production Pan-TE cluster-membership and summary tables are data-dependent outputs and are not included in this public repository.

The TE annotation and clustering scripts use generic input variables and do not embed assembly-specific FASTA paths. The clustering release copies retain the original thresholds and clustering logic.

The initial RepBase consensus screen and subsequent RepBase/Dfam rescue searches produce hit-level evidence. The `evidence_integration/` directory includes the Dfam, RepBase and TEsorter family-summary generators, reciprocal structural-evidence processing and the evidence-master builder. The evidence-integration directory now includes the portable, manifest-driven strict arbitration (`04_final_arbitration.py`), operational arbitration (`05_operational_arbitration.py`), v5.1 RepBase rescue (`06_relaxed_repbase_rescue_v5_1.py`) and v5.2 RepBase/Dfam rescue (`07_relaxed_homology_rescue_v5_2.py`) scripts. A synthetic end-to-end test is provided in `tests/test_classification_release.py`; it tests release logic with synthetic inputs and does not reproduce production classifications. For the manuscript workflow, the production sample set is BSB, TC, BTM and BTP; the manifest-based release scripts are configured with these four inputs.

**Four-assembly publication scope:** the manuscript-level Pan-TE analyses use BSB, TC, BT-M and BT-P only. BTB-A and BTB-B are excluded from the Unified v5.2 input, Pan-TE clustering and four-assembly comparisons. Numerical catalogue and cluster counts must be taken from the validated four-library production outputs.
