# TE evidence integration

## Scope

This subworkflow converts genome-specific RepeatModeler consensus evidence into family-level evidence tables and an integrated evidence master. It combines three evidence dimensions:

- **Homology:** RepBase and Dfam.
- **Domain:** TEsorter classification and conserved-domain matches.
- **Structure:** reciprocal LTR_retriever, MITE-Hunter and HelitronScanner matches.

All scripts use generic library identifiers and runtime input/output paths. Run the same scripts separately for each of the four focal genome libraries; no genome names, production paths or production data are embedded.

## Execution order

### 1. Filter reciprocal Class II evidence

After mapping the original RepeatModeler consensus sequences against the MITE-Hunter and HelitronScanner candidate databases, run:

```bash
export CLASSII_MAPPING_DIR=/path/to/classII_mapping_output
bash scripts/03_pan_TE_catalogue/evidence_integration/01_filter_classII_reciprocal.sh
```

The filter retains the source workflow's per-HSP `qcovhsp` values and calculates subject coverage as alignment length divided by subject length. Supportive thresholds are identity >=70%, alignment >=100 bp, E-value <=1e-10, query HSP coverage >=50% and subject coverage >=50%. High-confidence thresholds are identity >=80%, alignment >=100 bp, E-value <=1e-20, query HSP coverage >=80% and subject coverage >=80%. The best hit per query is selected by bit score, identity, qcovhsp and subject coverage.

### 2. Re-evaluate reciprocal LTR evidence

Prepare a TSV manifest with columns `sample_id` and `ltr_raw_blast_tsv`, then run:

```bash
python scripts/03_pan_TE_catalogue/evidence_integration/02_rethreshold_LTR_reciprocal.py \
  --manifest /path/to/ltr_manifest.tsv \
  --output-dir /path/to/ltr_reciprocal_output
```

HSP intervals are merged independently on query and subject sequences. Coverage is calculated from the merged interval unions; identity is alignment-length weighted. High-confidence matches require identity >=80%, query and subject coverage >=80%, at least 100 bp of reciprocal covered sequence and E-value <=1e-20. Supportive matches use identity >=70%, query and subject coverage >=50%, at least 100 bp and E-value <=1e-10. Evidence from direct LTR regions is distinguished from internal-region evidence.

### 3. Generate Dfam, RepBase and TEsorter family summaries

Each summary program is run once per consensus library.

**Dfam:** combine the Xenocyprididae and Cypriniformes BLAST tables:

```bash
python scripts/03_pan_TE_catalogue/evidence_integration/summarize_Dfam_v5.py \
  LIBRARY_ID OUTPUT_DIR XENOCYPRIDIDAE_RAW.tsv CYPRINIFORMES_RAW.tsv
```

**RepBase:** summarize the BLAST output from the original RepeatModeler consensus library:

```bash
python scripts/03_pan_TE_catalogue/evidence_integration/normalize_RepBase_v5.py \
  --sample-id LIBRARY_ID \
  --raw-tsv REPBASE_RAW.tsv \
  --output-dir OUTPUT_DIR
```

**TEsorter:** combine the `.cls.tsv` classification and `.dom.tsv` domain outputs:

```bash
python scripts/03_pan_TE_catalogue/evidence_integration/normalize_TEsorter_v5.py \
  --sample-id LIBRARY_ID \
  --cls-tsv TEsorter.cls.tsv \
  --dom-tsv TEsorter.dom.tsv \
  --output-dir OUTPUT_DIR
```

The Dfam and RepBase summaries classify query–subject pairs after merging overlapping HSP intervals. Query coverage is recalculated as the union of query-aligned bases divided by query length; subject coverage is calculated analogously. Although the BLAST search output includes `qcovs`, the normalization scripts calculate query coverage from the merged HSP coordinates rather than reading the `qcovs` output column.

For both Dfam and RepBase, STRONG evidence requires identity >=80%, query coverage >=80%, at least 80 query-covered bp and E-value <=1e-10. SUPPORTIVE evidence requires identity >=70%, query coverage >=50%, at least 80 query-covered bp and E-value <=1e-5. Family-level assignments are retained only when the qualifying broad-class evidence is not conflicting.

TEsorter domain evidence is STRONG at E-value <=1e-5 and domain coverage >=30%, and SUPPORTIVE at E-value <=1e-3 and coverage >=20%. The `.cls.tsv` classification is reconciled with domain-based class assignments; incompatible strong evidence is recorded as a conflict.

The programs produce pair-level (where applicable), family-level and per-library QC tables. These outputs are intermediate evidence tables; they are not themselves the final curated repeat library.

### 4. Build the integrated evidence master

Prepare `evidence_master_manifest.tsv` using the example template. Each row represents one genome-specific consensus library and points to its RepeatModeler FASTA, prior TEtrimmer/MCHelper evidence, Dfam/RepBase/TEsorter family summaries, reciprocal LTR family-best table, and four reciprocal Class II best-hit tables.

```bash
python scripts/03_pan_TE_catalogue/evidence_integration/03_build_evidence_master.py \
  --manifest /path/to/evidence_master_manifest.tsv \
  --output-dir /path/to/evidence_master \
  --review-dir /path/to/conflict_review
```

The master builder preserves the supplied evidence-combination rules. It combines RepBase and Dfam within the homology dimension, uses TEsorter as the domain dimension, and combines LTR/MITE/Helitron calls as the structural dimension. High-confidence Helitron evidence is conservatively capped at the supportive tier in the structural dimension. The output records cross-dimension status and flags records requiring review; final strict/operational classification is a separate downstream step.

## Required manifest columns

`sample_id`, `repeatmodeler_fasta`, `prior_evidence`, `dfam_family_summary`, `repbase_family_summary`, `tesorter_family_summary`, `ltr_reciprocal_family_best`, `mite_hc_best`, `mite_supportive_best`, `helitron_hc_best`, and `helitron_supportive_best`.

All paths may be absolute or relative to the manifest file. Required tables must exist; empty Class II best-hit files are allowed when no qualifying hit was observed.

## Remaining source-recovery item

A TC v5.1 candidate-classification table is available as a result, but the supplied `multi_sample_relaxed_repbase_v5_1.py` does not include TC in its configured sample list. The exact program or run-specific modification that generated the TC v5.1 table has not yet been identified. This downstream v5.1 step is therefore not claimed as fully reproducible here.

No production FASTAs, BLAST outputs, databases or large intermediate evidence tables are included in this directory.


## 5. Strict and operational classification; two-stage Unknown rescue

The evidence master is an intermediate evidence table, not the final repeat classification. Four portable, manifest-driven downstream scripts are now included in this directory:

| Stage | Release script | Historical source basis | Main output |
|---|---|---|---|
| Strict arbitration | `04_final_arbitration.py` | `final_arbitration_v5.py` | `*.v5_final_classification.tsv` plus a review table |
| Operational arbitration | `05_operational_arbitration.py` | `operational_arbitration_v5.py` | `*.v5_operational_classification.tsv` plus manual-review and summary tables |
| v5.1 RepBase rescue | `06_relaxed_repbase_rescue_v5_1.py` | `multi_sample_relaxed_repbase_v5_1.py` | `*.v5_1_candidate_classification.tsv`, query FASTA and v5.1 libraries |
| v5.2 RepBase/Dfam rescue | `07_relaxed_homology_rescue_v5_2.py` | `v52_relaxed_homology_rescue.py` | `*.v5_2_classification.tsv`, homology-call/QC tables and v5.2 libraries |

The release scripts use a runtime TSV manifest rather than hard-coded sample names or server paths. The production manuscript analysis uses four samples: BSB, TC, BTM (BT-M) and BTP (BT-P). BTBA and BTBB are excluded. The manifest-based release scripts should be configured with exactly these four production inputs.

### 5.1 Strict arbitration

Run `04_final_arbitration.py` after the evidence master has been generated. The script evaluates three independent evidence dimensions: homology (RepBase and Dfam combined), conserved domain (TEsorter), and structure (LTR/MITE/Helitron evidence). Prior/composite RepeatModeler, TEtrimmer and MCHelper classifications are not counted as independent dimensions for assigning a rescued broad class; compatible prior evidence may be used to retain a fine-level label after a broad class is supported externally. Conflicting or insufficient external evidence is retained as Unknown or sent to review according to the source rules.

Prepare a manifest with columns `sample_id` and `evidence_master`, then run:

```bash
python scripts/03_pan_TE_catalogue/evidence_integration/04_final_arbitration.py \
  --manifest /path/to/strict_manifest.tsv \
  --output-dir /path/to/strict_classification \
  --review-dir /path/to/strict_review
```

### 5.2 Operational arbitration

Run `05_operational_arbitration.py` using a manifest with `sample_id`, `evidence_master` and `strict_classification`. This preserves the operational source logic: prior/composite classifications are retained unless the external evidence rules support rescue or override; conflicts that cannot be resolved automatically are flagged for manual review. Strict evidence-based classification is retained as a separate confirmation layer rather than substituted for the operational classification.

```bash
python scripts/03_pan_TE_catalogue/evidence_integration/05_operational_arbitration.py \
  --manifest /path/to/operational_manifest.tsv \
  --output-dir /path/to/operational_classification \
  --review-dir /path/to/operational_review
```

### 5.3 v5.1 relaxed RepBase rescue

The v5.1 stage extracts families still classified as Unknown in the operational result. The historical BLAST search used `-evalue 1e-3`, `-perc_identity 60`, `-word_size 7`, `-max_target_seqs 100` and 6 threads. The family-level parser applies identity >=65%, merged query coverage >=30%, merged query-covered length >=100 bp and E-value <=1e-3. A broad-class rescue additionally requires at least two distinct qualifying RepBase subjects for the winning class and a winning-class bit-score fraction >=0.70.

Prepare the Unknown query FASTA:

```bash
python scripts/03_pan_TE_catalogue/evidence_integration/06_relaxed_repbase_rescue_v5_1.py \
  --manifest /path/to/classification_manifest.tsv --sample-id SAMPLE \
  --mode prepare --output-dir /path/to/v5_1
```

After the historical-equivalent relaxed BLAST search has been run and its 14-column output path is present in the manifest, summarize:

```bash
python scripts/03_pan_TE_catalogue/evidence_integration/06_relaxed_repbase_rescue_v5_1.py \
  --manifest /path/to/classification_manifest.tsv --sample-id SAMPLE \
  --mode summarize --output-dir /path/to/v5_1
```

The v5.1 parser retains operationally classified families unchanged. Its thresholds are stage-specific and are not interchangeable with the v5.2 thresholds below.

**Publication sample scope:** the production classification and downstream Pan-TE analysis use BSB, TC, BTM and BTP. The public release scripts use a runtime manifest to make this four-sample input set explicit; auxiliary helper-script arrays are not used to redefine the manuscript dataset.

### 5.4 v5.2 relaxed RepBase/Dfam rescue

The v5.2 stage starts from families remaining Unknown after v5.1, reuses the existing v5.1 RepBase BLAST output and incorporates a relaxed search against the curated Cyprinidae Dfam/FamDB TE-only reference. The historical Dfam resource contained 43,903 families before filtering, of which 34,402 classified TE references were retained across DNA, LTR, LINE, SINE and RC. The exact Dfam/FamDB release identifier has not been verified and is not assigned here.

The historical Dfam BLAST search used `-evalue 1e-3`, `-perc_identity 60`, `-word_size 7`, `-max_target_seqs 100` and 36 threads. For each database independently, the parser applies identity >=65%, merged query coverage >=20%, merged query-covered length >=80 bp and E-value <=1e-3; the top 10 qualifying subjects by bit score are considered. Rescue requires at least two distinct winning-class subjects and a winning-class bit-score fraction >=0.70. If RepBase and Dfam both pass and agree, the family is rescued; if they disagree, it remains Unknown. If only one database passes, that database's class is used. Families already classified in v5.1 remain unchanged.

```bash
python scripts/03_pan_TE_catalogue/evidence_integration/07_relaxed_homology_rescue_v5_2.py \
  --manifest /path/to/classification_manifest.tsv --sample-id SAMPLE \
  --mode prepare --output-dir /path/to/v5_2

# Run the Dfam BLAST search on the generated remaining-Unknown FASTA.
# Then use the existing v5.1 RepBase BLAST and new Dfam BLAST outputs:
python scripts/03_pan_TE_catalogue/evidence_integration/07_relaxed_homology_rescue_v5_2.py \
  --manifest /path/to/classification_manifest.tsv --sample-id SAMPLE \
  --mode summarize --output-dir /path/to/v5_2
```

The general Methods description (identity >=65%, query-covered length >=80 bp and E-value <=1e-3) describes the relaxed homology-rescue criteria; the exact query-coverage and minimum-length filters differ between the v5.1 RepBase stage and v5.2 database-specific stage and are documented separately above. Do not collapse these stage-specific filters into a single identical threshold set.

### 5.5 Classification test and provenance limitations

A synthetic end-to-end test is provided at `tests/test_classification_release.py`. Run it from the repository root:

```bash
python3 scripts/03_pan_TE_catalogue/tests/test_classification_release.py
```

The test exercises strict arbitration, operational arbitration, v5.1 rescue and v5.2 rescue using synthetic records; it does not rerun production data or validate the biological accuracy of every historical classification.

The first test attempt failed because the synthetic Dfam BLAST rows used shortened subject IDs (`DF1`, `DF2`), while `parse_dfam_map()` keys the Dfam reference map by the complete FASTA defline token (for example, `DF1#LTR/Gypsy|Broad=LTR`). The test fixture was corrected to use the full subject IDs; no production classification code or threshold was changed. GitHub Actions subsequently passed compilation, the synthetic end-to-end test and CLI checks on commit `a22ba563e752cc56731bf64de8e21c19197c2948` (run [41864679408](https://github.com/wangc-bio/BT_PanTE_Multiomics/actions/runs/41864679408)).

The publication analysis scope is the four assemblies BSB, TC, BT-M and BT-P. The manifest-driven release scripts support this exact input set and avoid hard-coded sample arrays.

The four-sample Pan-TE catalogue contains 19,132 classified consensus sequences and 11,950 refined clusters. The unified library contains 25,540 input consensus records and 25,130 exact-sequence-unique records before broad-class filtering. These are the catalogue counts used for the manuscript's four-sample analysis (BSB, TC, BTM and BTP).

No production FASTAs, BLAST outputs, databases or large intermediate evidence tables are included in this directory.
