# Module 03 source-script inventory

## Scope and release convention

Each analysis stage is represented by one generic, parameterized script. Genome FASTA files, identifiers, private server paths and per-run output directories are supplied through an external configuration file. The same scripts are reused across the four focal assemblies; no per-assembly script copies are maintained.

## Methods 2.5 — TE identification and genome-wide annotation

| Workflow stage | Generic release script | Source basis and retained parameters |
|---|---|---|
| De novo repeat discovery | `te_annotation/01_repeatmodeler2.sh` | RepeatModeler2 v2.0.5; `BuildDatabase -name DB_NAME GENOME.fa`; `RepeatModeler -database DB_NAME -threads 36 -LTRStruct` |
| Consensus curation | `te_annotation/02_tetrimmer.sh` | TEtrimmer v1.7.2 recorded command; `--classify_unknown`, Pfam-A HMM directory, log file and INFO log level |
| Consensus curation | `te_annotation/03_mchelper.sh` | MCHelper v1.7.0 recorded command; `-r A -t 24 -l LIB -g GENOME.fa -o OUT --input_type fasta -b BUSCO_HMM -a F` |
| LTR structural evidence | `te_annotation/04_ltr_structural_annotation.sh` | LTR_FINDER_parallel v1.2 using LTR_FINDER v1.0.7; GenomeTools v1.6.2; LTR_retriever v3.0.1; recorded LTR_FINDER_parallel, suffixerator, ltrharvest and LTR_retriever options retained |
| MITE structural evidence | `te_annotation/05_mite_hunter.sh` | MITE-Hunter manager command; `-P 1 -n 36 -c THREADS -S 12345678`; version not recorded |
| Helitron structural evidence | `te_annotation/06_helitronscanner.sh` | HelitronScanner v1.1; forward and reverse-complement scans, pairing and pure candidate export |
| Domain evidence | `te_annotation/07_tesorter.sh` | TEsorter v1.5.1; SeqKit version not recorded; Unknown consensus extraction; `-db rexdb-metazoa -p 4` |
| Class II evidence mapping | `te_annotation/08_classII_evidence_mapping.sh` | Original RepeatModeler consensus sequences are the BLASTn queries; MITE-Hunter and HelitronScanner candidate FASTAs are separate subject databases. BLASTn v2.15.0+ (confirmed); supportive: identity >=70%, qcov_hsp >=50%, alignment >=100 bp, E <=1e-10; high confidence: identity >=80%, qcov_hsp >=80%, alignment >=100 bp, E <=1e-20. The reciprocal filter additionally requires single-HSP subject coverage >=50% or >=80%, respectively, calculated as alignment length / subject length |
| Initial RepBase consensus screening | `te_annotation/09_repbase_consensus_mapping.sh` | Original RepeatModeler consensus libraries against the pre-built RepBase release 31.07 fish/vertebrate clean database; BLASTn v2.15.0+; supportive: identity >=70%, qcovs >=50%, alignment >=100 bp, E <=1e-10; high confidence: identity >=80%, qcovs >=80%, alignment >=100 bp, E <=1e-20; best subject per query selected by bit score, identity and qcovs |
| Homology database preparation | `te_annotation/10_build_homology_databases.sh` | BLAST+ v2.15.0+; `makeblastdb -dbtype nucl` build from supplied Xenocyprididae and Cypriniformes Dfam reference FASTAs; reference extraction/curation is not performed by this script |
| Homology rescue search | `te_annotation/11_homology_rescue_search.sh` | RepBase release 31.07 and supplied Dfam nucleotide databases; BLASTn v2.15.0+; identity >=65%, qcovs >=20%, qcovs-derived query-covered length >=80 bp, E <=1e-3 |
| Genome-wide annotation | `te_annotation/12_repeatmasker.sh` | RepeatMasker v4.1.8 and RMBlast v2.11.0+; `-pa 9 -lib LIB -gff -html -xsmall -e rmblast -no_is`; one script accepts either the genome-specific primary library or the unified sensitivity library |

The external configuration template is `te_annotation/te_annotation.env.example`. It contains generic placeholders only. Actual consensus FASTAs, reference FASTAs, database prefixes, tool paths and output directories are supplied at runtime. The RepBase database is supplied as a pre-built nucleotide BLAST database; the two Dfam databases are built by the generic database-preparation script.

## Unified library construction

| Generic release script | Function | Release changes |
|---|---|---|
| `unified_library/01_build_unified_library.py` | Merge one or more curated repeat libraries, deduplicate exact sequences by MD5, reconcile duplicate classifications, and write FASTA/provenance/summary/conflict outputs | Input libraries are supplied through a generic TSV manifest; no sample names, server paths or production inputs are hard-coded |

The implementation preserves uppercase normalization, exact-sequence MD5 grouping, broad-class assignment and fine-class fallback labels. When exactly one known TE fine-level classification conflicts only with Unknown, the known classification is retained (`KNOWN_VS_UNKNOWN_KEEP_KNOWN`); other unresolved broad-class conflicts are assigned Unknown. It uses only the Python standard library.

## Methods 2.6 — Pan-TE clustering

| Generic release script | Function | Release changes |
|---|---|---|
| `pante_clustering/01_split_unified_by_class.py` | Split the unified curated library into DNA, LTR, LINE, SINE and RC FASTA files | Input and output supplied by command-line arguments |
| `pante_clustering/02_run_pante_all_vs_all.sh` | Within-class all-versus-all BLASTn v2.15.0+ | Working directory supplied by `PANTE_WORKDIR`; job name matches script filename |
| `pante_clustering/03_build_pante_clusters.py` | Merge HSP intervals, apply high-confidence criteria and construct connected components | Working directory supplied by `--workdir` |
| `pante_clustering/04_audit_large_clusters.py` | Graph QC of large initial components | Working directory supplied by `--workdir`; audit does not alter membership |
| `pante_clustering/05_refine_pante_representative_clusters.py` | Representative-centered greedy refinement | Working directory supplied by `--workdir`; original representative selection and tie-breaking retained |

The clustering release copies remove hard-coded paths and environment activation commands. BLAST+ must be available in `PATH`.

## Methods 2.6 — genomic occupancy and parental/F1 comparison

| Generic release script | Function | Method retained |
|---|---|---|
| `genomic_occupancy/01_resolve_pante_abundance.py` | Map RepeatMasker consensus IDs to refined Pan-TE clusters and resolve all overlapping RepeatMasker annotations at base-pair resolution | All hits compete; non-`*` primary hits precede `*` secondary hits, followed by higher Smith–Waterman score and longer hit; each genomic base is assigned once |
| `genomic_occupancy/02_build_resolved_matrix.py` | Merge per-assembly Pan-TE occupancy tables and calculate assembly-normalized occupancy | Runtime manifest supplies sample ID, total assembly length and abundance-table path; `bp_per_Mb = occupied_bp / genome_size_bp × 1,000,000` |
| `genomic_occupancy/03_parental_f1_occupancy_contrast.py` | Identify parent-informative and parent-biased clusters and compare parental with corresponding F1 haplotype contrasts | Parent-informative: either parent has >=10 kb raw occupancy; parent bias: >=2-fold difference in bp/Mb; log2 contrasts use a 1 bp/Mb pseudocount; Spearman and sign concordance are calculated among parent-informative clusters; F1 haplotype bp/Mb direction is assessed and ties are recorded separately |

The occupancy resolver uses the RepeatMasker `.out` files generated with the unified v5.2 library for cross-assembly comparison. The matrix builder and parental/F1 comparison use external manifests; no assembly names or genome sizes are hard-coded in the release scripts. The matrix-builder script uses Python 3 standard-library modules; the parental/F1 comparison script uses Python 3, NumPy, pandas and SciPy. Exact package versions were not recorded in the source execution record.

## Evidence integration and final classification

The evidence-integration release includes the Dfam, RepBase and TEsorter family-summary generators, reciprocal LTR processing, Class II reciprocal filtering and evidence-master construction. RepBase and Dfam form one sequence-homology evidence dimension; TEsorter provides the domain dimension; LTR/MITE/Helitron provide structural evidence. The evidence master records the integrated dimension-level calls and conflicts but is not itself the final curated classification.

| Generic release script | Historical source basis | Function |
|---|---|---|
| `evidence_integration/01_filter_classII_reciprocal.sh` | Original Class II reciprocal filter | Applies query-HSP and subject-coverage criteria to MITE/Helitron evidence |
| `evidence_integration/02_rethreshold_LTR_reciprocal.py` | Original LTR reciprocal processing | Merges query/subject HSP intervals and applies reciprocal coverage thresholds |
| `evidence_integration/summarize_Dfam_v5.py` | Dfam family summary | Produces family-level Xenocyprididae/Cypriniformes homology calls |
| `evidence_integration/normalize_RepBase_v5.py` | RepBase family summary | Produces family-level RepBase homology calls |
| `evidence_integration/normalize_TEsorter_v5.py` | TEsorter class/domain integration | Produces family-level domain evidence |
| `evidence_integration/03_build_evidence_master.py` | Evidence-master construction | Combines homology, domain and structural dimensions |
| `evidence_integration/04_final_arbitration.py` | `final_arbitration_v5.py` | Strict evidence-based broad/fine classification and review outputs |
| `evidence_integration/05_operational_arbitration.py` | `operational_arbitration_v5.py` | Operational classification retaining prior labels unless source rules support rescue/override; retains strict result as a confirmation layer |
| `evidence_integration/06_relaxed_repbase_rescue_v5_1.py` | `multi_sample_relaxed_repbase_v5_1.py` | Prepares operational-Unknown queries and summarizes relaxed RepBase v5.1 results |
| `evidence_integration/07_relaxed_homology_rescue_v5_2.py` | `v52_relaxed_homology_rescue.py` | Prepares v5.1 remaining-Unknown queries and integrates relaxed RepBase/Dfam v5.2 results |

### Stage-specific rescue parameters

The historical v5.1 RepBase search used `-evalue 1e-3`, `-perc_identity 60`, `-word_size 7`, `-max_target_seqs 100` and 6 threads. Its family-level parser applies identity >=65%, merged query coverage >=30%, merged query-covered length >=100 bp and E-value <=1e-3. A broad-class assignment additionally requires at least two distinct winning-class subjects and a winning-class bit-score fraction >=0.70.

The historical v5.2 Dfam search used `-evalue 1e-3`, `-perc_identity 60`, `-word_size 7`, `-max_target_seqs 100` and 36 threads. Its database-specific parser applies identity >=65%, merged query coverage >=20%, merged query-covered length >=80 bp and E-value <=1e-3, considers the top 10 qualifying subjects by bit score, and requires at least two distinct winning-class subjects and a winning-class bit-score fraction >=0.70. RepBase and Dfam calls are reconciled independently: concordant passing calls are accepted; discordant passing calls remain Unknown; a single passing database can support rescue. Families already classified in v5.1 remain unchanged.

The manuscript's general description (identity >=65%, query-covered length >=80 bp and E-value <=1e-3) should not be read as saying the v5.1 and v5.2 parsers use identical query-coverage/minimum-length thresholds. These stage-specific operational filters are recorded separately above.

### Sample scope and production analysis

The manuscript Pan-TE analysis uses four focal assemblies: BSB, TC, BT-M (BTM) and BT-P (BTP). BTB-A (BTBA) and BTB-B (BTBB) are excluded from the manuscript analysis. This four-sample scope applies to the production classification, unified-library comparison, Pan-TE clustering and genomic-occupancy results described in the paper. Public release scripts use runtime manifests to specify these four inputs; auxiliary script variants or intermediate records do not redefine the production analysis scope.

The four-sample Pan-TE catalogue contains 19,132 classified consensus sequences and 11,950 refined clusters. The unified library contains 25,540 input consensus records and 25,130 exact-sequence-unique records before broad-class filtering. The genomic-occupancy comparison is based on the same four focal assemblies (BSB, TC, BTM and BTP).

A synthetic end-to-end test is provided at `tests/test_classification_release.py`. It exercises the four downstream classification/rescue scripts with synthetic inputs; it does not establish production-data reproducibility or biological classification accuracy.
