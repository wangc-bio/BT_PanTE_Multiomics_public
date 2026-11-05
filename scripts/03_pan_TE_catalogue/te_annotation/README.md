# TE annotation execution scripts

## Design

Each software stage has one reusable script. The same scripts are used for all four focal genome assemblies; only the external configuration changes. No assembly-specific script copies, hard-coded sample names, or private server paths are included.

Copy `te_annotation.env.example` to a private configuration file and edit the generic input/output paths. Export its location before submitting jobs:

```bash
export TE_CONFIG=/path/to/te_annotation.env
```

`GENOME_FASTA`, `GENOME_ID` and `WORK_DIR` identify the current run. Reference databases and their prefixes can be shared across genome runs. The configuration template deliberately uses generic filenames; set each input to the actual output produced by the preceding step.

## Execution order

Run de novo discovery and consensus curation in order:

```bash
sbatch 01_repeatmodeler2.sh
# After RepeatModeler2 finishes, set RM2_CLASSIFIED_FASTA in TE_CONFIG.
sbatch 02_tetrimmer.sh
# Set MCHelper input to the actual TEtrimmer consensus FASTA in TE_CONFIG.
sbatch 03_mchelper.sh
```

Structural evidence can be generated independently after the genome FASTA and tool paths have been checked:

```bash
sbatch 04_ltr_structural_annotation.sh
sbatch 05_mite_hunter.sh
sbatch 06_helitronscanner.sh
```

Domain and Class II evidence steps:

```bash
sbatch 07_tesorter.sh
# Run after the RepeatModeler2, MITE-Hunter and HelitronScanner outputs are available.
sbatch 08_classII_evidence_mapping.sh
```

Run the initial RepBase screen against the original RepeatModeler consensus library:

```bash
sbatch 09_repbase_consensus_mapping.sh
```

This step uses the pre-built RepBase release 31.07 fish/vertebrate clean database. It reports supportive and high-confidence matches separately.

Build the two Dfam nucleotide databases once from the supplied curated reference FASTAs; reuse these database prefixes across genome runs. Level 1 represents Xenocyprididae and level 2 represents Cypriniformes:

```bash
sbatch 10_build_homology_databases.sh
```

Set `HOMOLOGY_QUERY_FASTA` to the unresolved consensus set after prior/composite, structural and domain evidence has been reviewed. Then run the RepBase and first-level Dfam homology-rescue searches:

```bash
sbatch 11_homology_rescue_search.sh
```

The second Dfam reference level is run only when `DFAM_LEVEL2_QUERY_FASTA` is set to the sequences remaining unresolved after the first Dfam level. The script does not infer or generate this FASTA; it must come from the evidence-integration/query-removal step.

After the final genome-specific curated repeat library has been prepared, run primary RepeatMasker annotation:

```bash
sbatch 12_repeatmasker.sh
```

To perform the unified-library sensitivity annotation, reuse `12_repeatmasker.sh` and change only `REPEATMASKER_LIBRARY_FASTA` and `REPEATMASKER_OUT_DIR` in the external configuration.

Submit each SLURM script from this directory. `.out` and `.err` files are written beside the script. The scripts retain the analysis thread counts and software parameters, but do not prescribe cluster partitions, memory limits or wall-time limits. Software environment activation is left to the user's execution environment.

## Script inventory

| Script | Purpose | Main software / recorded parameters |
|---|---|---|
| `01_repeatmodeler2.sh` | De novo repeat discovery | RepeatModeler2 v2.0.5; `BuildDatabase`; `-threads 36 -LTRStruct` |
| `02_tetrimmer.sh` | Automated consensus curation | TEtrimmer v1.7.2; `--classify_unknown`; Pfam-A HMM database; `-ll INFO` |
| `03_mchelper.sh` | Additional consensus curation | MCHelper v1.7.0; `-r A -t 24 --input_type fasta -a F`; Actinopterygii odb10 HMM |
| `04_ltr_structural_annotation.sh` | LTR candidate discovery and curation | LTR_FINDER_parallel v1.2 (using LTR_FINDER v1.0.7); GenomeTools v1.6.2; LTR_retriever v3.0.1 |
| `05_mite_hunter.sh` | MITE candidate discovery | MITE-Hunter manager; `-P 1 -n 36 -c 36 -S 12345678`; version not recorded |
| `06_helitronscanner.sh` | Helitron candidate discovery | HelitronScanner v1.1; forward and reverse-complement scans |
| `07_tesorter.sh` | Domain-based classification of Unknown consensus sequences | TEsorter v1.5.1; SeqKit version not recorded; `rexdb-metazoa`; 4 threads |
| `08_classII_evidence_mapping.sh` | Query original RepeatModeler consensus sequences against separate MITE-Hunter and HelitronScanner candidate databases | BLASTn v2.15.0+ (confirmed); supportive and high-confidence `qcov_hsp` filters |
| `09_repbase_consensus_mapping.sh` | Search original RepeatModeler consensus sequences against RepBase release 31.07 | BLASTn v2.15.0+; supportive: identity >=70%, qcovs >=50%, alignment >=100 bp, E <=1e-10; high confidence: identity >=80%, qcovs >=80%, alignment >=100 bp, E <=1e-20 |
| `10_build_homology_databases.sh` | Build nucleotide BLAST databases from supplied Xenocyprididae and Cypriniformes Dfam reference FASTAs | BLAST+ v2.15.0+; `makeblastdb -dbtype nucl` |
| `11_homology_rescue_search.sh` | RepBase/Dfam homology-rescue searches and hit-level filtering | BLASTn v2.15.0+; identity >=65%, qcovs >=20%, qcovs-derived query-covered length >=80 bp, E <=1e-3 |
| `12_repeatmasker.sh` | Genome-wide repeat annotation | RepeatMasker v4.1.8; RMBlast v2.11.0+; `-pa 9 -lib -gff -html -xsmall -e rmblast -no_is` |

## Evidence integration and reproducibility boundary

The RepBase consensus-screening and homology-rescue search scripts produce hit-level BLAST evidence. The family-level normalizers in `evidence_integration/` summarize Dfam, RepBase and TEsorter evidence and record class-level conflicts. These summaries feed the integrated evidence master; they are not themselves the final curated classification. The downstream strict/operational arbitration and v5.1/v5.2 Unknown-rescue steps remain separate from the search and family-summary stages and still require generic release scripts.

The Dfam reference FASTAs are inputs to the database-build step. The exact Dfam/FamDB release and final reference-export provenance for the study have not been verified in the available BT project record; no release number is assigned in this repository.

The unified repeat library is constructed by 100% exact-sequence MD5 deduplication only; no near-identity clustering is used. The generic unified-library builder is provided in [the unified-library workflow](../unified_library/README.md). It accepts any number of input libraries through a runtime manifest and does not embed production sample identifiers or server paths.

Raw FASTQ, genome assemblies, consensus libraries, reference FASTAs, BLAST databases and large intermediate outputs remain outside GitHub.