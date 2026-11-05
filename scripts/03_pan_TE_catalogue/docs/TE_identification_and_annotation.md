# TE identification and genome-wide annotation (Methods 2.5)

The same multi-stage TE identification and annotation workflow was applied independently to each focal genome assembly. One generic script is used for each software stage; genome FASTA files, consensus libraries, reference databases and output directories are provided through an external configuration file.

## De novo repeat discovery and consensus curation

RepeatModeler2 v2.0.5 was run using a genome-specific database generated with `BuildDatabase`. The recorded commands were:

```bash
BuildDatabase -name DB_NAME GENOME.fa
RepeatModeler -database DB_NAME -threads 36 -LTRStruct
```

The resulting consensus library was curated using TEtrimmer v1.7.2. The recorded command used the RepeatModeler2 consensus FASTA, the corresponding genome FASTA, a dedicated output directory, the Pfam database, and the following options:

```bash
python3 TEtrimmer.py \
    -i consensi.fa.classified \
    -g GENOME.fa \
    -o TEtrimmer_OUT \
    -t THREADS \
    --classify_unknown \
    --pfam_dir PFAM_DIR \
    -l TEtrimmer.log \
    -ll INFO
```

MCHelper v1.7.0 was then run using the TEtrimmer consensus FASTA and the Actinopterygii odb10 BUSCO HMM dataset:

```bash
python3 MCHelper.py \
    -r A \
    -t 24 \
    -l TEtrimmer_curated_consensus.fa \
    -g GENOME.fa \
    -o MCHelper_OUT \
    --input_type fasta \
    -b actinopterygii_odb10.all.hmm \
    -a F
```

The actual TEtrimmer output FASTA selected for MCHelper is supplied explicitly in the configuration file; the script does not assume an output filename.

## Structural evidence

LTR retrotransposon candidates were identified using LTR_FINDER_parallel v1.2, which internally uses LTR_FINDER v1.0.7, and LTRharvest from GenomeTools v1.6.2. LTR_FINDER_parallel was run with `-threads 24 -harvest_out -size 1000000 -time 300`. GenomeTools `suffixerator` was run with `-tis -suf -lcp -des -ssp -sds -dna`. LTRharvest used `-minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes`. The LTR_FINDER_parallel and LTRharvest candidate outputs were combined and processed using LTR_retriever v3.0.1 with `-genome GENOME.fa -inharvest rawLTR.scn -threads THREADS`.

MITE candidates were identified using the MITE-Hunter manager with `-P 1 -n 36 -c THREADS -S 12345678`. The MITE-Hunter version was not recorded in the available execution record and is not specified.

Helitron candidates were identified using HelitronScanner v1.1. Both forward and reverse-complement scans were performed. The recorded scan parameters were `-bs 1000000 -ol 50 -th 1 -tl THREADS`; paired candidates were identified using `-ht 5 -tt 10 -hlr 200:20000`, followed by candidate export with `draw --pure`.

The original RepeatModeler consensus sequences were queried by BLASTn against separate databases built from the MITE-Hunter and HelitronScanner candidate sequences. Supportive matches required identity >=70%, query coverage (`qcov_hsp`) >=50%, alignment length >=100 bp and E-value <=1e-10. High-confidence matches required identity >=80%, `qcov_hsp` >=80%, alignment length >=100 bp and E-value <=1e-20. The best qualifying hit per query was selected by bit score, then identity and query coverage. Subject-coverage filtering was not part of this recorded mapping procedure.

## RepBase consensus screening

The original RepeatModeler consensus libraries were queried against the pre-built RepBase release 31.07 fish/vertebrate clean database using BLASTn v2.15.0+. Supportive matches required identity >=70%, subject-level query coverage (`qcovs`) >=50%, alignment length >=100 bp and E-value <=1e-10. High-confidence matches required identity >=80%, `qcovs` >=80%, alignment length >=100 bp and E-value <=1e-20. For each RepeatModeler query, the best subject was retained separately for the supportive and high-confidence sets, ranked by bit score, identity and `qcovs`. RepBase subject identifiers were used to extract the annotated RepBase class and species fields.

## Domain-based and homology-based classification

Unclassified RepeatModeler consensus sequences were subjected to domain-based classification using TEsorter v1.5.1 with the `rexdb-metazoa` database. Unknown sequences were extracted before running TEsorter with `-p 4`.

Homology-based classification incorporated RepBase release 31.07 and the supplied Dfam/FamDB reference FASTAs. RepBase and Dfam are two reference sources within the same sequence-homology evidence dimension; they are not treated as two independent biological evidence dimensions. The final Dfam/FamDB release and the original study-specific reference-export command have not been verified in the available BT project record, so no Dfam release number is assigned here.

Unclassified consensus sequences were subjected to BLASTn v2.15.0+ searches. The final study-specific hit filters were sequence identity >=65%, subject-level query coverage (`qcovs`) >=20%, qcovs-derived query-covered length >=80 bp and E-value <=1 x 10^-3. A broad TE-class assignment additionally required at least two qualifying reference subjects and a winning-class bit-score fraction >=70%. Concordant RepBase and Dfam broad-class assignments could support classification, whereas unresolved conflicts were retained as `Unknown`.

The Dfam search used two taxonomic reference levels sequentially. The second level was applied to sequences remaining unresolved after the first level. The search scripts generate raw and hit-level filtered subject matches; they do not by themselves make the final broad-class assignment.

Classifications derived from RepeatModeler2, TEtrimmer and MCHelper were treated as prior or composite evidence. Structural, conserved-domain and sequence-homology signals were integrated as complementary evidence. RepBase and Dfam constitute one homology evidence dimension, rather than two independent biological dimensions. Conflicting or insufficiently supported classifications were conservatively retained as `Unknown` during evidence arbitration.

The historical workflow then applied two relaxed Unknown-rescue stages. The v5.1 RepBase stage used identity >=65%, query coverage >=30%, at least 100 query-covered bp and E-value <=1e-3; a broad-class rescue additionally required at least two qualifying subjects and a winning-class bit-score fraction >=70%. The v5.2 stage evaluated v5.1 remaining-Unknown families against RepBase and the curated Cyprinidae Dfam TE-only library. Its database-specific filters were identity >=65%, query coverage >=20%, at least 80 query-covered bp and E-value <=1e-3, with the top 10 qualifying subjects considered per database, at least two winner-class subjects and a winning-class bit-score fraction >=70%. When both databases independently passed but disagreed on broad class, the family remained `Unknown`; a passing call from only one database could be used when the other did not pass. Existing v5.1 classifications were retained unchanged in v5.2.

The strict and operational arbitration programs, followed by the v5.1 and v5.2 rescue programs, are separate from the raw BLAST search steps. The exact Dfam/FamDB release identifier has not been verified and is not assigned here. The v5.1 source runner excludes TC; the historical TC v5.1 result exists, but its exact generating script/run modification remains unverified. The complete downstream source programs have been located but are not yet released here as manifest-driven scripts.

## Genome-wide annotation

Genome-specific curated repeat libraries were used for the primary genome-wide annotation with RepeatMasker v4.1.8 and RMBlast v2.11.0+. The recorded parameters were:

```bash
RepeatMasker \
    -pa 9 \
    -lib CURATED_LIBRARY.fa \
    -gff \
    -html \
    -xsmall \
    -e rmblast \
    -no_is \
    GENOME.fa
```

The same RepeatMasker script and configuration are used for primary and sensitivity annotations by changing the library and output-directory variables. Genome-specific curated libraries are used for the primary annotation. A unified curated repeat library is additionally used for cross-assembly quantitative comparison/sensitivity analysis and supplies the consensus sequences for Pan-TE clustering. The genome-specific primary annotation and unified-library comparative annotation serve different purposes and should be reported separately.

Do not add `-nolow`, `-norna` or `-s` to the recorded command. The execution record also indicates that `RepeatMasker -version` and `RepeatMasker --version` are unsupported by the installed build, so these should not be used as version-check commands.