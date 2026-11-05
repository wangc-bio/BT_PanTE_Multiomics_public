# Module 02 — Phased F1 genome and parental-origin validation

## Scope

This module contains the public analysis scripts used to characterize the two phased BT F1 haplotypes, assign parental-origin support, and validate the haplotype/gene relationships.

- BT-M: phased BT maternal haplotype
- BT-P: phased BT paternal haplotype
- BSB: *Megalobrama amblycephala*
- TC: *Culter alburnus*

## Analysis components

| Component | Software / method | Main purpose |
|---|---|---|
| 01 Whole-genome alignment | MUMmer 4.0.0rc1 | Compare each BT haplotype with the BSB and TC reference assemblies |
| 02 Synteny analysis | JCVI/MCscan | Compare collinear gene blocks among BSB, TC, BT-M and BT-P |
| 03 Phylogenetic validation | OrthoFinder, MUSCLE, Gblocks, SeqKit, ProtTest, RAxML | Infer chromosome-wise phylogenetic relationships |
| 04 Parent-specific k-mer | Jellyfish 2.3.1, Bowtie 1.3.1, BEDTools | Assign parental-source support in genomic windows and chromosomes |
| 05 Parental depth validation | BWA, SAMtools, BEDTools, Python | Independently compare normalized parental read depth with k-mer ancestry calls |
| 06 Homologous gene identification | MMseqs2 easy-rbh | Identify reciprocal-best-hit protein pairs between BT-M and BT-P |

## Key analysis settings

### 1. Whole-genome alignment

The same MUMmer commands and thresholds are used for the four comparisons (BT-M/BT-P against BSB/TC); only the reference and query FASTA inputs and output prefix change.

```bash
nucmer -c 90 -l 40 -t 30 -p PREFIX REF.fa QUERY.fa
delta-filter -1 -i 90 -o 90 -l 5000 PREFIX.delta > PREFIX.filtered.delta
mummerplot -p PREFIX.filtered -f PREFIX.filtered.delta -t postscript
```

The dot plots are qualitative. Chromosome-level parental-origin assessment also considers sequence identity, alignment coverage and uniquely aligned sequence statistics; the dot plots alone are not quantitative ancestry estimates.

### 2. Synteny analysis

Longest-isoform GFF3/CDS resources for BSB, TC, BT-M and BT-P are used. mRNA features are converted to BED with the ID attribute, duplicate BED entries are removed, and corresponding CDS sequences are extracted. Six pairwise JCVI/MCscan comparisons are run: BSB–BT-M, BSB–BT-P, BSB–TC, BT-M–BT-P, BT-M–TC and BT-P–TC. The MCscan C-score is 0.80; anchor filtering uses `--minspan=30 --simple`. Karyotype plots use the study-specific `seqids` and `layout` files.

### 3. Phylogenetic validation

The chromosome-wise analysis covers chr01–chr24 and seven taxa: Mamb, Ppek, Calb, Cery, Cide, BT-M and BT-P. OrthoFinder uses DIAMOND search (`-S diamond`) and six threads. MUSCLE v3 is used for protein alignment. Gblocks retains the original filtering settings `-b4=5 -b5=h -t=p`. SeqKit sorts the Gblocks alignments and converts them to single-line FASTA; the alignments are combined using the original paste/sed step. Taxon sequence IDs are then manually standardized and the resulting matrix is inspected before downstream analysis. ProtTest uses `-all-distributions -F -AIC -BIC -tc 0.5 -threads 36`. RAxML 8.2.8 uses `-f a -x 12345 -p 12345 -# 1000 -m PROTGAMMAJTTF`, with the chromosome-specific Cide outgroup.

### 4. Parent-specific k-mer assignment

Jellyfish counts canonical k-mers from the BSB and TC reference genome FASTA files (k=31 and k=21); the ancestry workflow uses k=31. K-mers with parental copy counts from 1 to 5 are retained, and set subtraction identifies parent-specific k-mers. Bowtie 1.3.1 maps them to the BT assembly with exact matching and unique-hit reporting (`-v 0 -a -m 1 --best --strata`). BEDTools counts mapped k-mers in 1-Mb windows. Windows with fewer than 50 total hits are Indistinct; BSB_source and TC_source calls require a parental fraction of at least 0.70; remaining supported windows are Mixed. Chromosome-level summaries classify chromosomes as BSBChr, TCChr, RecChr or IndistinctChr using the original mean-fraction and window-call thresholds. The 5-Mb analysis reuses the same mapped k-mer BED files and applies the same classification thresholds to primary chromosomes.

These k-mers are derived from parental reference genome FASTA files, not parental resequencing FASTQ files. Resequencing data are used in the separate depth-validation component.

### 5. Parental read-depth validation

Parental resequencing reads are mapped to the BT reference with BWA-MEM. BAM filtering retains MAPQ >=20 and properly paired reads; duplicate reads are retained by default, matching the original `REMOVE_DUP=0` setting. Mean depth is calculated in the selected 1-Mb or 5-Mb windows using SAMtools bedcov and normalized by each sample's length-weighted genome-wide mean depth.

The Python validation script calculates group mean ND, `log2((meanND_BSB + 1e-6)/(meanND_TC + 1e-6))`, a two-sided Welch t-test, and Benjamini–Hochberg FDR. The wrapper uses FDR 0.05 and absolute log2FC 0.3. Depth calls are compared with k-mer calls by exact chromosome/start/end coordinates; outputs include per-window calls, the merged comparison table, overall agreement and a confusion matrix. The same Python script accepts either 1-Mb or 5-Mb inputs.

### 6. Homologous gene pairs

MMseqs2 `easy-rbh` is run between BT-M and BT-P longest-isoform protein FASTA files with sensitivity 7.5, E-value 1e-5, minimum sequence identity 0.5, coverage 0.6 and `--cov-mode 0`. The output contains query/target IDs, identity, alignment and sequence lengths, coverage, E-value and bit score. RBH pairs should not automatically be equated with phylogenetically validated one-to-one orthologues.

## Repository contents and reproducibility

Each component has its own README and parameterized scripts. Set executable paths, input files and output directories through the documented environment variables before submission. Software versions are stated only where verified from the source material. Results already presented in the manuscript or supplementary materials are not duplicated in this repository.

See [SOURCE_SCRIPT_INVENTORY.md](SOURCE_SCRIPT_INVENTORY.md) for the source-to-public-script inventory. The analyses have already been completed by the study team and the corresponding results have been obtained. The scripts have been checked against the original commands and parameters to document those completed analyses; an independent end-to-end rerun is not required for this repository release.
