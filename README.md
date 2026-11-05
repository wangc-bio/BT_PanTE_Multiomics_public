# Pan-TE catalogue and multi-omics analysis workflows

This repository contains the computational workflows and supporting methods for constructing a transposable-element (TE) catalogue and analysing TE genomic occupancy, transcription, DNA methylation and chromatin accessibility in a diploid intergeneric hybrid fish system.

The repository is organized into six numbered modules (01–06), with Module 04 divided into 04A and 04B. It is intended to document the methods and provide reusable, parameterized scripts; it is not a data repository and does not contain raw sequencing reads, genome assemblies or large intermediate matrices.

## Study scope

The manuscript-level Pan-TE catalogue and comparative genomic analysis use four focal assemblies:

- **BSB** — blunt snout bream (*Megalobrama amblycephala*);
- **TC** — topmouth culter (*Culter alburnus*);
- **BT-M** — the maternal F1 haplotype;
- **BT-P** — the paternal F1 haplotype.

BT-M and BT-P are analysed as the two phased haplotypes of the F1 hybrid. The BC1 backcross assemblies are not part of the four-assembly Pan-TE analysis. Scripts that accept manifests should be supplied with the intended four inputs explicitly.

## Workflow modules

| Module | Main contents | Manuscript methods |
|---|---|---|
| 01 | Read-level QC, mapping summaries and assembly-quality workflows | Supplementary methods |
| 02 | Whole-genome alignment, synteny, phylogenetic validation, parent-specific k-mer and depth validation | Genome assembly and parental-origin analyses |
| 03 | TE discovery, evidence integration, consensus-library curation, Pan-TE clustering and genomic occupancy | 2.5–2.6 |
| 04A | RepeatMasker Kimura-divergence calculation, normalized landscapes and parental/F1 contrasts | TE divergence analyses |
| 04B | Pan-TE RNA-seq, WGBS and ATAC-seq processing, quantification and regulatory-state analyses | Transcriptomic and epigenomic analyses |
| 05 | Cross-assay integration and statistical modelling | 2.11–2.12 |
| 06 | Static validation, software records, release checks and SHA-256 manifest generation | Reproducibility |

See [docs/analysis_workflow.md](docs/analysis_workflow.md) for the module-level description and [docs/reproducibility.md](docs/reproducibility.md) for input, parameter and provenance requirements.

## Repository layout

```text
.github/workflows/   automated static and synthetic validation
docs/                workflow, data-release and reproducibility notes
metadata/            blank metadata and file-manifest templates
scripts/             module-specific scripts, examples and methods
```

## Important interpretation boundaries

- Pan-TE clusters are sequence-based consensus-homology groups. They are not genomic insertion loci and should not be interpreted as population-level insertion genotypes.
- Kimura divergence measures sequence divergence from a consensus; it is not a direct measurement of absolute insertion age.
- Reference-level RepeatMasker occupancy does not by itself estimate population insertion frequency.
- ATAC-derived accessibility values are relative measurements when no external spike-in normalization is used; they do not establish genome-wide absolute chromatin opening or closing.
- Integrated associations across genomic occupancy, RNA expression, DNA methylation and chromatin accessibility are not, by themselves, evidence of causality.
- Individual, tissue, sequencing library and sequencing run are different units. Technical libraries from one individual must not be treated as independent biological replicates.

## Inputs and data availability

Analysis inputs must be obtained from the relevant sequence archives or supplied by the user. Large FASTQ/FQ files, BAM/CRAM files, genome assemblies, reference databases and large result matrices are intentionally excluded from this repository. Each analysis should be run with a manifest or environment configuration that records the exact input files and assembly versions.

No sample-level QC tables or private analysis evidence bundles are distributed here. Archive accessions and data DOIs should be taken from the final manuscript and its associated data records; they are not guessed or inferred in this repository.

## Validation

GitHub Actions provides:

1. Python syntax checks across the released scripts;
2. Bash syntax checks for shell and Slurm scripts;
3. R parsing when the CI environment is available;
4. a synthetic end-to-end test for the Module 03 classification workflow;
5. a release validation report and SHA-256 manifest as a downloadable workflow artifact.

Passing static checks does not mean that every biological analysis has been rerun. The required input datasets and software environments must still be supplied to reproduce data-dependent results.

## License and citation

No software license has been assigned yet. Until a license is added, the code remains subject to the default copyright rules; public visibility does not itself grant permission to redistribute or reuse it. Please contact the corresponding author regarding reuse.

Please cite this software repository using the metadata in [CITATION.cff](CITATION.cff). The software citation follows the author order supplied for the associated manuscript. Chong Wang and Bo Zhang contributed equally; Min Tao and Jing Luo are co-corresponding authors. The manuscript's formal preferred-citation record (title, journal, year and DOI) will be added when those bibliographic details are finalized.
