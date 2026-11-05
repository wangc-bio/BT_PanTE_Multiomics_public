# Module 04A — Kimura-divergence analysis of TE landscapes

## Scope

This module documents the RepeatMasker-based Kimura-divergence analysis of transposable-element (TE) landscapes in four genome assemblies:

- BSB: blunt snout bream parental genome
- TC: topmouth culter parental genome
- BT-M (BTM): BSB-derived F1 haplotype
- BT-P (BTP): TC-derived F1 haplotype

All four genomes were annotated for this analysis with RepeatMasker v.4.1.8 using the same unified TE library (v5.2). The analysis uses the RepeatMasker alignment output to calculate CpG-adjusted Kimura divergence and compare TE coverage distributions across divergence bins.

This module covers **divergence-landscape comparisons**, not insertion-site presence/absence genotyping. It does not estimate population-level insertion frequencies or prove inheritance of individual TE insertion loci.

## Workflow

### 1. RepeatMasker alignment generation

Script: `run_RepeatMasker_unified_v52_Kimura_align.sh`

RepeatMasker v.4.1.8 is run independently for BSB, TC, BTM and BTP using the unified v5.2 TE library. The same annotation settings are retained; the `-a` option is included to generate the alignment file required by `calcDivergenceFromAlign.pl`.

Core RepeatMasker options used in the source workflow:

- `-pa 9`
- `-lib` unified v5.2 library
- `-gff -html -xsmall`
- `-e rmblast -no_is -a`

The script expects the four genome FASTA files and unified library to be available at the configured paths. Environment activation is configurable through `CONDA_SH`, `CONDA_ENV` and `RMBLAST_PROFILE`; the RepeatMasker executable can be specified with `REPEATMASKER_EXECUTABLE`. Set `TE_BASE` and, if needed, `TE_LIBRARY` and `KIMURA_RUN_DIR` before submission. Create the relative `logs/` directory before submitting the Slurm job because Slurm opens its output files before the script body executes.

### 2. Kimura-divergence calculation

The RepeatMasker utility `calcDivergenceFromAlign.pl` is used on the four `.align` files to generate the corresponding `.divsum` summaries. The divergence values are CpG-adjusted, as recorded in the generated summary files.

Script: `run_calc_Kimura_divergence.slurm`

The original calculation wrapper has now been added. It runs as a four-task Slurm array (BSB, TC, BTM and BTP), checks each input `.align` file and the RepeatMasker utility, and calls:

```bash
perl calcDivergenceFromAlign.pl -s SAMPLE.divsum SAMPLE.align
```

The public script parameterizes the analysis root and RepeatMasker utility directory through `KIMURA_BASE` and `REPEATMASKER_UTIL`. Before submission, export both variables and create the `logs/` directory under the submission working directory, because Slurm opens the log files before the script body runs. The output for each sample is written to `${KIMURA_BASE}/Kimura_divergence/SAMPLE/SAMPLE.divsum`.

Example submission (replace the two placeholder paths with the locations in the execution environment):

```bash
mkdir -p logs
export KIMURA_BASE=/path/to/RepeatMasker_unified_v5_2_Kimura_align
export REPEATMASKER_UTIL=/path/to/RepeatMasker/util
sbatch run_calc_Kimura_divergence.slurm
```

The array submits one task for each of BSB, TC, BT-M and BT-P. The four `.divsum` files are the direct outputs consumed by `parse_Kimura_divsum.py`.

### 3. Parse divsum summaries

Script: `parse_Kimura_divsum.py`

The parser reads the BSB, TC, BTM and BTP `.divsum` files, extracts family-level Kimura summaries and divergence-bin coverage matrices, maps RepeatMasker classes to the five broad classes (DNA, LTR, LINE, SINE and RC), and records whether each summary indicates CpG adjustment.

It produces parsed family-level, RepeatMasker-class, broad-class and main-five-class tables, together with parsing QC summaries. It also checks whether the four samples share the same Kimura-bin set.

### 4. Genome-normalized TE landscapes and parental/F1 contrasts

Script: `analyze_Kimura_landscape.py`

Coverage in each Kimura-divergence bin is normalized to the corresponding assembly size:

`Genome_pct = 100 × Coverage_bp_wellChar / Assembly_bp`

The five broad TE classes are analyzed separately and together. For each class and for the combined five-class landscape, the analysis calculates:

- BSB versus BT-M and TC versus BT-P landscape correlations (Spearman);
- correlation between parental and F1 haplotypic contrast profiles (Spearman and Pearson);
- the linear slope relating parental and F1 contrast profiles;
- sums of absolute parental and F1 contrast magnitudes.

The parental contrast is defined as BSB minus TC; the F1 haplotypic contrast is defined as BT-M minus BT-P. Positive and negative values therefore retain the direction of the corresponding contrast.

### 5. Non-zero-bin sensitivity analysis

Script: `analyze_Kimura_nonzero_sensitivity.py`

The full-bin correlations are recalculated after excluding bins where both members of a comparison have zero coverage. The script reports full and non-zero-bin Spearman/Pearson correlations, the number of bins removed, and the corresponding changes. It evaluates BSB versus BT-M, TC versus BT-P, and parental-contrast versus F1-contrast profiles for each of the five broad classes and the combined five-class landscape.

### 6. LTR parental-contrast contribution analysis

Script: `analyze_LTR_major_contrast_bins.py`

For LTR retrotransposons, Kimura bins are ranked by the absolute parental contrast. Each bin's contribution to the sum of absolute parental contrasts and its cumulative contribution are calculated. The script summarizes bins reaching cumulative contribution targets of 50%, 80% and 90%, and evaluates whether the parental and F1 contrast directions agree within those ranked bins.

### 7. Visualization

Script: `plot_Kimura_landscape_3panel_v5.py`

The figure summarizes:
- Panel A: five-class stacked TE coverage landscapes for BSB, BT-M, TC and BT-P;
- Panel B: parental (BSB minus TC) and F1 haplotypic (BT-M minus BT-P) coverage contrasts;
- Panel C: LTR divergence profiles for the two parental/haplotype pairs.

The plotting script exports PDF, SVG and 600-dpi PNG files. The displayed plotting range is Kimura divergence 0–50; the upstream landscape statistics are calculated from the parsed divergence-bin data.

## Interpretation

Kimura divergence is treated as a relative measure of sequence divergence from the RepeatMasker consensus, not as a direct estimate of absolute TE insertion age. Similarity between parental and corresponding F1-haplotype landscapes is interpreted as retention of divergence-distribution patterns, not as proof that individual insertion loci were inherited.

## Reproducibility notes

- Software/library: RepeatMasker v.4.1.8; unified TE library v5.2; RMBlast is invoked as the search engine.
- Inputs: BSB and TC parental genome FASTAs, BT-M and BT-P haplotype genome FASTAs, the unified TE library, four RepeatMasker `.align` files, four `.divsum` files, and assembly-size table.
- Python scripts use `KIMURA_BASE` or `KIMURA_DIR` environment variables to set analysis locations. The plotting font can be supplied through `FONT_FILE`; otherwise the script falls back to Arial.
- Python dependencies: NumPy, pandas and Matplotlib (for plotting).
- The manuscript and supplementary result tables/figures are not duplicated here.
- The original `run_calc_Kimura_divergence.slurm` wrapper is included in this directory. Set `KIMURA_BASE` to the directory containing the RepeatMasker Kimura run folder and `REPEATMASKER_UTIL` to the RepeatMasker `util` directory. Create `logs/` in the submission working directory before `sbatch`, because Slurm opens output files before the job script starts.
