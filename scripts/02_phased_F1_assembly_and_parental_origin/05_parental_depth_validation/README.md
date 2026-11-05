# Parental read-depth validation

The depth-based validation compares normalized read depth (ND) between BSB and TC parental samples in each genomic window, then evaluates agreement with the parent-specific k-mer ancestry calls. Window size is determined by the matched ND matrix and k-mer source table (e.g. 1 Mb or 5 Mb); the Python script itself is window-size agnostic.

## Statistical procedure

1. Read the ND matrix and obtain BSB/TC sample columns from the corresponding BAM-list basenames.
2. Optionally retain only primary chromosome identifiers matching `chrXXA` or `chrXXB` with `--only_chr`.
3. Calculate the mean ND for each parental group in every window.
4. Calculate `log2FC = log2((meanND_BSB + 1e-6) / (meanND_TC + 1e-6))`.
5. Perform a two-sided Welch two-sample t-test for BSB versus TC ND in each window (`equal_var=False`, missing values omitted).
6. Adjust the resulting p-values using the Benjamini–Hochberg procedure.
7. Assign `BSB` when FDR is below the selected cutoff and log2FC is at least the positive threshold; assign `TC` when FDR is below the cutoff and log2FC is at most the negative threshold; otherwise assign `Mixed`.
8. Join depth calls to k-mer calls by exact `chr`, `start` and `end` coordinates, then report the overall agreement and a depth-call versus k-mer-call count table.

## Parameters used in the study wrapper

The 1-Mb wrapper supplies `--fdr 0.05 --log2fc 0.3 --only_chr`. The Python CLI defaults are FDR 0.05 and absolute log2FC 0.5; the explicit command-line value overrides that default. The same Python implementation accepts either 1-Mb or 5-Mb matched input files.

## Outputs

- `<out_prefix>.depth_validate.tsv`: group mean ND, log2FC, p-value, BH-adjusted FDR and depth-based call for each window.
- `<out_prefix>.depth_vs_kmer.tsv`: depth results joined with the k-mer source call.
- `<out_prefix>.summary.txt`: sample counts, thresholds, number of matched windows, agreement rate and confusion-matrix counts.

## Dependencies

Python 3 with NumPy, pandas, SciPy and statsmodels. The shell wrapper activates the `py3.7` conda environment. The script uses the BAM-list file names to match sample columns in the ND matrix; the basenames (without `.bam`) must therefore match the matrix header exactly.

Run through the supplied `05_validate_depth_vs_kmer.sh` wrapper after setting the paths and environment variables for the ND matrix, BSB/TC BAM lists, k-mer source table and output prefix.