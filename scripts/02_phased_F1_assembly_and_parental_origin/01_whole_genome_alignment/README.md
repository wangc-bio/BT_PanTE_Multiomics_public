# Whole-genome alignment

MUMmer 4.0.0rc1 was used for whole-genome comparisons. The same command and thresholds are applied to all comparisons; only reference/query FASTA inputs and output prefix differ. Comparisons include BT-M and BT-P against the published BSB and TC references, plus the direct BT-M versus BT-P comparison.

Key settings:
- `nucmer -c 90 -l 40 -t 30`
- `delta-filter -1 -i 90 -o 90 -l 5000`
- `mummerplot` for dot-plot visualization; PostScript is converted to PDF with `ps2pdf`.

The `-1` option selects 1-to-1 alignments; it should not be described as a genome-wide guarantee that every aligned base is unique. Chromosome-level parental-origin assessment additionally considers sequence identity, alignment coverage and uniquely aligned sequence statistics. Dot plots are qualitative and are not, by themselves, quantitative ancestry estimates.

Run after loading MUMmer 4.0.0rc1 and Ghostscript:
```bash
export REF=/path/to/reference.fa
export QUERY=/path/to/query.fa
export PREFIX=comparison
export OUTDIR=/path/to/output
sbatch 01_MUMmer_alignment.sh
```
