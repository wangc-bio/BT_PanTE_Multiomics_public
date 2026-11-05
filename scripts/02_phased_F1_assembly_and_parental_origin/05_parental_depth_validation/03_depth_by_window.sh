#!/usr/bin/env bash
#SBATCH --job-name=03_depth_by_window
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --array=0-35
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
: "${SAMTOOLS:?}"; : "${WINBED:?}"; : "${BAM_LIST:?}"; : "${OUTDIR:?}"
WIN=${WIN:-1000000}
mkdir -p "$OUTDIR/depth_by_window" "$OUTDIR/genome_mean" "$OUTDIR/ND" logs
bam=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$BAM_LIST")
[[ -s "$bam" ]] || { echo "Missing BAM for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
prefix=$(basename "$bam" .filtered.bam)
mean="$OUTDIR/depth_by_window/$prefix.win$WIN.meanDepth.tsv"
gm="$OUTDIR/genome_mean/$prefix.genome_mean.txt"
nd="$OUTDIR/ND/$prefix.win$WIN.ND.tsv"
"$SAMTOOLS" bedcov "$WINBED" "$bam" | awk 'BEGIN{OFS="\t"}{len=$3-$2; print $1,$2,$3,(len>0?$4/len:0)}' > "$mean"
awk 'BEGIN{s=0;n=0}{l=$3-$2;s+=$4*l;n+=l}END{print(n>0?s/n:0)}' "$mean" > "$gm"
GM=$(cat "$gm")
awk -v gm="$GM" 'BEGIN{OFS="\t"}{print $1,$2,$3,($4/gm)}' "$mean" > "$nd"
