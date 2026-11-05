#!/usr/bin/env bash
#SBATCH --job-name=02_muscle_alignment
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-24%24
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
: "${MUSCLE:?Path to MUSCLE v3 executable}"; : "${ORTHO_DIR:?}"
: "${ALIGN_DIR:?}"
CHR=$(printf "chr%02d" "$SLURM_ARRAY_TASK_ID")
INDIR="$ORTHO_DIR/$CHR/Single_Copy_Orthologue_Sequences"
OUTDIR="$ALIGN_DIR/$CHR"
mkdir -p "$OUTDIR" logs
shopt -s nullglob
for fa in "$INDIR"/*.fa; do
  "$MUSCLE" -in "$fa" -out "$OUTDIR/$(basename "$fa").muscle"
done
