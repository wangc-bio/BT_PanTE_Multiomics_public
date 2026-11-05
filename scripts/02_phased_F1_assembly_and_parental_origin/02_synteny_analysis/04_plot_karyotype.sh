#!/usr/bin/env bash
#SBATCH --job-name=04_plot_karyotype
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
: "${CONDA_SH:?}"; : "${JCVI_ENV:?}"; : "${ANCHOR_DIR:?}"; : "${PLOT_CONFIG_DIR:?}"; : "${OUTDIR:?}"
mkdir -p "$OUTDIR" logs
source "$CONDA_SH"; conda activate "$JCVI_ENV"
cd "$ANCHOR_DIR"
for fmt in pdf svg; do
  python -m jcvi.graphics.karyotype "$PLOT_CONFIG_DIR/seqids" "$PLOT_CONFIG_DIR/layout" --notex -o "$OUTDIR/BT" --figsize 12x10 --format "$fmt" --dpi 600
done
