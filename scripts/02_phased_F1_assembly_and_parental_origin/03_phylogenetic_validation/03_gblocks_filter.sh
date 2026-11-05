#!/usr/bin/env bash
#SBATCH --job-name=03_gblocks_filter
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
: "${GBLOCKS:?Path to Gblocks executable}"; : "${ALIGN_DIR:?}"; : "${GBLOCK_DIR:?}"
mkdir -p logs
for n in $(seq -w 1 24); do
  chr=chr$n
  in="$ALIGN_DIR/$chr"
  out="$GBLOCK_DIR/$chr"
  mkdir -p "$out"
  shopt -s nullglob
  for fa in "$in"/*.muscle; do
    cp -f "$fa" "$out/"
    base=$(basename "$fa")
    (cd "$out" && "$GBLOCKS" "$base" -b4=5 -b5=h -t=p -e=.gblock)
  done
done
