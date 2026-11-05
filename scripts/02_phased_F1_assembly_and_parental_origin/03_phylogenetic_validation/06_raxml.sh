#!/usr/bin/env bash
#SBATCH --job-name=06_raxml
#SBATCH --cpus-per-task=36
#SBATCH --mem=100G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
: "${RAXML:?Path to RAxML 8.2.8 executable}"; : "${ALIGNMENT_DIR:?}"; : "${OUTGROUP_PREFIX:?Outgroup taxon prefix, e.g. Cide}"
mkdir -p logs
for n in $(seq -w 1 24); do
  chr=chr$n; indir="$ALIGNMENT_DIR/$chr"; out="$indir/trees"
  mkdir -p "$out"; cd "$out"
  "$RAXML" -T 36 -f a -x 12345 -p 12345 -# 1000 -m PROTGAMMAJTTF \
    -o "${OUTGROUP_PREFIX}${chr}" -s "$indir/all.phy" -n "BT${chr}"
done
