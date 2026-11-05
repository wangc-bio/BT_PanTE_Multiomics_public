#!/usr/bin/env bash
#SBATCH --job-name=01_window_definitions
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
: "${BT_FASTA:?}"; : "${SAMTOOLS:?}"; : "${BEDTOOLS:?}"; : "${OUTDIR:?}"
mkdir -p "$OUTDIR" logs; cd "$OUTDIR"
"$SAMTOOLS" faidx "$BT_FASTA"
cut -f1,2 "${BT_FASTA}.fai" > BT.genome
awk '$1 ~ /^chr[0-9][0-9][AB]$/' BT.genome > BT.primary.genome
for win in 1000000 5000000; do
  "$BEDTOOLS" makewindows -g BT.primary.genome -w "$win" > "win.$win.bed"
  "$BEDTOOLS" sort -g BT.primary.genome -i "win.$win.bed" > "win.$win.sorted.bed"
done
