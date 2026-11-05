#!/usr/bin/env bash
#SBATCH --job-name=01_orthofinder_chr
#SBATCH --cpus-per-task=6
#SBATCH --mem=60G
#SBATCH --array=1-24%24
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
: "${ORTHOFINDER:?Path to OrthoFinder executable}"
: "${PROTEOME_DIR:?Directory containing chr01–chr24 subdirectories with the seven proteomes}"
: "${OUT_BASE:?Output directory}"
CHR=$(printf "chr%02d" "$SLURM_ARRAY_TASK_ID")
THREADS=${SLURM_CPUS_PER_TASK:-6}
mkdir -p "$OUT_BASE/$CHR" logs
# Expected taxa: Mamb, Ppek, Calb, Cery, Cide, BT-M and BT-P.
# The chromosome-specific FASTA files in PROTEOME_DIR must use these prefixes.
cp "$PROTEOME_DIR/$CHR"/*.pep "$OUT_BASE/$CHR"/
"$ORTHOFINDER" -f "$OUT_BASE/$CHR" -S diamond -t "$THREADS"
