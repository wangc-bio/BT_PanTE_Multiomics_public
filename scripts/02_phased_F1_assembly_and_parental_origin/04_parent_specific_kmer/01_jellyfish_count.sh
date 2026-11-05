#!/usr/bin/env bash
#SBATCH --job-name=01_jellyfish_count
#SBATCH --cpus-per-task=36
#SBATCH --mem=160G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
: "${JELLYFISH:?Path to Jellyfish 2.3.1}"; : "${BSB_FASTA:?}"; : "${TC_FASTA:?}"; : "${OUTDIR:?}"
THREADS=${SLURM_CPUS_PER_TASK:-36}
mkdir -p "$OUTDIR" logs
cd "$OUTDIR"
for k in 31 21; do
  size=2G; [[ "$k" == 21 ]] && size=3G
  "$JELLYFISH" count -C -m "$k" -s "$size" -t "$THREADS" -o "BSB.k$k.jf" "$BSB_FASTA"
  "$JELLYFISH" count -C -m "$k" -s "$size" -t "$THREADS" -o "TC.k$k.jf" "$TC_FASTA"
  "$JELLYFISH" dump -c "BSB.k$k.jf" > "BSB.k$k.dump"
  "$JELLYFISH" dump -c "TC.k$k.jf" > "TC.k$k.dump"
  "$JELLYFISH" stats "BSB.k$k.jf" > "BSB.k$k.stats.txt"
  "$JELLYFISH" stats "TC.k$k.jf" > "TC.k$k.stats.txt"
done
