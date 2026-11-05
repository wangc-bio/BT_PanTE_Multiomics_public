#!/usr/bin/env bash
#SBATCH --job-name=01_MUMmer_alignment
#SBATCH --cpus-per-task=30
#SBATCH --mem=150G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail

# MUMmer 4.0.0rc1; same parameters are used for every reference/query pair.
# Set REF, QUERY, PREFIX and OUTDIR before submission, e.g. via an sbatch wrapper.
: "${REF:?Set REF to the reference FASTA}"
: "${QUERY:?Set QUERY to the query FASTA}"
: "${PREFIX:?Set PREFIX for this comparison}"
: "${OUTDIR:?Set OUTDIR for outputs}"
THREADS=${SLURM_CPUS_PER_TASK:-30}
mkdir -p "$OUTDIR" logs
cd "$OUTDIR"

nucmer -c 90 -l 40 -t "$THREADS" -p "$PREFIX" "$REF" "$QUERY"
delta-filter -1 -i 90 -o 90 -l 5000 "${PREFIX}.delta" > "${PREFIX}.filtered.delta"
mummerplot -p "${PREFIX}.filtered" -f "${PREFIX}.filtered.delta" -t postscript
ps2pdf "${PREFIX}.filtered.ps" "${PREFIX}.filtered.pdf"
