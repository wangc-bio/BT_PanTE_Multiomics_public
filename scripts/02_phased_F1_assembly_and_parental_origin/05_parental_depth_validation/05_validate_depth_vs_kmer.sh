#!/usr/bin/env bash
#SBATCH --job-name=05_validate_depth_vs_kmer
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
: "${PYTHON:?}"; : "${VALIDATION_PY:?Path to 05_depth_validate_nd.py}"
: "${ND_MATRIX:?}"; : "${BSB_LIST:?}"; : "${TC_LIST:?}"; : "${KMER_SOURCE:?}"; : "${OUT_PREFIX:?}"
"$PYTHON" "$VALIDATION_PY" --nd_matrix "$ND_MATRIX" --bsb_list "$BSB_LIST" --tc_list "$TC_LIST" --kmer_source "$KMER_SOURCE" --out_prefix "$OUT_PREFIX" --fdr 0.05 --log2fc 0.3 --only_chr
