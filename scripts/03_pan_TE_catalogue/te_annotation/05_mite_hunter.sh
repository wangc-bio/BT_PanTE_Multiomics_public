#!/bin/bash
# Software version: MITE-Hunter version not recorded in the source execution record.
#SBATCH --job-name=05_mite_hunter
#SBATCH --cpus-per-task=36
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

: "${TE_CONFIG:?Export TE_CONFIG=/path/to/te_annotation.env before submission}"
source "${TE_CONFIG}"

: "${GENOME_FASTA:?GENOME_FASTA is required}"
: "${WORK_DIR:?WORK_DIR is required}"
: "${MITE_HUNTER_MANAGER:?MITE_HUNTER_MANAGER is required}"

THREADS="${SLURM_CPUS_PER_TASK:-36}"
OUT_DIR="${WORK_DIR}/05_mite_hunter"
RUN_PREFIX="MITE_P1"

[[ -s "${GENOME_FASTA}" ]] || { echo "ERROR: genome FASTA is missing: ${GENOME_FASTA}" >&2; exit 1; }
[[ -s "${MITE_HUNTER_MANAGER}" ]] || { echo "ERROR: MITE_Hunter_manager.pl is missing: ${MITE_HUNTER_MANAGER}" >&2; exit 1; }

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

perl "${MITE_HUNTER_MANAGER}" \
    -i "${GENOME_FASTA}" \
    -g "${RUN_PREFIX}" \
    -P 1 \
    -n 36 \
    -c "${THREADS}" \
    -S 12345678

FINAL_FASTA="${OUT_DIR}/${RUN_PREFIX}_Step7.fa"
[[ -s "${FINAL_FASTA}" ]] || { echo "ERROR: expected MITE-Hunter Step7 FASTA not found: ${FINAL_FASTA}" >&2; exit 1; }

echo "MITE-Hunter candidate sequences: $(grep -c '^>' "${FINAL_FASTA}" || true)"
echo "MITE-Hunter finished: $(date)"
