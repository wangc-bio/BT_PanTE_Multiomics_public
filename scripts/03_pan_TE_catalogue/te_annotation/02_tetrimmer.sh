#!/bin/bash
# Software version: TEtrimmer v1.7.2
#SBATCH --job-name=02_tetrimmer
#SBATCH --cpus-per-task=36
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

: "${TE_CONFIG:?Export TE_CONFIG=/path/to/te_annotation.env before submission}"
source "${TE_CONFIG}"

: "${GENOME_FASTA:?GENOME_FASTA is required}"
: "${WORK_DIR:?WORK_DIR is required}"
: "${TETRIMMER_PY:?TETRIMMER_PY is required}"
: "${TETRIMMER_INPUT_FASTA:?TETRIMMER_INPUT_FASTA is required}"
: "${PFAM_DIR:?PFAM_DIR is required}"

THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-36}}"
OUT_DIR="${TETRIMMER_OUTPUT_DIR:-${WORK_DIR}/02_tetrimmer}"
LOG_FILE="${OUT_DIR}/TEtrimmer.log"

for F in "${GENOME_FASTA}" "${TETRIMMER_INPUT_FASTA}" "${TETRIMMER_PY}"; do
    [[ -s "${F}" ]] || { echo "ERROR: required input is missing or empty: ${F}" >&2; exit 1; }
done
[[ -d "${PFAM_DIR}" ]] || { echo "ERROR: Pfam directory not found: ${PFAM_DIR}" >&2; exit 1; }

mkdir -p "${OUT_DIR}"

python3 "${TETRIMMER_PY}" \
    -i "${TETRIMMER_INPUT_FASTA}" \
    -g "${GENOME_FASTA}" \
    -o "${OUT_DIR}" \
    -t "${THREADS}" \
    --classify_unknown \
    --pfam_dir "${PFAM_DIR}" \
    -l "${LOG_FILE}" \
    -ll INFO

echo "TEtrimmer finished: $(date)"
