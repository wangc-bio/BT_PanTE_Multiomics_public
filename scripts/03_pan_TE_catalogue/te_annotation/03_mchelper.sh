#!/bin/bash
# Software version: MCHelper v1.7.0
#SBATCH --job-name=03_mchelper
#SBATCH --cpus-per-task=24
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

: "${TE_CONFIG:?Export TE_CONFIG=/path/to/te_annotation.env before submission}"
source "${TE_CONFIG}"

: "${GENOME_FASTA:?GENOME_FASTA is required}"
: "${WORK_DIR:?WORK_DIR is required}"
: "${MCHELPER_PY:?MCHELPER_PY is required}"
: "${MCHELPER_INPUT_FASTA:?MCHELPER_INPUT_FASTA is required}"
: "${BUSCO_HMM:?BUSCO_HMM is required}"

THREADS="${MCHELPER_THREADS:-24}"
OUT_DIR="${MCHELPER_OUTPUT_DIR:-${WORK_DIR}/03_mchelper}"

for F in "${GENOME_FASTA}" "${MCHELPER_INPUT_FASTA}" "${MCHELPER_PY}" "${BUSCO_HMM}"; do
    [[ -s "${F}" ]] || { echo "ERROR: required input is missing or empty: ${F}" >&2; exit 1; }
done

mkdir -p "${OUT_DIR}"

python3 "${MCHELPER_PY}" \
    -r A \
    -t "${THREADS}" \
    -l "${MCHELPER_INPUT_FASTA}" \
    -g "${GENOME_FASTA}" \
    -o "${OUT_DIR}" \
    --input_type fasta \
    -b "${BUSCO_HMM}" \
    -a F

echo "MCHelper finished: $(date)"
