#!/bin/bash
# Software versions: TEsorter v1.5.1; SeqKit version not recorded in the source execution record.
#SBATCH --job-name=07_tesorter
#SBATCH --cpus-per-task=4
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

: "${TE_CONFIG:?Export TE_CONFIG=/path/to/te_annotation.env before submission}"
source "${TE_CONFIG}"

: "${TE_SORTER_INPUT_FASTA:?TE_SORTER_INPUT_FASTA is required}"
: "${WORK_DIR:?WORK_DIR is required}"

OUT_DIR="${WORK_DIR}/07_tesorter"
UNKNOWN_FASTA="${OUT_DIR}/Unknown.fa"
THREADS=4

[[ -s "${TE_SORTER_INPUT_FASTA}" ]] || { echo "ERROR: TEsorter input FASTA is missing: ${TE_SORTER_INPUT_FASTA}" >&2; exit 1; }
command -v seqkit >/dev/null 2>&1 || { echo "ERROR: seqkit not found in PATH" >&2; exit 1; }
command -v TEsorter >/dev/null 2>&1 || { echo "ERROR: TEsorter not found in PATH" >&2; exit 1; }

mkdir -p "${OUT_DIR}"

seqkit grep -r -p '#Unknown$' "${TE_SORTER_INPUT_FASTA}" > "${UNKNOWN_FASTA}"

N_UNKNOWN=$(grep -c '^>' "${UNKNOWN_FASTA}" || true)
echo "Unknown consensus sequences: ${N_UNKNOWN}"

if [[ "${N_UNKNOWN}" -gt 0 ]]; then
    TEsorter "${UNKNOWN_FASTA}" \
        -db rexdb-metazoa \
        -p "${THREADS}" \
        -pre "${OUT_DIR}/Unknown.rexdb"
else
    echo "No Unknown consensus sequences; TEsorter was not run."
fi

echo "TEsorter finished: $(date)"
