#!/bin/bash
# Software version: RepeatModeler2 v2.0.5
#SBATCH --job-name=01_repeatmodeler2
#SBATCH --cpus-per-task=36
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

: "${TE_CONFIG:?Export TE_CONFIG=/path/to/te_annotation.env before submission}"
source "${TE_CONFIG}"

: "${GENOME_ID:?GENOME_ID is required}"
: "${GENOME_FASTA:?GENOME_FASTA is required}"
: "${WORK_DIR:?WORK_DIR is required}"

THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-36}}"
OUT_DIR="${WORK_DIR}/01_repeatmodeler2"

if [[ ! -s "${GENOME_FASTA}" ]]; then
    echo "ERROR: genome FASTA is missing or empty: ${GENOME_FASTA}" >&2
    exit 1
fi

command -v BuildDatabase >/dev/null 2>&1 || { echo "ERROR: BuildDatabase not found in PATH" >&2; exit 1; }
command -v RepeatModeler >/dev/null 2>&1 || { echo "ERROR: RepeatModeler not found in PATH" >&2; exit 1; }

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

DB_NAME="${GENOME_ID}_RM2"

echo "Genome ID : ${GENOME_ID}"
echo "Genome    : ${GENOME_FASTA}"
echo "Database  : ${DB_NAME}"
echo "Threads   : ${THREADS}"
echo "Start     : $(date)"

BuildDatabase -name "${DB_NAME}" "${GENOME_FASTA}"

RepeatModeler \
    -database "${DB_NAME}" \
    -threads "${THREADS}" \
    -LTRStruct

echo "RepeatModeler2 finished: $(date)"
echo "Locate the resulting consensi.fa.classified and set RM2_CLASSIFIED_FASTA in TE_CONFIG."
