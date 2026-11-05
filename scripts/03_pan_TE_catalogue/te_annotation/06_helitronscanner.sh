#!/bin/bash
# Software version: HelitronScanner v1.1
#SBATCH --job-name=06_helitronscanner
#SBATCH --cpus-per-task=36
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

: "${TE_CONFIG:?Export TE_CONFIG=/path/to/te_annotation.env before submission}"
source "${TE_CONFIG}"

: "${GENOME_FASTA:?GENOME_FASTA is required}"
: "${WORK_DIR:?WORK_DIR is required}"
: "${HELITRONSCANNER_TRAINING_DIR:?HELITRONSCANNER_TRAINING_DIR is required}"

THREADS="${SLURM_CPUS_PER_TASK:-36}"
OUT_DIR="${WORK_DIR}/06_helitronscanner"
HEAD="${HELITRONSCANNER_TRAINING_DIR}/head.lcvs"
TAIL="${HELITRONSCANNER_TRAINING_DIR}/tail.lcvs"

for F in "${GENOME_FASTA}" "${HEAD}" "${TAIL}"; do
    [[ -s "${F}" ]] || { echo "ERROR: required input is missing or empty: ${F}" >&2; exit 1; }
done
command -v HelitronScanner >/dev/null 2>&1 || { echo "ERROR: HelitronScanner not found in PATH" >&2; exit 1; }

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# Forward strand
HelitronScanner scanHead \
    -lf "${HEAD}" -g "${GENOME_FASTA}" \
    -bs 1000000 -ol 50 -th 1 -tl "${THREADS}" \
    -o FWD.head.out

HelitronScanner scanTail \
    -lf "${TAIL}" -g "${GENOME_FASTA}" \
    -bs 1000000 -ol 50 -th 1 -tl "${THREADS}" \
    -o FWD.tail.out

HelitronScanner pairends \
    -hs FWD.head.out -ts FWD.tail.out \
    -ht 5 -tt 10 -hlr 200:20000 \
    -o FWD.paired.out

HelitronScanner draw \
    -p FWD.paired.out -g "${GENOME_FASTA}" \
    -o Helitron_FWD --pure

# Reverse-complement strand
HelitronScanner scanHead \
    -lf "${HEAD}" -g "${GENOME_FASTA}" \
    -bs 1000000 -ol 50 -th 1 -tl "${THREADS}" --rc \
    -o RC.head.out

HelitronScanner scanTail \
    -lf "${TAIL}" -g "${GENOME_FASTA}" \
    -bs 1000000 -ol 50 -th 1 -tl "${THREADS}" --rc \
    -o RC.tail.out

HelitronScanner pairends \
    -hs RC.head.out -ts RC.tail.out \
    -ht 5 -tt 10 -hlr 200:20000 --rc \
    -o RC.paired.out

HelitronScanner draw \
    -p RC.paired.out -g "${GENOME_FASTA}" \
    -o Helitron_RC --pure

ALL_FASTA="${OUT_DIR}/HelitronScanner.all.fa"
FWD_FASTA="${OUT_DIR}/Helitron_FWD.hel.fa"
RC_FASTA="${OUT_DIR}/Helitron_RC.hel.fa"

: > "${ALL_FASTA}"

if [[ -s "${FWD_FASTA}" ]]; then
    awk '/^>/ {print ">FWD_" substr($0,2); next} {print}' "${FWD_FASTA}" >> "${ALL_FASTA}"
fi
if [[ -s "${RC_FASTA}" ]]; then
    awk '/^>/ {print ">RC_" substr($0,2); next} {print}' "${RC_FASTA}" >> "${ALL_FASTA}"
fi

echo "Forward candidates : $(grep -c '^>' "${FWD_FASTA}" 2>/dev/null || true)"
echo "Reverse candidates : $(grep -c '^>' "${RC_FASTA}" 2>/dev/null || true)"
echo "Combined candidates: $(grep -c '^>' "${ALL_FASTA}" 2>/dev/null || true)"
echo "HelitronScanner finished: $(date)"
