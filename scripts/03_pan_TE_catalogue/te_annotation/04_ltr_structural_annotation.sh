#!/bin/bash
# Software versions: LTR_FINDER_parallel v1.2 (LTR_FINDER v1.0.7), GenomeTools v1.6.2, LTR_retriever v3.0.1
#SBATCH --job-name=04_ltr_structural_annotation
#SBATCH --cpus-per-task=24
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

: "${TE_CONFIG:?Export TE_CONFIG=/path/to/te_annotation.env before submission}"
source "${TE_CONFIG}"

: "${GENOME_FASTA:?GENOME_FASTA is required}"
: "${WORK_DIR:?WORK_DIR is required}"
: "${LTR_FINDER_PARALLEL:?LTR_FINDER_PARALLEL is required}"
: "${GENOMETOOLS_GT:?GENOMETOOLS_GT is required}"
: "${LTR_RETRIEVER:?LTR_RETRIEVER is required}"

THREADS="${SLURM_CPUS_PER_TASK:-24}"
OUT_DIR="${WORK_DIR}/04_ltr_structure"
GENOME_LINK="${OUT_DIR}/genome.fa"

for F in "${GENOME_FASTA}" "${LTR_FINDER_PARALLEL}" "${GENOMETOOLS_GT}" "${LTR_RETRIEVER}"; do
    [[ -s "${F}" ]] || { echo "ERROR: required file or executable is missing: ${F}" >&2; exit 1; }
done

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"
ln -sf "${GENOME_FASTA}" "${GENOME_LINK}"

"${LTR_FINDER_PARALLEL}" \
    -seq "${GENOME_LINK}" \
    -threads "${THREADS}" \
    -harvest_out \
    -size 1000000 \
    -time 300

"${GENOMETOOLS_GT}" suffixerator \
    -db "${GENOME_LINK}" \
    -indexname "${GENOME_LINK}" \
    -tis -suf -lcp -des -ssp -sds -dna

"${GENOMETOOLS_GT}" ltrharvest \
    -index "${GENOME_LINK}" \
    -minlenltr 100 \
    -maxlenltr 7000 \
    -mintsd 4 \
    -maxtsd 6 \
    -motif TGCA \
    -motifmis 1 \
    -similar 85 \
    -vic 10 \
    -seed 20 \
    -seqids yes \
    > "${GENOME_LINK}.harvest.scn"

FINDER_RESULT="${GENOME_LINK}.finder.combine.scn"
HARVEST_RESULT="${GENOME_LINK}.harvest.scn"
RAW_LTR="${GENOME_LINK}.rawLTR.scn"

[[ -s "${FINDER_RESULT}" ]] || { echo "ERROR: LTR_FINDER_parallel output not found: ${FINDER_RESULT}" >&2; exit 1; }
[[ -s "${HARVEST_RESULT}" ]] || { echo "ERROR: LTRharvest output not found: ${HARVEST_RESULT}" >&2; exit 1; }

cat "${HARVEST_RESULT}" "${FINDER_RESULT}" > "${RAW_LTR}"

"${LTR_RETRIEVER}" \
    -genome "${GENOME_LINK}" \
    -inharvest "${RAW_LTR}" \
    -threads "${THREADS}"

echo "LTR structural annotation finished: $(date)"
