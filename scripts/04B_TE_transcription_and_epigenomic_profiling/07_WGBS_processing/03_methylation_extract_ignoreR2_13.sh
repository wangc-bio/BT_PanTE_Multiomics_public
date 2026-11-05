#!/bin/bash
# Public-release note: original executed source; replace local absolute paths/environment activation with local configuration before running.
#SBATCH -J WGBS_extract_R2trim13
#SBATCH -p pms
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=80G
#SBATCH --array=1-36%4
#SBATCH -t 12:00:00
#SBATCH -o logs/WGBS_R2trim13.%A_%a.out
#SBATCH -e logs/WGBS_R2trim13.%A_%a.err

set -eo pipefail

# Ensure Bismark methylation-extractor is available in PATH before submitting this job.
# If your cluster requires environment activation, add it in your local job wrapper.

WORKDIR=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/WGBS/00_Data
META=${WORKDIR}/wgbs_metadata.tsv

INDEX_ROOT=${WORKDIR}/02_bismark_index

# ???? v2 ????? name-sorted BAM
NAME_BAM_DIR=${WORKDIR}/03_methylation_extract_v2/name_sorted_bam

OUTROOT=${WORKDIR}/03_methylation_extract_v3_ignoreR2_13

mkdir -p "${OUTROOT}"

EXTRACT_CORES=4

# ============================================================
# metadata
# ============================================================

line=$(awk -v id="${SLURM_ARRAY_TASK_ID}" \
    'BEGIN{FS=OFS="\t"} NR==id+1{print $0}' \
    "${META}" || true)

if [[ -z "${line}" ]]; then
    echo "[ERROR] Cannot find task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

IFS=$'\t' read -r \
sample_id individual_id species tissue tissue_code batch raw_R1 raw_R2 \
<<< "${line}"

echo "========================================"
echo "sample      : ${sample_id}"
echo "individual  : ${individual_id}"
echo "species     : ${species}"
echo "tissue      : ${tissue}"
echo "batch       : ${batch}"
echo "ignore_r2   : 13"
echo "========================================"

# ============================================================
# reference
# ============================================================

case "${species}" in
    BSB)
        GENOME_DIR=${INDEX_ROOT}/BSB
        ;;
    TC)
        GENOME_DIR=${INDEX_ROOT}/TC
        ;;
    BT)
        GENOME_DIR=${INDEX_ROOT}/BT
        ;;
    BTB)
        GENOME_DIR=${INDEX_ROOT}/BTB
        ;;
    *)
        echo "[ERROR] Unknown species: ${species}" >&2
        exit 1
        ;;
esac

if [[ ! -d "${GENOME_DIR}/Bisulfite_Genome" ]]; then
    echo "[ERROR] Missing Bismark genome: ${GENOME_DIR}" >&2
    exit 1
fi

# ============================================================
# existing name-sorted BAM
# ============================================================

NAME_BAM=${NAME_BAM_DIR}/${sample_id}.deduplicated.name_sorted.bam

if [[ ! -s "${NAME_BAM}" ]]; then
    echo "[ERROR] Missing name-sorted BAM:"
    echo "${NAME_BAM}"
    exit 1
fi

# ============================================================
# output
# ============================================================

SAMPLE_OUT=${OUTROOT}/${sample_id}

mkdir -p "${SAMPLE_OUT}"

echo "[INFO] BAM        : ${NAME_BAM}"
echo "[INFO] GENOME_DIR : ${GENOME_DIR}"
echo "[INFO] OUT        : ${SAMPLE_OUT}"
echo "[INFO] START      : $(date)"

# ============================================================
# methylation extraction
# ============================================================

bismark_methylation_extractor \
    --paired-end \
    --comprehensive \
    --merge_non_CpG \
    --bedGraph \
    --counts \
    --report \
    --gzip \
    --ignore_r2 13 \
    --multicore "${EXTRACT_CORES}" \
    --buffer_size 10G \
    --genome_folder "${GENOME_DIR}" \
    --output "${SAMPLE_OUT}" \
    "${NAME_BAM}"

echo
echo "===== output check ====="

COV=$(find "${SAMPLE_OUT}" \
    -maxdepth 1 \
    -type f \
    -name "*.bismark.cov.gz" \
    | head -1)

if [[ -z "${COV}" || ! -s "${COV}" ]]; then
    echo "[ERROR] Missing bismark cov output" >&2
    exit 1
fi

echo "[OK] COV: ${COV}"
echo "[OK] FINISHED: $(date)"