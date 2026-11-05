#!/bin/bash
# Public-release note: original executed source; replace local absolute paths/environment activation with local configuration before running.
#SBATCH --job-name=02b_bismark_align
#SBATCH --partition=pms
#SBATCH -N 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=80G
#SBATCH --array=1-36%4
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

# =========================
# 0. 激活 Bismark 环境
# =========================
# Ensure Bismark, Bowtie2 and samtools are available in PATH before submitting this job.
# If your cluster requires environment activation, add it in your local job wrapper.

echo "[INFO] conda env: ${CONDA_DEFAULT_ENV}"
echo "[INFO] bismark: $(which bismark)"
echo "[INFO] deduplicate_bismark: $(which deduplicate_bismark)"
echo "[INFO] bowtie2: $(which bowtie2)"
echo "[INFO] samtools: $(which samtools || true)"

BISMARK=bismark
DEDUP=deduplicate_bismark
SAMTOOLS=samtools

# =========================
# 1. 路径配置
# =========================
WORKDIR=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/WGBS/00_Data
META=${WORKDIR}/wgbs_metadata.tsv

CLEAN_DIR=${WORKDIR}/01_fastp/clean_fq
INDEX_ROOT=${WORKDIR}/02_bismark_index

OUTROOT=${WORKDIR}/02_bismark_align
ALIGN_DIR=${OUTROOT}/bam
DEDUP_DIR=${OUTROOT}/dedup_bam
REPORT_DIR=${OUTROOT}/reports
LOGDIR=${OUTROOT}/logs

mkdir -p "${ALIGN_DIR}" "${DEDUP_DIR}" "${REPORT_DIR}" "${LOGDIR}"

THREADS=${SLURM_CPUS_PER_TASK}

# =========================
# 2. 读取 metadata
# =========================
line=$(awk -v id="${SLURM_ARRAY_TASK_ID}" 'BEGIN{FS=OFS="\t"} NR==id+1{print $0}' "${META}" || true)

if [[ -z "${line}" ]]; then
    echo "[ERROR] Cannot find sample line for array id=${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

# metadata 列：
# sample_id individual_id species tissue tissue_code batch R1 R2
IFS=$'\t' read -r sample_id individual_id species tissue tissue_code batch raw_R1 raw_R2 <<< "${line}"

echo "======================================"
echo "[INFO] sample_id     : ${sample_id}"
echo "[INFO] individual_id : ${individual_id}"
echo "[INFO] species       : ${species}"
echo "[INFO] tissue        : ${tissue}"
echo "[INFO] batch         : ${batch}"
echo "[INFO] threads       : ${THREADS}"
echo "======================================"

# =========================
# 3. clean reads 路径
# =========================
R1=${CLEAN_DIR}/${sample_id}.clean_R1.fq.gz
R2=${CLEAN_DIR}/${sample_id}.clean_R2.fq.gz

if [[ ! -s "${R1}" ]]; then
    echo "[ERROR] clean R1 not found: ${R1}" >&2
    exit 1
fi

if [[ ! -s "${R2}" ]]; then
    echo "[ERROR] clean R2 not found: ${R2}" >&2
    exit 1
fi

# =========================
# 4. 根据物种选择参考
# =========================
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
    echo "[ERROR] Bismark genome index not found: ${GENOME_DIR}/Bisulfite_Genome" >&2
    exit 1
fi

echo "[INFO] genome dir: ${GENOME_DIR}"
echo "[INFO] R1: ${R1}"
echo "[INFO] R2: ${R2}"

# =========================
# 5. 每个样本独立输出目录
# =========================
SAMPLE_ALIGN_DIR=${ALIGN_DIR}/${sample_id}
SAMPLE_REPORT_DIR=${REPORT_DIR}/${sample_id}

mkdir -p "${SAMPLE_ALIGN_DIR}" "${SAMPLE_REPORT_DIR}"

# =========================
# 6. Bismark 比对
# =========================
# 默认 directional WGBS。
# 如果公司说明是 non-directional，需要加 --non_directional。

echo "[RUN] Bismark alignment: ${sample_id}"

"${BISMARK}" \
    --genome "${GENOME_DIR}" \
    --bowtie2 \
    -p "${THREADS}" \
    --score_min L,0,-0.6 \
    --gzip \
    --basename "${sample_id}" \
    -o "${SAMPLE_ALIGN_DIR}" \
    -1 "${R1}" \
    -2 "${R2}"

echo "[DONE] Bismark alignment: ${sample_id}"

# =========================
# 7. 定位 Bismark 输出 BAM
# =========================
bam=$(find "${SAMPLE_ALIGN_DIR}" -maxdepth 1 -type f -name "${sample_id}*.bam" | head -n 1)

if [[ -z "${bam}" || ! -s "${bam}" ]]; then
    echo "[ERROR] Bismark BAM not found for ${sample_id}" >&2
    exit 1
fi

echo "[INFO] Bismark BAM: ${bam}"

# 保存比对报告
find "${SAMPLE_ALIGN_DIR}" -maxdepth 1 -type f \( -name "*report.txt" -o -name "*splitting_report.txt" -o -name "*M-bias.txt" \) \
    -exec cp {} "${SAMPLE_REPORT_DIR}/" \; || true

# =========================
# 8. 去重复
# =========================
echo "[RUN] deduplicate_bismark: ${sample_id}"

"${DEDUP}" \
    --paired \
    --bam \
    --output_dir "${DEDUP_DIR}" \
    "${bam}"

echo "[DONE] deduplicate_bismark: ${sample_id}"

# =========================
# 9. 排序并建立索引
# =========================
dedup_bam=$(find "${DEDUP_DIR}" -maxdepth 1 -type f -name "${sample_id}*.deduplicated.bam" | head -n 1)

if [[ -z "${dedup_bam}" || ! -s "${dedup_bam}" ]]; then
    echo "[ERROR] deduplicated BAM not found for ${sample_id}" >&2
    exit 1
fi

final_bam=${DEDUP_DIR}/${sample_id}.deduplicated.sorted.bam

echo "[RUN] samtools sort/index: ${sample_id}"

"${SAMTOOLS}" sort \
    -@ "${THREADS}" \
    -o "${final_bam}" \
    "${dedup_bam}"

"${SAMTOOLS}" index \
    -@ "${THREADS}" \
    "${final_bam}"

if [[ "${dedup_bam}" != "${final_bam}" ]]; then
    rm -f "${dedup_bam}"
fi

echo "[OK] final BAM: ${final_bam}"
echo "[OK] final BAI: ${final_bam}.bai"
echo "[ALL DONE] ${sample_id}"