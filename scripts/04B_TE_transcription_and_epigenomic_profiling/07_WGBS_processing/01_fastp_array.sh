#!/bin/bash
# Public-release note: original executed source; replace local absolute paths/environment activation with local configuration before running.
#SBATCH --job-name=01_fastp_WGBS
#SBATCH --partition=pms
#SBATCH -N 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --array=1-36%6
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

# ====== 路径配置 ======
WORKDIR=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/WGBS/00_Data
FASTP=${FASTP:-fastp}
META=${WORKDIR}/wgbs_metadata.tsv   # 上一步生成的metadata

THREADS=${SLURM_CPUS_PER_TASK}

# ====== 输出目录 ======
OUTDIR=${WORKDIR}/01_fastp
CLEAN=${OUTDIR}/clean_fq
REPORT=${OUTDIR}/reports
LOGDIR=${OUTDIR}/logs

mkdir -p "${CLEAN}" "${REPORT}" "${LOGDIR}"

# ====== 读取第 N 个样本（跳过表头） ======
line=$(awk -v id="${SLURM_ARRAY_TASK_ID}" 'BEGIN{FS=OFS="\t"} NR==id+1{print $0}' "${META}" || true)
if [[ -z "${line}" ]]; then
  echo "[ERROR] Cannot find line for array id=${SLURM_ARRAY_TASK_ID} in ${META}" >&2
  exit 1
fi

# columns:
# sample_id individual_id species tissue tissue_code batch R1 R2
IFS=$'\t' read -r sample_id individual_id species tissue tissue_code batch R1 R2 <<< "${line}"

echo "==== WGBS fastp ===="
echo "sample_id     : ${sample_id}"
echo "individual_id : ${individual_id}"
echo "species       : ${species}"
echo "tissue        : ${tissue}"
echo "batch         : ${batch}"
echo "R1            : ${R1}"
echo "R2            : ${R2}"
echo "threads       : ${THREADS}"
echo "workdir       : ${WORKDIR}"
echo "===================="

# ====== 检查输入 ======
if [[ ! -s "${R1}" ]]; then echo "[ERROR] R1 not found: ${R1}" >&2; exit 1; fi
if [[ ! -s "${R2}" ]]; then echo "[ERROR] R2 not found: ${R2}" >&2; exit 1; fi

# ====== 输出文件名 ======
out1=${CLEAN}/${sample_id}.clean_R1.fq.gz
out2=${CLEAN}/${sample_id}.clean_R2.fq.gz
html=${REPORT}/${sample_id}.fastp.html
json=${REPORT}/${sample_id}.fastp.json
log=${LOGDIR}/${sample_id}.fastp.log

# ====== fastp 运行 ======
# 说明：
# --detect_adapter_for_pe：PE自动识别接头（WGBS常用）
# --trim_poly_g/poly_x：对NextSeq/NovaSeq等polyG/polyX做修剪（一般开着更稳）
# 质量参数你后续可按FastQC/fastp报告微调
"${FASTP}" \
  -i "${R1}" -I "${R2}" \
  -o "${out1}" -O "${out2}" \
  -w "${THREADS}" \
  --detect_adapter_for_pe \
  --trim_poly_g --trim_poly_x \
  -q 20 -u 30 -n 5 -l 30 \
  --compression 6 \
  --html "${html}" --json "${json}" \
  > "${log}" 2>&1

echo "[OK] fastp done: ${sample_id}"
echo "  clean: ${out1}"
echo "  clean: ${out2}"
echo "  report: ${html}"
echo "  report: ${json}"
echo "  log: ${log}"