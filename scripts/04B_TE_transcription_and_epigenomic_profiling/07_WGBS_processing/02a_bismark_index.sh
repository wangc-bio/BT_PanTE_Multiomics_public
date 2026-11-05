#!/bin/bash
# Public-release note: original executed source; replace local absolute paths/environment activation with local configuration before running.
#SBATCH --job-name=02a_bismark_index
#SBATCH --partition=pms
#SBATCH -N 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=80G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

# =========================
# 0. 激活 Bismark 环境
# =========================
# Ensure Bismark and Bowtie2 are available in PATH before submitting this job.
# If your cluster requires environment activation, add it in your local job wrapper.

echo "[INFO] conda env: ${CONDA_DEFAULT_ENV}"
echo "[INFO] bismark_genome_preparation: $(which bismark_genome_preparation)"
echo "[INFO] bowtie2: $(which bowtie2)"

BISMARK_GENOME_PREPARATION=bismark_genome_preparation

# =========================
# 1. 路径配置
# =========================
WORKDIR=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/WGBS/00_Data
INDEX_ROOT=${WORKDIR}/02_bismark_index
LOGDIR=${INDEX_ROOT}/logs

mkdir -p "${INDEX_ROOT}" "${LOGDIR}"

BSB_REF=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/ref_rename/Mamb/index/BSB.genome.fa
TC_REF=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/ref_rename/Calb/index/TC.genome.fa
BT_REF=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/Pb_asm/451/purge_hap/haphic/one_step/04.build/anno/EVM_v2/final_gff2/01_rename/index/BT.genome.rename.fa
BTB_REF=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/Pb_asm/477/purge_hap/haphic/hap1_hap2/anno/EVM_v3/merge_final/01_rename/index/BTB.final.rename.genome.fa

# =========================
# 2. 函数：构建索引
# =========================
prepare_one_ref () {
    species=$1
    ref=$2

    outdir=${INDEX_ROOT}/${species}
    mkdir -p "${outdir}"

    echo "======================================"
    echo "[INFO] species: ${species}"
    echo "[INFO] ref    : ${ref}"
    echo "[INFO] outdir : ${outdir}"
    echo "======================================"

    if [[ ! -s "${ref}" ]]; then
        echo "[ERROR] Reference fasta not found: ${ref}" >&2
        exit 1
    fi

    ln -sf "${ref}" "${outdir}/${species}.genome.fa"

    if [[ -d "${outdir}/Bisulfite_Genome" ]]; then
        echo "[SKIP] Bismark index already exists: ${outdir}/Bisulfite_Genome"
    else
        echo "[RUN] bismark_genome_preparation for ${species}"

        "${BISMARK_GENOME_PREPARATION}" \
            --bowtie2 \
            --parallel 12 \
            "${outdir}"
    fi

    echo "[DONE] ${species}"
}

# =========================
# 3. 分别构建索引
# =========================
prepare_one_ref BSB "${BSB_REF}"
prepare_one_ref TC  "${TC_REF}"
prepare_one_ref BT  "${BT_REF}"
prepare_one_ref BTB "${BTB_REF}"

echo "[ALL DONE] Bismark genome preparation finished."