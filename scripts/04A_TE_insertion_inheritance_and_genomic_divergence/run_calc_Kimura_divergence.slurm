#!/bin/bash
#SBATCH --job-name=run_calc_Kimura_divergence
#SBATCH --partition=pms
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=72:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/Kimura_%A_%a.out
#SBATCH --error=logs/Kimura_%A_%a.err

set -euo pipefail

# Set these variables in the submission environment:
# KIMURA_BASE: RepeatMasker_unified_v5_2_Kimura_align directory
# REPEATMASKER_UTIL: RepeatMasker util directory containing calcDivergenceFromAlign.pl
BASE=${KIMURA_BASE:?Set KIMURA_BASE to the RepeatMasker Kimura analysis directory}
RMUTIL=${REPEATMASKER_UTIL:?Set REPEATMASKER_UTIL to the RepeatMasker util directory}

NAMES=(BSB TC BTM BTP)

ALIGNS=(
    "${BASE}/results/BSB/BSB.genome.chr.fa.align"
    "${BASE}/results/TC/TC.genome.chr.fa.align"
    "${BASE}/results/BTM/BT.genome.subM.fa.align"
    "${BASE}/results/BTP/BT.genome.subP.fa.align"
)

NAME=${NAMES[$SLURM_ARRAY_TASK_ID]}
ALIGN=${ALIGNS[$SLURM_ARRAY_TASK_ID]}

OUTDIR="${BASE}/Kimura_divergence/${NAME}"
DIVSUM="${OUTDIR}/${NAME}.divsum"
CALC="${RMUTIL}/calcDivergenceFromAlign.pl"

mkdir -p "${OUTDIR}"

echo "=========================================="
echo "Sample      : ${NAME}"
echo "Input align : ${ALIGN}"
echo "Utility     : ${CALC}"
echo "Output      : ${DIVSUM}"
echo "Start       : $(date)"
echo "Host        : $(hostname)"
echo "=========================================="

if [[ ! -s "${ALIGN}" ]]; then
    echo "[ERROR] Missing or empty align file: ${ALIGN}" >&2
    exit 1
fi

if [[ ! -x "${CALC}" ]]; then
    echo "[ERROR] calcDivergenceFromAlign.pl not executable: ${CALC}" >&2
    exit 1
fi

perl "${CALC}" -s "${DIVSUM}" "${ALIGN}"

if [[ ! -s "${DIVSUM}" ]]; then
    echo "[ERROR] divsum was not generated: ${DIVSUM}" >&2
    exit 1
fi

echo
echo "=========================================="
echo "Completed   : ${NAME}"
echo "End         : $(date)"
echo "divsum size : $(du -h "${DIVSUM}" | cut -f1)"
echo "=========================================="
