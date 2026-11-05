#!/bin/bash
# Public-release note: original executed source; replace local absolute paths/environment activation with local configuration before running.
#SBATCH -J PanTE_mCG
#SBATCH -p pms
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH --array=1-27
#SBATCH -o logs/PanTE_mCG_%A_%a.out
#SBATCH -e logs/PanTE_mCG_%A_%a.err

set -eo pipefail

BASE=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/PanTE_v52
WGBS=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/WGBS/00_Data/03_methylation_extract_v3_ignoreR2_13
OUT=${BASE}/06_TE_WGBS

META=${OUT}/00_QC/WGBS_BSB_TC_BT_metadata.tsv

PY=${PYTHON:-python}


# ============================================================
# Robust metadata parsing
# awk treats tabs or spaces as whitespace
# ============================================================

ROW=$((SLURM_ARRAY_TASK_ID + 1))

SAMPLE=$(awk -v n=${ROW} 'NR==n {print $1}' ${META} | tr -d '\r')
GROUP=$(awk -v n=${ROW} 'NR==n {print $2}' ${META} | tr -d '\r')
TISSUE=$(awk -v n=${ROW} 'NR==n {print $3}' ${META} | tr -d '\r')


# ============================================================
# Sanity check
# ============================================================

if [[ -z "${SAMPLE}" || -z "${GROUP}" || -z "${TISSUE}" ]]; then
    echo "ERROR: failed to parse metadata row ${ROW}" >&2
    echo "Metadata file: ${META}" >&2
    awk -v n=${ROW} 'NR==n {print}' ${META} >&2
    exit 1
fi


if [[ "${GROUP}" != "BSB" && "${GROUP}" != "TC" && "${GROUP}" != "BT" ]]; then
    echo "ERROR: unexpected GROUP=${GROUP}" >&2
    exit 1
fi


if [[ "${TISSUE}" != "gill" && "${TISSUE}" != "kidney" && "${TISSUE}" != "muscle" ]]; then
    echo "ERROR: unexpected TISSUE=${TISSUE}" >&2
    exit 1
fi


# ============================================================
# Input COV
# ============================================================

COV=${WGBS}/${SAMPLE}_${TISSUE}/${SAMPLE}_${TISSUE}.deduplicated.name_sorted.bismark.cov.gz


# ============================================================
# PanTE intervals
# ============================================================

if [[ "${GROUP}" == "BT" ]]; then

    INTERVAL=${OUT}/01_annotation/BT.PanTE.WGBS_intervals.tsv

else

    INTERVAL=${OUT}/01_annotation/${GROUP}.PanTE.WGBS_intervals.tsv

fi


RESULT=${OUT}/02_quantification/${SAMPLE}_${TISSUE}.PanTE_weighted_mCG.tsv

QC=${OUT}/00_QC/${SAMPLE}_${TISSUE}.PanTE_mCG_QC.tsv


echo "============================================="
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Metadata row: ${ROW}"
echo "Sample: ${SAMPLE}"
echo "Group: ${GROUP}"
echo "Tissue: ${TISSUE}"
echo "COV: ${COV}"
echo "Intervals: ${INTERVAL}"
echo "Result: ${RESULT}"
echo "============================================="


# ============================================================
# Input checks
# ============================================================

if [[ ! -s "${COV}" ]]; then

    echo "ERROR: missing COV file:" >&2
    echo "${COV}" >&2

    echo "" >&2
    echo "Possible files for ${SAMPLE}:" >&2

    find ${WGBS} \
        -type f \
        -name "${SAMPLE}*.bismark.cov.gz" \
        >&2 || true

    exit 1
fi


if [[ ! -s "${INTERVAL}" ]]; then

    echo "ERROR: missing interval file:" >&2
    echo "${INTERVAL}" >&2

    exit 1
fi


# ============================================================
# Quantification
# ============================================================

${PY} \
${OUT}/02_quantification/quantify_PanTE_weighted_mCG.py \
--cov "${COV}" \
--intervals "${INTERVAL}" \
--sample "${SAMPLE}" \
--group "${GROUP}" \
--tissue "${TISSUE}" \
--out "${RESULT}" \
--qc "${QC}"


echo ""
echo "DONE: ${SAMPLE} ${GROUP} ${TISSUE}"