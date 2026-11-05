#!/bin/bash
# Public-release note: original executed source; replace local absolute paths/environment activation with your local configuration before running.
#SBATCH -J TEcount_multi
#SBATCH -p pms
#SBATCH -o logs/TEcount_multi.%A_%a.out
#SBATCH -e logs/TEcount_multi.%A_%a.err
#SBATCH --array=1-57%10
#SBATCH -c 1
#SBATCH --mem=80G
#SBATCH -t 12:00:00

set -eo pipefail

source "${CONDA_SH:?Set CONDA_SH to your conda.sh path}"

conda activate \
${CONDA_ENV:-TEtranscripts}

RNA=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/PanTE_v52/05_TE_RNA
JOBFILE=${RNA}/00_metadata/TEcount_jobs.tsv

LINE=$(awk -F'\t' -v id="${SLURM_ARRAY_TASK_ID}" \
    'NR>1 && $1==id {print; exit}' "${JOBFILE}")

if [[ -z "${LINE}" ]]; then
    echo "[ERROR] Cannot find task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

IFS=$'\t' read -r \
TASK_ID SAMPLE GROUP INDIVIDUAL TISSUE BAM GENE_GTF TE_GTF \
<<< "${LINE}"

OUT=${RNA}/02_TEcount_multi/${GROUP}/${TISSUE}/${SAMPLE}

mkdir -p "${OUT}"

echo "===== Task information ====="
echo "Task_ID     : ${TASK_ID}"
echo "Sample      : ${SAMPLE}"
echo "Group       : ${GROUP}"
echo "Individual  : ${INDIVIDUAL}"
echo "Tissue      : ${TISSUE}"
echo "BAM         : ${BAM}"
echo "Gene_GTF    : ${GENE_GTF}"
echo "TE_GTF      : ${TE_GTF}"
echo "Output      : ${OUT}"
echo "Start       : $(date)"
echo

for F in "${BAM}" "${GENE_GTF}" "${TE_GTF}"
do
    if [[ ! -s "${F}" ]]; then
        echo "[ERROR] Missing/empty input: ${F}" >&2
        exit 1
    fi
done

TEcount \
    -b "${BAM}" \
    --GTF "${GENE_GTF}" \
    --TE "${TE_GTF}" \
    --format BAM \
    --sortByPos \
    --stranded no \
    --mode multi \
    --project "${SAMPLE}" \
    --outdir "${OUT}" \
    --verbose 2

CNT=${OUT}/${SAMPLE}.cntTable

if [[ ! -s "${CNT}" ]]; then
    echo "[ERROR] Missing count table: ${CNT}" >&2
    exit 1
fi

echo
echo "===== PanTE summary ====="

awk -F'\t' '
$1 ~ /^PanTE_R/ {
    n++;
    if ($2 > 0) nz++;
    sum += $2
}
END {
    print "PanTE_rows\t" n;
    print "Nonzero_PanTE\t" nz;
    print "Total_PanTE_counts\t" sum;
}' "${CNT}"

echo
echo "Finished: $(date)"