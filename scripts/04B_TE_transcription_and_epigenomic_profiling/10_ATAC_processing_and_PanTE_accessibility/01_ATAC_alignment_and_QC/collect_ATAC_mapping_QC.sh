#!/bin/bash
# Public-release note: source script from the executed analysis. Replace the original local absolute paths and executable locations with paths in your environment before running.
#SBATCH -J ATAC_QCcollect
#SBATCH -p pms
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 00:30:00
#SBATCH -o logs/ATAC_QCcollect_%j.out
#SBATCH -e logs/ATAC_QCcollect_%j.err

set -eo pipefail

BASE=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/PanTE_v52
ATAC=${BASE}/07_TE_ATAC

META=${ATAC}/00_QC/ATAC_BSB_TC_BT_metadata.tsv
INDIR=${ATAC}/00_QC/mapq_per_sample
OUT=${ATAC}/00_QC/ATAC_mapping_quality_summary.tsv

echo -e \
"Sample\tGroup\tTissue\tProper_clean_fragments_Q0\tProper_clean_fragments_Q10\tProper_clean_fragments_Q20\tProper_clean_fragments_Q30\tQ10_fraction\tQ20_fraction\tQ30_fraction" \
> "${OUT}"

N=0

tail -n +2 "${META}" | \
while IFS=$'\t' read SAMPLE GROUP TISSUE BAM BAM_EXISTS INDEX_EXISTS
do
    F=${INDIR}/${SAMPLE}_${TISSUE}.mapq_qc.tsv

    if [[ ! -s "${F}" ]]; then
        echo "ERROR: missing result: ${F}" >&2
        exit 1
    fi

    tail -n +2 "${F}" >> "${OUT}"
done

N=$(($(wc -l < "${OUT}") - 1))

if [[ "${N}" -ne 27 ]]; then
    echo "ERROR: expected 27 records, found ${N}" >&2
    exit 1
fi

echo ""
echo "============================================"
echo "ATAC mapping QC summary completed"
echo "Records: ${N}"
echo "============================================"
echo ""

column -t "${OUT}"