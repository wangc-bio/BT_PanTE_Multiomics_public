#!/bin/bash
# Public-release note: source script from the executed analysis. Replace the original local absolute paths and executable locations with paths in your environment before running.
#SBATCH -J ATAC_mapQC
#SBATCH -p pms
#SBATCH -c 2
#SBATCH --mem=6G
#SBATCH -t 03:00:00
#SBATCH --array=1-27%27
#SBATCH -o logs/ATAC_mapQC_%A_%a.out
#SBATCH -e logs/ATAC_mapQC_%A_%a.err

set -eo pipefail

BASE=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/PanTE_v52
ATAC=${BASE}/07_TE_ATAC

META=${ATAC}/00_QC/ATAC_BSB_TC_BT_metadata.tsv
SAMTOOLS="${SAMTOOLS:-samtools}"
OUTDIR=${ATAC}/00_QC/mapq_per_sample

ROW=$((SLURM_ARRAY_TASK_ID + 1))

SAMPLE=$(awk -v n=${ROW} 'NR==n {print $1}' ${META} | tr -d '\r')
GROUP=$(awk -v n=${ROW} 'NR==n {print $2}' ${META} | tr -d '\r')
TISSUE=$(awk -v n=${ROW} 'NR==n {print $3}' ${META} | tr -d '\r')
BAM=$(awk -v n=${ROW} 'NR==n {print $4}' ${META} | tr -d '\r')

if [[ -z "${SAMPLE}" || -z "${GROUP}" || -z "${TISSUE}" || -z "${BAM}" ]]; then
    echo "ERROR: failed to parse metadata row ${ROW}" >&2
    exit 1
fi

if [[ ! -s "${BAM}" ]]; then
    echo "ERROR: missing BAM: ${BAM}" >&2
    exit 1
fi

OUT=${OUTDIR}/${SAMPLE}_${TISSUE}.mapq_qc.tsv

echo "============================================"
echo "Task: ${SLURM_ARRAY_TASK_ID}"
echo "Sample: ${SAMPLE}"
echo "Group: ${GROUP}"
echo "Tissue: ${TISSUE}"
echo "BAM: ${BAM}"
echo "============================================"

# ------------------------------------------------------------
# One BAM scan only.
#
# -f 66:
#   0x2  properly paired
#   0x40 first mate
#
# Thus each paired-end fragment contributes exactly one record.
#
# -F 3852 excludes:
#   0x4   read unmapped
#   0x8   mate unmapped
#   0x100 secondary
#   0x200 QC fail
#   0x400 duplicate
#   0x800 supplementary
#
# MAPQ is SAM field 5.
# Q0/Q10/Q20/Q30 are all counted in the SAME scan.
# ------------------------------------------------------------

${SAMTOOLS} view \
    -@ ${SLURM_CPUS_PER_TASK} \
    -f 66 \
    -F 3852 \
    "${BAM}" \
| awk \
    -v sample="${SAMPLE}" \
    -v group="${GROUP}" \
    -v tissue="${TISSUE}" '
BEGIN {
    q0=0;
    q10=0;
    q20=0;
    q30=0;
}
{
    q0++;

    mq=$5+0;

    if (mq >= 10) q10++;
    if (mq >= 20) q20++;
    if (mq >= 30) q30++;
}
END {
    if (q0 > 0) {
        f10=q10/q0;
        f20=q20/q0;
        f30=q30/q0;
    } else {
        f10=0;
        f20=0;
        f30=0;
    }

    print "Sample\tGroup\tTissue\tProper_clean_fragments_Q0\tProper_clean_fragments_Q10\tProper_clean_fragments_Q20\tProper_clean_fragments_Q30\tQ10_fraction\tQ20_fraction\tQ30_fraction";

    printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\n",
        sample,
        group,
        tissue,
        q0,
        q10,
        q20,
        q30,
        f10,
        f20,
        f30;
}' > "${OUT}"

echo ""
echo "Result:"
cat "${OUT}"

echo ""
echo "DONE: ${SAMPLE} ${GROUP} ${TISSUE}"