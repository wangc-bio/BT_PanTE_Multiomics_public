#!/bin/bash
# Public-release note: source script from the executed analysis. Replace the original local absolute paths and executable locations with paths in your environment before running.
#SBATCH -J ATAC_mapQC
#SBATCH -p pms
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH -t 06:00:00
#SBATCH -o logs/ATAC_mapQC_%j.out
#SBATCH -e logs/ATAC_mapQC_%j.err

set -eo pipefail

BASE=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/PanTE_v52
ATAC=${BASE}/07_TE_ATAC

META=${ATAC}/00_QC/ATAC_BSB_TC_BT_metadata.tsv

SAMTOOLS="${SAMTOOLS:-samtools}"
OUT=${ATAC}/00_QC/ATAC_mapping_quality_summary.tsv


echo -e \
"Sample\tGroup\tTissue\tProper_clean_fragments_Q0\tProper_clean_fragments_Q10\tProper_clean_fragments_Q20\tProper_clean_fragments_Q30\tQ10_fraction\tQ20_fraction\tQ30_fraction" \
> ${OUT}


tail -n +2 ${META} | \
while IFS=$'\t' read SAMPLE GROUP TISSUE BAM BAM_EXISTS INDEX_EXISTS
do

    echo "[INFO] ${SAMPLE} ${GROUP} ${TISSUE}"

    # --------------------------------------------------------
    # Count one alignment per paired-end fragment:
    #
    # -f 66 =
    #   0x2  properly paired
    #   0x40 first mate
    #
    # -F 3852 excludes:
    #   0x4   read unmapped
    #   0x8   mate unmapped
    #   0x100 secondary
    #   0x200 QC fail
    #   0x400 duplicate
    #   0x800 supplementary
    #
    # Thus one read1 record represents one clean proper fragment.
    # --------------------------------------------------------

    Q0=$(
        ${SAMTOOLS} view \
        -c \
        -f 66 \
        -F 3852 \
        "${BAM}"
    )

    Q10=$(
        ${SAMTOOLS} view \
        -c \
        -f 66 \
        -F 3852 \
        -q 10 \
        "${BAM}"
    )

    Q20=$(
        ${SAMTOOLS} view \
        -c \
        -f 66 \
        -F 3852 \
        -q 20 \
        "${BAM}"
    )

    Q30=$(
        ${SAMTOOLS} view \
        -c \
        -f 66 \
        -F 3852 \
        -q 30 \
        "${BAM}"
    )


    F10=$(
        awk -v a=${Q10} -v b=${Q0} \
        'BEGIN{if(b>0) printf "%.6f",a/b; else print "NA"}'
    )

    F20=$(
        awk -v a=${Q20} -v b=${Q0} \
        'BEGIN{if(b>0) printf "%.6f",a/b; else print "NA"}'
    )

    F30=$(
        awk -v a=${Q30} -v b=${Q0} \
        'BEGIN{if(b>0) printf "%.6f",a/b; else print "NA"}'
    )


    echo -e \
"${SAMPLE}\t${GROUP}\t${TISSUE}\t${Q0}\t${Q10}\t${Q20}\t${Q30}\t${F10}\t${F20}\t${F30}" \
>> ${OUT}

done


echo ""
echo "=========================================="
echo "ATAC mapping-quality QC completed"
echo "=========================================="

column -t ${OUT}