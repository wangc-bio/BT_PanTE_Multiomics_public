#!/bin/bash
# Public-release note: source script from the executed analysis. Replace the original local absolute paths and executable locations with paths in your environment before running.
#SBATCH -J ATAC_PanTE
#SBATCH -p pms
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 08:00:00
#SBATCH --array=1-27%27
#SBATCH -o logs/ATAC_PanTE_%A_%a.out
#SBATCH -e logs/ATAC_PanTE_%A_%a.err

set -eo pipefail


# ============================================================
# PATHS
# ============================================================

BASE=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/PanTE_v52
ATAC=${BASE}/07_TE_ATAC

META=${ATAC}/00_QC/ATAC_BSB_TC_BT_metadata.tsv

SAMTOOLS="${SAMTOOLS:-samtools}"
BEDTOOLS=${BEDTOOLS:-bedtools}
PY=${PYTHON:-python}


# ============================================================
# HARD DEPENDENCY CHECK
# ============================================================

for X in \
"${SAMTOOLS}" \
"${BEDTOOLS}" \
"${PY}"
do
    if [[ ! -x "${X}" ]]; then
        echo "ERROR: executable not found: ${X}" >&2
        exit 1
    fi
done


# ============================================================
# SAMPLE METADATA
# ============================================================

ROW=$((SLURM_ARRAY_TASK_ID + 1))

SAMPLE=$(awk -v n=${ROW} 'NR==n {print $1}' "${META}" | tr -d '\r')
GROUP=$(awk -v n=${ROW} 'NR==n {print $2}' "${META}" | tr -d '\r')
TISSUE=$(awk -v n=${ROW} 'NR==n {print $3}' "${META}" | tr -d '\r')
BAM=$(awk -v n=${ROW} 'NR==n {print $4}' "${META}" | tr -d '\r')


if [[ -z "${SAMPLE}" || -z "${GROUP}" || -z "${TISSUE}" || -z "${BAM}" ]]; then
    echo "ERROR: cannot parse metadata row ${ROW}" >&2
    exit 1
fi


if [[ ! -s "${BAM}" ]]; then
    echo "ERROR: BAM missing: ${BAM}" >&2
    exit 1
fi


if [[ "${GROUP}" == "BT" ]]; then

    SEG=${ATAC}/01_annotation/BT.PanTE.assignable_segments.bed

elif [[ "${GROUP}" == "BSB" ]]; then

    SEG=${ATAC}/01_annotation/BSB.PanTE.assignable_segments.bed

elif [[ "${GROUP}" == "TC" ]]; then

    SEG=${ATAC}/01_annotation/TC.PanTE.assignable_segments.bed

else

    echo "ERROR: unexpected group ${GROUP}" >&2
    exit 1

fi


if [[ ! -s "${SEG}" ]]; then
    echo "ERROR: PanTE segment file missing: ${SEG}" >&2
    exit 1
fi


PREFIX=${ATAC}/02_PanTE_count/${SAMPLE}_${TISSUE}

TOTAL=${PREFIX}.state_total_fragments.tsv
COUNT=${PREFIX}.PanTE_fragment_counts.tsv
QC=${PREFIX}.PanTE_overlap_QC.tsv
NORM=${PREFIX}.PanTE_accessibility.tsv
GENOME=${PREFIX}.main_chrom.genome


echo "============================================"
echo "Sample : ${SAMPLE}"
echo "Group  : ${GROUP}"
echo "Tissue : ${TISSUE}"
echo "BAM    : ${BAM}"
echo "SEG    : ${SEG}"
echo "============================================"


# ============================================================
# BUILD CHROMOSOME ORDER FILE DIRECTLY FROM BAM HEADER
#
# This is important because:
#
# TC:
# chr07 chr17 chr01 ...
#
# BT:
# chr01M chr01P chr02M chr02P ...
#
# are not ordinary lexicographic chromosome orders.
# bedtools -sorted therefore receives -g explicitly.
# ============================================================

${SAMTOOLS} view -H "${BAM}" \
| awk \
    -F'\t' \
    -v OFS='\t' \
    -v group="${GROUP}" '
$1=="@SQ" {

    sn="";
    ln="";

    for(i=2;i<=NF;i++) {

        if($i ~ /^SN:/)
            sn=substr($i,4);

        if($i ~ /^LN:/)
            ln=substr($i,4);
    }


    if(group=="BT") {

        if(sn ~ /^chr[0-9][0-9][MP]$/)
            print sn,ln;

    } else {

        if(sn ~ /^chr[0-9][0-9]$/)
            print sn,ln;
    }
}' > "${GENOME}"


if [[ ! -s "${GENOME}" ]]; then
    echo "ERROR: genome-order file is empty" >&2
    exit 1
fi


echo ""
echo "Main chromosomes:"
wc -l "${GENOME}"


# ============================================================
# ONE BAM SCAN
#
# Filters:
# - proper pair
# - mapped read + mapped mate
# - primary only
# - QC-pass
# - non-duplicate
# - non-supplementary
#
# Existing clean BAMs were already verified:
# all retained alignments have MAPQ >= 30.
#
# TLEN > 0 is used to represent one fragment per pair.
#
# No arbitrary maximum fragment-length threshold is imposed
# at this stage.
# ============================================================

${SAMTOOLS} view \
    -@ ${SLURM_CPUS_PER_TASK} \
    -f 2 \
    -F 3852 \
    "${BAM}" \
| awk \
    -v OFS='\t' \
    -v sample="${SAMPLE}" \
    -v group="${GROUP}" \
    -v totalfile="${TOTAL}" '
BEGIN {

    n=0;

}
{

    chr=$3;
    tlen=$9+0;


    # One fragment per pair
    if(tlen <= 0)
        next;


    if(group=="BT") {

        if(chr ~ /^chr[0-9][0-9]M$/)
            state="BTM";

        else if(chr ~ /^chr[0-9][0-9]P$/)
            state="BTP";

        else
            next;

    } else {

        if(chr !~ /^chr[0-9][0-9]$/)
            next;

        state=group;
    }


    start=$4-1;
    end=start+tlen;


    if(start < 0 || end <= start)
        next;


    n++;
    count[state]++;


    print \
        chr, \
        start, \
        end, \
        sample "_" n, \
        state;
}
END {

    print \
        "Genome_state\tTotal_fragments" \
        > totalfile;


    if(group=="BT") {

        print \
            "BTM\t" (("BTM" in count) ? count["BTM"] : 0) \
            >> totalfile;

        print \
            "BTP\t" (("BTP" in count) ? count["BTP"] : 0) \
            >> totalfile;

    } else {

        print \
            group "\t" ((group in count) ? count[group] : 0) \
            >> totalfile;
    }
}' \
| ${BEDTOOLS} intersect \
    -sorted \
    -g "${GENOME}" \
    -wa \
    -wb \
    -a - \
    -b "${SEG}" \
| awk \
    -F'\t' \
    -v OFS='\t' \
    -v countfile="${COUNT}" \
    -v qcfile="${QC}" '
function flush_fragment( key,a,pct ) {

    if(prev=="")
        return;


    overlap[current_state]++;


    if(nids==1) {

        unique[current_state]++;

        key=current_state SUBSEP one_id;

        pante_count[key]++;

    } else if(nids>1) {

        ambiguous[current_state]++;
    }
}


BEGIN {

    prev="";

}


{

    # A columns:
    # 1 chr
    # 2 start
    # 3 end
    # 4 fragment ID
    # 5 genome state
    #
    # B columns:
    # 6 chr
    # 7 start
    # 8 end
    # 9 PanTE_ID
    # 10 genome state

    frag=$4;
    state=$5;
    pante=$9;
    annotation_state=$10;


    if(state != annotation_state) {

        print \
            "ERROR: state mismatch: " state " vs " annotation_state \
            > "/dev/stderr";

        exit 2;
    }


    if(frag != prev) {

        flush_fragment();

        delete ids;

        nids=0;
        one_id="";
        current_state=state;

        prev=frag;
    }


    if(!(pante in ids)) {

        ids[pante]=1;

        nids++;

        one_id=pante;
    }
}


END {

    flush_fragment();


    print \
        "Genome_state\tPanTE_ID\tFragment_count" \
        > countfile;


    for(key in pante_count) {

        split(key,a,SUBSEP);

        print \
            a[1], \
            a[2], \
            pante_count[key] \
            >> countfile;
    }


    print \
        "Genome_state\tTE_overlapping_fragments\tUnique_assigned_fragments\tAmbiguous_fragments\tUnique_percent" \
        > qcfile;


    for(s in overlap) {

        pct=0;

        if(overlap[s] > 0)
            pct=100.0*unique[s]/overlap[s];


        printf \
            "%s\t%d\t%d\t%d\t%.6f\n", \
            s, \
            overlap[s], \
            unique[s], \
            ambiguous[s], \
            pct \
            >> qcfile;
    }
}'


# ============================================================
# NORMALIZE
# ============================================================

${PY} \
${ATAC}/02_PanTE_count/normalize_ATAC_PanTE_counts.py \
--counts "${COUNT}" \
--bp ${ATAC}/01_annotation/PanTE_assignable_bp.tsv \
--total "${TOTAL}" \
--sample "${SAMPLE}" \
--group "${GROUP}" \
--tissue "${TISSUE}" \
--out "${NORM}"


# ============================================================
# FINAL CHECK
# ============================================================

if [[ ! -s "${NORM}" ]]; then
    echo "ERROR: normalization output missing" >&2
    exit 1
fi


echo ""
echo "============================================"
echo "STATE TOTALS"
echo "============================================"

cat "${TOTAL}"


echo ""
echo "============================================"
echo "OVERLAP QC"
echo "============================================"

cat "${QC}"


echo ""
echo "============================================"
echo "DONE"
echo "============================================"