#!/usr/bin/env bash
# Public, path-portable command templates for read mapping and BAM-level QC.
#
# Recorded software versions:
#   BWA:        0.7.17-r1188
#   minimap2:   2.17-r941
#   STAR:       2.7.10b
#   samblaster: 0.1.26
#   samtools:   1.12
#   HTSlib:     1.12
#   HapHiC:     1.0.5
#
# RNA-seq:
#   STAR 2.7.10b was used to generate BAMs with the Aligned.sortedByCoord.out.bam
#   naming convention. The exact STAR command-line parameters were not present
#   in the supplied run records; no STAR alignment command is fabricated here.
#   Use 02_summarize_star_logs.sh to extract per-sample and aggregate statistics
#   from STAR Log.final.out files.
#
# This public template intentionally contains no cluster-specific paths,
# usernames, private directory names, or sample-specific FASTQ paths.
# Set the required variables in the environment before running.
#
# Usage:
#   MODE=reseq REF=... SAMPLE_ID=... R1=... R2=... OUTDIR=... bash 01_read_mapping_workflow_public.sh
#   MODE=hic REF=... HIC_R1=... HIC_R2=... HAPHIC_FILTER_BAM=... OUTDIR=... bash 01_read_mapping_workflow_public.sh
#   MODE=hifi REF=... HIFI_READS=... OUTDIR=... bash 01_read_mapping_workflow_public.sh
#   MODE=qc BAM=... LABEL=... OUTDIR=... bash 01_read_mapping_workflow_public.sh
#
# Default thread count for resequencing and HiFi: 36.
# Hi-C-specific thread counts reproduce the recorded command settings below.

set -euo pipefail

MODE="${MODE:-help}"
THREADS="${THREADS:-36}"
OUTDIR="${OUTDIR:-mapping_qc}"

require_var() {
    local name="$1"
    if [[ -z "${!name:-}" ]]; then
        echo "ERROR: required environment variable '$name' is not set." >&2
        exit 1
    fi
}

mkdir -p "$OUTDIR"

ensure_bwa_index() {
    require_var REF
    if [[ ! -s "${REF}.bwt" ]]; then
        bwa index "$REF"
    fi
}

run_resequencing() {
    require_var REF
    require_var SAMPLE_ID
    require_var R1
    require_var R2
    ensure_bwa_index

    bwa mem \
        -t "$THREADS" \
        -T 30 \
        -M \
        -R "@RG\tID:${SAMPLE_ID}\tLB:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:BGI" \
        "$REF" "$R1" "$R2" \
        | samtools sort -@ "$THREADS" -o "$OUTDIR/${SAMPLE_ID}.reseq.sorted.bam" -

    samtools index "$OUTDIR/${SAMPLE_ID}.reseq.sorted.bam"
    samtools flagstat "$OUTDIR/${SAMPLE_ID}.reseq.sorted.bam" > "$OUTDIR/${SAMPLE_ID}.reseq.flagstat.txt"
    samtools idxstats "$OUTDIR/${SAMPLE_ID}.reseq.sorted.bam" > "$OUTDIR/${SAMPLE_ID}.reseq.idxstats.txt"
}

run_hic() {
    require_var REF
    require_var HIC_R1
    require_var HIC_R2
    require_var HAPHIC_FILTER_BAM
    ensure_bwa_index

    local hic_threads="${HIC_THREADS:-24}"
    local hic_samtools_threads="${HIC_SAMTOOLS_THREADS:-12}"
    local filter_threads="${HAPHIC_FILTER_THREADS:-14}"
    local filter_samtools_threads="${HIC_FILTER_SAMTOOLS_THREADS:-22}"
    local hic_bam="${HIC_BAM:-$OUTDIR/HiC.bam}"
    local filtered_hic_bam="${FILTERED_HIC_BAM:-$OUTDIR/HiC.filtered.bam}"

    if [[ ! -s "$HIC_R1" || ! -s "$HIC_R2" ]]; then
        echo "ERROR: Hi-C R1/R2 input FASTQ is missing or empty." >&2
        exit 1
    fi
    if ! command -v samblaster >/dev/null 2>&1; then
        echo "ERROR: samblaster was not found in PATH." >&2
        exit 1
    fi
    if ! command -v "$HAPHIC_FILTER_BAM" >/dev/null 2>&1; then
        echo "ERROR: HapHiC filter_bam executable was not found: $HAPHIC_FILTER_BAM" >&2
        exit 1
    fi

    mkdir -p "$(dirname "$hic_bam")" "$(dirname "$filtered_hic_bam")"

    # The three Hi-C libraries belong to the same individual and tissue.
    # The supplied run used the merged paired-end Hi-C data for this mapping.
    bwa mem -5SP -t "$hic_threads" "$REF" "$HIC_R1" "$HIC_R2" \
        | samblaster \
        | samtools view -@ "$hic_samtools_threads" -S -h -b -F 3340 -o "$hic_bam" -

    # HapHiC filter_bam: MAPQ threshold 1 and NM < 3.
    "$HAPHIC_FILTER_BAM" "$hic_bam" 1 --nm 3 --threads "$filter_threads" \
        | samtools view -@ "$filter_samtools_threads" -b -o "$filtered_hic_bam" -

    echo "Hi-C BAM: $hic_bam"
    echo "Hi-C filtered BAM for HapHiC: $filtered_hic_bam"
    echo "Note: flagstat on the -F 3340 BAM is not an unfiltered raw-read mapping rate."
}

run_hifi() {
    require_var REF
    require_var HIFI_READS

    minimap2 -ax map-pb -t "$THREADS" "$REF" "$HIFI_READS" \
        | samtools sort -@ "$THREADS" -o "$OUTDIR/HiFi.sorted.bam" -

    samtools index "$OUTDIR/HiFi.sorted.bam"
    samtools flagstat "$OUTDIR/HiFi.sorted.bam" > "$OUTDIR/HiFi.flagstat.txt"
    samtools idxstats "$OUTDIR/HiFi.sorted.bam" > "$OUTDIR/HiFi.idxstats.txt"
}

run_bam_qc() {
    require_var BAM
    require_var LABEL
    samtools flagstat "$BAM" > "$OUTDIR/${LABEL}.flagstat.txt"
    samtools idxstats "$BAM" > "$OUTDIR/${LABEL}.idxstats.txt"
}

case "$MODE" in
    reseq)
        run_resequencing
        ;;
    hic)
        run_hic
        ;;
    hifi)
        run_hifi
        ;;
    qc)
        run_bam_qc
        ;;
    help)
        cat <<'USAGE'
Select MODE=reseq, hic, hifi, or qc and set the variables listed in the
header. For MODE=hic, also set HAPHIC_FILTER_BAM to the HapHiC filter_bam
executable. The Hi-C mode reproduces the recorded BWA/samblaster/samtools
filtering and HapHiC filter_bam steps.

The RNA-seq STAR invocation is intentionally not reconstructed: the exact
STAR command-line parameters were absent from the supplied records. Use
STAR Log.final.out to report the original input-read mapping rate; do not
interpret 100% flagstat values from mapped-only or filtered BAMs as the
overall raw-read mapping rate.
USAGE
        ;;
    *)
        echo "ERROR: MODE must be one of: reseq, hic, hifi, qc." >&2
        exit 2
        ;;
esac
