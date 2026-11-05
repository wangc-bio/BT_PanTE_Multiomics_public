#!/usr/bin/env bash
# Summarize STAR Log.final.out files into per-sample and weighted aggregate QC.
#
# Software:
#   STAR 2.7.10b
#
# Usage:
#   bash 02_summarize_star_logs.sh output.tsv \
#       SAMPLE_A=path/to/SAMPLE_A.Log.final.out \
#       SAMPLE_B=path/to/SAMPLE_B.Log.final.out
#
# Input paths are supplied by the user at runtime. No project-specific paths,
# usernames, sample FASTQ paths, or cluster configuration are embedded here.
#
# Overall mapping rate is calculated as:
#   (uniquely mapped + mapped to multiple loci + mapped to too many loci)
#   / number of input reads * 100
#
# The aggregate row uses summed read counts (a read-count-weighted rate),
# rather than the arithmetic mean of per-sample percentages.

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: bash $0 output.tsv SAMPLE_ID=STAR_Log.final.out [SAMPLE_ID=STAR_Log.final.out ...]" >&2
    exit 1
fi

OUTPUT_TSV="$1"
shift

mkdir -p "$(dirname "$OUTPUT_TSV")"

printf 'Sample_ID\tSTAR_version\tInput_reads\tUnique_mapped_reads\tUnique_mapped_pct\tMulti_locus_mapped_reads\tMulti_locus_mapped_pct\tToo_many_loci_reads\tTotal_mapped_reads\tOverall_mapping_rate_pct\tUnmapped_mismatch_reads\tUnmapped_too_short_reads\tUnmapped_other_reads\tTotal_unmapped_reads\tUnmapped_rate_pct\n' > "$OUTPUT_TSV"

total_input=0
total_unique=0
total_multi=0
total_too_many=0
total_unmapped_mismatch=0
total_unmapped_short=0
total_unmapped_other=0

for item in "$@"; do
    if [[ "$item" != *=* ]]; then
        echo "ERROR: expected SAMPLE_ID=STAR_Log.final.out, received: $item" >&2
        exit 2
    fi

    sample_id="${item%%=*}"
    star_log="${item#*=}"

    if [[ -z "$sample_id" || ! -s "$star_log" ]]; then
        echo "ERROR: sample ID is empty or STAR log is missing/empty: $item" >&2
        exit 2
    fi

    values=$(awk -F'|' '
        function trim(s) {
            sub(/^[[:space:]]+/, "", s)
            sub(/[[:space:]]+$/, "", s)
            return s
        }
        {
            key = trim($1)
            value = trim($2)
            gsub(/%/, "", value)

            if (key == "Number of input reads") input = value + 0
            else if (key == "Uniquely mapped reads number") unique = value + 0
            else if (key == "Uniquely mapped reads %") unique_pct = value + 0
            else if (key == "Number of reads mapped to multiple loci") multi = value + 0
            else if (key == "% of reads mapped to multiple loci") multi_pct = value + 0
            else if (key == "Number of reads mapped to too many loci") too_many = value + 0
            else if (key == "Number of reads unmapped: too many mismatches") unmapped_mismatch = value + 0
            else if (key == "Number of reads unmapped: too short") unmapped_short = value + 0
            else if (key == "Number of reads unmapped: other") unmapped_other = value + 0
        }
        END {
            printf "%d\t%d\t%.2f\t%d\t%.2f\t%d\t%d\t%d\t%d\n",
                input, unique, unique_pct, multi, multi_pct, too_many,
                unmapped_mismatch, unmapped_short, unmapped_other
        }
    ' "$star_log")

    read -r input unique unique_pct multi multi_pct too_many unmapped_mismatch unmapped_short unmapped_other <<< "$values"

    if [[ "$input" -le 0 ]]; then
        echo "ERROR: could not extract a valid 'Number of input reads' from: $star_log" >&2
        exit 3
    fi

    mapped=$((unique + multi + too_many))
    unmapped=$((unmapped_mismatch + unmapped_short + unmapped_other))
    classified=$((mapped + unmapped))
    if [[ "$classified" -ne "$input" ]]; then
        echo "ERROR: STAR mapping categories ($classified) do not sum to input reads ($input) for $sample_id." >&2
        exit 4
    fi

    overall_pct=$(awk -v n="$mapped" -v d="$input" 'BEGIN {printf "%.4f", n/d*100}')
    unmapped_pct=$(awk -v n="$unmapped" -v d="$input" 'BEGIN {printf "%.4f", n/d*100}')

    printf '%s\t2.7.10b\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\n' \
        "$sample_id" "$input" "$unique" "$unique_pct" "$multi" "$multi_pct" \
        "$too_many" "$mapped" "$overall_pct" "$unmapped_mismatch" \
        "$unmapped_short" "$unmapped_other" "$unmapped" "$unmapped_pct" \
        >> "$OUTPUT_TSV"

    total_input=$((total_input + input))
    total_unique=$((total_unique + unique))
    total_multi=$((total_multi + multi))
    total_too_many=$((total_too_many + too_many))
    total_unmapped_mismatch=$((total_unmapped_mismatch + unmapped_mismatch))
    total_unmapped_short=$((total_unmapped_short + unmapped_short))
    total_unmapped_other=$((total_unmapped_other + unmapped_other))
done

total_mapped=$((total_unique + total_multi + total_too_many))
total_unmapped=$((total_unmapped_mismatch + total_unmapped_short + total_unmapped_other))

aggregate_unique_pct=$(awk -v n="$total_unique" -v d="$total_input" 'BEGIN {printf "%.2f", n/d*100}')
aggregate_multi_pct=$(awk -v n="$total_multi" -v d="$total_input" 'BEGIN {printf "%.2f", n/d*100}')
aggregate_overall_pct=$(awk -v n="$total_mapped" -v d="$total_input" 'BEGIN {printf "%.4f", n/d*100}')
aggregate_unmapped_pct=$(awk -v n="$total_unmapped" -v d="$total_input" 'BEGIN {printf "%.4f", n/d*100}')

printf 'ALL_SAMPLES\t2.7.10b\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\n' \
    "$total_input" "$total_unique" "$aggregate_unique_pct" \
    "$total_multi" "$aggregate_multi_pct" "$total_too_many" \
    "$total_mapped" "$aggregate_overall_pct" "$total_unmapped_mismatch" \
    "$total_unmapped_short" "$total_unmapped_other" "$total_unmapped" \
    "$aggregate_unmapped_pct" >> "$OUTPUT_TSV"

echo "STAR mapping summary written to: $OUTPUT_TSV"
