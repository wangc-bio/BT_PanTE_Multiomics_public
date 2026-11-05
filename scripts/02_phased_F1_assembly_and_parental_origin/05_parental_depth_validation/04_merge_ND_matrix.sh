#!/usr/bin/env bash
#SBATCH --job-name=04_merge_ND_matrix
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
: "${WINBED:?}"; : "${ND_DIR:?}"; : "${BSB_LIST:?}"; : "${TC_LIST:?}"; : "${OUTDIR:?}"
WIN=${WIN:-1000000}
mkdir -p "$OUTDIR/columns" logs; cd "$OUTDIR"
cut -f1-3 "$WINBED" > win.coords
: > columns.order; : > samples.order
i=0
for list in "$BSB_LIST" "$TC_LIST"; do
  while IFS= read -r bam || [[ -n "$bam" ]]; do
    [[ -z "$bam" || "$bam" == \#* ]] && continue
    sample=$(basename "$bam" | sed -E 's/\.filtered\.bam$//I; s/\.(bam|sam|fq|fastq)(\.gz)?$//I')
    nd="$ND_DIR/$sample.win$WIN.ND.tsv"
    [[ -s "$nd" ]] || { echo "Missing normalized-depth file: $nd" >&2; exit 1; }
    i=$((i+1)); col=$(printf "columns/%03d_%s.col" "$i" "$sample")
    cut -f4 "$nd" > "$col"
    printf '%s\n' "$col" >> columns.order
    printf '%s\n' "$sample" >> samples.order
  done < "$list"
done
{ printf 'chr\tstart\tend'; while IFS= read -r s; do printf '\t%s' "$s"; done < samples.order; printf '\n'; } > "ND_matrix.win$WIN.tsv"
paste win.coords $(cat columns.order) >> "ND_matrix.win$WIN.tsv"
