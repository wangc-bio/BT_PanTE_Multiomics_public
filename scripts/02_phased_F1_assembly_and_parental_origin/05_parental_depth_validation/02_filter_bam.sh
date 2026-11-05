#!/usr/bin/env bash
#SBATCH --job-name=02_filter_bam
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --array=0-23
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
: "${SAMTOOLS:?}"; : "${BAM_LIST:?}"; : "${OUTDIR:?}"
mkdir -p "$OUTDIR" logs
bam=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$BAM_LIST")
[[ -s "$bam" ]] || { echo "Missing BAM: $bam" >&2; exit 1; }
sample=$(basename "$bam" .bam)
FLAG_FILTER=0x904
if [[ "${REMOVE_DUP:-0}" == "1" ]]; then FLAG_FILTER=0xD04; fi
"$SAMTOOLS" view -@ "${SLURM_CPUS_PER_TASK:-8}" -b -q 20 -f 2 -F "$FLAG_FILTER" "$bam" > "$OUTDIR/$sample.filtered.bam"
"$SAMTOOLS" index "$OUTDIR/$sample.filtered.bam"
