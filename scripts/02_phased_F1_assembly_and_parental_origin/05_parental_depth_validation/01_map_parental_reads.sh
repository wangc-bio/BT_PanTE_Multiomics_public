#!/usr/bin/env bash
#SBATCH --job-name=01_map_parental_reads
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --array=0-23
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
: "${BWA:?}"; : "${SAMTOOLS:?}"; : "${BT_REFERENCE:?}"; : "${READS_LIST:?}"; : "${OUTDIR:?}"
mkdir -p "$OUTDIR" logs
r1=$(sed -n "$((2*SLURM_ARRAY_TASK_ID+1))p" "$READS_LIST")
r2=$(sed -n "$((2*SLURM_ARRAY_TASK_ID+2))p" "$READS_LIST")
[[ -s "$r1" && -s "$r2" ]] || { echo "Missing paired FASTQ for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
sample=$(basename "$r1" | sed 's/_R1.fastp.fq.gz$//')
"$BWA" mem -t "${SLURM_CPUS_PER_TASK:-8}" "$BT_REFERENCE" "$r1" "$r2" | "$SAMTOOLS" sort -@ 4 -o "$OUTDIR/$sample.bam"
