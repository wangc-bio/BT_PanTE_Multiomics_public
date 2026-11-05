#!/bin/bash
# Software version: BLAST+ v2.15.0+ (makeblastdb, blastn)
#SBATCH --job-name=02_run_pante_all_vs_all
#SBATCH --array=0-4%5
#SBATCH --cpus-per-task=36
#SBATCH --output=02_run_pante_all_vs_all_%A_%a.out
#SBATCH --error=02_run_pante_all_vs_all_%A_%a.err

set -euo pipefail

: "${PANTE_WORKDIR:?Set PANTE_WORKDIR to the Pan-TE working directory before submission}"
ROOT="${PANTE_WORKDIR}"
mkdir -p "${ROOT}/01_blast" "${ROOT}/02_cluster"
command -v makeblastdb >/dev/null 2>&1 || { echo "ERROR: makeblastdb not found in PATH" >&2; exit 1; }
command -v blastn >/dev/null 2>&1 || { echo "ERROR: blastn not found in PATH" >&2; exit 1; }

CLASSES=(
DNA
LTR
LINE
SINE
RC
)

TASK_ID="${SLURM_ARRAY_TASK_ID:?Submit this script as a SLURM array job}"
C=${CLASSES[$TASK_ID]}

FA=${ROOT}/00_input/${C}.fa

DB=${ROOT}/01_blast/${C}.db

OUT=${ROOT}/01_blast/${C}.all_vs_all.tsv

echo "===================================="
echo "Class : ${C}"
echo "Input : ${FA}"
echo "===================================="

makeblastdb \
    -in "${FA}" \
    -dbtype nucl \
    -parse_seqids \
    -out "${DB}"

blastn \
    -task blastn \
    -query "${FA}" \
    -db "${DB}" \
    -perc_identity 70 \
    -evalue 1e-10 \
    -dust no \
    -soft_masking false \
    -max_hsps 50 \
    -num_threads "${SLURM_CPUS_PER_TASK:-36}" \
    -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen' \
    > "${OUT}"

echo
echo "Finished ${C}"
wc -l "${OUT}"