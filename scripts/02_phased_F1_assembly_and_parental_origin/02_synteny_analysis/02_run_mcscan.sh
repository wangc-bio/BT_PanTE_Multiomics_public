#!/usr/bin/env bash
#SBATCH --job-name=02_run_mcscan
#SBATCH --cpus-per-task=36
#SBATCH --mem=160G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
: "${WORKDIR:?Set WORKDIR containing *.bed and *.cds files}"
: "${CONDA_SH:?}"; : "${JCVI_ENV:?}"
mkdir -p logs
cd "$WORKDIR"
source "$CONDA_SH"
conda activate "$JCVI_ENV"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-36}
# Six pairwise comparisons among BSB, TC, BT-M and BT-P.
for pair in BSB.BTP BSB.BTM BSB.TC BTM.BTP BTM.TC BTP.TC; do
  python -m jcvi.compara.catalog ortholog ${pair%%.*} ${pair#*.} --cscore=.80 --cpus 36
done
