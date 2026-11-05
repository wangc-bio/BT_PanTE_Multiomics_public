#!/usr/bin/env bash
#SBATCH --job-name=03_filter_anchors
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
: "${WORKDIR:?Set WORKDIR containing *.anchors}"; : "${CONDA_SH:?}"; : "${JCVI_ENV:?}"
mkdir -p logs
cd "$WORKDIR"
source "$CONDA_SH"; conda activate "$JCVI_ENV"
for pair in BSB.BTP BSB.BTM BSB.TC BTM.BTP BTM.TC BTP.TC; do
  python -m jcvi.compara.synteny screen --minspan=30 --simple "${pair}.anchors" "${pair}.anchors.simple"
done
