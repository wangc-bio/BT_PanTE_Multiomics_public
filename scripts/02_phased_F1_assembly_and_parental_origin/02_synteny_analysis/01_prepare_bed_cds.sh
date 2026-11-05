#!/usr/bin/env bash
#SBATCH --job-name=01_prepare_bed_cds
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail

# Required environment variables: BSB_GFF, TC_GFF, BTM_GFF, BTP_GFF,
# BSB_CDS, TC_CDS, BTM_CDS, BTP_CDS, OUTDIR.
: "${BSB_GFF:?}"; : "${TC_GFF:?}"; : "${BTM_GFF:?}"; : "${BTP_GFF:?}"
: "${BSB_CDS:?}"; : "${TC_CDS:?}"; : "${BTM_CDS:?}"; : "${BTP_CDS:?}"
: "${OUTDIR:?}"
mkdir -p "$OUTDIR" logs
cd "$OUTDIR"

source "${CONDA_SH:?Set CONDA_SH to conda.sh}"
conda activate "${JCVI_ENV:?Set JCVI_ENV}"
for spec in "BSB:$BSB_GFF" "TC:$TC_GFF" "BTM:$BTM_GFF" "BTP:$BTP_GFF"; do
  name=${spec%%:*}; gff=${spec#*:}
  python -m jcvi.formats.gff bed --type=mRNA --key=ID "$gff" -o "${name}.bed"
  python -m jcvi.formats.bed uniq "${name}.bed"
done
conda deactivate
source "${SEQKIT_CONDA_SH:?Set SEQKIT_CONDA_SH}"
conda activate "${SEQKIT_ENV:?Set SEQKIT_ENV}"
for name in BSB TC BTM BTP; do
  case "$name" in
    BSB) cds=$BSB_CDS;; TC) cds=$TC_CDS;; BTM) cds=$BTM_CDS;; BTP) cds=$BTP_CDS;;
  esac
  seqkit grep -f <(cut -f4 "${name}.uniq.bed") "$cds" > "${name}.cds"
done
