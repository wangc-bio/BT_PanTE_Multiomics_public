#!/usr/bin/env bash
#SBATCH --job-name=01_MMseqs2_RBH
#SBATCH --cpus-per-task=36
#SBATCH --mem=160G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
: "${MMSEQS:?Path to MMseqs2}"; : "${BTM_PROTEINS:?}"; : "${BTP_PROTEINS:?}"; : "${OUTDIR:?}"
THREADS=${SLURM_CPUS_PER_TASK:-36}
mkdir -p "$OUTDIR" logs
cd "$OUTDIR"
"$MMSEQS" easy-rbh "$BTM_PROTEINS" "$BTP_PROTEINS" BT_subM-P_RBH.tsv tmp_rbh \
  --threads "$THREADS" -s 7.5 -e 1e-5 --min-seq-id 0.5 -c 0.6 --cov-mode 0 \
  --format-output "query,target,pident,alnlen,qlen,tlen,qcov,tcov,evalue,bits"
echo "# RBH pairs: $(wc -l < BT_subM-P_RBH.tsv)"
