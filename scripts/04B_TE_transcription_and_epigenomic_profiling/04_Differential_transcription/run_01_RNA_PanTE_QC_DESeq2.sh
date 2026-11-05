#!/bin/bash
# Public-release note: original executed source; replace local absolute paths/environment activation with your local configuration before running.
#SBATCH -J PanTE_RNA_DE
#SBATCH -p pms
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH -t 04:00:00
#SBATCH -o logs/PanTE_RNA_DE.%j.out
#SBATCH -e logs/PanTE_RNA_DE.%j.err

set -eo pipefail

source "${CONDA_SH:?Set CONDA_SH to your conda.sh path}"

conda activate "${CONDA_ENV:-R}"

RSCRIPT=${RSCRIPT:-Rscript}

echo "===== Environment ====="
echo "Rscript: ${RSCRIPT}"
"${RSCRIPT}" --version

echo
echo "===== Package check ====="

"${RSCRIPT}" -e '
cat("DESeq2:", as.character(packageVersion("DESeq2")), "\n")
cat("ggplot2:", as.character(packageVersion("ggplot2")), "\n")
cat("pheatmap:", as.character(packageVersion("pheatmap")), "\n")
'

echo
echo "===== Start PanTE RNA analysis ====="

"${RSCRIPT}" \
${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/PanTE_v52/05_TE_RNA/05_DE/01_RNA_PanTE_QC_DESeq2.R

echo
echo "===== Finished ====="
date