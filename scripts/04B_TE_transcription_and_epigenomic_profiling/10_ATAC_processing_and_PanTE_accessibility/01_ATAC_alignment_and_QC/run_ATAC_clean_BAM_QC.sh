#!/bin/bash
# Public-release note: source script from the executed analysis. Replace the original local absolute paths and executable locations with paths in your environment before running.
#SBATCH -J ATAC_PanTE_QC
#SBATCH -p pms
#SBATCH -c 12
#SBATCH --mem=30G
#SBATCH -t 04:00:00
#SBATCH -o logs/ATAC_PanTE_QC.%j.out
#SBATCH -e logs/ATAC_PanTE_QC.%j.err

set -eo pipefail

${PYTHON:-python} \
${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/PanTE_v52/07_TE_ATAC/00_QC/qc_ATAC_clean_BAMs.py