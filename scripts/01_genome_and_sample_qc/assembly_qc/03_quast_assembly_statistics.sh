#!/usr/bin/env bash
#SBATCH --job-name=03_quast_assembly_statistics.sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# QUAST v5.2.0
# Calculate assembly statistics for the final BT assembly.
#
# Required environment variable:
#   ASSEMBLY_FASTA  - input final assembly FASTA
#
# Optional environment variables:
#   OUTPUT_DIR      - QUAST output directory (default: quast_out)
#   THREADS         - number of threads (default: 5)
#   QUAST_BIN       - QUAST executable (default: quast.py)
#
# Example:
#   sbatch --export=ALL,ASSEMBLY_FASTA=/path/to/final_assembly.fa,OUTPUT_DIR=quast_out \
#       03_quast_assembly_statistics.sh
#
# Run from the directory where the Slurm log files should be written.
# QUAST results are summarized in Supplementary Table S2; this script is
# provided as a path-portable record of the command and is not a new analysis.

set -euo pipefail

: "${ASSEMBLY_FASTA:?Please set ASSEMBLY_FASTA to the final assembly FASTA}"

OUTPUT_DIR="${OUTPUT_DIR:-quast_out}"
THREADS="${THREADS:-5}"
QUAST_BIN="${QUAST_BIN:-quast.py}"

if [[ ! -s "$ASSEMBLY_FASTA" ]]; then
    echo "ERROR: assembly FASTA does not exist or is empty: $ASSEMBLY_FASTA" >&2
    exit 1
fi

if ! command -v "$QUAST_BIN" >/dev/null 2>&1; then
    echo "ERROR: QUAST executable was not found in PATH: $QUAST_BIN" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT_DIR")"

"$QUAST_BIN" -t "$THREADS" -o "$OUTPUT_DIR" "$ASSEMBLY_FASTA"
