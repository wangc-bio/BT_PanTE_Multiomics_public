#!/usr/bin/env bash
#SBATCH --job-name=02_merqury_assembly_evaluation.sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Merqury version: 1.4.1
# Run Merqury using a Meryl k-mer database and the target assembly FASTA.
#
# Required environment variables:
#   MERYL_DB       - Meryl database generated from HiFi reads
#   ASSEMBLY_FASTA - final assembly FASTA to be evaluated
#
# Optional environment variable:
#   OUTPUT_PREFIX  - output prefix (default: result)
#
# Example:
#   sbatch --export=ALL,MERYL_DB=/path/to/read.meryl,ASSEMBLY_FASTA=/path/to/assembly.fa,OUTPUT_PREFIX=result 02_merqury_assembly_evaluation.sh
#
# Run from the directory where the Merqury output should be written.
# Ensure Merqury 1.4.1 is available in PATH on the compute node.

set -euo pipefail

: "${MERYL_DB:?Please set MERYL_DB to the Meryl k-mer database}"
: "${ASSEMBLY_FASTA:?Please set ASSEMBLY_FASTA to the assembly FASTA}"

OUTPUT_PREFIX="${OUTPUT_PREFIX:-result}"

if ! command -v merqury.sh >/dev/null 2>&1; then
    echo "ERROR: merqury.sh was not found in PATH. Please load Merqury 1.4.1." >&2
    exit 1
fi

if [[ ! -e "$MERYL_DB" ]]; then
    echo "ERROR: Meryl database was not found: $MERYL_DB" >&2
    exit 1
fi

if [[ ! -s "$ASSEMBLY_FASTA" ]]; then
    echo "ERROR: assembly FASTA does not exist or is empty: $ASSEMBLY_FASTA" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

merqury.sh "$MERYL_DB" "$ASSEMBLY_FASTA" "$OUTPUT_PREFIX"
