#!/usr/bin/env bash
#SBATCH --job-name=01_meryl_k21_count.sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=100G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Meryl version: 1.4.1
# Generate a k=21 count database from PacBio HiFi reads.
#
# Required environment variable:
#   HIFI_READS  - input HiFi FASTQ/FASTQ.GZ file
#
# Optional environment variables:
#   MERYL_DB    - output Meryl database (default: read.meryl)
#   THREADS     - Meryl thread count (default: 36)
#
# Example:
#   sbatch --export=ALL,HIFI_READS=/path/to/hifi_reads.fastq.gz,MERYL_DB=/path/to/read.meryl 01_meryl_k21_count.sh
#
# Run from the directory where the Slurm log files should be written.
# Ensure Meryl 1.4.1 is available in PATH on the compute node.

set -euo pipefail

: "${HIFI_READS:?Please set HIFI_READS to the input PacBio HiFi reads}"

MERYL_DB="${MERYL_DB:-read.meryl}"
THREADS="${THREADS:-36}"

if ! command -v meryl >/dev/null 2>&1; then
    echo "ERROR: meryl was not found in PATH. Please load Meryl 1.4.1." >&2
    exit 1
fi

if [[ ! -s "$HIFI_READS" ]]; then
    echo "ERROR: HiFi input file does not exist or is empty: $HIFI_READS" >&2
    exit 1
fi

mkdir -p "$(dirname "$MERYL_DB")"

meryl k=21 count threads="$THREADS" output "$MERYL_DB" "$HIFI_READS"
