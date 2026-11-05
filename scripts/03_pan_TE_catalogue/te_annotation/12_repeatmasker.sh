#!/bin/bash
# Software versions: RepeatMasker v4.1.8, RMBlast v2.11.0+
#SBATCH --job-name=12_repeatmasker
#SBATCH --cpus-per-task=9
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -eo pipefail
if [[ -z "$TE_CONFIG" ]]; then
    echo "ERROR: export TE_CONFIG=/path/to/te_annotation.env before submission" >&2
    exit 1
fi
source "$TE_CONFIG"

if [[ -z "$GENOME_FASTA" || -z "$WORK_DIR" || -z "$REPEATMASKER_LIBRARY_FASTA" ]]; then
    echo "ERROR: GENOME_FASTA, WORK_DIR and REPEATMASKER_LIBRARY_FASTA are required" >&2
    exit 1
fi
set -u

OUT_DIR="${REPEATMASKER_OUT_DIR:-}"
if [[ -z "$OUT_DIR" ]]; then OUT_DIR="$WORK_DIR/12_repeatmasker"; fi

for F in "$GENOME_FASTA" "$REPEATMASKER_LIBRARY_FASTA"; do
    [[ -s "$F" ]] || { echo "ERROR: required input is missing or empty: $F" >&2; exit 1; }
done
command -v RepeatMasker >/dev/null 2>&1 || { echo "ERROR: RepeatMasker not found in PATH" >&2; exit 1; }

mkdir -p "$OUT_DIR"

RepeatMasker \
    -pa 9 \
    -lib "$REPEATMASKER_LIBRARY_FASTA" \
    -gff \
    -html \
    -xsmall \
    -e rmblast \
    -no_is \
    -dir "$OUT_DIR" \
    "$GENOME_FASTA"

echo "RepeatMasker finished: $(date)"
