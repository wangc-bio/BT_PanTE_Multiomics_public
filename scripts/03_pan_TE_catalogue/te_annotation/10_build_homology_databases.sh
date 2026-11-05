#!/bin/bash
# Software version: BLAST+ v2.15.0+ (makeblastdb, blastdbcmd)
#SBATCH --job-name=10_build_homology_databases
#SBATCH --cpus-per-task=4
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -eo pipefail
if [[ -z "$TE_CONFIG" ]]; then
    echo "ERROR: export TE_CONFIG=/path/to/te_annotation.env before submission" >&2
    exit 1
fi
source "$TE_CONFIG"

if [[ -z "$DFAM_REFERENCE_FASTA_1" || -z "$DFAM_REFERENCE_FASTA_2" || -z "$DFAM_DB_PREFIX_1" || -z "$DFAM_DB_PREFIX_2" ]]; then
    echo "ERROR: both Dfam reference FASTAs and database prefixes are required in TE_CONFIG" >&2
    exit 1
fi
set -u

for F in "$DFAM_REFERENCE_FASTA_1" "$DFAM_REFERENCE_FASTA_2"; do
    [[ -s "$F" ]] || { echo "ERROR: Dfam reference FASTA is missing or empty: $F" >&2; exit 1; }
done

command -v makeblastdb >/dev/null 2>&1 || { echo "ERROR: makeblastdb not found in PATH" >&2; exit 1; }
command -v blastdbcmd >/dev/null 2>&1 || { echo "ERROR: blastdbcmd not found in PATH" >&2; exit 1; }

build_database() {
    LABEL="$1"
    FASTA="$2"
    PREFIX="$3"

    mkdir -p "$(dirname "$PREFIX")"

    echo "Building $LABEL"
    echo "Input FASTA : $FASTA"
    echo "DB prefix   : $PREFIX"

    makeblastdb \
        -in "$FASTA" \
        -dbtype nucl \
        -out "$PREFIX" \
        > "$PREFIX.makeblastdb.log" 2>&1

    blastdbcmd -db "$PREFIX" -info > "$PREFIX.database_info.txt"
}

build_database "Dfam reference level 1" "$DFAM_REFERENCE_FASTA_1" "$DFAM_DB_PREFIX_1"
build_database "Dfam reference level 2" "$DFAM_REFERENCE_FASTA_2" "$DFAM_DB_PREFIX_2"

echo "Homology reference databases finished: $(date)"
