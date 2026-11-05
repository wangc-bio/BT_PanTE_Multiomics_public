#!/bin/bash
# Software version: BLAST+ v2.15.0+ (blastn, blastdbcmd)
#SBATCH --job-name=11_homology_rescue_search
#SBATCH --cpus-per-task=36
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -eo pipefail
if [[ -z "$TE_CONFIG" ]]; then
    echo "ERROR: export TE_CONFIG=/path/to/te_annotation.env before submission" >&2
    exit 1
fi
source "$TE_CONFIG"

if [[ -z "$HOMOLOGY_QUERY_FASTA" || -z "$REPBASE_DB_PREFIX" || -z "$DFAM_DB_PREFIX_1" || -z "$WORK_DIR" ]]; then
    echo "ERROR: HOMOLOGY_QUERY_FASTA, REPBASE_DB_PREFIX, DFAM_DB_PREFIX_1 and WORK_DIR are required" >&2
    exit 1
fi
export LC_ALL=C

THREADS="$SLURM_CPUS_PER_TASK"
if [[ -z "$THREADS" ]]; then THREADS=36; fi
set -u
OUT_DIR="$WORK_DIR/11_homology_rescue"

[[ -s "$HOMOLOGY_QUERY_FASTA" ]] || { echo "ERROR: homology query FASTA is missing: $HOMOLOGY_QUERY_FASTA" >&2; exit 1; }
command -v blastn >/dev/null 2>&1 || { echo "ERROR: blastn not found in PATH" >&2; exit 1; }
command -v blastdbcmd >/dev/null 2>&1 || { echo "ERROR: blastdbcmd not found in PATH" >&2; exit 1; }

mkdir -p "$OUT_DIR"

run_search() {
    LABEL="$1"
    QUERY_FASTA="$2"
    DB_PREFIX="$3"
    DB_DIR="$OUT_DIR/$LABEL"

    [[ -s "$QUERY_FASTA" ]] || { echo "ERROR: query FASTA for $LABEL is missing: $QUERY_FASTA" >&2; exit 1; }
    blastdbcmd -db "$DB_PREFIX" -info >/dev/null

    mkdir -p "$DB_DIR"

    blastn \
        -query "$QUERY_FASTA" \
        -db "$DB_PREFIX" \
        -task blastn \
        -evalue 1e-3 \
        -perc_identity 65 \
        -dust no \
        -soft_masking false \
        -num_threads "$THREADS" \
        -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore qcovs" \
        -out "$DB_DIR/$LABEL.raw.tsv"

    awk -F'\t' 'BEGIN{
        OFS="\t"
        print "Query_ID","Subject_ID","Identity","Alignment_length","Query_start","Query_end","Query_length","Subject_start","Subject_end","Subject_length","Evalue","Bit_score","Query_coverage","Query_covered_bp_from_qcovs"
    }
    NF>=13 {
        covered_bp=($13*$7)/100
    }
    NF>=13 && $3>=65 && $13>=20 && covered_bp>=80 && $11<=1e-3 {
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,covered_bp
    }' "$DB_DIR/$LABEL.raw.tsv" > "$DB_DIR/$LABEL.qualifying_hits.tsv"

    RAW_N=$(wc -l < "$DB_DIR/$LABEL.raw.tsv")
    QUAL_N=$(($(wc -l < "$DB_DIR/$LABEL.qualifying_hits.tsv") - 1))
    echo "$LABEL: raw hits=$RAW_N; qualifying hits=$QUAL_N"
}

# RepBase and the first Dfam reference level use the same unresolved query set.
run_search "RepBase" "$HOMOLOGY_QUERY_FASTA" "$REPBASE_DB_PREFIX"
run_search "Dfam_level_1" "$HOMOLOGY_QUERY_FASTA" "$DFAM_DB_PREFIX_1"

# Level 2 is run only when the unresolved-after-level-1 FASTA is supplied.
if [[ -n "$DFAM_LEVEL2_QUERY_FASTA" ]]; then
    if [[ -z "$DFAM_DB_PREFIX_2" ]]; then
        echo "ERROR: DFAM_DB_PREFIX_2 is required when DFAM_LEVEL2_QUERY_FASTA is set" >&2
        exit 1
    fi
    run_search "Dfam_level_2" "$DFAM_LEVEL2_QUERY_FASTA" "$DFAM_DB_PREFIX_2"
else
    echo "Dfam level 2 skipped: set DFAM_LEVEL2_QUERY_FASTA after first-level evidence integration."
fi

echo "Homology searches and hit-level filtering finished: $(date)"
echo "Broad-class winner selection, minimum subject support, bit-score fraction calculation, and RepBase/Dfam conflict resolution are not performed by this search script."
