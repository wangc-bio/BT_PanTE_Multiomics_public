#!/bin/bash
# Software version: BLAST+ v2.15.0+ (confirmed from the execution environment).
#SBATCH --job-name=08_classII_evidence_mapping
#SBATCH --cpus-per-task=36
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -eo pipefail
if [[ -z "$TE_CONFIG" ]]; then
    echo "ERROR: export TE_CONFIG=/path/to/te_annotation.env before submission" >&2
    exit 1
fi
source "$TE_CONFIG"

if [[ -z "$RM2_CLASSIFIED_FASTA" || -z "$MITE_FASTA" || -z "$HELITRON_FASTA" || -z "$WORK_DIR" ]]; then
    echo "ERROR: RM2_CLASSIFIED_FASTA, MITE_FASTA, HELITRON_FASTA and WORK_DIR are required in TE_CONFIG" >&2
    exit 1
fi
export LC_ALL=C

THREADS="$SLURM_CPUS_PER_TASK"
if [[ -z "$THREADS" ]]; then THREADS=36; fi
set -u
OUT_DIR="$WORK_DIR/08_classII_mapping"

for F in "$RM2_CLASSIFIED_FASTA" "$MITE_FASTA" "$HELITRON_FASTA"; do
    [[ -s "$F" ]] || { echo "ERROR: required FASTA is missing or empty: $F" >&2; exit 1; }
done
command -v makeblastdb >/dev/null 2>&1 || { echo "ERROR: makeblastdb not found in PATH" >&2; exit 1; }
command -v blastn >/dev/null 2>&1 || { echo "ERROR: blastn not found in PATH" >&2; exit 1; }

mkdir -p "$OUT_DIR"

run_mapping() {
    TYPE="$1"
    CANDIDATE_FASTA="$2"
    TYPE_DIR="$OUT_DIR/$TYPE"
    DB_PREFIX="$TYPE_DIR/$TYPE.db"
    RAW="$TYPE_DIR/$TYPE.raw.tsv"

    mkdir -p "$TYPE_DIR"

    makeblastdb \
        -in "$CANDIDATE_FASTA" \
        -dbtype nucl \
        -out "$DB_PREFIX" \
        > "$TYPE_DIR/makeblastdb.log" 2>&1

    blastn \
        -query "$RM2_CLASSIFIED_FASTA" \
        -db "$DB_PREFIX" \
        -task blastn \
        -evalue 1e-10 \
        -perc_identity 70 \
        -qcov_hsp_perc 50 \
        -max_target_seqs 100 \
        -num_threads "$THREADS" \
        -dust no \
        -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp" \
        -out "$RAW"

    # Supportive: identity >=70%, qcov_hsp >=50%, alignment >=100 bp, E <=1e-10.
    awk -F'\t' 'BEGIN{OFS="\t"}
        $3 >= 70 && $4 >= 100 && $11 <= 1e-10 && $13 >= 50
    ' "$RAW" > "$TYPE_DIR/$TYPE.supportive.tsv"

    # High confidence: identity >=80%, qcov_hsp >=80%, alignment >=100 bp, E <=1e-20.
    awk -F'\t' 'BEGIN{OFS="\t"}
        $3 >= 80 && $4 >= 100 && $11 <= 1e-20 && $13 >= 80
    ' "$RAW" > "$TYPE_DIR/$TYPE.HC.tsv"

    for LEVEL in supportive HC; do
        INPUT="$TYPE_DIR/$TYPE.$LEVEL.tsv"
        OUTPUT="$TYPE_DIR/$TYPE.$LEVEL.best.tsv"

        if [[ -s "$INPUT" ]]; then
            sort -t $'\t' \
                -k1,1 \
                -k12,12gr \
                -k3,3gr \
                -k13,13gr \
                "$INPUT" \
            | awk -F'\t' '!seen[$1]++' > "$OUTPUT"
        else
            : > "$OUTPUT"
        fi
    done

    printf "%s\t%s\t%s\t%s\n" \
        "$TYPE" \
        "$(grep -c '^>' "$CANDIDATE_FASTA" || true)" \
        "$(wc -l < "$TYPE_DIR/$TYPE.supportive.best.tsv")" \
        "$(wc -l < "$TYPE_DIR/$TYPE.HC.best.tsv")"
}

SUMMARY="$OUT_DIR/ClassII_mapping.summary.tsv"
printf "Evidence\tCandidate_sequences\tSupportive_RM_families\tHC_RM_families\n" > "$SUMMARY"

run_mapping MITE "$MITE_FASTA" >> "$SUMMARY"
run_mapping Helitron "$HELITRON_FASTA" >> "$SUMMARY"

cat "$SUMMARY"
echo "Class II evidence mapping finished: $(date)"
