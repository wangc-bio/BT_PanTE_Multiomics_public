#!/usr/bin/env bash
# Reciprocal subject-coverage filtering for MITE and Helitron Class II evidence.
# Input columns: qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp

set -euo pipefail
export LC_ALL=C

: "${CLASSII_MAPPING_DIR:?Set CLASSII_MAPPING_DIR to the 08_classII_mapping output directory}"

for TYPE in MITE Helitron; do
    TYPE_DIR="${CLASSII_MAPPING_DIR}/${TYPE}"
    RAW="${TYPE_DIR}/${TYPE}.raw.tsv"
    SUP="${TYPE_DIR}/${TYPE}.recip.supportive.tsv"
    HC="${TYPE_DIR}/${TYPE}.recip.HC.tsv"
    BESTSUP="${TYPE_DIR}/${TYPE}.recip.supportive.best.tsv"
    BESTHC="${TYPE_DIR}/${TYPE}.recip.HC.best.tsv"

    [[ -s "$RAW" ]] || { echo "ERROR: missing or empty raw table: $RAW" >&2; exit 1; }

    # Supportive: identity >=70%, alignment >=100 bp, E <=1e-10,
    # qcovhsp >=50%, and single-HSP subject coverage >=50%.
    awk -F'\t' '
        BEGIN{OFS="\t"}
        NF < 13 { next }
        { scov=($6>0 ? 100*$4/$6 : 0) }
        $3 >= 70 && $4 >= 100 && $11 <= 1e-10 && $13 >= 50 && scov >= 50 {
            print $0,scov
        }
    ' "$RAW" > "$SUP"

    # High confidence: identity >=80%, alignment >=100 bp, E <=1e-20,
    # qcovhsp >=80%, and single-HSP subject coverage >=80%.
    awk -F'\t' '
        BEGIN{OFS="\t"}
        NF < 13 { next }
        { scov=($6>0 ? 100*$4/$6 : 0) }
        $3 >= 80 && $4 >= 100 && $11 <= 1e-20 && $13 >= 80 && scov >= 80 {
            print $0,scov
        }
    ' "$RAW" > "$HC"

    for LEVEL in supportive HC; do
        if [[ "$LEVEL" == "supportive" ]]; then
            INPUT="$SUP"; OUTPUT="$BESTSUP"
        else
            INPUT="$HC"; OUTPUT="$BESTHC"
        fi
        if [[ -s "$INPUT" ]]; then
            sort -t $'\t' -k1,1 -k12,12gr -k3,3gr -k13,13gr -k14,14gr "$INPUT" \
                | awk -F'\t' '!seen[$1]++' > "$OUTPUT"
        else
            : > "$OUTPUT"
        fi
    done

    printf '%s\t%s\t%s\n' "$TYPE" \
        "$(wc -l < "$BESTSUP")" "$(wc -l < "$BESTHC")"
done
