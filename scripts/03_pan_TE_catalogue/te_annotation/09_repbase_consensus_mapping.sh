#!/bin/bash
# Software version: BLAST+ v2.15.0+ (blastn, blastdbcmd)
#SBATCH --job-name=09_repbase_consensus_mapping
#SBATCH --cpus-per-task=36
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail
export LC_ALL=C

: "${TE_CONFIG:?Export TE_CONFIG=/path/to/te_annotation.env before submission}"
source "${TE_CONFIG}"

: "${RM2_CLASSIFIED_FASTA:?RM2_CLASSIFIED_FASTA is required}"
: "${REPBASE_DB_PREFIX:?REPBASE_DB_PREFIX is required}"
: "${WORK_DIR:?WORK_DIR is required}"

THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-36}}"
OUT_DIR="${REPBASE_OUTPUT_DIR:-${WORK_DIR}/09_repbase_consensus_mapping}"

[[ -s "${RM2_CLASSIFIED_FASTA}" ]] || {
    echo "ERROR: RepeatModeler consensus FASTA is missing: ${RM2_CLASSIFIED_FASTA}" >&2
    exit 1
}
command -v blastn >/dev/null 2>&1 || {
    echo "ERROR: blastn not found in PATH" >&2
    exit 1
}
command -v blastdbcmd >/dev/null 2>&1 || {
    echo "ERROR: blastdbcmd not found in PATH" >&2
    exit 1
}

blastdbcmd -db "${REPBASE_DB_PREFIX}" -info >/dev/null
mkdir -p "${OUT_DIR}"

RAW="${OUT_DIR}/RepBase.raw.tsv"
SUP="${OUT_DIR}/RepBase.supportive.tsv"
HC="${OUT_DIR}/RepBase.HC.tsv"
BEST_SUP="${OUT_DIR}/RepBase.supportive.best.tsv"
BEST_HC="${OUT_DIR}/RepBase.HC.best.tsv"

echo "RepeatModeler query: ${RM2_CLASSIFIED_FASTA}"
echo "RepBase database   : ${REPBASE_DB_PREFIX}"
echo "Threads            : ${THREADS}"
echo "Output             : ${OUT_DIR}"

# Query the original RepeatModeler consensus library against RepBase.
# qcovs is subject-level query coverage and is filtered below.
blastn \
    -query "${RM2_CLASSIFIED_FASTA}" \
    -db "${REPBASE_DB_PREFIX}" \
    -task blastn \
    -evalue 1e-10 \
    -perc_identity 70 \
    -max_target_seqs 100 \
    -num_threads "${THREADS}" \
    -dust no \
    -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovs" \
    -out "${RAW}"

# Supportive evidence: identity >=70%, qcovs >=50%,
# alignment length >=100 bp and E-value <=1e-10.
awk -F'\t' 'BEGIN{OFS="\t"}
    $3 >= 70 && $4 >= 100 && $11 <= 1e-10 && $13 >= 50
' "${RAW}" > "${SUP}"

# High-confidence evidence: identity >=80%, qcovs >=80%,
# alignment length >=100 bp and E-value <=1e-20.
awk -F'\t' 'BEGIN{OFS="\t"}
    $3 >= 80 && $4 >= 100 && $11 <= 1e-20 && $13 >= 80
' "${RAW}" > "${HC}"

# Retain the best subject per RepeatModeler query:
# bit score, then identity, then qcovs.
for TYPE in supportive HC; do
    if [[ "${TYPE}" == "supportive" ]]; then
        INPUT="${SUP}"
        OUTPUT="${BEST_SUP}"
    else
        INPUT="${HC}"
        OUTPUT="${BEST_HC}"
    fi

    if [[ -s "${INPUT}" ]]; then
        sort -t $'\t' -k1,1 -k12,12gr -k3,3gr -k13,13gr "${INPUT}" \
            | awk -F'\t' '!seen[$1]++' > "${OUTPUT}"
    else
        : > "${OUTPUT}"
    fi
done

# Extract RepBase class and species annotations from subject identifiers.
for TYPE in supportive HC; do
    if [[ "${TYPE}" == "supportive" ]]; then
        INPUT="${BEST_SUP}"
    else
        INPUT="${BEST_HC}"
    fi

    awk -F'\t' 'BEGIN{
        OFS="\t"
        print "RM_family","RepBase_subject","pident","aln_len","qlen","slen","qstart","qend","sstart","send","evalue","bitscore","qcovs","RepBase_class","RepBase_species"
    }
    NF >= 13 {
        cls="NA"
        spp="NA"
        n=split($2,a,/[|]/)
        for(i=1;i<=n;i++){
            if(a[i] ~ /^RBclass=/){
                cls=a[i]
                sub(/^RBclass=/,"",cls)
            }
            if(a[i] ~ /^Species=/){
                spp=a[i]
                sub(/^Species=/,"",spp)
            }
        }
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,cls,spp
    }' "${INPUT}" > "${OUT_DIR}/RepBase.${TYPE}.best.annotated.tsv"
done

# Summarize RepBase classes among best hits.
for TYPE in supportive HC; do
    awk -F'\t' 'BEGIN{OFS="\t"}
        NR > 1 {count[$14]++}
        END {
            print "RepBase_class","RM_families"
            for (c in count) print c,count[c]
        }
    ' "${OUT_DIR}/RepBase.${TYPE}.best.annotated.tsv" \
      | { read -r header; printf '%s\n' "${header}"; sort -t $'\t' -k2,2nr; } \
      > "${OUT_DIR}/RepBase.${TYPE}.class_summary.tsv"
done

awk -F'\t' '$1 ~ /#Unknown/' "${BEST_SUP}" \
    > "${OUT_DIR}/RepBase.supportive.best.Unknown.tsv"
awk -F'\t' '$1 ~ /#Unknown/' "${BEST_HC}" \
    > "${OUT_DIR}/RepBase.HC.best.Unknown.tsv"

{
    printf "Metric\tValue\n"
    printf "Total_RM_families\t%s\n" "$(grep -c '^>' "${RM2_CLASSIFIED_FASTA}" || true)"
    printf "Original_Unknown\t%s\n" "$(grep '^>' "${RM2_CLASSIFIED_FASTA}" | grep -c '#Unknown' || true)"
    printf "Raw_BLAST_hits\t%s\n" "$(wc -l < "${RAW}")"
    printf "Supportive_hits\t%s\n" "$(wc -l < "${SUP}")"
    printf "Supportive_families\t%s\n" "$(wc -l < "${BEST_SUP}")"
    printf "High_confidence_hits\t%s\n" "$(wc -l < "${HC}")"
    printf "High_confidence_families\t%s\n" "$(wc -l < "${BEST_HC}")"
    printf "Supportive_Unknown_families\t%s\n" "$(wc -l < "${OUT_DIR}/RepBase.supportive.best.Unknown.tsv")"
    printf "High_confidence_Unknown_families\t%s\n" "$(wc -l < "${OUT_DIR}/RepBase.HC.best.Unknown.tsv")"
} > "${OUT_DIR}/RepBase.summary.tsv"

echo "RepBase consensus mapping finished: $(date)"
echo "Summary: ${OUT_DIR}/RepBase.summary.tsv"
