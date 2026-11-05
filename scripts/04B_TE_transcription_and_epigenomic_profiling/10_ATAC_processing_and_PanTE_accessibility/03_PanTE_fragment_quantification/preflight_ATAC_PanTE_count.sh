#!/bin/bash
# Public-release note: source script from the executed analysis. Replace the original local absolute paths and executable locations with paths in your environment before running.

set -eo pipefail

BASE=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/PanTE_v52
ATAC=${BASE}/07_TE_ATAC

SAMTOOLS="${SAMTOOLS:-samtools}"
BEDTOOLS=${BEDTOOLS:-bedtools}
PY=${PYTHON:-python}

META=${ATAC}/00_QC/ATAC_BSB_TC_BT_metadata.tsv

echo "=========================================="
echo "1. EXECUTABLE CHECK"
echo "=========================================="

for X in \
"${SAMTOOLS}" \
"${BEDTOOLS}" \
"${PY}"
do
    if [[ ! -x "${X}" ]]; then
        echo "FAIL: not executable: ${X}"
        exit 1
    else
        echo "PASS: ${X}"
    fi
done

echo ""
echo "samtools:"
${SAMTOOLS} --version | head -1

echo ""
echo "bedtools:"
${BEDTOOLS} --version

echo ""
echo "Python:"
${PY} --version


echo ""
echo "=========================================="
echo "2. PYTHON MODULE CHECK"
echo "=========================================="

${PY} - <<'PY'
import numpy
import pandas

print("numpy :", numpy.__version__)
print("pandas:", pandas.__version__)
print("PASS: Python modules")
PY


echo ""
echo "=========================================="
echo "3. INPUT FILE CHECK"
echo "=========================================="

FILES=(
"${META}"
"${ATAC}/01_annotation/BSB.PanTE.assignable_segments.bed"
"${ATAC}/01_annotation/TC.PanTE.assignable_segments.bed"
"${ATAC}/01_annotation/BT.PanTE.assignable_segments.bed"
"${ATAC}/01_annotation/PanTE_assignable_bp.tsv"
"${ATAC}/02_PanTE_count/normalize_ATAC_PanTE_counts.py"
)

for F in "${FILES[@]}"
do
    if [[ ! -s "${F}" ]]; then
        echo "FAIL: missing/empty: ${F}"
        exit 1
    else
        echo "PASS: ${F}"
    fi
done


echo ""
echo "=========================================="
echo "4. METADATA CHECK"
echo "=========================================="

N=$(($(wc -l < "${META}") - 1))

echo "Samples = ${N}"

if [[ "${N}" -ne 27 ]]; then
    echo "FAIL: expected 27 samples"
    exit 1
fi

echo "PASS: 27 samples"


echo ""
echo "=========================================="
echo "5. BAM INTEGRITY + SORT ORDER CHECK"
echo "=========================================="

tail -n +2 "${META}" | \
while IFS=$'\t' read SAMPLE GROUP TISSUE BAM BAM_EXISTS INDEX_EXISTS
do

    echo -n "${SAMPLE} ${GROUP} ${TISSUE}: "

    ${SAMTOOLS} quickcheck "${BAM}"

    SO=$(
        ${SAMTOOLS} view -H "${BAM}" \
        | awk -F'\t' '
        $1=="@HD" {
            for(i=2;i<=NF;i++)
                if($i ~ /^SO:/) {
                    sub(/^SO:/,"",$i);
                    print $i
                }
        }'
    )

    if [[ "${SO}" != "coordinate" ]]; then
        echo "FAIL: BAM sort order = ${SO}"
        exit 1
    fi

    echo "PASS"
done


echo ""
echo "=========================================="
echo "ALL PREFLIGHT CHECKS PASSED"
echo "=========================================="