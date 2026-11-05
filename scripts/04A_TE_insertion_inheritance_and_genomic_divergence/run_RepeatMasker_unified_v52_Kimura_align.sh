#!/bin/bash
#SBATCH --job-name=RM_KIMURA_v52
#SBATCH --partition=pms
#SBATCH --array=0-3%4
#SBATCH --cpus-per-task=36
#SBATCH --mem=160G
#SBATCH --time=72:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

# ============================================================
# Environment
# ============================================================

# Activate the RepeatMasker environment and RMBlast before submission, or adapt these two lines to the local installation.
if [[ -n "${CONDA_SH:-}" ]]; then source "${CONDA_SH}"; conda activate "${CONDA_ENV:-RepeatMasker}"; fi
if [[ -n "${RMBLAST_PROFILE:-}" ]]; then source "${RMBLAST_PROFILE}"; fi


# ============================================================
# Paths
# ============================================================

BASE=${TE_BASE:?Set TE_BASE to the directory containing the four genome FASTA files and unified library}

ROOT=${KIMURA_RUN_DIR:-${BASE}/RepeatMasker_unified_v5_2_Kimura_align}

LIB=${TE_LIBRARY:-${BASE}/Final_evidence_v5/13_unified_v52/Unified_v5_2.curated_repeat.fa}

RM=${REPEATMASKER_EXECUTABLE:-$(command -v RepeatMasker)}


# ============================================================
# Results 2 genomes only
#
# BSB  = parental BSB genome
# TC   = parental TC genome
# BTM  = BSB-derived F1 haplotype
# BTP  = TC-derived F1 haplotype
# ============================================================

SAMPLES=(
BSB
TC
BTM
BTP
)

GENOMES=(
"${BASE}/RepeatModeler/BSB.genome.chr.fa"
"${BASE}/RepeatModeler/TC.genome.chr.fa"
"${BASE}/RepeatModeler/BT.genome.subM.fa"
"${BASE}/RepeatModeler/BT.genome.subP.fa"
)


# ============================================================
# Current sample
# ============================================================

S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
GENOME=${GENOMES[$SLURM_ARRAY_TASK_ID]}

OUTDIR=${ROOT}/results/${S}

mkdir -p "${OUTDIR}"


# ============================================================
# Information
# ============================================================

echo "============================================================"
echo "Unified v5.2 RepeatMasker for Kimura divergence"
echo "============================================================"
echo "Date     : $(date)"
echo "Hostname : $(hostname)"
echo "Task ID  : ${SLURM_ARRAY_TASK_ID}"
echo "Sample   : ${S}"
echo "Genome   : ${GENOME}"
echo "Library  : ${LIB}"
echo "RM       : ${RM}"
echo "Output   : ${OUTDIR}"
echo "============================================================"


# ============================================================
# Input checks
# ============================================================

if [[ ! -s "${GENOME}" ]]; then
    echo "ERROR: genome missing or empty:"
    echo "${GENOME}"
    exit 1
fi

if [[ ! -s "${LIB}" ]]; then
    echo "ERROR: unified TE library missing or empty:"
    echo "${LIB}"
    exit 1
fi

if [[ ! -x "${RM}" ]]; then
    echo "ERROR: RepeatMasker executable missing:"
    echo "${RM}"
    exit 1
fi


# ============================================================
# Basic QC
# ============================================================

echo
echo "===== Input QC ====="

echo -n "Genome sequences : "
grep -c '^>' "${GENOME}"

echo -n "Unified families : "
grep -c '^>' "${LIB}"

echo -n "Genome size (bp) : "
awk '
    /^>/ {next}
    {
        gsub(/[[:space:]]/, "")
        n += length($0)
    }
    END {print n}
' "${GENOME}"

echo
echo -n "Library MD5 : "
md5sum "${LIB}"

echo


# ============================================================
# RepeatMasker
#
# IMPORTANT:
# - Same unified library as the frozen Results 2 annotation
# - Same core parameters as previous run
# - ONLY new analytical option is "-a"
#
# -a writes the alignment file required by:
# calcDivergenceFromAlign.pl
# ============================================================

echo "===== Starting RepeatMasker ====="
echo "Start time: $(date)"
echo

"${RM}" \
    -pa 9 \
    -lib "${LIB}" \
    -gff \
    -html \
    -xsmall \
    -e rmblast \
    -no_is \
    -a \
    -dir "${OUTDIR}" \
    "${GENOME}"

echo
echo "RepeatMasker finished: $(date)"


# ============================================================
# Output checks
# ============================================================

echo
echo "===== Output check ====="


TBL=$(find "${OUTDIR}" \
    -maxdepth 1 \
    -type f \
    -name '*.tbl' \
    | head -1)

OUT=$(find "${OUTDIR}" \
    -maxdepth 1 \
    -type f \
    -name '*.out' \
    | head -1)

ALIGN=$(find "${OUTDIR}" \
    -maxdepth 1 \
    -type f \
    -name '*.align' \
    | head -1)


if [[ -z "${TBL}" || ! -s "${TBL}" ]]; then
    echo "ERROR: .tbl missing for ${S}"
    exit 1
fi

if [[ -z "${OUT}" || ! -s "${OUT}" ]]; then
    echo "ERROR: RepeatMasker .out missing for ${S}"
    exit 1
fi

if [[ -z "${ALIGN}" || ! -s "${ALIGN}" ]]; then
    echo "ERROR: .align missing for ${S}"
    echo "Kimura divergence cannot be calculated without .align"
    exit 1
fi


echo "[OK] TBL   : ${TBL}"
echo "[OK] OUT   : ${OUT}"
echo "[OK] ALIGN : ${ALIGN}"

echo
echo "===== Output sizes ====="

ls -lh \
    "${TBL}" \
    "${OUT}" \
    "${ALIGN}"


# ============================================================
# Show RepeatMasker summary
# ============================================================

echo
echo "===== RepeatMasker summary ====="
cat "${TBL}"


# ============================================================
# Alignment QC
# ============================================================

echo
echo "===== Alignment QC ====="

echo -n "ALIGN size (bytes): "
stat -c '%s' "${ALIGN}"

echo -n "ALIGN lines       : "
wc -l < "${ALIGN}"

echo
echo "First 15 lines:"
head -15 "${ALIGN}"


# ============================================================
# Finished
# ============================================================

echo
echo "============================================================"
echo "Finished sample: ${S}"
echo "Finish time     : $(date)"
echo "ALIGN file      : ${ALIGN}"
echo "============================================================"