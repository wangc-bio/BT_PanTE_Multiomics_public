#!/bin/bash
# Public-release note: original executed source; replace local absolute paths/environment activation with local configuration before running.
#SBATCH -J WGBS_TE_coordQC
#SBATCH -p pms
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH -t 01:00:00
#SBATCH -o logs/coordQC.%j.out
#SBATCH -e logs/coordQC.%j.err

set -eo pipefail

OUTROOT=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/PanTE_v52/06_TE_WGBS/00_QC
mkdir -p "${OUTROOT}"

# ============================================================
# WGBS reference FASTA
# ============================================================

BSB_REF=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/ref_rename/Mamb/index/BSB.genome.fa

TC_REF=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/ref_rename/Calb/index/TC.genome.fa

BT_REF=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/Pb_asm/451/purge_hap/haphic/one_step/04.build/anno/EVM_v2/final_gff2/01_rename/index/BT.genome.rename.fa

# ============================================================
# Unified RepeatMasker GFF
# ============================================================

BSB_GFF=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/RepeatMasker_unified_v5_2_fixed/results/BSB/BSB.genome.chr.fa.out.gff

TC_GFF=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/RepeatMasker_unified_v5_2_fixed/results/TC/TC.genome.chr.fa.out.gff

BTM_GFF=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/RepeatMasker_unified_v5_2_fixed/results/BTM/BT.genome.subM.fa.out.gff

BTP_GFF=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/RepeatMasker_unified_v5_2_fixed/results/BTP/BT.genome.subP.fa.out.gff


echo "===== START ====="
date

# ============================================================
# Check input
# ============================================================

for F in \
    "${BSB_REF}" \
    "${TC_REF}" \
    "${BT_REF}" \
    "${BSB_GFF}" \
    "${TC_GFF}" \
    "${BTM_GFF}" \
    "${BTP_GFF}"
do
    if [[ ! -s "${F}" ]]; then
        echo "[ERROR] Missing input: ${F}" >&2
        exit 1
    fi
done

# ============================================================
# FASTA lengths
# ============================================================

for REF in "${BSB_REF}" "${TC_REF}" "${BT_REF}"
do
    if [[ ! -s "${REF}.fai" ]]; then
        echo "[INFO] Building fai: ${REF}"
        samtools faidx "${REF}"
    fi
done

awk 'BEGIN{OFS="\t"} $1~/^chr[0-9][0-9]$/ {print $1,$2}' \
    "${BSB_REF}.fai" \
    | sort -k1,1 \
    > "${OUTROOT}/BSB.WGBS_ref.chrom_lengths.tsv"

awk 'BEGIN{OFS="\t"} $1~/^chr[0-9][0-9]$/ {print $1,$2}' \
    "${TC_REF}.fai" \
    | sort -k1,1 \
    > "${OUTROOT}/TC.WGBS_ref.chrom_lengths.tsv"

awk 'BEGIN{OFS="\t"} $1~/^chr[0-9][0-9][MP]$/ {print $1,$2}' \
    "${BT_REF}.fai" \
    | sort -k1,1 \
    > "${OUTROOT}/BT.WGBS_ref.chrom_lengths.tsv"

# ============================================================
# RepeatMasker sequence-region lengths
# ============================================================

awk '
BEGIN{OFS="\t"}
$1=="##sequence-region" {
    if($2~/^chr[0-9][0-9]$/)
        print $2,$4
}' "${BSB_GFF}" \
| sort -k1,1 \
> "${OUTROOT}/BSB.PanTE.chrom_lengths.tsv"

awk '
BEGIN{OFS="\t"}
$1=="##sequence-region" {
    if($2~/^chr[0-9][0-9]$/)
        print $2,$4
}' "${TC_GFF}" \
| sort -k1,1 \
> "${OUTROOT}/TC.PanTE.chrom_lengths.tsv"

{
    awk '
    BEGIN{OFS="\t"}
    $1=="##sequence-region" {
        if($2~/^chr[0-9][0-9]M$/)
            print $2,$4
    }' "${BTM_GFF}"

    awk '
    BEGIN{OFS="\t"}
    $1=="##sequence-region" {
        if($2~/^chr[0-9][0-9]P$/)
            print $2,$4
    }' "${BTP_GFF}"
} \
| sort -k1,1 \
> "${OUTROOT}/BT.PanTE.chrom_lengths.tsv"

# ============================================================
# Python exact comparison
# ============================================================

${PYTHON:-python} - <<'PY'
import os
from pathlib import Path
import csv

root = Path(
    os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52", "06_TE_WGBS", "00_QC")
)

groups = ["BSB", "TC", "BT"]

summary = []

def read_lengths(path):
    d = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            chrom, length = line.rstrip().split("\t")[:2]
            d[chrom] = int(length)
    return d

for group in groups:

    ref_file = root / f"{group}.WGBS_ref.chrom_lengths.tsv"
    te_file  = root / f"{group}.PanTE.chrom_lengths.tsv"

    ref = read_lengths(ref_file)
    te  = read_lengths(te_file)

    chromosomes = sorted(set(ref) | set(te))

    out = root / f"{group}.coordinate_comparison.tsv"

    n_match = 0
    n_mismatch = 0
    n_missing_ref = 0
    n_missing_te = 0

    with open(out, "w") as fh:
        fh.write(
            "Chromosome\tWGBS_ref_length\tPanTE_length\t"
            "Difference\tStatus\n"
        )

        for chrom in chromosomes:

            r = ref.get(chrom)
            t = te.get(chrom)

            if r is None:
                status = "MISSING_IN_WGBS_REF"
                diff = "NA"
                n_missing_ref += 1

            elif t is None:
                status = "MISSING_IN_PanTE"
                diff = "NA"
                n_missing_te += 1

            else:
                diff = r - t

                if diff == 0:
                    status = "MATCH"
                    n_match += 1
                else:
                    status = "LENGTH_MISMATCH"
                    n_mismatch += 1

            fh.write(
                f"{chrom}\t"
                f"{r if r is not None else 'NA'}\t"
                f"{t if t is not None else 'NA'}\t"
                f"{diff}\t{status}\n"
            )

    expected = 48 if group == "BT" else 24

    if (
        n_match == expected
        and n_mismatch == 0
        and n_missing_ref == 0
        and n_missing_te == 0
    ):
        overall = "PASS"
    else:
        overall = "FAIL"

    summary.append([
        group,
        expected,
        len(ref),
        len(te),
        n_match,
        n_mismatch,
        n_missing_ref,
        n_missing_te,
        overall
    ])


summary_file = root / "coordinate_QC_summary.tsv"

with open(summary_file, "w") as fh:

    fh.write(
        "Group\tExpected_chr\tWGBS_ref_chr\tPanTE_chr\t"
        "Matched_chr\tLength_mismatch\t"
        "Missing_in_WGBS_ref\tMissing_in_PanTE\tStatus\n"
    )

    for row in summary:
        fh.write("\t".join(map(str,row)) + "\n")

print("")
print("===== Coordinate QC summary =====")

for row in summary:
    print("\t".join(map(str,row)))

print("")
print("Output:", summary_file)

if any(row[-1] != "PASS" for row in summary):
    raise SystemExit(
        "ERROR: At least one genome failed coordinate QC."
    )

print("")
print("ALL GENOMES PASSED.")
PY

echo
echo "===== FINAL TABLE ====="

column -t "${OUTROOT}/coordinate_QC_summary.tsv"

echo
echo "===== END ====="
date