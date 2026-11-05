#!/bin/bash
# Public-release note: original executed source; replace local absolute paths/environment activation with your local configuration before running.
#SBATCH -J STAR_QC
#SBATCH -p pms
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH -o logs/STAR_QC.%j.out
#SBATCH -e logs/STAR_QC.%j.err

set -eo pipefail

RNA=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/PanTE_v52/05_TE_RNA

mkdir -p "${RNA}/logs"
mkdir -p "${RNA}/00_metadata"

${PYTHON:-python} - <<'PY'
import os
import csv
from pathlib import Path

meta = Path(
    os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52", "05_TE_RNA", "00_metadata", "sample_metadata.tsv")
)

out = Path(
    os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52", "05_TE_RNA", "00_metadata", "STAR_mapping_QC.tsv")
)

# STAR Log.final.out ?? -> ?????
wanted = {
    "Number of input reads": "Input_reads",
    "Average input read length": "Average_read_length",
    "Uniquely mapped reads number": "Unique_reads",
    "Uniquely mapped reads %": "Unique_percent",
    "Average mapped length": "Average_mapped_length",
    "Mismatch rate per base, %": "Mismatch_percent",
    "Number of reads mapped to multiple loci": "Multi_reads",
    "% of reads mapped to multiple loci": "Multi_percent",
    "Number of reads mapped to too many loci": "Too_many_loci_reads",
    "% of reads mapped to too many loci": "Too_many_loci_percent",
    "% of reads unmapped: too many mismatches": "Unmapped_mismatch_percent",
    "% of reads unmapped: too short": "Unmapped_short_percent",
    "% of reads unmapped: other": "Unmapped_other_percent",
}

records = []

with meta.open() as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:

        bam = Path(row["BAM"])

        log = Path(
            str(bam).replace(
                ".Aligned.sortedByCoord.out.bam",
                ".Log.final.out"
            )
        )

        rec = {
            "Sample": row["Sample"],
            "Group": row["Group"],
            "Individual": row["Individual"],
            "Tissue": row["Tissue"],
            "Input_reads": "NA",
            "Average_read_length": "NA",
            "Average_mapped_length": "NA",
            "Unique_reads": "NA",
            "Unique_percent": "NA",
            "Multi_reads": "NA",
            "Multi_percent": "NA",
            "Too_many_loci_reads": "NA",
            "Too_many_loci_percent": "NA",
            "Mismatch_percent": "NA",
            "Unmapped_mismatch_percent": "NA",
            "Unmapped_short_percent": "NA",
            "Unmapped_other_percent": "NA",
            "QC_status": "NA",
            "Log_final": str(log),
        }

        if not log.exists():

            rec["QC_status"] = "MISSING_LOG"
            records.append(rec)
            continue

        with log.open() as lf:

            for line in lf:

                if "|" not in line:
                    continue

                key, value = line.split("|", 1)

                key = key.strip()
                value = value.strip()

                if key in wanted:
                    rec[wanted[key]] = value

        rec["QC_status"] = "OK"

        records.append(rec)


cols = [
    "Sample",
    "Group",
    "Individual",
    "Tissue",
    "Input_reads",
    "Average_read_length",
    "Average_mapped_length",
    "Unique_reads",
    "Unique_percent",
    "Multi_reads",
    "Multi_percent",
    "Too_many_loci_reads",
    "Too_many_loci_percent",
    "Mismatch_percent",
    "Unmapped_mismatch_percent",
    "Unmapped_short_percent",
    "Unmapped_other_percent",
    "QC_status",
    "Log_final",
]


with out.open("w", newline="") as oh:

    writer = csv.DictWriter(
        oh,
        fieldnames=cols,
        delimiter="\t",
        extrasaction="raise"
    )

    writer.writeheader()
    writer.writerows(records)


print("===== STAR QC SUMMARY =====")
print("Total samples :", len(records))

n_ok = sum(
    x["QC_status"] == "OK"
    for x in records
)

n_missing = sum(
    x["QC_status"] != "OK"
    for x in records
)

print("Logs found    :", n_ok)
print("Missing logs  :", n_missing)
print("Output        :", out)

if len(records) != 57:
    print(
        "[WARNING] Expected 57 RNA-seq samples, "
        "but found {}".format(len(records))
    )

PY

echo
echo "===== Group x Tissue sample numbers ====="

awk -F'\t' '
NR>1 {
    a[$2 FS $4]++
}
END {
    for(k in a)
        print k,a[k]
}' "${RNA}/00_metadata/STAR_mapping_QC.tsv" \
| sort

echo
echo "===== Missing STAR logs ====="

awk -F'\t' '
NR>1 && $18!="OK" {
    print $1,$2,$4,$18,$19
}' OFS="\t" \
"${RNA}/00_metadata/STAR_mapping_QC.tsv" || true

echo
echo "===== DONE ====="