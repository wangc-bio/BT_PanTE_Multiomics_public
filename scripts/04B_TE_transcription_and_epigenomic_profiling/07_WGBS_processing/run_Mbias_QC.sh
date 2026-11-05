#!/bin/bash
# Public-release note: original executed source; replace local absolute paths/environment activation with local configuration before running.
#SBATCH -J WGBS_Mbias_QC
#SBATCH -p pms
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH -o logs/Mbias_QC.%j.out
#SBATCH -e logs/Mbias_QC.%j.err

set -eo pipefail

WORKDIR=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/WGBS/00_Data
META=${WORKDIR}/wgbs_metadata.tsv
METH=${WORKDIR}/03_methylation_extract_v2

OUTROOT=${FB_ROOT:?Set FB_ROOT to the parent FB directory}/TE/PanTE_v52/06_TE_WGBS/00_QC

mkdir -p "${OUTROOT}"

${PYTHON:-python} - <<'PY'
import os
from pathlib import Path
import csv
import glob
import statistics

WORKDIR = Path(os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "WGBS", "00_Data"))
META = WORKDIR / "wgbs_metadata.tsv"
METH = WORKDIR / "03_methylation_extract_v2"

OUTROOT = Path(
    os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52", "06_TE_WGBS", "00_QC")
)

OUT_ALL = OUTROOT / "M_bias_CpG_all_positions.tsv"
OUT_SUM = OUTROOT / "M_bias_CpG_sample_summary.tsv"


def read_metadata():
    rows = []
    with META.open() as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            rows.append(r)
    return rows


def parse_mbias(path):
    """
    Returns:
      {
        "R1": [(pos,meth,unmeth,pct,cov), ...],
        "R2": [...]
      }

    Only CpG sections are retained.
    """
    data = {"R1": [], "R2": []}
    current = None

    with open(path, "r", errors="ignore") as f:
        for raw in f:
            line = raw.strip()

            if line.startswith("CpG context (R1)"):
                current = "R1"
                continue

            if line.startswith("CpG context (R2)"):
                current = "R2"
                continue

            # Once another context starts, stop collecting CpG
            if (
                "context (R1)" in line
                or "context (R2)" in line
            ) and not line.startswith("CpG context"):
                current = None
                continue

            if current is None:
                continue

            if not line or line.startswith("="):
                continue

            if line.startswith("position"):
                continue

            p = line.split("\t")

            if len(p) < 5:
                continue

            try:
                pos = int(p[0])
                meth = int(p[1])
                unmeth = int(p[2])
                pct = float(p[3])
                cov = int(p[4])
            except ValueError:
                continue

            data[current].append(
                (pos, meth, unmeth, pct, cov)
            )

    return data


def weighted_pct(rows):
    m = sum(x[1] for x in rows)
    u = sum(x[2] for x in rows)

    if m + u == 0:
        return None

    return 100.0 * m / (m + u)


def mean_pct(rows):
    vals = [x[3] for x in rows if x[4] > 0]

    if not vals:
        return None

    return statistics.mean(vals)


meta = read_metadata()

all_rows = []
summary = []

for r in meta:

    sample = r["sample_id"]
    species = r["species"]
    tissue = r["tissue"]
    individual = r["individual_id"]

    hits = glob.glob(
        str(METH / sample / "*M-bias.txt")
    )

    if not hits:

        summary.append({
            "sample_id": sample,
            "individual_id": individual,
            "species": species,
            "tissue": tissue,
            "read": "NA",
            "read_length": "NA",
            "central_weighted_mCG": "NA",
            "first1_mCG": "NA",
            "first2_weighted_mCG": "NA",
            "first5_weighted_mCG": "NA",
            "last5_weighted_mCG": "NA",
            "first2_diff_pp": "NA",
            "first5_diff_pp": "NA",
            "last5_diff_pp": "NA",
            "heuristic_flag": "MISSING_MBIAS",
            "mbias_file": "NA"
        })

        continue

    mbias = hits[0]
    dat = parse_mbias(mbias)

    for read in ["R1", "R2"]:

        rows = sorted(
            dat.get(read, []),
            key=lambda x: x[0]
        )

        if not rows:
            continue

        n = max(x[0] for x in rows)

        # central 80% of the read:
        # exclude first 10% and last 10%
        lo = max(1, int(n * 0.10) + 1)
        hi = max(lo, int(n * 0.90))

        central = [
            x for x in rows
            if lo <= x[0] <= hi
        ]

        first1 = [
            x for x in rows
            if x[0] == 1
        ]

        first2 = [
            x for x in rows
            if x[0] <= 2
        ]

        first5 = [
            x for x in rows
            if x[0] <= 5
        ]

        last5 = [
            x for x in rows
            if x[0] > n - 5
        ]

        base = weighted_pct(central)
        f1 = weighted_pct(first1)
        f2 = weighted_pct(first2)
        f5 = weighted_pct(first5)
        l5 = weighted_pct(last5)

        def diff(x):
            if x is None or base is None:
                return None
            return x - base

        d2 = diff(f2)
        d5 = diff(f5)
        dl5 = diff(l5)

        # This is ONLY a screening heuristic.
        # It is not a Bismark-defined cutoff.
        extreme = []

        if d2 is not None and abs(d2) >= 5:
            extreme.append("FIRST2")

        if d5 is not None and abs(d5) >= 5:
            extreme.append("FIRST5")

        if dl5 is not None and abs(dl5) >= 5:
            extreme.append("LAST5")

        flag = (
            ",".join(extreme)
            if extreme
            else "NONE"
        )

        summary.append({
            "sample_id": sample,
            "individual_id": individual,
            "species": species,
            "tissue": tissue,
            "read": read,
            "read_length": n,
            "central_weighted_mCG": (
                "%.3f" % base if base is not None else "NA"
            ),
            "first1_mCG": (
                "%.3f" % f1 if f1 is not None else "NA"
            ),
            "first2_weighted_mCG": (
                "%.3f" % f2 if f2 is not None else "NA"
            ),
            "first5_weighted_mCG": (
                "%.3f" % f5 if f5 is not None else "NA"
            ),
            "last5_weighted_mCG": (
                "%.3f" % l5 if l5 is not None else "NA"
            ),
            "first2_diff_pp": (
                "%.3f" % d2 if d2 is not None else "NA"
            ),
            "first5_diff_pp": (
                "%.3f" % d5 if d5 is not None else "NA"
            ),
            "last5_diff_pp": (
                "%.3f" % dl5 if dl5 is not None else "NA"
            ),
            "heuristic_flag": flag,
            "mbias_file": mbias
        })

        for pos, m, u, pct, cov in rows:
            all_rows.append([
                sample,
                individual,
                species,
                tissue,
                read,
                pos,
                m,
                u,
                pct,
                cov
            ])


with OUT_ALL.open("w") as out:

    out.write(
        "sample_id\tindividual_id\tspecies\ttissue\t"
        "read\tposition\tmethylated\tunmethylated\t"
        "mCG_percent\tcoverage\n"
    )

    for x in all_rows:
        out.write(
            "\t".join(map(str, x)) + "\n"
        )


fields = [
    "sample_id",
    "individual_id",
    "species",
    "tissue",
    "read",
    "read_length",
    "central_weighted_mCG",
    "first1_mCG",
    "first2_weighted_mCG",
    "first5_weighted_mCG",
    "last5_weighted_mCG",
    "first2_diff_pp",
    "first5_diff_pp",
    "last5_diff_pp",
    "heuristic_flag",
    "mbias_file",
]

with OUT_SUM.open("w", newline="") as out:

    wr = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t"
    )

    wr.writeheader()
    wr.writerows(summary)


print("===== M-bias QC =====")
print("metadata samples :", len(meta))
print("position rows    :", len(all_rows))
print("summary rows     :", len(summary))
print("output           :", OUT_SUM)

missing = [
    x for x in summary
    if x["heuristic_flag"] == "MISSING_MBIAS"
]

print("missing M-bias   :", len(missing))

PY

echo
echo "===== Sample summary ====="

column -t "${OUTROOT}/M_bias_CpG_sample_summary.tsv" \
| head -30

echo
echo "===== Flag counts ====="

awk -F'\t' '
NR>1 {
    a[$3 FS $5 FS $15]++
}
END {
    for(k in a)
        print k,a[k]
}' \
"${OUTROOT}/M_bias_CpG_sample_summary.tsv" \
| sort

echo
echo "===== DONE ====="