#!/usr/bin/env python3
import os
# Public-release note: original executed source; replace local absolute paths with paths from your environment before running.

import csv
import re
import sys
from pathlib import Path
from collections import defaultdict

RNA = Path(
    os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52", "05_TE_RNA")
)

META = RNA / "00_metadata/sample_metadata.tsv"

MEMBERSHIP = Path(
    os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52", "02_cluster", "PanTE_refined_cluster_members.tsv")
)

TE_GTF = {
    "BSB": RNA / "01_annotation/BSB.PanTE.TEtranscripts.gtf",
    "TC":  RNA / "01_annotation/TC.PanTE.TEtranscripts.gtf",
    "BT":  RNA / "01_annotation/BT.PanTE.TEtranscripts.gtf",
}

OUTDIR = RNA / "04_matrix"
OUTDIR.mkdir(parents=True, exist_ok=True)

QC_OUT = OUTDIR / "TEcount_sample_QC.tsv"
MATRIX_OUT = OUTDIR / "PanTE_raw_counts.tsv"
PRES_OUT = OUTDIR / "PanTE_genome_presence.tsv"


# ============================================================
# 1. metadata
# ============================================================

with META.open() as f:
    metadata = list(csv.DictReader(f, delimiter="\t"))

if len(metadata) != 57:
    raise RuntimeError(
        "Expected 57 samples, found {}".format(len(metadata))
    )

samples = [x["Sample"] for x in metadata]

if len(samples) != len(set(samples)):
    raise RuntimeError("Duplicated sample names in metadata")


# ============================================================
# 2. all 11,950 PanTE IDs
# ============================================================

all_pante = set()

with MEMBERSHIP.open() as f:
    rd = csv.DictReader(f, delimiter="\t")
    for r in rd:
        all_pante.add(r["PanTE_ID"])

all_pante = sorted(all_pante)

print("Total homologous PanTE clusters:", len(all_pante))

if len(all_pante) != 11950:
    print(
        "[WARNING] Expected 11,950 PanTE clusters, found",
        len(all_pante)
    )


# ============================================================
# 3. genomic annotation presence from TE GTF
# ============================================================

gene_re = re.compile(r'gene_id\s+"([^"]+)"')

genome_presence = {}

for group, gtf in TE_GTF.items():

    present = set()

    with gtf.open() as f:
        for line in f:

            if not line.strip() or line.startswith("#"):
                continue

            m = gene_re.search(line)

            if not m:
                continue

            x = m.group(1)

            pm = re.match(r"(PanTE_R\d+)", x)

            if pm:
                present.add(pm.group(1))

    genome_presence[group] = present

    print(
        group,
        "annotated PanTE:",
        len(present)
    )


# ============================================================
# 4. read 57 TEcount tables
# ============================================================

matrix = {
    p: {s: 0 for s in samples}
    for p in all_pante
}

qc = []

failed = []

noninteger_values = []

for r in metadata:

    sample = r["Sample"]
    group = r["Group"]
    tissue = r["Tissue"]

    cnt = (
        RNA /
        "02_TEcount_multi" /
        group /
        tissue /
        sample /
        "{}.cntTable".format(sample)
    )

    if not cnt.exists() or cnt.stat().st_size == 0:

        failed.append(sample)

        qc.append({
            "Sample": sample,
            "Group": group,
            "Tissue": tissue,
            "Status": "MISSING",
            "PanTE_rows": 0,
            "Nonzero_PanTE": 0,
            "Total_PanTE_counts": 0,
            "Unexpected_PanTE": 0,
            "Count_table": str(cnt)
        })

        continue

    counts = defaultdict(float)

    with cnt.open() as f:

        header = f.readline().rstrip("\n").split("\t")

        if len(header) < 2:
            raise RuntimeError(
                "Unexpected cntTable header: {}".format(cnt)
            )

        for line in f:

            if not line.strip():
                continue

            p = line.rstrip("\n").split("\t")

            if len(p) < 2:
                continue

            feature = p[0]

            # Handles either:
            # PanTE_R000001
            # or PanTE_R000001:family:class
            m = re.match(r"(PanTE_R\d+)", feature)

            if not m:
                continue

            pan = m.group(1)

            try:
                val = float(p[1])
            except ValueError:
                raise RuntimeError(
                    "Non-numeric count in {}: {}".format(
                        cnt, line.rstrip()
                    )
                )

            counts[pan] += val

    unexpected = [
        p for p in counts
        if p not in genome_presence[group]
    ]

    nonzero = sum(
        1 for v in counts.values()
        if v > 0
    )

    total = sum(counts.values())

    # TEcount output used for DESeq2 should be integer-valued
    for p, v in counts.items():
        if abs(v - round(v)) > 1e-8:
            noninteger_values.append(
                (sample, p, v)
            )

    for p, v in counts.items():

        if p not in matrix:
            raise RuntimeError(
                "Unknown PanTE ID {} in {}".format(
                    p, sample
                )
            )

        matrix[p][sample] = v

    qc.append({
        "Sample": sample,
        "Group": group,
        "Tissue": tissue,
        "Status": "OK",
        "PanTE_rows": len(counts),
        "Nonzero_PanTE": nonzero,
        "Total_PanTE_counts": total,
        "Unexpected_PanTE": len(unexpected),
        "Count_table": str(cnt)
    })


# ============================================================
# 5. QC table
# ============================================================

qc_cols = [
    "Sample",
    "Group",
    "Tissue",
    "Status",
    "PanTE_rows",
    "Nonzero_PanTE",
    "Total_PanTE_counts",
    "Unexpected_PanTE",
    "Count_table"
]

with QC_OUT.open("w", newline="") as f:

    wr = csv.DictWriter(
        f,
        fieldnames=qc_cols,
        delimiter="\t"
    )

    wr.writeheader()
    wr.writerows(qc)


if failed:

    print("\nERROR: missing TEcount results:")
    for s in failed:
        print(" ", s)

    sys.exit(2)


# ============================================================
# 6. verify integer counts
# ============================================================

if noninteger_values:

    print(
        "\nERROR: non-integer TEcount values detected."
    )

    print("First 20:")

    for x in noninteger_values[:20]:
        print(x)

    print(
        "Do NOT round automatically before checking."
    )

    sys.exit(3)


# ============================================================
# 7. raw count matrix
# ============================================================

with MATRIX_OUT.open("w") as f:

    f.write(
        "PanTE_ID\t" +
        "\t".join(samples) +
        "\n"
    )

    for pan in all_pante:

        vals = [
            str(int(round(matrix[pan][s])))
            for s in samples
        ]

        f.write(
            pan + "\t" +
            "\t".join(vals) +
            "\n"
        )


# ============================================================
# 8. genomic presence matrix
# ============================================================

with PRES_OUT.open("w") as f:

    f.write(
        "PanTE_ID\tBSB_present\tTC_present\tBT_present\n"
    )

    for pan in all_pante:

        f.write(
            "{}\t{}\t{}\t{}\n".format(
                pan,
                1 if pan in genome_presence["BSB"] else 0,
                1 if pan in genome_presence["TC"] else 0,
                1 if pan in genome_presence["BT"] else 0
            )
        )


# ============================================================
# 9. summary
# ============================================================

print("")
print("===== TEcount merge completed =====")
print("Samples             :", len(samples))
print("PanTE clusters      :", len(all_pante))
print("Missing samples     :", len(failed))
print("Noninteger values   :", len(noninteger_values))

print("")
print("Outputs:")
print(QC_OUT)
print(MATRIX_OUT)
print(PRES_OUT)