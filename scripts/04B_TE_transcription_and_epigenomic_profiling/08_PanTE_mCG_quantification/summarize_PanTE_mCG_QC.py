#!/usr/bin/env python
# Public-release note: original executed source; replace local absolute paths with local paths before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import glob
import os
import pandas as pd


BASE = os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52")
OUT = os.path.join(BASE, "06_TE_WGBS")

files = sorted(
    glob.glob(
        os.path.join(
            OUT,
            "00_QC/*.PanTE_mCG_QC.tsv"
        )
    )
)

if len(files) == 0:
    raise RuntimeError(
        "No *.PanTE_mCG_QC.tsv files found."
    )

rows = []

for f in files:

    d = {}

    with open(f, "r") as fh:

        # skip header
        fh.readline()

        for line in fh:

            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 2:
                continue

            d[parts[0]] = parts[1]

    rows.append({
        "Sample":
            d.get("Sample", ""),

        "Group":
            d.get("Group", ""),

        "Tissue":
            d.get("Tissue", ""),

        "Total_cov_records":
            int(d.get(
                "Total_cov_records",
                0
            )),

        "CpG_overlapping_PanTE":
            int(d.get(
                "CpG_records_overlapping_PanTE",
                0
            )),

        "Unique_assigned_CpG":
            int(d.get(
                "Unique_assigned_CpG_records",
                0
            )),

        "Ambiguous_CpG":
            int(d.get(
                "Ambiguous_CpG_records",
                0
            )),

        "Unique_assignment_percent":
            float(d.get(
                "Unique_assignment_percent_of_TE_CpGs",
                0
            )),

        "Ambiguous_percent":
            float(d.get(
                "Ambiguous_percent_of_TE_CpGs",
                0
            )),

        "Quantified_PanTE_rows":
            int(d.get(
                "Quantified_PanTE_state_rows",
                0
            ))
    })


df = pd.DataFrame(rows)

df = df.sort_values(
    [
        "Group",
        "Sample",
        "Tissue"
    ]
)


outfile = os.path.join(
    OUT,
    "00_QC/PanTE_mCG_QC_summary.tsv"
)

df.to_csv(
    outfile,
    sep="\t",
    index=False
)


print("")
print("===== QC FILE COUNT =====")
print(len(files))


print("")
print("===== PER-SAMPLE QC =====")
print(
    df.to_string(
        index=False
    )
)


print("")
print("===== OVERALL SUMMARY =====")

for col in [
    "Unique_assignment_percent",
    "Ambiguous_percent",
    "Quantified_PanTE_rows"
]:

    print(
        "{}: min={:.3f}, median={:.3f}, max={:.3f}".format(
            col,
            df[col].min(),
            df[col].median(),
            df[col].max()
        )
    )


print("")
print("===== BY GROUP =====")

print(
    df.groupby("Group")[
        [
            "Unique_assignment_percent",
            "Ambiguous_percent",
            "Quantified_PanTE_rows"
        ]
    ].median()
)


print("")
print("===== BY TISSUE =====")

print(
    df.groupby("Tissue")[
        [
            "Unique_assignment_percent",
            "Ambiguous_percent",
            "Quantified_PanTE_rows"
        ]
    ].median()
)


print("")
print("Saved:")
print(outfile)