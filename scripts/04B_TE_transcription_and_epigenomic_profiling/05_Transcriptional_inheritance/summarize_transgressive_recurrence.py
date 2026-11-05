#!/usr/bin/env python
# Public-release note: original executed source; replace local absolute paths with paths from your environment before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import pandas as pd

BASE = os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52")
INFILE = os.path.join(
    BASE,
    "05_TE_RNA/05_DE/04_hybrid_pattern/PanTE_inheritance_all_tissues.tsv"
)

OUTDIR = (
    os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52", "05_TE_RNA", "05_DE", "04_hybrid_pattern")
)

df = pd.read_csv(
    INFILE,
    sep="\t"
)

# ============================================================
# Strict statistically supported transgressive events
# ============================================================

x = df[
    (df["Comparable_all3"] == 1) &
    (df["Transgressive"] == "YES")
].copy()

print("")
print("Strict transgressive tissue-PanTE events:")
print(len(x))

print("")
print("Unique transgressive PanTE:")
print(x["PanTE_ID"].nunique())


# ============================================================
# Per-PanTE recurrence
# ============================================================

rows = []

for pante, g in x.groupby("PanTE_ID"):

    tissues = sorted(
        g["Tissue"]
        .astype(str)
        .unique()
    )

    directions = sorted(
        g["Transgressive_direction"]
        .astype(str)
        .unique()
    )

    rows.append({
        "PanTE_ID":
            pante,

        "Broad_class":
            ";".join(
                sorted(
                    g["Broad_class"]
                    .dropna()
                    .astype(str)
                    .unique()
                )
            ),

        "N_tissues":
            len(tissues),

        "Tissues":
            ";".join(tissues),

        "Directions":
            ";".join(directions),

        "Direction_consistent":
            "YES"
            if len(directions) == 1
            else "NO"
    })


summary = pd.DataFrame(
    rows
)

summary.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_transgressive_cross_tissue.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Recurrence distribution
# ============================================================

rec = (
    summary[
        "N_tissues"
    ]
    .value_counts()
    .sort_index()
    .rename_axis(
        "N_tissues"
    )
    .reset_index(
        name="N_PanTE"
    )
)

rec[
    "Percent"
] = (
    rec["N_PanTE"]
    /
    summary["PanTE_ID"].nunique()
    * 100.0
)

rec.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_transgressive_cross_tissue_recurrence.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Exact tissue combinations for UpSet plot
# ============================================================

combo = (
    summary
    .groupby(
        [
            "Tissues",
            "Directions"
        ]
    )
    .size()
    .reset_index(
        name="N_PanTE"
    )
    .sort_values(
        "N_PanTE",
        ascending=False
    )
)

combo.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_transgressive_tissue_combinations.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Up and Down separately
# ============================================================

for direction in [
    "Up",
    "Down"
]:

    y = x[
        x[
            "Transgressive_direction"
        ] == direction
    ].copy()

    rows2 = []

    for pante, g in y.groupby(
        "PanTE_ID"
    ):

        tissues = sorted(
            g["Tissue"]
            .astype(str)
            .unique()
        )

        rows2.append({
            "PanTE_ID":
                pante,

            "N_tissues":
                len(tissues),

            "Tissues":
                ";".join(tissues)
        })


    z = pd.DataFrame(
        rows2
    )

    z.to_csv(
        os.path.join(
            OUTDIR,
            "PanTE_transgressive_{}_cross_tissue.tsv".format(
                direction.lower()
            )
        ),
        sep="\t",
        index=False
    )


    r = (
        z[
            "N_tissues"
        ]
        .value_counts()
        .sort_index()
        .rename_axis(
            "N_tissues"
        )
        .reset_index(
            name="N_PanTE"
        )
    )

    if len(z) > 0:

        r["Percent"] = (
            r["N_PanTE"]
            /
            len(z)
            * 100.0
        )

    r.to_csv(
        os.path.join(
            OUTDIR,
            "PanTE_transgressive_{}_recurrence.tsv".format(
                direction.lower()
            )
        ),
        sep="\t",
        index=False
    )


# ============================================================
# Direction consistency
# ============================================================

dc = (
    summary[
        "Direction_consistent"
    ]
    .value_counts()
)

print("")
print("Cross-tissue recurrence:")
print(
    rec.to_string(
        index=False
    )
)

print("")
print("Direction consistency:")
print(
    dc.to_string()
)

print("")
print("Top tissue combinations:")
print(
    combo.head(30).to_string(
        index=False
    )
)

print("")
print("Completed successfully.")