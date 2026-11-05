#!/usr/bin/env python
# Public-release note: original executed source; replace local absolute paths with paths from your environment before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import pandas as pd


BASE = os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52")
OUT = os.path.join(
    BASE,
    "05_TE_RNA/05_DE/05_genomic_RNA_integration"
)

INFILE = os.path.join(
    OUT,
    "PanTE_genomic_RNA_integrated_all_tissues.tsv"
)


df = pd.read_csv(
    INFILE,
    sep="\t"
)


# ============================================================
# Candidate:
#
# significant RNA transgressive
# but without corresponding genomic abundance change
# ============================================================

contexts = [
    "RNA_up_without_genomic_abundance_increase",
    "RNA_down_without_genomic_abundance_decrease"
]


cand = df[
    df[
        "Transgressive_genomic_context"
    ].isin(
        contexts
    )
].copy()


# ============================================================
# Save better-named full candidate file
# ============================================================

cand.to_csv(
    os.path.join(
        OUT,
        "PanTE_transgressive_without_corresponding_genomic_change.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Summarize each PanTE across tissues
# ============================================================

rows = []


for pante, x in cand.groupby(
    "PanTE_ID"
):

    tissues = sorted(
        x["Tissue"]
        .astype(str)
        .unique()
    )

    directions = sorted(
        x["Transgressive_direction"]
        .astype(str)
        .unique()
    )

    strong_n = int(
        (
            x[
                "Strong_transgressive"
            ].astype(str)
            ==
            "YES"
        ).sum()
    )


    # tissues with ATAC/WGBS available
    epi_tissues = sorted(
        set(tissues)
        &
        set([
            "Gill",
            "Kidney",
            "Muscle"
        ])
    )


    # direction consistency
    direction_consistent = (
        len(directions) == 1
    )


    # simple priority
    if (
        strong_n > 0
        and len(epi_tissues) > 0
        and len(tissues) >= 2
        and direction_consistent
    ):
        priority = "Tier1"

    elif (
        strong_n > 0
        and len(epi_tissues) > 0
    ):
        priority = "Tier2"

    elif (
        len(epi_tissues) > 0
    ):
        priority = "Tier3"

    else:
        priority = "Tier4"


    rows.append({
        "PanTE_ID":
            pante,

        "Broad_class":
            ";".join(
                sorted(
                    x[
                        "Broad_class"
                    ]
                    .dropna()
                    .astype(str)
                    .unique()
                )
            ),

        "Parent_genomic_bias":
            ";".join(
                sorted(
                    x[
                        "Parent_genomic_bias"
                    ]
                    .dropna()
                    .astype(str)
                    .unique()
                )
            ),

        "N_tissues":
            len(tissues),

        "Tissues":
            ";".join(tissues),

        "Direction":
            ";".join(directions),

        "Direction_consistent":
            "YES"
            if direction_consistent
            else "NO",

        "Strong_event_N":
            strong_n,

        "Epigenomic_tissue_N":
            len(epi_tissues),

        "Epigenomic_tissues":
            ";".join(
                epi_tissues
            ),

        "Priority":
            priority
    })


summary = pd.DataFrame(
    rows
)


priority_order = {
    "Tier1": 1,
    "Tier2": 2,
    "Tier3": 3,
    "Tier4": 4
}


summary["_priority"] = (
    summary[
        "Priority"
    ].map(
        priority_order
    )
)


summary = summary.sort_values(
    [
        "_priority",
        "Strong_event_N",
        "N_tissues"
    ],
    ascending=[
        True,
        False,
        False
    ]
)


summary = summary.drop(
    columns=[
        "_priority"
    ]
)


summary.to_csv(
    os.path.join(
        OUT,
        "PanTE_transgressive_candidate_priority.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Priority summary
# ============================================================

p = (
    summary[
        "Priority"
    ]
    .value_counts()
    .rename_axis(
        "Priority"
    )
    .reset_index(
        name="N_PanTE"
    )
)


p.to_csv(
    os.path.join(
        OUT,
        "PanTE_transgressive_candidate_priority_summary.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Recurrence summary
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


rec.to_csv(
    os.path.join(
        OUT,
        "PanTE_transgressive_recurrence_summary.tsv"
    ),
    sep="\t",
    index=False
)


print("")
print("Tissue-PanTE candidate records:")
print(len(cand))

print("")
print("Unique PanTE candidates:")
print(summary["PanTE_ID"].nunique())

print("")
print("Recurrence:")
print(rec.to_string(index=False))

print("")
print("Priority:")
print(p.to_string(index=False))

print("")
print("Top candidates:")
print(
    summary.head(30).to_string(
        index=False
    )
)

print("")
print("Priority analysis completed successfully.")