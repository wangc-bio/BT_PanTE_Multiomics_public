#!/usr/bin/env python
# Public-release note: source script from the executed analysis. Replace the original local absolute paths and executable locations with paths in your environment before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import numpy as np
import pandas as pd


BASE = os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52")
ATAC = os.path.join(BASE, "07_TE_ATAC")
INDIR = os.path.join(ATAC, "03_matrix")

INFILE = os.path.join(
    INDIR,
    "PanTE_ATAC_accessibility_long.tsv"
)

OUTPREFIX = os.path.join(
    INDIR,
    "ATAC_parental_state_retention"
)


# ============================================================
# MAIN FILTER
# ============================================================

MIN_BP = 500
MIN_FRAG = 3


dat = pd.read_csv(
    INFILE,
    sep="\t"
)


states = [
    "BSB",
    "TC",
    "BTM",
    "BTP"
]

tissues = [
    "gill",
    "kidney",
    "muscle"
]


# ============================================================
# Require:
# 1. assignable bp >= 500
# 2. all 3 replicates >= 3 fragments
# ============================================================

qual_rows = []


for (state, tissue, pante), x in dat.groupby(
    [
        "Genome_state",
        "Tissue",
        "PanTE_ID"
    ]
):

    if len(x) != 3:
        continue

    bp = x["Assignable_bp"].iloc[0]

    if bp < MIN_BP:
        continue

    if (
        x["Fragment_count"] < MIN_FRAG
    ).any():
        continue

    qual_rows.append(
        {
            "Genome_state":
                state,

            "Tissue":
                tissue,

            "PanTE_ID":
                pante,

            "Assignable_bp":
                bp,

            "Mean_FPKM":
                x[
                    "Accessibility_FPKM"
                ].mean(),

            "Mean_log2FPKM1":
                x[
                    "log2_Accessibility_FPKM1"
                ].mean(),

            "Min_rep_fragments":
                x[
                    "Fragment_count"
                ].min(),

            "Mean_fragments":
                x[
                    "Fragment_count"
                ].mean()
        }
    )


qual = pd.DataFrame(
    qual_rows
)


qual.to_csv(
    OUTPREFIX + ".qualified_long.tsv",
    sep="\t",
    index=False
)


# ============================================================
# BUILD FOUR-STATE COMPARABLE MATRIX
# ============================================================

matrix_rows = []


for tissue in tissues:

    sub = qual[
        qual["Tissue"] == tissue
    ]

    sets = {}

    for state in states:

        sets[state] = set(
            sub.loc[
                sub["Genome_state"] == state,
                "PanTE_ID"
            ]
        )

    common = (
        sets["BSB"]
        &
        sets["TC"]
        &
        sets["BTM"]
        &
        sets["BTP"]
    )


    for pante in sorted(common):

        row = {
            "Tissue":
                tissue,

            "PanTE_ID":
                pante
        }


        for state in states:

            x = sub[
                (sub["Genome_state"] == state)
                &
                (sub["PanTE_ID"] == pante)
            ].iloc[0]

            row[
                state + "_FPKM"
            ] = x["Mean_FPKM"]

            row[
                state + "_log2FPKM1"
            ] = x["Mean_log2FPKM1"]

            row[
                state + "_Assignable_bp"
            ] = x["Assignable_bp"]


        # Parent -> hybrid change
        row[
            "BSB_to_BTM_delta"
        ] = (
            row["BTM_log2FPKM1"]
            -
            row["BSB_log2FPKM1"]
        )

        row[
            "TC_to_BTP_delta"
        ] = (
            row["BTP_log2FPKM1"]
            -
            row["TC_log2FPKM1"]
        )


        # Parent-specific contrast
        row[
            "Parent_delta"
        ] = (
            row["BSB_log2FPKM1"]
            -
            row["TC_log2FPKM1"]
        )


        # Hybrid haplotype contrast
        row[
            "BT_delta"
        ] = (
            row["BTM_log2FPKM1"]
            -
            row["BTP_log2FPKM1"]
        )


        matrix_rows.append(
            row
        )


mat = pd.DataFrame(
    matrix_rows
)


mat.to_csv(
    OUTPREFIX + ".four_state_matrix.tsv",
    sep="\t",
    index=False
)


# ============================================================
# ABSOLUTE STATE RETENTION
# ============================================================

retention_rows = []


for tissue in tissues:

    x = mat[
        mat["Tissue"] == tissue
    ].copy()


    comparisons = [
        (
            "BSB_to_BTM",
            "BSB_log2FPKM1",
            "BTM_log2FPKM1"
        ),

        (
            "TC_to_BTP",
            "TC_log2FPKM1",
            "BTP_log2FPKM1"
        )
    ]


    for name, a, b in comparisons:

        rho = x[a].corr(
            x[b],
            method="spearman"
        )

        delta = (
            x[b] - x[a]
        )


        retention_rows.append(
            {
                "Tissue":
                    tissue,

                "Comparison":
                    name,

                "N_PanTE":
                    len(x),

                "Spearman_rho":
                    rho,

                "Median_delta_log2FPKM1":
                    delta.median(),

                "Mean_delta_log2FPKM1":
                    delta.mean()
            }
        )


retention = pd.DataFrame(
    retention_rows
)


retention.to_csv(
    OUTPREFIX + ".absolute_state_summary.tsv",
    sep="\t",
    index=False
)


# ============================================================
# PARENTAL CONTRAST RETENTION
# ============================================================

contrast_rows = []


for tissue in tissues:

    x = mat[
        mat["Tissue"] == tissue
    ].copy()


    rho = x[
        "Parent_delta"
    ].corr(
        x["BT_delta"],
        method="spearman"
    )


    contrast_rows.append(
        {
            "Tissue":
                tissue,

            "N_PanTE":
                len(x),

            "Spearman_rho":
                rho,

            "Median_parent_delta":
                x[
                    "Parent_delta"
                ].median(),

            "Median_BT_delta":
                x[
                    "BT_delta"
                ].median()
        }
    )


contrast = pd.DataFrame(
    contrast_rows
)


contrast.to_csv(
    OUTPREFIX + ".parental_contrast_summary.tsv",
    sep="\t",
    index=False
)


# ============================================================
# RETAINED / ATTENUATED / REVERSED
#
# Sensitivity:
# parental contrast threshold = 0.5 / 1.0 / 1.5
#
# This is an effect-size framework, NOT statistical significance.
# ============================================================

class_rows = []


for cutoff in [
    0.5,
    1.0,
    1.5
]:

    for tissue in tissues:

        x = mat[
            mat["Tissue"] == tissue
        ].copy()


        x = x[
            x["Parent_delta"].abs()
            >=
            cutoff
        ].copy()


        n = len(x)


        if n == 0:
            continue


        same_direction = (
            np.sign(
                x["Parent_delta"]
            )
            ==
            np.sign(
                x["BT_delta"]
            )
        )


        strong_retained = (
            same_direction
            &
            (
                x["BT_delta"].abs()
                >=
                cutoff
            )
        )


        attenuated = (
            same_direction
            &
            (
                x["BT_delta"].abs()
                <
                cutoff
            )
        )


        reversed_state = (
            ~same_direction
        )


        categories = {
            "Retained_magnitude":
                strong_retained,

            "Attenuated":
                attenuated,

            "Reversed":
                reversed_state
        }


        for category, mask in categories.items():

            count = int(
                mask.sum()
            )


            class_rows.append(
                {
                    "Contrast_cutoff":
                        cutoff,

                    "Tissue":
                        tissue,

                    "Category":
                        category,

                    "N_parent_differential":
                        n,

                    "N":
                        count,

                    "Percent":
                        100.0
                        *
                        count
                        /
                        n
                }
            )


classes = pd.DataFrame(
    class_rows
)


classes.to_csv(
    OUTPREFIX + ".contrast_classification.tsv",
    sep="\t",
    index=False
)


# ============================================================
# PRINT SUMMARY
# ============================================================

print("")
print("============================================")
print("FILTER")
print("============================================")

print(
    "Assignable_bp >= {}".format(
        MIN_BP
    )
)

print(
    "Each of 3 replicates fragments >= {}".format(
        MIN_FRAG
    )
)


print("")
print("============================================")
print("FOUR-STATE COMPARABLE PanTE")
print("============================================")

print(
    mat.groupby(
        "Tissue"
    ).size()
)


print("")
print("============================================")
print("ABSOLUTE STATE RETENTION")
print("============================================")

print(
    retention.to_string(
        index=False
    )
)


print("")
print("============================================")
print("PARENTAL CONTRAST RETENTION")
print("============================================")

print(
    contrast.to_string(
        index=False
    )
)


print("")
print("============================================")
print("CONTRAST CLASSIFICATION")
print("============================================")

print(
    classes.to_string(
        index=False
    )
)


print("")
print("Completed successfully.")