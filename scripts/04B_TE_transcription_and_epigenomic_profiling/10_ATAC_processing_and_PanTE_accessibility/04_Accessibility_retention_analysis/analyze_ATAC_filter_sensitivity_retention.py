#!/usr/bin/env python
# Public-release note: source script from the executed analysis. Replace the original local absolute paths and executable locations with paths in your environment before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import pandas as pd


BASE = os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52")
ATAC = os.path.join(BASE, "07_TE_ATAC")

INFILE = os.path.join(
    ATAC,
    "03_matrix/PanTE_ATAC_accessibility_long.tsv"
)

OUT = os.path.join(
    ATAC,
    "03_matrix/ATAC_filter_sensitivity_retention.tsv"
)


FILTERS = [
    ("Relaxed",   100,  1),
    ("Main",      500,  3),
    ("Stringent", 1000, 5)
]

STATES = [
    "BSB",
    "TC",
    "BTM",
    "BTP"
]

TISSUES = [
    "gill",
    "kidney",
    "muscle"
]


dat = pd.read_csv(
    INFILE,
    sep="\t"
)


results = []


for filter_name, min_bp, min_frag in FILTERS:

    print("")
    print("============================================")
    print(
        "{}: bp >= {}, each replicate fragments >= {}".format(
            filter_name,
            min_bp,
            min_frag
        )
    )
    print("============================================")


    # ========================================================
    # QUALIFY EACH STATE / TISSUE / PanTE
    # ========================================================

    qrows = []


    for (
        state,
        tissue,
        pante
    ), x in dat.groupby(
        [
            "Genome_state",
            "Tissue",
            "PanTE_ID"
        ]
    ):

        # Need exactly 3 biological replicates
        if len(x) != 3:
            continue


        bp = x[
            "Assignable_bp"
        ].iloc[0]


        if bp < min_bp:
            continue


        if (
            x[
                "Fragment_count"
            ]
            <
            min_frag
        ).any():

            continue


        qrows.append(
            {
                "Genome_state":
                    state,

                "Tissue":
                    tissue,

                "PanTE_ID":
                    pante,

                "Mean_log2FPKM1":
                    x[
                        "log2_Accessibility_FPKM1"
                    ].mean()
            }
        )


    qual = pd.DataFrame(
        qrows
    )


    # ========================================================
    # FOUR-STATE COMMON PanTE, BY TISSUE
    # ========================================================

    for tissue in TISSUES:

        sub = qual[
            qual["Tissue"] == tissue
        ]


        state_tables = {}


        for state in STATES:

            x = sub[
                sub["Genome_state"] == state
            ][
                [
                    "PanTE_ID",
                    "Mean_log2FPKM1"
                ]
            ].copy()


            x = x.rename(
                columns={
                    "Mean_log2FPKM1":
                    state
                }
            )


            state_tables[state] = x


        # Inner join = must pass filter in all 4 states
        mat = state_tables["BSB"]

        for state in [
            "TC",
            "BTM",
            "BTP"
        ]:

            mat = pd.merge(
                mat,
                state_tables[state],
                on="PanTE_ID",
                how="inner"
            )


        n = len(mat)


        # ====================================================
        # ABSOLUTE STATE RETENTION
        # ====================================================

        rho_bsb_btm = mat[
            "BSB"
        ].corr(
            mat["BTM"],
            method="spearman"
        )


        rho_tc_btp = mat[
            "TC"
        ].corr(
            mat["BTP"],
            method="spearman"
        )


        # ====================================================
        # PARENTAL CONTRAST RETENTION
        # ====================================================

        parent_delta = (
            mat["BSB"]
            -
            mat["TC"]
        )


        bt_delta = (
            mat["BTM"]
            -
            mat["BTP"]
        )


        rho_contrast = parent_delta.corr(
            bt_delta,
            method="spearman"
        )


        results.append(
            {
                "Filter":
                    filter_name,

                "Min_assignable_bp":
                    min_bp,

                "Min_fragments_each_rep":
                    min_frag,

                "Tissue":
                    tissue,

                "N_four_state_PanTE":
                    n,

                "Rho_BSB_BTM":
                    rho_bsb_btm,

                "Rho_TC_BTP":
                    rho_tc_btp,

                "Rho_parent_vs_BT_contrast":
                    rho_contrast
            }
        )


        print("")
        print(tissue)

        print(
            "N = {}".format(
                n
            )
        )

        print(
            "BSB -> BTM rho = {:.6f}".format(
                rho_bsb_btm
            )
        )

        print(
            "TC  -> BTP rho = {:.6f}".format(
                rho_tc_btp
            )
        )

        print(
            "Parent contrast -> BT contrast rho = {:.6f}".format(
                rho_contrast
            )
        )


res = pd.DataFrame(
    results
)


res.to_csv(
    OUT,
    sep="\t",
    index=False
)


print("")
print("============================================")
print("FINAL SUMMARY")
print("============================================")

print(
    res.to_string(
        index=False
    )
)


print("")
print("Saved:")
print(OUT)