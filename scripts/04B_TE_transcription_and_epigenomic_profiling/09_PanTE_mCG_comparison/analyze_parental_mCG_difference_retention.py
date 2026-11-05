#!/usr/bin/env python
# Public-release note: original executed source; replace local absolute paths with local paths before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import pandas as pd
import numpy as np
from scipy.stats import spearmanr


BASE = os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52")
OUTDIR = os.path.join(
    BASE,
    "06_TE_WGBS/03_matrix"
)

TISSUES = [
    "gill",
    "kidney",
    "muscle"
]


summary_rows = []
sensitivity_rows = []


for tissue in TISSUES:

    # --------------------------------------------------------
    # Main-filter matrices
    # --------------------------------------------------------

    f_bsb = os.path.join(
        OUTDIR,
        "main.{}.BSB_to_BTM.PanTE_mCG.tsv".format(
            tissue
        )
    )

    f_tc = os.path.join(
        OUTDIR,
        "main.{}.TC_to_BTP.PanTE_mCG.tsv".format(
            tissue
        )
    )


    a = pd.read_csv(
        f_bsb,
        sep="\t"
    )

    b = pd.read_csv(
        f_tc,
        sep="\t"
    )


    # Only retain columns needed here
    a = a[
        [
            "PanTE_ID",
            "BSB_pooled_mCG",
            "BTM_pooled_mCG"
        ]
    ].copy()


    b = b[
        [
            "PanTE_ID",
            "TC_pooled_mCG",
            "BTP_pooled_mCG"
        ]
    ].copy()


    x = pd.merge(
        a,
        b,
        on="PanTE_ID",
        how="inner"
    )


    # --------------------------------------------------------
    # Parental contrast and hybrid haplotype contrast
    #
    # positive = BSB / BTM more methylated
    # negative = TC / BTP more methylated
    # --------------------------------------------------------

    x["Parent_delta_mCG"] = (
        x["BSB_pooled_mCG"]
        -
        x["TC_pooled_mCG"]
    )


    x["BT_haplotype_delta_mCG"] = (
        x["BTM_pooled_mCG"]
        -
        x["BTP_pooled_mCG"]
    )


    x["Retention_residual"] = (
        x["BT_haplotype_delta_mCG"]
        -
        x["Parent_delta_mCG"]
    )


    x["Abs_parent_delta_mCG"] = np.abs(
        x["Parent_delta_mCG"]
    )


    x["Abs_BT_delta_mCG"] = np.abs(
        x["BT_haplotype_delta_mCG"]
    )


    # --------------------------------------------------------
    # Overall correlation
    # --------------------------------------------------------

    rho, pvalue = spearmanr(
        x["Parent_delta_mCG"],
        x["BT_haplotype_delta_mCG"]
    )


    summary_rows.append({
        "Tissue":
            tissue,

        "N_all_four_state_PanTE":
            len(x),

        "Spearman_rho_parent_vs_BT_delta":
            rho,

        "Spearman_p":
            pvalue,

        "Median_parent_delta":
            x[
                "Parent_delta_mCG"
            ].median(),

        "Median_BT_delta":
            x[
                "BT_haplotype_delta_mCG"
            ].median(),

        "Median_retention_residual":
            x[
                "Retention_residual"
            ].median()
    })


    # --------------------------------------------------------
    # Sensitivity:
    # define a parent-differential subset by effect size
    # --------------------------------------------------------

    for cutoff in [
        0.05,
        0.10,
        0.15
    ]:

        z = x[
            x[
                "Abs_parent_delta_mCG"
            ]
            >=
            cutoff
        ].copy()


        if len(z) == 0:
            continue


        z[
            "Parent_direction"
        ] = np.where(
            z[
                "Parent_delta_mCG"
            ] > 0,
            "BSB_higher",
            "TC_higher"
        )


        z[
            "BT_direction"
        ] = np.where(
            z[
                "BT_haplotype_delta_mCG"
            ] > 0,
            "BTM_higher",
            "BTP_higher"
        )


        # Direction retention:
        #
        # parent BSB higher -> BTM higher
        # parent TC higher  -> BTP higher

        z[
            "Direction_retained"
        ] = (
            (
                (
                    z[
                        "Parent_delta_mCG"
                    ] > 0
                )
                &
                (
                    z[
                        "BT_haplotype_delta_mCG"
                    ] > 0
                )
            )
            |
            (
                (
                    z[
                        "Parent_delta_mCG"
                    ] < 0
                )
                &
                (
                    z[
                        "BT_haplotype_delta_mCG"
                    ] < 0
                )
            )
        )


        retained = int(
            z[
                "Direction_retained"
            ].sum()
        )

        reversed_n = (
            len(z)
            -
            retained
        )


        # Strong preservation:
        # same direction AND BT haplotype difference
        # still reaches the same cutoff
        strong_retained = (
            z[
                "Direction_retained"
            ]
            &
            (
                z[
                    "Abs_BT_delta_mCG"
                ]
                >=
                cutoff
            )
        ).sum()


        sensitivity_rows.append({
            "Tissue":
                tissue,

            "Parent_difference_cutoff":
                cutoff,

            "Parent_differential_N":
                len(z),

            "Direction_retained_N":
                retained,

            "Direction_reversed_N":
                reversed_n,

            "Direction_retained_percent":
                100.0
                *
                retained
                /
                len(z),

            "Strong_retained_N":
                int(
                    strong_retained
                ),

            "Strong_retained_percent":
                100.0
                *
                strong_retained
                /
                len(z)
        })


    outfile = os.path.join(
        OUTDIR,
        "{}.parental_mCG_difference_retention.tsv".format(
            tissue
        )
    )

    x.to_csv(
        outfile,
        sep="\t",
        index=False
    )


# ============================================================
# Save summaries
# ============================================================

summary = pd.DataFrame(
    summary_rows
)

sensitivity = pd.DataFrame(
    sensitivity_rows
)


summary.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_parental_mCG_difference_retention_summary.tsv"
    ),
    sep="\t",
    index=False
)


sensitivity.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_parental_mCG_direction_retention_sensitivity.tsv"
    ),
    sep="\t",
    index=False
)


print("")
print(
    "============================================"
)
print(
    "PARENTAL mCG DIFFERENCE RETENTION"
)
print(
    "============================================"
)

print("")
print(
    summary.to_string(
        index=False
    )
)


print("")
print(
    "============================================"
)
print(
    "DIRECTION RETENTION SENSITIVITY"
)
print(
    "============================================"
)

print("")
print(
    sensitivity.to_string(
        index=False
    )
)


print("")
print(
    "Completed successfully."
)