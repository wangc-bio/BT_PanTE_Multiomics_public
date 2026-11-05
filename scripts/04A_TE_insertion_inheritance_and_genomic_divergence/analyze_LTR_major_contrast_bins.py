#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import numpy as np
import pandas as pd


KIM = os.environ.get("KIMURA_DIR", os.path.join(os.environ.get("KIMURA_BASE", "PATH_TO_RepeatMasker_unified_v5_2_Kimura_align"), "Kimura_divergence"))

ANA = os.path.join(
    KIM,
    "02_landscape_analysis"
)

INFILE = os.path.join(
    ANA,
    "Kimura_LTR_parental_contrast_by_bin.tsv"
)

OUT1 = os.path.join(
    ANA,
    "LTR_parental_contrast_ranked_bins.tsv"
)

OUT2 = os.path.join(
    ANA,
    "LTR_major_contrast_summary.tsv"
)


d = pd.read_csv(
    INFILE,
    sep="\t"
)

d = d.sort_values(
    "Kimura_bin"
).copy()


# ============================================================
# absolute parental difference
# ============================================================

d[
    "Abs_parent_delta"
] = np.abs(
    d[
        "Parent_delta_pct"
    ]
)


d[
    "Abs_BT_delta"
] = np.abs(
    d[
        "BT_delta_pct"
    ]
)


total_parent_abs = d[
    "Abs_parent_delta"
].sum()


# ============================================================
# sign concordance
# ============================================================

def sign_value(x):

    if x > 0:
        return 1

    if x < 0:
        return -1

    return 0


d[
    "Parent_sign"
] = d[
    "Parent_delta_pct"
].map(
    sign_value
)


d[
    "BT_sign"
] = d[
    "BT_delta_pct"
].map(
    sign_value
)


d[
    "Sign_concordant"
] = np.where(
    (
        d[
            "Parent_sign"
        ]
        ==
        d[
            "BT_sign"
        ]
    )
    &
    (
        d[
            "Parent_sign"
        ] != 0
    ),
    "YES",
    "NO"
)


# ============================================================
# rank by parental contrast magnitude
# ============================================================

ranked = d.sort_values(
    "Abs_parent_delta",
    ascending=False
).copy()


ranked[
    "Parent_contrast_contribution_pct"
] = (
    100.0
    *
    ranked[
        "Abs_parent_delta"
    ]
    /
    total_parent_abs
)


ranked[
    "Cumulative_parent_contrast_pct"
] = ranked[
    "Parent_contrast_contribution_pct"
].cumsum()


ranked[
    "Rank"
] = np.arange(
    1,
    len(
        ranked
    ) + 1
)


# ============================================================
# summary for top contributors
# ============================================================

rows = []


for cutoff in [
    50,
    80,
    90
]:

    x = ranked[
        ranked[
            "Cumulative_parent_contrast_pct"
        ] <= cutoff
    ].copy()


    # include first bin crossing the cutoff
    if len(x) < len(ranked):

        next_row = ranked.iloc[
            [
                len(x)
            ]
        ]

        x = pd.concat(
            [
                x,
                next_row
            ],
            ignore_index=True
        )


    n = len(
        x
    )


    sign_ok = int(
        (
            x[
                "Sign_concordant"
            ] == "YES"
        ).sum()
    )


    rows.append({
        "Cumulative_parent_contrast_target_pct":
            cutoff,

        "N_bins":
            n,

        "Kimura_bins":
            ",".join(
                map(
                    str,
                    x[
                        "Kimura_bin"
                    ].astype(
                        int
                    ).tolist()
                )
            ),

        "Observed_parent_contrast_pct":
            float(
                x[
                    "Parent_contrast_contribution_pct"
                ].sum()
            ),

        "Sign_concordant_bins":
            sign_ok,

        "Sign_concordance_fraction":
            (
                float(
                    sign_ok
                )
                /
                n
                if n > 0
                else np.nan
            ),

        "Parent_abs_delta_sum":
            float(
                x[
                    "Abs_parent_delta"
                ].sum()
            ),

        "BT_abs_delta_sum":
            float(
                x[
                    "Abs_BT_delta"
                ].sum()
            )
    })


summary = pd.DataFrame(
    rows
)


ranked.to_csv(
    OUT1,
    sep="\t",
    index=False
)


summary.to_csv(
    OUT2,
    sep="\t",
    index=False
)


print("")
print("============================================")
print("Top 15 LTR parental-contrast Kimura bins")
print("============================================")
print("")

print(
    ranked[
        [
            "Rank",
            "Kimura_bin",
            "BSB",
            "TC",
            "BTM",
            "BTP",
            "Parent_delta_pct",
            "BT_delta_pct",
            "Parent_contrast_contribution_pct",
            "Cumulative_parent_contrast_pct",
            "Sign_concordant"
        ]
    ]
    .head(
        15
    )
    .to_string(
        index=False
    )
)

print("")
print("============================================")
print("Cumulative contrast summary")
print("============================================")
print("")

print(
    summary.to_string(
        index=False
    )
)

print("")
print("Outputs:")
print(OUT1)
print(OUT2)
