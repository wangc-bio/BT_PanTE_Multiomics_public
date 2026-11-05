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
    "Kimura_landscape_main5.genome_normalized.wide.tsv"
)

OUTFILE = os.path.join(
    ANA,
    "Kimura_nonzero_bin_sensitivity.tsv"
)

DETAIL_OUT = os.path.join(
    ANA,
    "Kimura_nonzero_bin_detail.tsv"
)


CLASSES = [
    "DNA",
    "LTR",
    "LINE",
    "SINE",
    "RC"
]


def spearman(x, y):

    x = pd.Series(
        x,
        dtype=float
    )

    y = pd.Series(
        y,
        dtype=float
    )

    valid = (
        x.notna()
        &
        y.notna()
    )

    x = x[
        valid
    ]

    y = y[
        valid
    ]

    if len(x) < 3:
        return np.nan

    rx = x.rank(
        method="average"
    )

    ry = y.rank(
        method="average"
    )

    if (
        rx.std() == 0
        or
        ry.std() == 0
    ):
        return np.nan

    return float(
        np.corrcoef(
            rx,
            ry
        )[0, 1]
    )


def pearson(x, y):

    x = np.asarray(
        x,
        dtype=float
    )

    y = np.asarray(
        y,
        dtype=float
    )

    valid = (
        np.isfinite(x)
        &
        np.isfinite(y)
    )

    x = x[
        valid
    ]

    y = y[
        valid
    ]

    if len(x) < 3:
        return np.nan

    if (
        np.std(x) == 0
        or
        np.std(y) == 0
    ):
        return np.nan

    return float(
        np.corrcoef(
            x,
            y
        )[0, 1]
    )


d = pd.read_csv(
    INFILE,
    sep="\t"
)


results = []
details = []


def analyze_pair(
    cls,
    comparison,
    x,
    y,
    bins
):

    x = np.asarray(
        x,
        dtype=float
    )

    y = np.asarray(
        y,
        dtype=float
    )

    bins = np.asarray(
        bins
    )

    valid = (
        np.isfinite(x)
        &
        np.isfinite(y)
    )

    x = x[
        valid
    ]

    y = y[
        valid
    ]

    bins = bins[
        valid
    ]


    # ----------------------------------------------
    # full correlation
    # ----------------------------------------------

    rho_full = spearman(
        x,
        y
    )

    r_full = pearson(
        x,
        y
    )


    # ----------------------------------------------
    # sensitivity:
    # remove bins where BOTH sides equal zero
    # ----------------------------------------------

    keep = (
        (x > 0)
        |
        (y > 0)
    )


    x2 = x[
        keep
    ]

    y2 = y[
        keep
    ]

    bins2 = bins[
        keep
    ]


    rho_nonzero = spearman(
        x2,
        y2
    )

    r_nonzero = pearson(
        x2,
        y2
    )


    results.append({
        "Broad_class":
            cls,

        "Comparison":
            comparison,

        "N_bins_full":
            len(x),

        "N_bins_nonzero":
            len(x2),

        "N_both_zero_removed":
            len(x) - len(x2),

        "Spearman_full":
            rho_full,

        "Spearman_nonzero":
            rho_nonzero,

        "Delta_Spearman":
            (
                rho_nonzero
                -
                rho_full
            ),

        "Pearson_full":
            r_full,

        "Pearson_nonzero":
            r_nonzero,

        "Delta_Pearson":
            (
                r_nonzero
                -
                r_full
            )
    })


    for b, xv, yv in zip(
        bins2,
        x2,
        y2
    ):

        details.append({
            "Broad_class":
                cls,

            "Comparison":
                comparison,

            "Kimura_bin":
                int(
                    b
                ),

            "X":
                xv,

            "Y":
                yv
        })


# ============================================================
# CLASS-SPECIFIC
# ============================================================

for cls in CLASSES:

    x = d[
        d[
            "Broad_class"
        ] == cls
    ].copy()

    x = x.sort_values(
        "Kimura_bin"
    )


    # BSB -> BT-M
    analyze_pair(
        cls,
        "BSB_vs_BTM",
        x[
            "BSB"
        ],
        x[
            "BTM"
        ],
        x[
            "Kimura_bin"
        ]
    )


    # TC -> BT-P
    analyze_pair(
        cls,
        "TC_vs_BTP",
        x[
            "TC"
        ],
        x[
            "BTP"
        ],
        x[
            "Kimura_bin"
        ]
    )


    # parental contrast -> BT contrast
    #
    # Here zero means the contrast itself is exactly zero.
    analyze_pair(
        cls,
        "ParentContrast_vs_BTContrast",
        x[
            "Parent_delta_pct"
        ],
        x[
            "BT_delta_pct"
        ],
        x[
            "Kimura_bin"
        ]
    )


# ============================================================
# ALL MAIN5
#
# First sum five classes within each Kimura bin.
# ============================================================

all5 = (
    d
    .groupby(
        "Kimura_bin",
        as_index=False
    )[
        [
            "BSB",
            "TC",
            "BTM",
            "BTP"
        ]
    ]
    .sum()
)


all5[
    "Parent_delta_pct"
] = (
    all5[
        "BSB"
    ]
    -
    all5[
        "TC"
    ]
)


all5[
    "BT_delta_pct"
] = (
    all5[
        "BTM"
    ]
    -
    all5[
        "BTP"
    ]
)


analyze_pair(
    "All_main5",
    "BSB_vs_BTM",
    all5[
        "BSB"
    ],
    all5[
        "BTM"
    ],
    all5[
        "Kimura_bin"
    ]
)


analyze_pair(
    "All_main5",
    "TC_vs_BTP",
    all5[
        "TC"
    ],
    all5[
        "BTP"
    ],
    all5[
        "Kimura_bin"
    ]
)


analyze_pair(
    "All_main5",
    "ParentContrast_vs_BTContrast",
    all5[
        "Parent_delta_pct"
    ],
    all5[
        "BT_delta_pct"
    ],
    all5[
        "Kimura_bin"
    ]
)


# ============================================================
# OUTPUT
# ============================================================

res = pd.DataFrame(
    results
)

detail = pd.DataFrame(
    details
)


res.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)


detail.to_csv(
    DETAIL_OUT,
    sep="\t",
    index=False
)


print("")
print("============================================")
print("Non-zero-bin sensitivity completed")
print("============================================")
print("")

print(
    res.to_string(
        index=False
    )
)

print("")
print("Output:")
print(
    OUTFILE
)
