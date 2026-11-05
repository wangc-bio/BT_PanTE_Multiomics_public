#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import numpy as np
import pandas as pd

from scipy.stats import spearmanr


# ============================================================
# PATHS
# ============================================================

INTEGRATION_DIR = os.environ.get("PANTE_INTEGRATION_DIR", os.path.join(os.getcwd(), "integration"))

MASTER = os.environ.get("PANTE_MASTER_MATRIX", os.path.join(INTEGRATION_DIR, "01_matrix", "PanTE_multiomics_master.RNA_testable.tsv"))

OUTDIR = os.environ.get("PANTE_OUTPUT_DIR", os.path.join(os.getcwd(), "results", "05_genomic_divergence_vs_epigenetic_deviation"))

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# SETTINGS
# ============================================================

TISSUES = [
    "Gill",
    "Kidney",
    "Muscle"
]

ASSAYS = {
    "mCG": {
        "metric": "D_mCG_mean",
        "available": "WGBS_available"
    },

    "ATAC": {
        "metric": "D_ATAC_mean",
        "available": "ATAC_available"
    }
}

BACKGROUNDS = [
    "All_assay_available",
    "Parent_informative"
]


# ============================================================
# FUNCTIONS
# ============================================================

def bh_adjust(pvalues):

    p = np.asarray(
        pvalues,
        dtype=float
    )

    result = np.full(
        len(p),
        np.nan,
        dtype=float
    )

    valid = np.isfinite(
        p
    )

    if valid.sum() == 0:
        return result


    pv = p[
        valid
    ]

    order = np.argsort(
        pv
    )

    ranked = pv[
        order
    ]

    n = len(
        ranked
    )

    adjusted_ranked = np.empty(
        n,
        dtype=float
    )

    running = 1.0


    for i in range(
        n - 1,
        -1,
        -1
    ):

        rank = i + 1

        q = (
            ranked[i]
            *
            n
            /
            float(rank)
        )

        running = min(
            running,
            q
        )

        adjusted_ranked[
            i
        ] = min(
            running,
            1.0
        )


    tmp = np.empty(
        n,
        dtype=float
    )

    tmp[
        order
    ] = adjusted_ranked


    result[
        np.where(
            valid
        )[0]
    ] = tmp


    return result


def available_mask(series):

    """
    Robust handling of availability columns:
    1/0, True/False, YES/NO, etc.
    """

    if pd.api.types.is_bool_dtype(
        series
    ):

        return series.fillna(
            False
        )


    if pd.api.types.is_numeric_dtype(
        series
    ):

        return (
            series.fillna(
                0
            )
            >
            0
        )


    s = (
        series
        .astype(str)
        .str.strip()
        .str.lower()
    )


    return s.isin(
        [
            "1",
            "true",
            "yes",
            "y"
        ]
    )


# ============================================================
# LOAD
# ============================================================

if not os.path.isfile(
    MASTER
):

    raise RuntimeError(
        "Missing master matrix: {}".format(
            MASTER
        )
    )


df = pd.read_csv(
    MASTER,
    sep="\t"
)


required = [
    "PanTE_ID",
    "Tissue",
    "D_genomic",
    "Parent_max_bp",
    "D_mCG_mean",
    "D_ATAC_mean",
    "WGBS_available",
    "ATAC_available"
]


for col in required:

    if col not in df.columns:

        raise RuntimeError(
            "Missing required column: {}".format(
                col
            )
        )


# ============================================================
# BASIC QC
# ============================================================

df = df[
    df[
        "Tissue"
    ].isin(
        TISSUES
    )
].copy()


if df.duplicated(
    subset=[
        "PanTE_ID",
        "Tissue"
    ]
).any():

    raise RuntimeError(
        "Duplicated PanTE_ID x Tissue detected."
    )


# Confirm D_genomic range
finite_D = df.loc[
    df[
        "D_genomic"
    ].notna(),
    "D_genomic"
].astype(float)


if len(
    finite_D
) == 0:

    raise RuntimeError(
        "No finite D_genomic values."
    )


if (
    finite_D.min() < -1e-12
    or
    finite_D.max() > 1.0 + 1e-12
):

    raise RuntimeError(
        "D_genomic outside expected range [0,1]: min={}, max={}".format(
            finite_D.min(),
            finite_D.max()
        )
    )


print("")
print("==============================================")
print("D_genomic QC")
print("==============================================")
print("")

print(
    "N finite = {:,}".format(
        len(
            finite_D
        )
    )
)

print(
    "Range = {:.6f} - {:.6f}".format(
        finite_D.min(),
        finite_D.max()
    )
)

print(
    "Median = {:.6f}".format(
        finite_D.median()
    )
)


# ============================================================
# CORRELATION ANALYSIS
# ============================================================

rows = []
used_rows = []


for background in BACKGROUNDS:

    for assay in [
        "mCG",
        "ATAC"
    ]:

        metric = ASSAYS[
            assay
        ][
            "metric"
        ]

        available_col = ASSAYS[
            assay
        ][
            "available"
        ]


        for tissue in TISSUES:

            x = df[
                df[
                    "Tissue"
                ] == tissue
            ].copy()


            # ------------------------------------------------
            # assay-specific availability
            # ------------------------------------------------

            x = x[
                available_mask(
                    x[
                        available_col
                    ]
                )
            ].copy()


            # ------------------------------------------------
            # Parent-informative sensitivity
            #
            # Same occupancy criterion as Result 3.3:
            # max parental raw occupancy >=10 kb
            # ------------------------------------------------

            if background == "Parent_informative":

                x = x[
                    x[
                        "Parent_max_bp"
                    ]
                    >=
                    10000
                ].copy()


            # ------------------------------------------------
            # Complete pair for correlation
            # ------------------------------------------------

            x = x[
                x[
                    "D_genomic"
                ].notna()
                &
                x[
                    metric
                ].notna()
            ].copy()


            x[
                "D_genomic"
            ] = x[
                "D_genomic"
            ].astype(
                float
            )


            x[
                metric
            ] = x[
                metric
            ].astype(
                float
            )


            n = len(
                x
            )


            if n < 3:

                rho = np.nan
                p = np.nan

            else:

                rho, p = spearmanr(
                    x[
                        "D_genomic"
                    ].values,
                    x[
                        metric
                    ].values
                )

                rho = float(
                    rho
                )

                p = float(
                    p
                )


            rows.append({

                "Background":
                    background,

                "Assay":
                    assay,

                "Metric":
                    metric,

                "Tissue":
                    tissue,

                "N":
                    n,

                "Median_D_genomic":
                    float(
                        x[
                            "D_genomic"
                        ].median()
                    )
                    if n
                    else np.nan,

                "Q25_D_genomic":
                    float(
                        x[
                            "D_genomic"
                        ].quantile(
                            0.25
                        )
                    )
                    if n
                    else np.nan,

                "Q75_D_genomic":
                    float(
                        x[
                            "D_genomic"
                        ].quantile(
                            0.75
                        )
                    )
                    if n
                    else np.nan,

                "Median_epigenetic_deviation":
                    float(
                        x[
                            metric
                        ].median()
                    )
                    if n
                    else np.nan,

                "Spearman_rho":
                    rho,

                "P":
                    p
            })


            # save row-level data for plotting
            temp = x[
                [
                    "PanTE_ID",
                    "Tissue",
                    "D_genomic",
                    metric
                ]
            ].copy()


            temp[
                "Background"
            ] = background


            temp[
                "Assay"
            ] = assay


            temp[
                "Metric"
            ] = metric


            temp = temp.rename(
                columns={
                    metric:
                        "Epigenetic_deviation"
                }
            )


            used_rows.append(
                temp
            )


res = pd.DataFrame(
    rows
)


# ============================================================
# MULTIPLE TESTING
#
# 12 total tests:
# 2 backgrounds x 2 assays x 3 tissues
# ============================================================

res[
    "FDR_all12"
] = bh_adjust(
    res[
        "P"
    ].values
)


# 6 tests separately within each background
res[
    "FDR_within_background"
] = np.nan


for background in BACKGROUNDS:

    idx = (
        res[
            "Background"
        ]
        ==
        background
    )


    res.loc[
        idx,
        "FDR_within_background"
    ] = bh_adjust(
        res.loc[
            idx,
            "P"
        ].values
    )


# 3 tissues separately within each background + assay
res[
    "FDR_within_background_assay"
] = np.nan


for background in BACKGROUNDS:

    for assay in [
        "mCG",
        "ATAC"
    ]:

        idx = (
            (
                res[
                    "Background"
                ]
                ==
                background
            )
            &
            (
                res[
                    "Assay"
                ]
                ==
                assay
            )
        )


        res.loc[
            idx,
            "FDR_within_background_assay"
        ] = bh_adjust(
            res.loc[
                idx,
                "P"
            ].values
        )


# ============================================================
# DECILE SUMMARIES FOR POSSIBLE FIGURE
#
# This is descriptive only, not inferential.
# ============================================================

used = pd.concat(
    used_rows,
    ignore_index=True
)


decile_rows = []


for (
    background,
    assay,
    tissue
), x in used.groupby(
    [
        "Background",
        "Assay",
        "Tissue"
    ]
):

    x = x.copy()


    if len(x) < 20:
        continue


    # rank first so ties in D_genomic do not make qcut fail
    ranks = x[
        "D_genomic"
    ].rank(
        method="first"
    )


    x[
        "D_genomic_decile"
    ] = pd.qcut(
        ranks,
        q=10,
        labels=False
    ) + 1


    for decile, y in x.groupby(
        "D_genomic_decile"
    ):

        decile_rows.append({

            "Background":
                background,

            "Assay":
                assay,

            "Tissue":
                tissue,

            "D_genomic_decile":
                int(
                    decile
                ),

            "N":
                len(
                    y
                ),

            "Median_D_genomic":
                float(
                    y[
                        "D_genomic"
                    ].median()
                ),

            "Median_epigenetic_deviation":
                float(
                    y[
                        "Epigenetic_deviation"
                    ].median()
                ),

            "Q25_epigenetic_deviation":
                float(
                    y[
                        "Epigenetic_deviation"
                    ].quantile(
                        0.25
                    )
                ),

            "Q75_epigenetic_deviation":
                float(
                    y[
                        "Epigenetic_deviation"
                    ].quantile(
                        0.75
                    )
                )
        })


deciles = pd.DataFrame(
    decile_rows
)


# ============================================================
# OUTPUT
# ============================================================

RESULT_OUT = os.path.join(
    OUTDIR,
    "Dgenomic_vs_epigenetic_deviation.Spearman.tsv"
)

DATA_OUT = os.path.join(
    OUTDIR,
    "Dgenomic_vs_epigenetic_deviation.data.tsv"
)

DECILE_OUT = os.path.join(
    OUTDIR,
    "Dgenomic_vs_epigenetic_deviation.decile_summary.tsv"
)

SUMMARY_OUT = os.path.join(
    OUTDIR,
    "Dgenomic_vs_epigenetic_deviation.summary.txt"
)


res.to_csv(
    RESULT_OUT,
    sep="\t",
    index=False
)


used.to_csv(
    DATA_OUT,
    sep="\t",
    index=False,
    na_rep="NA"
)


deciles.to_csv(
    DECILE_OUT,
    sep="\t",
    index=False
)


with open(
    SUMMARY_OUT,
    "w"
) as out:

    out.write(
        "Parental genomic divergence vs matched-parent epigenetic deviation\n"
    )

    out.write(
        "=================================================================\n\n"
    )

    out.write(
        "D_genomic = |BSB_bp_per_Mb - TC_bp_per_Mb| / "
        "(BSB_bp_per_Mb + TC_bp_per_Mb)\n\n"
    )

    out.write(
        "mCG deviation = D_mCG_mean\n"
    )

    out.write(
        "ATAC deviation = D_ATAC_mean\n\n"
    )

    out.write(
        "Backgrounds:\n"
    )

    out.write(
        "1. All assay-available PanTE x tissue observations\n"
    )

    out.write(
        "2. Parent-informative subset: Parent_max_bp >= 10,000\n\n"
    )

    out.write(
        res.to_string(
            index=False
        )
    )

    out.write(
        "\n"
    )


# ============================================================
# SCREEN
# ============================================================

print("")
print("==============================================")
print("D_genomic vs epigenetic deviation")
print("==============================================")
print("")

print(
    res.to_string(
        index=False
    )
)

print("")
print("Output directory:")
print(OUTDIR)
