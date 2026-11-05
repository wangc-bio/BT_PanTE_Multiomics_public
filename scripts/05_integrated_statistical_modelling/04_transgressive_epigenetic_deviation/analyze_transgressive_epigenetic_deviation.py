#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu


# ============================================================
# PATHS
# ============================================================

INTEGRATION_DIR = os.environ.get("PANTE_INTEGRATION_DIR", os.path.join(os.getcwd(), "integration"))

MASTER = os.environ.get("PANTE_MASTER_MATRIX", os.path.join(INTEGRATION_DIR, "01_matrix", "PanTE_multiomics_master.RNA_testable.tsv"))

OUTDIR = os.environ.get("PANTE_OUTPUT_DIR", os.path.join(os.getcwd(), "results", "02_transgressive_vs_epigenetic_deviation"))

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# SETTINGS
# ============================================================

TISSUES = [
    "Gill",
    "Kidney",
    "Muscle"
]

METRICS = [
    ("D_mCG_mean", "mCG"),
    ("D_ATAC_mean", "ATAC")
]


# ============================================================
# FUNCTIONS
# ============================================================

def bh_adjust(pvalues):

    p = np.asarray(
        pvalues,
        dtype=float
    )

    n = len(p)

    order = np.argsort(
        p
    )

    ranked = p[
        order
    ]

    adjusted = np.empty(
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

        value = (
            ranked[i]
            *
            n
            /
            rank
        )

        running = min(
            running,
            value
        )

        adjusted[
            order[i]
        ] = min(
            running,
            1.0
        )

    return adjusted


def cliffs_delta_from_u(
    u,
    n1,
    n0
):

    """
    Cliff's delta calculated from Mann-Whitney U.

    Group 1 = strict transgressive
    Group 0 = non-transgressive

    Positive delta:
        transgressive group tends to have larger deviation.

    Negative delta:
        transgressive group tends to have smaller deviation.
    """

    if n1 == 0 or n0 == 0:
        return np.nan

    return (
        2.0
        *
        float(u)
        /
        (
            float(n1)
            *
            float(n0)
        )
        -
        1.0
    )


def effect_label(delta):

    if pd.isna(delta):
        return "NA"

    x = abs(
        delta
    )

    if x < 0.147:
        return "negligible"

    if x < 0.330:
        return "small"

    if x < 0.474:
        return "medium"

    return "large"


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
    "RNA_testable",
    "Strict_transgressive",
    "D_mCG_mean",
    "D_ATAC_mean"
]


for col in required:

    if col not in df.columns:

        raise RuntimeError(
            "Missing required column: {}".format(
                col
            )
        )


# ============================================================
# BASIC FILTER
# ============================================================

df = df[
    df[
        "Tissue"
    ].isin(
        TISSUES
    )
].copy()


df = df[
    df[
        "RNA_testable"
    ] == 1
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


if not df[
    "Strict_transgressive"
].isin(
    [
        "YES",
        "NO"
    ]
).all():

    raise RuntimeError(
        "Unexpected Strict_transgressive values."
    )


# ============================================================
# QC
# ============================================================

print("")
print("==============================================")
print("RNA-testable input")
print("==============================================")
print("")

print(
    "Rows: {:,}".format(
        len(df)
    )
)

print(
    "Unique PanTE IDs: {:,}".format(
        df[
            "PanTE_ID"
        ].nunique()
    )
)

print("")

for tissue in TISSUES:

    x = df[
        df[
            "Tissue"
        ] == tissue
    ]

    print(
        "{}: rows={:,}, strict_transgressive={:,}".format(
            tissue,
            len(x),
            int(
                (
                    x[
                        "Strict_transgressive"
                    ] == "YES"
                ).sum()
            )
        )
    )


# ============================================================
# STATISTICS
# ============================================================

results = []
group_summary = []


for metric, assay in METRICS:

    for tissue in TISSUES:

        x = df[
            df[
                "Tissue"
            ] == tissue
        ][
            [
                "PanTE_ID",
                "Strict_transgressive",
                metric
            ]
        ].copy()


        # metric-specific complete observations
        x = x[
            x[
                metric
            ].notna()
        ].copy()


        trans = x[
            x[
                "Strict_transgressive"
            ] == "YES"
        ][
            metric
        ].astype(
            float
        ).values


        nontrans = x[
            x[
                "Strict_transgressive"
            ] == "NO"
        ][
            metric
        ].astype(
            float
        ).values


        n_trans = len(
            trans
        )

        n_nontrans = len(
            nontrans
        )


        if n_trans == 0 or n_nontrans == 0:

            raise RuntimeError(
                "Empty group for {} / {}".format(
                    tissue,
                    metric
                )
            )


        # ----------------------------------------------------
        # Mann-Whitney U
        #
        # U is returned for first sample = transgressive
        # ----------------------------------------------------

        test = mannwhitneyu(
            trans,
            nontrans,
            alternative="two-sided"
        )


        u = float(
            test.statistic
        )

        p = float(
            test.pvalue
        )


        delta = cliffs_delta_from_u(
            u,
            n_trans,
            n_nontrans
        )


        med_trans = float(
            np.median(
                trans
            )
        )

        med_nontrans = float(
            np.median(
                nontrans
            )
        )


        mean_trans = float(
            np.mean(
                trans
            )
        )

        mean_nontrans = float(
            np.mean(
                nontrans
            )
        )


        median_difference = (
            med_trans
            -
            med_nontrans
        )


        if med_nontrans != 0:

            median_ratio = (
                med_trans
                /
                med_nontrans
            )

        else:

            median_ratio = np.nan


        results.append({

            "Assay":
                assay,

            "Metric":
                metric,

            "Tissue":
                tissue,

            "N_transgressive":
                n_trans,

            "N_nontransgressive":
                n_nontrans,

            "Median_transgressive":
                med_trans,

            "Median_nontransgressive":
                med_nontrans,

            "Median_difference":
                median_difference,

            "Median_ratio":
                median_ratio,

            "Mean_transgressive":
                mean_trans,

            "Mean_nontransgressive":
                mean_nontrans,

            "MannWhitney_U":
                u,

            "Cliffs_delta":
                delta,

            "Effect_size":
                effect_label(
                    delta
                ),

            "P":
                p
        })


        group_summary.append({

            "Assay":
                assay,

            "Metric":
                metric,

            "Tissue":
                tissue,

            "Group":
                "Strict_transgressive",

            "N":
                n_trans,

            "Median":
                med_trans,

            "Mean":
                mean_trans,

            "Q25":
                float(
                    np.percentile(
                        trans,
                        25
                    )
                ),

            "Q75":
                float(
                    np.percentile(
                        trans,
                        75
                    )
                )
        })


        group_summary.append({

            "Assay":
                assay,

            "Metric":
                metric,

            "Tissue":
                tissue,

            "Group":
                "Non_transgressive",

            "N":
                n_nontrans,

            "Median":
                med_nontrans,

            "Mean":
                mean_nontrans,

            "Q25":
                float(
                    np.percentile(
                        nontrans,
                        25
                    )
                ),

            "Q75":
                float(
                    np.percentile(
                        nontrans,
                        75
                    )
                )
        })


res = pd.DataFrame(
    results
)


# ============================================================
# MULTIPLE TESTING
# ============================================================

# all six tests together
res[
    "FDR_all6"
] = bh_adjust(
    res[
        "P"
    ].values
)


# within assay: 3 tissues separately
res[
    "FDR_within_assay"
] = np.nan


for assay in [
    "mCG",
    "ATAC"
]:

    idx = res[
        "Assay"
    ] == assay

    res.loc[
        idx,
        "FDR_within_assay"
    ] = bh_adjust(
        res.loc[
            idx,
            "P"
        ].values
    )


# ============================================================
# SORT
# ============================================================

assay_order = {
    "mCG": 0,
    "ATAC": 1
}

tissue_order = {
    "Gill": 0,
    "Kidney": 1,
    "Muscle": 2
}


res[
    "_assay_order"
] = res[
    "Assay"
].map(
    assay_order
)


res[
    "_tissue_order"
] = res[
    "Tissue"
].map(
    tissue_order
)


res = res.sort_values(
    [
        "_assay_order",
        "_tissue_order"
    ]
).drop(
    columns=[
        "_assay_order",
        "_tissue_order"
    ]
)


summary_df = pd.DataFrame(
    group_summary
)


# ============================================================
# ROW-LEVEL DATA USED
# ============================================================

used_cols = [
    "PanTE_ID",
    "Tissue",
    "Broad_class",
    "Parent_divergent",
    "Parent_genomic_bias",
    "Strict_transgressive",
    "Strict_transgressive_direction",
    "D_mCG_M",
    "D_mCG_P",
    "D_mCG_mean",
    "D_mCG_max",
    "D_ATAC_M",
    "D_ATAC_P",
    "D_ATAC_mean",
    "D_ATAC_max"
]


used_cols = [
    c
    for c in used_cols
    if c in df.columns
]


used = df[
    used_cols
].copy()


# ============================================================
# OUTPUT
# ============================================================

RESULT_OUT = os.path.join(
    OUTDIR,
    "strict_transgressive_vs_epigenetic_deviation.tsv"
)

GROUP_OUT = os.path.join(
    OUTDIR,
    "strict_transgressive_vs_epigenetic_deviation.group_summary.tsv"
)

USED_OUT = os.path.join(
    OUTDIR,
    "strict_transgressive_vs_epigenetic_deviation.data.tsv"
)

SUMMARY_OUT = os.path.join(
    OUTDIR,
    "strict_transgressive_vs_epigenetic_deviation.summary.txt"
)


res.to_csv(
    RESULT_OUT,
    sep="\t",
    index=False
)


summary_df.to_csv(
    GROUP_OUT,
    sep="\t",
    index=False
)


used.to_csv(
    USED_OUT,
    sep="\t",
    index=False,
    na_rep="NA"
)


with open(
    SUMMARY_OUT,
    "w"
) as out:

    out.write(
        "Strict transgressive PanTEs vs matched-parent epigenetic deviation\n"
    )

    out.write(
        "=================================================================\n\n"
    )

    out.write(
        "Primary metrics:\n"
    )

    out.write(
        "  D_mCG_mean  = mean absolute matched-parent mCG deviation\n"
    )

    out.write(
        "  D_ATAC_mean = mean absolute matched-parent ATAC deviation\n\n"
    )

    out.write(
        "Tests: two-sided Mann-Whitney U, separately by tissue.\n"
    )

    out.write(
        "Positive Cliff's delta means larger deviation in strict transgressive PanTEs.\n\n"
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
print("Strict transgressive vs epigenetic deviation")
print("==============================================")
print("")

print(
    res.to_string(
        index=False
    )
)

print("")
print(
    "Output directory:"
)

print(
    OUTDIR
)
