#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import numpy as np
import pandas as pd

from scipy.stats import spearmanr, pearsonr
from scipy.stats import t as t_dist


# ============================================================
# PATHS
# ============================================================

INTEGRATION_DIR = os.environ.get("PANTE_INTEGRATION_DIR", os.path.join(os.getcwd(), "integration"))

MASTER = os.environ.get("PANTE_MASTER_MATRIX", os.path.join(INTEGRATION_DIR, "01_matrix", "PanTE_multiomics_master.RNA_testable.tsv"))

OUTDIR = os.environ.get("PANTE_OUTPUT_DIR", os.path.join(os.getcwd(), "results", "06_TEclass_abundance_controlled_sensitivity"))

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# SETTINGS
# ============================================================

TISSUES = [
    "Gill",
    "Kidney",
    "Muscle"
]

CLASSES = [
    "DNA",
    "LTR",
    "LINE",
    "SINE",
    "RC"
]

ASSAYS = {
    "mCG": "D_mCG_mean",
    "ATAC": "D_ATAC_mean"
}

MIN_N = 30


# ============================================================
# FUNCTIONS
# ============================================================

def bh_adjust(pvalues):

    p = np.asarray(
        pvalues,
        dtype=float
    )

    out = np.full(
        len(p),
        np.nan,
        dtype=float
    )

    valid = np.isfinite(
        p
    )

    if valid.sum() == 0:
        return out

    pv = p[valid]

    order = np.argsort(
        pv
    )

    ranked = pv[
        order
    ]

    n = len(
        ranked
    )

    adj_ranked = np.empty(
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

        adj_ranked[i] = min(
            running,
            1.0
        )

    tmp = np.empty(
        n,
        dtype=float
    )

    tmp[order] = adj_ranked

    out[
        np.where(
            valid
        )[0]
    ] = tmp

    return out


def residualize(
    y,
    z
):

    """
    Residualize y on one covariate z with intercept.
    """

    y = np.asarray(
        y,
        dtype=float
    )

    z = np.asarray(
        z,
        dtype=float
    )

    X = np.column_stack(
        [
            np.ones(
                len(z)
            ),
            z
        ]
    )

    beta = np.linalg.lstsq(
        X,
        y,
        rcond=None
    )[0]

    fitted = np.dot(
        X,
        beta
    )

    return (
        y
        -
        fitted
    )


def partial_spearman(
    x,
    y,
    covariate
):

    """
    Partial Spearman correlation controlling one covariate.

    Procedure:
      1. rank-transform x, y, covariate;
      2. residualize rank(x) and rank(y) on rank(covariate);
      3. Pearson correlation between residuals.

    P value uses t distribution with df = n - 3.
    """

    dat = pd.DataFrame({
        "x": x,
        "y": y,
        "z": covariate
    }).dropna()


    n = len(
        dat
    )

    if n < 4:

        return (
            np.nan,
            np.nan,
            n
        )


    rx = dat[
        "x"
    ].rank(
        method="average"
    ).values


    ry = dat[
        "y"
    ].rank(
        method="average"
    ).values


    rz = dat[
        "z"
    ].rank(
        method="average"
    ).values


    if (
        np.std(rx) == 0
        or
        np.std(ry) == 0
        or
        np.std(rz) == 0
    ):

        return (
            np.nan,
            np.nan,
            n
        )


    res_x = residualize(
        rx,
        rz
    )

    res_y = residualize(
        ry,
        rz
    )


    r = pearsonr(
        res_x,
        res_y
    )[0]


    r = float(
        r
    )


    dfree = (
        n
        -
        3
    )


    if (
        dfree <= 0
        or
        abs(r) >= 1
    ):

        p = np.nan

    else:

        tstat = (
            r
            *
            np.sqrt(
                dfree
                /
                (
                    1.0
                    -
                    r * r
                )
            )
        )


        p = float(
            2.0
            *
            t_dist.sf(
                abs(
                    tstat
                ),
                df=dfree
            )
        )


    return (
        r,
        p,
        n
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
    "Broad_class",

    "BSB_bp",
    "TC_bp",

    "BSB_bp_per_Mb",
    "TC_bp_per_Mb",

    "Parent_max_bp",
    "D_genomic",

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


if df.duplicated(
    subset=[
        "PanTE_ID",
        "Tissue"
    ]
).any():

    raise RuntimeError(
        "Duplicated PanTE_ID x Tissue detected."
    )


# ------------------------------------------------------------
# Parent-informative:
# same occupancy filter as Result 3.3
# ------------------------------------------------------------

df = df[
    df[
        "Parent_max_bp"
    ]
    >=
    10000
].copy()


# ------------------------------------------------------------
# Total parental normalized abundance
# ------------------------------------------------------------

df[
    "Parent_total_bp_per_Mb"
] = (
    df[
        "BSB_bp_per_Mb"
    ]
    +
    df[
        "TC_bp_per_Mb"
    ]
)


df[
    "log10_Parent_total_bp_per_Mb_plus1"
] = np.log10(
    df[
        "Parent_total_bp_per_Mb"
    ]
    +
    1.0
)


# ============================================================
# QC
# ============================================================

print("")
print("==============================================")
print("Parent-informative background QC")
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
        "{}: {:,} rows".format(
            tissue,
            len(x)
        )
    )


# ============================================================
# CLASS COUNTS
# ============================================================

count_rows = []


for tissue in TISSUES:

    for cls in CLASSES:

        x = df[
            (
                df[
                    "Tissue"
                ] == tissue
            )
            &
            (
                df[
                    "Broad_class"
                ] == cls
            )
        ]

        count_rows.append({

            "Tissue":
                tissue,

            "Broad_class":
                cls,

            "N_parent_informative":
                len(x),

            "N_mCG_available":
                int(
                    x[
                        "D_mCG_mean"
                    ].notna().sum()
                ),

            "N_ATAC_available":
                int(
                    x[
                        "D_ATAC_mean"
                    ].notna().sum()
                )
        })


class_counts = pd.DataFrame(
    count_rows
)


# ============================================================
# ANALYSIS
# ============================================================

rows = []
used_rows = []


for assay in [
    "mCG",
    "ATAC"
]:

    metric = ASSAYS[
        assay
    ]


    for tissue in TISSUES:

        for cls in CLASSES:

            x = df[
                (
                    df[
                        "Tissue"
                    ] == tissue
                )
                &
                (
                    df[
                        "Broad_class"
                    ] == cls
                )
            ][
                [
                    "PanTE_ID",
                    "D_genomic",
                    metric,
                    "Parent_total_bp_per_Mb",
                    "log10_Parent_total_bp_per_Mb_plus1"
                ]
            ].dropna().copy()


            n = len(
                x
            )


            # ------------------------------------------------
            # Raw within-class Spearman
            # ------------------------------------------------

            if n >= MIN_N:

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

            else:

                rho = np.nan
                p = np.nan


            # ------------------------------------------------
            # Partial Spearman controlling abundance
            # ------------------------------------------------

            if n >= MIN_N:

                partial_rho, partial_p, partial_n = partial_spearman(
                    x[
                        "D_genomic"
                    ].values,
                    x[
                        metric
                    ].values,
                    x[
                        "log10_Parent_total_bp_per_Mb_plus1"
                    ].values
                )

            else:

                partial_rho = np.nan
                partial_p = np.nan
                partial_n = n


            rows.append({

                "Assay":
                    assay,

                "Metric":
                    metric,

                "Tissue":
                    tissue,

                "Broad_class":
                    cls,

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

                "Median_parent_total_bp_per_Mb":
                    float(
                        x[
                            "Parent_total_bp_per_Mb"
                        ].median()
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

                "Spearman_P":
                    p,

                "Partial_Spearman_rho_control_abundance":
                    partial_rho,

                "Partial_Spearman_P_control_abundance":
                    partial_p,

                "Partial_N":
                    partial_n
            })


            temp = x.copy()

            temp[
                "Assay"
            ] = assay

            temp[
                "Metric"
            ] = metric

            temp[
                "Tissue"
            ] = tissue

            temp[
                "Broad_class"
            ] = cls

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
# 30 tests:
# 2 assays x 3 tissues x 5 classes
# ============================================================

res[
    "Spearman_FDR_all30"
] = bh_adjust(
    res[
        "Spearman_P"
    ].values
)


res[
    "Partial_FDR_all30"
] = bh_adjust(
    res[
        "Partial_Spearman_P_control_abundance"
    ].values
)


# ------------------------------------------------------------
# FDR within assay:
# 15 tests each
# ------------------------------------------------------------

res[
    "Spearman_FDR_within_assay"
] = np.nan

res[
    "Partial_FDR_within_assay"
] = np.nan


for assay in [
    "mCG",
    "ATAC"
]:

    idx = (
        res[
            "Assay"
        ]
        ==
        assay
    )


    res.loc[
        idx,
        "Spearman_FDR_within_assay"
    ] = bh_adjust(
        res.loc[
            idx,
            "Spearman_P"
        ].values
    )


    res.loc[
        idx,
        "Partial_FDR_within_assay"
    ] = bh_adjust(
        res.loc[
            idx,
            "Partial_Spearman_P_control_abundance"
        ].values
    )


# ------------------------------------------------------------
# FDR within tissue + assay:
# 5 classes each
# ------------------------------------------------------------

res[
    "Spearman_FDR_within_tissue_assay"
] = np.nan

res[
    "Partial_FDR_within_tissue_assay"
] = np.nan


for assay in [
    "mCG",
    "ATAC"
]:

    for tissue in TISSUES:

        idx = (
            (
                res[
                    "Assay"
                ]
                ==
                assay
            )
            &
            (
                res[
                    "Tissue"
                ]
                ==
                tissue
            )
        )


        res.loc[
            idx,
            "Spearman_FDR_within_tissue_assay"
        ] = bh_adjust(
            res.loc[
                idx,
                "Spearman_P"
            ].values
        )


        res.loc[
            idx,
            "Partial_FDR_within_tissue_assay"
        ] = bh_adjust(
            res.loc[
                idx,
                "Partial_Spearman_P_control_abundance"
            ].values
        )


# ============================================================
# DIRECTION CONSISTENCY SUMMARY
# ============================================================

summary_rows = []


for assay in [
    "mCG",
    "ATAC"
]:

    for tissue in TISSUES:

        x = res[
            (
                res[
                    "Assay"
                ] == assay
            )
            &
            (
                res[
                    "Tissue"
                ] == tissue
            )
            &
            (
                res[
                    "N"
                ] >= MIN_N
            )
        ].copy()


        summary_rows.append({

            "Assay":
                assay,

            "Tissue":
                tissue,

            "N_classes_tested":
                len(x),

            "N_raw_positive":
                int(
                    (
                        x[
                            "Spearman_rho"
                        ] > 0
                    ).sum()
                ),

            "N_raw_FDR05":
                int(
                    (
                        x[
                            "Spearman_FDR_within_tissue_assay"
                        ] < 0.05
                    ).sum()
                ),

            "N_partial_positive":
                int(
                    (
                        x[
                            "Partial_Spearman_rho_control_abundance"
                        ] > 0
                    ).sum()
                ),

            "N_partial_FDR05":
                int(
                    (
                        x[
                            "Partial_FDR_within_tissue_assay"
                        ] < 0.05
                    ).sum()
                ),

            "Median_raw_rho":
                float(
                    x[
                        "Spearman_rho"
                    ].median()
                )
                if len(x)
                else np.nan,

            "Median_partial_rho":
                float(
                    x[
                        "Partial_Spearman_rho_control_abundance"
                    ].median()
                )
                if len(x)
                else np.nan
        })


direction_summary = pd.DataFrame(
    summary_rows
)


# ============================================================
# OUTPUT
# ============================================================

RESULT_OUT = os.path.join(
    OUTDIR,
    "TEclass_Dgenomic_vs_epigenetic_deviation.partialSpearman.tsv"
)

COUNT_OUT = os.path.join(
    OUTDIR,
    "TEclass_parent_informative_counts.tsv"
)

DIRECTION_OUT = os.path.join(
    OUTDIR,
    "TEclass_direction_consistency_summary.tsv"
)

DATA_OUT = os.path.join(
    OUTDIR,
    "TEclass_Dgenomic_vs_epigenetic_deviation.data.tsv"
)

SUMMARY_OUT = os.path.join(
    OUTDIR,
    "TEclass_abundance_controlled_sensitivity.summary.txt"
)


res.to_csv(
    RESULT_OUT,
    sep="\t",
    index=False
)


class_counts.to_csv(
    COUNT_OUT,
    sep="\t",
    index=False
)


direction_summary.to_csv(
    DIRECTION_OUT,
    sep="\t",
    index=False
)


if len(
    used_rows
):

    used = pd.concat(
        used_rows,
        ignore_index=True
    )

else:

    used = pd.DataFrame()


used.to_csv(
    DATA_OUT,
    sep="\t",
    index=False,
    na_rep="NA"
)


with open(
    SUMMARY_OUT,
    "w"
) as out:

    out.write(
        "TE-class-stratified and abundance-controlled sensitivity analysis\n"
    )

    out.write(
        "===============================================================\n\n"
    )

    out.write(
        "Background: parent-informative PanTEs only "
        "(Parent_max_bp >= 10,000).\n\n"
    )

    out.write(
        "D_genomic = |BSB_bp_per_Mb - TC_bp_per_Mb| / "
        "(BSB_bp_per_Mb + TC_bp_per_Mb)\n"
    )

    out.write(
        "Covariate = log10(BSB_bp_per_Mb + TC_bp_per_Mb + 1)\n\n"
    )

    out.write(
        "Partial Spearman = rank residual correlation controlling "
        "parental total genomic abundance.\n\n"
    )

    out.write(
        "RESULTS\n"
    )

    out.write(
        "-------\n"
    )

    out.write(
        res.to_string(
            index=False
        )
    )

    out.write(
        "\n\n"
    )

    out.write(
        "DIRECTION CONSISTENCY\n"
    )

    out.write(
        "---------------------\n"
    )

    out.write(
        direction_summary.to_string(
            index=False
        )
    )

    out.write(
        "\n"
    )


# ============================================================
# SCREEN OUTPUT
# ============================================================

print("")
print("==============================================")
print("TE-class counts")
print("==============================================")
print("")

print(
    class_counts.to_string(
        index=False
    )
)


print("")
print("==============================================")
print("Class-stratified correlations")
print("==============================================")
print("")

print(
    res.to_string(
        index=False
    )
)


print("")
print("==============================================")
print("Direction consistency")
print("==============================================")
print("")

print(
    direction_summary.to_string(
        index=False
    )
)


print("")
print("Output directory:")
print(OUTDIR)
