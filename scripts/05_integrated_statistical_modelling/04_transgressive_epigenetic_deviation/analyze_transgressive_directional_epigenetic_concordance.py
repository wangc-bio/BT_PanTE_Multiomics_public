#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu, wilcoxon

# SciPy >=1.12 replaced scipy.stats.binom_test with binomtest.
# Return the scalar P value expected by the original analysis code.
try:
    from scipy.stats import binom_test
except ImportError:
    from scipy.stats import binomtest as _binomtest

    def binom_test(k, n, p=0.5, alternative="two-sided"):
        return _binomtest(k, n, p=p, alternative=alternative).pvalue


# ============================================================
# PATHS
# ============================================================

INTEGRATION_DIR = os.environ.get("PANTE_INTEGRATION_DIR", os.path.join(os.getcwd(), "integration"))

MASTER = os.environ.get("PANTE_MASTER_MATRIX", os.path.join(INTEGRATION_DIR, "01_matrix", "PanTE_multiomics_master.RNA_testable.tsv"))

OUTDIR = os.environ.get("PANTE_OUTPUT_DIR", os.path.join(os.getcwd(), "results", "03_transgressive_directional_epigenetic_concordance"))

os.makedirs(OUTDIR, exist_ok=True)

TISSUES = [
    "Gill",
    "Kidney",
    "Muscle"
]


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

        adj_ranked[
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
    ] = adj_ranked

    out[
        np.where(
            valid
        )[0]
    ] = tmp

    return out


def cliffs_delta_from_u(
    u,
    n1,
    n0
):

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


def effect_label(
    delta
):

    if not np.isfinite(
        delta
    ):
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


def safe_wilcoxon_against_zero(
    values
):

    values = np.asarray(
        values,
        dtype=float
    )

    values = values[
        np.isfinite(
            values
        )
    ]

    values = values[
        values != 0
    ]

    if len(values) == 0:

        return (
            np.nan,
            np.nan,
            0
        )

    try:

        stat, p = wilcoxon(
            values,
            alternative="two-sided"
        )

    except TypeError:

        stat, p = wilcoxon(
            values
        )

    return (
        float(stat),
        float(p),
        len(values)
    )


def mean_if_both(
    a,
    b
):

    out = (
        a
        +
        b
    ) / 2.0

    return out.where(
        a.notna()
        &
        b.notna()
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
    "RNA_testable",

    "Strict_transgressive",
    "Strict_transgressive_direction",

    "mCG_BSB",
    "mCG_BTM",
    "mCG_TC",
    "mCG_BTP",

    "mCG_BSB_to_BTM_delta",
    "mCG_TC_to_BTP_delta",

    "ATAC_BSB_log2FPKM1",
    "ATAC_BTM_log2FPKM1",

    "ATAC_TC_log2FPKM1",
    "ATAC_BTP_log2FPKM1",

    "ATAC_BSB_to_BTM_delta",
    "ATAC_TC_to_BTP_delta"
]


for col in required:

    if col not in df.columns:

        raise RuntimeError(
            "Missing required column: {}".format(
                col
            )
        )


# ============================================================
# KEEP STRICT TRANSGRESSIVE ONLY
# ============================================================

df = df[
    df[
        "Tissue"
    ].isin(
        TISSUES
    )
    &
    (
        df[
            "RNA_testable"
        ] == 1
    )
    &
    (
        df[
            "Strict_transgressive"
        ] == "YES"
    )
].copy()


# ============================================================
# RNA DIRECTION
# ============================================================

direction = (
    df[
        "Strict_transgressive_direction"
    ]
    .astype(str)
    .str.strip()
    .str.lower()
)


df[
    "RNA_direction"
] = np.where(
    direction.eq(
        "up"
    ),
    "Up",
    np.where(
        direction.eq(
            "down"
        ),
        "Down",
        "Other"
    )
)


if (
    df[
        "RNA_direction"
    ] == "Other"
).any():

    bad = sorted(
        df.loc[
            df[
                "RNA_direction"
            ] == "Other",
            "Strict_transgressive_direction"
        ]
        .astype(str)
        .unique()
    )

    raise RuntimeError(
        "Unexpected strict transgressive directions: {}".format(
            ",".join(
                bad
            )
        )
    )


df[
    "RNA_sign"
] = np.where(
    df[
        "RNA_direction"
    ] == "Up",
    1.0,
    -1.0
)


# ============================================================
# SIGNED EPIGENETIC SHIFTS
#
# Positive shift:
# hybrid haplotype > matched parent
# ============================================================

# ------------------------------------------------------------
# mCG
# ------------------------------------------------------------

df[
    "mCG_shift_M"
] = (
    df[
        "mCG_BTM"
    ]
    -
    df[
        "mCG_BSB"
    ]
)


df[
    "mCG_shift_P"
] = (
    df[
        "mCG_BTP"
    ]
    -
    df[
        "mCG_TC"
    ]
)


df[
    "mCG_shift_mean"
] = mean_if_both(
    df[
        "mCG_shift_M"
    ],
    df[
        "mCG_shift_P"
    ]
)


# ------------------------------------------------------------
# ATAC
#
# use log2(FPKM + 1)
# ------------------------------------------------------------

df[
    "ATAC_shift_M"
] = (
    df[
        "ATAC_BTM_log2FPKM1"
    ]
    -
    df[
        "ATAC_BSB_log2FPKM1"
    ]
)


df[
    "ATAC_shift_P"
] = (
    df[
        "ATAC_BTP_log2FPKM1"
    ]
    -
    df[
        "ATAC_TC_log2FPKM1"
    ]
)


df[
    "ATAC_shift_mean"
] = mean_if_both(
    df[
        "ATAC_shift_M"
    ],
    df[
        "ATAC_shift_P"
    ]
)


# ============================================================
# RNA-EPIGENETIC CONCORDANCE SCORE
#
# Positive score = direction consistent with RNA
#
# RNA up:
#   ATAC > 0
#   mCG  < 0
#
# RNA down:
#   ATAC < 0
#   mCG  > 0
# ============================================================

df[
    "mCG_concordance_score"
] = (
    -df[
        "RNA_sign"
    ]
    *
    df[
        "mCG_shift_mean"
    ]
)


df[
    "ATAC_concordance_score"
] = (
    df[
        "RNA_sign"
    ]
    *
    df[
        "ATAC_shift_mean"
    ]
)


# ============================================================
# HAPLOTYPE-SPECIFIC CONCORDANCE
# ============================================================

df[
    "mCG_M_concordant"
] = (
    -df[
        "RNA_sign"
    ]
    *
    df[
        "mCG_shift_M"
    ]
    >
    0
)


df[
    "mCG_P_concordant"
] = (
    -df[
        "RNA_sign"
    ]
    *
    df[
        "mCG_shift_P"
    ]
    >
    0
)


df[
    "ATAC_M_concordant"
] = (
    df[
        "RNA_sign"
    ]
    *
    df[
        "ATAC_shift_M"
    ]
    >
    0
)


df[
    "ATAC_P_concordant"
] = (
    df[
        "RNA_sign"
    ]
    *
    df[
        "ATAC_shift_P"
    ]
    >
    0
)


df[
    "mCG_both_haplotypes_concordant"
] = (
    df[
        "mCG_M_concordant"
    ]
    &
    df[
        "mCG_P_concordant"
    ]
    &
    df[
        "mCG_shift_M"
    ].notna()
    &
    df[
        "mCG_shift_P"
    ].notna()
)


df[
    "ATAC_both_haplotypes_concordant"
] = (
    df[
        "ATAC_M_concordant"
    ]
    &
    df[
        "ATAC_P_concordant"
    ]
    &
    df[
        "ATAC_shift_M"
    ].notna()
    &
    df[
        "ATAC_shift_P"
    ].notna()
)


# ============================================================
# QC:
# compare recalculated shifts with existing delta columns
# ============================================================

qc_rows = []


for assay, calc_m, calc_p, old_m, old_p in [

    (
        "mCG",
        "mCG_shift_M",
        "mCG_shift_P",
        "mCG_BSB_to_BTM_delta",
        "mCG_TC_to_BTP_delta"
    ),

    (
        "ATAC",
        "ATAC_shift_M",
        "ATAC_shift_P",
        "ATAC_BSB_to_BTM_delta",
        "ATAC_TC_to_BTP_delta"
    )

]:

    for hap, calc_col, old_col in [

        (
            "M",
            calc_m,
            old_m
        ),

        (
            "P",
            calc_p,
            old_p
        )

    ]:

        x = df[
            [
                calc_col,
                old_col
            ]
        ].dropna()


        if len(x) == 0:

            mae_same = np.nan
            mae_reverse = np.nan
            max_abs_same = np.nan

        else:

            mae_same = float(
                np.mean(
                    np.abs(
                        x[
                            calc_col
                        ]
                        -
                        x[
                            old_col
                        ]
                    )
                )
            )


            mae_reverse = float(
                np.mean(
                    np.abs(
                        x[
                            calc_col
                        ]
                        +
                        x[
                            old_col
                        ]
                    )
                )
            )


            max_abs_same = float(
                np.max(
                    np.abs(
                        x[
                            calc_col
                        ]
                        -
                        x[
                            old_col
                        ]
                    )
                )
            )


        if (
            np.isfinite(
                mae_same
            )
            and
            np.isfinite(
                mae_reverse
            )
            and
            mae_same <= mae_reverse
        ):

            match = "same_sign"

        else:

            match = "reverse_or_unclear"


        qc_rows.append({

            "Assay":
                assay,

            "Haplotype":
                hap,

            "Calculated_shift":
                calc_col,

            "Existing_delta":
                old_col,

            "N_compared":
                len(x),

            "MAE_same_sign":
                mae_same,

            "MAE_reverse_sign":
                mae_reverse,

            "Max_abs_difference_same_sign":
                max_abs_same,

            "Existing_delta_matches":
                match
        })


qc = pd.DataFrame(
    qc_rows
)


# ============================================================
# ANALYSIS 1:
# TRANSGRESSIVE UP VS DOWN
# ============================================================

ud_rows = []


for assay, metric, expected in [

    (
        "mCG",
        "mCG_shift_mean",
        "Up < Down"
    ),

    (
        "ATAC",
        "ATAC_shift_mean",
        "Up > Down"
    )

]:

    for tissue in TISSUES:

        x = df[
            df[
                "Tissue"
            ] == tissue
        ][
            [
                "RNA_direction",
                metric
            ]
        ].dropna()


        up = x.loc[
            x[
                "RNA_direction"
            ] == "Up",
            metric
        ].astype(
            float
        ).values


        down = x.loc[
            x[
                "RNA_direction"
            ] == "Down",
            metric
        ].astype(
            float
        ).values


        if (
            len(up) == 0
            or
            len(down) == 0
        ):

            u = np.nan
            p = np.nan
            delta = np.nan

        else:

            test = mannwhitneyu(
                up,
                down,
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
                len(up),
                len(down)
            )


        ud_rows.append({

            "Assay":
                assay,

            "Metric":
                metric,

            "Tissue":
                tissue,

            "Expected_direction":
                expected,

            "N_up":
                len(up),

            "N_down":
                len(down),

            "Median_up":
                float(
                    np.median(
                        up
                    )
                )
                if len(up)
                else np.nan,

            "Median_down":
                float(
                    np.median(
                        down
                    )
                )
                if len(down)
                else np.nan,

            "Median_up_minus_down":
                (
                    float(
                        np.median(
                            up
                        )
                        -
                        np.median(
                            down
                        )
                    )
                    if (
                        len(up)
                        and
                        len(down)
                    )
                    else np.nan
                ),

            "MannWhitney_U":
                u,

            "Cliffs_delta_Up_vs_Down":
                delta,

            "Effect_size":
                effect_label(
                    delta
                ),

            "P":
                p
        })


ud = pd.DataFrame(
    ud_rows
)


ud[
    "FDR_all6"
] = bh_adjust(
    ud[
        "P"
    ].values
)


ud[
    "FDR_within_assay"
] = np.nan


for assay in [
    "mCG",
    "ATAC"
]:

    idx = (
        ud[
            "Assay"
        ] == assay
    )

    ud.loc[
        idx,
        "FDR_within_assay"
    ] = bh_adjust(
        ud.loc[
            idx,
            "P"
        ].values
    )


# ============================================================
# ANALYSIS 2:
# CONCORDANCE SCORE > 0 ?
# ============================================================

score_rows = []


for assay, score_col in [

    (
        "mCG",
        "mCG_concordance_score"
    ),

    (
        "ATAC",
        "ATAC_concordance_score"
    )

]:

    for tissue in TISSUES:

        vals = (
            df.loc[
                df[
                    "Tissue"
                ] == tissue,
                score_col
            ]
            .dropna()
            .astype(float)
            .values
        )


        nonzero = vals[
            vals != 0
        ]


        n_pos = int(
            (
                nonzero > 0
            ).sum()
        )


        n_neg = int(
            (
                nonzero < 0
            ).sum()
        )


        n_nonzero = len(
            nonzero
        )


        if n_nonzero > 0:

            p_binom = float(
                binom_test(
                    n_pos,
                    n=n_nonzero,
                    p=0.5,
                    alternative="two-sided"
                )
            )

        else:

            p_binom = np.nan


        w_stat, w_p, n_w = safe_wilcoxon_against_zero(
            vals
        )


        score_rows.append({

            "Assay":
                assay,

            "Score":
                score_col,

            "Tissue":
                tissue,

            "N_available":
                len(vals),

            "Median_score":
                float(
                    np.median(
                        vals
                    )
                )
                if len(vals)
                else np.nan,

            "Mean_score":
                float(
                    np.mean(
                        vals
                    )
                )
                if len(vals)
                else np.nan,

            "N_positive_concordant":
                n_pos,

            "N_negative_discordant":
                n_neg,

            "N_nonzero":
                n_nonzero,

            "Concordant_pct_nonzero":
                (
                    100.0
                    *
                    n_pos
                    /
                    n_nonzero
                    if n_nonzero
                    else np.nan
                ),

            "Wilcoxon_signed_rank_stat":
                w_stat,

            "Wilcoxon_P_vs_zero":
                w_p,

            "Binomial_P_concordance_vs_50pct":
                p_binom
        })


score = pd.DataFrame(
    score_rows
)


score[
    "Wilcoxon_FDR_all6"
] = bh_adjust(
    score[
        "Wilcoxon_P_vs_zero"
    ].values
)


score[
    "Binomial_FDR_all6"
] = bh_adjust(
    score[
        "Binomial_P_concordance_vs_50pct"
    ].values
)


# ============================================================
# ANALYSIS 3:
# BOTH HAPLOTYPES CONCORDANT
# ============================================================

both_rows = []


for assay, col_m, col_p, both_col in [

    (
        "mCG",
        "mCG_shift_M",
        "mCG_shift_P",
        "mCG_both_haplotypes_concordant"
    ),

    (
        "ATAC",
        "ATAC_shift_M",
        "ATAC_shift_P",
        "ATAC_both_haplotypes_concordant"
    )

]:

    for tissue in TISSUES:

        x = df[
            df[
                "Tissue"
            ] == tissue
        ].copy()


        x = x[
            x[
                col_m
            ].notna()
            &
            x[
                col_p
            ].notna()
        ]


        n = len(
            x
        )


        n_both = int(
            x[
                both_col
            ].sum()
        )


        both_rows.append({

            "Assay":
                assay,

            "Tissue":
                tissue,

            "N_both_haplotypes_available":
                n,

            "N_both_haplotypes_concordant":
                n_both,

            "Pct_both_haplotypes_concordant":
                (
                    100.0
                    *
                    n_both
                    /
                    n
                    if n
                    else np.nan
                )
        })


both = pd.DataFrame(
    both_rows
)


# ============================================================
# OUTPUT DATA
# ============================================================

qc.to_csv(
    os.path.join(
        OUTDIR,
        "signed_shift_definition_QC.tsv"
    ),
    sep="\t",
    index=False
)


ud.to_csv(
    os.path.join(
        OUTDIR,
        "transgressive_up_vs_down.signed_epigenetic_shift.tsv"
    ),
    sep="\t",
    index=False
)


score.to_csv(
    os.path.join(
        OUTDIR,
        "RNA_epigenetic_directional_concordance.tsv"
    ),
    sep="\t",
    index=False
)


both.to_csv(
    os.path.join(
        OUTDIR,
        "both_haplotypes_directional_concordance.tsv"
    ),
    sep="\t",
    index=False
)


used_cols = [

    "PanTE_ID",
    "Tissue",
    "Broad_class",

    "Parent_divergent",
    "Parent_genomic_bias",

    "Strict_transgressive_direction",

    "RNA_direction",
    "RNA_sign",

    "mCG_shift_M",
    "mCG_shift_P",
    "mCG_shift_mean",

    "mCG_concordance_score",

    "mCG_M_concordant",
    "mCG_P_concordant",

    "mCG_both_haplotypes_concordant",

    "ATAC_shift_M",
    "ATAC_shift_P",
    "ATAC_shift_mean",

    "ATAC_concordance_score",

    "ATAC_M_concordant",
    "ATAC_P_concordant",

    "ATAC_both_haplotypes_concordant"
]


used_cols = [
    c
    for c in used_cols
    if c in df.columns
]


df[
    used_cols
].to_csv(
    os.path.join(
        OUTDIR,
        "strict_transgressive.directional_multiomics_data.tsv"
    ),
    sep="\t",
    index=False,
    na_rep="NA"
)


# ============================================================
# SUMMARY
# ============================================================

SUMMARY = os.path.join(
    OUTDIR,
    "directional_epigenetic_concordance.summary.txt"
)


with open(
    SUMMARY,
    "w"
) as out:

    out.write(
        "Strict transgressive PanTE directional epigenetic concordance\n"
    )

    out.write(
        "===========================================================\n\n"
    )


    out.write(
        "Signed-shift definition:\n"
    )

    out.write(
        "  mCG shift M = BT-M - BSB\n"
    )

    out.write(
        "  mCG shift P = BT-P - TC\n"
    )

    out.write(
        "  ATAC shift M = log2FPKM1(BT-M) - log2FPKM1(BSB)\n"
    )

    out.write(
        "  ATAC shift P = log2FPKM1(BT-P) - log2FPKM1(TC)\n\n"
    )


    out.write(
        "Positive concordance score means epigenetic change agrees with RNA direction:\n"
    )

    out.write(
        "  RNA up:   mCG down and/or ATAC up\n"
    )

    out.write(
        "  RNA down: mCG up and/or ATAC down\n\n"
    )


    out.write(
        "UP VS DOWN SIGNED SHIFTS\n"
    )

    out.write(
        "------------------------\n"
    )

    out.write(
        ud.to_string(
            index=False
        )
    )

    out.write(
        "\n\n"
    )


    out.write(
        "DIRECTIONAL CONCORDANCE SCORES\n"
    )

    out.write(
        "------------------------------\n"
    )

    out.write(
        score.to_string(
            index=False
        )
    )

    out.write(
        "\n\n"
    )


    out.write(
        "BOTH-HAPLOTYPE CONCORDANCE\n"
    )

    out.write(
        "--------------------------\n"
    )

    out.write(
        both.to_string(
            index=False
        )
    )

    out.write(
        "\n\n"
    )


    out.write(
        "DELTA-COLUMN QC\n"
    )

    out.write(
        "---------------\n"
    )

    out.write(
        qc.to_string(
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
print("Signed-shift definition QC")
print("==============================================")
print("")

print(
    qc.to_string(
        index=False
    )
)


print("")
print("==============================================")
print("Transgressive Up vs Down")
print("==============================================")
print("")

print(
    ud.to_string(
        index=False
    )
)


print("")
print("==============================================")
print("Directional concordance")
print("==============================================")
print("")

print(
    score.to_string(
        index=False
    )
)


print("")
print("==============================================")
print("Both-haplotype concordance")
print("==============================================")
print("")

print(
    both.to_string(
        index=False
    )
)


print("")
print("Output directory:")
print(OUTDIR)
