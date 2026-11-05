#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import numpy as np
import pandas as pd

from scipy.stats import fisher_exact, mannwhitneyu

try:
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
except ImportError:
    raise RuntimeError(
        "statsmodels is required for abundance-adjusted logistic regression."
    )


# ============================================================
# PATHS
# ============================================================

INTEGRATION_DIR = os.environ.get("PANTE_INTEGRATION_DIR", os.path.join(os.getcwd(), "integration"))

MASTER = os.environ.get("PANTE_MASTER_MATRIX", os.path.join(INTEGRATION_DIR, "01_matrix", "PanTE_multiomics_master.RNA_testable.tsv"))

OUTDIR = os.environ.get("PANTE_OUTPUT_DIR", os.path.join(os.getcwd(), "results", "07_Figure6C_abundance_adjusted_reanalysis"))

os.makedirs(OUTDIR, exist_ok=True)

TISSUES = [
    "Gill",
    "Kidney",
    "Muscle"
]


# ============================================================
# HELPERS
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

    valid = np.isfinite(p)

    if valid.sum() == 0:
        return out

    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]

    n = len(ranked)

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
        np.where(valid)[0]
    ] = tmp

    return out


def pct(
    num,
    den
):

    if den == 0:
        return np.nan

    return (
        100.0
        *
        float(num)
        /
        float(den)
    )


def fit_logistic(
    formula,
    data,
    target_term="ParentDivergent_num",
    cov_type="HC0",
    groups=None
):

    model = smf.glm(
        formula=formula,
        data=data,
        family=sm.families.Binomial()
    )

    if cov_type == "cluster":

        fit = model.fit(
            cov_type="cluster",
            cov_kwds={
                "groups": groups
            }
        )

    else:

        fit = model.fit(
            cov_type=cov_type
        )

    beta = float(
        fit.params[
            target_term
        ]
    )

    se = float(
        fit.bse[
            target_term
        ]
    )

    p = float(
        fit.pvalues[
            target_term
        ]
    )

    ci = fit.conf_int().loc[
        target_term
    ]

    ci_low_beta = float(
        ci.iloc[0]
    )

    ci_high_beta = float(
        ci.iloc[1]
    )

    return {
        "beta":
            beta,

        "SE":
            se,

        "OR":
            float(
                np.exp(beta)
            ),

        "CI_low":
            float(
                np.exp(ci_low_beta)
            ),

        "CI_high":
            float(
                np.exp(ci_high_beta)
            ),

        "P":
            p,

        "AIC":
            float(
                fit.aic
            )
    }


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
    "Parent_divergent",

    "Parent_max_bp",

    "BSB_bp_per_Mb",
    "TC_bp_per_Mb",

    "Broad_class"
]


for col in required:

    if col not in df.columns:

        raise RuntimeError(
            "Missing required column: {}".format(
                col
            )
        )


# ============================================================
# BACKGROUND
#
# IMPORTANT:
# RNA-testable AND parent-informative only.
# ============================================================

df = df[
    (
        df[
            "Tissue"
        ].isin(
            TISSUES
        )
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
            "Parent_max_bp"
        ] >= 10000
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


if not df[
    "Parent_divergent"
].isin(
    [
        "YES",
        "NO"
    ]
).all():

    raise RuntimeError(
        "Unexpected Parent_divergent values."
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
# VARIABLES
# ============================================================

df[
    "ParentDivergent_num"
] = np.where(
    df[
        "Parent_divergent"
    ] == "YES",
    1,
    0
)


df[
    "StrictTransgressive_num"
] = np.where(
    df[
        "Strict_transgressive"
    ] == "YES",
    1,
    0
)


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
    "log10_parent_total_abundance"
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
print("===================================================")
print("Figure 6C reanalysis background")
print("===================================================")
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
        "{}: N={:,}; divergent={:,}; non-divergent={:,}; "
        "strict-transgressive={:,}".format(
            tissue,
            len(x),
            int(
                x[
                    "ParentDivergent_num"
                ].sum()
            ),
            int(
                (
                    x[
                        "ParentDivergent_num"
                    ] == 0
                ).sum()
            ),
            int(
                x[
                    "StrictTransgressive_num"
                ].sum()
            )
        )
    )


# ============================================================
# ABUNDANCE IMBALANCE QC
# ============================================================

abundance_rows = []


for tissue in TISSUES:

    x = df[
        df[
            "Tissue"
        ] == tissue
    ].copy()


    div = x.loc[
        x[
            "ParentDivergent_num"
        ] == 1,
        "Parent_total_bp_per_Mb"
    ].astype(float)


    nondiv = x.loc[
        x[
            "ParentDivergent_num"
        ] == 0,
        "Parent_total_bp_per_Mb"
    ].astype(float)


    test = mannwhitneyu(
        div,
        nondiv,
        alternative="two-sided"
    )


    abundance_rows.append({

        "Tissue":
            tissue,

        "N_divergent":
            len(div),

        "N_nondivergent":
            len(nondiv),

        "Median_total_abundance_divergent":
            float(
                np.median(
                    div
                )
            ),

        "Median_total_abundance_nondivergent":
            float(
                np.median(
                    nondiv
                )
            ),

        "Median_log10_abundance_divergent":
            float(
                np.median(
                    np.log10(
                        div
                        +
                        1.0
                    )
                )
            ),

        "Median_log10_abundance_nondivergent":
            float(
                np.median(
                    np.log10(
                        nondiv
                        +
                        1.0
                    )
                )
            ),

        "MannWhitney_U":
            float(
                test.statistic
            ),

        "MannWhitney_P":
            float(
                test.pvalue
            )
    })


abundance_qc = pd.DataFrame(
    abundance_rows
)


abundance_qc[
    "FDR_BH"
] = bh_adjust(
    abundance_qc[
        "MannWhitney_P"
    ].values
)


# ============================================================
# PER-TISSUE ANALYSIS
# ============================================================

rows = []
decile_rows = []


for tissue in TISSUES:

    x = df[
        df[
            "Tissue"
        ] == tissue
    ].copy()


    # --------------------------------------------------------
    # Abundance deciles
    #
    # rank first to avoid qcut problems with ties
    # --------------------------------------------------------

    abundance_rank = x[
        "log10_parent_total_abundance"
    ].rank(
        method="first"
    )


    x[
        "Abundance_decile"
    ] = (
        pd.qcut(
            abundance_rank,
            q=10,
            labels=False
        )
        +
        1
    ).astype(int)


    # --------------------------------------------------------
    # Crude informative-only 2 x 2 table
    # --------------------------------------------------------

    div = (
        x[
            "ParentDivergent_num"
        ] == 1
    )

    trans = (
        x[
            "StrictTransgressive_num"
        ] == 1
    )


    a = int(
        (
            div
            &
            trans
        ).sum()
    )

    b = int(
        (
            div
            &
            ~trans
        ).sum()
    )

    c = int(
        (
            ~div
            &
            trans
        ).sum()
    )

    d = int(
        (
            ~div
            &
            ~trans
        ).sum()
    )


    fisher_or, fisher_p = fisher_exact(
        [
            [
                a,
                b
            ],
            [
                c,
                d
            ]
        ],
        alternative="two-sided"
    )


    # --------------------------------------------------------
    # Model 1:
    # continuous abundance adjustment
    #
    # OR for ParentDivergent_num is the quantity of interest
    # --------------------------------------------------------

    m_cont = fit_logistic(

        formula=(
            "StrictTransgressive_num ~ "
            "ParentDivergent_num + "
            "log10_parent_total_abundance"
        ),

        data=x,

        cov_type="HC0"
    )


    # --------------------------------------------------------
    # Model 2:
    # non-parametric abundance adjustment
    #
    # abundance decile as categorical covariate
    # --------------------------------------------------------

    m_decile = fit_logistic(

        formula=(
            "StrictTransgressive_num ~ "
            "ParentDivergent_num + "
            "C(Abundance_decile)"
        ),

        data=x,

        cov_type="HC0"
    )


    # --------------------------------------------------------
    # Abundance coefficient from Model 1
    # --------------------------------------------------------

    model_abundance = smf.glm(

        formula=(
            "StrictTransgressive_num ~ "
            "ParentDivergent_num + "
            "log10_parent_total_abundance"
        ),

        data=x,

        family=sm.families.Binomial()

    ).fit(
        cov_type="HC0"
    )


    abundance_beta = float(
        model_abundance.params[
            "log10_parent_total_abundance"
        ]
    )


    abundance_ci = model_abundance.conf_int().loc[
        "log10_parent_total_abundance"
    ]


    abundance_or = float(
        np.exp(
            abundance_beta
        )
    )


    abundance_ci_low = float(
        np.exp(
            abundance_ci.iloc[0]
        )
    )


    abundance_ci_high = float(
        np.exp(
            abundance_ci.iloc[1]
        )
    )


    abundance_p = float(
        model_abundance.pvalues[
            "log10_parent_total_abundance"
        ]
    )


    rows.append({

        "Tissue":
            tissue,

        "N_parent_informative_RNA_testable":
            len(x),

        "N_parent_divergent":
            a + b,

        "N_parent_nondivergent":
            c + d,

        "Strict_transgressive_divergent":
            a,

        "Strict_transgressive_nondivergent":
            c,

        "Strict_transgressive_pct_divergent":
            pct(
                a,
                a + b
            ),

        "Strict_transgressive_pct_nondivergent":
            pct(
                c,
                c + d
            ),

        "Informative_only_Fisher_OR":
            float(
                fisher_or
            ),

        "Informative_only_Fisher_P":
            float(
                fisher_p
            ),

        "Abundance_adjusted_OR":
            m_cont[
                "OR"
            ],

        "Abundance_adjusted_95CI_low":
            m_cont[
                "CI_low"
            ],

        "Abundance_adjusted_95CI_high":
            m_cont[
                "CI_high"
            ],

        "Abundance_adjusted_P":
            m_cont[
                "P"
            ],

        "Abundance_decile_adjusted_OR":
            m_decile[
                "OR"
            ],

        "Abundance_decile_adjusted_95CI_low":
            m_decile[
                "CI_low"
            ],

        "Abundance_decile_adjusted_95CI_high":
            m_decile[
                "CI_high"
            ],

        "Abundance_decile_adjusted_P":
            m_decile[
                "P"
            ],

        "Abundance_effect_OR_per_1log10":
            abundance_or,

        "Abundance_effect_95CI_low":
            abundance_ci_low,

        "Abundance_effect_95CI_high":
            abundance_ci_high,

        "Abundance_effect_P":
            abundance_p
    })


    # --------------------------------------------------------
    # Descriptive outcome rate by abundance decile
    # --------------------------------------------------------

    for decile, y in x.groupby(
        "Abundance_decile"
    ):

        for status_value, status_name in [
            (
                1,
                "Parent_divergent"
            ),
            (
                0,
                "Parent_nondivergent"
            )
        ]:

            z = y[
                y[
                    "ParentDivergent_num"
                ] == status_value
            ]


            if len(z) == 0:

                continue


            decile_rows.append({

                "Tissue":
                    tissue,

                "Abundance_decile":
                    int(
                        decile
                    ),

                "Parent_status":
                    status_name,

                "N":
                    len(z),

                "Median_parent_total_bp_per_Mb":
                    float(
                        z[
                            "Parent_total_bp_per_Mb"
                        ].median()
                    ),

                "N_strict_transgressive":
                    int(
                        z[
                            "StrictTransgressive_num"
                        ].sum()
                    ),

                "Pct_strict_transgressive":
                    100.0
                    *
                    float(
                        z[
                            "StrictTransgressive_num"
                        ].mean()
                    )
            })


res = pd.DataFrame(
    rows
)


# ============================================================
# MULTIPLE-TESTING CORRECTION
# ============================================================

res[
    "Informative_only_Fisher_FDR"
] = bh_adjust(
    res[
        "Informative_only_Fisher_P"
    ].values
)


res[
    "Abundance_adjusted_FDR"
] = bh_adjust(
    res[
        "Abundance_adjusted_P"
    ].values
)


res[
    "Abundance_decile_adjusted_FDR"
] = bh_adjust(
    res[
        "Abundance_decile_adjusted_P"
    ].values
)


res[
    "Abundance_effect_FDR"
] = bh_adjust(
    res[
        "Abundance_effect_P"
    ].values
)


decile_summary = pd.DataFrame(
    decile_rows
)


# ============================================================
# POOLED ANALYSIS
#
# Repeated PanTE_ID across tissues:
# cluster-robust SE by PanTE_ID.
# ============================================================

pooled = df.copy()


# ------------------------------------------------------------
# Continuous abundance model
# ------------------------------------------------------------

pooled_cont = fit_logistic(

    formula=(
        "StrictTransgressive_num ~ "
        "ParentDivergent_num + "
        "log10_parent_total_abundance + "
        "C(Tissue)"
    ),

    data=pooled,

    cov_type="cluster",

    groups=pooled[
        "PanTE_ID"
    ].values
)


# ------------------------------------------------------------
# Tissue-specific abundance deciles
# ------------------------------------------------------------

pooled[
    "Abundance_decile"
] = np.nan


for tissue in TISSUES:

    idx = (
        pooled[
            "Tissue"
        ] == tissue
    )


    ranks = pooled.loc[
        idx,
        "log10_parent_total_abundance"
    ].rank(
        method="first"
    )


    pooled.loc[
        idx,
        "Abundance_decile"
    ] = (
        pd.qcut(
            ranks,
            q=10,
            labels=False
        )
        +
        1
    ).values


pooled[
    "Abundance_decile"
] = pooled[
    "Abundance_decile"
].astype(int)


pooled[
    "Tissue_abundance_stratum"
] = (
    pooled[
        "Tissue"
    ].astype(str)
    +
    "_D"
    +
    pooled[
        "Abundance_decile"
    ].astype(str)
)


pooled_strata = fit_logistic(

    formula=(
        "StrictTransgressive_num ~ "
        "ParentDivergent_num + "
        "C(Tissue_abundance_stratum)"
    ),

    data=pooled,

    cov_type="cluster",

    groups=pooled[
        "PanTE_ID"
    ].values
)


pooled_result = pd.DataFrame([

    {
        "Model":
            "Continuous_abundance_plus_tissue_clustered_by_PanTE",

        "N":
            len(
                pooled
            ),

        "N_unique_PanTE":
            pooled[
                "PanTE_ID"
            ].nunique(),

        "Adjusted_OR_parent_divergent":
            pooled_cont[
                "OR"
            ],

        "Adjusted_95CI_low":
            pooled_cont[
                "CI_low"
            ],

        "Adjusted_95CI_high":
            pooled_cont[
                "CI_high"
            ],

        "P":
            pooled_cont[
                "P"
            ]
    },

    {
        "Model":
            "Tissue_specific_abundance_decile_strata_clustered_by_PanTE",

        "N":
            len(
                pooled
            ),

        "N_unique_PanTE":
            pooled[
                "PanTE_ID"
            ].nunique(),

        "Adjusted_OR_parent_divergent":
            pooled_strata[
                "OR"
            ],

        "Adjusted_95CI_low":
            pooled_strata[
                "CI_low"
            ],

        "Adjusted_95CI_high":
            pooled_strata[
                "CI_high"
            ],

        "P":
            pooled_strata[
                "P"
            ]
    }
])


# ============================================================
# SAVE DATA USED
# ============================================================

used_cols = [
    "PanTE_ID",
    "Tissue",
    "Broad_class",

    "Parent_divergent",
    "ParentDivergent_num",

    "Strict_transgressive",
    "StrictTransgressive_num",

    "Parent_max_bp",

    "BSB_bp_per_Mb",
    "TC_bp_per_Mb",

    "Parent_total_bp_per_Mb",
    "log10_parent_total_abundance"
]


df[
    used_cols
].to_csv(

    os.path.join(
        OUTDIR,
        "Figure6C_reanalysis.data.tsv"
    ),

    sep="\t",
    index=False,
    na_rep="NA"
)


# ============================================================
# OUTPUT
# ============================================================

res.to_csv(

    os.path.join(
        OUTDIR,
        "Figure6C_reanalysis.per_tissue.tsv"
    ),

    sep="\t",
    index=False
)


pooled_result.to_csv(

    os.path.join(
        OUTDIR,
        "Figure6C_reanalysis.pooled_clustered.tsv"
    ),

    sep="\t",
    index=False
)


abundance_qc.to_csv(

    os.path.join(
        OUTDIR,
        "Figure6C_parent_abundance_balance_QC.tsv"
    ),

    sep="\t",
    index=False
)


decile_summary.to_csv(

    os.path.join(
        OUTDIR,
        "Figure6C_abundance_decile_outcome_summary.tsv"
    ),

    sep="\t",
    index=False
)


# ============================================================
# SUMMARY
# ============================================================

SUMMARY = os.path.join(
    OUTDIR,
    "Figure6C_reanalysis.summary.txt"
)


with open(
    SUMMARY,
    "w"
) as out:

    out.write(
        "Figure 6C abundance-adjusted reanalysis\n"
    )

    out.write(
        "=======================================\n\n"
    )

    out.write(
        "Background:\n"
    )

    out.write(
        "RNA-testable AND parent-informative PanTEs only "
        "(Parent_max_bp >= 10,000).\n\n"
    )

    out.write(
        "Parent total genomic abundance:\n"
    )

    out.write(
        "BSB_bp_per_Mb + TC_bp_per_Mb\n\n"
    )

    out.write(
        "Continuous covariate:\n"
    )

    out.write(
        "log10(parent total genomic abundance + 1)\n\n"
    )


    out.write(
        "PER-TISSUE RESULTS\n"
    )

    out.write(
        "------------------\n"
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
        "POOLED CLUSTER-ROBUST RESULTS\n"
    )

    out.write(
        "-----------------------------\n"
    )

    out.write(
        pooled_result.to_string(
            index=False
        )
    )

    out.write(
        "\n\n"
    )


    out.write(
        "ABUNDANCE BALANCE QC\n"
    )

    out.write(
        "--------------------\n"
    )

    out.write(
        abundance_qc.to_string(
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
print("===================================================")
print("Per-tissue reanalysis")
print("===================================================")
print("")

print(
    res.to_string(
        index=False
    )
)


print("")
print("===================================================")
print("Pooled cluster-robust reanalysis")
print("===================================================")
print("")

print(
    pooled_result.to_string(
        index=False
    )
)


print("")
print("===================================================")
print("Parental abundance balance QC")
print("===================================================")
print("")

print(
    abundance_qc.to_string(
        index=False
    )
)


print("")
print("Output:")
print(OUTDIR)