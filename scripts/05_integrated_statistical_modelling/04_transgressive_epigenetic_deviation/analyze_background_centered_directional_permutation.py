#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import numpy as np
import pandas as pd


INTEGRATION_DIR = os.environ.get("PANTE_INTEGRATION_DIR", os.path.join(os.getcwd(), "integration"))

MASTER = os.environ.get("PANTE_MASTER_MATRIX", os.path.join(INTEGRATION_DIR, "01_matrix", "PanTE_multiomics_master.RNA_testable.tsv"))

OUTDIR = os.environ.get("PANTE_OUTPUT_DIR", os.path.join(os.getcwd(), "results", "04_background_centered_directional_permutation"))

os.makedirs(OUTDIR, exist_ok=True)

TISSUES = [
    "Gill",
    "Kidney",
    "Muscle"
]

N_PERM = 10000
SEED = 20260909

rng = np.random.RandomState(
    SEED
)


# ============================================================
# LOAD
# ============================================================

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

    "ATAC_BSB_log2FPKM1",
    "ATAC_BTM_log2FPKM1",
    "ATAC_TC_log2FPKM1",
    "ATAC_BTP_log2FPKM1"
]

for c in required:

    if c not in df.columns:

        raise RuntimeError(
            "Missing required column: {}".format(c)
        )


df = df[
    (
        df["RNA_testable"] == 1
    )
    &
    (
        df["Tissue"].isin(TISSUES)
    )
].copy()


# ============================================================
# RECALCULATE SIGNED SHIFTS
# ============================================================

df["mCG_shift_M"] = (
    df["mCG_BTM"]
    -
    df["mCG_BSB"]
)

df["mCG_shift_P"] = (
    df["mCG_BTP"]
    -
    df["mCG_TC"]
)

df["mCG_shift_mean"] = (
    df["mCG_shift_M"]
    +
    df["mCG_shift_P"]
) / 2.0

df.loc[
    df["mCG_shift_M"].isna()
    |
    df["mCG_shift_P"].isna(),
    "mCG_shift_mean"
] = np.nan


df["ATAC_shift_M"] = (
    df["ATAC_BTM_log2FPKM1"]
    -
    df["ATAC_BSB_log2FPKM1"]
)

df["ATAC_shift_P"] = (
    df["ATAC_BTP_log2FPKM1"]
    -
    df["ATAC_TC_log2FPKM1"]
)

df["ATAC_shift_mean"] = (
    df["ATAC_shift_M"]
    +
    df["ATAC_shift_P"]
) / 2.0

df.loc[
    df["ATAC_shift_M"].isna()
    |
    df["ATAC_shift_P"].isna(),
    "ATAC_shift_mean"
] = np.nan


# ============================================================
# FUNCTIONS
# ============================================================

def empirical_p_greater(
    observed,
    null
):

    null = np.asarray(
        null,
        dtype=float
    )

    return (
        1.0
        +
        np.sum(
            null >= observed
        )
    ) / (
        len(null)
        +
        1.0
    )


def empirical_p_less(
    observed,
    null
):

    null = np.asarray(
        null,
        dtype=float
    )

    return (
        1.0
        +
        np.sum(
            null <= observed
        )
    ) / (
        len(null)
        +
        1.0
    )


# ============================================================
# ANALYSIS
# ============================================================

results = []
background_rows = []
data_rows = []


for tissue in TISSUES:

    tissue_df = df[
        df["Tissue"] == tissue
    ].copy()


    for assay in [
        "mCG",
        "ATAC"
    ]:

        shift_mean = (
            assay
            +
            "_shift_mean"
        )

        shift_M = (
            assay
            +
            "_shift_M"
        )

        shift_P = (
            assay
            +
            "_shift_P"
        )


        # ----------------------------------------------------
        # Non-transgressive background
        # ----------------------------------------------------

        bg = tissue_df[
            tissue_df["Strict_transgressive"] == "NO"
        ].copy()


        bg_mean = bg[
            shift_mean
        ].dropna()


        bg_M = bg[
            shift_M
        ].dropna()


        bg_P = bg[
            shift_P
        ].dropna()


        if (
            len(bg_mean) == 0
            or
            len(bg_M) == 0
            or
            len(bg_P) == 0
        ):

            raise RuntimeError(
                "Missing non-transgressive background: {} {}".format(
                    tissue,
                    assay
                )
            )


        median_bg_mean = float(
            np.median(
                bg_mean
            )
        )

        median_bg_M = float(
            np.median(
                bg_M
            )
        )

        median_bg_P = float(
            np.median(
                bg_P
            )
        )


        background_rows.append({

            "Tissue":
                tissue,

            "Assay":
                assay,

            "N_background_mean":
                len(bg_mean),

            "Median_background_shift_mean":
                median_bg_mean,

            "Median_background_shift_M":
                median_bg_M,

            "Median_background_shift_P":
                median_bg_P
        })


        # ----------------------------------------------------
        # Strict-transgressive observations
        # ----------------------------------------------------

        x = tissue_df[
            tissue_df["Strict_transgressive"] == "YES"
        ].copy()


        x = x[
            x[shift_mean].notna()
            &
            x[shift_M].notna()
            &
            x[shift_P].notna()
        ].copy()


        direction = (
            x[
                "Strict_transgressive_direction"
            ]
            .astype(str)
            .str.strip()
            .str.lower()
        )


        if not direction.isin(
            [
                "up",
                "down"
            ]
        ).all():

            raise RuntimeError(
                "Unexpected RNA direction in {} {}".format(
                    tissue,
                    assay
                )
            )


        sign = np.where(
            direction.values == "up",
            1.0,
            -1.0
        )


        # ----------------------------------------------------
        # Background-centered shifts
        # ----------------------------------------------------

        centered_mean = (
            x[
                shift_mean
            ].values
            -
            median_bg_mean
        )


        centered_M = (
            x[
                shift_M
            ].values
            -
            median_bg_M
        )


        centered_P = (
            x[
                shift_P
            ].values
            -
            median_bg_P
        )


        if assay == "ATAC":

            score = (
                sign
                *
                centered_mean
            )

            score_M = (
                sign
                *
                centered_M
            )

            score_P = (
                sign
                *
                centered_P
            )

        else:

            score = (
                -sign
                *
                centered_mean
            )

            score_M = (
                -sign
                *
                centered_M
            )

            score_P = (
                -sign
                *
                centered_P
            )


        obs_median = float(
            np.median(
                score
            )
        )

        obs_mean = float(
            np.mean(
                score
            )
        )


        nonzero = score[
            score != 0
        ]


        if len(nonzero):

            obs_prop = float(
                np.mean(
                    nonzero > 0
                )
            )

        else:

            obs_prop = np.nan


        obs_both = float(
            np.mean(
                (
                    score_M > 0
                )
                &
                (
                    score_P > 0
                )
            )
        )


        # ----------------------------------------------------
        # Permutation:
        # shuffle RNA direction labels within tissue/assay
        #
        # This preserves:
        #   - actual epigenetic shifts
        #   - number of Up and Down PanTEs
        # ----------------------------------------------------

        null_median = np.empty(
            N_PERM,
            dtype=float
        )

        null_mean = np.empty(
            N_PERM,
            dtype=float
        )

        null_prop = np.empty(
            N_PERM,
            dtype=float
        )

        null_both = np.empty(
            N_PERM,
            dtype=float
        )


        for i in range(
            N_PERM
        ):

            perm_sign = rng.permutation(
                sign
            )


            if assay == "ATAC":

                perm_score = (
                    perm_sign
                    *
                    centered_mean
                )

                perm_M = (
                    perm_sign
                    *
                    centered_M
                )

                perm_P = (
                    perm_sign
                    *
                    centered_P
                )

            else:

                perm_score = (
                    -perm_sign
                    *
                    centered_mean
                )

                perm_M = (
                    -perm_sign
                    *
                    centered_M
                )

                perm_P = (
                    -perm_sign
                    *
                    centered_P
                )


            null_median[i] = np.median(
                perm_score
            )

            null_mean[i] = np.mean(
                perm_score
            )


            nz = perm_score[
                perm_score != 0
            ]


            null_prop[i] = (
                np.mean(
                    nz > 0
                )
                if len(nz)
                else np.nan
            )


            null_both[i] = np.mean(
                (
                    perm_M > 0
                )
                &
                (
                    perm_P > 0
                )
            )


        # ----------------------------------------------------
        # output
        # ----------------------------------------------------

        results.append({

            "Tissue":
                tissue,

            "Assay":
                assay,

            "N_strict_transgressive":
                len(x),

            "N_up":
                int(
                    (
                        sign == 1
                    ).sum()
                ),

            "N_down":
                int(
                    (
                        sign == -1
                    ).sum()
                ),

            "Background_median_shift":
                median_bg_mean,

            "Observed_median_centered_score":
                obs_median,

            "Null_median_centered_score":
                float(
                    np.median(
                        null_median
                    )
                ),

            "P_enrichment_median_score":
                empirical_p_greater(
                    obs_median,
                    null_median
                ),

            "P_depletion_median_score":
                empirical_p_less(
                    obs_median,
                    null_median
                ),

            "Observed_mean_centered_score":
                obs_mean,

            "Null_mean_centered_score":
                float(
                    np.mean(
                        null_mean
                    )
                ),

            "P_enrichment_mean_score":
                empirical_p_greater(
                    obs_mean,
                    null_mean
                ),

            "P_depletion_mean_score":
                empirical_p_less(
                    obs_mean,
                    null_mean
                ),

            "Observed_concordant_fraction":
                obs_prop,

            "Null_concordant_fraction":
                float(
                    np.mean(
                        null_prop
                    )
                ),

            "P_enrichment_concordant_fraction":
                empirical_p_greater(
                    obs_prop,
                    null_prop
                ),

            "P_depletion_concordant_fraction":
                empirical_p_less(
                    obs_prop,
                    null_prop
                ),

            "Observed_both_haplotype_concordant_fraction":
                obs_both,

            "Null_both_haplotype_concordant_fraction":
                float(
                    np.mean(
                        null_both
                    )
                ),

            "P_enrichment_both_haplotypes":
                empirical_p_greater(
                    obs_both,
                    null_both
                ),

            "P_depletion_both_haplotypes":
                empirical_p_less(
                    obs_both,
                    null_both
                )
        })


        for j in range(
            len(x)
        ):

            data_rows.append({

                "PanTE_ID":
                    x.iloc[j][
                        "PanTE_ID"
                    ],

                "Tissue":
                    tissue,

                "Assay":
                    assay,

                "RNA_direction":
                    (
                        "Up"
                        if sign[j] == 1
                        else "Down"
                    ),

                "Raw_shift_mean":
                    float(
                        x.iloc[j][
                            shift_mean
                        ]
                    ),

                "Centered_shift_mean":
                    float(
                        centered_mean[j]
                    ),

                "Centered_concordance_score":
                    float(
                        score[j]
                    ),

                "Centered_M_score":
                    float(
                        score_M[j]
                    ),

                "Centered_P_score":
                    float(
                        score_P[j]
                    )
            })


# ============================================================
# OUTPUT
# ============================================================

res = pd.DataFrame(
    results
)

background = pd.DataFrame(
    background_rows
)

used = pd.DataFrame(
    data_rows
)


res.to_csv(
    os.path.join(
        OUTDIR,
        "background_centered_directional_permutation.tsv"
    ),
    sep="\t",
    index=False
)


background.to_csv(
    os.path.join(
        OUTDIR,
        "background_shift_medians.tsv"
    ),
    sep="\t",
    index=False
)


used.to_csv(
    os.path.join(
        OUTDIR,
        "background_centered_directional_data.tsv"
    ),
    sep="\t",
    index=False
)


with open(
    os.path.join(
        OUTDIR,
        "background_centered_directional_permutation.summary.txt"
    ),
    "w"
) as out:

    out.write(
        "Background-centered directional permutation analysis\n"
    )

    out.write(
        "====================================================\n\n"
    )

    out.write(
        "Background = median shift among RNA-testable non-transgressive PanTEs\n"
    )

    out.write(
        "within each tissue and assay.\n"
    )

    out.write(
        "RNA Up/Down labels were permuted within each tissue and assay.\n"
    )

    out.write(
        "N permutations = {:,}\n\n".format(
            N_PERM
        )
    )

    out.write(
        res.to_string(
            index=False
        )
    )

    out.write(
        "\n"
    )


print("")
print("==============================================")
print("Background shift medians")
print("==============================================")
print("")

print(
    background.to_string(
        index=False
    )
)


print("")
print("==============================================")
print("Background-centered permutation")
print("==============================================")
print("")

print(
    res.to_string(
        index=False
    )
)

print("")
print("Output:")
print(OUTDIR)
