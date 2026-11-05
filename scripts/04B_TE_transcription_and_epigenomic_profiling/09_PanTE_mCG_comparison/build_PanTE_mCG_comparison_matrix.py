#!/usr/bin/env python
# Public-release note: original executed source; replace local absolute paths with local paths before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import glob
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


BASE = os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52")
ROOT = os.path.join(BASE, "06_TE_WGBS")

INDIR = os.path.join(
    ROOT,
    "02_quantification"
)

OUTDIR = os.path.join(
    ROOT,
    "03_matrix"
)

os.makedirs(
    OUTDIR,
    exist_ok=True
)


# ============================================================
# Threshold schemes
# ============================================================

SCHEMES = {
    "relaxed": {
        "min_cpg": 3,
        "min_mean_depth": 2.0
    },

    "main": {
        "min_cpg": 5,
        "min_mean_depth": 3.0
    },

    "stringent": {
        "min_cpg": 10,
        "min_mean_depth": 5.0
    }
}


TISSUES = [
    "gill",
    "kidney",
    "muscle"
]


COMPARISONS = [
    (
        "BSB",
        "BTM",
        "BSB_to_BTM"
    ),

    (
        "TC",
        "BTP",
        "TC_to_BTP"
    )
]


# ============================================================
# Read all 27 sample files
# ============================================================

files = sorted(
    glob.glob(
        os.path.join(
            INDIR,
            "*.PanTE_weighted_mCG.tsv"
        )
    )
)

if len(files) != 27:

    raise RuntimeError(
        "Expected 27 files, found {}".format(
            len(files)
        )
    )


dfs = []

for f in files:

    x = pd.read_csv(
        f,
        sep="\t"
    )

    x["Mean_CpG_depth"] = (
        x["Total_C_reads"].astype(float)
        /
        x["N_CpG_sites"].astype(float)
    )

    dfs.append(
        x
    )


df = pd.concat(
    dfs,
    ignore_index=True
)


print("")
print("Total input records:")
print(len(df))


# ============================================================
# Summarize one genome state
# ============================================================

def summarize_state(
    x,
    tissue,
    state,
    min_cpg,
    min_depth
):

    z = x[
        (
            x["Tissue"] == tissue
        )
        &
        (
            x["Genome_state"] == state
        )
    ].copy()


    z["PASS"] = (
        (
            z["N_CpG_sites"]
            >=
            min_cpg
        )
        &
        (
            z["Mean_CpG_depth"]
            >=
            min_depth
        )
    )


    z = z[
        z["PASS"]
    ].copy()


    rows = []

    for pante, g in z.groupby(
        "PanTE_ID"
    ):

        nrep = g["Sample"].nunique()

        if nrep != 3:
            continue


        mc = g["mC"].sum()
        total = g["Total_C_reads"].sum()

        if total <= 0:
            continue


        pooled = (
            float(mc)
            /
            float(total)
        )


        replicate_values = (
            g["weighted_mCG"]
            .astype(float)
            .values
        )


        rows.append({
            "PanTE_ID":
                pante,

            "{}_Nrep".format(state):
                nrep,

            "{}_pooled_mCG".format(state):
                pooled,

            "{}_mean_mCG".format(state):
                np.mean(
                    replicate_values
                ),

            "{}_sd_mCG".format(state):
                np.std(
                    replicate_values,
                    ddof=1
                ),

            "{}_total_mC".format(state):
                mc,

            "{}_total_reads".format(state):
                total,

            "{}_median_CpG".format(state):
                g[
                    "N_CpG_sites"
                ].median(),

            "{}_median_depth".format(state):
                g[
                    "Mean_CpG_depth"
                ].median()
        })


    return pd.DataFrame(
        rows
    )


# ============================================================
# Main analysis
# ============================================================

qc_rows = []
corr_rows = []
delta_rows = []


for scheme_name, pars in SCHEMES.items():

    min_cpg = pars[
        "min_cpg"
    ]

    min_depth = pars[
        "min_mean_depth"
    ]


    print("")
    print(
        "=========================================="
    )

    print(
        "SCHEME: {}".format(
            scheme_name
        )
    )

    print(
        "N_CpG >= {}, mean depth >= {}".format(
            min_cpg,
            min_depth
        )
    )

    print(
        "=========================================="
    )


    for tissue in TISSUES:

        state_tables = {}


        for state in [
            "BSB",
            "TC",
            "BTM",
            "BTP"
        ]:

            s = summarize_state(
                df,
                tissue,
                state,
                min_cpg,
                min_depth
            )

            state_tables[
                state
            ] = s


            qc_rows.append({
                "Scheme":
                    scheme_name,

                "Tissue":
                    tissue,

                "Genome_state":
                    state,

                "N_PanTE_3of3":
                    len(s)
            })


        for parent, hybrid, comparison in COMPARISONS:

            a = state_tables[
                parent
            ]

            b = state_tables[
                hybrid
            ]


            merged = pd.merge(
                a,
                b,
                on="PanTE_ID",
                how="inner"
            )


            pcol = (
                "{}_pooled_mCG".format(
                    parent
                )
            )

            hcol = (
                "{}_pooled_mCG".format(
                    hybrid
                )
            )


            merged[
                "Delta_mCG"
            ] = (
                merged[hcol]
                -
                merged[pcol]
            )


            merged[
                "Abs_Delta_mCG"
            ] = np.abs(
                merged[
                    "Delta_mCG"
                ]
            )


            merged.insert(
                0,
                "Scheme",
                scheme_name
            )

            merged.insert(
                1,
                "Tissue",
                tissue
            )

            merged.insert(
                2,
                "Comparison",
                comparison
            )


            outfile = os.path.join(
                OUTDIR,
                "{}.{}.{}.PanTE_mCG.tsv".format(
                    scheme_name,
                    tissue,
                    comparison
                )
            )

            merged.to_csv(
                outfile,
                sep="\t",
                index=False
            )


            rho, pvalue = spearmanr(
                merged[pcol],
                merged[hcol]
            )


            corr_rows.append({
                "Scheme":
                    scheme_name,

                "Tissue":
                    tissue,

                "Comparison":
                    comparison,

                "N_PanTE":
                    len(merged),

                "Spearman_rho":
                    rho,

                "Spearman_p":
                    pvalue
            })


            for cutoff in [
                0.05,
                0.10,
                0.15
            ]:

                gain = (
                    merged[
                        "Delta_mCG"
                    ]
                    >=
                    cutoff
                ).sum()

                loss = (
                    merged[
                        "Delta_mCG"
                    ]
                    <=
                    -cutoff
                ).sum()

                stable = (
                    merged[
                        "Abs_Delta_mCG"
                    ]
                    <
                    cutoff
                ).sum()


                delta_rows.append({
                    "Scheme":
                        scheme_name,

                    "Tissue":
                        tissue,

                    "Comparison":
                        comparison,

                    "Delta_cutoff":
                        cutoff,

                    "N_PanTE":
                        len(merged),

                    "Stable_N":
                        stable,

                    "Gain_N":
                        gain,

                    "Loss_N":
                        loss,

                    "Stable_percent":
                        100.0
                        *
                        stable
                        /
                        len(merged),

                    "Gain_percent":
                        100.0
                        *
                        gain
                        /
                        len(merged),

                    "Loss_percent":
                        100.0
                        *
                        loss
                        /
                        len(merged)
                })


            print(
                "{} {}: N={}, rho={:.4f}, median_delta={:.4f}".format(
                    tissue,
                    comparison,
                    len(merged),
                    rho,
                    merged[
                        "Delta_mCG"
                    ].median()
                )
            )


# ============================================================
# Save summaries
# ============================================================

qc = pd.DataFrame(
    qc_rows
)

corr = pd.DataFrame(
    corr_rows
)

delta = pd.DataFrame(
    delta_rows
)


qc.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_mCG_filtering_QC.tsv"
    ),
    sep="\t",
    index=False
)


corr.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_mCG_parent_hybrid_correlations.tsv"
    ),
    sep="\t",
    index=False
)


delta.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_mCG_delta_sensitivity.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Print main scheme
# ============================================================

print("")
print(
    "=========================================="
)

print(
    "MAIN SCHEME: 3/3 replicates"
)

print(
    "N_CpG >= 5 and mean depth >= 3"
)

print(
    "=========================================="
)


print("")
print("Comparable PanTE counts:")

print(
    corr[
        corr["Scheme"] == "main"
    ].to_string(
        index=False
    )
)


print("")
print(
    "Delta sensitivity:"
)

print(
    delta[
        delta["Scheme"] == "main"
    ].to_string(
        index=False
    )
)


print("")
print(
    "Completed successfully."
)