#!/usr/bin/env python
# Public-release note: original executed source; replace local absolute paths with paths from your environment before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import numpy as np
import pandas as pd


BASE = os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52")

GENOMIC = os.path.join(
    BASE,
    "03_abundance_resolved/matrix/"
    "PanTE_resolved_abundance_matrix.tsv"
)

RNA = os.path.join(
    BASE,
    "05_TE_RNA/05_DE/04_hybrid_pattern/"
    "PanTE_inheritance_all_tissues.tsv"
)

OUTDIR = os.path.join(
    BASE,
    "05_TE_RNA/05_DE/05_genomic_RNA_integration"
)

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# Genome sizes
# ============================================================

GSIZE = {
    "BSB": 1078381124.0,
    "TC":  1051140989.0,
    "BTM": 1083690569.0,
    "BTP": 1021101091.0
}


# ============================================================
# Detect raw assigned-bp columns
# ============================================================

def find_bp_col(df, genome):

    cols = list(df.columns)

    preferred = [
        genome + "_raw_bp",
        genome + "_assigned_bp",
        genome + "_resolved_bp",
        genome + "_bp",
        genome
    ]

    for c in preferred:
        if c in cols:
            return c

    candidates = []

    for c in cols:

        low = c.lower()

        if genome.lower() not in low:
            continue

        if "per_mb" in low:
            continue

        if "bp_per" in low:
            continue

        if (
            low.endswith("_bp")
            or "assigned_bp" in low
            or "resolved_bp" in low
            or "raw_bp" in low
        ):
            candidates.append(c)

    if len(candidates) == 1:
        return candidates[0]

    print("")
    print("Cannot uniquely detect raw bp column for:", genome)
    print("Candidates:", candidates)
    print("")
    print("Available genomic matrix columns:")
    for c in cols:
        print(c)

    raise RuntimeError(
        "Cannot detect genomic bp column for {}".format(genome)
    )


# ============================================================
# Spearman without scipy dependency
# ============================================================

def spearman(x, y):

    tmp = pd.DataFrame({
        "x": x,
        "y": y
    })

    tmp = tmp.replace(
        [np.inf, -np.inf],
        np.nan
    ).dropna()

    if len(tmp) < 3:
        return np.nan, len(tmp)

    rx = tmp["x"].rank()
    ry = tmp["y"].rank()

    rho = rx.corr(ry)

    return rho, len(tmp)


# ============================================================
# Load
# ============================================================

g = pd.read_csv(
    GENOMIC,
    sep="\t"
)

r = pd.read_csv(
    RNA,
    sep="\t"
)


if "PanTE_ID" not in g.columns:
    raise RuntimeError(
        "PanTE_ID missing from genomic matrix."
    )

if "PanTE_ID" not in r.columns:
    raise RuntimeError(
        "PanTE_ID missing from RNA table."
    )


# ============================================================
# Detect genomic bp columns
# ============================================================

bp_cols = {}

for genome in [
    "BSB",
    "TC",
    "BTM",
    "BTP"
]:

    bp_cols[genome] = find_bp_col(
        g,
        genome
    )


print("")
print("Detected genomic bp columns:")
for k, v in bp_cols.items():
    print(
        "{} -> {}".format(
            k,
            v
        )
    )
print("")


# ============================================================
# Keep genomic variables
# ============================================================

gg = g[
    [
        "PanTE_ID",
        bp_cols["BSB"],
        bp_cols["TC"],
        bp_cols["BTM"],
        bp_cols["BTP"]
    ]
].copy()


gg.columns = [
    "PanTE_ID",
    "BSB_TE_bp",
    "TC_TE_bp",
    "BTM_TE_bp",
    "BTP_TE_bp"
]


for c in [
    "BSB_TE_bp",
    "TC_TE_bp",
    "BTM_TE_bp",
    "BTP_TE_bp"
]:

    gg[c] = pd.to_numeric(
        gg[c],
        errors="coerce"
    ).fillna(0)


# ============================================================
# Genomic abundance: bp per Mb
# ============================================================

gg["BSB_genomic_bp_per_Mb"] = (
    gg["BSB_TE_bp"]
    / GSIZE["BSB"]
    * 1e6
)

gg["TC_genomic_bp_per_Mb"] = (
    gg["TC_TE_bp"]
    / GSIZE["TC"]
    * 1e6
)


gg["BT_TE_bp"] = (
    gg["BTM_TE_bp"]
    +
    gg["BTP_TE_bp"]
)


gg["BT_genomic_bp_per_Mb"] = (
    gg["BT_TE_bp"]
    /
    (
        GSIZE["BTM"]
        +
        GSIZE["BTP"]
    )
    * 1e6
)


# ============================================================
# Parent genomic bias
#
# same frozen informative threshold:
# max parental raw bp >= 10 kb
#
# >=2-fold = parental biased
# ============================================================

gg["Parent_max_bp"] = (
    gg[
        [
            "BSB_TE_bp",
            "TC_TE_bp"
        ]
    ].max(axis=1)
)


gg["Parent_genomic_log2ratio"] = np.log2(
    (
        gg["BSB_genomic_bp_per_Mb"]
        + 1e-6
    )
    /
    (
        gg["TC_genomic_bp_per_Mb"]
        + 1e-6
    )
)


def genomic_bias(row):

    if row["Parent_max_bp"] < 10000:
        return "Low_abundance"

    x = row[
        "Parent_genomic_log2ratio"
    ]

    if x >= 1:
        return "BSB_biased"

    if x <= -1:
        return "TC_biased"

    return "Nonbiased"


gg["Parent_genomic_bias"] = gg.apply(
    genomic_bias,
    axis=1
)


# ============================================================
# Merge RNA
# ============================================================

df = r.merge(
    gg,
    on="PanTE_ID",
    how="left"
)


# ============================================================
# RNA log values
# ============================================================

df["BSB_RNA_log2"] = np.log2(
    df["BSB_mean"] + 1
)

df["TC_RNA_log2"] = np.log2(
    df["TC_mean"] + 1
)

df["BT_RNA_log2"] = np.log2(
    df["BT_mean"] + 1
)


df["Parent_RNA_log2ratio"] = np.log2(
    (
        df["BSB_mean"] + 1
    )
    /
    (
        df["TC_mean"] + 1
    )
)


# ============================================================
# BT genomic position relative to parental range
# ============================================================

df["Parent_genomic_max"] = (
    df[
        [
            "BSB_genomic_bp_per_Mb",
            "TC_genomic_bp_per_Mb"
        ]
    ].max(axis=1)
)


df["Parent_genomic_min"] = (
    df[
        [
            "BSB_genomic_bp_per_Mb",
            "TC_genomic_bp_per_Mb"
        ]
    ].min(axis=1)
)


df["BT_genomic_vs_parent_max_log2"] = np.log2(
    (
        df["BT_genomic_bp_per_Mb"]
        + 1e-6
    )
    /
    (
        df["Parent_genomic_max"]
        + 1e-6
    )
)


df["BT_genomic_vs_parent_min_log2"] = np.log2(
    (
        df["BT_genomic_bp_per_Mb"]
        + 1e-6
    )
    /
    (
        df["Parent_genomic_min"]
        + 1e-6
    )
)


# ============================================================
# Genomic context of significant transgressive RNA
#
# NO arbitrary genomic fold cutoff.
# ============================================================

def transgressive_context(row):

    if row["Transgressive"] != "YES":
        return "Not_transgressive"

    direction = row[
        "Transgressive_direction"
    ]

    btg = row[
        "BT_genomic_bp_per_Mb"
    ]

    gmax = row[
        "Parent_genomic_max"
    ]

    gmin = row[
        "Parent_genomic_min"
    ]


    if direction == "Up":

        if btg > gmax:
            return (
                "RNA_up_with_genomic_abundance_increase"
            )
        else:
            return (
                "RNA_up_without_genomic_abundance_increase"
            )


    if direction == "Down":

        if btg < gmin:
            return (
                "RNA_down_with_genomic_abundance_decrease"
            )
        else:
            return (
                "RNA_down_without_genomic_abundance_decrease"
            )


    return "Unknown"


df["Transgressive_genomic_context"] = df.apply(
    transgressive_context,
    axis=1
)


# ============================================================
# Priority regulatory candidate
#
# Up:
# RNA transgressive up but genomic abundance not above parents
#
# Down:
# RNA transgressive down but genomic abundance not below parents
# ============================================================

df[
    "Copy_number_independent_transgressive_candidate"
] = np.where(
    df["Transgressive_genomic_context"].isin([
        "RNA_up_without_genomic_abundance_increase",
        "RNA_down_without_genomic_abundance_decrease"
    ]),
    "YES",
    "NO"
)


# ============================================================
# Save full integrated table
# ============================================================

df.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_genomic_RNA_integrated_all_tissues.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Correlation analyses
# ============================================================

corr_rows = []


for tissue in sorted(
    df["Tissue"].unique()
):

    x = df[
        (
            df["Tissue"] == tissue
        )
        &
        (
            df["Comparable_all3"] == 1
        )
    ].copy()


    tests = [
        (
            "BSB_genomic_vs_RNA",
            x["BSB_genomic_bp_per_Mb"],
            x["BSB_RNA_log2"]
        ),
        (
            "TC_genomic_vs_RNA",
            x["TC_genomic_bp_per_Mb"],
            x["TC_RNA_log2"]
        ),
        (
            "BT_genomic_vs_RNA",
            x["BT_genomic_bp_per_Mb"],
            x["BT_RNA_log2"]
        ),
        (
            "Parent_genomic_logratio_vs_RNA_logratio",
            x["Parent_genomic_log2ratio"],
            x["Parent_RNA_log2ratio"]
        )
    ]


    for name, a, b in tests:

        rho, n = spearman(
            a,
            b
        )

        corr_rows.append({
            "Tissue":
                tissue,

            "Comparison":
                name,

            "Spearman_rho":
                rho,

            "N":
                n
        })


corr = pd.DataFrame(
    corr_rows
)


corr.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_genomic_RNA_correlations.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Transgressive genomic-context summary
# ============================================================

trans = df[
    (
        df["Comparable_all3"] == 1
    )
    &
    (
        df["Transgressive"] == "YES"
    )
].copy()


trans_summary = (
    trans
    .groupby(
        [
            "Tissue",
            "Transgressive_direction",
            "Transgressive_genomic_context"
        ]
    )
    .size()
    .reset_index(
        name="N"
    )
)


tot = (
    trans
    .groupby(
        [
            "Tissue",
            "Transgressive_direction"
        ]
    )
    .size()
    .reset_index(
        name="Direction_total"
    )
)


trans_summary = trans_summary.merge(
    tot,
    on=[
        "Tissue",
        "Transgressive_direction"
    ],
    how="left"
)


trans_summary[
    "Percent_within_direction"
] = (
    trans_summary["N"]
    /
    trans_summary[
        "Direction_total"
    ]
    * 100.0
)


trans_summary.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_transgressive_genomic_context_summary.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Candidate list for ATAC/WGBS
# ============================================================

candidate = df[
    (
        df[
            "Copy_number_independent_transgressive_candidate"
        ] == "YES"
    )
].copy()


candidate.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_copy_number_independent_transgressive_candidates.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Parent-biased PanTE counts
# ============================================================

bias_counts = (
    gg[
        "Parent_genomic_bias"
    ]
    .value_counts()
)


print("")
print(
    "========================================"
)
print(
    "Parent genomic bias counts"
)
print(
    "========================================"
)

print(
    bias_counts.to_string()
)


# ============================================================
# Final
# ============================================================

print("")
print(
    "========================================"
)
print(
    "Genomic-RNA correlations"
)
print(
    "========================================"
)

print(
    corr.to_string(
        index=False
    )
)


print("")
print(
    "========================================"
)
print(
    "Transgressive genomic context"
)
print(
    "========================================"
)

print(
    trans_summary.to_string(
        index=False
    )
)


print("")
print(
    "Copy-number-independent transgressive candidates:"
)

print(
    len(candidate)
)

print("")
print(
    "Integration completed successfully."
)