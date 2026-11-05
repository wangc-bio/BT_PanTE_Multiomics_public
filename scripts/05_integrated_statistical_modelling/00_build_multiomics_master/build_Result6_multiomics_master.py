#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import numpy as np
import pandas as pd


# ============================================================
# PATHS
# ============================================================

BASE = os.environ.get("PANTE_BASE_DIR", os.getcwd())
INT = os.environ.get("PANTE_INTEGRATION_DIR", os.path.join(BASE, "09_Result6_multiomics_integration"))
OUT_MATRIX = os.environ.get("PANTE_MASTER_MATRIX_DIR", os.path.join(INT, "01_matrix"))
OUT_QC = os.environ.get("PANTE_MASTER_QC_DIR", os.path.join(INT, "00_input_check"))

GENOMIC_FILE = os.environ.get("PANTE_GENOMIC_ABUNDANCE_MATRIX", os.path.join(BASE, "03_abundance_resolved/matrix/PanTE_resolved_abundance_matrix.tsv"))
RNA_FILE = os.environ.get("PANTE_RNA_INHERITANCE_TABLE", os.path.join(BASE, "05_TE_RNA/05_DE/04_hybrid_pattern/PanTE_inheritance_all_tissues.tsv"))
ATAC_FILE = os.environ.get("PANTE_ATAC_RETENTION_MATRIX", os.path.join(BASE, "07_TE_ATAC/03_matrix/ATAC_parental_state_retention.four_state_matrix.tsv"))

WGBS_FILES = {
    "Gill": os.environ.get("PANTE_WGBS_GILL_RETENTION", os.path.join(BASE, "06_TE_WGBS/03_matrix/gill.parental_mCG_difference_retention.tsv")),
    "Kidney": os.environ.get("PANTE_WGBS_KIDNEY_RETENTION", os.path.join(BASE, "06_TE_WGBS/03_matrix/kidney.parental_mCG_difference_retention.tsv")),
    "Muscle": os.environ.get("PANTE_WGBS_MUSCLE_RETENTION", os.path.join(BASE, "06_TE_WGBS/03_matrix/muscle.parental_mCG_difference_retention.tsv"))
}

TISSUES = ["Gill", "Kidney", "Muscle"]

# ============================================================
# CHECK FILES
# ============================================================

all_files = [
    GENOMIC_FILE,
    RNA_FILE,
    ATAC_FILE
] + list(WGBS_FILES.values())

for f in all_files:

    if not os.path.exists(f):
        raise RuntimeError(
            "Missing input file: {}".format(f)
        )


# ============================================================
# HELPERS
# ============================================================

def normalize_tissue(x):

    x = str(x).strip().lower()

    return {
        "gill": "Gill",
        "kidney": "Kidney",
        "muscle": "Muscle"
    }.get(
        x,
        str(x)
    )


def yes_no(x):

    x = str(x).strip().upper()

    if x in [
        "YES",
        "TRUE",
        "1"
    ]:
        return "YES"

    if x in [
        "NO",
        "FALSE",
        "0"
    ]:
        return "NO"

    return "NA"


def safe_numeric(
    df,
    columns
):

    for c in columns:

        if c in df.columns:

            df[c] = pd.to_numeric(
                df[c],
                errors="coerce"
            )

    return df


# ============================================================
# 1. GENOMIC
# ============================================================

print("")
print("============================================")
print("Loading genomic abundance matrix")
print("============================================")

g = pd.read_csv(
    GENOMIC_FILE,
    sep="\t"
)

if g["PanTE_ID"].duplicated().any():

    raise RuntimeError(
        "Duplicated PanTE_ID found in genomic matrix."
    )


genomic_numeric = [
    "BSB_bp",
    "BSB_bp_per_Mb",
    "TC_bp",
    "TC_bp_per_Mb",
    "BTM_bp",
    "BTM_bp_per_Mb",
    "BTP_bp",
    "BTP_bp_per_Mb"
]

g = safe_numeric(
    g,
    genomic_numeric
)


# ------------------------------------------------------------
# Frozen parent-biased definition from Result 3.3:
#
# max parental raw bp >= 10 kb
# AND >= 2-fold parental difference
# ------------------------------------------------------------

g[
    "Parent_max_bp"
] = g[
    [
        "BSB_bp",
        "TC_bp"
    ]
].max(
    axis=1
)


g[
    "Parent_min_bp"
] = g[
    [
        "BSB_bp",
        "TC_bp"
    ]
].min(
    axis=1
)


def raw_fold_ratio(row):

    mx = row[
        "Parent_max_bp"
    ]

    mn = row[
        "Parent_min_bp"
    ]

    if pd.isna(mx) or pd.isna(mn):
        return np.nan

    if mx == 0:
        return np.nan

    if mn == 0:
        return np.inf

    return mx / mn


g[
    "Parent_raw_fold_ratio"
] = g.apply(
    raw_fold_ratio,
    axis=1
)


# ------------------------------------------------------------
# Frozen Result 3.3 definition:
#
# Parent-informative:
#   max parental RAW genomic occupancy >= 10 kb
#
# Parent-biased:
#   >= 2-fold difference in assembly-normalized genomic
#   abundance (bp/Mb)
#
# Direction is also defined from bp/Mb.
# ------------------------------------------------------------

g[
    "Parent_max_bp_per_Mb"
] = g[
    [
        "BSB_bp_per_Mb",
        "TC_bp_per_Mb"
    ]
].max(
    axis=1
)


g[
    "Parent_min_bp_per_Mb"
] = g[
    [
        "BSB_bp_per_Mb",
        "TC_bp_per_Mb"
    ]
].min(
    axis=1
)


def normalized_fold_ratio(row):

    mx = row[
        "Parent_max_bp_per_Mb"
    ]

    mn = row[
        "Parent_min_bp_per_Mb"
    ]

    if pd.isna(mx) or pd.isna(mn):
        return np.nan

    if mx == 0:
        return np.nan

    if mn == 0:
        return np.inf

    return mx / mn


g[
    "Parent_normalized_fold_ratio"
] = g.apply(
    normalized_fold_ratio,
    axis=1
)


def parent_bias(row):

    raw_max = row[
        "Parent_max_bp"
    ]

    norm_max = row[
        "Parent_max_bp_per_Mb"
    ]

    norm_min = row[
        "Parent_min_bp_per_Mb"
    ]

    bsb_norm = row[
        "BSB_bp_per_Mb"
    ]

    tc_norm = row[
        "TC_bp_per_Mb"
    ]

    if pd.isna(raw_max):
        return "NA"

    # parent-informative filter
    if raw_max < 10000:
        return "Nonbiased"

    if pd.isna(norm_max) or pd.isna(norm_min):
        return "NA"

    if norm_max == 0:
        return "Nonbiased"

    if norm_min == 0:

        if bsb_norm > tc_norm:
            return "BSB_biased"

        if tc_norm > bsb_norm:
            return "TC_biased"

        return "Nonbiased"

    if norm_max / norm_min < 2.0:
        return "Nonbiased"

    if bsb_norm > tc_norm:
        return "BSB_biased"

    if tc_norm > bsb_norm:
        return "TC_biased"

    return "Nonbiased"


g[
    "Parent_genomic_bias"
] = g.apply(
    parent_bias,
    axis=1
)


g[
    "Parent_divergent"
] = np.where(
    g[
        "Parent_genomic_bias"
    ].isin(
        [
            "BSB_biased",
            "TC_biased"
        ]
    ),
    "YES",
    "NO"
)


# ------------------------------------------------------------
# Continuous parental genomic divergence
#
# Use bp/Mb rather than raw bp, because this controls for
# genome-size differences.
#
# D = |BSB - TC| / (BSB + TC)
#
# range: 0-1
# ------------------------------------------------------------

gen_sum = (
    g[
        "BSB_bp_per_Mb"
    ]
    +
    g[
        "TC_bp_per_Mb"
    ]
)


g[
    "D_genomic"
] = np.where(
    gen_sum > 0,

    np.abs(
        g[
            "BSB_bp_per_Mb"
        ]
        -
        g[
            "TC_bp_per_Mb"
        ]
    )
    /
    gen_sum,

    np.nan
)


gen_bias_counts = (
    g[
        "Parent_genomic_bias"
    ]
    .value_counts(
        dropna=False
    )
)


n_parent_divergent_full = int(
    (
        g[
            "Parent_divergent"
        ] == "YES"
    ).sum()
)


print(
    "Genomic PanTE IDs: {:,}".format(
        len(g)
    )
)

print(
    "Parent-divergent PanTEs: {:,}".format(
        n_parent_divergent_full
    )
)

print(
    gen_bias_counts
)


if n_parent_divergent_full != 507:

    print(
        "[WARNING] Parent-divergent count is {}, "
        "not the expected frozen value 507.".format(
            n_parent_divergent_full
        )
    )

else:

    print(
        "[OK] Reproduced frozen parent-divergent count: 507"
    )


# retain useful genomic columns

g_keep = g[
    [
        "PanTE_ID",
        "Broad_class",
        "Cluster_size",
        "Representative_ID",

        "BSB_bp",
        "BSB_bp_per_Mb",

        "TC_bp",
        "TC_bp_per_Mb",

        "BTM_bp",
        "BTM_bp_per_Mb",

        "BTP_bp",
        "BTP_bp_per_Mb",

        "Parent_max_bp",
        "Parent_min_bp",
        "Parent_raw_fold_ratio",

        "Parent_max_bp_per_Mb",
        "Parent_min_bp_per_Mb",
        "Parent_normalized_fold_ratio",

        "Parent_genomic_bias",
        "Parent_divergent",
        "D_genomic"
    ]
].copy()


# ============================================================
# 2. RNA
# ============================================================

print("")
print("============================================")
print("Loading RNA inheritance matrix")
print("============================================")

r = pd.read_csv(
    RNA_FILE,
    sep="\t"
)


r[
    "Tissue"
] = r[
    "Tissue"
].map(
    normalize_tissue
)


r = r[
    r[
        "Tissue"
    ].isin(
        TISSUES
    )
].copy()


if r.duplicated(
    subset=[
        "PanTE_ID",
        "Tissue"
    ]
).any():

    dup = r[
        r.duplicated(
            subset=[
                "PanTE_ID",
                "Tissue"
            ],
            keep=False
        )
    ]

    raise RuntimeError(
        "Duplicated PanTE_ID x Tissue in RNA matrix: {}".format(
            len(dup)
        )
    )


rna_numeric = [
    "BSB_mean",
    "TC_mean",
    "BT_mean",
    "APE",

    "BSB_vs_TC_log2FC",
    "BSB_vs_TC_padj",

    "BT_vs_BSB_log2FC",
    "BT_vs_BSB_padj",

    "BT_vs_TC_log2FC",
    "BT_vs_TC_padj",

    "Comparable_all3"
]

r = safe_numeric(
    r,
    rna_numeric
)


# ------------------------------------------------------------
# RNA strict-testable background
#
# Must:
#   Comparable_all3 == 1
#   BT-vs-BSB padj available
#   BT-vs-TC padj available
#
# This avoids classifying untestable PanTEs as non-transgressive.
# ------------------------------------------------------------

r[
    "RNA_testable"
] = (
    (
        r[
            "Comparable_all3"
        ] == 1
    )
    &
    r[
        "BT_vs_BSB_padj"
    ].notna()
    &
    r[
        "BT_vs_TC_padj"
    ].notna()
).astype(int)


r[
    "Transgressive_original"
] = r[
    "Transgressive"
].map(
    yes_no
)


# In this project, "Transgressive" is the frozen strict
# transgressive definition:
#
# BT outside parental range
# AND BT-vs-BSB padj < 0.05
# AND BT-vs-TC padj < 0.05
#
# No extra fold-change cutoff.

r[
    "Strict_transgressive"
] = "NA"


mask_testable = (
    r[
        "RNA_testable"
    ] == 1
)


r.loc[
    mask_testable,
    "Strict_transgressive"
] = r.loc[
    mask_testable,
    "Transgressive_original"
]


r[
    "Strict_transgressive_direction"
] = "NA"


mask_trans = (
    r[
        "Strict_transgressive"
    ] == "YES"
)


r.loc[
    mask_trans,
    "Strict_transgressive_direction"
] = r.loc[
    mask_trans,
    "Transgressive_direction"
].fillna(
    "NA"
)


# ------------------------------------------------------------
# Independent QC:
# recalculate the frozen strict-transgressive rule
# ------------------------------------------------------------

parent_max_rna = r[
    [
        "BSB_mean",
        "TC_mean"
    ]
].max(
    axis=1
)


parent_min_rna = r[
    [
        "BSB_mean",
        "TC_mean"
    ]
].min(
    axis=1
)


calc_up = (
    mask_testable
    &
    (
        r[
            "BT_mean"
        ] > parent_max_rna
    )
    &
    (
        r[
            "BT_vs_BSB_padj"
        ] < 0.05
    )
    &
    (
        r[
            "BT_vs_TC_padj"
        ] < 0.05
    )
)


calc_down = (
    mask_testable
    &
    (
        r[
            "BT_mean"
        ] < parent_min_rna
    )
    &
    (
        r[
            "BT_vs_BSB_padj"
        ] < 0.05
    )
    &
    (
        r[
            "BT_vs_TC_padj"
        ] < 0.05
    )
)


r[
    "Strict_transgressive_recalculated"
] = np.where(
    calc_up | calc_down,
    "YES",
    "NO"
)


r[
    "Strict_transgressive_direction_recalculated"
] = np.where(
    calc_up,
    "Up",
    np.where(
        calc_down,
        "Down",
        "NA"
    )
)


rule_mismatch = (
    mask_testable
    &
    (
        r[
            "Strict_transgressive"
        ]
        !=
        r[
            "Strict_transgressive_recalculated"
        ]
    )
)


n_rule_mismatch = int(
    rule_mismatch.sum()
)


print(
    "RNA rows (Gill/Kidney/Muscle): {:,}".format(
        len(r)
    )
)

print(
    "Strict-testable RNA rows: {:,}".format(
        int(
            r[
                "RNA_testable"
            ].sum()
        )
    )
)

print(
    "Strict-transgressive rule mismatches: {:,}".format(
        n_rule_mismatch
    )
)


if n_rule_mismatch > 0:

    print(
        "[WARNING] Original and independently recalculated "
        "transgressive calls differ for {} rows.".format(
            n_rule_mismatch
        )
    )

    r[
        rule_mismatch
    ].to_csv(
        os.path.join(
            OUT_QC,
            "RNA_strict_transgressive_rule_mismatches.tsv"
        ),
        sep="\t",
        index=False
    )

else:

    print(
        "[OK] RNA strict-transgressive rule reproduced exactly."
    )


# Preserve original RNA Broad_class separately for QC
r = r.rename(
    columns={
        "Broad_class":
            "RNA_Broad_class"
    }
)


# ============================================================
# 3. WGBS
# ============================================================

print("")
print("============================================")
print("Loading WGBS matrices")
print("============================================")

w_list = []


for tissue in TISSUES:

    w = pd.read_csv(
        WGBS_FILES[
            tissue
        ],
        sep="\t"
    )


    if w[
        "PanTE_ID"
    ].duplicated().any():

        raise RuntimeError(
            "Duplicated WGBS PanTE_ID in {}".format(
                tissue
            )
        )


    w[
        "Tissue"
    ] = tissue


    w = w.rename(
        columns={
            "BSB_pooled_mCG":
                "mCG_BSB",

            "BTM_pooled_mCG":
                "mCG_BTM",

            "TC_pooled_mCG":
                "mCG_TC",

            "BTP_pooled_mCG":
                "mCG_BTP",

            "Parent_delta_mCG":
                "mCG_Parent_delta",

            "BT_haplotype_delta_mCG":
                "mCG_BT_delta",

            "Retention_residual":
                "mCG_Retention_residual",

            "Abs_parent_delta_mCG":
                "mCG_Abs_parent_delta",

            "Abs_BT_delta_mCG":
                "mCG_Abs_BT_delta"
        }
    )


    w_list.append(
        w
    )


w = pd.concat(
    w_list,
    ignore_index=True
)


w = safe_numeric(
    w,
    [
        "mCG_BSB",
        "mCG_BTM",
        "mCG_TC",
        "mCG_BTP",

        "mCG_Parent_delta",
        "mCG_BT_delta",
        "mCG_Retention_residual",
        "mCG_Abs_parent_delta",
        "mCG_Abs_BT_delta"
    ]
)


if w.duplicated(
    subset=[
        "PanTE_ID",
        "Tissue"
    ]
).any():

    raise RuntimeError(
        "Duplicated PanTE_ID x Tissue in WGBS."
    )


# matched-parent signed deviation

w[
    "mCG_BSB_to_BTM_delta"
] = (
    w[
        "mCG_BTM"
    ]
    -
    w[
        "mCG_BSB"
    ]
)


w[
    "mCG_TC_to_BTP_delta"
] = (
    w[
        "mCG_BTP"
    ]
    -
    w[
        "mCG_TC"
    ]
)


# absolute deviations

w[
    "D_mCG_M"
] = np.abs(
    w[
        "mCG_BSB_to_BTM_delta"
    ]
)


w[
    "D_mCG_P"
] = np.abs(
    w[
        "mCG_TC_to_BTP_delta"
    ]
)


# main epigenetic deviation metric:
# mean of matched-parent deviations

w[
    "D_mCG_mean"
] = (
    w[
        "D_mCG_M"
    ]
    +
    w[
        "D_mCG_P"
    ]
) / 2.0


# sensitivity metric

w[
    "D_mCG_max"
] = w[
    [
        "D_mCG_M",
        "D_mCG_P"
    ]
].max(
    axis=1
)


print(
    "WGBS PanTE x tissue rows: {:,}".format(
        len(w)
    )
)


# ============================================================
# 4. ATAC
# ============================================================

print("")
print("============================================")
print("Loading ATAC matrix")
print("============================================")

a = pd.read_csv(
    ATAC_FILE,
    sep="\t"
)


a[
    "Tissue"
] = a[
    "Tissue"
].map(
    normalize_tissue
)


a = a[
    a[
        "Tissue"
    ].isin(
        TISSUES
    )
].copy()


if a.duplicated(
    subset=[
        "PanTE_ID",
        "Tissue"
    ]
).any():

    raise RuntimeError(
        "Duplicated PanTE_ID x Tissue in ATAC matrix."
    )


a = a.rename(
    columns={
        "BSB_FPKM":
            "ATAC_BSB_FPKM",

        "BSB_log2FPKM1":
            "ATAC_BSB_log2FPKM1",

        "BSB_Assignable_bp":
            "ATAC_BSB_Assignable_bp",

        "TC_FPKM":
            "ATAC_TC_FPKM",

        "TC_log2FPKM1":
            "ATAC_TC_log2FPKM1",

        "TC_Assignable_bp":
            "ATAC_TC_Assignable_bp",

        "BTM_FPKM":
            "ATAC_BTM_FPKM",

        "BTM_log2FPKM1":
            "ATAC_BTM_log2FPKM1",

        "BTM_Assignable_bp":
            "ATAC_BTM_Assignable_bp",

        "BTP_FPKM":
            "ATAC_BTP_FPKM",

        "BTP_log2FPKM1":
            "ATAC_BTP_log2FPKM1",

        "BTP_Assignable_bp":
            "ATAC_BTP_Assignable_bp",

        "BSB_to_BTM_delta":
            "ATAC_BSB_to_BTM_delta",

        "TC_to_BTP_delta":
            "ATAC_TC_to_BTP_delta",

        "Parent_delta":
            "ATAC_Parent_delta",

        "BT_delta":
            "ATAC_BT_delta"
    }
)


a = safe_numeric(
    a,
    [
        "ATAC_BSB_FPKM",
        "ATAC_BSB_log2FPKM1",
        "ATAC_BSB_Assignable_bp",

        "ATAC_TC_FPKM",
        "ATAC_TC_log2FPKM1",
        "ATAC_TC_Assignable_bp",

        "ATAC_BTM_FPKM",
        "ATAC_BTM_log2FPKM1",
        "ATAC_BTM_Assignable_bp",

        "ATAC_BTP_FPKM",
        "ATAC_BTP_log2FPKM1",
        "ATAC_BTP_Assignable_bp",

        "ATAC_BSB_to_BTM_delta",
        "ATAC_TC_to_BTP_delta",

        "ATAC_Parent_delta",
        "ATAC_BT_delta"
    ]
)


a[
    "D_ATAC_M"
] = np.abs(
    a[
        "ATAC_BSB_to_BTM_delta"
    ]
)


a[
    "D_ATAC_P"
] = np.abs(
    a[
        "ATAC_TC_to_BTP_delta"
    ]
)


a[
    "D_ATAC_mean"
] = (
    a[
        "D_ATAC_M"
    ]
    +
    a[
        "D_ATAC_P"
    ]
) / 2.0


a[
    "D_ATAC_max"
] = a[
    [
        "D_ATAC_M",
        "D_ATAC_P"
    ]
].max(
    axis=1
)


print(
    "ATAC PanTE x tissue rows: {:,}".format(
        len(a)
    )
)


# ============================================================
# 5. MERGE
#
# RNA is the backbone because Result 3.6 asks whether
# transcriptionally remodeled PanTEs differ from the
# RNA-testable background.
# ============================================================

print("")
print("============================================")
print("Merging multi-omic master matrix")
print("============================================")


master = r.merge(
    g_keep,
    on="PanTE_ID",
    how="left",
    validate="many_to_one"
)


master = master.merge(
    w,
    on=[
        "PanTE_ID",
        "Tissue"
    ],
    how="left",
    validate="one_to_one"
)


master = master.merge(
    a,
    on=[
        "PanTE_ID",
        "Tissue"
    ],
    how="left",
    validate="one_to_one"
)


# ============================================================
# 6. QC FLAGS
# ============================================================

master[
    "WGBS_available"
] = (
    master[
        [
            "mCG_BSB",
            "mCG_BTM",
            "mCG_TC",
            "mCG_BTP"
        ]
    ]
    .notna()
    .all(
        axis=1
    )
).astype(int)


master[
    "ATAC_available"
] = (
    master[
        [
            "ATAC_BSB_log2FPKM1",
            "ATAC_BTM_log2FPKM1",
            "ATAC_TC_log2FPKM1",
            "ATAC_BTP_log2FPKM1"
        ]
    ]
    .notna()
    .all(
        axis=1
    )
).astype(int)


master[
    "Complete_multiomics"
] = (
    (
        master[
            "RNA_testable"
        ] == 1
    )
    &
    (
        master[
            "WGBS_available"
        ] == 1
    )
    &
    (
        master[
            "ATAC_available"
        ] == 1
    )
).astype(int)


# ============================================================
# Broad class consistency
# ============================================================

class_comparable = (
    master[
        "Broad_class"
    ].notna()
    &
    master[
        "RNA_Broad_class"
    ].notna()
)


class_mismatch = (
    class_comparable
    &
    (
        master[
            "Broad_class"
        ]
        !=
        master[
            "RNA_Broad_class"
        ]
    )
)


n_class_mismatch = int(
    class_mismatch.sum()
)


if n_class_mismatch > 0:

    master[
        class_mismatch
    ].to_csv(
        os.path.join(
            OUT_QC,
            "Broad_class_mismatches.tsv"
        ),
        sep="\t",
        index=False
    )


# ============================================================
# 7. OUTPUT MASTER MATRICES
# ============================================================

MASTER_ALL = os.path.join(
    OUT_MATRIX,
    "PanTE_multiomics_master.Gill_Kidney_Muscle.tsv"
)


master.to_csv(
    MASTER_ALL,
    sep="\t",
    index=False,
    na_rep="NA"
)


MASTER_TESTABLE = os.path.join(
    OUT_MATRIX,
    "PanTE_multiomics_master.RNA_testable.tsv"
)


master[
    master[
        "RNA_testable"
    ] == 1
].to_csv(
    MASTER_TESTABLE,
    sep="\t",
    index=False,
    na_rep="NA"
)


MASTER_COMPLETE = os.path.join(
    OUT_MATRIX,
    "PanTE_multiomics_master.complete_multiomics.tsv"
)


master[
    master[
        "Complete_multiomics"
    ] == 1
].to_csv(
    MASTER_COMPLETE,
    sep="\t",
    index=False,
    na_rep="NA"
)


# ============================================================
# 8. QC SUMMARY
# ============================================================

qc_rows = []


for tissue in TISSUES:

    d = master[
        master[
            "Tissue"
        ] == tissue
    ]


    test = d[
        d[
            "RNA_testable"
        ] == 1
    ]


    trans = test[
        test[
            "Strict_transgressive"
        ] == "YES"
    ]


    qc_rows.append({
        "Tissue":
            tissue,

        "RNA_rows":
            len(d),

        "RNA_testable":
            len(test),

        "Strict_transgressive":
            len(trans),

        "Transgressive_up":
            int(
                (
                    trans[
                        "Strict_transgressive_direction"
                    ] == "Up"
                ).sum()
            ),

        "Transgressive_down":
            int(
                (
                    trans[
                        "Strict_transgressive_direction"
                    ] == "Down"
                ).sum()
            ),

        "WGBS_available":
            int(
                d[
                    "WGBS_available"
                ].sum()
            ),

        "ATAC_available":
            int(
                d[
                    "ATAC_available"
                ].sum()
            ),

        "Complete_multiomics":
            int(
                d[
                    "Complete_multiomics"
                ].sum()
            ),

        "Parent_divergent_in_RNA_rows":
            int(
                (
                    d[
                        "Parent_divergent"
                    ] == "YES"
                ).sum()
            )
    })


qc_table = pd.DataFrame(
    qc_rows
)


QC_TABLE = os.path.join(
    OUT_QC,
    "Result6_master_matrix_counts.tsv"
)


qc_table.to_csv(
    QC_TABLE,
    sep="\t",
    index=False
)


QC_LOG = os.path.join(
    OUT_QC,
    "Result6_master_matrix_QC.txt"
)


with open(
    QC_LOG,
    "w"
) as out:

    out.write(
        "Result 3.6 multi-omics master matrix QC\n"
    )

    out.write(
        "=========================================\n\n"
    )


    out.write(
        "Genomic PanTE IDs: {:,}\n".format(
            len(g)
        )
    )


    out.write(
        "Frozen parent-divergent PanTEs: {:,}\n".format(
            n_parent_divergent_full
        )
    )


    out.write(
        "Expected frozen value: 507\n\n"
    )


    out.write(
        "RNA rows, Gill/Kidney/Muscle: {:,}\n".format(
            len(r)
        )
    )


    out.write(
        "RNA strict-testable rows: {:,}\n".format(
            int(
                r[
                    "RNA_testable"
                ].sum()
            )
        )
    )


    out.write(
        "Strict transgressive rule mismatches: {:,}\n".format(
            n_rule_mismatch
        )
    )


    out.write(
        "Broad-class mismatches: {:,}\n\n".format(
            n_class_mismatch
        )
    )


    out.write(
        "Master rows: {:,}\n".format(
            len(master)
        )
    )


    out.write(
        "Complete multi-omics rows: {:,}\n\n".format(
            int(
                master[
                    "Complete_multiomics"
                ].sum()
            )
        )
    )


    out.write(
        "Per-tissue counts:\n"
    )


    out.write(
        qc_table.to_string(
            index=False
        )
    )


    out.write(
        "\n"
    )


# ============================================================
# DONE
# ============================================================

print("")
print("============================================")
print("Master matrix completed successfully")
print("============================================")

print("")
print("All RNA-background rows:")
print(MASTER_ALL)

print("")
print("Strict RNA-testable background:")
print(MASTER_TESTABLE)

print("")
print("Complete RNA + WGBS + ATAC:")
print(MASTER_COMPLETE)

print("")
print("QC counts:")
print(QC_TABLE)

print("")
print("QC log:")
print(QC_LOG)

print("")
print(qc_table.to_string(index=False))
