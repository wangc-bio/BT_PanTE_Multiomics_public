#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import numpy as np
import pandas as pd


# ============================================================
# Settings
# ============================================================

PANTE_ROOT = os.environ.get(
    "PANTE_ROOT",
    "PATH_TO_PanTE_v52"
)

ROOT = os.path.join(
    PANTE_ROOT,
    "05_TE_RNA"
)

DE_ROOT = os.path.join(
    ROOT,
    "05_DE/03_DESeq2"
)

OUTDIR = os.path.join(
    ROOT,
    "05_DE/04_hybrid_pattern"
)

MEMBERSHIP = os.path.join(
    PANTE_ROOT,
    "02_cluster",
    "PanTE_refined_cluster_members.tsv"
)

TISSUES = [
    "Brain",
    "Gill",
    "Kidney",
    "Liver",
    "Muscle"
]

FDR = 0.05

# strong-effect annotation only
STRONG_LFC = 1.0

os.makedirs(
    OUTDIR,
    exist_ok=True
)


# ============================================================
# Helpers
# ============================================================

def sig(p):
    if pd.isna(p):
        return False
    return float(p) < FDR


def nonsig(p):
    if pd.isna(p):
        return False
    return float(p) >= FDR


def strong(lfc, p):
    if pd.isna(lfc) or pd.isna(p):
        return False

    return (
        float(p) < FDR
        and abs(float(lfc)) >= STRONG_LFC
    )


# ============================================================
# PanTE broad class
# ============================================================

membership = pd.read_csv(
    MEMBERSHIP,
    sep="\t"
)

pc = (
    membership[
        [
            "PanTE_ID",
            "Broad_class"
        ]
    ]
    .drop_duplicates()
)

check = (
    pc.groupby(
        "PanTE_ID"
    )["Broad_class"]
    .nunique()
)

if (check > 1).any():
    raise RuntimeError(
        "Some PanTE IDs have multiple Broad_class values."
    )

pante_class = (
    pc
    .drop_duplicates(
        "PanTE_ID"
    )
    .set_index(
        "PanTE_ID"
    )["Broad_class"]
    .to_dict()
)


# ============================================================
# Load one tissue
# ============================================================

def load_tissue(tissue):

    d = os.path.join(
        DE_ROOT,
        tissue
    )

    parent_file = os.path.join(
        d,
        "BSB_vs_TC.DESeq2.tsv"
    )

    hb_file = os.path.join(
        d,
        "BT_vs_BSB.DESeq2.tsv"
    )

    ht_file = os.path.join(
        d,
        "BT_vs_TC.DESeq2.tsv"
    )

    means_file = os.path.join(
        d,
        "PanTE_normalized_group_means.tsv"
    )

    for f in [
        parent_file,
        hb_file,
        ht_file,
        means_file
    ]:
        if not os.path.exists(f):
            raise RuntimeError(
                "Missing input: {}".format(f)
            )


    parent = pd.read_csv(
        parent_file,
        sep="\t"
    )[
        [
            "PanTE_ID",
            "log2FoldChange",
            "padj"
        ]
    ]

    parent.columns = [
        "PanTE_ID",
        "BSB_vs_TC_log2FC",
        "BSB_vs_TC_padj"
    ]


    hb = pd.read_csv(
        hb_file,
        sep="\t"
    )[
        [
            "PanTE_ID",
            "log2FoldChange",
            "padj"
        ]
    ]

    hb.columns = [
        "PanTE_ID",
        "BT_vs_BSB_log2FC",
        "BT_vs_BSB_padj"
    ]


    ht = pd.read_csv(
        ht_file,
        sep="\t"
    )[
        [
            "PanTE_ID",
            "log2FoldChange",
            "padj"
        ]
    ]

    ht.columns = [
        "PanTE_ID",
        "BT_vs_TC_log2FC",
        "BT_vs_TC_padj"
    ]


    means = pd.read_csv(
        means_file,
        sep="\t"
    )


    df = (
        means
        .merge(
            parent,
            on="PanTE_ID"
        )
        .merge(
            hb,
            on="PanTE_ID"
        )
        .merge(
            ht,
            on="PanTE_ID"
        )
    )


    df["Tissue"] = tissue

    df["Broad_class"] = (
        df["PanTE_ID"]
        .map(
            pante_class
        )
    )

    df["APE"] = (
        df["BSB_mean"]
        +
        df["TC_mean"]
    ) / 2.0


    return df


# ============================================================
# Final classifier
# ============================================================

def classify(row):

    b = float(
        row["BSB_mean"]
    )

    t = float(
        row["TC_mean"]
    )

    h = float(
        row["BT_mean"]
    )


    parent_p = row[
        "BSB_vs_TC_padj"
    ]

    hb_p = row[
        "BT_vs_BSB_padj"
    ]

    ht_p = row[
        "BT_vs_TC_padj"
    ]


    parent_sig = sig(
        parent_p
    )

    hb_sig = sig(
        hb_p
    )

    ht_sig = sig(
        ht_p
    )


    parent_ns = nonsig(
        parent_p
    )

    hb_ns = nonsig(
        hb_p
    )

    ht_ns = nonsig(
        ht_p
    )


    # --------------------------------------------------------
    # No-change label
    # --------------------------------------------------------

    no_change = (
        parent_ns
        and hb_ns
        and ht_ns
    )


    # --------------------------------------------------------
    # Strong effects
    # --------------------------------------------------------

    strong_hb = strong(
        row["BT_vs_BSB_log2FC"],
        hb_p
    )

    strong_ht = strong(
        row["BT_vs_TC_log2FC"],
        ht_p
    )


    # --------------------------------------------------------
    # Non-comparable
    # --------------------------------------------------------

    if int(
        row["Comparable_all3"]
    ) != 1:

        return {
            "Inheritance_class":
                "Genomic_noncomparable",

            "Pattern":
                "NA",

            "No_change":
                "NA",

            "Transgressive":
                "NA",

            "Transgressive_direction":
                "NA",

            "Strong_transgressive":
                "NA",

            "BSB_TC_significant":
                parent_sig,

            "BT_BSB_significant":
                hb_sig,

            "BT_TC_significant":
                ht_sig
        }


    # ========================================================
    # 1. TC expression dominance
    #
    # BT similar to TC,
    # but significantly different from BSB.
    # ========================================================

    if (
        ht_ns
        and hb_sig
    ):

        if t < b:
            pattern = "II"
        else:
            pattern = "XI"

        return {
            "Inheritance_class":
                "TC_expression_dominance",

            "Pattern":
                pattern,

            "No_change":
                "NO",

            "Transgressive":
                "NO",

            "Transgressive_direction":
                "NA",

            "Strong_transgressive":
                "NO",

            "BSB_TC_significant":
                parent_sig,

            "BT_BSB_significant":
                hb_sig,

            "BT_TC_significant":
                ht_sig
        }


    # ========================================================
    # 2. BSB expression dominance
    # ========================================================

    if (
        hb_ns
        and ht_sig
    ):

        if b < t:
            pattern = "IV"
        else:
            pattern = "IX"

        return {
            "Inheritance_class":
                "BSB_expression_dominance",

            "Pattern":
                pattern,

            "No_change":
                "NO",

            "Transgressive":
                "NO",

            "Transgressive_direction":
                "NA",

            "Strong_transgressive":
                "NO",

            "BSB_TC_significant":
                parent_sig,

            "BT_BSB_significant":
                hb_sig,

            "BT_TC_significant":
                ht_sig
        }


    # ========================================================
    # 3. Additive
    #
    # BT lies between parental means.
    #
    # Significance is not required.
    # ========================================================

    if (
        b < h < t
    ):

        return {
            "Inheritance_class":
                "Additive",

            "Pattern":
                "I",

            "No_change":
                "YES" if no_change else "NO",

            "Transgressive":
                "NO",

            "Transgressive_direction":
                "NA",

            "Strong_transgressive":
                "NO",

            "BSB_TC_significant":
                parent_sig,

            "BT_BSB_significant":
                hb_sig,

            "BT_TC_significant":
                ht_sig
        }


    if (
        t < h < b
    ):

        return {
            "Inheritance_class":
                "Additive",

            "Pattern":
                "XII",

            "No_change":
                "YES" if no_change else "NO",

            "Transgressive":
                "NO",

            "Transgressive_direction":
                "NA",

            "Strong_transgressive":
                "NO",

            "BSB_TC_significant":
                parent_sig,

            "BT_BSB_significant":
                hb_sig,

            "BT_TC_significant":
                ht_sig
        }


    # ========================================================
    # 4. Above parental range
    #
    # Main class depends ONLY on position.
    #
    # Transgressive=YES only when BT is significantly
    # higher than BOTH parents.
    # ========================================================

    if (
        h > b
        and h > t
    ):

        trans = (
            hb_sig
            and ht_sig
        )

        strong_trans = (
            trans
            and strong_hb
            and strong_ht
        )


        if parent_ns:

            pattern = "VIII"

        elif b < t:

            pattern = "V"

        else:

            pattern = "VI"


        return {
            "Inheritance_class":
                "Above_parental_range",

            "Pattern":
                pattern,

            "No_change":
                "YES" if no_change else "NO",

            "Transgressive":
                "YES" if trans else "NO",

            "Transgressive_direction":
                "Up" if trans else "NA",

            "Strong_transgressive":
                "YES" if strong_trans else "NO",

            "BSB_TC_significant":
                parent_sig,

            "BT_BSB_significant":
                hb_sig,

            "BT_TC_significant":
                ht_sig
        }


    # ========================================================
    # 5. Below parental range
    #
    # Transgressive=YES only when BT is significantly
    # lower than BOTH parents.
    # ========================================================

    if (
        h < b
        and h < t
    ):

        trans = (
            hb_sig
            and ht_sig
        )

        strong_trans = (
            trans
            and strong_hb
            and strong_ht
        )


        if parent_ns:

            pattern = "VII"

        elif t < b:

            pattern = "III"

        else:

            pattern = "X"


        return {
            "Inheritance_class":
                "Below_parental_range",

            "Pattern":
                pattern,

            "No_change":
                "YES" if no_change else "NO",

            "Transgressive":
                "YES" if trans else "NO",

            "Transgressive_direction":
                "Down" if trans else "NA",

            "Strong_transgressive":
                "YES" if strong_trans else "NO",

            "BSB_TC_significant":
                parent_sig,

            "BT_BSB_significant":
                hb_sig,

            "BT_TC_significant":
                ht_sig
        }


    # ========================================================
    # Exact numerical ties
    # ========================================================

    return {
        "Inheritance_class":
            "Unresolved_tie",

        "Pattern":
            "Unresolved_tie",

        "No_change":
            "YES" if no_change else "NO",

        "Transgressive":
            "NO",

        "Transgressive_direction":
            "NA",

        "Strong_transgressive":
            "NO",

        "BSB_TC_significant":
            parent_sig,

        "BT_BSB_significant":
            hb_sig,

        "BT_TC_significant":
            ht_sig
    }


# ============================================================
# Run
# ============================================================

all_tables = []

class_rows = []
pattern_rows = []
trans_rows = []
nochange_rows = []
qc_rows = []


print("")
print(
    "=============================================="
)
print(
    "FINAL PanTE inheritance-state classification"
)
print(
    "=============================================="
)

print(
    "Main classes:"
)
print(
    "Additive / TC dominance / BSB dominance / "
    "Above parental range / Below parental range"
)

print(
    "Transgressive requires BT vs BOTH parents padj < 0.05"
)

print(
    "No main-class log2FC cutoff."
)

print("")


for tissue in TISSUES:

    print(
        "===== {} =====".format(
            tissue
        )
    )


    df = load_tissue(
        tissue
    )


    result = pd.DataFrame(
        [
            classify(row)
            for _, row in df.iterrows()
        ]
    )


    out = pd.concat(
        [
            df.reset_index(drop=True),
            result.reset_index(drop=True)
        ],
        axis=1
    )


    # overwrite old tissue results
    out.to_csv(
        os.path.join(
            OUTDIR,
            "{}.PanTE_inheritance.tsv".format(
                tissue
            )
        ),
        sep="\t",
        index=False
    )


    out.to_csv(
        os.path.join(
            OUTDIR,
            "{}.PanTE_12pattern.tsv".format(
                tissue
            )
        ),
        sep="\t",
        index=False
    )


    all_tables.append(
        out
    )


    comp = out[
        out[
            "Comparable_all3"
        ] == 1
    ].copy()


    total = len(
        comp
    )


    # ========================================================
    # Main class summary
    # ========================================================

    ct = (
        comp[
            "Inheritance_class"
        ]
        .value_counts()
    )


    for cls, n in ct.items():

        class_rows.append({
            "Tissue":
                tissue,

            "Inheritance_class":
                cls,

            "N":
                int(n),

            "Comparable_total":
                total,

            "Percent":
                100.0 * n / total
        })


    # ========================================================
    # Pattern summary
    # ========================================================

    pt = (
        comp[
            "Pattern"
        ]
        .value_counts()
    )


    for pat, n in pt.items():

        pattern_rows.append({
            "Tissue":
                tissue,

            "Pattern":
                pat,

            "N":
                int(n),

            "Comparable_total":
                total,

            "Percent":
                100.0 * n / total
        })


    # ========================================================
    # Transgressive summary
    # ========================================================

    up = comp[
        (
            comp["Transgressive"] == "YES"
        )
        &
        (
            comp["Transgressive_direction"] == "Up"
        )
    ]


    down = comp[
        (
            comp["Transgressive"] == "YES"
        )
        &
        (
            comp["Transgressive_direction"] == "Down"
        )
    ]


    strong_trans = comp[
        comp[
            "Strong_transgressive"
        ] == "YES"
    ]


    trans_rows.append({
        "Tissue":
            tissue,

        "Transgressive_up_N":
            len(up),

        "Transgressive_up_percent":
            100.0 * len(up) / total,

        "Transgressive_down_N":
            len(down),

        "Transgressive_down_percent":
            100.0 * len(down) / total,

        "Transgressive_total_N":
            len(up) + len(down),

        "Transgressive_total_percent":
            100.0 * (
                len(up) + len(down)
            ) / total,

        "Strong_transgressive_N":
            len(strong_trans)
    })


    # ========================================================
    # No-change
    # ========================================================

    nc = (
        comp[
            "No_change"
        ]
        .value_counts()
    )


    for state, n in nc.items():

        nochange_rows.append({
            "Tissue":
                tissue,

            "No_change":
                state,

            "N":
                int(n),

            "Comparable_total":
                total,

            "Percent":
                100.0 * n / total
        })


    # ========================================================
    # QC
    # ========================================================

    trans = comp[
        comp["Transgressive"] == "YES"
    ]


    trans_fail = int(
        (
            ~trans[
                "BT_BSB_significant"
            ]
            |
            ~trans[
                "BT_TC_significant"
            ]
        ).sum()
    )


    unresolved = int(
        (
            comp[
                "Inheritance_class"
            ]
            ==
            "Unresolved_tie"
        ).sum()
    )


    status = (
        "PASS"
        if trans_fail == 0
        else "FAIL"
    )


    qc_rows.append({
        "Tissue":
            tissue,

        "Comparable_total":
            total,

        "Assigned_total":
            len(comp),

        "Transgressive_total":
            len(trans),

        "Transgressive_rule_fail":
            trans_fail,

        "Unresolved_tie_N":
            unresolved,

        "QC_status":
            status
    })


    print("")
    print(
        ct.to_string()
    )

    print("")
    print(
        "Transgressive up   = {}".format(
            len(up)
        )
    )

    print(
        "Transgressive down = {}".format(
            len(down)
        )
    )

    print(
        "QC = {}".format(
            status
        )
    )

    print("")


# ============================================================
# Combined table
# ============================================================

combined = pd.concat(
    all_tables,
    ignore_index=True
)


combined.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_inheritance_all_tissues.tsv"
    ),
    sep="\t",
    index=False
)


combined.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_12pattern_all_tissues.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Main summaries
# ============================================================

class_summary = pd.DataFrame(
    class_rows
)

class_summary.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_inheritance_summary.tsv"
    ),
    sep="\t",
    index=False
)


pattern_summary = pd.DataFrame(
    pattern_rows
)

pattern_summary.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_12pattern_summary.tsv"
    ),
    sep="\t",
    index=False
)


trans_summary = pd.DataFrame(
    trans_rows
)

trans_summary.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_transgressive_summary.tsv"
    ),
    sep="\t",
    index=False
)


# keep old filename compatible
trans_summary.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_strong_transgressive_summary.tsv"
    ),
    sep="\t",
    index=False
)


nochange_summary = pd.DataFrame(
    nochange_rows
)

nochange_summary.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_no_change_summary.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Main class x No_change
# ============================================================

comp_all = combined[
    combined[
        "Comparable_all3"
    ] == 1
].copy()


class_nc = (
    comp_all
    .groupby(
        [
            "Tissue",
            "Inheritance_class",
            "No_change"
        ]
    )
    .size()
    .reset_index(
        name="N"
    )
)


class_nc.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_inheritance_by_no_change.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Broad class
# ============================================================

broad = (
    comp_all
    .groupby(
        [
            "Tissue",
            "Broad_class",
            "Inheritance_class"
        ]
    )
    .size()
    .reset_index(
        name="N"
    )
)


broad.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_broadclass_inheritance_summary.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Transgressive PanTE list
# ============================================================

transgressive_list = comp_all[
    comp_all[
        "Transgressive"
    ] == "YES"
].copy()


transgressive_list.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_transgressive_supported.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# Settings
# ============================================================

settings = pd.DataFrame([
    {
        "Parameter":
            "Main_classes",

        "Value":
            (
                "Additive;TC_expression_dominance;"
                "BSB_expression_dominance;"
                "Above_parental_range;"
                "Below_parental_range"
            )
    },
    {
        "Parameter":
            "Dominance_FDR",

        "Value":
            FDR
    },
    {
        "Parameter":
            "Transgressive_rule",

        "Value":
            (
                "BT outside parental range AND "
                "BT_vs_BSB padj<0.05 AND "
                "BT_vs_TC padj<0.05"
            )
    },
    {
        "Parameter":
            "Main_log2FC_cutoff",

        "Value":
            "NONE"
    },
    {
        "Parameter":
            "Strong_transgressive",

        "Value":
            (
                "Transgressive=YES and both "
                "|BT-parent log2FC|>=1"
            )
    },
    {
        "Parameter":
            "No_change",

        "Value":
            (
                "all three pairwise padj>=0.05"
            )
    }
])


settings.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_inheritance_settings.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# QC
# ============================================================

qc = pd.DataFrame(
    qc_rows
)


qc.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_inheritance_classification_QC.tsv"
    ),
    sep="\t",
    index=False
)


qc.to_csv(
    os.path.join(
        OUTDIR,
        "PanTE_12pattern_classification_QC.tsv"
    ),
    sep="\t",
    index=False
)


if (
    qc[
        "QC_status"
    ] != "PASS"
).any():

    raise RuntimeError(
        "Classification QC failed."
    )


print("")
print(
    "=============================================="
)
print(
    "FINAL MAIN CLASS SUMMARY"
)
print(
    "=============================================="
)

print(
    class_summary.to_string(
        index=False
    )
)


print("")
print(
    "=============================================="
)
print(
    "TRANSGRESSIVE SUMMARY"
)
print(
    "=============================================="
)

print(
    trans_summary.to_string(
        index=False
    )
)


print("")
print(
    "=============================================="
)
print(
    "QC"
)
print(
    "=============================================="
)

print(
    qc.to_string(
        index=False
    )
)


print("")
print(
    "FINAL classification completed successfully."
)