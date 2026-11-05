#!/usr/bin/env python
# Public-release note: source script from the executed analysis. Replace the original local absolute paths and executable locations with paths in your environment before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import glob
import os
import pandas as pd


BASE = os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52")
ATAC = os.path.join(
    BASE,
    "07_TE_ATAC"
)

INDIR = os.path.join(
    ATAC,
    "00_QC/mapq_per_sample"
)

OUT = os.path.join(
    ATAC,
    "00_QC/ATAC_mapping_quality_summary.tsv"
)


files = sorted(
    glob.glob(
        os.path.join(
            INDIR,
            "*.mapq_qc.tsv"
        )
    )
)


if len(files) != 27:
    raise RuntimeError(
        "Expected 27 QC files, found {}".format(
            len(files)
        )
    )


dfs = []

for f in files:

    x = pd.read_csv(
        f,
        sep="\t"
    )

    if len(x) != 1:
        raise RuntimeError(
            "Expected exactly 1 data row in {}, found {}".format(
                f,
                len(x)
            )
        )

    dfs.append(
        x
    )


df = pd.concat(
    dfs,
    ignore_index=True
)


tissue_order = {
    "gill": 0,
    "kidney": 1,
    "muscle": 2
}

group_order = {
    "BSB": 0,
    "TC": 1,
    "BT": 2
}


df["_group_order"] = df["Group"].map(
    group_order
)

df["_tissue_order"] = df["Tissue"].map(
    tissue_order
)


df = df.sort_values(
    [
        "_group_order",
        "Sample",
        "_tissue_order"
    ]
).drop(
    columns=[
        "_group_order",
        "_tissue_order"
    ]
)


df.to_csv(
    OUT,
    sep="\t",
    index=False
)


print("")
print("===== FILE COUNT =====")
print(len(files))

print("")
print("===== RECORD COUNT =====")
print(len(df))

print("")
print("===== FULL SUMMARY =====")
print(
    df.to_string(
        index=False
    )
)


print("")
print("===== Q30 FRACTION SUMMARY =====")

print(
    df.groupby(
        [
            "Group",
            "Tissue"
        ]
    )[
        "Q30_fraction"
    ].agg(
        [
            "min",
            "median",
            "max"
        ]
    )
)


print("")
print("===== SAMPLES WITH Q30 < Q0 =====")

different = df[
    df[
        "Proper_clean_fragments_Q30"
    ]
    <
    df[
        "Proper_clean_fragments_Q0"
    ]
]

if len(different) == 0:

    print(
        "NONE"
    )

else:

    print(
        different[
            [
                "Sample",
                "Group",
                "Tissue",
                "Proper_clean_fragments_Q0",
                "Proper_clean_fragments_Q30",
                "Q30_fraction"
            ]
        ].to_string(
            index=False
        )
    )


print("")
print("Saved:")
print(OUT)