#!/usr/bin/env python
# Public-release note: original executed source; replace local absolute paths with local paths before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import glob
import os
import pandas as pd


BASE = os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52")
OUT = os.path.join(BASE, "06_TE_WGBS")

files = sorted(
    glob.glob(
        os.path.join(
            OUT,
            "02_quantification/*.PanTE_weighted_mCG.tsv"
        )
    )
)

if len(files) == 0:
    raise RuntimeError(
        "No PanTE weighted mCG files found."
    )


all_df = []

for f in files:

    df = pd.read_csv(
        f,
        sep="\t"
    )

    all_df.append(
        df
    )


x = pd.concat(
    all_df,
    ignore_index=True
)


print("")
print("===== RESULT FILE COUNT =====")
print(len(files))


print("")
print("===== TOTAL SAMPLE-PANTE-STATE RECORDS =====")
print(len(x))


print("")
print("===== N_CpG_sites DISTRIBUTION =====")

for q in [
    0,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99,
    1.0
]:

    print(
        "q={:.2f}: {:.1f}".format(
            q,
            x["N_CpG_sites"].quantile(q)
        )
    )


print("")
print("===== Total_C_reads DISTRIBUTION =====")

for q in [
    0,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99,
    1.0
]:

    print(
        "q={:.2f}: {:.1f}".format(
            q,
            x["Total_C_reads"].quantile(q)
        )
    )


print("")
print("===== RECORDS PASSING CpG THRESHOLDS =====")

for n in [
    1,
    2,
    3,
    5,
    10,
    20
]:

    k = (
        x["N_CpG_sites"]
        >=
        n
    ).sum()

    print(
        "N_CpG >= {:2d}: {:8d} ({:.2f}%)".format(
            n,
            k,
            100.0 * k / len(x)
        )
    )


print("")
print("===== RECORDS PASSING TOTAL-READ THRESHOLDS =====")

for n in [
    1,
    5,
    10,
    20,
    50,
    100
]:

    k = (
        x["Total_C_reads"]
        >=
        n
    ).sum()

    print(
        "Total_C_reads >= {:3d}: {:8d} ({:.2f}%)".format(
            n,
            k,
            100.0 * k / len(x)
        )
    )


print("")
print("===== UNIQUE PanTE BY GENOME STATE =====")

print(
    x.groupby(
        "Genome_state"
    )["PanTE_ID"]
    .nunique()
)


print("")
print("===== MEDIAN N_CpG BY GENOME STATE =====")

print(
    x.groupby(
        "Genome_state"
    )["N_CpG_sites"]
    .median()
)


print("")
print("===== MEDIAN TOTAL READS BY GENOME STATE =====")

print(
    x.groupby(
        "Genome_state"
    )["Total_C_reads"]
    .median()
)


print("")
print("===== RECORDS BY GROUP / TISSUE / STATE =====")

print(
    x.groupby(
        [
            "Group",
            "Tissue",
            "Genome_state"
        ]
    ).size()
)