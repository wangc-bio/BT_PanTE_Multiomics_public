#!/usr/bin/env python
# Public-release note: source script from the executed analysis. Replace the original local absolute paths and executable locations with paths in your environment before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import argparse
import numpy as np
import pandas as pd


parser = argparse.ArgumentParser()

parser.add_argument("--counts", required=True)
parser.add_argument("--bp", required=True)
parser.add_argument("--total", required=True)

parser.add_argument("--sample", required=True)
parser.add_argument("--group", required=True)
parser.add_argument("--tissue", required=True)

parser.add_argument("--out", required=True)

args = parser.parse_args()


bp = pd.read_csv(
    args.bp,
    sep="\t"
)

counts = pd.read_csv(
    args.counts,
    sep="\t"
)

total = pd.read_csv(
    args.total,
    sep="\t"
)


states = {
    "BSB": ["BSB"],
    "TC": ["TC"],
    "BT": ["BTM", "BTP"]
}[args.group]


rows = []


for state in states:

    b = bp[
        bp["Genome_state"] == state
    ].copy()

    c = counts[
        counts["Genome_state"] == state
    ].copy()

    t = total[
        total["Genome_state"] == state
    ]


    if len(t) != 1:

        raise RuntimeError(
            "Expected exactly one total-fragment row for {}, got {}".format(
                state,
                len(t)
            )
        )


    libsize = int(
        t.iloc[0][
            "Total_fragments"
        ]
    )


    c = c[
        [
            "PanTE_ID",
            "Fragment_count"
        ]
    ]


    x = pd.merge(
        b,
        c,
        on="PanTE_ID",
        how="left"
    )


    x[
        "Fragment_count"
    ] = (
        x["Fragment_count"]
        .fillna(0)
        .astype(int)
    )


    x[
        "Total_state_fragments"
    ] = libsize


    # fragments / kb / million fragments
    x[
        "Accessibility_FPKM"
    ] = (
        x[
            "Fragment_count"
        ].astype(float)
        *
        1.0e9
        /
        (
            x[
                "Assignable_bp"
            ].astype(float)
            *
            float(libsize)
        )
    )


    x[
        "log2_Accessibility_FPKM1"
    ] = np.log2(
        x[
            "Accessibility_FPKM"
        ]
        +
        1.0
    )


    x.insert(
        0,
        "Sample",
        args.sample
    )

    x.insert(
        1,
        "Group",
        args.group
    )

    x.insert(
        2,
        "Tissue",
        args.tissue
    )


    rows.append(
        x
    )


result = pd.concat(
    rows,
    ignore_index=True
)


result.to_csv(
    args.out,
    sep="\t",
    index=False
)


print(
    "Sample:",
    args.sample
)

print(
    "Group:",
    args.group
)

print(
    "Tissue:",
    args.tissue
)

print(
    "Rows:",
    len(result)
)