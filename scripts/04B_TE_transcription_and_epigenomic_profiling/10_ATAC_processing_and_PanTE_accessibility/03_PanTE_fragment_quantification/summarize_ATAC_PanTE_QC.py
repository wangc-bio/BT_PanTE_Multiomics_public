#!/usr/bin/env python
# Public-release note: source script from the executed analysis. Replace the original local absolute paths and executable locations with paths in your environment before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import glob
import numpy as np
import pandas as pd


BASE = os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52")
ATAC = os.path.join(BASE, "07_TE_ATAC")

COUNTDIR = os.path.join(
    ATAC,
    "02_PanTE_count"
)

OUTDIR = os.path.join(
    ATAC,
    "03_matrix"
)

os.makedirs(
    OUTDIR,
    exist_ok=True
)


# ============================================================
# 1. ACCESSIBILITY FILES
# ============================================================

files = sorted(
    glob.glob(
        os.path.join(
            COUNTDIR,
            "*.PanTE_accessibility.tsv"
        )
    )
)

print("")
print("============================================")
print("ACCESSIBILITY FILE CHECK")
print("============================================")
print("Files:", len(files))

if len(files) != 27:
    raise RuntimeError(
        "Expected 27 accessibility files, found {}".format(
            len(files)
        )
    )


dfs = []

for f in files:

    x = pd.read_csv(
        f,
        sep="\t"
    )

    dfs.append(x)


dat = pd.concat(
    dfs,
    ignore_index=True
)


LONG_OUT = os.path.join(
    OUTDIR,
    "PanTE_ATAC_accessibility_long.tsv"
)

dat.to_csv(
    LONG_OUT,
    sep="\t",
    index=False
)


print("Long records:", len(dat))
print("Samples:", dat["Sample"].nunique())


# ============================================================
# 2. SAMPLE / STATE CHECK
# ============================================================

rep = (
    dat[
        [
            "Sample",
            "Group",
            "Tissue",
            "Genome_state"
        ]
    ]
    .drop_duplicates()
    .groupby(
        [
            "Genome_state",
            "Tissue"
        ]
    )
    .size()
    .reset_index(
        name="N_replicates"
    )
)

print("")
print("============================================")
print("REPLICATES PER STATE / TISSUE")
print("============================================")
print(rep.to_string(index=False))


# ============================================================
# 3. MAIN-CHROMOSOME FRACTION
# ============================================================

mapq = pd.read_csv(
    os.path.join(
        ATAC,
        "00_QC/ATAC_mapping_quality_summary.tsv"
    ),
    sep="\t"
)


total_files = sorted(
    glob.glob(
        os.path.join(
            COUNTDIR,
            "*.state_total_fragments.tsv"
        )
    )
)


total_rows = []

for f in total_files:

    name = os.path.basename(f)

    prefix = name.replace(
        ".state_total_fragments.tsv",
        ""
    )

    sample, tissue = prefix.split(
        "_",
        1
    )

    x = pd.read_csv(
        f,
        sep="\t"
    )

    x["Sample"] = sample
    x["Tissue"] = tissue

    total_rows.append(x)


tot = pd.concat(
    total_rows,
    ignore_index=True
)


tot_sum = (
    tot.groupby(
        [
            "Sample",
            "Tissue"
        ],
        as_index=False
    )[
        "Total_fragments"
    ]
    .sum()
    .rename(
        columns={
            "Total_fragments":
            "Main_chr_fragments"
        }
    )
)


main = pd.merge(
    mapq,
    tot_sum,
    on=[
        "Sample",
        "Tissue"
    ],
    how="left"
)


main[
    "Main_chr_fraction"
] = (
    main["Main_chr_fragments"]
    /
    main["Proper_clean_fragments_Q0"]
)


MAIN_OUT = os.path.join(
    OUTDIR,
    "ATAC_main_chromosome_fragment_QC.tsv"
)

main.to_csv(
    MAIN_OUT,
    sep="\t",
    index=False
)


print("")
print("============================================")
print("MAIN-CHROMOSOME FRAGMENT FRACTION")
print("============================================")

print(
    main.groupby(
        [
            "Group",
            "Tissue"
        ]
    )[
        "Main_chr_fraction"
    ].agg(
        [
            "min",
            "median",
            "max"
        ]
    ).to_string()
)


# ============================================================
# 4. UNIQUE ASSIGNMENT QC
# ============================================================

qc_files = sorted(
    glob.glob(
        os.path.join(
            COUNTDIR,
            "*.PanTE_overlap_QC.tsv"
        )
    )
)


qc_rows = []

for f in qc_files:

    name = os.path.basename(f)

    prefix = name.replace(
        ".PanTE_overlap_QC.tsv",
        ""
    )

    sample, tissue = prefix.split(
        "_",
        1
    )

    x = pd.read_csv(
        f,
        sep="\t"
    )

    x["Sample"] = sample
    x["Tissue"] = tissue

    qc_rows.append(x)


qc = pd.concat(
    qc_rows,
    ignore_index=True
)


print("")
print("============================================")
print("UNIQUE PanTE ASSIGNMENT")
print("============================================")

print(
    qc.groupby(
        "Genome_state"
    )[
        "Unique_percent"
    ].agg(
        [
            "min",
            "median",
            "max"
        ]
    ).to_string()
)


# ============================================================
# 5. ASSIGNABLE BP DISTRIBUTION
# ============================================================

bp = (
    dat[
        [
            "Genome_state",
            "PanTE_ID",
            "Assignable_bp",
            "N_segments"
        ]
    ]
    .drop_duplicates()
)


quantiles = [
    0,
    0.01,
    0.05,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99,
    1.00
]


print("")
print("============================================")
print("ASSIGNABLE BP QUANTILES")
print("============================================")

for state in [
    "BSB",
    "TC",
    "BTM",
    "BTP"
]:

    x = bp.loc[
        bp["Genome_state"] == state,
        "Assignable_bp"
    ].values

    print("")
    print("----- {} -----".format(state))
    print("PanTE:", len(x))

    for q in quantiles:

        print(
            "q{:g}\t{:.1f}".format(
                q,
                np.quantile(x, q)
            )
        )


# ============================================================
# 6. POOLED FRAGMENT COUNTS ACROSS 3 REPLICATES
# ============================================================

pooled = (
    dat.groupby(
        [
            "Genome_state",
            "Tissue",
            "PanTE_ID"
        ],
        as_index=False
    )[
        "Fragment_count"
    ]
    .sum()
    .rename(
        columns={
            "Fragment_count":
            "Pooled_fragment_count"
        }
    )
)


print("")
print("============================================")
print("POOLED FRAGMENT COUNT QUANTILES")
print("============================================")

for (state, tissue), x in pooled.groupby(
    [
        "Genome_state",
        "Tissue"
    ]
):

    v = x[
        "Pooled_fragment_count"
    ].values

    print("")
    print(
        "{} {}".format(
            state,
            tissue
        )
    )

    for q in [
        0,
        0.10,
        0.25,
        0.50,
        0.75,
        0.90,
        0.95,
        0.99,
        1.00
    ]:

        print(
            "q{:g}\t{:.1f}".format(
                q,
                np.quantile(v, q)
            )
        )


# ============================================================
# 7. REPLICATE SPEARMAN CORRELATIONS
# ============================================================

corr_rows = []


for (
    state,
    tissue
), x in dat.groupby(
    [
        "Genome_state",
        "Tissue"
    ]
):

    mat = x.pivot_table(
        index="PanTE_ID",
        columns="Sample",
        values="log2_Accessibility_FPKM1"
    )


    samples = list(
        mat.columns
    )


    for i in range(
        len(samples)
    ):

        for j in range(
            i + 1,
            len(samples)
        ):

            a = samples[i]
            b = samples[j]

            tmp = mat[
                [
                    a,
                    b
                ]
            ].dropna()


            rho = tmp[a].corr(
                tmp[b],
                method="spearman"
            )


            corr_rows.append({
                "Genome_state":
                    state,

                "Tissue":
                    tissue,

                "Sample1":
                    a,

                "Sample2":
                    b,

                "N_PanTE":
                    len(tmp),

                "Spearman_rho":
                    rho
            })


corr = pd.DataFrame(
    corr_rows
)


CORR_OUT = os.path.join(
    OUTDIR,
    "ATAC_replicate_Spearman.tsv"
)

corr.to_csv(
    CORR_OUT,
    sep="\t",
    index=False
)


print("")
print("============================================")
print("REPLICATE SPEARMAN")
print("============================================")

print(
    corr.groupby(
        [
            "Genome_state",
            "Tissue"
        ]
    )[
        "Spearman_rho"
    ].agg(
        [
            "min",
            "median",
            "max"
        ]
    ).to_string()
)


# ============================================================
# 8. COVERAGE SENSITIVITY
# ============================================================

# minimum fragment count across the 3 replicates
rep_counts = dat.pivot_table(
    index=[
        "Genome_state",
        "Tissue",
        "PanTE_ID",
        "Assignable_bp"
    ],
    columns="Sample",
    values="Fragment_count"
).reset_index()


meta_cols = [
    "Genome_state",
    "Tissue",
    "PanTE_ID",
    "Assignable_bp"
]

sample_cols = [
    c for c in rep_counts.columns
    if c not in meta_cols
]


rep_counts[
    "Min_rep_fragment_count"
] = rep_counts[
    sample_cols
].min(
    axis=1,
    skipna=True
)


bp_cutoffs = [
    0,
    100,
    500,
    1000,
    2000,
    5000,
    10000
]

count_cutoffs = [
    0,
    1,
    3,
    5,
    10
]


sensitivity_rows = []


for state in [
    "BSB",
    "TC",
    "BTM",
    "BTP"
]:

    for tissue in [
        "gill",
        "kidney",
        "muscle"
    ]:

        x = rep_counts[
            (rep_counts["Genome_state"] == state)
            &
            (rep_counts["Tissue"] == tissue)
        ]

        total_n = len(x)


        for bpc in bp_cutoffs:

            for cc in count_cutoffs:

                n = (
                    (
                        x["Assignable_bp"] >= bpc
                    )
                    &
                    (
                        x[
                            "Min_rep_fragment_count"
                        ] >= cc
                    )
                ).sum()


                sensitivity_rows.append({
                    "Genome_state":
                        state,

                    "Tissue":
                        tissue,

                    "Assignable_bp_cutoff":
                        bpc,

                    "Min_fragment_each_rep_cutoff":
                        cc,

                    "Total_PanTE":
                        total_n,

                    "Pass_PanTE":
                        n,

                    "Pass_percent":
                        100.0 * n / total_n
                        if total_n > 0
                        else np.nan
                })


sens = pd.DataFrame(
    sensitivity_rows
)


SENS_OUT = os.path.join(
    OUTDIR,
    "ATAC_filter_sensitivity.tsv"
)

sens.to_csv(
    SENS_OUT,
    sep="\t",
    index=False
)


# ============================================================
# 9. COMMON PanTE ACROSS FOUR GENOME STATES
# ============================================================

sets = {}

for state in [
    "BSB",
    "TC",
    "BTM",
    "BTP"
]:

    sets[state] = set(
        bp.loc[
            bp["Genome_state"] == state,
            "PanTE_ID"
        ]
    )


common = (
    sets["BSB"]
    &
    sets["TC"]
    &
    sets["BTM"]
    &
    sets["BTP"]
)


print("")
print("============================================")
print("COMMON PanTE ACROSS FOUR STATES")
print("============================================")

for state in [
    "BSB",
    "TC",
    "BTM",
    "BTP"
]:

    print(
        "{}\t{}".format(
            state,
            len(sets[state])
        )
    )

print(
    "Four-state common\t{}".format(
        len(common)
    )
)


with open(
    os.path.join(
        OUTDIR,
        "PanTE_four_state_common.txt"
    ),
    "w"
) as out:

    for p in sorted(common):
        out.write(
            p + "\n"
        )


print("")
print("============================================")
print("OUTPUT")
print("============================================")

print(LONG_OUT)
print(MAIN_OUT)
print(CORR_OUT)
print(SENS_OUT)

print("")
print("Completed successfully.")