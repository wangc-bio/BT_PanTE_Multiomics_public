#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import numpy as np
import pandas as pd


# ============================================================
# PATHS
# ============================================================

KIM = os.environ.get(
    "KIMURA_DIR",
    os.path.join(
        os.environ.get(
            "KIMURA_BASE",
            "PATH_TO_RepeatMasker_unified_v5_2_Kimura_align"
        ),
        "Kimura_divergence"
    )
)

PARSED = os.path.join(
    KIM,
    "01_parsed"
)

OUT = os.path.join(
    KIM,
    "02_landscape_analysis"
)

if not os.path.exists(OUT):
    os.makedirs(OUT)


COVERAGE_FILE = os.path.join(
    PARSED,
    "Kimura_coverage_main5.long.tsv"
)

SIZE_FILE = os.path.join(
    PARSED,
    "assembly_sizes.tsv"
)


# ============================================================
# SETTINGS
# ============================================================

SAMPLES = [
    "BSB",
    "TC",
    "BTM",
    "BTP"
]

CLASSES = [
    "DNA",
    "LTR",
    "LINE",
    "SINE",
    "RC"
]


# ============================================================
# HELPERS
# ============================================================

def spearman_no_scipy(x, y):

    x = pd.Series(
        x,
        dtype=float
    )

    y = pd.Series(
        y,
        dtype=float
    )

    valid = (
        x.notna()
        &
        y.notna()
    )

    x = x[
        valid
    ]

    y = y[
        valid
    ]

    if len(x) < 3:
        return np.nan

    rx = x.rank(
        method="average"
    )

    ry = y.rank(
        method="average"
    )

    if (
        rx.std() == 0
        or
        ry.std() == 0
    ):
        return np.nan

    return float(
        np.corrcoef(
            rx,
            ry
        )[0, 1]
    )


def pearson_no_scipy(x, y):

    x = pd.Series(
        x,
        dtype=float
    )

    y = pd.Series(
        y,
        dtype=float
    )

    valid = (
        x.notna()
        &
        y.notna()
    )

    x = x[
        valid
    ].values

    y = y[
        valid
    ].values

    if len(x) < 3:
        return np.nan

    if (
        np.std(x) == 0
        or
        np.std(y) == 0
    ):
        return np.nan

    return float(
        np.corrcoef(
            x,
            y
        )[0, 1]
    )


def linear_slope(x, y):

    x = np.asarray(
        x,
        dtype=float
    )

    y = np.asarray(
        y,
        dtype=float
    )

    valid = (
        np.isfinite(x)
        &
        np.isfinite(y)
    )

    x = x[
        valid
    ]

    y = y[
        valid
    ]

    if len(x) < 3:
        return np.nan

    if np.std(x) == 0:
        return np.nan

    return float(
        np.polyfit(
            x,
            y,
            1
        )[0]
    )


# ============================================================
# LOAD DATA
# ============================================================

print("")
print("============================================")
print("Loading Kimura coverage")
print("============================================")

d = pd.read_csv(
    COVERAGE_FILE,
    sep="\t"
)

sizes = pd.read_csv(
    SIZE_FILE,
    sep="\t"
)


sizes[
    "Assembly_bp"
] = pd.to_numeric(
    sizes[
        "Assembly_bp"
    ],
    errors="raise"
)


size_map = dict(
    zip(
        sizes[
            "Sample"
        ],
        sizes[
            "Assembly_bp"
        ]
    )
)


for sample in SAMPLES:

    if sample not in size_map:

        raise RuntimeError(
            "Missing assembly size for {}".format(
                sample
            )
        )


# ============================================================
# QC INPUT
# ============================================================

expected_rows = (
    len(
        SAMPLES
    )
    *
    len(
        CLASSES
    )
    *
    71
)


if len(d) != expected_rows:

    print(
        "[WARNING] main5 matrix has {} rows; expected {}.".format(
            len(d),
            expected_rows
        )
    )


# ============================================================
# GENOME NORMALIZATION
#
# % assembly =
# 100 * well-characterized TE bp / assembly length
# ============================================================

d[
    "Assembly_bp"
] = d[
    "Sample"
].map(
    size_map
)


d[
    "Genome_pct"
] = (
    100.0
    *
    d[
        "Coverage_bp_wellChar"
    ]
    /
    d[
        "Assembly_bp"
    ]
)


# ============================================================
# TOTAL COVERAGE BY CLASS
# ============================================================

class_total = (
    d
    .groupby(
        [
            "Sample",
            "Broad_class"
        ],
        as_index=False
    )
    .agg(
        Coverage_bp_wellChar=(
            "Coverage_bp_wellChar",
            "sum"
        ),
        Genome_pct=(
            "Genome_pct",
            "sum"
        )
    )
)


# ============================================================
# TOTAL MAIN5 TE COVERAGE
# ============================================================

main5_total = (
    d
    .groupby(
        "Sample",
        as_index=False
    )
    .agg(
        Main5_TE_bp=(
            "Coverage_bp_wellChar",
            "sum"
        ),
        Main5_TE_genome_pct=(
            "Genome_pct",
            "sum"
        )
    )
)


# ============================================================
# MAKE WIDE MATRIX
# ============================================================

wide = (
    d
    .pivot_table(
        index=[
            "Broad_class",
            "Kimura_bin"
        ],
        columns="Sample",
        values="Genome_pct",
        aggfunc="sum"
    )
    .reset_index()
)


for sample in SAMPLES:

    if sample not in wide.columns:

        raise RuntimeError(
            "Missing {} in wide matrix.".format(
                sample
            )
        )


# ============================================================
# PARENTAL AND BT CONTRASTS
# ============================================================

wide[
    "Parent_delta_pct"
] = (
    wide[
        "BSB"
    ]
    -
    wide[
        "TC"
    ]
)


wide[
    "BT_delta_pct"
] = (
    wide[
        "BTM"
    ]
    -
    wide[
        "BTP"
    ]
)


wide[
    "BSB_to_BTM_delta_pct"
] = (
    wide[
        "BTM"
    ]
    -
    wide[
        "BSB"
    ]
)


wide[
    "TC_to_BTP_delta_pct"
] = (
    wide[
        "BTP"
    ]
    -
    wide[
        "TC"
    ]
)


wide[
    "Abs_parent_delta_pct"
] = np.abs(
    wide[
        "Parent_delta_pct"
    ]
)


wide[
    "Abs_BT_delta_pct"
] = np.abs(
    wide[
        "BT_delta_pct"
    ]
)


# ============================================================
# CLASS-SPECIFIC RETENTION STATISTICS
# ============================================================

stats_rows = []


for cls in CLASSES:

    x = wide[
        wide[
            "Broad_class"
        ] == cls
    ].copy()


    x = x.sort_values(
        "Kimura_bin"
    )


    stats_rows.append({
        "Broad_class":
            cls,

        "N_bins":
            len(x),

        "BSB_vs_BTM_Spearman":
            spearman_no_scipy(
                x[
                    "BSB"
                ],
                x[
                    "BTM"
                ]
            ),

        "TC_vs_BTP_Spearman":
            spearman_no_scipy(
                x[
                    "TC"
                ],
                x[
                    "BTP"
                ]
            ),

        "ParentContrast_vs_BTContrast_Spearman":
            spearman_no_scipy(
                x[
                    "Parent_delta_pct"
                ],
                x[
                    "BT_delta_pct"
                ]
            ),

        "ParentContrast_vs_BTContrast_Pearson":
            pearson_no_scipy(
                x[
                    "Parent_delta_pct"
                ],
                x[
                    "BT_delta_pct"
                ]
            ),

        "ParentContrast_to_BTContrast_slope":
            linear_slope(
                x[
                    "Parent_delta_pct"
                ],
                x[
                    "BT_delta_pct"
                ]
            ),

        "Parent_absContrast_sum":
            float(
                x[
                    "Abs_parent_delta_pct"
                ].sum()
            ),

        "BT_absContrast_sum":
            float(
                x[
                    "Abs_BT_delta_pct"
                ].sum()
            )
    })


# ============================================================
# ALL MAIN5, SUMMED BY KIMURA BIN
#
# Important:
# This is not pooled 5*71 observations.
# Classes are first summed within each Kimura bin.
# ============================================================

all_bins = (
    wide
    .groupby(
        "Kimura_bin",
        as_index=False
    )[
        [
            "BSB",
            "TC",
            "BTM",
            "BTP"
        ]
    ]
    .sum()
)


all_bins[
    "Parent_delta_pct"
] = (
    all_bins[
        "BSB"
    ]
    -
    all_bins[
        "TC"
    ]
)


all_bins[
    "BT_delta_pct"
] = (
    all_bins[
        "BTM"
    ]
    -
    all_bins[
        "BTP"
    ]
)


all_bins[
    "Abs_parent_delta_pct"
] = np.abs(
    all_bins[
        "Parent_delta_pct"
    ]
)


all_bins[
    "Abs_BT_delta_pct"
] = np.abs(
    all_bins[
        "BT_delta_pct"
    ]
)


stats_rows.append({
    "Broad_class":
        "All_main5",

    "N_bins":
        len(
            all_bins
        ),

    "BSB_vs_BTM_Spearman":
        spearman_no_scipy(
            all_bins[
                "BSB"
            ],
            all_bins[
                "BTM"
            ]
        ),

    "TC_vs_BTP_Spearman":
        spearman_no_scipy(
            all_bins[
                "TC"
            ],
            all_bins[
                "BTP"
            ]
        ),

    "ParentContrast_vs_BTContrast_Spearman":
        spearman_no_scipy(
            all_bins[
                "Parent_delta_pct"
            ],
            all_bins[
                "BT_delta_pct"
            ]
        ),

    "ParentContrast_vs_BTContrast_Pearson":
        pearson_no_scipy(
            all_bins[
                "Parent_delta_pct"
            ],
            all_bins[
                "BT_delta_pct"
            ]
        ),

    "ParentContrast_to_BTContrast_slope":
        linear_slope(
            all_bins[
                "Parent_delta_pct"
            ],
            all_bins[
                "BT_delta_pct"
            ]
        ),

    "Parent_absContrast_sum":
        float(
            all_bins[
                "Abs_parent_delta_pct"
            ].sum()
        ),

    "BT_absContrast_sum":
        float(
            all_bins[
                "Abs_BT_delta_pct"
            ].sum()
        )
})


stats = pd.DataFrame(
    stats_rows
)


stats[
    "AbsContrast_ratio_BT_over_parent"
] = (
    stats[
        "BT_absContrast_sum"
    ]
    /
    stats[
        "Parent_absContrast_sum"
    ]
)


# ============================================================
# DOMINANT KIMURA PEAKS
#
# Purely descriptive:
# identify bin with maximum genome coverage for each
# sample x broad class.
# ============================================================

peak_rows = []


for sample in SAMPLES:

    for cls in CLASSES:

        x = d[
            (
                d[
                    "Sample"
                ] == sample
            )
            &
            (
                d[
                    "Broad_class"
                ] == cls
            )
        ].copy()


        if x.empty:
            continue


        idx = x[
            "Genome_pct"
        ].idxmax()


        row = x.loc[
            idx
        ]


        peak_rows.append({
            "Sample":
                sample,

            "Broad_class":
                cls,

            "Peak_Kimura_bin":
                int(
                    row[
                        "Kimura_bin"
                    ]
                ),

            "Peak_genome_pct":
                float(
                    row[
                        "Genome_pct"
                    ]
                )
        })


peaks = pd.DataFrame(
    peak_rows
)


# ============================================================
# LTR-SPECIFIC TABLE
# ============================================================

ltr = wide[
    wide[
        "Broad_class"
    ] == "LTR"
].copy()


ltr = ltr.sort_values(
    "Kimura_bin"
)


# ============================================================
# OUTPUT
# ============================================================

LONG_OUT = os.path.join(
    OUT,
    "Kimura_landscape_main5.genome_normalized.long.tsv"
)

WIDE_OUT = os.path.join(
    OUT,
    "Kimura_landscape_main5.genome_normalized.wide.tsv"
)

TOTAL_OUT = os.path.join(
    OUT,
    "Kimura_main5_total_genome_coverage.tsv"
)

CLASS_OUT = os.path.join(
    OUT,
    "Kimura_class_total_genome_coverage.tsv"
)

STATS_OUT = os.path.join(
    OUT,
    "Kimura_parental_landscape_retention_statistics.tsv"
)

ALLBIN_OUT = os.path.join(
    OUT,
    "Kimura_all_main5_by_bin.tsv"
)

LTR_OUT = os.path.join(
    OUT,
    "Kimura_LTR_parental_contrast_by_bin.tsv"
)

PEAK_OUT = os.path.join(
    OUT,
    "Kimura_peak_bins.tsv"
)

LOG_OUT = os.path.join(
    OUT,
    "Kimura_landscape_analysis_summary.txt"
)


d.to_csv(
    LONG_OUT,
    sep="\t",
    index=False
)


wide.to_csv(
    WIDE_OUT,
    sep="\t",
    index=False
)


main5_total.to_csv(
    TOTAL_OUT,
    sep="\t",
    index=False
)


class_total.to_csv(
    CLASS_OUT,
    sep="\t",
    index=False
)


stats.to_csv(
    STATS_OUT,
    sep="\t",
    index=False
)


all_bins.to_csv(
    ALLBIN_OUT,
    sep="\t",
    index=False
)


ltr.to_csv(
    LTR_OUT,
    sep="\t",
    index=False
)


peaks.to_csv(
    PEAK_OUT,
    sep="\t",
    index=False
)


with open(
    LOG_OUT,
    "w"
) as out:

    out.write(
        "Kimura landscape analysis summary\n"
    )

    out.write(
        "=================================\n\n"
    )

    out.write(
        "Assembly sizes\n"
    )

    out.write(
        "--------------\n"
    )

    out.write(
        sizes.to_string(
            index=False
        )
    )

    out.write(
        "\n\n"
    )


    out.write(
        "Main five TE classes: total genome coverage\n"
    )

    out.write(
        "-------------------------------------------\n"
    )

    out.write(
        main5_total.to_string(
            index=False
        )
    )

    out.write(
        "\n\n"
    )


    out.write(
        "Class-specific genome coverage\n"
    )

    out.write(
        "------------------------------\n"
    )

    out.write(
        class_total.to_string(
            index=False
        )
    )

    out.write(
        "\n\n"
    )


    out.write(
        "Parental landscape retention statistics\n"
    )

    out.write(
        "---------------------------------------\n"
    )

    out.write(
        stats.to_string(
            index=False
        )
    )

    out.write(
        "\n\n"
    )


    out.write(
        "Peak Kimura bins\n"
    )

    out.write(
        "----------------\n"
    )

    out.write(
        peaks.to_string(
            index=False
        )
    )

    out.write(
        "\n"
    )


print("")
print("============================================")
print("Kimura landscape analysis completed")
print("============================================")

print("")
print("Main5 genome coverage:")
print(
    main5_total.to_string(
        index=False
    )
)

print("")
print("Retention statistics:")
print(
    stats.to_string(
        index=False
    )
)

print("")
print("Peak bins:")
print(
    peaks.to_string(
        index=False
    )
)

print("")
print("Summary:")
print(
    LOG_OUT
)
