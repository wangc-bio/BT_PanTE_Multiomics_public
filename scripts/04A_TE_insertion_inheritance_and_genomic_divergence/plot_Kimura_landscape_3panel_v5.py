#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch


# ============================================================
# PATHS
# ============================================================

KIM = os.environ.get("KIMURA_DIR", os.path.join(os.environ.get("KIMURA_BASE", "PATH_TO_RepeatMasker_unified_v5_2_Kimura_align"), "Kimura_divergence"))

ANA = os.path.join(
    KIM,
    "02_landscape_analysis"
)

FIG = os.path.join(
    KIM,
    "03_figure_v5"
)

os.makedirs(
    FIG,
    exist_ok=True
)


LONG_FILE = os.path.join(
    ANA,
    "Kimura_landscape_main5.genome_normalized.long.tsv"
)

PEAK_FILE = os.path.join(
    ANA,
    "Kimura_peak_bins.tsv"
)


for f in [
    LONG_FILE,
    PEAK_FILE
]:
    if not os.path.isfile(f):
        raise RuntimeError(
            "Missing input file: {}".format(f)
        )


# ============================================================
# SETTINGS
# ============================================================

XMIN = 0
XMAX = 50

BINS = np.arange(
    XMIN,
    XMAX + 1
)


CLASSES = [
    "DNA",
    "LTR",
    "LINE",
    "SINE",
    "RC"
]


# Panel A arrangement
PANEL_A_SAMPLES = [
    ["BSB", "BTM"],
    ["TC", "BTP"]
]


DISPLAY = {
    "BSB": "BSB",
    "BTM": "BT-M",
    "TC": "TC",
    "BTP": "BT-P"
}


# ============================================================
# FONT
# ============================================================

FONT_FILE = os.environ.get("FONT_FILE", "")


if os.path.isfile(
    FONT_FILE
):

    font_manager.fontManager.addfont(
        FONT_FILE
    )

    FONT_NAME = font_manager.FontProperties(
        fname=FONT_FILE
    ).get_name()

else:

    FONT_NAME = "Arial"


plt.rcParams.update({

    "font.family": FONT_NAME,
    "font.weight": "bold",

    "axes.labelweight": "bold",
    "axes.titleweight": "bold",

    "axes.linewidth": 1.25,

    "font.size": 10,

    "xtick.major.width": 1.15,
    "ytick.major.width": 1.15,

    "xtick.major.size": 4.5,
    "ytick.major.size": 4.5,

    "pdf.fonttype": 42,
    "ps.fonttype": 42,

    "svg.fonttype": "none"
})


# ============================================================
# COLORS
# ============================================================

CLASS_COLORS = {

    "DNA": "#4C78A8",

    "LTR": "#E45756",

    "LINE": "#72B7B2",

    "SINE": "#E3B341",

    "RC": "#B279A2"
}


# Panel C:
# same parental-origin color,
# line style distinguishes parent and BT haplotype.

ORIGIN_COLOR = {

    "BSB": "#B86500",
    "BTM": "#B86500",

    "TC": "#337AB7",
    "BTP": "#337AB7"
}


LINESTYLE = {

    "BSB": "-",
    "BTM": "--",

    "TC": "-",
    "BTP": "--"
}


# ============================================================
# LOAD DATA
# ============================================================

long = pd.read_csv(
    LONG_FILE,
    sep="\t"
)

peaks = pd.read_csv(
    PEAK_FILE,
    sep="\t"
)


required_long = [
    "Sample",
    "Kimura_bin",
    "Broad_class",
    "Genome_pct"
]


for c in required_long:

    if c not in long.columns:

        raise RuntimeError(
            "Missing column in long matrix: {}".format(
                c
            )
        )


required_peak = [
    "Sample",
    "Broad_class",
    "Peak_Kimura_bin"
]


for c in required_peak:

    if c not in peaks.columns:

        raise RuntimeError(
            "Missing column in peak table: {}".format(
                c
            )
        )


# ============================================================
# PLOTTING RANGE
# ============================================================

plot_long = long[
    (
        long[
            "Kimura_bin"
        ] >= XMIN
    )
    &
    (
        long[
            "Kimura_bin"
        ] <= XMAX
    )
].copy()


# ============================================================
# HELPERS
# ============================================================

def style_axis(
    ax,
    show_grid=True
):

    ax.spines[
        "top"
    ].set_visible(
        False
    )

    ax.spines[
        "right"
    ].set_visible(
        False
    )

    ax.spines[
        "left"
    ].set_linewidth(
        1.25
    )

    ax.spines[
        "bottom"
    ].set_linewidth(
        1.25
    )


    ax.tick_params(
        axis="both",
        direction="out",
        labelsize=9
    )


    if show_grid:

        ax.grid(
            axis="y",
            linewidth=0.55,
            alpha=0.28,
            zorder=0
        )


    ax.set_axisbelow(
        True
    )


    for t in (
        ax.get_xticklabels()
        +
        ax.get_yticklabels()
    ):

        t.set_fontweight(
            "bold"
        )


def save_figure(
    fig,
    basename
):

    pdf = os.path.join(
        FIG,
        basename + ".pdf"
    )

    svg = os.path.join(
        FIG,
        basename + ".svg"
    )

    png = os.path.join(
        FIG,
        basename + ".png"
    )


    fig.savefig(
        pdf,
        bbox_inches="tight"
    )


    fig.savefig(
        svg,
        bbox_inches="tight"
    )


    fig.savefig(
        png,
        dpi=600,
        bbox_inches="tight"
    )


def get_curve(
    sample,
    broad_class
):

    x = plot_long[
        (
            plot_long[
                "Sample"
            ] == sample
        )
        &
        (
            plot_long[
                "Broad_class"
            ] == broad_class
        )
    ][
        [
            "Kimura_bin",
            "Genome_pct"
        ]
    ].copy()


    s = (
        x
        .groupby(
            "Kimura_bin"
        )[
            "Genome_pct"
        ]
        .sum()
        .reindex(
            BINS,
            fill_value=0.0
        )
    )


    return s


def get_all_class_matrix(
    sample
):

    data = {}


    for cls in CLASSES:

        data[
            cls
        ] = get_curve(
            sample,
            cls
        ).values


    return pd.DataFrame(
        data,
        index=BINS
    )


def get_peak(
    sample
):

    x = peaks[
        (
            peaks[
                "Sample"
            ] == sample
        )
        &
        (
            peaks[
                "Broad_class"
            ] == "LTR"
        )
    ]


    if len(x) != 1:

        raise RuntimeError(
            "Cannot uniquely find LTR peak for {}".format(
                sample
            )
        )


    return int(
        x.iloc[
            0
        ][
            "Peak_Kimura_bin"
        ]
    )


# ============================================================
# PRECALCULATE
# ============================================================

STACK = {

    s:
        get_all_class_matrix(
            s
        )

    for s in [
        "BSB",
        "BTM",
        "TC",
        "BTP"
    ]
}


# ------------------------------------------------------------
# Panel A common y-axis maximum
# ------------------------------------------------------------

PANEL_A_YMAX = 0.0


for sample in STACK:

    m = STACK[
        sample
    ].sum(
        axis=1
    ).max()


    PANEL_A_YMAX = max(
        PANEL_A_YMAX,
        float(
            m
        )
    )


PANEL_A_YMAX *= 1.07


# ------------------------------------------------------------
# Panel B contrasts
# ------------------------------------------------------------

PARENT_DIFF = (
    STACK[
        "BSB"
    ]
    -
    STACK[
        "TC"
    ]
)


BT_DIFF = (
    STACK[
        "BTM"
    ]
    -
    STACK[
        "BTP"
    ]
)


def signed_stack_limit(
    matrix
):

    positive = (
        matrix
        .clip(
            lower=0
        )
        .sum(
            axis=1
        )
    )


    negative = (
        matrix
        .clip(
            upper=0
        )
        .sum(
            axis=1
        )
    )


    return max(
        float(
            positive.max()
        ),
        abs(
            float(
                negative.min()
            )
        )
    )


PANEL_B_YMAX = (
    max(
        signed_stack_limit(
            PARENT_DIFF
        ),
        signed_stack_limit(
            BT_DIFF
        )
    )
    *
    1.12
)


# ------------------------------------------------------------
# LTR peaks
# ------------------------------------------------------------

PEAKS = {

    "BSB":
        get_peak(
            "BSB"
        ),

    "BTM":
        get_peak(
            "BTM"
        ),

    "TC":
        get_peak(
            "TC"
        ),

    "BTP":
        get_peak(
            "BTP"
        )
}


# ============================================================
# PANEL A
# ============================================================

def draw_panel_A(
    fig,
    spec,
    panel_letter=False
):

    outer = fig.add_subplot(
        spec
    )


    outer.axis(
        "off"
    )


    gs = spec.subgridspec(
        2,
        2,

        wspace=0.18,
        hspace=0.24
    )


    for row in range(
        2
    ):

        for col in range(
            2
        ):

            sample = PANEL_A_SAMPLES[
                row
            ][
                col
            ]


            ax = fig.add_subplot(
                gs[
                    row,
                    col
                ]
            )


            matrix = STACK[
                sample
            ]


            bottom = np.zeros(
                len(
                    BINS
                )
            )


            for cls in CLASSES:

                values = matrix[
                    cls
                ].values


                ax.bar(
                    BINS,
                    values,

                    bottom=bottom,

                    width=0.92,

                    color=CLASS_COLORS[
                        cls
                    ],

                    edgecolor="none",
                    linewidth=0,

                    zorder=2
                )


                bottom += values


            ax.set_xlim(
                -0.5,
                XMAX + 0.5
            )


            ax.set_ylim(
                0,
                PANEL_A_YMAX
            )


            ax.set_xticks(
                np.arange(
                    0,
                    XMAX + 1,
                    10
                )
            )


            # ------------------------------------------------
            # Sample title
            # deliberately outside the plotting data region
            # ------------------------------------------------

            ax.set_title(
                DISPLAY[
                    sample
                ],

                loc="left",

                pad=5,

                fontsize=12,

                fontweight="bold"
            )


            # top row:
            # no duplicated x labels
            if row == 0:

                ax.tick_params(
                    axis="x",
                    labelbottom=False
                )


            # right column:
            # no duplicated y labels
            if col == 1:

                ax.tick_params(
                    axis="y",
                    labelleft=False
                )


            style_axis(
                ax,
                show_grid=True
            )


    # --------------------------------------------------------
    # Shared axis labels
    # --------------------------------------------------------

    outer.text(
        0.50,
        -0.073,

        "Kimura divergence (%)",

        transform=outer.transAxes,

        ha="center",
        va="top",

        fontsize=11,
        fontweight="bold"
    )


    outer.text(
        -0.060,
        0.50,

        "TE coverage (% assembly)",

        transform=outer.transAxes,

        ha="center",
        va="center",

        rotation=90,

        fontsize=11,
        fontweight="bold"
    )


    # --------------------------------------------------------
    # One shared TE-class legend
    #
    # moved upward relative to v4:
    # 1.055 -> 1.095
    # --------------------------------------------------------

    handles = [

        Patch(
            facecolor=CLASS_COLORS[
                cls
            ],

            edgecolor="none",

            label=cls
        )

        for cls in CLASSES
    ]


    outer.legend(
        handles=handles,
        labels=CLASSES,

        frameon=False,

        ncol=5,

        fontsize=9,

        loc="upper center",

        bbox_to_anchor=(
            0.50,
            1.095
        ),

        handlelength=1.15,
        handletextpad=0.45,

        columnspacing=1.00,

        borderaxespad=0
    )


    if panel_letter:

        outer.text(
            -0.075,
            1.045,

            "A",

            transform=outer.transAxes,

            fontsize=24,
            fontweight="bold",

            ha="left",
            va="top"
        )


    return outer


# ============================================================
# PANEL B
# ============================================================

def draw_signed_stack(
    ax,
    matrix,
    sample_label,
    show_x=False
):

    # --------------------------------------------------------
    # Highlight Kimura 0-2
    # --------------------------------------------------------

    ax.axvspan(
        -0.5,
        2.5,

        facecolor="#E6E6E6",
        alpha=0.55,

        linewidth=0,

        zorder=0
    )


    positive_bottom = np.zeros(
        len(
            BINS
        )
    )


    negative_bottom = np.zeros(
        len(
            BINS
        )
    )


    for cls in CLASSES:

        values = matrix[
            cls
        ].values


        positive = np.where(
            values > 0,
            values,
            0.0
        )


        negative = np.where(
            values < 0,
            values,
            0.0
        )


        ax.bar(
            BINS,
            positive,

            bottom=positive_bottom,

            width=0.92,

            color=CLASS_COLORS[
                cls
            ],

            edgecolor="none",
            linewidth=0,

            zorder=2
        )


        ax.bar(
            BINS,
            negative,

            bottom=negative_bottom,

            width=0.92,

            color=CLASS_COLORS[
                cls
            ],

            edgecolor="none",
            linewidth=0,

            zorder=2
        )


        positive_bottom += positive
        negative_bottom += negative