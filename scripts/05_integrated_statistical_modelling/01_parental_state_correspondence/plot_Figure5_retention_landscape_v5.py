#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib as mpl

from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap


# ============================================================
# CONFIG
# ============================================================

WGBS_DIR = os.environ.get("PANTE_WGBS_RETENTION_DIR", os.path.join(os.getcwd(), "data", "wgbs_retention"))
ATAC_DIR = os.environ.get("PANTE_ATAC_RETENTION_DIR", os.path.join(os.getcwd(), "data", "atac_retention"))
FIG = os.environ.get("PANTE_FIG5_SUMMARY_DIR", os.path.join(os.getcwd(), "results", "Figure5_summary"))
OUT = os.environ.get("PANTE_OUTPUT_DIR", os.path.join(os.getcwd(), "results", "Figure5_retention_landscape_v5"))
FONT_FILE = os.environ.get("PANTE_FONT_FILE", "").strip()

ABS_SUMMARY = os.environ.get("PANTE_FIG5_ABSOLUTE_SUMMARY", os.path.join(FIG, "Figure5_absolute_state_retention.tsv"))
CON_SUMMARY = os.environ.get("PANTE_FIG5_CONTRAST_SUMMARY", os.path.join(FIG, "Figure5_parental_contrast_retention.tsv"))
ATAC_MATRIX = os.environ.get("PANTE_ATAC_RETENTION_MATRIX", os.path.join(ATAC_DIR, "ATAC_parental_state_retention.four_state_matrix.tsv"))

WGBS_FILES = {
    "Gill": os.environ.get("PANTE_WGBS_GILL_RETENTION", os.path.join(WGBS_DIR, "gill.parental_mCG_difference_retention.tsv")),
    "Kidney": os.environ.get("PANTE_WGBS_KIDNEY_RETENTION", os.path.join(WGBS_DIR, "kidney.parental_mCG_difference_retention.tsv")),
    "Muscle": os.environ.get("PANTE_WGBS_MUSCLE_RETENTION", os.path.join(WGBS_DIR, "muscle.parental_mCG_difference_retention.tsv"))
}

N_BOOT = 1000
RANDOM_SEED = 20260908

os.makedirs(OUT, exist_ok=True)
required_files = [ABS_SUMMARY, CON_SUMMARY, ATAC_MATRIX] + list(WGBS_FILES.values())
if FONT_FILE:
    required_files.append(FONT_FILE)
for f in required_files:
    if not os.path.exists(f):
        raise RuntimeError("Missing input file: {}".format(f))

FP = FontProperties(fname=FONT_FILE) if FONT_FILE else FontProperties()


# ============================================================
# TYPOGRAPHY
# ============================================================

FP = FontProperties(fname=FONT_FILE) if FONT_FILE else FontProperties()


mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"

mpl.rcParams["axes.linewidth"] = 1.0

mpl.rcParams["xtick.major.width"] = 1.0
mpl.rcParams["ytick.major.width"] = 1.0

mpl.rcParams["xtick.direction"] = "out"
mpl.rcParams["ytick.direction"] = "out"


# ============================================================
# COLORS
# ============================================================

MCG_COLOR = "#70439A"
ATAC_COLOR = "#168C80"

MCG_REG = "#3F205C"
ATAC_REG = "#075C51"

GREY = "#747474"


MCG_CMAP = LinearSegmentedColormap.from_list(
    "mCG_density",
    [
        "#FFFFFF",
        "#F1ECF6",
        "#DACAE8",
        "#B696D0",
        "#7D55A5",
        "#402060"
    ]
)


ATAC_CMAP = LinearSegmentedColormap.from_list(
    "ATAC_density",
    [
        "#FFFFFF",
        "#E9F6F3",
        "#C7E8E1",
        "#84CBBE",
        "#349E91",
        "#075C51"
    ]
)


TISSUES = [
    "Gill",
    "Kidney",
    "Muscle"
]


TISSUE_MARKERS = {
    "Gill": "o",
    "Kidney": "s",
    "Muscle": "D"
}


# ============================================================
# BASIC HELPERS
# ============================================================

def clean_layer(x):

    x = str(x)

    if x == "DNA_methylation":
        return "DNA methylation"

    if x == "Chromatin_accessibility":
        return "Chromatin accessibility"

    return x


def clean_tissue(x):

    return {
        "gill": "Gill",
        "kidney": "Kidney",
        "muscle": "Muscle"
    }.get(
        str(x).lower(),
        str(x)
    )


def font_ticks(
    ax,
    size=9
):

    for lab in ax.get_xticklabels():

        lab.set_fontproperties(FP)
        lab.set_fontsize(size)


    for lab in ax.get_yticklabels():

        lab.set_fontproperties(FP)
        lab.set_fontsize(size)


def style_axis(
    ax
):

    ax.spines[
        "top"
    ].set_visible(False)

    ax.spines[
        "right"
    ].set_visible(False)


    ax.spines[
        "left"
    ].set_linewidth(1.0)

    ax.spines[
        "bottom"
    ].set_linewidth(1.0)


    ax.tick_params(
        width=1.0,
        length=4
    )


    font_ticks(
        ax,
        9
    )


def save_fig(
    fig,
    prefix
):

    fig.savefig(
        prefix + ".pdf",
        bbox_inches="tight",
        facecolor="white"
    )


    fig.savefig(
        prefix + ".svg",
        bbox_inches="tight",
        facecolor="white"
    )


    fig.savefig(
        prefix + ".png",
        dpi=600,
        bbox_inches="tight",
        facecolor="white"
    )


def spearman_rho(
    x,
    y
):

    d = pd.DataFrame({
        "x": x,
        "y": y
    }).dropna()


    rho = d[
        "x"
    ].corr(
        d[
            "y"
        ],
        method="spearman"
    )


    return (
        float(rho),
        len(d)
    )


def symmetric_limit(
    x,
    y,
    q=0.995
):

    v = np.concatenate(
        [
            np.asarray(x),
            np.asarray(y)
        ]
    )


    v = v[
        np.isfinite(v)
    ]


    lim = np.quantile(
        np.abs(v),
        q
    )


    return float(
        lim * 1.07
    )


# ============================================================
# FIND PanTE ID COLUMN
# ============================================================

def find_id_column(df):

    candidates = [
        "PanTE_ID",
        "PanTE",
        "PanTE_cluster",
        "PanTE_cluster_ID",
        "Cluster_ID",
        "cluster_id",
        "panTE_ID",
        "pante_id"
    ]


    lower_map = {
        str(c).lower(): c
        for c in df.columns
    }


    for candidate in candidates:

        key = candidate.lower()

        if key in lower_map:

            return lower_map[
                key
            ]


    return None


# ============================================================
# OLS
# ============================================================

def fit_weighted_ols(
    x,
    y,
    w=None
):

    x = np.asarray(
        x,
        dtype=float
    )

    y = np.asarray(
        y,
        dtype=float
    )


    if w is None:

        w = np.ones(
            len(x),
            dtype=float
        )

    else:

        w = np.asarray(
            w,
            dtype=float
        )


    sw = np.sum(
        w
    )


    if sw <= 0:

        return (
            np.nan,
            np.nan
        )


    xbar = np.sum(
        w * x
    ) / sw


    ybar = np.sum(
        w * y
    ) / sw


    denom = np.sum(
        w
        *
        (
            x - xbar
        ) ** 2
    )


    if denom <= 0:

        return (
            np.nan,
            np.nan
        )


    beta = (
        np.sum(
            w
            *
            (
                x - xbar
            )
            *
            (
                y - ybar
            )
        )
        /
        denom
    )


    intercept = (
        ybar
        -
        beta * xbar
    )


    return (
        float(intercept),
        float(beta)
    )


# ============================================================
# BOOTSTRAP OLS
# ============================================================

def bootstrap_regression(
    data,
    n_boot=N_BOOT,
    seed=RANDOM_SEED
):

    """
    Preferred:
      cluster bootstrap by PanTE_ID

    Fallback:
      row bootstrap

    Returns:
      intercept
      beta
      beta CI
      intercept CI
      bootstrap arrays
      bootstrap unit
    """


    d = data.copy()


    d = d.dropna(
        subset=[
            "Parent_delta",
            "BT_delta"
        ]
    )


    x = d[
        "Parent_delta"
    ].values.astype(float)


    y = d[
        "BT_delta"
    ].values.astype(float)


    intercept, beta = fit_weighted_ols(
        x,
        y
    )


    rng = np.random.RandomState(
        seed
    )


    boot_intercept = np.empty(
        n_boot,
        dtype=float
    )


    boot_beta = np.empty(
        n_boot,
        dtype=float
    )


    use_cluster = (
        "PanTE_ID" in d.columns
        and
        d[
            "PanTE_ID"
        ].notna().all()
        and
        d[
            "PanTE_ID"
        ].nunique() >= 100
    )


    if use_cluster:

        codes, unique_ids = pd.factorize(
            d[
                "PanTE_ID"
            ].astype(str)
        )


        n_cluster = len(
            unique_ids
        )


        probability = np.ones(
            n_cluster,
            dtype=float
        ) / float(
            n_cluster
        )


        for b in range(
            n_boot
        ):

            cluster_count = rng.multinomial(
                n_cluster,
                probability
            )


            row_weights = cluster_count[
                codes
            ].astype(float)


            ai, bi = fit_weighted_ols(
                x,
                y,
                row_weights
            )


            boot_intercept[
                b
            ] = ai


            boot_beta[
                b
            ] = bi


        bootstrap_unit = (
            "PanTE_ID cluster"
        )


    else:

        n = len(
            x
        )


        for b in range(
            n_boot
        ):

            idx = rng.randint(
                0,
                n,
                size=n
            )


            ai, bi = fit_weighted_ols(
                x[
                    idx
                ],
                y[
                    idx
                ]
            )


            boot_intercept[
                b
            ] = ai


            boot_beta[
                b
            ] = bi


        bootstrap_unit = (
            "row"
        )


    valid = (
        np.isfinite(
            boot_intercept
        )
        &
        np.isfinite(
            boot_beta
        )
    )


    boot_intercept = boot_intercept[
        valid
    ]


    boot_beta = boot_beta[
        valid
    ]


    beta_low = np.quantile(
        boot_beta,
        0.025
    )


    beta_high = np.quantile(
        boot_beta,
        0.975
    )


    int_low = np.quantile(
        boot_intercept,
        0.025
    )


    int_high = np.quantile(
        boot_intercept,
        0.975
    )


    return {
        "intercept":
            intercept,

        "beta":
            beta,

        "beta_low":
            float(beta_low),

        "beta_high":
            float(beta_high),

        "intercept_low":
            float(int_low),

        "intercept_high":
            float(int_high),

        "boot_intercept":
            boot_intercept,

        "boot_beta":
            boot_beta,

        "bootstrap_unit":
            bootstrap_unit,

        "n_boot_valid":
            len(boot_beta)
    }


# ============================================================
# PREDICTION CI
# ============================================================

def bootstrap_prediction_band(
    regression,
    xgrid
):

    intercepts = regression[
        "boot_intercept"
    ]


    betas = regression[
        "boot_beta"
    ]


    pred = (
        intercepts[
            :,
            None
        ]
        +
        betas[
            :,
            None
        ]
        *
        xgrid[
            None,
            :
        ]
    )


    low = np.quantile(
        pred,
        0.025,
        axis=0
    )


    high = np.quantile(
        pred,
        0.975,
        axis=0
    )


    return (
        low,
        high
    )


# ============================================================
# SUMMARY DATA
# ============================================================

ABS = pd.read_csv(
    ABS_SUMMARY,
    sep="\t"
)


CON = pd.read_csv(
    CON_SUMMARY,
    sep="\t"
)


ABS[
    "Layer"
] = ABS[
    "Layer"
].map(
    clean_layer
)


ABS[
    "Tissue"
] = ABS[
    "Tissue"
].map(
    clean_tissue
)


CON[
    "Layer"
] = CON[
    "Layer"
].map(
    clean_layer
)


CON[
    "Tissue"
] = CON[
    "Tissue"
].map(
    clean_tissue
)


# ============================================================
# PANEL A DATA
# ============================================================

phase_rows = []


for layer in [
    "DNA methylation",
    "Chromatin accessibility"
]:

    for tissue in TISSUES:

        a = ABS[
            (ABS["Layer"] == layer)
            &
            (ABS["Tissue"] == tissue)
        ]


        c = CON[
            (CON["Layer"] == layer)
            &
            (CON["Tissue"] == tissue)
        ]


        if len(c) != 1:

            raise RuntimeError(
                "Expected one contrast row: {} {}".format(
                    layer,
                    tissue
                )
            )


        y = float(
            c.iloc[0][
                "Spearman_rho"
            ]
        )


        for comparison in [
            "BSB_to_BTM",
            "TC_to_BTP"
        ]:

            q = a[
                a[
                    "Comparison"
                ] == comparison
            ]


            if len(q) != 1:

                raise RuntimeError(
                    "Expected one state row: {} {} {}".format(
                        layer,
                        tissue,
                        comparison
                    )
                )


            phase_rows.append({
                "Layer":
                    layer,

                "Tissue":
                    tissue,

                "Comparison":
                    comparison,

                "State_rho":
                    float(
                        q.iloc[0][
                            "Spearman_rho"
                        ]
                    ),

                "Contrast_rho":
                    y
            })


PHASE = pd.DataFrame(
    phase_rows
)


PHASE.to_csv(
    os.path.join(
        OUT,
        "Figure5_v5_phase_map_data.tsv"
    ),
    sep="\t",
    index=False
)


# ============================================================
# LOAD WGBS
# ============================================================

wgbs_rows = []


for tissue in TISSUES:

    raw = pd.read_csv(
        WGBS_FILES[
            tissue
        ],
        sep="\t"
    )


    id_col = find_id_column(
        raw
    )


    d = pd.DataFrame({
        "Tissue":
            tissue,

        "Parent_delta":
            pd.to_numeric(
                raw[
                    "Parent_delta_mCG"
                ],
                errors="coerce"
            ),

        "BT_delta":
            pd.to_numeric(
                raw[
                    "BT_haplotype_delta_mCG"
                ],
                errors="coerce"
            )
    })


    if id_col is not None:

        d[
            "PanTE_ID"
        ] = raw[
            id_col
        ].astype(str).values


    d = d.dropna(
        subset=[
            "Parent_delta",
            "BT_delta"
        ]
    )


    wgbs_rows.append(
        d
    )


WGBS_DAT = pd.concat(
    wgbs_rows,
    ignore_index=True
)


# ============================================================
# LOAD ATAC
# ============================================================

ATAC_RAW = pd.read_csv(
    ATAC_MATRIX,
    sep="\t"
)


atac_id_col = find_id_column(
    ATAC_RAW
)


ATAC_DAT = pd.DataFrame({
    "Tissue":
        ATAC_RAW[
            "Tissue"
        ].map(
            clean_tissue
        ),

    "Parent_delta":
        pd.to_numeric(
            ATAC_RAW[
                "Parent_delta"
            ],
            errors="coerce"
        ),

    "BT_delta":
        pd.to_numeric(
            ATAC_RAW[
                "BT_delta"
            ],
            errors="coerce"
        )
})


if atac_id_col is not None:

    ATAC_DAT[
        "PanTE_ID"
    ] = ATAC_RAW[
        atac_id_col
    ].astype(str).values


ATAC_DAT = ATAC_DAT.dropna(
    subset=[
        "Parent_delta",
        "BT_delta"
    ]
)


# ============================================================
# RUN REGRESSIONS
# ============================================================

print("")
print("Running DNA methylation bootstrap regression...")


REG_WGBS = bootstrap_regression(
    WGBS_DAT
)


print(
    "  beta = {:.6f} [{:.6f}, {:.6f}]".format(
        REG_WGBS[
            "beta"
        ],
        REG_WGBS[
            "beta_low"
        ],
        REG_WGBS[
            "beta_high"
        ]
    )
)


print(
    "  bootstrap unit = {}".format(
        REG_WGBS[
            "bootstrap_unit"
        ]
    )
)


print("")
print("Running ATAC bootstrap regression...")


REG_ATAC = bootstrap_regression(
    ATAC_DAT,
    seed=RANDOM_SEED + 1
)


print(
    "  beta = {:.6f} [{:.6f}, {:.6f}]".format(
        REG_ATAC[
            "beta"
        ],
        REG_ATAC[
            "beta_low"
        ],
        REG_ATAC[
            "beta_high"
        ]
    )
)


print(
    "  bootstrap unit = {}".format(
        REG_ATAC[
            "bootstrap_unit"
        ]
    )
)


# ============================================================
# PANEL A
# ============================================================

def draw_phase_map(
    ax
):

    style_axis(
        ax
    )


    # subtle high-state-retention region
    ax.axvspan(
        0.83,
        0.872,
        color="#F5F7F8",
        alpha=0.80,
        zorder=0
    )


    for _, r in PHASE.iterrows():

        layer = r[
            "Layer"
        ]


        tissue = r[
            "Tissue"
        ]


        comparison = r[
            "Comparison"
        ]


        color = (
            MCG_COLOR
            if layer == "DNA methylation"
            else ATAC_COLOR
        )


        x = r[
            "State_rho"
        ]


        # EXACT same y for both matched-parent comparisons
        y = r[
            "Contrast_rho"
        ]


        # deliberately longer local segment
        seg = 0.0052


        if comparison == "BSB_to_BTM":

            # solid line + filled marker
            ax.plot(
                [
                    x - seg,
                    x + seg
                ],
                [
                    y,
                    y
                ],
                linestyle="-",
                color=color,
                linewidth=2.35,
                solid_capstyle="round",
                zorder=2
            )


            ax.scatter(
                x,
                y,
                s=105,
                marker=TISSUE_MARKERS[
                    tissue
                ],
                facecolor=color,
                edgecolor="white",
                linewidth=1.25,
                zorder=4
            )


        else:

            # strong custom dash + OPEN marker
            ax.plot(
                [
                    x - seg,
                    x + seg
                ],
                [
                    y,
                    y
                ],
                linestyle=(
                    0,
                    (
                        5,
                        2
                    )
                ),
                color=color,
                linewidth=2.25,
                dash_capstyle="butt",
                zorder=2
            )


            ax.scatter(
                x,
                y,
                s=105,
                marker=TISSUE_MARKERS[
                    tissue
                ],
                facecolor="white",
                edgecolor=color,
                linewidth=1.9,
                zorder=5
            )


    ax.set_xlim(
        0.790,
        0.872
    )


    ax.set_ylim(
        0.280,
        0.452
    )


    ax.set_xlabel(
        "Matched parental-state retention\n(Spearman $\\rho$)",
        fontproperties=FP,
        fontsize=10.5
    )


    ax.set_ylabel(
        "Parental-contrast retention\n(Spearman $\\rho$)",
        fontproperties=FP,
        fontsize=10.5
    )


    ax.set_title(
        "PanTE epigenetic retention landscape",
        loc="left",
        fontproperties=FP,
        fontsize=13.5,
        pad=10
    )


    # ========================================================
    # LEGENDS
    # ========================================================

    layer_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=MCG_COLOR,
            markeredgecolor="white",
            markersize=8,
            label="DNA methylation"
        ),

        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=ATAC_COLOR,
            markeredgecolor="white",
            markersize=8,
            label="Chromatin accessibility"
        )
    ]


    tissue_handles = [
        Line2D(
            [0],
            [0],
            marker=TISSUE_MARKERS[
                t
            ],
            linestyle="none",
            markerfacecolor="#606060",
            markeredgecolor="white",
            markersize=7,
            label=t
        )
        for t in TISSUES
    ]


    comparison_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="-",
            color="#4D4D4D",
            markerfacecolor="#4D4D4D",
            markeredgecolor="white",
            linewidth=2.1,
            markersize=6.5,
            label="BSB -> BT-M"
        ),

        Line2D(
            [0],
            [0],
            marker="o",
            linestyle=(
                0,
                (
                    5,
                    2
                )
            ),
            color="#4D4D4D",
            markerfacecolor="white",
            markeredgecolor="#4D4D4D",
            markeredgewidth=1.5,
            linewidth=2.1,
            markersize=6.5,
            label="TC -> BT-P"
        )
    ]


    leg1 = ax.legend(
        handles=layer_handles,
        loc="upper left",
        frameon=False,
        fontsize=8.0
    )


    for t in leg1.get_texts():
        t.set_fontproperties(FP)


    ax.add_artist(
        leg1
    )


    leg2 = ax.legend(
        handles=tissue_handles,
        loc="lower right",
        frameon=False,
        fontsize=7.5
    )


    for t in leg2.get_texts():
        t.set_fontproperties(FP)


    ax.add_artist(
        leg2
    )


    leg3 = ax.legend(
        handles=comparison_handles,
        loc="center left",
        bbox_to_anchor=(
            0.01,
            0.55
        ),
        frameon=False,
        fontsize=7.4,
        handlelength=2.6
    )


    for t in leg3.get_texts():
        t.set_fontproperties(FP)


# ============================================================
# B/C JOINT REGRESSION LANDSCAPE
# ============================================================

def draw_joint_panel(
    fig,
    subspec,
    data,
    regression,
    layer_name,
    layer_color,
    regression_color,
    cmap,
    x_label,
    y_label
):

    gs = GridSpecFromSubplotSpec(
        2,
        2,
        subplot_spec=subspec,
        width_ratios=[
            5.0,
            1.05
        ],
        height_ratios=[
            1.00,
            4.7
        ],
        wspace=0.025,
        hspace=0.025
    )


    ax_top = fig.add_subplot(
        gs[
            0,
            0
        ]
    )


    ax_main = fig.add_subplot(
        gs[
            1,
            0
        ]
    )


    ax_right = fig.add_subplot(
        gs[
            1,
            1
        ]
    )


    x = data[
        "Parent_delta"
    ].values


    y = data[
        "BT_delta"
    ].values


    lim = symmetric_limit(
        x,
        y
    )


    # ========================================================
    # DENSITY
    # ========================================================

    ax_main.hexbin(
        x,
        y,
        gridsize=70,
        extent=[
            -lim,
            lim,
            -lim,
            lim
        ],
        bins="log",
        mincnt=1,
        cmap=cmap,
        linewidths=0,
        rasterized=True
    )


    ax_main.axhline(
        0,
        color="#B7B7B7",
        linewidth=0.8
    )


    ax_main.axvline(
        0,
        color="#B7B7B7",
        linewidth=0.8
    )


    # y = x reference
    ax_main.plot(
        [
            -lim,
            lim
        ],
        [
            -lim,
            lim
        ],
        linestyle="--",
        linewidth=1.05,
        color="#555555",
        alpha=0.78,
        zorder=4
    )


    # ========================================================
    # LINEAR REGRESSION
    # ========================================================

    # only central x range for regression display
    x_low = np.quantile(
        x,
        0.01
    )


    x_high = np.quantile(
        x,
        0.99
    )


    xgrid = np.linspace(
        x_low,
        x_high,
        220
    )


    yhat = (
        regression[
            "intercept"
        ]
        +
        regression[
            "beta"
        ]
        *
        xgrid
    )


    ci_low, ci_high = bootstrap_prediction_band(
        regression,
        xgrid
    )


    # 95% bootstrap CI
    ax_main.fill_between(
        xgrid,
        ci_low,
        ci_high,
        color=layer_color,
        alpha=0.20,
        linewidth=0,
        zorder=5
    )


    # fitted line
    ax_main.plot(
        xgrid,
        yhat,
        color=regression_color,
        linewidth=2.8,
        zorder=7
    )


    # ========================================================
    # AXIS
    # ========================================================

    ax_main.set_xlim(
        -lim,
        lim
    )


    ax_main.set_ylim(
        -lim,
        lim
    )


    ax_main.set_aspect(
        "equal",
        adjustable="box"
    )


    ax_main.set_xlabel(
        x_label,
        fontproperties=FP,
        fontsize=9.7
    )


    ax_main.set_ylabel(
        y_label,
        fontproperties=FP,
        fontsize=9.7
    )


    style_axis(
        ax_main
    )


    # ========================================================
    # STATS
    # ========================================================

    rho, n = spearman_rho(
        x,
        y
    )


    med_parent = np.median(
        np.abs(x)
    )


    med_bt = np.median(
        np.abs(y)
    )


    ratio = (
        med_bt
        /
        med_parent
    )


    annotation = (
        "$\\rho$={:.2f}\n"
        "$\\beta$={:.2f} [{:.2f}, {:.2f}]"
    ).format(
        rho,
        regression[
            "beta"
        ],
        regression[
            "beta_low"
        ],
        regression[
            "beta_high"
        ]
    )


    ax_main.text(
        0.035,
        0.965,
        annotation,
        transform=ax_main.transAxes,
        ha="left",
        va="top",
        fontproperties=FP,
        fontsize=8.1,
        color="#303030",
        bbox=dict(
            boxstyle="round,pad=0.32",
            facecolor="white",
            edgecolor="#D5D5D5",
            alpha=0.95
        )
    )


    # ========================================================
    # TOP MARGINAL
    # ========================================================

    bins = np.linspace(
        -lim,
        lim,
        72
    )


    ax_top.hist(
        x,
        bins=bins,
        density=True,
        color=layer_color,
        alpha=0.66,
        edgecolor="none"
    )


    ax_top.axvline(
        0,
        color="#B0B0B0",
        linewidth=0.7
    )


    ax_top.set_xlim(
        -lim,
        lim
    )


    ax_top.set_xticks([])
    ax_top.set_yticks([])


    for s in [
        "left",
        "right",
        "top",
        "bottom"
    ]:

        ax_top.spines[
            s
        ].set_visible(False)


    ax_top.text(
        0.0,
        1.28,
        layer_name,
        transform=ax_top.transAxes,
        ha="left",
        va="top",
        fontproperties=FP,
        fontsize=13.5,
        color=layer_color
    )


    # ========================================================
    # RIGHT MARGINAL
    # ========================================================

    ax_right.hist(
        y,
        bins=bins,
        density=True,
        orientation="horizontal",
        color=layer_color,
        alpha=0.66,
        edgecolor="none"
    )


    ax_right.axhline(
        0,
        color="#B0B0B0",
        linewidth=0.7
    )


    ax_right.set_ylim(
        -lim,
        lim
    )


    ax_right.set_xticks([])
    ax_right.set_yticks([])


    for s in [
        "left",
        "right",
        "top",
        "bottom"
    ]:

        ax_right.spines[
            s
        ].set_visible(False)


    return {
        "rho":
            rho,

        "N":
            n,

        "R_absDelta":
            ratio,

        "intercept":
            regression[
                "intercept"
            ],

        "intercept_low":
            regression[
                "intercept_low"
            ],

        "intercept_high":
            regression[
                "intercept_high"
            ],

        "beta":
            regression[
                "beta"
            ],

        "beta_low":
            regression[
                "beta_low"
            ],

        "beta_high":
            regression[
                "beta_high"
            ],

        "bootstrap_unit":
            regression[
                "bootstrap_unit"
            ]
    }


# ============================================================
# INDIVIDUAL A
# ============================================================

figA, axA = plt.subplots(
    figsize=(
        5.2,
        5.0
    ),
    facecolor="white"
)


draw_phase_map(
    axA
)


figA.subplots_adjust(
    left=0.19,
    right=0.97,
    bottom=0.16,
    top=0.90
)


save_fig(
    figA,
    os.path.join(
        OUT,
        "Figure5_v5_panel_A_retention_map"
    )
)


plt.close(
    figA
)


# ============================================================
# INDIVIDUAL B
# ============================================================

figB = plt.figure(
    figsize=(
        6.7,
        5.6
    ),
    facecolor="white"
)


gB = GridSpec(
    1,
    1,
    figure=figB
)


statB_single = draw_joint_panel(
    figB,
    gB[
        0
    ],
    WGBS_DAT,
    REG_WGBS,
    "DNA methylation",
    MCG_COLOR,
    MCG_REG,
    MCG_CMAP,
    "Parent $\\Delta$mCG",
    "BT haplotype $\\Delta$mCG"
)


figB.subplots_adjust(
    left=0.12,
    right=0.96,
    bottom=0.13,
    top=0.88
)


save_fig(
    figB,
    os.path.join(
        OUT,
        "Figure5_v5_panel_B_mCG_regression"
    )
)


plt.close(
    figB
)


# ============================================================
# INDIVIDUAL C
# ============================================================

figC = plt.figure(
    figsize=(
        6.7,
        5.6
    ),
    facecolor="white"
)


gC = GridSpec(
    1,
    1,
    figure=figC
)


statC_single = draw_joint_panel(
    figC,
    gC[
        0
    ],
    ATAC_DAT,
    REG_ATAC,
    "Chromatin accessibility",
    ATAC_COLOR,
    ATAC_REG,
    ATAC_CMAP,
    "Parent $\\Delta$accessibility",
    "BT haplotype $\\Delta$accessibility"
)


figC.subplots_adjust(
    left=0.12,
    right=0.96,
    bottom=0.13,
    top=0.88
)


save_fig(
    figC,
    os.path.join(
        OUT,
        "Figure5_v5_panel_C_ATAC_regression"
    )
)


plt.close(
    figC
)


# ============================================================
# COMBINED
#
#     A | B
#     A | C
# ============================================================

fig = plt.figure(
    figsize=(
        13.8,
        9.2
    ),
    facecolor="white"
)


outer = GridSpec(
    2,
    2,
    figure=fig,
    width_ratios=[
        0.37,
        0.63
    ],
    height_ratios=[
        1.0,
        1.0
    ],
    wspace=0.27,
    hspace=0.24
)


axA = fig.add_subplot(
    outer[
        :,
        0
    ]
)


draw_phase_map(
    axA
)


statB = draw_joint_panel(
    fig,
    outer[
        0,
        1
    ],
    WGBS_DAT,
    REG_WGBS,
    "DNA methylation",
    MCG_COLOR,
    MCG_REG,
    MCG_CMAP,
    "Parent $\\Delta$mCG",
    "BT haplotype $\\Delta$mCG"
)


statC = draw_joint_panel(
    fig,
    outer[
        1,
        1
    ],
    ATAC_DAT,
    REG_ATAC,
    "Chromatin accessibility",
    ATAC_COLOR,
    ATAC_REG,
    ATAC_CMAP,
    "Parent $\\Delta$accessibility",
    "BT haplotype $\\Delta$accessibility"
)


# ============================================================
# ALIGN PANEL LETTERS
# ============================================================

fig.subplots_adjust(
    left=0.070,
    right=0.985,
    bottom=0.072,
    top=0.945
)


# force geometry realization
fig.canvas.draw()


A_pos = axA.get_position()


# obtain top histogram / main axes positions
# axes created after A:
all_axes = fig.axes


# Find B/C left boundaries by inspecting axes positions.
# B and C are on the right column and share same x0.
right_axes = [
    a
    for a in all_axes
    if a.get_position().x0 > 0.40
]


right_x0 = min(
    [
        a.get_position().x0
        for a in right_axes
    ]
)


# B top = highest right-side axis
right_top = max(
    [
        a.get_position().y1
        for a in right_axes
    ]
)


# A and B panel letters share exact top baseline
letter_top = max(
    A_pos.y1,
    right_top
) + 0.012


# C begins around lower right group
lower_right_axes = [
    a
    for a in right_axes
    if a.get_position().y0 < 0.50
]


C_top = max(
    [
        a.get_position().y1
        for a in lower_right_axes
    ]
) + 0.012


fig.text(
    A_pos.x0 - 0.032,
    letter_top,
    "A",
    ha="left",
    va="top",
    fontproperties=FP,
    fontsize=25
)


fig.text(
    right_x0 - 0.032,
    letter_top,
    "B",
    ha="left",
    va="top",
    fontproperties=FP,
    fontsize=25
)


fig.text(
    right_x0 - 0.032,
    C_top,
    "C",
    ha="left",
    va="top",
    fontproperties=FP,
    fontsize=25
)


PREFIX = os.path.join(
    OUT,
    "Figure5_v5_journal_layout"
)


save_fig(
    fig,
    PREFIX
)


plt.close(
    fig
)


# ============================================================
# STATISTICS TABLE
# ============================================================

stats = pd.DataFrame([
    {
        "Layer":
            "DNA_methylation",

        "Spearman_rho":
            statB[
                "rho"
            ],

        "N":
            statB[
                "N"
            ],

        "R_absDelta":
            statB[
                "R_absDelta"
            ],

        "Intercept":
            statB[
                "intercept"
            ],

        "Intercept_CI_low":
            statB[
                "intercept_low"
            ],

        "Intercept_CI_high":
            statB[
                "intercept_high"
            ],

        "Beta":
            statB[
                "beta"
            ],

        "Beta_CI_low":
            statB[
                "beta_low"
            ],

        "Beta_CI_high":
            statB[
                "beta_high"
            ],

        "Bootstrap_unit":
            statB[
                "bootstrap_unit"
            ],

        "N_boot":
            N_BOOT
    },

    {
        "Layer":
            "Chromatin_accessibility",

        "Spearman_rho":
            statC[
                "rho"
            ],

        "N":
            statC[
                "N"
            ],

        "R_absDelta":
            statC[
                "R_absDelta"
            ],

        "Intercept":
            statC[
                "intercept"
            ],

        "Intercept_CI_low":
            statC[
                "intercept_low"
            ],

        "Intercept_CI_high":
            statC[
                "intercept_high"
            ],

        "Beta":
            statC[
                "beta"
            ],

        "Beta_CI_low":
            statC[
                "beta_low"
            ],

        "Beta_CI_high":
            statC[
                "beta_high"
            ],

        "Bootstrap_unit":
            statC[
                "bootstrap_unit"
            ],

        "N_boot":
            N_BOOT
    }
])


STATS_OUT = os.path.join(
    OUT,
    "Figure5_v5_regression_statistics.tsv"
)


stats.to_csv(
    STATS_OUT,
    sep="\t",
    index=False
)


# ============================================================
# LOG
# ============================================================

LOG = os.path.join(
    OUT,
    "Figure5_v5_journal_layout.log"
)


with open(
    LOG,
    "w"
) as out:

    out.write(
        "Figure 5 v5 completed.\n\n"
    )


    out.write(
        "IMPORTANT:\n"
    )


    out.write(
        "Regression is descriptive pooled OLS.\n"
    )


    out.write(
        "95% CI is based on {} bootstrap replicates.\n".format(
            N_BOOT
        )
    )


    out.write(
        "Bootstrap unit is PanTE_ID when available; "
        "otherwise row-level resampling is used.\n\n"
    )


    for _, r in stats.iterrows():

        out.write(
            "{}\n".format(
                r[
                    "Layer"
                ]
            )
        )


        out.write(
            "Spearman rho = {:.6f}\n".format(
                r[
                    "Spearman_rho"
                ]
            )
        )


        out.write(
            "N = {}\n".format(
                int(
                    r[
                        "N"
                    ]
                )
            )
        )


        out.write(
            "R_absDelta = {:.6f}\n".format(
                r[
                    "R_absDelta"
                ]
            )
        )


        out.write(
            "Intercept = {:.6f} "
            "[{:.6f}, {:.6f}]\n".format(
                r[
                    "Intercept"
                ],
                r[
                    "Intercept_CI_low"
                ],
                r[
                    "Intercept_CI_high"
                ]
            )
        )


        out.write(
            "Beta = {:.6f} "
            "[{:.6f}, {:.6f}]\n".format(
                r[
                    "Beta"
                ],
                r[
                    "Beta_CI_low"
                ],
                r[
                    "Beta_CI_high"
                ]
            )
        )


        out.write(
            "Bootstrap unit = {}\n\n".format(
                r[
                    "Bootstrap_unit"
                ]
            )
        )


print("")
print("============================================")
print("Figure 5 v5 completed successfully.")
print("============================================")

print("")
print("Combined:")
print(
    PREFIX + ".png"
)

print("")
print("Individual panels:")

print(
    os.path.join(
        OUT,
        "Figure5_v5_panel_A_retention_map.png"
    )
)

print(
    os.path.join(
        OUT,
        "Figure5_v5_panel_B_mCG_regression.png"
    )
)

print(
    os.path.join(
        OUT,
        "Figure5_v5_panel_C_ATAC_regression.png"
    )
)

print("")
print("Regression statistics:")
print(
    STATS_OUT
)

print("")
print("Log:")
print(
    LOG
)