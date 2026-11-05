#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Prepare Figure 5 matched-parent state and parental-contrast summaries.

This release version calculates WGBS statistics from the six main-filter
matrices rather than embedding fixed summary values.
"""

from __future__ import print_function

import os
import pandas as pd
from scipy.stats import spearmanr

WGBS_DIR = os.environ.get("PANTE_WGBS_MCG_MATRIX_DIR")
ATAC_MATRIX = os.environ.get("PANTE_ATAC_RETENTION_MATRIX")
OUTDIR = os.environ.get("PANTE_OUTPUT_DIR", os.path.join(os.getcwd(), "results", "Figure5_summary"))

if not WGBS_DIR:
    raise RuntimeError("Set PANTE_WGBS_MCG_MATRIX_DIR to the directory containing the six main-filter WGBS matrices.")
if not ATAC_MATRIX:
    raise RuntimeError("Set PANTE_ATAC_RETENTION_MATRIX to ATAC_parental_state_retention.four_state_matrix.tsv.")
os.makedirs(OUTDIR, exist_ok=True)

TISSUES = ["Gill", "Kidney", "Muscle"]


def read_tsv(path):
    if not os.path.isfile(path):
        raise RuntimeError("Missing input file: {}".format(path))
    return pd.read_csv(path, sep="\t")


def rho_and_n(x, y):
    d = pd.DataFrame({
        "x": pd.to_numeric(x, errors="coerce"),
        "y": pd.to_numeric(y, errors="coerce")
    }).dropna()
    if len(d) < 2:
        return len(d), float("nan")
    return len(d), float(spearmanr(d["x"].values, d["y"].values)[0])


absolute_rows = []
contrast_rows = []

for tissue in TISSUES:
    stem = tissue.lower()
    bsb_file = os.path.join(WGBS_DIR, "main.{}.BSB_to_BTM.PanTE_mCG.tsv".format(stem))
    tc_file = os.path.join(WGBS_DIR, "main.{}.TC_to_BTP.PanTE_mCG.tsv".format(stem))
    bsb = read_tsv(bsb_file)
    tc = read_tsv(tc_file)

    for col in ["PanTE_ID", "BSB_pooled_mCG", "BTM_pooled_mCG"]:
        if col not in bsb.columns:
            raise RuntimeError("{} missing column {}".format(bsb_file, col))
    for col in ["PanTE_ID", "TC_pooled_mCG", "BTP_pooled_mCG"]:
        if col not in tc.columns:
            raise RuntimeError("{} missing column {}".format(tc_file, col))
    if bsb["PanTE_ID"].duplicated().any() or tc["PanTE_ID"].duplicated().any():
        raise RuntimeError("Duplicated PanTE_ID found in WGBS input for {}".format(tissue))

    n, rho = rho_and_n(bsb["BSB_pooled_mCG"], bsb["BTM_pooled_mCG"])
    absolute_rows.append({"Layer": "DNA_methylation", "Tissue": tissue,
                          "Comparison": "BSB_to_BTM", "N_PanTE": n, "Spearman_rho": rho})
    n, rho = rho_and_n(tc["TC_pooled_mCG"], tc["BTP_pooled_mCG"])
    absolute_rows.append({"Layer": "DNA_methylation", "Tissue": tissue,
                          "Comparison": "TC_to_BTP", "N_PanTE": n, "Spearman_rho": rho})

    joined = bsb[["PanTE_ID", "BSB_pooled_mCG", "BTM_pooled_mCG"]].merge(
        tc[["PanTE_ID", "TC_pooled_mCG", "BTP_pooled_mCG"]],
        on="PanTE_ID", how="inner", validate="one_to_one"
    )
    parent_delta = joined["BSB_pooled_mCG"] - joined["TC_pooled_mCG"]
    bt_delta = joined["BTM_pooled_mCG"] - joined["BTP_pooled_mCG"]
    n, rho = rho_and_n(parent_delta, bt_delta)
    contrast_rows.append({"Layer": "DNA_methylation", "Tissue": tissue,
                          "N_PanTE": n, "Spearman_rho": rho})

atac = read_tsv(ATAC_MATRIX)
required_atac = [
    "Tissue", "PanTE_ID", "BSB_log2FPKM1", "BTM_log2FPKM1",
    "TC_log2FPKM1", "BTP_log2FPKM1", "Parent_delta", "BT_delta"
]
for col in required_atac:
    if col not in atac.columns:
        raise RuntimeError("{} missing column {}".format(ATAC_MATRIX, col))
atac["Tissue_clean"] = atac["Tissue"].astype(str).str.strip().str.lower().map(
    {"gill": "Gill", "kidney": "Kidney", "muscle": "Muscle"}
)

for tissue in TISSUES:
    x = atac[atac["Tissue_clean"] == tissue].copy()
    if x["PanTE_ID"].duplicated().any():
        raise RuntimeError("Duplicated PanTE_ID in ATAC matrix for {}".format(tissue))

    n, rho = rho_and_n(x["BSB_log2FPKM1"], x["BTM_log2FPKM1"])
    absolute_rows.append({"Layer": "Chromatin_accessibility", "Tissue": tissue,
                          "Comparison": "BSB_to_BTM", "N_PanTE": n, "Spearman_rho": rho})
    n, rho = rho_and_n(x["TC_log2FPKM1"], x["BTP_log2FPKM1"])
    absolute_rows.append({"Layer": "Chromatin_accessibility", "Tissue": tissue,
                          "Comparison": "TC_to_BTP", "N_PanTE": n, "Spearman_rho": rho})
    n, rho = rho_and_n(x["Parent_delta"], x["BT_delta"])
    contrast_rows.append({"Layer": "Chromatin_accessibility", "Tissue": tissue,
                          "N_PanTE": n, "Spearman_rho": rho})

absolute = pd.DataFrame(absolute_rows)
contrast = pd.DataFrame(contrast_rows)

summary_rows = []
for layer in ["DNA_methylation", "Chromatin_accessibility"]:
    a = absolute.loc[absolute["Layer"] == layer, "Spearman_rho"]
    c = contrast.loc[contrast["Layer"] == layer, "Spearman_rho"]
    summary_rows.append({
        "Layer": layer,
        "Mean_absolute_rho": a.mean(),
        "Min_absolute_rho": a.min(),
        "Max_absolute_rho": a.max(),
        "Mean_contrast_rho": c.mean(),
        "Min_contrast_rho": c.min(),
        "Max_contrast_rho": c.max()
    })
summary = pd.DataFrame(summary_rows)

absolute.to_csv(os.path.join(OUTDIR, "Figure5_absolute_state_retention.tsv"), sep="\t", index=False)
contrast.to_csv(os.path.join(OUTDIR, "Figure5_parental_contrast_retention.tsv"), sep="\t", index=False)
summary.to_csv(os.path.join(OUTDIR, "Figure5_layer_summary.tsv"), sep="\t", index=False)

print("Figure 5 summary tables written to: {}".format(OUTDIR))
print("\nAbsolute-state retention:\n{}".format(absolute.to_string(index=False)))
print("\nParental-contrast retention:\n{}".format(contrast.to_string(index=False)))
