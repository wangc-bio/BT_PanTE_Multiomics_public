#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
import sys
import numpy as np
import pandas as pd

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


def sample_names_from_bam_list(path: str):
    """Extract sample prefix from bam list: basename without .bam"""
    names = []
    with open(path, "r") as f:
        for line in f:
            bam = line.strip()
            if not bam:
                continue
            base = os.path.basename(bam)
            if base.endswith(".bam"):
                base = base[:-4]
            names.append(base)
    return names


def kmer_call_to_simple(call: str) -> str:
    if call == "BSB_source":
        return "BSB"
    if call == "TC_source":
        return "TC"
    return "Mixed"


def main():
    ap = argparse.ArgumentParser(
        description="Validate k-mer genome-source calls using normalized window depth (ND matrix)."
    )
    ap.add_argument("--nd_matrix", required=True, help="ND_matrix.winXXXX.tsv")
    ap.add_argument("--bsb_list", required=True, help="BSB.bams.list (absolute paths, one per line)")
    ap.add_argument("--tc_list", required=True, help="TC.bams.list (absolute paths, one per line)")
    ap.add_argument("--kmer_source", required=True, help="kmer window source TSV (5Mb or 1Mb)")
    ap.add_argument("--out_prefix", required=True, help="output prefix")
    ap.add_argument("--fdr", type=float, default=0.05, help="FDR cutoff (default 0.05)")
    ap.add_argument("--log2fc", type=float, default=0.5, help="abs(log2FC) cutoff (default 0.5)")
    ap.add_argument("--only_chr", action="store_true", help="keep only chrXXA/chrXXB")
    args = ap.parse_args()

    # -----------------------
    # Load ND matrix
    # -----------------------
    df = pd.read_csv(args.nd_matrix, sep="\t", dtype={"chr": str})
    if not {"chr", "start", "end"}.issubset(df.columns[:3]):
        # allow alternate header
        pass

    # Optional: only chr
    if args.only_chr:
        df = df[df["chr"].str.match(r"^chr\d{2}[AB]$")].copy()

    # Determine sample columns by bam lists
    bsb_samples = sample_names_from_bam_list(args.bsb_list)
    tc_samples = sample_names_from_bam_list(args.tc_list)

    # ND matrix columns are sample names; keep only those present
    nd_cols = [c for c in df.columns if c not in ["chr", "start", "end"]]
    bsb_cols = [s for s in bsb_samples if s in nd_cols]
    tc_cols = [s for s in tc_samples if s in nd_cols]

    if len(bsb_cols) == 0 or len(tc_cols) == 0:
        print("[ERROR] Could not find BSB/TC sample columns in ND matrix header.", file=sys.stderr)
        print(f"  BSB columns found: {len(bsb_cols)} / expected {len(bsb_samples)}", file=sys.stderr)
        print(f"  TC  columns found: {len(tc_cols)} / expected {len(tc_samples)}", file=sys.stderr)
        print("  Hint: check ND_matrix header sample names match your bam basenames.", file=sys.stderr)
        sys.exit(1)

    # Warn if some missing
    miss_bsb = [s for s in bsb_samples if s not in bsb_cols]
    miss_tc = [s for s in tc_samples if s not in tc_cols]
    if miss_bsb:
        print(f"[WARN] {len(miss_bsb)} BSB samples not found in matrix header, e.g. {miss_bsb[:5]}", file=sys.stderr)
    if miss_tc:
        print(f"[WARN] {len(miss_tc)} TC samples not found in matrix header, e.g. {miss_tc[:5]}", file=sys.stderr)

    # -----------------------
    # Compute group means, log2FC
    # -----------------------
    Xb = df[bsb_cols].to_numpy(dtype=float)
    Xt = df[tc_cols].to_numpy(dtype=float)

    meanB = np.nanmean(Xb, axis=1)
    meanT = np.nanmean(Xt, axis=1)

    eps = 1e-6
    log2fc = np.log2((meanB + eps) / (meanT + eps))

    # Welch t-test per row (vectorized)
    # Note: scipy ttest_ind supports axis
    t_res = ttest_ind(Xb, Xt, axis=1, equal_var=False, nan_policy="omit")
    pval = t_res.pvalue

    # BH-FDR
    fdr = multipletests(pval, alpha=args.fdr, method="fdr_bh")[1]

    # Call by depth
    call = np.full(df.shape[0], "Mixed", dtype=object)
    call[(fdr < args.fdr) & (log2fc >= args.log2fc)] = "BSB"
    call[(fdr < args.fdr) & (log2fc <= -args.log2fc)] = "TC"

    out_depth = df[["chr", "start", "end"]].copy()
    out_depth["meanND_BSB"] = meanB
    out_depth["meanND_TC"] = meanT
    out_depth["log2FC"] = log2fc
    out_depth["pvalue"] = pval
    out_depth["FDR"] = fdr
    out_depth["depth_call"] = call

    depth_tsv = f"{args.out_prefix}.depth_validate.tsv"
    out_depth.to_csv(depth_tsv, sep="\t", index=False)
    print(f"[DONE] wrote {depth_tsv}", file=sys.stderr)

    # -----------------------
    # Merge with k-mer window calls
    # -----------------------
    kmer = pd.read_csv(args.kmer_source, sep="\t", dtype={"chr": str})
    if args.only_chr:
        kmer = kmer[kmer["chr"].str.match(r"^chr\d{2}[AB]$")].copy()

    # robustly find kmer call column (usually win_call)
    kcol = None
    for cand in ["win_call", "kmer_call", "call"]:
        if cand in kmer.columns:
            kcol = cand
            break
    if kcol is None:
        # try last column
        kcol = kmer.columns[-1]

    kmer_slim = kmer[["chr", "start", "end", kcol]].copy()
    kmer_slim["kmer_call"] = kmer_slim[kcol].astype(str).map(kmer_call_to_simple)

    merged = out_depth.merge(kmer_slim[["chr", "start", "end", "kmer_call"]],
                             on=["chr", "start", "end"], how="inner")

    merged_tsv = f"{args.out_prefix}.depth_vs_kmer.tsv"
    merged.to_csv(merged_tsv, sep="\t", index=False)
    print(f"[DONE] wrote {merged_tsv}", file=sys.stderr)

    # -----------------------
    # Summary: agreement + confusion matrix
    # -----------------------
    same = (merged["depth_call"] == merged["kmer_call"]).sum()
    tot = merged.shape[0]
    rate = same / tot if tot else 0.0

    # confusion
    cm = (merged
          .groupby(["depth_call", "kmer_call"])
          .size()
          .reset_index(name="count")
          .sort_values(["depth_call", "kmer_call"]))

    summary_txt = f"{args.out_prefix}.summary.txt"
    with open(summary_txt, "w") as o:
        o.write(f"ND matrix: {args.nd_matrix}\n")
        o.write(f"Kmer source: {args.kmer_source}\n")
        o.write(f"BSB n={len(bsb_cols)}; TC n={len(tc_cols)}\n")
        o.write(f"FDR<{args.fdr}, |log2FC|>={args.log2fc}\n")
        o.write(f"Merged windows: {tot}\n")
        o.write(f"Agreement: {same}/{tot} = {rate:.4f}\n\n")
        o.write("Confusion matrix (depth_call vs kmer_call):\n")
        for _, r in cm.iterrows():
            o.write(f"{r['depth_call']}\t{r['kmer_call']}\t{int(r['count'])}\n")

    print(f"[DONE] wrote {summary_txt}", file=sys.stderr)


if __name__ == "__main__":
    main()
