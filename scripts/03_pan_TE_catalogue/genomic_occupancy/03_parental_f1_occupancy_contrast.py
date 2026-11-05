#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Compare parental Pan-TE occupancy contrasts with the corresponding F1 haplotypes."""

import argparse
import csv
import math
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


REQUIRED_ROLES = {
    "parent_A",
    "parent_B",
    "f1_haplotype_A",
    "f1_haplotype_B",
}


def read_role_manifest(path):
    roles = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"role", "sample_id"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                "Role manifest must contain tab-separated columns: role, sample_id"
            )

        for line_number, row in enumerate(reader, 2):
            role = (row.get("role") or "").strip()
            sample = (row.get("sample_id") or "").strip()
            if role not in REQUIRED_ROLES:
                raise ValueError(
                    "Unknown role '%s' at manifest line %d" % (role, line_number)
                )
            if not sample:
                raise ValueError(
                    "Empty sample_id at manifest line %d" % line_number
                )
            if role in roles:
                raise ValueError("Duplicate role: %s" % role)
            if not re.match(r"^[A-Za-z0-9_.-]+$", sample):
                raise ValueError(
                    "sample_id must contain only letters, numbers, '.', '_' or '-': "
                    + sample
                )
            roles[role] = sample

    missing = REQUIRED_ROLES - set(roles)
    if missing:
        raise ValueError(
            "Role manifest is missing roles: %s" % ", ".join(sorted(missing))
        )
    if len(set(roles.values())) != 4:
        raise ValueError("The four analysis roles must map to four distinct sample IDs.")
    return roles


def normalized_name(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def resolve_matrix_column(df, sample_id, suffix):
    target = normalized_name(sample_id + suffix)
    matches = [
        column for column in df.columns
        if normalized_name(column) == target
    ]
    if len(matches) != 1:
        raise KeyError(
            "Cannot uniquely resolve matrix column for sample '%s', suffix '%s'. "
            "Candidates: %s" % (sample_id, suffix, matches)
        )
    return matches[0]


def safe_spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]

    if len(x) < 3 or np.all(x == x[0]) or np.all(y == y[0]):
        return np.nan, np.nan, len(x)

    result = spearmanr(x, y)
    rho = getattr(result, "statistic", None)
    if rho is None:
        rho = getattr(result, "correlation", np.nan)
    pvalue = getattr(result, "pvalue", np.nan)
    return float(rho), float(pvalue), len(x)


def ratio_without_pseudocount(numerator, denominator):
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    ratio = np.full(numerator.shape, np.nan, dtype=float)

    both_positive = (numerator > 0) & (denominator > 0)
    ratio[both_positive] = numerator[both_positive] / denominator[both_positive]
    ratio[(numerator > 0) & (denominator == 0)] = np.inf
    ratio[(numerator == 0) & (denominator > 0)] = 0.0
    return ratio


def classify_parent_bias(ratio, fold_threshold):
    labels = np.full(len(ratio), "similar", dtype=object)
    labels[ratio >= fold_threshold] = "parent_A-biased"
    labels[ratio <= 1.0 / fold_threshold] = "parent_B-biased"
    labels[np.isnan(ratio)] = "uninformative"
    return labels


def directional_retention(parent_bias, f1_a, f1_b):
    if parent_bias == "parent_A-biased":
        if f1_a > f1_b:
            return "concordant"
        if f1_a < f1_b:
            return "discordant"
        return "tie"

    if parent_bias == "parent_B-biased":
        if f1_b > f1_a:
            return "concordant"
        if f1_b < f1_a:
            return "discordant"
        return "tie"

    return "not_assessed"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Classify parent-informative and parent-biased Pan-TE clusters, "
            "then compare parental and F1 haplotypic occupancy contrasts."
        )
    )
    parser.add_argument("--matrix", required=True,
                        help="Pan-TE resolved-abundance matrix TSV")
    parser.add_argument("--role-manifest", required=True,
                        help="TSV mapping parent_A, parent_B and the two F1 haplotypes to matrix sample IDs")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--presence-bp", type=float, default=10000.0,
                        help="Raw occupied-bp threshold for parent-informative clusters (default: 10000)")
    parser.add_argument("--parent-fold", type=float, default=2.0,
                        help="Fold threshold for parental bias based on bp/Mb (default: 2)")
    parser.add_argument("--pseudocount", type=float, default=1.0,
                        help="Pseudocount in bp/Mb for log2 contrasts only (default: 1)")
    args = parser.parse_args()

    if args.presence_bp < 0:
        raise ValueError("--presence-bp must be non-negative")
    if args.parent_fold <= 1:
        raise ValueError("--parent-fold must be greater than 1")
    if args.pseudocount <= 0:
        raise ValueError("--pseudocount must be positive")

    roles = read_role_manifest(args.role_manifest)
    df = pd.read_csv(args.matrix, sep="\t", low_memory=False)

    id_candidates = [
        column for column in df.columns
        if normalized_name(column) == "panteid"
    ]
    if len(id_candidates) != 1:
        raise KeyError("Matrix must contain exactly one PanTE_ID column.")
    id_col = id_candidates[0]

    if df[id_col].duplicated().any():
        raise ValueError("PanTE_ID values are not unique in the abundance matrix.")

    core = pd.DataFrame({"PanTE_ID": df[id_col].astype(str)})

    for metadata_col in [
        "Broad_class", "Cluster_size", "Representative_ID"
    ]:
        if metadata_col in df.columns:
            core[metadata_col] = df[metadata_col]

    role_columns = {}
    for role, sample in roles.items():
        bp_col = resolve_matrix_column(df, sample, "_bp")
        bpm_col = resolve_matrix_column(df, sample, "_bp_per_Mb")

        bp_name = role + "_occupied_bp"
        bpm_name = role + "_bp_per_Mb"

        core[bp_name] = pd.to_numeric(df[bp_col], errors="coerce")
        core[bpm_name] = pd.to_numeric(df[bpm_col], errors="coerce")

        if core[bp_name].isna().any() or core[bpm_name].isna().any():
            raise ValueError(
                "Missing or non-numeric occupancy values for role %s (%s)"
                % (role, sample)
            )
        if (core[bp_name] < 0).any() or (core[bpm_name] < 0).any():
            raise ValueError("Negative occupancy values found for role %s" % role)

        role_columns[role] = (bp_name, bpm_name)

    parent_a_bp, parent_a_bpm = role_columns["parent_A"]
    parent_b_bp, parent_b_bpm = role_columns["parent_B"]
    f1_a_bp, f1_a_bpm = role_columns["f1_haplotype_A"]
    f1_b_bp, f1_b_bpm = role_columns["f1_haplotype_B"]

    # Parent-informative status uses raw parental occupied base pairs.
    core["parent_informative"] = (
        (core[parent_a_bp] >= args.presence_bp)
        | (core[parent_b_bp] >= args.presence_bp)
    )

    # Parental bias uses normalized bp/Mb values without a pseudocount.
    parent_ratio = ratio_without_pseudocount(
        core[parent_a_bpm].to_numpy(),
        core[parent_b_bpm].to_numpy()
    )
    core["parent_A_over_parent_B_ratio_no_pc"] = parent_ratio
    parent_fold_difference = np.full(parent_ratio.shape, np.nan, dtype=float)
    positive_ratio = parent_ratio > 0
    parent_fold_difference[positive_ratio] = np.maximum(
        parent_ratio[positive_ratio],
        1.0 / parent_ratio[positive_ratio]
    )
    parent_fold_difference[(parent_ratio == 0) | np.isinf(parent_ratio)] = np.inf
    core["parent_fold_difference"] = parent_fold_difference
    core["parent_bias"] = classify_parent_bias(parent_ratio, args.parent_fold)
    core.loc[~core["parent_informative"], "parent_bias"] = (
        "low_abundance_unclassified"
    )

    # Log2 contrasts use the stated pseudocount in bp/Mb.
    pc = args.pseudocount
    core["X_log2_parent_A_over_parent_B"] = np.log2(
        (core[parent_a_bpm] + pc) / (core[parent_b_bpm] + pc)
    )
    core["Y_log2_F1_A_over_F1_B"] = np.log2(
        (core[f1_a_bpm] + pc) / (core[f1_b_bpm] + pc)
    )

    core["F1_directional_retention"] = [
        directional_retention(bias, a, b)
        for bias, a, b in zip(
            core["parent_bias"],
            core[f1_a_bpm],
            core[f1_b_bpm]
        )
    ]

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    metadata_cols = [
        column for column in [
            "PanTE_ID", "Broad_class", "Cluster_size", "Representative_ID"
        ] if column in core.columns
    ]
    value_cols = [
        parent_a_bp, parent_b_bp, f1_a_bp, f1_b_bp,
        parent_a_bpm, parent_b_bpm, f1_a_bpm, f1_b_bpm,
        "parent_informative",
        "parent_A_over_parent_B_ratio_no_pc",
        "parent_fold_difference",
        "parent_bias",
        "X_log2_parent_A_over_parent_B",
        "Y_log2_F1_A_over_F1_B",
        "F1_directional_retention",
    ]
    core[metadata_cols + value_cols].to_csv(
        os.path.join(outdir, "parent_f1_contrast_per_PanTE.tsv"),
        sep="\t", index=False
    )

    # Association and sign concordance are evaluated among parent-informative clusters.
    informative = core["parent_informative"]
    x = core.loc[informative, "X_log2_parent_A_over_parent_B"]
    y = core.loc[informative, "Y_log2_F1_A_over_F1_B"]
    rho, pvalue, n = safe_spearman(x, y)

    sign_keep = x.notna() & y.notna() & (x != 0) & (y != 0)
    sign_concordance = np.nan
    if sign_keep.any():
        sign_concordance = float(
            (np.sign(x[sign_keep]) == np.sign(y[sign_keep])).mean()
        )

    stats_rows = [{
        "analysis_set": "parent_informative",
        "n_PanTE": n,
        "pseudocount_bp_per_Mb": pc,
        "spearman_rho": rho,
        "spearman_pvalue_descriptive": pvalue,
        "n_nonzero_sign_pairs": int(sign_keep.sum()),
        "directional_sign_concordance": sign_concordance,
    }]

    # Directional retention among parent-biased clusters; ties are reported separately.
    biased = core[
        core["parent_informative"]
        & core["parent_bias"].isin(["parent_A-biased", "parent_B-biased"])
    ]
    direction_rows = []
    for label in ["parent_A-biased", "parent_B-biased", "combined_parent_biased"]:
        subset = biased if label == "combined_parent_biased" else biased[
            biased["parent_bias"] == label
        ]
        concordant = int((subset["F1_directional_retention"] == "concordant").sum())
        discordant = int((subset["F1_directional_retention"] == "discordant").sum())
        ties = int((subset["F1_directional_retention"] == "tie").sum())
        assessed = concordant + discordant
        direction_rows.append({
            "parent_bias": label,
            "n_PanTE": len(subset),
            "concordant": concordant,
            "discordant": discordant,
            "tie": ties,
            "directional_concordance_fraction_excluding_ties": (
                concordant / assessed if assessed else np.nan
            ),
        })

    pd.DataFrame(direction_rows).to_csv(
        os.path.join(outdir, "parent_biased_directional_retention.tsv"),
        sep="\t", index=False
    )

    parameters = [
        ("parent_A_sample_id", roles["parent_A"]),
        ("parent_B_sample_id", roles["parent_B"]),
        ("f1_haplotype_A_sample_id", roles["f1_haplotype_A"]),
        ("f1_haplotype_B_sample_id", roles["f1_haplotype_B"]),
        ("parent_informative_threshold_bp", args.presence_bp),
        ("parent_bias_fold_threshold", args.parent_fold),
        ("log2_contrast_pseudocount_bp_per_Mb", args.pseudocount),
    ]
    with open(os.path.join(outdir, "run_parameters.tsv"), "w",
              encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["parameter", "value"])
        writer.writerows(parameters)

    print("Pan-TE clusters:", len(core))
    print("Parent-informative clusters:", int(informative.sum()))
    print("Parent-biased clusters:", len(biased))
    print("Spearman rho (parent-informative):", rho)
    print("Output directory:", outdir)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("ERROR: %s" % error, file=sys.stderr)
        sys.exit(1)
