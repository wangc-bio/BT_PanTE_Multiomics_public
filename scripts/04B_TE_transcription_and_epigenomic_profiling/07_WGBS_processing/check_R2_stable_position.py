#!/usr/bin/env python3
import os
# Public-release note: original executed source; replace local absolute paths with local paths before running.

import pandas as pd
import numpy as np

ROOT = os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52", "06_TE_WGBS", "00_QC")

POS = f"{ROOT}/M_bias_CpG_all_positions.tsv"
SUM = f"{ROOT}/M_bias_CpG_sample_summary.tsv"

OUT = f"{ROOT}/R2_position_stability.tsv"
OUT_SAMPLE = f"{ROOT}/R2_sample_stable_position.tsv"

pos = pd.read_csv(POS, sep="\t")
sm  = pd.read_csv(SUM, sep="\t")

# ??? R2
pos = pos[pos["read"] == "R2"].copy()
sm  = sm[sm["read"] == "R2"].copy()

central = dict(
    zip(
        sm["sample_id"],
        pd.to_numeric(sm["central_weighted_mCG"])
    )
)

pos["central_mCG"] = pos["sample_id"].map(central)
pos["mCG_percent"] = pd.to_numeric(pos["mCG_percent"])
pos["coverage"] = pd.to_numeric(pos["coverage"])

pos["diff_pp"] = (
    pos["mCG_percent"] -
    pos["central_mCG"]
)

pos["abs_diff_pp"] = pos["diff_pp"].abs()

# ------------------------------------------------------------
# 1. ????? position ??
#    ??? 20 bp
# ------------------------------------------------------------

x = pos[
    (pos["position"] >= 1) &
    (pos["position"] <= 20)
].copy()

rows = []

for p, g in x.groupby("position"):

    vals = g["diff_pp"].dropna()
    av = vals.abs()

    rows.append({
        "position": int(p),
        "n_samples": len(vals),
        "median_diff_pp": np.median(vals),
        "median_abs_diff_pp": np.median(av),
        "P90_abs_diff_pp": np.percentile(av, 90),
        "max_abs_diff_pp": np.max(av),
        "n_abs_le_2pp": int((av <= 2).sum()),
        "pct_abs_le_2pp": 100.0 * (av <= 2).mean(),
        "n_abs_le_5pp": int((av <= 5).sum()),
        "pct_abs_le_5pp": 100.0 * (av <= 5).mean(),
    })

out = pd.DataFrame(rows).sort_values("position")
out.to_csv(OUT, sep="\t", index=False)


# ------------------------------------------------------------
# 2. ??????"???????"
#
# ????? QC heuristic?
# ????????? 5 bp?
# ??? mCG ??? <= 2 percentage points
# ------------------------------------------------------------

sample_rows = []

for sample, g in x.groupby("sample_id"):

    g = g.sort_values("position")

    d = dict(
        zip(
            g["position"].astype(int),
            g["abs_diff_pp"]
        )
    )

    stable_start = None

    for start in range(1, 17):

        vals = [
            d.get(p, np.nan)
            for p in range(start, start + 5)
        ]

        if (
            all(np.isfinite(v) for v in vals)
            and all(v <= 2 for v in vals)
        ):
            stable_start = start
            break

    sample_rows.append({
        "sample_id": sample,
        "species": g["species"].iloc[0],
        "tissue": g["tissue"].iloc[0],
        "first_stable_position_5bp": (
            stable_start
            if stable_start is not None
            else "NOT_STABLE_BY_20BP"
        )
    })

sample_out = pd.DataFrame(sample_rows)
sample_out.to_csv(
    OUT_SAMPLE,
    sep="\t",
    index=False
)

print("===== R2 position stability =====")
print(out.to_string(index=False))

print("")
print("===== Stable-start distribution =====")
print(
    sample_out[
        "first_stable_position_5bp"
    ].value_counts(dropna=False).sort_index()
)

print("")
print("Outputs:")
print(OUT)
print(OUT_SAMPLE)