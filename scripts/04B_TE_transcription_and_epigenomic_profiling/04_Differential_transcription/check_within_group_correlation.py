#!/usr/bin/env python
# Public-release note: original executed source; replace local absolute paths with paths from your environment before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import pandas as pd
import numpy as np

ROOT = os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52", "05_TE_RNA")

META = os.path.join(
    ROOT,
    "00_metadata/sample_metadata.tsv"
)

meta = pd.read_csv(META, sep="\t")
meta = meta.set_index("Sample")

tissues = ["Brain", "Gill", "Kidney", "Liver", "Muscle"]

all_rows = []

for tissue in tissues:

    corfile = os.path.join(
        ROOT,
        "05_DE/03_DESeq2",
        tissue,
        tissue + ".sample_correlation.tsv"
    )

    cor = pd.read_csv(
        corfile,
        sep="\t",
        index_col=0
    )

    samples = list(cor.columns)

    for sample in samples:

        group = meta.loc[sample, "Group"]

        peers = [
            x for x in samples
            if x != sample
            and meta.loc[x, "Group"] == group
        ]

        if len(peers) == 0:
            mean_cor = np.nan
            min_cor = np.nan
        else:
            vals = cor.loc[sample, peers].astype(float)
            mean_cor = vals.mean()
            min_cor = vals.min()

        all_rows.append({
            "Tissue": tissue,
            "Sample": sample,
            "Group": group,
            "N_peers": len(peers),
            "Mean_within_group_Pearson": mean_cor,
            "Min_within_group_Pearson": min_cor
        })

out = pd.DataFrame(all_rows)

outfile = os.path.join(
    ROOT,
    "05_DE/00_QC/within_group_sample_correlation.tsv"
)

out.to_csv(
    outfile,
    sep="\t",
    index=False
)

print(out.to_string(index=False))
print("")
print("Output:", outfile)