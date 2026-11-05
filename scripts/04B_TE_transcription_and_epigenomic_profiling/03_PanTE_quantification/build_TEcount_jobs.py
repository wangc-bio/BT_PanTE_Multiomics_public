#!/usr/bin/env python3
import os
# Public-release note: original executed source; replace local absolute paths with paths from your environment before running.

import pandas as pd

meta_file = (
    os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52", "05_TE_RNA", "00_metadata", "sample_metadata.tsv")
)

out_file = (
    os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52", "05_TE_RNA", "00_metadata", "TEcount_jobs.tsv")
)

gene_gtf = {
    "BSB": os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "RNAseq", "BSB", "01_index", "BSB_longest_isoform.gtf"),
    "TC":  os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "RNAseq", "TC", "01_index", "TC_longest_isoform.gtf"),
    "BT":  os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "RNAseq", "BT", "index", "BT.rename.longest_isoform.gtf"),
}

te_gtf = {
    "BSB": os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52", "05_TE_RNA", "01_annotation", "BSB.PanTE.TEtranscripts.gtf"),
    "TC":  os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52", "05_TE_RNA", "01_annotation", "TC.PanTE.TEtranscripts.gtf"),
    "BT":  os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52", "05_TE_RNA", "01_annotation", "BT.PanTE.TEtranscripts.gtf"),
}

df = pd.read_csv(meta_file, sep="\t")

df["Gene_GTF"] = df["Group"].map(gene_gtf)
df["TE_GTF"] = df["Group"].map(te_gtf)

if df["Gene_GTF"].isna().any():
    raise ValueError("Missing Gene_GTF mapping")

if df["TE_GTF"].isna().any():
    raise ValueError("Missing TE_GTF mapping")

# SLURM array ? 1-based Task_ID
df.insert(0, "Task_ID", range(1, len(df) + 1))

cols = [
    "Task_ID",
    "Sample",
    "Group",
    "Individual",
    "Tissue",
    "BAM",
    "Gene_GTF",
    "TE_GTF",
]

df[cols].to_csv(out_file, sep="\t", index=False)

print("Total jobs:", len(df))
print("")
print(df.groupby(["Group", "Tissue"]).size())
print("")
print("Output:", out_file)