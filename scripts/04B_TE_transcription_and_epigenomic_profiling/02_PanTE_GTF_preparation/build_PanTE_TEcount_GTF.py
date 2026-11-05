#!/usr/bin/env python3
import os
# Public-release note: original executed source; replace local absolute paths with paths from your environment before running.
# -*- coding: utf-8 -*-

import re
import sys
from pathlib import Path
from collections import Counter, defaultdict

import pandas as pd


# ============================================================
# PATHS
# ============================================================

BASE = Path(os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE"))
PANTE = BASE / "PanTE_v52"

MEMBERSHIP = (
    PANTE /
    "02_cluster/PanTE_refined_cluster_members.tsv"
)

UNIFIED_FA = (
    BASE /
    "Final_evidence_v5/13_unified_v52/"
    "Unified_v5_2.curated_repeat.fa"
)

GFF_FILES = {
    "BSB": (
        BASE /
        "RepeatMasker_unified_v5_2_fixed/results/"
        "BSB/BSB.genome.chr.fa.out.gff"
    ),

    "TC": (
        BASE /
        "RepeatMasker_unified_v5_2_fixed/results/"
        "TC/TC.genome.chr.fa.out.gff"
    ),

    "BTM": (
        BASE /
        "RepeatMasker_unified_v5_2_fixed/results/"
        "BTM/BT.genome.subM.fa.out.gff"
    ),

    "BTP": (
        BASE /
        "RepeatMasker_unified_v5_2_fixed/results/"
        "BTP/BT.genome.subP.fa.out.gff"
    ),
}

OUTDIR = (
    PANTE /
    "05_TE_RNA/01_annotation"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# HELPERS
# ============================================================

def parse_unified_headers(fasta):
    """
    Expected header examples:
      >UNI_v52_xxxxxx#LTR/Gypsy
      >UNI_v52_xxxxxx#DNA/Unknown

    Only the first whitespace-delimited token is used.
    """

    anno = {}

    with open(fasta) as fh:

        for line in fh:

            if not line.startswith(">"):
                continue

            header = line[1:].strip()
            token = header.split()[0]

            if "#" in token:

                uid, cls = token.split("#", 1)

            else:

                uid = token
                cls = "Unknown"

            cls = cls.strip()

            if "/" in cls:

                broad, subclass = cls.split("/", 1)

            else:

                broad = cls
                subclass = "Unknown"

            anno[uid] = {
                "header_class": cls,
                "header_broad": broad,
                "header_subclass": subclass,
            }

    return anno


def load_membership(path):

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str
    )

    required = [
        "PanTE_ID",
        "Broad_class",
        "Representative_ID",
        "Unified_ID"
    ]

    for c in required:

        if c not in df.columns:
            raise ValueError(
                "Missing membership column: {}".format(c)
            )

    # Every Unified_ID should belong to one PanTE
    if df["Unified_ID"].duplicated().any():

        bad = df.loc[
            df["Unified_ID"].duplicated(False),
            "Unified_ID"
        ].head(20).tolist()

        raise ValueError(
            "Duplicated Unified_ID in membership: {}".format(
                ",".join(bad)
            )
        )

    return df


def extract_motif(attributes):
    """
    RepeatMasker GFF example:
    ID=1;Target "Motif:UNI_v52_011087" 5902 6146
    """

    m = re.search(
        r'Target\s+"Motif:([^"]+)"',
        attributes
    )

    if m:
        return m.group(1)

    return None


def clean_value(x):

    x = str(x)

    x = x.replace('"', "")
    x = x.replace(";", "_")
    x = x.replace(" ", "_")

    return x


# ============================================================
# BUILD ANNOTATION TABLES
# ============================================================

print("===== Loading PanTE membership =====")

mem = load_membership(
    MEMBERSHIP
)

print("Membership rows :", len(mem))
print("Unique PanTE    :", mem["PanTE_ID"].nunique())
print("Unified IDs     :", mem["Unified_ID"].nunique())


print("")
print("===== Parsing Unified FASTA headers =====")

header_anno = parse_unified_headers(
    UNIFIED_FA
)

print("Unified FASTA headers:", len(header_anno))


# Unified ID -> PanTE metadata

uid_to_meta = {}

# PanTE representative annotation
pante_rep_anno = {}


for _, r in mem.iterrows():

    uid = r["Unified_ID"]
    pante = r["PanTE_ID"]
    broad = r["Broad_class"]
    rep = r["Representative_ID"]

    if rep in header_anno:

        rep_broad = header_anno[rep]["header_broad"]
        rep_subclass = header_anno[rep]["header_subclass"]
        rep_full = header_anno[rep]["header_class"]

    else:

        rep_broad = broad
        rep_subclass = "Unknown"
        rep_full = broad + "/Unknown"

    pante_rep_anno[pante] = {
        "Broad_class": broad,
        "Representative_ID": rep,
        "Representative_broad": rep_broad,
        "Subclass": rep_subclass,
        "Full_class": rep_full,
    }

    uid_to_meta[uid] = {
        "PanTE_ID": pante,
        "Broad_class": broad,
        "Representative_ID": rep,
        "Subclass": rep_subclass,
        "Full_class": rep_full,
    }


# ============================================================
# WRITE PanTE master annotation
# ============================================================

master_rows = []

for pante in sorted(pante_rep_anno):

    x = pante_rep_anno[pante]

    master_rows.append({
        "PanTE_ID": pante,
        "Broad_class": x["Broad_class"],
        "Subclass": x["Subclass"],
        "Full_class": x["Full_class"],
        "Representative_ID": x["Representative_ID"],
    })

master = pd.DataFrame(
    master_rows
)

master.to_csv(
    OUTDIR / "PanTE_master_annotation.tsv",
    sep="\t",
    index=False
)

print("")
print("PanTE master annotation written:")
print(OUTDIR / "PanTE_master_annotation.tsv")


# ============================================================
# GFF -> TEcount GTF
# ============================================================

def convert_one(
    label,
    gff,
    output_gtf,
    locus_prefix,
    origin=None,
    append=False,
    origin_fh=None
):

    mode = "a" if append else "w"

    stats = Counter()

    pantes_seen = set()
    chromosomes = set()
    subclass_counter = Counter()
    broad_counter = Counter()

    locus_index = 0

    with open(gff) as inp, open(output_gtf, mode) as out:

        for line in inp:

            if line.startswith("#"):
                continue

            f = line.rstrip("\n").split("\t")

            if len(f) < 9:
                stats["malformed"] += 1
                continue

            stats["gff_records"] += 1

            chrom = f[0]
            feature = f[2]
            start = f[3]
            end = f[4]
            score = f[5]
            strand = f[6]
            attrs = f[8]

            motif = extract_motif(attrs)

            if motif is None:

                stats["no_motif"] += 1
                continue

            # Keep only unified classified TE consensuses
            if not motif.startswith("UNI_v52_"):

                stats["non_UNI"] += 1
                continue

            stats["UNI_records"] += 1

            if motif not in uid_to_meta:

                stats["UNI_not_in_PanTE"] += 1
                continue

            meta = uid_to_meta[motif]

            pante = meta["PanTE_ID"]
            broad = meta["Broad_class"]
            subclass = meta["Subclass"]

            locus_index += 1

            locus_id = "{}_{:09d}".format(
                locus_prefix,
                locus_index
            )

            family = clean_value(
                subclass
            )

            broad_clean = clean_value(
                broad
            )

            attrs_gtf = (
                'gene_id "{}"; '
                'transcript_id "{}"; '
                'family_id "{}"; '
                'class_id "{}";'
            ).format(
                pante,
                locus_id,
                family,
                broad_clean
            )

            gtf_fields = [
                chrom,
                "PanTE",
                "exon",
                start,
                end,
                score,
                strand,
                ".",
                attrs_gtf,
            ]

            out.write(
                "\t".join(gtf_fields) + "\n"
            )

            stats["kept"] += 1
            chromosomes.add(chrom)
            pantes_seen.add(pante)

            broad_counter[broad] += 1
            subclass_counter[
                broad + "/" + subclass
            ] += 1

            if origin_fh is not None:

                origin_fh.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        locus_id,
                        origin if origin else label,
                        chrom,
                        start,
                        end,
                        motif,
                        pante,
                        broad + "/" + subclass
                    )
                )

    return {
        "stats": stats,
        "pantes": pantes_seen,
        "chromosomes": chromosomes,
        "broad": broad_counter,
        "subclass": subclass_counter,
    }


# ============================================================
# BSB
# ============================================================

print("")
print("===== Building BSB TE GTF =====")

bsb_gtf = (
    OUTDIR /
    "BSB.PanTE.TEtranscripts.gtf"
)

bsb_origin = (
    OUTDIR /
    "BSB.PanTE.locus_metadata.tsv"
)

with open(bsb_origin, "w") as md:

    md.write(
        "Locus_ID\tOrigin\tChr\tStart\tEnd\tUnified_ID\t"
        "PanTE_ID\tClassification\n"
    )

    bsb_res = convert_one(
        "BSB",
        GFF_FILES["BSB"],
        bsb_gtf,
        "BSB_TE",
        origin="BSB",
        append=False,
        origin_fh=md
    )


# ============================================================
# TC
# ============================================================

print("===== Building TC TE GTF =====")

tc_gtf = (
    OUTDIR /
    "TC.PanTE.TEtranscripts.gtf"
)

tc_origin = (
    OUTDIR /
    "TC.PanTE.locus_metadata.tsv"
)

with open(tc_origin, "w") as md:

    md.write(
        "Locus_ID\tOrigin\tChr\tStart\tEnd\tUnified_ID\t"
        "PanTE_ID\tClassification\n"
    )

    tc_res = convert_one(
        "TC",
        GFF_FILES["TC"],
        tc_gtf,
        "TC_TE",
        origin="TC",
        append=False,
        origin_fh=md
    )


# ============================================================
# BT = BTM + BTP
# ============================================================

print("===== Building BT merged TE GTF =====")

bt_gtf = (
    OUTDIR /
    "BT.PanTE.TEtranscripts.gtf"
)

bt_origin = (
    OUTDIR /
    "BT.PanTE.locus_metadata.tsv"
)

with open(bt_origin, "w") as md:

    md.write(
        "Locus_ID\tOrigin\tChr\tStart\tEnd\tUnified_ID\t"
        "PanTE_ID\tClassification\n"
    )

    btm_res = convert_one(
        "BTM",
        GFF_FILES["BTM"],
        bt_gtf,
        "BTM_TE",
        origin="M",
        append=False,
        origin_fh=md
    )

    btp_res = convert_one(
        "BTP",
        GFF_FILES["BTP"],
        bt_gtf,
        "BTP_TE",
        origin="P",
        append=True,
        origin_fh=md
    )


# ============================================================
# QC
# ============================================================

def write_qc(name, res, outfile):

    with open(outfile, "w") as out:

        s = res["stats"]

        out.write("Sample\t{}\n".format(name))

        for k in sorted(s):

            out.write(
                "{}\t{}\n".format(
                    k,
                    s[k]
                )
            )

        out.write(
            "unique_chromosomes\t{}\n".format(
                len(res["chromosomes"])
            )
        )

        out.write(
            "unique_PanTE\t{}\n".format(
                len(res["pantes"])
            )
        )

        out.write("\n# Broad_class loci\n")

        for k, v in sorted(
            res["broad"].items()
        ):

            out.write(
                "{}\t{}\n".format(k, v)
            )

        out.write("\n# Subclass loci\n")

        for k, v in sorted(
            res["subclass"].items(),
            key=lambda x: (-x[1], x[0])
        ):

            out.write(
                "{}\t{}\n".format(k, v)
            )


write_qc(
    "BSB",
    bsb_res,
    OUTDIR / "BSB.PanTE.GTF_QC.txt"
)

write_qc(
    "TC",
    tc_res,
    OUTDIR / "TC.PanTE.GTF_QC.txt"
)


# Merge BT stats for report

bt_res = {
    "stats": btm_res["stats"] + btp_res["stats"],
    "pantes": btm_res["pantes"] | btp_res["pantes"],
    "chromosomes": (
        btm_res["chromosomes"]
        | btp_res["chromosomes"]
    ),
    "broad": btm_res["broad"] + btp_res["broad"],
    "subclass": (
        btm_res["subclass"]
        + btp_res["subclass"]
    ),
}

write_qc(
    "BT",
    bt_res,
    OUTDIR / "BT.PanTE.GTF_QC.txt"
)


# ============================================================
# FINAL SUMMARY
# ============================================================

print("")
print("===== FINAL SUMMARY =====")

for name, res, gtf in [
    ("BSB", bsb_res, bsb_gtf),
    ("TC", tc_res, tc_gtf),
    ("BT", bt_res, bt_gtf),
]:

    print(
        "{}: loci={}, PanTE={}, chromosomes={}".format(
            name,
            res["stats"]["kept"],
            len(res["pantes"]),
            len(res["chromosomes"])
        )
    )

    print("  GTF:", gtf)

print("")
print("Done.")