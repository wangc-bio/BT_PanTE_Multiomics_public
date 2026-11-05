#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Generalized from the supplied v5 reciprocal LTR evidence script.

from __future__ import print_function

import argparse
import os
import csv
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description="Re-evaluate RepeatModeler-to-LTR_retriever matches using reciprocal interval coverage.")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def read_manifest(path):
    base = os.path.dirname(os.path.abspath(path))
    entries, seen = [], set()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "ltr_raw_blast_tsv"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Manifest must contain sample_id and ltr_raw_blast_tsv columns")
        for line_no, row in enumerate(reader, 2):
            sample = (row.get("sample_id") or "").strip()
            raw = (row.get("ltr_raw_blast_tsv") or "").strip()
            if not sample or not raw:
                raise ValueError("Missing sample_id or ltr_raw_blast_tsv at line %d" % line_no)
            if sample in seen:
                raise ValueError("Duplicate sample_id: %s" % sample)
            seen.add(sample)
            if not os.path.isabs(raw):
                raw = os.path.join(base, raw)
            raw = os.path.abspath(raw)
            if not os.path.isfile(raw):
                raise FileNotFoundError("LTR raw BLAST table not found: %s" % raw)
            entries.append((sample, raw))
    if not entries:
        raise ValueError("The LTR manifest contains no entries")
    return entries


ARGS = parse_args()
ENTRIES = read_manifest(ARGS.manifest)
OUT = os.path.abspath(ARGS.output_dir)
os.makedirs(OUT, exist_ok=True)


def family_id(x):
    return x.strip().lstrip(">").split("#", 1)[0]


def target_class(x):
    if "#" in x:
        return x.split("#", 1)[1]
    return "LTR/Unknown"


def target_region(x):
    name = x.split("#", 1)[0]

    if name.endswith("_LTR"):
        return "LTR"

    if name.endswith("_INT"):
        return "INT"

    return "OTHER"


def union_length(intervals):
    if not intervals:
        return 0

    norm = []

    for a, b in intervals:
        lo = min(a, b)
        hi = max(a, b)
        norm.append((lo, hi))

    norm.sort()

    total = 0
    cur_s, cur_e = norm[0]

    for s, e in norm[1:]:

        if s <= cur_e + 1:
            if e > cur_e:
                cur_e = e
        else:
            total += cur_e - cur_s + 1
            cur_s, cur_e = s, e

    total += cur_e - cur_s + 1

    return total


def reciprocal_tier(pident, covered_bp, qcov, scov, evalue):

    if (
        pident >= 80.0
        and qcov >= 80.0
        and scov >= 80.0
        and covered_bp >= 100
        and evalue <= 1e-20
    ):
        return "HC"

    if (
        pident >= 70.0
        and qcov >= 50.0
        and scov >= 50.0
        and covered_bp >= 100
        and evalue <= 1e-10
    ):
        return "SUPPORTIVE"

    return "WEAK"


def evidence_level(region, tier):

    if region == "LTR":

        if tier == "HC":
            return "DIRECT_LTR_HC"

        if tier == "SUPPORTIVE":
            return "DIRECT_LTR_SUPPORTIVE"

        return "WEAK"

    if region == "INT":

        if tier == "HC":
            return "INTERNAL_HC"

        if tier == "SUPPORTIVE":
            return "INTERNAL_SUPPORTIVE"

        return "WEAK"

    return "WEAK"


EVIDENCE_RANK = {
    "DIRECT_LTR_HC": 6,
    "DIRECT_LTR_SUPPORTIVE": 5,
    "INTERNAL_HC": 4,
    "INTERNAL_SUPPORTIVE": 3,
    "WEAK": 1
}


summary_rows = []


for sample, infile in ENTRIES:

    # group by query-target pair
    pair_hits = defaultdict(list)

    raw_n = 0

    with open(infile, "r") as fh:

        for line_no, line in enumerate(fh, 1):

            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            f = line.split("\t")

            if len(f) < 13:
                raise RuntimeError(
                    "%s line %d: expected >=13 columns"
                    % (infile, line_no)
                )

            qseqid = f[0]
            sseqid = f[1]

            hit = {
                "qseqid": qseqid,
                "sseqid": sseqid,
                "pident": float(f[2]),
                "aln_len": int(float(f[3])),
                "qstart": int(float(f[4])),
                "qend": int(float(f[5])),
                "qlen": int(float(f[6])),
                "sstart": int(float(f[7])),
                "send": int(float(f[8])),
                "slen": int(float(f[9])),
                "evalue": float(f[10]),
                "bitscore": float(f[11])
            }

            pair_hits[(qseqid, sseqid)].append(hit)
            raw_n += 1

    pair_rows = []

    for (qseqid, sseqid), hits in pair_hits.items():

        qlen = hits[0]["qlen"]
        slen = hits[0]["slen"]

        q_intervals = [
            (h["qstart"], h["qend"])
            for h in hits
        ]

        s_intervals = [
            (h["sstart"], h["send"])
            for h in hits
        ]

        q_union_bp = union_length(q_intervals)
        s_union_bp = union_length(s_intervals)

        qcov = (
            100.0 * q_union_bp / float(qlen)
            if qlen > 0 else 0.0
        )

        scov = (
            100.0 * s_union_bp / float(slen)
            if slen > 0 else 0.0
        )

        qcov = min(qcov, 100.0)
        scov = min(scov, 100.0)

        # approximate weighted identity across HSPs
        total_aln = sum(h["aln_len"] for h in hits)

        if total_aln > 0:
            weighted_pident = sum(
                h["pident"] * h["aln_len"]
                for h in hits
            ) / float(total_aln)
        else:
            weighted_pident = 0.0

        # conservative aligned coverage length
        covered_bp = min(
            q_union_bp,
            s_union_bp
        )

        best_evalue = min(
            h["evalue"]
            for h in hits
        )

        bitscore_sum = sum(
            h["bitscore"]
            for h in hits
        )

        region = target_region(sseqid)

        tier = reciprocal_tier(
            weighted_pident,
            covered_bp,
            qcov,
            scov,
            best_evalue
        )

        evidence = evidence_level(
            region,
            tier
        )

        pair_rows.append({
            "RepeatModeler_ID":
                family_id(qseqid),

            "Query_raw":
                qseqid,

            "LTR_target":
                sseqid,

            "LTR_target_class":
                target_class(sseqid),

            "Target_region":
                region,

            "HSP_count":
                len(hits),

            "weighted_pident":
                round(weighted_pident, 3),

            "qlen":
                qlen,

            "q_union_bp":
                q_union_bp,

            "qcov":
                round(qcov, 3),

            "slen":
                slen,

            "s_union_bp":
                s_union_bp,

            "scov":
                round(scov, 3),

            "covered_bp":
                covered_bp,

            "best_evalue":
                best_evalue,

            "bitscore_sum":
                round(bitscore_sum, 3),

            "Reciprocal_tier":
                tier,

            "LTR_evidence_level":
                evidence
        })

    # best target for each RM family
    best_family = {}

    for r in pair_rows:

        fid = r["RepeatModeler_ID"]

        key = (
            EVIDENCE_RANK[
                r["LTR_evidence_level"]
            ],
            min(
                float(r["qcov"]),
                float(r["scov"])
            ),
            float(r["weighted_pident"]),
            int(r["covered_bp"]),
            float(r["bitscore_sum"])
        )

        if fid not in best_family:
            best_family[fid] = (key, r)

        elif key > best_family[fid][0]:
            best_family[fid] = (key, r)

    family_rows = [
        x[1]
        for x in best_family.values()
    ]

    family_rows.sort(
        key=lambda x: x["RepeatModeler_ID"]
    )

    fields = [
        "RepeatModeler_ID",
        "Query_raw",
        "LTR_target",
        "LTR_target_class",
        "Target_region",
        "HSP_count",
        "weighted_pident",
        "qlen",
        "q_union_bp",
        "qcov",
        "slen",
        "s_union_bp",
        "scov",
        "covered_bp",
        "best_evalue",
        "bitscore_sum",
        "Reciprocal_tier",
        "LTR_evidence_level"
    ]

    pair_out = os.path.join(
        OUT,
        "%s.LTR.reciprocal.pairs.tsv" % sample
    )

    with open(pair_out, "w") as fo:

        writer = csv.DictWriter(
            fo,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writeheader()
        writer.writerows(pair_rows)

    family_out = os.path.join(
        OUT,
        "%s.LTR.reciprocal.family_best.tsv"
        % sample
    )

    with open(family_out, "w") as fo:

        writer = csv.DictWriter(
            fo,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writeheader()
        writer.writerows(family_rows)

    counts = defaultdict(int)

    for r in family_rows:
        counts[r["LTR_evidence_level"]] += 1

    summary_rows.append({
        "Sample": sample,
        "Raw_HSPs": raw_n,
        "Query_target_pairs": len(pair_rows),
        "RM_families_with_LTRretriever_hit":
            len(family_rows),
        "DIRECT_LTR_HC":
            counts["DIRECT_LTR_HC"],
        "DIRECT_LTR_SUPPORTIVE":
            counts["DIRECT_LTR_SUPPORTIVE"],
        "INTERNAL_HC":
            counts["INTERNAL_HC"],
        "INTERNAL_SUPPORTIVE":
            counts["INTERNAL_SUPPORTIVE"],
        "WEAK_ONLY":
            counts["WEAK"]
    })

    print(
        "%s\tfamilies=%d\tDIRECT_HC=%d\tDIRECT_SUP=%d\t"
        "INT_HC=%d\tINT_SUP=%d\tWEAK=%d"
        % (
            sample,
            len(family_rows),
            counts["DIRECT_LTR_HC"],
            counts["DIRECT_LTR_SUPPORTIVE"],
            counts["INTERNAL_HC"],
            counts["INTERNAL_SUPPORTIVE"],
            counts["WEAK"]
        )
    )


summary_file = os.path.join(
    OUT,
    "all_libraries.LTR_reciprocal.summary.tsv"
)

summary_fields = [
    "Sample",
    "Raw_HSPs",
    "Query_target_pairs",
    "RM_families_with_LTRretriever_hit",
    "DIRECT_LTR_HC",
    "DIRECT_LTR_SUPPORTIVE",
    "INTERNAL_HC",
    "INTERNAL_SUPPORTIVE",
    "WEAK_ONLY"
]

with open(summary_file, "w") as fo:

    writer = csv.DictWriter(
        fo,
        fieldnames=summary_fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(summary_rows)

print("")
print("Finished; libraries processed:", len(ENTRIES))
print(summary_file)