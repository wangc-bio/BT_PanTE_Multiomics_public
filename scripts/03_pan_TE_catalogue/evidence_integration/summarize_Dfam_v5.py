#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import sys
import csv
from collections import defaultdict

if len(sys.argv) != 5:
    raise RuntimeError(
        "Usage: summarize_Dfam_v5.py SAMPLE OUTDIR XENO_RAW CYPRINO_RAW"
    )

SAMPLE = sys.argv[1]
OUTDIR = os.path.abspath(sys.argv[2])
XENO_RAW = os.path.abspath(sys.argv[3])
CYPRINO_RAW = os.path.abspath(sys.argv[4])
os.makedirs(OUTDIR, exist_ok=True)


def family_id(x):
    return x.strip().lstrip(">").split("#", 1)[0]


def target_class(x):
    if "#" in x:
        c = x.split("#", 1)[1].strip()
        if c:
            return c
    return "Unknown"


def broad_class(c):
    s = c.lower()

    if "sine" in s:
        return "SINE"

    if (
        "ltr" in s or "gypsy" in s or "copia" in s or
        "bel" in s or "pao" in s or "erv" in s or "dirs" in s
    ):
        return "LTR"

    if (
        "line" in s or "l1" in s or "l2" in s or "cr1" in s or
        "rex" in s or "jockey" in s or "rte" in s or "bov-b" in s
    ):
        return "LINE"

    if "penelope" in s:
        return "PENELOPE"

    if "helitron" in s or s.startswith("rc/") or s == "rc":
        return "RC"

    if (
        "dna" in s or "tir" in s or "mite" in s or
        "maverick" in s or "polinton" in s
    ):
        return "DNA"

    if (
        "simple_repeat" in s or "simple repeat" in s or
        "satellite" in s or "low_complexity" in s or
        "low complexity" in s or s == "trna" or s == "rrna"
    ):
        return "nonTE"

    return "Unknown"


def union_length(intervals):
    if not intervals:
        return 0

    norm = []
    for a, b in intervals:
        lo = min(a, b)
        hi = max(a, b)
        norm.append((lo, hi))

    norm.sort()
    cur_s, cur_e = norm[0]
    total = 0

    for s, e in norm[1:]:
        if s <= cur_e + 1:
            if e > cur_e:
                cur_e = e
        else:
            total += cur_e - cur_s + 1
            cur_s, cur_e = s, e

    total += cur_e - cur_s + 1
    return total


def strength(pident, qcov, qbp, evalue):
    if (
        pident >= 80.0 and
        qcov >= 80.0 and
        qbp >= 80 and
        evalue <= 1e-10
    ):
        return "STRONG"

    if (
        pident >= 70.0 and
        qcov >= 50.0 and
        qbp >= 80 and
        evalue <= 1e-5
    ):
        return "SUPPORTIVE"

    return "WEAK"


def read_raw(path, database):
    groups = defaultdict(list)

    if not os.path.exists(path):
        raise RuntimeError("Missing raw file: %s" % path)

    raw_n = 0

    with open(path, "r") as fh:
        for line_no, line in enumerate(fh, 1):
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            f = line.split("\t")

            if len(f) < 13:
                raise RuntimeError(
                    "%s line %d has %d columns; expected >=13"
                    % (path, line_no, len(f))
                )

            hit = {
                "qseqid": f[0],
                "sseqid": f[1],
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

            groups[(f[0], f[1])].append(hit)
            raw_n += 1

    rows = []

    for (qseqid, sseqid), hits in groups.items():
        qlen = hits[0]["qlen"]
        slen = hits[0]["slen"]

        q_union_bp = union_length([
            (x["qstart"], x["qend"]) for x in hits
        ])

        s_union_bp = union_length([
            (x["sstart"], x["send"]) for x in hits
        ])

        qcov = 100.0 * q_union_bp / float(qlen) if qlen > 0 else 0.0
        scov = 100.0 * s_union_bp / float(slen) if slen > 0 else 0.0

        qcov = min(qcov, 100.0)
        scov = min(scov, 100.0)

        total_aln = sum(x["aln_len"] for x in hits)

        if total_aln > 0:
            pident = sum(
                x["pident"] * x["aln_len"] for x in hits
            ) / float(total_aln)
        else:
            pident = 0.0

        best_evalue = min(x["evalue"] for x in hits)
        bitscore_sum = sum(x["bitscore"] for x in hits)

        fine = target_class(sseqid)
        broad = broad_class(fine)

        st = strength(
            pident,
            qcov,
            q_union_bp,
            best_evalue
        )

        rows.append({
            "Sample": SAMPLE,
            "Database": database,
            "RepeatModeler_ID": family_id(qseqid),
            "Query_raw": qseqid,
            "Dfam_target": sseqid,
            "Dfam_fine_class": fine,
            "Dfam_broad_class": broad,
            "HSP_count": len(hits),
            "weighted_pident": round(pident, 3),
            "qlen": qlen,
            "query_covered_bp": q_union_bp,
            "qcov": round(qcov, 3),
            "slen": slen,
            "subject_covered_bp": s_union_bp,
            "scov": round(scov, 3),
            "best_evalue": best_evalue,
            "bitscore_sum": round(bitscore_sum, 3),
            "Dfam_strength": st
        })

    return rows, raw_n


xeno_rows, xeno_raw_n = read_raw(
    XENO_RAW,
    "Xenocyprididae"
)

cyprino_rows, cyprino_raw_n = read_raw(
    CYPRINO_RAW,
    "Cypriniformes"
)

all_rows = xeno_rows + cyprino_rows


PAIR_FIELDS = [
    "Sample",
    "Database",
    "RepeatModeler_ID",
    "Query_raw",
    "Dfam_target",
    "Dfam_fine_class",
    "Dfam_broad_class",
    "HSP_count",
    "weighted_pident",
    "qlen",
    "query_covered_bp",
    "qcov",
    "slen",
    "subject_covered_bp",
    "scov",
    "best_evalue",
    "bitscore_sum",
    "Dfam_strength"
]

pair_out = os.path.join(
    OUTDIR,
    "%s.Dfam_v5.pair_summary.tsv" % SAMPLE
)

with open(pair_out, "w") as fo:
    writer = csv.DictWriter(
        fo,
        fieldnames=PAIR_FIELDS,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(all_rows)


rank = {
    "STRONG": 3,
    "SUPPORTIVE": 2,
    "WEAK": 1
}


best_by_db = {}

for r in all_rows:
    key_id = (
        r["RepeatModeler_ID"],
        r["Database"]
    )

    key_rank = (
        rank[r["Dfam_strength"]],
        float(r["qcov"]),
        float(r["weighted_pident"]),
        int(r["query_covered_bp"]),
        float(r["bitscore_sum"])
    )

    if (
        key_id not in best_by_db or
        key_rank > best_by_db[key_id][0]
    ):
        best_by_db[key_id] = (
            key_rank,
            r
        )


families = defaultdict(list)

for r in all_rows:
    families[r["RepeatModeler_ID"]].append(r)


family_rows = []

for fid in sorted(families.keys()):
    rows = families[fid]

    strong_rows = [
        r for r in rows
        if (
            r["Dfam_strength"] == "STRONG" and
            r["Dfam_broad_class"] not in ("", "Unknown")
        )
    ]

    supportive_rows = [
        r for r in rows
        if (
            r["Dfam_strength"] == "SUPPORTIVE" and
            r["Dfam_broad_class"] not in ("", "Unknown")
        )
    ]

    strong_classes = sorted(set(
        r["Dfam_broad_class"]
        for r in strong_rows
    ))

    supportive_classes = sorted(set(
        r["Dfam_broad_class"]
        for r in supportive_rows
    ))

    if len(strong_classes) > 1:
        final_class = "CONFLICT"
        final_strength = "STRONG_CONFLICT"
        conflict = "YES"

    elif len(strong_classes) == 1:
        final_class = strong_classes[0]

        opposite_support = [
            x for x in supportive_classes
            if x != final_class
        ]

        if opposite_support:
            final_strength = "STRONG_WITH_WEAK_CONFLICT"
            conflict = "WEAK"
        else:
            final_strength = "STRONG"
            conflict = "NO"

    elif len(supportive_classes) > 1:
        final_class = "CONFLICT"
        final_strength = "SUPPORTIVE_CONFLICT"
        conflict = "YES"

    elif len(supportive_classes) == 1:
        final_class = supportive_classes[0]
        final_strength = "SUPPORTIVE"
        conflict = "NO"

    else:
        final_class = ""
        final_strength = "NONE"
        conflict = "NO"

    xeno = best_by_db.get(
        (fid, "Xenocyprididae")
    )

    cypr = best_by_db.get(
        (fid, "Cypriniformes")
    )

    xr = xeno[1] if xeno else {}
    cr = cypr[1] if cypr else {}

    family_rows.append({
        "Sample": SAMPLE,
        "RepeatModeler_ID": fid,

        "Xeno_best_broad":
            xr.get("Dfam_broad_class", ""),
        "Xeno_best_fine":
            xr.get("Dfam_fine_class", ""),
        "Xeno_best_strength":
            xr.get("Dfam_strength", ""),
        "Xeno_best_pident":
            xr.get("weighted_pident", ""),
        "Xeno_best_qcov":
            xr.get("qcov", ""),
        "Xeno_best_scov":
            xr.get("scov", ""),
        "Xeno_best_evalue":
            xr.get("best_evalue", ""),

        "Cyprino_best_broad":
            cr.get("Dfam_broad_class", ""),
        "Cyprino_best_fine":
            cr.get("Dfam_fine_class", ""),
        "Cyprino_best_strength":
            cr.get("Dfam_strength", ""),
        "Cyprino_best_pident":
            cr.get("weighted_pident", ""),
        "Cyprino_best_qcov":
            cr.get("qcov", ""),
        "Cyprino_best_scov":
            cr.get("scov", ""),
        "Cyprino_best_evalue":
            cr.get("best_evalue", ""),

        "Strong_broad_classes":
            ",".join(strong_classes),

        "Supportive_broad_classes":
            ",".join(supportive_classes),

        "Strong_hit_count":
            len(strong_rows),

        "Supportive_hit_count":
            len(supportive_rows),

        "Dfam_v5_class":
            final_class,

        "Dfam_v5_strength":
            final_strength,

        "Dfam_v5_conflict":
            conflict
    })


FAMILY_FIELDS = [
    "Sample",
    "RepeatModeler_ID",

    "Xeno_best_broad",
    "Xeno_best_fine",
    "Xeno_best_strength",
    "Xeno_best_pident",
    "Xeno_best_qcov",
    "Xeno_best_scov",
    "Xeno_best_evalue",

    "Cyprino_best_broad",
    "Cyprino_best_fine",
    "Cyprino_best_strength",
    "Cyprino_best_pident",
    "Cyprino_best_qcov",
    "Cyprino_best_scov",
    "Cyprino_best_evalue",

    "Strong_broad_classes",
    "Supportive_broad_classes",

    "Strong_hit_count",
    "Supportive_hit_count",

    "Dfam_v5_class",
    "Dfam_v5_strength",
    "Dfam_v5_conflict"
]


family_out = os.path.join(
    OUTDIR,
    "%s.Dfam_v5.family_summary.tsv" % SAMPLE
)

with open(family_out, "w") as fo:
    writer = csv.DictWriter(
        fo,
        fieldnames=FAMILY_FIELDS,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(family_rows)


counts = defaultdict(int)

for r in family_rows:
    counts[r["Dfam_v5_strength"]] += 1


qc_out = os.path.join(
    OUTDIR,
    "%s.Dfam_v5.summary.tsv" % SAMPLE
)

with open(qc_out, "w") as fo:
    fo.write(
        "Sample\tXeno_raw_HSPs\tCyprino_raw_HSPs\t"
        "Families_with_any_Dfam_hit\t"
        "STRONG\tSTRONG_WITH_WEAK_CONFLICT\t"
        "STRONG_CONFLICT\tSUPPORTIVE\t"
        "SUPPORTIVE_CONFLICT\n"
    )

    fo.write(
        "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            SAMPLE,
            xeno_raw_n,
            cyprino_raw_n,
            len(family_rows),
            counts["STRONG"],
            counts["STRONG_WITH_WEAK_CONFLICT"],
            counts["STRONG_CONFLICT"],
            counts["SUPPORTIVE"],
            counts["SUPPORTIVE_CONFLICT"]
        )
    )


print("Sample:", SAMPLE)
print(
    "Dfam families with any hit:",
    len(family_rows)
)
print(
    "STRONG:",
    counts["STRONG"]
)
print(
    "STRONG_WITH_WEAK_CONFLICT:",
    counts["STRONG_WITH_WEAK_CONFLICT"]
)
print(
    "STRONG_CONFLICT:",
    counts["STRONG_CONFLICT"]
)
print(
    "SUPPORTIVE:",
    counts["SUPPORTIVE"]
)
print(
    "SUPPORTIVE_CONFLICT:",
    counts["SUPPORTIVE_CONFLICT"]
)
print(
    "Finished:",
    family_out
)