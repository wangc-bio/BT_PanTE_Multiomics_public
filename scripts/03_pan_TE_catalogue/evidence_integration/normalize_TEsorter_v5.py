#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

import argparse
import os
import csv
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize TEsorter classification and domain evidence for one consensus library."
    )
    parser.add_argument("--sample-id", required=True,
                        help="Generic identifier for this consensus library.")
    parser.add_argument("--cls-tsv", required=True,
                        help="TEsorter .cls.tsv output for this library.")
    parser.add_argument("--dom-tsv", required=True,
                        help="TEsorter .dom.tsv output for this library.")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for family-level and QC summaries.")
    return parser.parse_args()


ARGS = parse_args()
sample = ARGS.sample_id
work = os.path.abspath(ARGS.output_dir)
cls_file = os.path.abspath(ARGS.cls_tsv)
dom_file = os.path.abspath(ARGS.dom_tsv)
os.makedirs(work, exist_ok=True)

for required_file in (cls_file, dom_file):
    if not os.path.isfile(required_file):
        raise FileNotFoundError("Missing TEsorter input: %s" % required_file)


def family_id(x):
    x = x.strip().lstrip(">")
    x = x.split("|", 1)[0]
    return x.split("#", 1)[0]


def broad_class(text):
    s = (text or "").strip().lower()

    if not s or s in ("unknown", "na", "none"):
        return "Unknown"

    if "mixture" in s or "mixed" in s:
        return "Unknown"

    if "maverick" in s or "polinton" in s:
        return "DNA"

    if "helitron" in s:
        return "RC"

    if "sine" in s:
        return "SINE"

    if (
        "penelope" in s or
        s == "ple" or
        "/ple" in s
    ):
        return "LINE"

    if (
        "line" in s or
        "non-ltr" in s or
        "nonltr" in s
    ):
        return "LINE"

    if (
        "ltr" in s or
        "gypsy" in s or
        "copia" in s or
        "bel" in s or
        "pao" in s or
        "erv" in s or
        "dirs" in s
    ):
        return "LTR"

    if (
        "class_ii" in s or
        "class ii" in s or
        "dna" in s or
        "tir" in s or
        "transposase" in s
    ):
        return "DNA"

    return "Unknown"


def domain_strength(evalue, coverage):

    if (
        evalue <= 1e-5 and
        coverage >= 30.0
    ):
        return "STRONG"

    if (
        evalue <= 1e-3 and
        coverage >= 20.0
    ):
        return "SUPPORTIVE"

    return "WEAK"


def read_cls(path):
    data = {}

    with open(path, "r") as fh:
        reader = csv.reader(
            fh,
            delimiter="\t"
        )

        for row in reader:

            if not row:
                continue

            if row[0].startswith("#"):
                continue

            if len(row) < 7:
                continue

            fid = family_id(row[0])

            order = row[1].strip()
            superfamily = row[2].strip()
            clade = row[3].strip()
            complete = row[4].strip()
            strand = row[5].strip()
            domains = row[6].strip()

            cls_text = order

            if (
                superfamily and
                superfamily.lower() != "unknown"
            ):
                cls_text += "/" + superfamily

            data[fid] = {
                "TEsorter_order":
                    order,

                "TEsorter_superfamily":
                    superfamily,

                "TEsorter_clade":
                    clade,

                "TEsorter_complete":
                    complete,

                "TEsorter_strand":
                    strand,

                "TEsorter_domains_from_cls":
                    domains,

                "TEsorter_cls_broad":
                    broad_class(cls_text)
            }

    return data


def read_dom(path):
    hits = defaultdict(list)

    with open(path, "r") as fh:
        reader = csv.reader(
            fh,
            delimiter="\t"
        )

        for row in reader:

            if not row:
                continue

            if row[0].startswith("#"):
                continue

            if len(row) < 6:
                continue

            raw_id = row[0].strip()
            fid = family_id(raw_id)

            try:
                length = int(float(row[1]))
                evalue = float(row[2])
                coverage = float(row[3])
                probability = float(row[4])
                score = float(row[5])

            except ValueError:
                continue

            annotation = ""

            if "|" in raw_id:
                annotation = raw_id.split(
                    "|",
                    1
                )[1]

            hits[fid].append({
                "raw_id":
                    raw_id,

                "annotation":
                    annotation,

                "length":
                    length,

                "evalue":
                    evalue,

                "coverage":
                    coverage,

                "probability":
                    probability,

                "score":
                    score,

                "broad":
                    broad_class(annotation),

                "strength":
                    domain_strength(
                        evalue,
                        coverage
                    )
            })

    return hits


def best_domain(rows):

    rank = {
        "STRONG": 3,
        "SUPPORTIVE": 2,
        "WEAK": 1
    }

    return max(
        rows,
        key=lambda x: (
            rank[x["strength"]],
            x["coverage"],
            -x["evalue"],
            x["length"],
            x["score"]
        )
    )


summary_rows = []


cls = read_cls(
    cls_file
)

dom = read_dom(
    dom_file
)

all_ids = sorted(
    set(cls.keys()) |
    set(dom.keys())
)

out_rows = []

for fid in all_ids:

    c = cls.get(
        fid,
        {}
    )

    dh = dom.get(
        fid,
        []
    )

    cls_broad = c.get(
        "TEsorter_cls_broad",
        "Unknown"
    )

    strong_classes = sorted(
        set(
            x["broad"]
            for x in dh
            if (
                x["strength"] == "STRONG" and
                x["broad"] != "Unknown"
            )
        )
    )

    support_classes = sorted(
        set(
            x["broad"]
            for x in dh
            if (
                x["strength"] == "SUPPORTIVE" and
                x["broad"] != "Unknown"
            )
        )
    )

    strongest = (
        best_domain(dh)
        if dh else None
    )

    best_strength = (
        strongest["strength"]
        if strongest
        else "NONE"
    )

    final_class = "Unknown"
    final_strength = "NONE"
    conflict = "NO"
    reason = (
        "No usable domain evidence"
    )

    if len(strong_classes) > 1:

        final_class = "CONFLICT"
        final_strength = (
            "STRONG_CONFLICT"
        )
        conflict = "YES"
        reason = (
            "Strong domains support "
            "multiple broad classes"
        )

    elif len(strong_classes) == 1:

        domain_class = (
            strong_classes[0]
        )

        if (
            cls_broad != "Unknown" and
            cls_broad != domain_class
        ):

            final_class = "CONFLICT"
            final_strength = (
                "STRONG_CONFLICT"
            )
            conflict = "YES"

            reason = (
                "TEsorter cls broad "
                "conflicts with strong "
                "domain class"
            )

        else:

            final_class = domain_class
            final_strength = "STRONG"

            reason = (
                "Strong domain evidence"
            )

            if cls_broad != "Unknown":

                final_class = cls_broad

                reason = (
                    "Strong domain evidence "
                    "consistent with "
                    "TEsorter cls"
                )

            if any(
                x != final_class
                for x in support_classes
            ):

                final_strength = (
                    "STRONG_WITH_WEAK_CONFLICT"
                )

                conflict = "WEAK"

                reason += (
                    "; supportive domain "
                    "conflict present"
                )

    elif len(support_classes) > 1:

        final_class = "CONFLICT"

        final_strength = (
            "SUPPORTIVE_CONFLICT"
        )

        conflict = "YES"

        reason = (
            "Supportive domains support "
            "multiple broad classes"
        )

    elif len(support_classes) == 1:

        domain_class = (
            support_classes[0]
        )

        if (
            cls_broad != "Unknown" and
            cls_broad != domain_class
        ):

            final_class = "CONFLICT"

            final_strength = (
                "SUPPORTIVE_CONFLICT"
            )

            conflict = "YES"

            reason = (
                "TEsorter cls broad "
                "conflicts with supportive "
                "domain class"
            )

        else:

            final_class = (
                cls_broad
                if cls_broad != "Unknown"
                else domain_class
            )

            final_strength = (
                "SUPPORTIVE"
            )

            reason = (
                "Supportive domain evidence"
            )

    elif cls_broad != "Unknown":

        final_class = cls_broad

        final_strength = (
            "CLS_ONLY_WEAK"
        )

        reason = (
            "TEsorter cls exists but no "
            "domain passed v5 threshold"
        )

    out_rows.append({
        "Sample":
            sample,

        "RepeatModeler_ID":
            fid,

        "TEsorter_order":
            c.get(
                "TEsorter_order",
                ""
            ),

        "TEsorter_superfamily":
            c.get(
                "TEsorter_superfamily",
                ""
            ),

        "TEsorter_clade":
            c.get(
                "TEsorter_clade",
                ""
            ),

        "TEsorter_complete":
            c.get(
                "TEsorter_complete",
                ""
            ),

        "TEsorter_strand":
            c.get(
                "TEsorter_strand",
                ""
            ),

        "TEsorter_domains_from_cls":
            c.get(
                "TEsorter_domains_from_cls",
                ""
            ),

        "TEsorter_cls_broad":
            cls_broad,

        "Domain_record_count":
            len(dh),

        "Strong_domain_classes":
            ",".join(
                strong_classes
            ),

        "Supportive_domain_classes":
            ",".join(
                support_classes
            ),

        "Best_domain_annotation":
            (
                strongest["annotation"]
                if strongest
                else ""
            ),

        "Best_domain_length":
            (
                strongest["length"]
                if strongest
                else ""
            ),

        "Best_domain_evalue":
            (
                strongest["evalue"]
                if strongest
                else ""
            ),

        "Best_domain_coverage":
            (
                strongest["coverage"]
                if strongest
                else ""
            ),

        "Best_domain_probability":
            (
                strongest["probability"]
                if strongest
                else ""
            ),

        "Best_domain_score":
            (
                strongest["score"]
                if strongest
                else ""
            ),

        "Best_domain_strength":
            best_strength,

        "TEsorter_v5_class":
            final_class,

        "TEsorter_v5_strength":
            final_strength,

        "TEsorter_v5_conflict":
            conflict,

        "TEsorter_v5_reason":
            reason
    })

fields = [
    "Sample",
    "RepeatModeler_ID",
    "TEsorter_order",
    "TEsorter_superfamily",
    "TEsorter_clade",
    "TEsorter_complete",
    "TEsorter_strand",
    "TEsorter_domains_from_cls",
    "TEsorter_cls_broad",
    "Domain_record_count",
    "Strong_domain_classes",
    "Supportive_domain_classes",
    "Best_domain_annotation",
    "Best_domain_length",
    "Best_domain_evalue",
    "Best_domain_coverage",
    "Best_domain_probability",
    "Best_domain_score",
    "Best_domain_strength",
    "TEsorter_v5_class",
    "TEsorter_v5_strength",
    "TEsorter_v5_conflict",
    "TEsorter_v5_reason"
]

outfile = os.path.join(
    work,
    "%s.TEsorter_v5.family_summary.tsv"
    % sample
)

with open(outfile, "w") as fo:

    writer = csv.DictWriter(
        fo,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        out_rows
    )

counts = defaultdict(int)
class_counts = defaultdict(int)

for r in out_rows:

    counts[
        r["TEsorter_v5_strength"]
    ] += 1

    if r["TEsorter_v5_strength"] in (
        "STRONG",
        "STRONG_WITH_WEAK_CONFLICT",
        "SUPPORTIVE"
    ):

        class_counts[
            r["TEsorter_v5_class"]
        ] += 1

summary_rows.append({
    "Sample":
        sample,

    "Families_in_cls":
        len(cls),

    "Families_with_domains":
        len(dom),

    "STRONG":
        counts["STRONG"],

    "STRONG_WITH_WEAK_CONFLICT":
        counts[
            "STRONG_WITH_WEAK_CONFLICT"
        ],

    "STRONG_CONFLICT":
        counts["STRONG_CONFLICT"],

    "SUPPORTIVE":
        counts["SUPPORTIVE"],

    "SUPPORTIVE_CONFLICT":
        counts[
            "SUPPORTIVE_CONFLICT"
        ],

    "CLS_ONLY_WEAK":
        counts["CLS_ONLY_WEAK"],

    "DNA":
        class_counts["DNA"],

    "LTR":
        class_counts["LTR"],

    "LINE":
        class_counts["LINE"],

    "SINE":
        class_counts["SINE"],

    "RC":
        class_counts["RC"]
})

print(
    "%s cls=%d dom_families=%d "
    "strong=%d strong_conflict=%d "
    "supportive=%d"
    % (
        sample,
        len(cls),
        len(dom),

        counts["STRONG"] +
        counts[
            "STRONG_WITH_WEAK_CONFLICT"
        ],

        counts[
            "STRONG_CONFLICT"
        ],

        counts[
            "SUPPORTIVE"
        ]
    )
)


summary_file = os.path.join(
    work,
    "%s.TEsorter_v5.summary.tsv" % sample
)

summary_fields = [
    "Sample",
    "Families_in_cls",
    "Families_with_domains",
    "STRONG",
    "STRONG_WITH_WEAK_CONFLICT",
    "STRONG_CONFLICT",
    "SUPPORTIVE",
    "SUPPORTIVE_CONFLICT",
    "CLS_ONLY_WEAK",
    "DNA",
    "LTR",
    "LINE",
    "SINE",
    "RC"
]

with open(summary_file, "w") as fo:

    writer = csv.DictWriter(
        fo,
        fieldnames=summary_fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        summary_rows
    )

print("")
print("Finished:")
print(summary_file)