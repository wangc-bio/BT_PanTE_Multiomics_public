#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

import argparse
import os
import csv
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize RepBase BLAST evidence for one consensus library."
    )
    parser.add_argument("--sample-id", required=True,
                        help="Generic identifier for this consensus library.")
    parser.add_argument("--raw-tsv", required=True,
                        help="RepBase BLAST tabular output for this library.")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for pair-level and family-level summaries.")
    return parser.parse_args()


ARGS = parse_args()
sample = ARGS.sample_id
work = os.path.abspath(ARGS.output_dir)
raw = os.path.abspath(ARGS.raw_tsv)
os.makedirs(work, exist_ok=True)

if not os.path.isfile(raw):
    raise FileNotFoundError("Missing RepBase BLAST table: %s" % raw)


def family_id(x):
    return x.strip().lstrip(">").split("#", 1)[0]


def parse_subject(subject):
    parts = subject.split("|")

    name = parts[0]
    rbclass = ""
    species = ""

    for p in parts[1:]:

        if p.startswith("RBclass="):
            rbclass = p.split("=", 1)[1]

        elif p.startswith("Species="):
            species = p.split("=", 1)[1]

    return name, rbclass, species


def broad_class(rbclass, subject_name):

    rc = (rbclass or "").strip().lower()
    nm = (subject_name or "").strip().lower()

    text = rc + " " + nm

    # ==================================================
    # RepBase explicit v5 broad-class map
    # Conservative: only taxonomically unambiguous calls
    # ==================================================

    explicit_map = {
        # DNA transposons
        "dada": "DNA",
        "p": "DNA",
        "sola1": "DNA",
        "sola2": "DNA",
        "isl2eu": "DNA",

        # LTR retrotransposons
        "troyka": "LTR",

        # Penelope / non-LTR broad bin
        "neptune": "LINE",
        "nimb": "LINE",
        "crack": "LINE",
        "daphne": "LINE",
        "hero": "LINE",
        "proto2": "LINE",
        "r2": "LINE",

        # explicit non-TE repeats
        "sat": "nonTE",
        "msat": "nonTE",
        "satellite": "nonTE",
        "trna": "nonTE",
        "rp5s": "nonTE",
        "multicopy_gene": "nonTE"
    }

    if rc in explicit_map:
        return explicit_map[rc]

    if not rc and not nm:
        return "Unknown"

    # --------------------------------------------------
    # SINE
    # --------------------------------------------------

    if (
        "sine" in text or
        rc.startswith("sine") or
        rc in ("alu", "mir")
    ):
        return "SINE"

    # --------------------------------------------------
    # Rolling-circle
    # --------------------------------------------------

    if (
        "helitron" in text or
        rc in (
            "helitron",
            "rolling_circle",
            "rolling-circle"
        )
    ):
        return "RC"

    # --------------------------------------------------
    # LTR retrotransposons
    # --------------------------------------------------

    ltr_tokens = (
        "gypsy",
        "copia",
        "bel",
        "pao",
        "erv",
        "retrovirus",
        "ltr",
        "dirs"
    )

    if any(
        token in text
        for token in ltr_tokens
    ):
        return "LTR"

    # --------------------------------------------------
    # LINE / non-LTR
    # --------------------------------------------------

    line_tokens = (
        "line",
        "non-ltr",
        "nonltr",
        "cr1",
        "rex",
        "rte",
        "jockey",
        "tx1",
        "penelope",
        "ingi"
    )

    if any(
        token in text
        for token in line_tokens
    ):
        return "LINE"

    if rc in (
        "l1",
        "l2"
    ):
        return "LINE"

    # --------------------------------------------------
    # DNA transposons
    # --------------------------------------------------

    dna_tokens = (
        "dna_transposon",
        "dna transposon",
        "hat",
        "harbinger",
        "mariner",
        "tc1",
        "piggybac",
        "mutator",
        "mudr",
        "cacta",
        "enspm",
        "kolobok",
        "p-element",
        "p element",
        "maverick",
        "polinton",
        "transib",
        "zisupton",
        "academ",
        "ginger",
        "is3eu",
        "piggybat",
        "merlin",
        "crypton"
    )

    if any(
        token in text
        for token in dna_tokens
    ):
        return "DNA"

    return "Unknown"


def union_length(intervals):

    if not intervals:
        return 0

    norm = []

    for a, b in intervals:

        lo = min(a, b)
        hi = max(a, b)

        norm.append(
            (lo, hi)
        )

    norm.sort()

    cur_s, cur_e = norm[0]

    total = 0

    for s, e in norm[1:]:

        if s <= cur_e + 1:

            if e > cur_e:
                cur_e = e

        else:

            total += (
                cur_e - cur_s + 1
            )

            cur_s, cur_e = s, e

    total += (
        cur_e - cur_s + 1
    )

    return total


def strength(
    pident,
    qcov,
    qbp,
    evalue
):

    # v5 STRONG
    if (
        pident >= 80.0 and
        qcov >= 80.0 and
        qbp >= 80 and
        evalue <= 1e-10
    ):
        return "STRONG"

    # v5 SUPPORTIVE
    if (
        pident >= 70.0 and
        qcov >= 50.0 and
        qbp >= 80 and
        evalue <= 1e-5
    ):
        return "SUPPORTIVE"

    return "WEAK"


def read_raw(path, sample):

    groups = defaultdict(list)

    raw_n = 0

    with open(path, "r") as fh:

        for line_no, line in enumerate(
            fh,
            1
        ):

            line = line.rstrip("\n")

            if (
                not line or
                line.startswith("#")
            ):
                continue

            f = line.split("\t")

            if len(f) < 13:

                raise RuntimeError(
                    "%s line %d has %d columns; expected >=13"
                    % (
                        path,
                        line_no,
                        len(f)
                    )
                )

            hit = {
                "qseqid":
                    f[0],

                "sseqid":
                    f[1],

                "pident":
                    float(f[2]),

                "aln_len":
                    int(float(f[3])),

                "qstart":
                    int(float(f[4])),

                "qend":
                    int(float(f[5])),

                "qlen":
                    int(float(f[6])),

                "sstart":
                    int(float(f[7])),

                "send":
                    int(float(f[8])),

                "slen":
                    int(float(f[9])),

                "evalue":
                    float(f[10]),

                "bitscore":
                    float(f[11])
            }

            groups[
                (
                    f[0],
                    f[1]
                )
            ].append(hit)

            raw_n += 1

    pair_rows = []

    for (
        qseqid,
        sseqid
    ), hits in groups.items():

        qlen = hits[0]["qlen"]
        slen = hits[0]["slen"]

        q_union_bp = union_length([
            (
                x["qstart"],
                x["qend"]
            )
            for x in hits
        ])

        s_union_bp = union_length([
            (
                x["sstart"],
                x["send"]
            )
            for x in hits
        ])

        if qlen > 0:
            qcov = (
                100.0 *
                q_union_bp /
                float(qlen)
            )
        else:
            qcov = 0.0

        if slen > 0:
            scov = (
                100.0 *
                s_union_bp /
                float(slen)
            )
        else:
            scov = 0.0

        qcov = min(
            qcov,
            100.0
        )

        scov = min(
            scov,
            100.0
        )

        total_aln = sum(
            x["aln_len"]
            for x in hits
        )

        if total_aln > 0:

            pident = sum(
                x["pident"] *
                x["aln_len"]
                for x in hits
            ) / float(total_aln)

        else:
            pident = 0.0

        best_evalue = min(
            x["evalue"]
            for x in hits
        )

        bitscore_sum = sum(
            x["bitscore"]
            for x in hits
        )

        (
            subject_name,
            rbclass,
            species
        ) = parse_subject(
            sseqid
        )

        broad = broad_class(
            rbclass,
            subject_name
        )

        st = strength(
            pident,
            qcov,
            q_union_bp,
            best_evalue
        )

        pair_rows.append({

            "Sample":
                sample,

            "RepeatModeler_ID":
                family_id(qseqid),

            "Query_raw":
                qseqid,

            "RepBase_subject_raw":
                sseqid,

            "RepBase_subject":
                subject_name,

            "RepBase_RBclass":
                rbclass,

            "RepBase_species":
                species,

            "RepBase_broad":
                broad,

            "HSP_count":
                len(hits),

            "weighted_pident":
                round(
                    pident,
                    3
                ),

            "qlen":
                qlen,

            "query_covered_bp":
                q_union_bp,

            "qcov":
                round(
                    qcov,
                    3
                ),

            "slen":
                slen,

            "subject_covered_bp":
                s_union_bp,

            "scov":
                round(
                    scov,
                    3
                ),

            "best_evalue":
                best_evalue,

            "bitscore_sum":
                round(
                    bitscore_sum,
                    3
                ),

            "RepBase_strength":
                st
        })

    return (
        pair_rows,
        raw_n
    )


PAIR_FIELDS = [
    "Sample",
    "RepeatModeler_ID",
    "Query_raw",
    "RepBase_subject_raw",
    "RepBase_subject",
    "RepBase_RBclass",
    "RepBase_species",
    "RepBase_broad",
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
    "RepBase_strength"
]


FAMILY_FIELDS = [
    "Sample",
    "RepeatModeler_ID",

    "Best_subject",
    "Best_RBclass",
    "Best_species",
    "Best_broad",
    "Best_strength",

    "Best_pident",
    "Best_qcov",
    "Best_scov",
    "Best_query_covered_bp",
    "Best_evalue",
    "Best_bitscore_sum",

    "Strong_broad_classes",
    "Supportive_broad_classes",

    "Strong_hit_count",
    "Supportive_hit_count",

    "RepBase_v5_class",
    "RepBase_v5_strength",
    "RepBase_v5_conflict"
]


RANK = {
    "STRONG": 3,
    "SUPPORTIVE": 2,
    "WEAK": 1
}


summary_rows = []

unknown_rbclasses_global = (
    defaultdict(int)
)


pair_rows, raw_n = read_raw(
    raw,
    sample
)

pair_out = os.path.join(
    work,
    "%s.RepBase_v5.pair_summary.tsv"
    % sample
)

with open(
    pair_out,
    "w"
) as fo:

    writer = csv.DictWriter(
        fo,
        fieldnames=PAIR_FIELDS,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()

    writer.writerows(
        pair_rows
    )

by_family = defaultdict(list)

for r in pair_rows:

    by_family[
        r["RepeatModeler_ID"]
    ].append(r)

    if (
        r["RepBase_broad"]
        == "Unknown"
    ):

        key = (
            r["RepBase_RBclass"]
            or
            "EMPTY_RBCLASS"
        )

        unknown_rbclasses_global[
            key
        ] += 1

family_rows = []

for fid in sorted(
    by_family.keys()
):

    rows = by_family[fid]

    best = max(
        rows,
        key=lambda r: (
            RANK[
                r[
                    "RepBase_strength"
                ]
            ],
            float(r["qcov"]),
            float(
                r[
                    "weighted_pident"
                ]
            ),
            int(
                r[
                    "query_covered_bp"
                ]
            ),
            float(
                r[
                    "bitscore_sum"
                ]
            )
        )
    )

    strong_rows = [
        r
        for r in rows
        if (
            r[
                "RepBase_strength"
            ] == "STRONG"
            and
            r[
                "RepBase_broad"
            ] != "Unknown"
        )
    ]

    supportive_rows = [
        r
        for r in rows
        if (
            r[
                "RepBase_strength"
            ] == "SUPPORTIVE"
            and
            r[
                "RepBase_broad"
            ] != "Unknown"
        )
    ]

    strong_classes = sorted(
        set(
            r[
                "RepBase_broad"
            ]
            for r
            in strong_rows
        )
    )

    supportive_classes = sorted(
        set(
            r[
                "RepBase_broad"
            ]
            for r
            in supportive_rows
        )
    )

    # --------------------------------------------------
    # family-level classification
    # --------------------------------------------------

    if len(
        strong_classes
    ) > 1:

        final_class = (
            "CONFLICT"
        )

        final_strength = (
            "STRONG_CONFLICT"
        )

        conflict = "YES"

    elif len(
        strong_classes
    ) == 1:

        final_class = (
            strong_classes[0]
        )

        opposite_support = [
            x
            for x
            in supportive_classes
            if x != final_class
        ]

        if opposite_support:

            final_strength = (
                "STRONG_WITH_WEAK_CONFLICT"
            )

            conflict = "WEAK"

        else:

            final_strength = (
                "STRONG"
            )

            conflict = "NO"

    elif len(
        supportive_classes
    ) > 1:

        final_class = (
            "CONFLICT"
        )

        final_strength = (
            "SUPPORTIVE_CONFLICT"
        )

        conflict = "YES"

    elif len(
        supportive_classes
    ) == 1:

        final_class = (
            supportive_classes[0]
        )

        final_strength = (
            "SUPPORTIVE"
        )

        conflict = "NO"

    else:

        final_class = ""

        final_strength = (
            "NONE"
        )

        conflict = "NO"

    family_rows.append({

        "Sample":
            sample,

        "RepeatModeler_ID":
            fid,

        "Best_subject":
            best[
                "RepBase_subject"
            ],

        "Best_RBclass":
            best[
                "RepBase_RBclass"
            ],

        "Best_species":
            best[
                "RepBase_species"
            ],

        "Best_broad":
            best[
                "RepBase_broad"
            ],

        "Best_strength":
            best[
                "RepBase_strength"
            ],

        "Best_pident":
            best[
                "weighted_pident"
            ],

        "Best_qcov":
            best[
                "qcov"
            ],

        "Best_scov":
            best[
                "scov"
            ],

        "Best_query_covered_bp":
            best[
                "query_covered_bp"
            ],

        "Best_evalue":
            best[
                "best_evalue"
            ],

        "Best_bitscore_sum":
            best[
                "bitscore_sum"
            ],

        "Strong_broad_classes":
            ",".join(
                strong_classes
            ),

        "Supportive_broad_classes":
            ",".join(
                supportive_classes
            ),

        "Strong_hit_count":
            len(
                strong_rows
            ),

        "Supportive_hit_count":
            len(
                supportive_rows
            ),

        "RepBase_v5_class":
            final_class,

        "RepBase_v5_strength":
            final_strength,

        "RepBase_v5_conflict":
            conflict
    })

family_out = os.path.join(
    work,
    "%s.RepBase_v5.family_summary.tsv"
    % sample
)

with open(
    family_out,
    "w"
) as fo:

    writer = csv.DictWriter(
        fo,
        fieldnames=FAMILY_FIELDS,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()

    writer.writerows(
        family_rows
    )

strength_counts = (
    defaultdict(int)
)

class_counts = (
    defaultdict(int)
)

for r in family_rows:

    strength_counts[
        r[
            "RepBase_v5_strength"
        ]
    ] += 1

    if r[
        "RepBase_v5_strength"
    ] in (
        "STRONG",
        "STRONG_WITH_WEAK_CONFLICT",
        "SUPPORTIVE"
    ):

        class_counts[
            r[
                "RepBase_v5_class"
            ]
        ] += 1

summary_rows.append({

    "Sample":
        sample,

    "Raw_HSPs":
        raw_n,

    "Families_with_any_hit":
        len(
            family_rows
        ),

    "STRONG":
        strength_counts[
            "STRONG"
        ],

    "STRONG_WITH_WEAK_CONFLICT":
        strength_counts[
            "STRONG_WITH_WEAK_CONFLICT"
        ],

    "STRONG_CONFLICT":
        strength_counts[
            "STRONG_CONFLICT"
        ],

    "SUPPORTIVE":
        strength_counts[
            "SUPPORTIVE"
        ],

    "SUPPORTIVE_CONFLICT":
        strength_counts[
            "SUPPORTIVE_CONFLICT"
        ],

    "NONE":
        strength_counts[
            "NONE"
        ],

    "DNA":
        class_counts[
            "DNA"
        ],

    "LTR":
        class_counts[
            "LTR"
        ],

    "LINE":
        class_counts[
            "LINE"
        ],

    "SINE":
        class_counts[
            "SINE"
        ],

    "RC":
        class_counts[
            "RC"
        ]
})

print(
    "%s hit=%d strong=%d "
    "strong_conflict=%d supportive=%d"
    % (
        sample,

        len(
            family_rows
        ),

        strength_counts[
            "STRONG"
        ]
        +
        strength_counts[
            "STRONG_WITH_WEAK_CONFLICT"
        ],

        strength_counts[
            "STRONG_CONFLICT"
        ],

        strength_counts[
            "SUPPORTIVE"
        ]
    )
)


summary_file = os.path.join(
    work,
    "%s.RepBase_v5.summary.tsv" % sample
)


summary_fields = [
    "Sample",
    "Raw_HSPs",
    "Families_with_any_hit",
    "STRONG",
    "STRONG_WITH_WEAK_CONFLICT",
    "STRONG_CONFLICT",
    "SUPPORTIVE",
    "SUPPORTIVE_CONFLICT",
    "NONE",
    "DNA",
    "LTR",
    "LINE",
    "SINE",
    "RC"
]


with open(
    summary_file,
    "w"
) as fo:

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


unknown_file = os.path.join(
    work,
    "RepBase_v5.unmapped_RBclass.tsv"
)


with open(
    unknown_file,
    "w"
) as fo:

    fo.write(
        "RBclass\tPair_count\n"
    )

    for rbclass, n in sorted(
        unknown_rbclasses_global.items(),
        key=lambda x: (
            -x[1],
            x[0]
        )
    ):

        fo.write(
            "%s\t%d\n"
            % (
                rbclass,
                n
            )
        )


print("")
print("Finished:")
print(summary_file)
print(unknown_file)