#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Generalized from the first complete program in the supplied concatenated source file.

from __future__ import print_function

import argparse
import os
import csv
from collections import defaultdict

REQUIRED_MANIFEST_COLUMNS = {
    "sample_id", "repeatmodeler_fasta", "prior_evidence",
    "dfam_family_summary", "repbase_family_summary",
    "tesorter_family_summary", "ltr_reciprocal_family_best",
    "mite_hc_best", "mite_supportive_best",
    "helitron_hc_best", "helitron_supportive_best",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Integrate homology, domain and structural TE evidence for each consensus library.")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--review-dir", required=True)
    return parser.parse_args()


def read_manifest(path):
    base = os.path.dirname(os.path.abspath(path))
    entries, seen = [], set()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not REQUIRED_MANIFEST_COLUMNS.issubset(reader.fieldnames):
            missing = sorted(REQUIRED_MANIFEST_COLUMNS - set(reader.fieldnames or []))
            raise ValueError("Manifest is missing required columns: %s" % ", ".join(missing))
        for line_no, row in enumerate(reader, 2):
            sample = (row.get("sample_id") or "").strip()
            if not sample:
                raise ValueError("Empty sample_id at manifest line %d" % line_no)
            if sample in seen:
                raise ValueError("Duplicate sample_id: %s" % sample)
            seen.add(sample)
            item = {"sample_id": sample}
            for key in REQUIRED_MANIFEST_COLUMNS - {"sample_id"}:
                value = (row.get(key) or "").strip()
                if not value:
                    raise ValueError("Missing %s for %s" % (key, sample))
                if not os.path.isabs(value):
                    value = os.path.join(base, value)
                value = os.path.abspath(value)
                if not os.path.isfile(value):
                    raise FileNotFoundError("Missing %s for %s: %s" % (key, sample, value))
                item[key] = value
            entries.append(item)
    if not entries:
        raise ValueError("The manifest contains no genome/library entries.")
    return entries


ARGS = parse_args()
ENTRIES = read_manifest(ARGS.manifest)
OUT = os.path.abspath(ARGS.output_dir)
REVIEW = os.path.abspath(ARGS.review_dir)
os.makedirs(OUT, exist_ok=True)
os.makedirs(REVIEW, exist_ok=True)


def fam_id(x):
    return (
        x.strip()
        .lstrip(">")
        .split()[0]
        .split("#", 1)[0]
    )


def broad_from_fine(x):
    s = (x or "").lower()

    if (
        "simple_repeat" in s or
        "satellite" in s or
        "low_complexity" in s or
        s in ("trna", "rrna")
    ):
        return "nonTE"

    if "sine" in s:
        return "SINE"

    if (
        "helitron" in s or
        s.startswith("rc/")
    ):
        return "RC"

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
        "line" in s or
        "non-ltr" in s or
        "nonltr" in s or
        "cr1" in s or
        "rex" in s or
        "rte" in s or
        "jockey" in s or
        "penelope" in s
    ):
        return "LINE"

    if (
        "dna" in s or
        "tir" in s or
        "mite" in s or
        "maverick" in s or
        "polinton" in s
    ):
        return "DNA"

    return "Unknown"


def read_fasta(path):
    data = {}

    cur = None
    cur_class = "Unknown"
    seq_len = 0

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")

            if line.startswith(">"):

                if cur is not None:
                    data[cur] = (
                        cur_class,
                        broad_from_fine(cur_class),
                        seq_len
                    )

                token = line[1:].split()[0]

                cur = token.split("#", 1)[0]

                if "#" in token:
                    cur_class = token.split("#", 1)[1]
                else:
                    cur_class = "Unknown"

                seq_len = 0

            else:
                seq_len += len(line.strip())

    if cur is not None:
        data[cur] = (
            cur_class,
            broad_from_fine(cur_class),
            seq_len
        )

    return data


def read_dict(path, key="RepeatModeler_ID"):
    out = {}

    if not os.path.isfile(path):
        raise FileNotFoundError("Required evidence table is missing: %s" % path)

    with open(path) as fh:
        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )
        if not reader.fieldnames or key not in reader.fieldnames:
            raise ValueError("Evidence table %s must contain column %s" % (path, key))

        for row in reader:
            fid = fam_id(
                row.get(key, "")
            )

            if fid:
                out[fid] = row

    return out


def read_id_set(path):
    result = set()

    if not os.path.isfile(path):
        raise FileNotFoundError("Required Class II evidence table is missing: %s" % path)

    with open(path) as fh:
        for line in fh:

            if (
                not line.strip() or
                line.startswith("#")
            ):
                continue

            token = (
                line.rstrip("\n")
                .split("\t")[0]
            )

            low = token.lower()

            if (
                "query" in low or
                "repeatmodeler" in low
            ):
                continue

            fid = fam_id(token)

            if fid:
                result.add(fid)

    return result


def src_call(
    row,
    class_col,
    strength_col
):
    if not row:
        return "", "NONE", False

    cls = (
        row.get(class_col, "")
        or ""
    )

    strength = (
        row.get(strength_col, "")
        or "NONE"
    )

    conflict = (
        "CONFLICT" in strength or
        cls == "CONFLICT"
    )

    if cls in (
        "",
        "Unknown",
        "CONFLICT"
    ):
        usable_class = ""
    else:
        usable_class = cls

    return (
        usable_class,
        strength,
        conflict
    )


def combine_sources(calls):

    strong = []
    supportive = []
    any_conflict = False

    for cls, strength, conflict in calls:

        any_conflict = (
            any_conflict or conflict
        )

        if not cls:
            continue

        if strength.startswith("STRONG"):
            strong.append(cls)

        elif strength.startswith(
            "SUPPORTIVE"
        ):
            supportive.append(cls)

    strong_set = sorted(
        set(strong)
    )

    support_set = sorted(
        set(supportive)
    )

    if len(strong_set) > 1:

        return (
            "CONFLICT",
            "STRONG_CONFLICT",
            "YES",
            ",".join(strong_set),
            ",".join(support_set)
        )

    if len(strong_set) == 1:

        cls = strong_set[0]

        opposite_support = [
            x
            for x in support_set
            if x != cls
        ]

        if any_conflict:

            strength = (
                "STRONG_WITH_SOURCE_CONFLICT"
            )

            flag = "YES"

        elif opposite_support:

            strength = (
                "STRONG_WITH_WEAK_CONFLICT"
            )

            flag = "WEAK"

        else:

            strength = "STRONG"
            flag = "NO"

        return (
            cls,
            strength,
            flag,
            ",".join(strong_set),
            ",".join(support_set)
        )

    if len(support_set) > 1:

        return (
            "CONFLICT",
            "SUPPORTIVE_CONFLICT",
            "YES",
            "",
            ",".join(support_set)
        )

    if len(support_set) == 1:

        if any_conflict:
            strength = (
                "SUPPORTIVE_WITH_SOURCE_CONFLICT"
            )
            flag = "YES"
        else:
            strength = "SUPPORTIVE"
            flag = "NO"

        return (
            support_set[0],
            strength,
            flag,
            "",
            support_set[0]
        )

    if any_conflict:

        return (
            "",
            "SOURCE_CONFLICT_ONLY",
            "YES",
            "",
            ""
        )

    return (
        "",
        "NONE",
        "NO",
        "",
        ""
    )


def dimension_tier(
    cls,
    strength
):
    if cls in (
        "",
        "Unknown",
        "CONFLICT"
    ):
        return "NONE"

    if strength.startswith("STRONG"):
        return "STRONG"

    if strength.startswith(
        "SUPPORTIVE"
    ):
        return "SUPPORTIVE"

    return "NONE"


summary_rows = []
status_long = []


for entry in ENTRIES:

    sample = entry["sample_id"]
    universe = read_fasta(entry["repeatmodeler_fasta"])

    # --------------------------------------------------
    # Prior information only
    # --------------------------------------------------

    prior = read_dict(entry["prior_evidence"])

    # --------------------------------------------------
    # New v5 evidence
    # --------------------------------------------------

    dfam = read_dict(entry["dfam_family_summary"])

    repbase = read_dict(entry["repbase_family_summary"])

    tesorter = read_dict(entry["tesorter_family_summary"])

    ltr = read_dict(entry["ltr_reciprocal_family_best"])

    # --------------------------------------------------
    # MITE / Helitron reciprocal evidence
    # --------------------------------------------------

    mite_hc = read_id_set(entry["mite_hc_best"])
    mite_sup_all = read_id_set(entry["mite_supportive_best"])
    helitron_hc = read_id_set(entry["helitron_hc_best"])
    helitron_sup_all = read_id_set(entry["helitron_supportive_best"])

    # HC precedence
    mite_sup = (
        mite_sup_all - mite_hc
    )

    helitron_sup = (
        helitron_sup_all - helitron_hc
    )

    rows = []
    status_counts = defaultdict(int)

    for fid in sorted(universe):

        (
            original_class,
            original_broad,
            seq_len
        ) = universe[fid]

        p = prior.get(fid, {})
        d = dfam.get(fid, {})
        r = repbase.get(fid, {})
        t = tesorter.get(fid, {})
        l = ltr.get(fid, {})

        # ==================================================
        # HOMOLOGY DIMENSION
        # ==================================================

        dfam_call = src_call(
            d,
            "Dfam_v5_class",
            "Dfam_v5_strength"
        )

        repbase_call = src_call(
            r,
            "RepBase_v5_class",
            "RepBase_v5_strength"
        )

        (
            hom_class,
            hom_strength,
            hom_conflict,
            hom_strong_classes,
            hom_support_classes
        ) = combine_sources(
            [
                dfam_call,
                repbase_call
            ]
        )

        # ==================================================
        # DOMAIN DIMENSION
        # ==================================================

        (
            domain_class,
            domain_strength,
            domain_source_conflict
        ) = src_call(
            t,
            "TEsorter_v5_class",
            "TEsorter_v5_strength"
        )

        domain_conflict = (
            "YES"
            if domain_source_conflict
            else "NO"
        )

        # ==================================================
        # STRUCTURE DIMENSION
        # ==================================================

        structure_calls = []

        ltr_level = (
            l.get(
                "LTR_evidence_level",
                ""
            )
            if l
            else ""
        )

        if ltr_level == "DIRECT_LTR_HC":

            structure_calls.append(
                (
                    "LTR",
                    "STRONG",
                    False
                )
            )

        elif ltr_level in (
            "DIRECT_LTR_SUPPORTIVE",
            "INTERNAL_HC"
        ):

            structure_calls.append(
                (
                    "LTR",
                    "SUPPORTIVE",
                    False
                )
            )

        # INTERNAL_SUPPORTIVE is deliberately weak-only.

        if fid in mite_hc:

            structure_calls.append(
                (
                    "DNA",
                    "STRONG",
                    False
                )
            )

        elif fid in mite_sup:

            structure_calls.append(
                (
                    "DNA",
                    "SUPPORTIVE",
                    False
                )
            )

        # Helitron HC is capped at supportive.
        if fid in helitron_hc:

            structure_calls.append(
                (
                    "RC",
                    "SUPPORTIVE",
                    False
                )
            )

        # Helitron supportive = weak-only.

        (
            struct_class,
            struct_strength,
            struct_conflict,
            struct_strong_classes,
            struct_support_classes
        ) = combine_sources(
            structure_calls
        )

        # ==================================================
        # CROSS-DIMENSION DIAGNOSTIC
        # ==================================================

        dimensions = [
            (
                "HOMOLOGY",
                hom_class,
                hom_strength
            ),
            (
                "DOMAIN",
                domain_class,
                domain_strength
            ),
            (
                "STRUCTURE",
                struct_class,
                struct_strength
            )
        ]

        strong_calls = [
            (name, cls)
            for name, cls, strength
            in dimensions
            if dimension_tier(
                cls,
                strength
            ) == "STRONG"
        ]

        support_calls = [
            (name, cls)
            for name, cls, strength
            in dimensions
            if dimension_tier(
                cls,
                strength
            ) == "SUPPORTIVE"
        ]

        strong_classes = sorted(
            set(
                cls
                for _, cls
                in strong_calls
            )
        )

        support_classes = sorted(
            set(
                cls
                for _, cls
                in support_calls
            )
        )

        internal_review = (
            hom_conflict == "YES" or
            domain_conflict == "YES" or
            struct_conflict == "YES"
        )

        candidate = ""

        if len(strong_classes) > 1:

            status = (
                "STRONG_DIMENSION_CONFLICT"
            )

        elif len(strong_classes) == 1:

            candidate = strong_classes[0]

            nstrong = sum(
                1
                for _, cls
                in strong_calls
                if cls == candidate
            )

            opposite_support = [
                cls
                for _, cls
                in support_calls
                if cls != candidate
            ]

            if opposite_support:

                status = (
                    "STRONG_WITH_CROSS_SUPPORT_CONFLICT"
                )

            elif nstrong >= 3:

                status = (
                    "3D_STRONG_CONCORDANT"
                )

            elif nstrong == 2:

                status = (
                    "2D_STRONG_CONCORDANT"
                )

            else:

                status = (
                    "1D_STRONG_ONLY"
                )

        elif len(support_classes) > 1:

            status = (
                "SUPPORTIVE_DIMENSION_CONFLICT"
            )

        elif len(support_classes) == 1:

            candidate = support_classes[0]

            nsupport = sum(
                1
                for _, cls
                in support_calls
                if cls == candidate
            )

            if nsupport >= 2:
                status = (
                    "MULTI_SUPPORT_CONCORDANT"
                )
            else:
                status = "SUPPORT_ONLY"

        else:

            status = (
                "NO_EXTERNAL_EVIDENCE"
            )

        needs_review = (
            "CONFLICT" in status or
            internal_review
        )

        row = {
            "Sample":
                sample,

            "RepeatModeler_ID":
                fid,

            "Consensus_length":
                seq_len,

            "Original_class":
                original_class,

            "Original_broad":
                original_broad,

            # ------------------------------
            # PRIOR / COMPOSITE
            # ------------------------------

            "TEtrimmer_class":
                p.get(
                    "TEtrimmer_class",
                    ""
                ),

            "TEtrimmer_identity":
                p.get(
                    "TEtrimmer_identity",
                    ""
                ),

            "TEtrimmer_qcov":
                p.get(
                    "TEtrimmer_qcov",
                    ""
                ),

            "TEtrimmer_alignment_length":
                p.get(
                    "TEtrimmer_alignment_length",
                    ""
                ),

            "MCHelper_class":
                p.get(
                    "MCHelper_class",
                    ""
                ),

            "MCHelper_identity":
                p.get(
                    "MCHelper_identity",
                    ""
                ),

            "MCHelper_qcov":
                p.get(
                    "MCHelper_qcov",
                    ""
                ),

            "MCHelper_alignment_length":
                p.get(
                    "MCHelper_alignment_length",
                    ""
                ),

            # ------------------------------
            # HOMOLOGY
            # ------------------------------

            "Dfam_class":
                d.get(
                    "Dfam_v5_class",
                    ""
                ),

            "Dfam_strength":
                d.get(
                    "Dfam_v5_strength",
                    "NONE"
                ),

            "Dfam_conflict":
                d.get(
                    "Dfam_v5_conflict",
                    "NO"
                ),

            "RepBase_class":
                r.get(
                    "RepBase_v5_class",
                    ""
                ),

            "RepBase_strength":
                r.get(
                    "RepBase_v5_strength",
                    "NONE"
                ),

            "RepBase_conflict":
                r.get(
                    "RepBase_v5_conflict",
                    "NO"
                ),

            "Homology_class":
                hom_class,

            "Homology_strength":
                hom_strength,

            "Homology_conflict":
                hom_conflict,

            "Homology_strong_classes":
                hom_strong_classes,

            "Homology_supportive_classes":
                hom_support_classes,

            # ------------------------------
            # DOMAIN
            # ------------------------------

            "TEsorter_class":
                t.get(
                    "TEsorter_v5_class",
                    ""
                ),

            "TEsorter_strength":
                t.get(
                    "TEsorter_v5_strength",
                    "NONE"
                ),

            "TEsorter_conflict":
                t.get(
                    "TEsorter_v5_conflict",
                    "NO"
                ),

            "Domain_class":
                domain_class,

            "Domain_strength":
                domain_strength,

            "Domain_conflict":
                domain_conflict,

            # ------------------------------
            # STRUCTURE
            # ------------------------------

            "LTR_evidence_level":
                ltr_level,

            "LTR_target_class":
                l.get(
                    "LTR_target_class",
                    ""
                ),

            "LTR_qcov":
                l.get(
                    "qcov",
                    ""
                ),

            "LTR_scov":
                l.get(
                    "scov",
                    ""
                ),

            "MITE_HC":
                (
                    "YES"
                    if fid in mite_hc
                    else "NO"
                ),

            "MITE_supportive_only":
                (
                    "YES"
                    if fid in mite_sup
                    else "NO"
                ),

            "Helitron_HC":
                (
                    "YES"
                    if fid in helitron_hc
                    else "NO"
                ),

            "Helitron_supportive_only":
                (
                    "YES"
                    if fid in helitron_sup
                    else "NO"
                ),

            "Structure_class":
                struct_class,

            "Structure_strength":
                struct_strength,

            "Structure_conflict":
                struct_conflict,

            "Structure_strong_classes":
                struct_strong_classes,

            "Structure_supportive_classes":
                struct_support_classes,

            # ------------------------------
            # CROSS DIMENSIONS
            # ------------------------------

            "Strong_dimension_classes":
                ",".join(
                    strong_classes
                ),

            "Supportive_dimension_classes":
                ",".join(
                    support_classes
                ),

            "Cross_candidate_class":
                candidate,

            "Cross_dimension_status":
                status,

            "Needs_review":
                (
                    "YES"
                    if needs_review
                    else "NO"
                )
        }

        rows.append(row)
        status_counts[status] += 1

    # ==================================================
    # MASTER
    # ==================================================

    fields = list(
        rows[0].keys()
    )

    master_file = os.path.join(
        OUT,
        "%s.v5_evidence_master.tsv"
        % sample
    )

    with open(
        master_file,
        "w"
    ) as fo:

        writer = csv.DictWriter(
            fo,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writeheader()
        writer.writerows(rows)

    # ==================================================
    # REVIEW FILE
    # ==================================================

    review_file = os.path.join(
        REVIEW,
        "%s.v5_conflict_review.tsv"
        % sample
    )

    with open(
        review_file,
        "w"
    ) as fo:

        writer = csv.DictWriter(
            fo,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writeheader()

        writer.writerows([
            x
            for x in rows
            if x["Needs_review"] == "YES"
        ])

    def n_dim(
        prefix,
        kind
    ):
        if kind == "STRONG":

            return sum(
                1
                for x in rows
                if (
                    x[
                        prefix +
                        "_strength"
                    ].startswith(
                        "STRONG"
                    )
                    and
                    x[
                        prefix +
                        "_class"
                    ]
                    not in (
                        "",
                        "Unknown",
                        "CONFLICT"
                    )
                )
            )

        if kind == "SUPPORTIVE":

            return sum(
                1
                for x in rows
                if (
                    x[
                        prefix +
                        "_strength"
                    ].startswith(
                        "SUPPORTIVE"
                    )
                    and
                    x[
                        prefix +
                        "_class"
                    ]
                    not in (
                        "",
                        "Unknown",
                        "CONFLICT"
                    )
                )
            )

        return 0

    summary_rows.append({
        "Sample":
            sample,

        "Total_families":
            len(rows),

        "Homology_strong":
            n_dim(
                "Homology",
                "STRONG"
            ),

        "Homology_supportive":
            n_dim(
                "Homology",
                "SUPPORTIVE"
            ),

        "Homology_conflict":
            sum(
                1
                for x in rows
                if x[
                    "Homology_conflict"
                ] == "YES"
            ),

        "Domain_strong":
            n_dim(
                "Domain",
                "STRONG"
            ),

        "Domain_supportive":
            n_dim(
                "Domain",
                "SUPPORTIVE"
            ),

        "Domain_conflict":
            sum(
                1
                for x in rows
                if x[
                    "Domain_conflict"
                ] == "YES"
            ),

        "Structure_strong":
            n_dim(
                "Structure",
                "STRONG"
            ),

        "Structure_supportive":
            n_dim(
                "Structure",
                "SUPPORTIVE"
            ),

        "Structure_conflict":
            sum(
                1
                for x in rows
                if x[
                    "Structure_conflict"
                ] == "YES"
            ),

        "3D_STRONG_CONCORDANT":
            status_counts[
                "3D_STRONG_CONCORDANT"
            ],

        "2D_STRONG_CONCORDANT":
            status_counts[
                "2D_STRONG_CONCORDANT"
            ],

        "1D_STRONG_ONLY":
            status_counts[
                "1D_STRONG_ONLY"
            ],

        "STRONG_WITH_CROSS_SUPPORT_CONFLICT":
            status_counts[
                "STRONG_WITH_CROSS_SUPPORT_CONFLICT"
            ],

        "STRONG_DIMENSION_CONFLICT":
            status_counts[
                "STRONG_DIMENSION_CONFLICT"
            ],

        "MULTI_SUPPORT_CONCORDANT":
            status_counts[
                "MULTI_SUPPORT_CONCORDANT"
            ],

        "SUPPORT_ONLY":
            status_counts[
                "SUPPORT_ONLY"
            ],

        "SUPPORTIVE_DIMENSION_CONFLICT":
            status_counts[
                "SUPPORTIVE_DIMENSION_CONFLICT"
            ],

        "NO_EXTERNAL_EVIDENCE":
            status_counts[
                "NO_EXTERNAL_EVIDENCE"
            ],

        "Needs_review":
            sum(
                1
                for x in rows
                if x[
                    "Needs_review"
                ] == "YES"
            )
    })

    for status, count in sorted(
        status_counts.items()
    ):

        status_long.append({
            "Sample":
                sample,

            "Status":
                status,

            "Count":
                count
        })


# ======================================================
# Aggregate summary across manifest entries
# ======================================================

summary_file = os.path.join(
    OUT,
    "all_libraries.v5_evidence_master.summary.tsv"
)

with open(
    summary_file,
    "w"
) as fo:

    fields = list(
        summary_rows[0].keys()
    )

    writer = csv.DictWriter(
        fo,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        summary_rows
    )


status_file = os.path.join(
    OUT,
    "all_libraries.v5_cross_status.tsv"
)

with open(
    status_file,
    "w"
) as fo:

    writer = csv.DictWriter(
        fo,
        fieldnames=[
            "Sample",
            "Status",
            "Count"
        ],
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        status_long
    )


print("Finished")
print(summary_file)
