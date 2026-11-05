#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import csv
import statistics
from collections import defaultdict

parser = argparse.ArgumentParser(description="Refine initial Pan-TE components using representative-centered greedy clustering.")
parser.add_argument("--workdir", required=True, help="Pan-TE working directory")
args = parser.parse_args()
ROOT = os.path.abspath(args.workdir)

MEMBERS = os.path.join(
    ROOT,
    "02_cluster",
    "PanTE_cluster_members.tsv"
)

CLASSES = [
    "DNA",
    "LTR",
    "LINE",
    "SINE",
    "RC"
]

CLASS_ORDER = {
    "DNA": 1,
    "LTR": 2,
    "LINE": 3,
    "SINE": 4,
    "RC": 5
}

OUT_MEMBERS = os.path.join(
    ROOT,
    "02_cluster",
    "PanTE_refined_cluster_members.tsv"
)

OUT_SUMMARY = os.path.join(
    ROOT,
    "02_cluster",
    "PanTE_refined_cluster_summary.tsv"
)

OUT_AUDIT = os.path.join(
    ROOT,
    "02_cluster",
    "PanTE_refined_component_audit.tsv"
)

OUT_CLUSTER_QC = os.path.join(
    ROOT,
    "02_cluster",
    "PanTE_refined_cluster_qc.tsv"
)


# ==========================================================
# 1. Read original connected components
# ==========================================================

old_members = defaultdict(set)
old_class = {}

with open(MEMBERS) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:

        pid = row["PanTE_ID"]
        broad = row["Broad_class"]
        uid = row["Unified_ID"].split("#", 1)[0]

        old_members[pid].add(uid)
        old_class[pid] = broad


uid_to_old = {}

for pid, members in old_members.items():

    for uid in members:
        uid_to_old[uid] = pid


# ==========================================================
# 2. Read and deduplicate high-confidence edges
# ==========================================================

# adjacency[node][neighbor] = edge information
adjacency = defaultdict(dict)

for broad in CLASSES:

    edgefile = os.path.join(
        ROOT,
        "02_cluster",
        "%s.high_confidence_edges.tsv" % broad
    )

    if not os.path.exists(edgefile):
        continue

    pair_best = {}

    with open(edgefile) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        for row in reader:

            q = row["Query"].split("#", 1)[0]
            s = row["Subject"].split("#", 1)[0]

            if q == s:
                continue

            # Only relationships belonging to the same
            # original connected component are relevant.
            oq = uid_to_old.get(q)
            os_ = uid_to_old.get(s)

            if oq is None or oq != os_:
                continue

            pair = tuple(sorted((q, s)))

            rec = {
                "identity": float(row["Identity"]),
                "qcov": float(row["Qcov"]),
                "scov": float(row["Scov"]),
                "qbp": int(row["Qbp"]),
                "sbp": int(row["Sbp"]),
                "evalue": float(row["Best_Evalue"]),
                "bitscore": float(row["Bitscore"])
            }

            old = pair_best.get(pair)

            # q->s and s->q may both occur.
            # Retain the record with the greater bit score.
            if old is None or rec["bitscore"] > old["bitscore"]:
                pair_best[pair] = rec

    for (a, b), rec in pair_best.items():

        adjacency[a][b] = rec
        adjacency[b][a] = rec


# ==========================================================
# 3. Representative-centered greedy clustering
#
# Within each original connected component:
#
# - choose the node with the largest number of direct
#   high-confidence neighbours among remaining nodes;
#
# - ties are resolved using summed bit score and then ID;
#
# - representative + all directly connected remaining
#   neighbours form one refined PanTE cluster;
#
# - remove these sequences and repeat.
#
# Therefore every non-representative member has a direct
# 80/80 high-confidence relationship to its representative.
# ==========================================================

refined_records = []
component_audit = []
cluster_qc = []

counter = 1


def old_pid_sort_key(pid):

    broad = old_class[pid]

    try:
        number = int(
            pid.replace("PanTE", "")
        )
    except Exception:
        number = 10**12

    return (
        CLASS_ORDER.get(broad, 99),
        number,
        pid
    )


for old_pid in sorted(
    old_members,
    key=old_pid_sort_key
):

    broad = old_class[old_pid]

    original_nodes = set(
        old_members[old_pid]
    )

    remaining = set(
        original_nodes
    )

    subclusters = []

    while remaining:

        # ----------------------------------------------
        # Find best representative among remaining nodes
        # ----------------------------------------------

        candidates = []

        for node in remaining:

            neighbours = (
                set(adjacency.get(node, {}))
                & remaining
            )

            degree = len(neighbours)

            summed_bitscore = sum(
                adjacency[node][x]["bitscore"]
                for x in neighbours
            )

            candidates.append(
                (
                    -degree,
                    -summed_bitscore,
                    node
                )
            )

        candidates.sort()

        representative = candidates[0][2]

        direct_neighbours = (
            set(
                adjacency.get(
                    representative,
                    {}
                )
            )
            & remaining
        )

        direct_neighbours.discard(
            representative
        )

        members = {
            representative
        } | direct_neighbours

        # ----------------------------------------------
        # QC statistics for representative relationships
        # ----------------------------------------------

        identities = []
        reciprocal_cov = []
        bitscores = []

        for member in sorted(members):

            if member == representative:
                continue

            rec = adjacency[
                representative
            ][member]

            identities.append(
                rec["identity"]
            )

            reciprocal_cov.append(
                min(
                    rec["qcov"],
                    rec["scov"]
                )
            )

            bitscores.append(
                rec["bitscore"]
            )

        refined_id = "PanTE_R%06d" % counter
        counter += 1

        subclusters.append(
            (
                refined_id,
                representative,
                sorted(members)
            )
        )

        for member in sorted(members):

            refined_records.append({
                "PanTE_ID": refined_id,
                "Broad_class": broad,
                "Cluster_size": len(members),
                "Representative_ID": representative,
                "Original_component": old_pid,
                "Unified_ID": member
            })

        cluster_qc.append({
            "PanTE_ID": refined_id,
            "Broad_class": broad,
            "Cluster_size": len(members),
            "Representative_ID": representative,
            "Original_component": old_pid,
            "Direct_members": max(
                0,
                len(members) - 1
            ),
            "Min_rep_identity":
                min(identities)
                if identities else 100.0,
            "Median_rep_identity":
                statistics.median(identities)
                if identities else 100.0,
            "Min_rep_recip_cov":
                min(reciprocal_cov)
                if reciprocal_cov else 100.0,
            "Median_rep_recip_cov":
                statistics.median(reciprocal_cov)
                if reciprocal_cov else 100.0
        })

        remaining -= members


    sizes = sorted(
        [
            len(x[2])
            for x in subclusters
        ],
        reverse=True
    )

    component_audit.append({
        "Original_component": old_pid,
        "Broad_class": broad,
        "Original_size": len(original_nodes),
        "Refined_cluster_n": len(subclusters),
        "Refined_sizes": ",".join(
            map(str, sizes)
        ),
        "Largest_refined_cluster":
            max(sizes) if sizes else 0,
        "Was_split":
            "YES"
            if len(subclusters) > 1
            else "NO"
    })


# ==========================================================
# 4. Write membership table
# ==========================================================

with open(OUT_MEMBERS, "w") as fo:

    fields = [
        "PanTE_ID",
        "Broad_class",
        "Cluster_size",
        "Representative_ID",
        "Original_component",
        "Unified_ID"
    ]

    writer = csv.DictWriter(
        fo,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()

    for row in refined_records:
        writer.writerow(row)


# ==========================================================
# 5. Component-level split audit
# ==========================================================

with open(OUT_AUDIT, "w") as fo:

    fields = [
        "Original_component",
        "Broad_class",
        "Original_size",
        "Refined_cluster_n",
        "Refined_sizes",
        "Largest_refined_cluster",
        "Was_split"
    ]

    writer = csv.DictWriter(
        fo,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()

    for row in component_audit:
        writer.writerow(row)


# ==========================================================
# 6. Refined cluster QC
# ==========================================================

with open(OUT_CLUSTER_QC, "w") as fo:

    fields = [
        "PanTE_ID",
        "Broad_class",
        "Cluster_size",
        "Representative_ID",
        "Original_component",
        "Direct_members",
        "Min_rep_identity",
        "Median_rep_identity",
        "Min_rep_recip_cov",
        "Median_rep_recip_cov"
    ]

    writer = csv.DictWriter(
        fo,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()

    for row in cluster_qc:
        writer.writerow(row)


# ==========================================================
# 7. Summary by TE class
# ==========================================================

cluster_info = defaultdict(dict)

for row in refined_records:

    pid = row["PanTE_ID"]

    cluster_info[pid]["broad"] = row[
        "Broad_class"
    ]

    cluster_info[pid]["size"] = row[
        "Cluster_size"
    ]


with open(OUT_SUMMARY, "w") as fo:

    fo.write(
        "Broad_class\t"
        "Consensus_n\t"
        "Original_components\t"
        "Refined_clusters\t"
        "Singleton_clusters\t"
        "Multi_member_clusters\t"
        "Median_cluster_size\t"
        "Max_cluster_size\n"
    )

    for broad in CLASSES:

        consensus_n = sum(
            1
            for row in refined_records
            if row["Broad_class"] == broad
        )

        original_n = sum(
            1
            for pid in old_members
            if old_class[pid] == broad
        )

        sizes = [
            info["size"]
            for info in cluster_info.values()
            if info["broad"] == broad
        ]

        refined_n = len(sizes)

        singleton_n = sum(
            1
            for x in sizes
            if x == 1
        )

        multi_n = sum(
            1
            for x in sizes
            if x > 1
        )

        median_size = (
            statistics.median(sizes)
            if sizes
            else 0
        )

        max_size = (
            max(sizes)
            if sizes
            else 0
        )

        fo.write(
            "%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\n"
            % (
                broad,
                consensus_n,
                original_n,
                refined_n,
                singleton_n,
                multi_n,
                median_size,
                max_size
            )
        )


print("==========================================")
print("Representative-centered PanTE refinement")
print("==========================================")
print("Consensus records :", len(refined_records))
print("Final clusters    :", len(cluster_info))
print("")
print("Membership :", OUT_MEMBERS)
print("Summary    :", OUT_SUMMARY)
print("Audit      :", OUT_AUDIT)
print("QC         :", OUT_CLUSTER_QC)