#!/usr/bin/env python3

import os
import argparse
import csv
import math
import statistics
from collections import defaultdict, deque

parser = argparse.ArgumentParser(description="Audit graph structure of large initial Pan-TE components.")
parser.add_argument("--workdir", required=True, help="Pan-TE working directory")
args = parser.parse_args()
ROOT = os.path.abspath(args.workdir)

MEMBERS = os.path.join(
    ROOT,
    "02_cluster",
    "PanTE_cluster_members.tsv"
)

OUT = os.path.join(
    ROOT,
    "02_cluster",
    "PanTE_large_cluster_graph_qc.tsv"
)

CLASSES = [
    "DNA",
    "LTR",
    "LINE",
    "SINE",
    "RC"
]

MIN_CLUSTER_SIZE = 10


# ==========================================================
# Read cluster membership
# ==========================================================

cluster_members = defaultdict(set)
cluster_class = {}

with open(MEMBERS) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:

        pid = row["PanTE_ID"]
        broad = row["Broad_class"]
        uid = row["Unified_ID"].split("#", 1)[0]

        cluster_members[pid].add(uid)
        cluster_class[pid] = broad


uid_to_cluster = {}

for pid, members in cluster_members.items():

    for uid in members:
        uid_to_cluster[uid] = pid


large_clusters = {
    pid
    for pid, members in cluster_members.items()
    if len(members) >= MIN_CLUSTER_SIZE
}


# ==========================================================
# Read accepted high-confidence edges
#
# all-vs-all BLAST contains q->s and s->q.
# Deduplicate them into one undirected edge.
# ==========================================================

edges = defaultdict(dict)

for broad in CLASSES:

    edgefile = os.path.join(
        ROOT,
        "02_cluster",
        f"{broad}.high_confidence_edges.tsv"
    )

    if not os.path.exists(edgefile):
        continue

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

            pid_q = uid_to_cluster.get(q)
            pid_s = uid_to_cluster.get(s)

            if pid_q is None or pid_q != pid_s:
                continue

            pid = pid_q

            if pid not in large_clusters:
                continue

            pair = tuple(sorted((q, s)))

            rec = {
                "identity": float(row["Identity"]),
                "qcov": float(row["Qcov"]),
                "scov": float(row["Scov"]),
                "bitscore": float(row["Bitscore"])
            }

            # q->s and s->q may both exist.
            # Keep the record with the highest bitscore.
            old = edges[pid].get(pair)

            if old is None or rec["bitscore"] > old["bitscore"]:
                edges[pid][pair] = rec


# ==========================================================
# Graph utilities
# ==========================================================

def graph_stats(nodes, pair_records):

    adj = {
        x: set()
        for x in nodes
    }

    for a, b in pair_records:

        adj[a].add(b)
        adj[b].add(a)

    degrees = [
        len(adj[x])
        for x in nodes
    ]

    diameter = 0
    all_distances = []

    for source in nodes:

        dist = {
            source: 0
        }

        q = deque([source])

        while q:

            u = q.popleft()

            for v in adj[u]:

                if v not in dist:

                    dist[v] = dist[u] + 1
                    q.append(v)

        # Cluster should be connected.
        if len(dist) != len(nodes):

            return {
                "connected": 0,
                "diameter": -1,
                "mean_path": -1,
                "degrees": degrees
            }

        diameter = max(
            diameter,
            max(dist.values())
        )

        for target, d in dist.items():

            if target != source:
                all_distances.append(d)

    mean_path = (
        sum(all_distances) / len(all_distances)
        if all_distances
        else 0
    )

    return {
        "connected": 1,
        "diameter": diameter,
        "mean_path": mean_path,
        "degrees": degrees
    }


# ==========================================================
# Produce QC
# ==========================================================

rows = []

for pid in sorted(
    large_clusters,
    key=lambda x: (
        -len(cluster_members[x]),
        x
    )
):

    nodes = sorted(cluster_members[pid])

    n = len(nodes)

    pair_records = edges.get(pid, {})

    e = len(pair_records)

    possible = n * (n - 1) // 2

    density = (
        e / possible
        if possible > 0
        else 0
    )

    stats = graph_stats(
        nodes,
        pair_records
    )

    degrees = stats["degrees"]

    identities = [
        x["identity"]
        for x in pair_records.values()
    ]

    recip_cov = [
        min(
            x["qcov"],
            x["scov"]
        )
        for x in pair_records.values()
    ]

    row = {
        "PanTE_ID": pid,
        "Broad_class": cluster_class[pid],
        "Cluster_size": n,
        "Unique_direct_edges": e,
        "Possible_pairs": possible,
        "Edge_density": density,
        "Connected": stats["connected"],
        "Diameter": stats["diameter"],
        "Mean_shortest_path": stats["mean_path"],
        "Min_degree": min(degrees) if degrees else 0,
        "Median_degree": statistics.median(degrees) if degrees else 0,
        "Mean_degree": statistics.mean(degrees) if degrees else 0,
        "Max_degree": max(degrees) if degrees else 0,
        "Mean_identity": statistics.mean(identities) if identities else 0,
        "Median_identity": statistics.median(identities) if identities else 0,
        "Min_identity": min(identities) if identities else 0,
        "Mean_min_recip_cov": statistics.mean(recip_cov) if recip_cov else 0,
        "Median_min_recip_cov": statistics.median(recip_cov) if recip_cov else 0,
        "Min_min_recip_cov": min(recip_cov) if recip_cov else 0
    }

    rows.append(row)


fields = [
    "PanTE_ID",
    "Broad_class",
    "Cluster_size",
    "Unique_direct_edges",
    "Possible_pairs",
    "Edge_density",
    "Connected",
    "Diameter",
    "Mean_shortest_path",
    "Min_degree",
    "Median_degree",
    "Mean_degree",
    "Max_degree",
    "Mean_identity",
    "Median_identity",
    "Min_identity",
    "Mean_min_recip_cov",
    "Median_min_recip_cov",
    "Min_min_recip_cov"
]

with open(OUT, "w") as fo:

    writer = csv.DictWriter(
        fo,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()

    for row in rows:

        writer.writerow(row)


print("Large clusters:", len(rows))
print("Output:", OUT)

print("")
print("Top clusters:")
print(
    "PanTE_ID\tClass\tN\tEdges\tDensity\t"
    "Diameter\tMedianDegree\tMedianIdentity\tMedianRecipCov"
)

for r in rows[:30]:

    print(
        "%s\t%s\t%d\t%d\t%.4f\t%d\t%.1f\t%.2f\t%.2f"
        % (
            r["PanTE_ID"],
            r["Broad_class"],
            r["Cluster_size"],
            r["Unique_direct_edges"],
            r["Edge_density"],
            r["Diameter"],
            r["Median_degree"],
            r["Median_identity"],
            r["Median_min_recip_cov"]
        )
    )