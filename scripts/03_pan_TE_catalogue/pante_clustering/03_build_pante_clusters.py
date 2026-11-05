#!/usr/bin/env python3

import os
import argparse
import csv
from collections import defaultdict

parser = argparse.ArgumentParser(description="Build initial Pan-TE connected components from BLASTn results.")
parser.add_argument("--workdir", required=True, help="Pan-TE working directory")
args = parser.parse_args()
ROOT = os.path.abspath(args.workdir)
os.makedirs(os.path.join(ROOT, "02_cluster"), exist_ok=True)

CLASSES = [
    "DNA",
    "LTR",
    "LINE",
    "SINE",
    "RC"
]

MIN_ID = 80.0
MIN_QCOV = 80.0
MIN_SCOV = 80.0
MIN_BP = 100
MAX_E = 1e-20


def union_length(intervals):

    if not intervals:
        return 0

    a = sorted(intervals)

    s, e = a[0]
    total = 0

    for x, y in a[1:]:

        if x <= e + 1:
            e = max(e, y)

        else:
            total += e - s + 1
            s, e = x, y

    total += e - s + 1

    return total


def read_ids(fasta):

    ids = []

    with open(fasta) as fh:

        for line in fh:

            if line.startswith(">"):

                token = line[1:].split()[0]

                uid = token.split("#", 1)[0]

                ids.append(uid)

    return ids


class DSU:

    def __init__(self, nodes):

        self.parent = {
            x: x for x in nodes
        }

        self.rank = {
            x: 0 for x in nodes
        }

    def find(self, x):

        while self.parent[x] != x:

            self.parent[x] = self.parent[
                self.parent[x]
            ]

            x = self.parent[x]

        return x

    def union(self, a, b):

        ra = self.find(a)
        rb = self.find(b)

        if ra == rb:
            return

        if self.rank[ra] < self.rank[rb]:

            self.parent[ra] = rb

        elif self.rank[ra] > self.rank[rb]:

            self.parent[rb] = ra

        else:

            self.parent[rb] = ra
            self.rank[ra] += 1


all_clusters = []

edge_summary = []


for broad in CLASSES:

    fasta = os.path.join(
        ROOT,
        "00_input",
        "%s.fa" % broad
    )

    blast = os.path.join(
        ROOT,
        "01_blast",
        "%s.all_vs_all.tsv" % broad
    )

    ids = read_ids(fasta)

    hsps = defaultdict(list)

    with open(blast) as fh:

        for line in fh:

            f = line.rstrip().split("\t")

            if len(f) != 12:
                continue

            # BLAST qseqid/sseqid retain the RepeatMasker-style
            # classification suffix, e.g.
            # UNI_v52_000002#DNA/Unknown.
            # PanTE node IDs use only UNI_v52_000002, so normalize here.
            q = f[0].split("#", 1)[0]
            s = f[1].split("#", 1)[0]

            if q == s:
                continue

            pid = float(f[2])
            aln = int(f[3])

            qs = min(
                int(f[4]),
                int(f[5])
            )

            qe = max(
                int(f[4]),
                int(f[5])
            )

            ss = min(
                int(f[6]),
                int(f[7])
            )

            se = max(
                int(f[6]),
                int(f[7])
            )

            ev = float(f[8])
            bit = float(f[9])

            qlen = int(f[10])
            slen = int(f[11])

            hsps[(q, s)].append({
                "pid": pid,
                "aln": aln,
                "qs": qs,
                "qe": qe,
                "ss": ss,
                "se": se,
                "ev": ev,
                "bit": bit,
                "qlen": qlen,
                "slen": slen
            })


    dsu = DSU(ids)

    n_pairs = 0
    n_edges = 0

    edge_file = os.path.join(
        ROOT,
        "02_cluster",
        "%s.high_confidence_edges.tsv"
        % broad
    )

    with open(edge_file, "w") as fo:

        fo.write(
            "Query\tSubject\tIdentity\t"
            "Qcov\tScov\tQbp\tSbp\t"
            "Best_Evalue\tBitscore\n"
        )

        for (q, s), arr in hsps.items():

            n_pairs += 1

            qbp = union_length([
                (x["qs"], x["qe"])
                for x in arr
            ])

            sbp = union_length([
                (x["ss"], x["se"])
                for x in arr
            ])

            qlen = arr[0]["qlen"]
            slen = arr[0]["slen"]

            qcov = (
                100.0 * qbp / qlen
                if qlen else 0
            )

            scov = (
                100.0 * sbp / slen
                if slen else 0
            )

            total_aln = sum(
                x["aln"]
                for x in arr
            )

            ident = (
                sum(
                    x["pid"] * x["aln"]
                    for x in arr
                )
                / total_aln
                if total_aln
                else 0
            )

            best_e = min(
                x["ev"] for x in arr
            )

            bitscore = sum(
                x["bit"] for x in arr
            )

            if (
                ident >= MIN_ID
                and qcov >= MIN_QCOV
                and scov >= MIN_SCOV
                and min(qbp, sbp) >= MIN_BP
                and best_e <= MAX_E
            ):

                n_edges += 1

                dsu.union(q, s)

                fo.write(
                    "%s\t%s\t%.3f\t%.3f\t%.3f\t"
                    "%d\t%d\t%.3e\t%.3f\n"
                    % (
                        q,
                        s,
                        ident,
                        qcov,
                        scov,
                        qbp,
                        sbp,
                        best_e,
                        bitscore
                    )
                )


    groups = defaultdict(list)

    for uid in ids:

        groups[
            dsu.find(uid)
        ].append(uid)


    groups = sorted(
        groups.values(),
        key=lambda x: (
            -len(x),
            sorted(x)[0]
        )
    )


    for g in groups:

        all_clusters.append(
            (
                broad,
                sorted(g)
            )
        )


    edge_summary.append(
        (
            broad,
            len(ids),
            n_pairs,
            n_edges,
            len(groups),
            sum(
                1
                for g in groups
                if len(g) == 1
            )
        )
    )


# =========================================================
# Stable PanTE IDs
# =========================================================

outfile = os.path.join(
    ROOT,
    "02_cluster",
    "PanTE_cluster_members.tsv"
)

summary = os.path.join(
    ROOT,
    "02_cluster",
    "PanTE_cluster_summary.tsv"
)

with open(outfile, "w") as fo:

    fo.write(
        "PanTE_ID\tBroad_class\t"
        "Cluster_size\tUnified_ID\n"
    )

    counter = 1

    for broad, members in all_clusters:

        pid = "PanTE%06d" % counter

        for uid in members:

            fo.write(
                "%s\t%s\t%d\t%s\n"
                % (
                    pid,
                    broad,
                    len(members),
                    uid
                )
            )

        counter += 1


with open(summary, "w") as fo:

    fo.write(
        "Broad_class\tConsensus_n\t"
        "BLAST_pairs\tAccepted_edges\t"
        "PanTE_clusters\tSingleton_clusters\n"
    )

    for row in edge_summary:

        fo.write(
            "\t".join(
                map(str, row)
            )
            + "\n"
        )


print("Finished.")
print(outfile)
print(summary)