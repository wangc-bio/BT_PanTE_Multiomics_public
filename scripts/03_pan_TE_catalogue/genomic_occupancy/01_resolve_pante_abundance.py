#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import csv
import heapq
from collections import defaultdict

if len(sys.argv) != 5:
    sys.stderr.write(
        "Usage: 01_resolve_pante_abundance.py "
        "<PanTE_members.tsv> <RepeatMasker.out> "
        "<sample> <output_prefix>\n"
    )
    sys.exit(1)

MEMBERS = sys.argv[1]
RMOUT   = sys.argv[2]
SAMPLE  = sys.argv[3]
PREFIX  = sys.argv[4]

OUT_AB = PREFIX + ".PanTE_resolved_abundance.tsv"
OUT_QC = PREFIX + ".resolution_qc.tsv"


# ==========================================================
# Read PanTE membership
# ==========================================================

uid_to_pante = {}
pante_class = {}

with open(MEMBERS) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:

        uid = row["Unified_ID"].split("#", 1)[0]
        pid = row["PanTE_ID"]

        uid_to_pante[uid] = pid
        pante_class[pid] = row["Broad_class"]


# ==========================================================
# Read RepeatMasker hits
#
# Store ALL annotations because Unknown / simple repeats /
# nonTE annotations must also participate in overlap
# resolution.
#
# tuple:
# start0, end, primary, score, length, PanTE_ID/OTHER
# ==========================================================

by_chr = defaultdict(list)

total_hits = 0
primary_hits = 0
star_hits = 0
classified_hits = 0
other_hits = 0

with open(RMOUT) as fh:

    for line in fh:

        line = line.strip()

        if not line:
            continue

        f = line.split()

        if len(f) < 14:
            continue

        try:
            score = float(f[0])
        except ValueError:
            continue

        total_hits += 1

        chrom = f[4]

        try:
            s = int(f[5])
            e = int(f[6])
        except ValueError:
            continue

        start0 = min(s, e) - 1
        end = max(s, e)

        if end <= start0:
            continue

        length = end - start0

        repeat_id = f[9].split("#", 1)[0]

        pid = uid_to_pante.get(
            repeat_id,
            "OTHER"
        )

        if pid == "OTHER":
            other_hits += 1
        else:
            classified_hits += 1

        is_star = (
            len(f) > 14
            and f[-1] == "*"
        )

        if is_star:
            primary = 0
            star_hits += 1
        else:
            primary = 1
            primary_hits += 1

        by_chr[chrom].append(
            (
                start0,
                end,
                primary,
                score,
                length,
                pid
            )
        )


# ==========================================================
# Sweep-line overlap resolution
#
# For every genomic segment:
#
# primary/non-* > secondary/*
# then greater SW score
# then greater hit length
#
# Every base is therefore assigned at most once.
# ==========================================================

occupied_bp = defaultdict(int)
resolved_segments = defaultdict(int)

total_repeat_union_bp = 0
other_resolved_bp = 0


for chrom in sorted(by_chr):

    intervals = by_chr[chrom]

    intervals.sort(
        key=lambda x: (
            x[0],
            x[1]
        )
    )

    n = len(intervals)

    if n == 0:
        continue

    active = set()

    # (end, interval_index)
    end_heap = []

    # (-primary, -score, -length, interval_index)
    winner_heap = []

    i = 0
    pos = intervals[0][0]

    while i < n or active:

        # ----------------------------------------------
        # Remove intervals that have already ended
        # ----------------------------------------------

        while (
            end_heap
            and end_heap[0][0] <= pos
        ):

            _, idx = heapq.heappop(
                end_heap
            )

            active.discard(idx)


        # ----------------------------------------------
        # Add intervals starting at current position
        # ----------------------------------------------

        while (
            i < n
            and intervals[i][0] <= pos
        ):

            start0, end, primary, score, length, pid = \
                intervals[i]

            active.add(i)

            heapq.heappush(
                end_heap,
                (
                    end,
                    i
                )
            )

            heapq.heappush(
                winner_heap,
                (
                    -primary,
                    -score,
                    -length,
                    i
                )
            )

            i += 1


        # ----------------------------------------------
        # Clean lazy winner heap
        # ----------------------------------------------

        while (
            winner_heap
            and winner_heap[0][3] not in active
        ):

            heapq.heappop(
                winner_heap
            )


        # ----------------------------------------------
        # No active intervals: jump to next hit
        # ----------------------------------------------

        if not active:

            if i < n:

                pos = intervals[i][0]
                continue

            break


        # ----------------------------------------------
        # Next coordinate where annotation can change
        # ----------------------------------------------

        while (
            end_heap
            and end_heap[0][1] not in active
        ):

            heapq.heappop(
                end_heap
            )

        next_end = end_heap[0][0]

        if i < n:
            next_start = intervals[i][0]
        else:
            next_start = float("inf")

        next_pos = min(
            next_end,
            next_start
        )

        if next_pos <= pos:
            pos += 1
            continue


        # ----------------------------------------------
        # Current winning annotation
        # ----------------------------------------------

        winner_idx = winner_heap[0][3]

        (
            start0,
            end,
            primary,
            score,
            length,
            pid
        ) = intervals[winner_idx]

        span = next_pos - pos

        total_repeat_union_bp += span

        if pid == "OTHER":

            other_resolved_bp += span

        else:

            occupied_bp[pid] += span
            resolved_segments[pid] += 1

        pos = next_pos


# ==========================================================
# Write PanTE abundance
# ==========================================================

with open(OUT_AB, "w") as fo:

    fo.write(
        "PanTE_ID\tBroad_class\t"
        "Occupied_bp\tResolved_segments\tSample\n"
    )

    for pid in sorted(occupied_bp):

        fo.write(
            "%s\t%s\t%d\t%d\t%s\n"
            % (
                pid,
                pante_class[pid],
                occupied_bp[pid],
                resolved_segments[pid],
                SAMPLE
            )
        )


pante_total_bp = sum(
    occupied_bp.values()
)


# ==========================================================
# QC
# ==========================================================

with open(OUT_QC, "w") as fo:

    fo.write(
        "Sample\tTotal_RM_hits\tPrimary_hits\tStar_hits\t"
        "Classified_hits\tOther_hits\t"
        "Repeat_union_bp\tPanTE_resolved_bp\t"
        "Other_resolved_bp\tDetected_PanTE\n"
    )

    fo.write(
        "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
        % (
            SAMPLE,
            total_hits,
            primary_hits,
            star_hits,
            classified_hits,
            other_hits,
            total_repeat_union_bp,
            pante_total_bp,
            other_resolved_bp,
            len(occupied_bp)
        )
    )


print("============================================")
print("Sample              :", SAMPLE)
print("Total RM hits       :", total_hits)
print("Primary hits        :", primary_hits)
print("Star hits           :", star_hits)
print("Detected PanTE      :", len(occupied_bp))
print("Repeat union bp     :", total_repeat_union_bp)
print("PanTE resolved bp   :", pante_total_bp)
print("Other resolved bp   :", other_resolved_bp)
print("")
print("Abundance:", OUT_AB)
print("QC       :", OUT_QC)