#!/usr/bin/env python
# Public-release note: source script from the executed analysis. Replace the original local absolute paths and executable locations with paths in your environment before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import re
import shutil
from collections import Counter, defaultdict


BASE = os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52")

WGBS_ANNO = os.path.join(
    BASE,
    "06_TE_WGBS/01_annotation"
)

ATAC = os.path.join(
    BASE,
    "07_TE_ATAC"
)

OUTDIR = os.path.join(
    ATAC,
    "01_annotation"
)

TMPDIR = os.path.join(
    OUTDIR,
    "tmp_chr_segments"
)

os.makedirs(OUTDIR, exist_ok=True)

if os.path.exists(TMPDIR):
    shutil.rmtree(TMPDIR)

os.makedirs(TMPDIR)


ORDER_FILE = os.path.join(
    ATAC,
    "00_QC/ATAC_reference_chromosomes.tsv"
)


INPUTS = {
    "BSB": os.path.join(
        WGBS_ANNO,
        "BSB.PanTE.WGBS_intervals.tsv"
    ),

    "TC": os.path.join(
        WGBS_ANNO,
        "TC.PanTE.WGBS_intervals.tsv"
    ),

    "BT": os.path.join(
        WGBS_ANNO,
        "BT.PanTE.WGBS_intervals.tsv"
    )
}


# ============================================================
# Reference chromosome order
# ============================================================

order = defaultdict(list)

with open(ORDER_FILE) as fh:

    header = fh.readline().rstrip("\n").split("\t")

    idx_group = header.index("Group")
    idx_chr = header.index("Chromosome")

    for line in fh:

        p = line.rstrip("\n").split("\t")

        order[p[idx_group]].append(
            p[idx_chr]
        )


def is_main_chr(group, chrom):

    if group in ("BSB", "TC"):

        return re.match(
            r"^chr[0-9][0-9]$",
            chrom
        ) is not None

    if group == "BT":

        return re.match(
            r"^chr[0-9][0-9][MP]$",
            chrom
        ) is not None

    return False


def state_from_chr(group, chrom):

    if group == "BSB":
        return "BSB"

    if group == "TC":
        return "TC"

    if group == "BT":

        if chrom.endswith("M"):
            return "BTM"

        if chrom.endswith("P"):
            return "BTP"

    raise RuntimeError(
        "Cannot determine state: {} {}".format(
            group,
            chrom
        )
    )


# ============================================================
# Sweep one chromosome
# ============================================================

def process_chr(
    group,
    chrom,
    intervals,
    outfile,
    bp_stat,
    seg_stat
):

    events = []

    for start1, end1, pante in intervals:

        # GTF / WGBS interval:
        # 1-based inclusive
        #
        # BED:
        # 0-based half-open

        start0 = start1 - 1
        end0 = end1

        if end0 <= start0:
            continue

        events.append(
            (
                start0,
                1,
                pante
            )
        )

        events.append(
            (
                end0,
                -1,
                pante
            )
        )


    if not events:
        return


    events.sort(
        key=lambda x: x[0]
    )


    active = Counter()

    prev = None
    i = 0

    state = state_from_chr(
        group,
        chrom
    )


    with open(outfile, "w") as out:

        while i < len(events):

            pos = events[i][0]


            # region between prev and current position
            if (
                prev is not None
                and
                pos > prev
                and
                len(active) == 1
            ):

                pante = next(
                    iter(active.keys())
                )

                out.write(
                    "{}\t{}\t{}\t{}\t{}\n".format(
                        chrom,
                        prev,
                        pos,
                        pante,
                        state
                    )
                )

                bp = pos - prev

                bp_stat[
                    (
                        state,
                        pante
                    )
                ] += bp

                seg_stat[
                    (
                        state,
                        pante
                    )
                ] += 1


            # apply all events at this coordinate
            while (
                i < len(events)
                and
                events[i][0] == pos
            ):

                _, delta, pante = events[i]

                active[pante] += delta

                if active[pante] <= 0:
                    del active[pante]

                i += 1


            prev = pos


# ============================================================
# Process each genome
# ============================================================

bp_stat = Counter()
seg_stat = Counter()


for group in [
    "BSB",
    "TC",
    "BT"
]:

    infile = INPUTS[group]

    print("")
    print(
        "Processing {}".format(
            group
        )
    )


    # input was chromosome-sorted; process one chromosome at a time
    current_chr = None
    intervals = []

    produced = {}


    def flush_chr(chrom, intervals):

        if chrom is None:
            return

        if not is_main_chr(
            group,
            chrom
        ):
            return

        tmp = os.path.join(
            TMPDIR,
            "{}.{}.bed".format(
                group,
                chrom
            )
        )

        process_chr(
            group,
            chrom,
            intervals,
            tmp,
            bp_stat,
            seg_stat
        )

        produced[chrom] = tmp


    with open(infile) as fh:

        fh.readline()

        for line in fh:

            p = line.rstrip("\n").split("\t")

            chrom = p[0]
            start = int(p[1])
            end = int(p[2])
            pante = p[3]

            if current_chr is None:
                current_chr = chrom

            if chrom != current_chr:

                flush_chr(
                    current_chr,
                    intervals
                )

                current_chr = chrom
                intervals = []

            intervals.append(
                (
                    start,
                    end,
                    pante
                )
            )


        flush_chr(
            current_chr,
            intervals
        )


    # --------------------------------------------------------
    # concatenate in BAM reference order
    # --------------------------------------------------------

    outfile = os.path.join(
        OUTDIR,
        "{}.PanTE.assignable_segments.bed".format(
            group
        )
    )


    with open(outfile, "w") as out:

        for chrom in order[group]:

            if not is_main_chr(
                group,
                chrom
            ):
                continue

            tmp = produced.get(
                chrom
            )

            if tmp is None:
                continue

            with open(tmp) as fh:

                for line in fh:
                    out.write(line)


    print(
        "Written:",
        outfile
    )


# ============================================================
# bp summary
# ============================================================

SUMMARY = os.path.join(
    OUTDIR,
    "PanTE_assignable_bp.tsv"
)

with open(SUMMARY, "w") as out:

    out.write(
        "Genome_state\tPanTE_ID\tAssignable_bp\tN_segments\n"
    )

    for key in sorted(bp_stat):

        state, pante = key

        out.write(
            "{}\t{}\t{}\t{}\n".format(
                state,
                pante,
                bp_stat[key],
                seg_stat[key]
            )
        )


print("")
print("Assignable bp summary:")
print(SUMMARY)

shutil.rmtree(TMPDIR)

print("")
print("Completed successfully.")