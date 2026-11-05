#!/usr/bin/env python
# Public-release note: original executed source; replace local absolute paths with local paths before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import argparse
import gzip
import heapq
import os
import re
from collections import Counter, defaultdict


def load_intervals(filename):

    by_chr = defaultdict(
        list
    )

    pante_state = {}

    with open(
        filename,
        "r"
    ) as fh:

        header = fh.readline()

        for line in fh:

            if not line.strip():
                continue

            chrom, start, end, pante, state = (
                line.rstrip("\n").split("\t")
            )

            start = int(start)
            end = int(end)

            by_chr[
                chrom
            ].append(
                (
                    start,
                    end,
                    pante,
                    state
                )
            )

            pante_state[
                (
                    pante,
                    state
                )
            ] = 1

    for chrom in by_chr:

        by_chr[
            chrom
        ].sort(
            key=lambda x: (
                x[0],
                x[1],
                x[2]
            )
        )

    return by_chr


def process_cov(
    covfile,
    intervals,
    sample,
    group,
    tissue,
    outfile,
    qcfile
):

    # key:
    # (PanTE_ID, Genome_state)
    #
    # values:
    # N_CpG_sites
    # mC
    # uC

    stat = defaultdict(
        lambda: [
            0,
            0,
            0
        ]
    )


    n_cov = 0
    n_on_te = 0
    n_unique = 0
    n_ambiguous = 0
    n_no_interval = 0


    current_chr = None

    chr_intervals = []

    interval_idx = 0

    active_heap = []

    active_counter = Counter()

    serial = 0


    with gzip.open(
        covfile,
        "rt"
    ) as fh:

        for line in fh:

            if not line.strip():
                continue

            parts = line.rstrip(
                "\n"
            ).split(
                "\t"
            )

            if len(parts) < 6:
                continue

            chrom = parts[0]

            try:

                pos = int(
                    parts[1]
                )

                meth = int(
                    parts[4]
                )

                unmeth = int(
                    parts[5]
                )

            except ValueError:

                continue


            n_cov += 1


            # New chromosome
            if chrom != current_chr:

                current_chr = chrom

                chr_intervals = intervals.get(
                    chrom,
                    []
                )

                interval_idx = 0

                active_heap = []

                active_counter = Counter()

                serial = 0


            if not chr_intervals:

                n_no_interval += 1
                continue


            # Add intervals that have started
            while (
                interval_idx
                <
                len(chr_intervals)
                and
                chr_intervals[
                    interval_idx
                ][0]
                <=
                pos
            ):

                start, end, pante, state = (
                    chr_intervals[
                        interval_idx
                    ]
                )

                serial += 1

                heapq.heappush(
                    active_heap,
                    (
                        end,
                        serial,
                        pante,
                        state
                    )
                )

                active_counter[
                    (
                        pante,
                        state
                    )
                ] += 1

                interval_idx += 1


            # Remove intervals that ended before this CpG
            while (
                active_heap
                and
                active_heap[0][0]
                <
                pos
            ):

                end, sid, pante, state = (
                    heapq.heappop(
                        active_heap
                    )
                )

                key = (
                    pante,
                    state
                )

                active_counter[
                    key
                ] -= 1

                if active_counter[
                    key
                ] <= 0:

                    del active_counter[
                        key
                    ]


            if not active_counter:

                n_no_interval += 1
                continue


            n_on_te += 1


            # Unique PanTE assignment.
            #
            # Multiple overlapping intervals belonging to the SAME
            # PanTE_ID/state count as one assignment.
            #
            # If different PanTE_IDs overlap the CpG, discard it.
            keys = list(
                active_counter.keys()
            )


            unique_pante_ids = set(
                x[0]
                for x in keys
            )


            if len(
                unique_pante_ids
            ) != 1:

                n_ambiguous += 1
                continue


            # In BT, M/P chromosome identity guarantees a single state.
            # For BSB/TC there is likewise one state.
            #
            # Take the unique key.
            if len(keys) > 1:

                states = set(
                    x[1]
                    for x in keys
                )

                if len(states) != 1:

                    n_ambiguous += 1
                    continue


            key = keys[0]


            stat[
                key
            ][0] += 1

            stat[
                key
            ][1] += meth

            stat[
                key
            ][2] += unmeth

            n_unique += 1


    # --------------------------------------------------------
    # Output
    # --------------------------------------------------------

    rows = []

    for (
        pante,
        state
    ), values in stat.items():

        n_cpg, mc, uc = values

        total = mc + uc

        if total > 0:

            weighted = (
                float(mc)
                /
                float(total)
            )

        else:

            weighted = float(
                "nan"
            )


        rows.append(
            (
                sample,
                group,
                tissue,
                state,
                pante,
                n_cpg,
                mc,
                uc,
                total,
                weighted
            )
        )


    rows.sort(
        key=lambda x: (
            x[3],
            x[4]
        )
    )


    with open(
        outfile,
        "w"
    ) as out:

        out.write(
            "Sample\tGroup\tTissue\tGenome_state\t"
            "PanTE_ID\tN_CpG_sites\tmC\tuC\t"
            "Total_C_reads\tweighted_mCG\n"
        )

        for row in rows:

            out.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.8f}\n".format(
                    *row
                )
            )


    ambiguous_pct = (
        100.0
        *
        n_ambiguous
        /
        n_on_te
        if n_on_te > 0
        else 0.0
    )


    assigned_pct = (
        100.0
        *
        n_unique
        /
        n_on_te
        if n_on_te > 0
        else 0.0
    )


    with open(
        qcfile,
        "w"
    ) as q:

        q.write(
            "Metric\tValue\n"
        )

        q.write(
            "Sample\t{}\n".format(
                sample
            )
        )

        q.write(
            "Group\t{}\n".format(
                group
            )
        )

        q.write(
            "Tissue\t{}\n".format(
                tissue
            )
        )

        q.write(
            "Total_cov_records\t{}\n".format(
                n_cov
            )
        )

        q.write(
            "CpG_records_overlapping_PanTE\t{}\n".format(
                n_on_te
            )
        )

        q.write(
            "Unique_assigned_CpG_records\t{}\n".format(
                n_unique
            )
        )

        q.write(
            "Ambiguous_CpG_records\t{}\n".format(
                n_ambiguous
            )
        )

        q.write(
            "Unique_assignment_percent_of_TE_CpGs\t{:.6f}\n".format(
                assigned_pct
            )
        )

        q.write(
            "Ambiguous_percent_of_TE_CpGs\t{:.6f}\n".format(
                ambiguous_pct
            )
        )

        q.write(
            "Quantified_PanTE_state_rows\t{}\n".format(
                len(rows)
            )
        )


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--cov",
        required=True
    )

    parser.add_argument(
        "--intervals",
        required=True
    )

    parser.add_argument(
        "--sample",
        required=True
    )

    parser.add_argument(
        "--group",
        required=True
    )

    parser.add_argument(
        "--tissue",
        required=True
    )

    parser.add_argument(
        "--out",
        required=True
    )

    parser.add_argument(
        "--qc",
        required=True
    )

    args = parser.parse_args()


    intervals = load_intervals(
        args.intervals
    )


    process_cov(
        args.cov,
        intervals,
        args.sample,
        args.group,
        args.tissue,
        args.out,
        args.qc
    )


if __name__ == "__main__":
    main()