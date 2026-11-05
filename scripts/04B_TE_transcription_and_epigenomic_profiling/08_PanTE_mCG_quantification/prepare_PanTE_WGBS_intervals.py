#!/usr/bin/env python
# Public-release note: original executed source; replace local absolute paths with local paths before running.
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import re
import sys
from collections import defaultdict


BASE = os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52")

OUTDIR = os.path.join(
    BASE,
    "06_TE_WGBS/01_annotation"
)

os.makedirs(
    OUTDIR,
    exist_ok=True
)


GTFS = {
    "BSB": os.path.join(
        BASE,
        "05_TE_RNA/01_annotation/BSB.PanTE.TEtranscripts.gtf"
    ),

    "TC": os.path.join(
        BASE,
        "05_TE_RNA/01_annotation/TC.PanTE.TEtranscripts.gtf"
    ),

    "BT": os.path.join(
        BASE,
        "05_TE_RNA/01_annotation/BT.PanTE.TEtranscripts.gtf"
    )
}


PANTE_RE = re.compile(
    r"(PanTE_[A-Za-z0-9_]+)"
)


def classify_haplotype(genome, chrom):

    if genome == "BSB":
        return "BSB"

    if genome == "TC":
        return "TC"

    if genome == "BT":

        if re.match(
            r"^chr\d+M$",
            chrom
        ):
            return "BTM"

        if re.match(
            r"^chr\d+P$",
            chrom
        ):
            return "BTP"

        return "BT_other"

    return "NA"


def parse_gtf(genome, filename):

    data = defaultdict(
        list
    )

    n_line = 0
    n_pante = 0

    with open(filename, "r") as fh:

        for line in fh:

            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            n_line += 1

            parts = line.rstrip(
                "\n"
            ).split(
                "\t"
            )

            if len(parts) < 9:
                continue

            chrom = parts[0]

            try:
                start = int(parts[3])
                end = int(parts[4])

            except ValueError:
                continue

            m = PANTE_RE.search(
                line
            )

            if m is None:
                continue

            pante = m.group(1)

            hap = classify_haplotype(
                genome,
                chrom
            )

            # Main chromosomes only
            if genome == "BT":

                if hap not in (
                    "BTM",
                    "BTP"
                ):
                    continue

            else:

                if not re.match(
                    r"^chr\d+$",
                    chrom
                ):
                    continue

            data[
                (
                    chrom,
                    pante,
                    hap
                )
            ].append(
                (
                    start,
                    end
                )
            )

            n_pante += 1

    print(
        "{}: parsed GTF lines={}, PanTE records={}".format(
            genome,
            n_line,
            n_pante
        )
    )

    return data


def merge_intervals(intervals):

    intervals = sorted(
        intervals
    )

    merged = []

    for start, end in intervals:

        if not merged:

            merged.append(
                [
                    start,
                    end
                ]
            )

            continue

        last_start, last_end = merged[-1]

        if start <= last_end + 1:

            if end > last_end:
                merged[-1][1] = end

        else:

            merged.append(
                [
                    start,
                    end
                ]
            )

    return merged


def write_genome(genome, filename):

    data = parse_gtf(
        genome,
        filename
    )

    rows = []

    for (
        chrom,
        pante,
        hap
    ), intervals in data.items():

        merged = merge_intervals(
            intervals
        )

        for start, end in merged:

            rows.append(
                (
                    chrom,
                    start,
                    end,
                    pante,
                    hap
                )
            )

    # chromosome lexical sorting is fine here because names are zero-padded
    rows.sort(
        key=lambda x: (
            x[0],
            x[1],
            x[2],
            x[3]
        )
    )

    outfile = os.path.join(
        OUTDIR,
        "{}.PanTE.WGBS_intervals.tsv".format(
            genome
        )
    )

    with open(
        outfile,
        "w"
    ) as out:

        out.write(
            "Chrom\tStart\tEnd\tPanTE_ID\tGenome_state\n"
        )

        for row in rows:

            out.write(
                "{}\t{}\t{}\t{}\t{}\n".format(
                    *row
                )
            )

    n_pante = len(
        set(
            r[3]
            for r in rows
        )
    )

    print(
        "{}: merged intervals={}, PanTE_IDs={}".format(
            genome,
            len(rows),
            n_pante
        )
    )

    print(
        "  {}".format(
            outfile
        )
    )


def main():

    for genome in [
        "BSB",
        "TC",
        "BT"
    ]:

        f = GTFS[genome]

        if not os.path.exists(f):

            raise RuntimeError(
                "Missing GTF: {}".format(
                    f
                )
            )

        write_genome(
            genome,
            f
        )


if __name__ == "__main__":
    main()