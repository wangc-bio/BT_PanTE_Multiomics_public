#!/usr/bin/env python3

import os
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="Split the unified curated TE library into broad-class FASTA files.")
parser.add_argument("--input", required=True, help="Unified curated repeat FASTA")
parser.add_argument("--output-dir", required=True, help="Output directory for class-specific FASTA files")
args = parser.parse_args()
FA = args.input
OUT = args.output_dir
os.makedirs(OUT, exist_ok=True)
if not os.path.isfile(FA):
    raise FileNotFoundError(FA)

KEEP = {"DNA", "LTR", "LINE", "SINE", "RC"}


def get_broad(cls):

    s = cls.lower()

    if s.startswith("dna"):
        return "DNA"

    if s.startswith("ltr"):
        return "LTR"

    if s.startswith("line"):
        return "LINE"

    if s.startswith("sine"):
        return "SINE"

    if s.startswith("rc") or "helitron" in s:
        return "RC"

    return "OTHER"


records = defaultdict(list)

header = None
seq = []

with open(FA) as fh:

    for line in fh:

        line = line.rstrip()

        if line.startswith(">"):

            if header is not None:

                token = header.split()[0]

                if "#" in token:
                    uid, cls = token.split("#", 1)
                else:
                    uid = token
                    cls = "Unknown"

                b = get_broad(cls)

                if b in KEEP:
                    records[b].append(
                        (uid, cls, "".join(seq))
                    )

            header = line[1:]
            seq = []

        else:
            seq.append(line.strip())


if header is not None:

    token = header.split()[0]

    if "#" in token:
        uid, cls = token.split("#", 1)
    else:
        uid = token
        cls = "Unknown"

    b = get_broad(cls)

    if b in KEEP:
        records[b].append(
            (uid, cls, "".join(seq))
        )


for b in ["DNA", "LTR", "LINE", "SINE", "RC"]:

    outfile = os.path.join(
        OUT,
        "%s.fa" % b
    )

    with open(outfile, "w") as fo:

        for uid, cls, s in records[b]:

            fo.write(
                ">%s#%s\n%s\n"
                % (
                    uid,
                    cls,
                    s
                )
            )

    print(
        "%s\t%d"
        % (
            b,
            len(records[b])
        )
    )