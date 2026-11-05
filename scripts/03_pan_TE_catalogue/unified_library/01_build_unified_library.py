#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build a unified curated repeat library from one or more genome-specific libraries.

Method retained from the supplied v5.2 production script:
  1. Read all input FASTA libraries listed in a TSV manifest.
  2. Uppercase sequences and deduplicate by exact-sequence MD5 only.
  3. Reconcile classifications for exact duplicate sequences.
  4. Write the unified FASTA, sequence provenance map, summary and conflict table.

No sample names, server paths or production input files are hard-coded.
"""

from __future__ import print_function

import argparse
import csv
import hashlib
import os
import re
import sys
from collections import defaultdict


def read_fasta(path):
    header = None
    seq = []

    with open(path, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")

            if not line.strip():
                continue

            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(seq)
                    if not sequence:
                        raise ValueError(
                            "Empty FASTA sequence before line %d in %s"
                            % (line_number, path)
                        )
                    yield header, sequence

                header = line[1:].strip()
                if not header:
                    raise ValueError(
                        "Empty FASTA header at line %d in %s"
                        % (line_number, path)
                    )
                seq = []
            else:
                if header is None:
                    raise ValueError(
                        "Sequence data before first FASTA header in %s"
                        % path
                    )
                seq.append("".join(line.split()))

    if header is not None:
        sequence = "".join(seq)
        if not sequence:
            raise ValueError("Empty final FASTA sequence in %s" % path)
        yield header, sequence


def parse_header(header):
    token = header.split()[0]

    if "#" in token:
        record_id, classification = token.split("#", 1)
    else:
        record_id = token
        classification = "Unknown"

    return record_id, classification


def broad(classification):
    value = (classification or "").strip().lower()

    if not value:
        return "Unknown"

    # Non-TE annotations
    if (
        "simple_repeat" in value
        or "low_complexity" in value
        or "satellite" in value
        or value in (
            "trna",
            "rrna",
            "scrna",
            "snrna",
            "srprna",
            "ncrna",
        )
    ):
        return "nonTE"

    # RC / Helitron
    if value.startswith("rc") or "helitron" in value:
        return "RC"

    # SINE
    if value.startswith("sine") or "sine" in value:
        return "SINE"

    # LTR
    if (
        value.startswith("ltr")
        or "gypsy" in value
        or "copia" in value
        or "bel" in value
        or "pao" in value
        or "erv" in value
        or "dirs" in value
    ):
        return "LTR"

    # LINE
    if (
        value.startswith("line")
        or "non-ltr" in value
        or "cr1" in value
        or "rex" in value
        or "jockey" in value
        or "rte" in value
        or "penelope" in value
    ):
        return "LINE"

    # DNA transposons
    if (
        value.startswith("dna")
        or "tc1" in value
        or "mariner" in value
        or "hat" in value
        or "harbinger" in value
        or "piggybac" in value
        or "mudr" in value
        or "mule" in value
        or "transib" in value
        or "kolobok" in value
        or "cacta" in value
        or "crypton" in value
        or "sola" in value
        or "tir" in value
    ):
        return "DNA"

    if "unknown" in value:
        return "Unknown"

    return "Other"


def broad_fallback(broad_class):
    if broad_class == "DNA":
        return "DNA/Unknown"
    if broad_class == "LTR":
        return "LTR/Unknown"
    if broad_class == "LINE":
        return "LINE/Unknown"
    if broad_class == "SINE":
        return "SINE/Unknown"
    if broad_class == "RC":
        return "RC/Helitron"
    if broad_class == "nonTE":
        return "nonTE"
    return "Unknown"


def read_manifest(manifest_path):
    manifest_dir = os.path.dirname(os.path.abspath(manifest_path))
    entries = []
    seen_source_ids = set()

    with open(manifest_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"source_id", "library_fasta"}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                "Manifest must contain tab-separated columns: "
                "source_id and library_fasta"
            )

        for line_number, row in enumerate(reader, 2):
            source_id = (row.get("source_id") or "").strip()
            library_fasta = (row.get("library_fasta") or "").strip()

            if not source_id or not library_fasta:
                raise ValueError(
                    "Missing source_id or library_fasta at manifest line %d"
                    % line_number
                )

            if not re.match(r"^[A-Za-z0-9_.-]+$", source_id):
                raise ValueError(
                    "source_id must contain only letters, numbers, '.', '_' "
                    "or '-': %s" % source_id
                )

            if source_id in seen_source_ids:
                raise ValueError("Duplicate source_id in manifest: %s" % source_id)
            seen_source_ids.add(source_id)

            if not os.path.isabs(library_fasta):
                library_fasta = os.path.join(manifest_dir, library_fasta)

            library_fasta = os.path.abspath(library_fasta)
            if not os.path.isfile(library_fasta):
                raise FileNotFoundError(
                    "Input library does not exist: %s" % library_fasta
                )

            entries.append((source_id, library_fasta))

    if not entries:
        raise ValueError("The input manifest contains no library entries.")

    return entries


def write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writerow(header)
        writer.writerows(rows)


def build_unified_library(manifest, output_dir, output_prefix, id_prefix):
    entries = read_manifest(manifest)
    os.makedirs(output_dir, exist_ok=True)

    output_fasta = os.path.join(
        output_dir, output_prefix + ".curated_repeat.fa"
    )
    map_file = os.path.join(
        output_dir, output_prefix + ".sequence_map.tsv"
    )
    summary_file = os.path.join(
        output_dir, output_prefix + ".summary.tsv"
    )
    conflict_file = os.path.join(
        output_dir, output_prefix + ".class_conflicts.tsv"
    )

    groups = {}
    input_counts = defaultdict(int)

    for source_id, library_path in entries:
        for header, sequence in read_fasta(library_path):
            record_id, classification = parse_header(header)
            sequence = sequence.upper()

            try:
                sequence_bytes = sequence.encode("ascii")
            except UnicodeEncodeError:
                raise ValueError(
                    "Non-ASCII sequence encountered: %s in %s"
                    % (record_id, library_path)
                )

            sequence_md5 = hashlib.md5(sequence_bytes).hexdigest()
            input_counts[source_id] += 1

            if sequence_md5 not in groups:
                groups[sequence_md5] = {
                    "sequence": sequence,
                    "members": [],
                }

            groups[sequence_md5]["members"].append({
                "source_id": source_id,
                "record_id": record_id,
                "class": classification,
                "broad": broad(classification),
            })

    if not groups:
        raise ValueError("No FASTA records were found in the input libraries.")

    records = []
    conflicts = []
    decision_counts = defaultdict(int)
    broad_counts = defaultdict(int)

    for sequence_md5, data in groups.items():
        members = data["members"]

        classes = sorted(set(member["class"] for member in members))
        broad_classes = sorted(set(member["broad"] for member in members))

        # Exact duplicate records with identical fine-level classification.
        if len(classes) == 1:
            final_class = classes[0]
            final_broad = broad(final_class)
            decision = "SAME_CLASS"

        # Fine-level labels differ, but all records agree on a TE broad class.
        elif (
            len(broad_classes) == 1
            and broad_classes[0] in ("DNA", "LTR", "LINE", "SINE", "RC")
        ):
            final_broad = broad_classes[0]
            final_class = broad_fallback(final_broad)
            decision = "SAME_BROAD_FINE_CONFLICT"

        elif len(broad_classes) == 1 and broad_classes[0] == "nonTE":
            final_broad = "nonTE"
            final_class = "nonTE"
            decision = "NONTE_AGREEMENT"

        elif len(broad_classes) == 1 and broad_classes[0] == "Unknown":
            final_broad = "Unknown"
            final_class = "Unknown"
            decision = "ALL_UNKNOWN"

        # If exactly one known TE classification conflicts only with
        # Unknown, retain the known classification. This is the
        # KNOWN_VS_UNKNOWN_KEEP_KNOWN rule used in the final v5.2 library.
        elif (
            "Unknown" in broad_classes
            and len(set(broad_classes) - {"Unknown"}) == 1
            and len(set(classes) - {"Unknown"}) == 1
            and next(iter(set(broad_classes) - {"Unknown"}))
                in ("DNA", "LTR", "LINE", "SINE", "RC")
        ):
            final_class = next(iter(set(classes) - {"Unknown"}))
            final_broad = broad(final_class)
            decision = "KNOWN_VS_UNKNOWN_KEEP_KNOWN"

            # Preserve the original disagreement in the conflict audit;
            # the sequence map records how it was resolved.
            conflicts.append([
                sequence_md5,
                ";".join(broad_classes),
                ";".join(classes),
                ";".join(
                    "%s:%s#%s" % (
                        member["source_id"],
                        member["record_id"],
                        member["class"],
                    )
                    for member in members
                ),
            ])

        # Other broad-class conflicts remain unresolved and are assigned
        # Unknown; their source records are retained in the conflict audit.
        else:
            final_broad = "Unknown"
            final_class = "Unknown"
            decision = "BROAD_CONFLICT_TO_UNKNOWN"

            conflicts.append([
                sequence_md5,
                ";".join(broad_classes),
                ";".join(classes),
                ";".join(
                    "%s:%s#%s" % (
                        member["source_id"],
                        member["record_id"],
                        member["class"],
                    )
                    for member in members
                ),
            ])

        decision_counts[decision] += 1
        broad_counts[final_broad] += 1

        records.append({
            "md5": sequence_md5,
            "sequence": data["sequence"],
            "members": members,
            "final_class": final_class,
            "final_broad": final_broad,
            "decision": decision,
        })

    # Stable ordering makes IDs reproducible across runs with the same inputs.
    records.sort(key=lambda record: record["md5"])

    with open(output_fasta, "w", encoding="utf-8") as fasta_handle, \
            open(map_file, "w", encoding="utf-8", newline="") as map_handle:
        map_writer = csv.writer(
            map_handle,
            delimiter="\t",
            lineterminator="\n",
        )
        map_writer.writerow([
            "Unified_ID",
            "MD5",
            "Final_class",
            "Final_broad",
            "Length",
            "Source_count",
            "Decision",
            "Source_members",
        ])

        for index, record in enumerate(records, 1):
            unified_id = "%s_%06d" % (id_prefix, index)
            fasta_handle.write(
                ">%s#%s\n" % (unified_id, record["final_class"])
            )

            sequence = record["sequence"]
            for start in range(0, len(sequence), 60):
                fasta_handle.write(sequence[start:start + 60] + "\n")

            map_writer.writerow([
                unified_id,
                record["md5"],
                record["final_class"],
                record["final_broad"],
                len(sequence),
                len(record["members"]),
                record["decision"],
                ";".join(
                    "%s:%s#%s" % (
                        member["source_id"],
                        member["record_id"],
                        member["class"],
                    )
                    for member in record["members"]
                ),
            ])

    write_tsv(
        conflict_file,
        ["MD5", "Broads", "Classes", "Members"],
        conflicts,
    )

    summary_rows = []
    for source_id, _ in entries:
        summary_rows.append(["Input_%s" % source_id, input_counts[source_id]])

    summary_rows.extend([
        ["Input_total", sum(input_counts.values())],
        ["Unified_unique_sequences", len(records)],
        [
            "Exact_duplicate_groups",
            sum(1 for record in records if len(record["members"]) > 1),
        ],
        [
            "Known_vs_Unknown_keep_known",
            decision_counts["KNOWN_VS_UNKNOWN_KEEP_KNOWN"],
        ],
        [
            "Broad_conflicts_to_Unknown",
            decision_counts["BROAD_CONFLICT_TO_UNKNOWN"],
        ],
        [
            "Same_broad_fine_conflicts",
            decision_counts["SAME_BROAD_FINE_CONFLICT"],
        ],
        ["", ""],
        ["Final_broad", "Families"],
    ])

    for broad_class in (
        "DNA", "LTR", "LINE", "SINE", "RC", "Unknown", "nonTE", "Other"
    ):
        summary_rows.append([broad_class, broad_counts[broad_class]])

    write_tsv(summary_file, ["Metric", "Count"], summary_rows)

    print("Unified curated repeat library")
    print("========================================")
    for source_id, _ in entries:
        print("%-20s : %d" % (source_id, input_counts[source_id]))
    print("")
    print("Input total              = %d" % sum(input_counts.values()))
    print("Unified unique sequences = %d" % len(records))
    print(
        "Broad conflicts          = %d"
        % decision_counts["BROAD_CONFLICT_TO_UNKNOWN"]
    )
    print("")
    print("Library :", output_fasta)
    print("Map     :", map_file)
    print("Summary :", summary_file)
    print("Conflict:", conflict_file)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build a unified curated repeat library using exact-sequence "
            "MD5 deduplication and classification-conflict reconciliation."
        )
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="TSV input manifest with columns: source_id and library_fasta",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for the unified FASTA and result tables",
    )
    parser.add_argument(
        "--output-prefix",
        default="Unified_v5_2",
        help="Output filename prefix (default: Unified_v5_2)",
    )
    parser.add_argument(
        "--id-prefix",
        default="UNI_v52",
        help="Unified FASTA sequence ID prefix (default: UNI_v52)",
    )
    args = parser.parse_args()

    build_unified_library(
        manifest=args.manifest,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        id_prefix=args.id_prefix,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("ERROR: %s" % error, file=sys.stderr)
        sys.exit(1)
