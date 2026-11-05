#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Build a genome-normalized Pan-TE abundance matrix from per-assembly tables."""

import argparse
import csv
import os
import re
import sys


def read_manifest(path):
    manifest_dir = os.path.dirname(os.path.abspath(path))
    entries = []
    seen = set()

    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "genome_size_bp", "abundance_tsv"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                "Manifest must contain tab-separated columns: "
                "sample_id, genome_size_bp, abundance_tsv"
            )

        for line_number, row in enumerate(reader, 2):
            sample = (row.get("sample_id") or "").strip()
            size_text = (row.get("genome_size_bp") or "").strip()
            abundance_path = (row.get("abundance_tsv") or "").strip()

            if not sample or not size_text or not abundance_path:
                raise ValueError(
                    "Missing sample_id, genome_size_bp or abundance_tsv "
                    "at manifest line %d" % line_number
                )
            if not re.match(r"^[A-Za-z0-9_.-]+$", sample):
                raise ValueError("Invalid sample_id: %s" % sample)
            if sample in seen:
                raise ValueError("Duplicate sample_id: %s" % sample)
            seen.add(sample)

            genome_size = int(size_text)
            if genome_size <= 0:
                raise ValueError(
                    "genome_size_bp must be positive for %s" % sample
                )

            if not os.path.isabs(abundance_path):
                abundance_path = os.path.join(manifest_dir, abundance_path)
            abundance_path = os.path.abspath(abundance_path)
            if not os.path.isfile(abundance_path):
                raise FileNotFoundError(
                    "Abundance table not found: %s" % abundance_path
                )

            entries.append({
                "sample_id": sample,
                "genome_size_bp": genome_size,
                "abundance_tsv": abundance_path,
            })

    if not entries:
        raise ValueError("The manifest contains no sample entries.")
    return entries


def read_cluster_metadata(path):
    metadata = {}

    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "PanTE_ID", "Broad_class", "Cluster_size",
            "Representative_ID", "Unified_ID"
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                "Pan-TE membership table must contain: %s"
                % ", ".join(sorted(required))
            )

        for row in reader:
            pid = row["PanTE_ID"]
            record = {
                "Broad_class": row["Broad_class"],
                "Cluster_size": int(row["Cluster_size"]),
                "Representative_ID": row["Representative_ID"],
            }

            if pid in metadata and metadata[pid] != record:
                raise ValueError(
                    "Inconsistent metadata within Pan-TE cluster %s" % pid
                )
            metadata[pid] = record

    if not metadata:
        raise ValueError("No Pan-TE clusters found in membership table.")
    return metadata


def read_abundance(path, sample_id, known_clusters):
    abundance = {}
    segments = {}

    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "PanTE_ID", "Occupied_bp", "Resolved_segments", "Sample"
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                "Abundance table %s must contain: %s"
                % (path, ", ".join(sorted(required)))
            )

        for row in reader:
            pid = row["PanTE_ID"]
            if row["Sample"] != sample_id:
                raise ValueError(
                    "Sample label mismatch in %s: expected %s, found %s"
                    % (path, sample_id, row["Sample"])
                )
            if pid not in known_clusters:
                raise ValueError(
                    "Pan-TE ID %s in %s is absent from the membership table"
                    % (pid, path)
                )
            if pid in abundance:
                raise ValueError(
                    "Duplicate Pan-TE ID %s in %s" % (pid, path)
                )

            abundance[pid] = int(row["Occupied_bp"])
            segments[pid] = int(row["Resolved_segments"])

    return abundance, segments


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Merge per-assembly Pan-TE resolved-abundance tables and "
            "normalize occupied base pairs by assembly size."
        )
    )
    parser.add_argument("--members", required=True,
                        help="PanTE_refined_cluster_members.tsv")
    parser.add_argument("--manifest", required=True,
                        help="TSV with sample_id, genome_size_bp, abundance_tsv")
    parser.add_argument("--out", required=True,
                        help="Output Pan-TE abundance matrix TSV")
    parser.add_argument("--presence-bp", type=int, default=10000,
                        help="Occupied-bp threshold for summary counts (default: 10000)")
    args = parser.parse_args()

    if args.presence_bp < 0:
        raise ValueError("--presence-bp must be non-negative")

    samples = read_manifest(args.manifest)
    metadata = read_cluster_metadata(args.members)
    known_clusters = set(metadata)

    abundance_by_sample = {}
    segments_by_sample = {}

    for entry in samples:
        sample = entry["sample_id"]
        abundance, segments = read_abundance(
            entry["abundance_tsv"], sample, known_clusters
        )
        abundance_by_sample[sample] = abundance
        segments_by_sample[sample] = segments

    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)
    summary_path = os.path.splitext(args.out)[0] + ".summary.tsv"

    fields = [
        "PanTE_ID", "Broad_class", "Cluster_size", "Representative_ID"
    ]
    for entry in samples:
        sample = entry["sample_id"]
        fields.extend([
            "%s_bp" % sample,
            "%s_bp_per_Mb" % sample,
            "%s_segments" % sample,
        ])

    with open(args.out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t",
            lineterminator="\n"
        )
        writer.writeheader()

        for pid in sorted(metadata):
            row = {
                "PanTE_ID": pid,
                "Broad_class": metadata[pid]["Broad_class"],
                "Cluster_size": metadata[pid]["Cluster_size"],
                "Representative_ID": metadata[pid]["Representative_ID"],
            }

            for entry in samples:
                sample = entry["sample_id"]
                genome_size = entry["genome_size_bp"]
                bp = abundance_by_sample[sample].get(pid, 0)
                n_segments = segments_by_sample[sample].get(pid, 0)
                row["%s_bp" % sample] = bp
                row["%s_bp_per_Mb" % sample] = "%.6f" % (
                    bp / genome_size * 1000000.0
                )
                row["%s_segments" % sample] = n_segments

            writer.writerow(row)

    with open(summary_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(
            handle, delimiter="\t", lineterminator="\n"
        )
        writer.writerow([
            "Sample", "Genome_size_bp", "Occupied_bp",
            "Detected_PanTE", "PanTE_ge_presence_threshold_bp"
        ])

        for entry in samples:
            sample = entry["sample_id"]
            abundance = abundance_by_sample[sample]
            total_bp = sum(abundance.values())
            detected = sum(1 for bp in abundance.values() if bp > 0)
            above_threshold = sum(
                1 for bp in abundance.values() if bp >= args.presence_bp
            )
            writer.writerow([
                sample,
                entry["genome_size_bp"],
                total_bp,
                detected,
                above_threshold,
            ])

    print("Pan-TE clusters:", len(metadata))
    print("Assemblies:", len(samples))
    print("Matrix:", args.out)
    print("Summary:", summary_path)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("ERROR: %s" % error, file=sys.stderr)
        sys.exit(1)
