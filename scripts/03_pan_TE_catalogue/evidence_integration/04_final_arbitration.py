#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Strict evidence-based TE classification arbitration.

Generic release derived from the historical final_arbitration_v5.py.
The three independent evidence dimensions are Homology, Domain and Structure.
Prior/composite labels are used only to refine the fine class after broad-class
assignment from external evidence.
"""
from __future__ import print_function
import argparse
import csv
import os
from collections import defaultdict

TE_BROAD = {"DNA", "LTR", "LINE", "SINE", "RC"}
CLASSIFIABLE_CONTEXTS = {
    "TWO_STRONG_PLUS_SUPPORT", "TWO_STRONG",
    "ONE_STRONG_PLUS_2SUPPORT", "ONE_STRONG_PLUS_SUPPORT",
    "ONE_STRONG_ONLY", "MULTI_SUPPORT"
}

def parse_args():
    p = argparse.ArgumentParser(description="Strict evidence arbitration for TE evidence-master tables.")
    p.add_argument("--manifest", required=True, help="TSV with sample_id and evidence_master columns.")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--review-dir", required=True)
    return p.parse_args()

def resolve_path(value, base):
    value = (value or "").strip()
    if not value:
        raise ValueError("A required manifest path is blank.")
    return os.path.abspath(value if os.path.isabs(value) else os.path.join(base, value))

def read_manifest(path):
    base = os.path.dirname(os.path.abspath(path))
    out, seen = [], set()
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"sample_id", "evidence_master"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Manifest missing columns: " + ", ".join(sorted(missing)))
        for line_no, row in enumerate(reader, 2):
            sid = (row.get("sample_id") or "").strip()
            if not sid or sid in seen:
                raise ValueError("Blank or duplicate sample_id at manifest line %d" % line_no)
            seen.add(sid)
            path_in = resolve_path(row["evidence_master"], base)
            if not os.path.isfile(path_in):
                raise FileNotFoundError("Evidence master not found: " + path_in)
            out.append({"sample_id": sid, "evidence_master": path_in})
    if not out:
        raise ValueError("Manifest contains no sample rows.")
    return out

def broad_from_fine(value):
    s = (value or "").strip().lower()
    if not s:
        return "Unknown"
    if any(x in s for x in ("simple_repeat", "simple repeat", "satellite",
                            "low_complexity", "low complexity")) or s in ("trna", "rrna"):
        return "nonTE"
    if "sine" in s:
        return "SINE"
    if "helitron" in s or s.startswith("rc/"):
        return "RC"
    if any(x in s for x in ("ltr", "gypsy", "copia", "bel", "pao", "erv", "dirs")):
        return "LTR"
    if any(x in s for x in ("line", "non-ltr", "nonltr", "cr1", "rex", "rte", "jockey", "penelope")):
        return "LINE"
    if any(x in s for x in ("dna", "tir", "mite", "maverick", "polinton")):
        return "DNA"
    return "Unknown"

def tier(cls, strength):
    if not cls or cls in ("Unknown", "CONFLICT"):
        return "NONE"
    strength = strength or ""
    if strength.startswith("STRONG"):
        return "STRONG"
    if strength.startswith("SUPPORTIVE"):
        return "SUPPORTIVE"
    return "NONE"

def evidence_context(row):
    dims = [
        ("H", row.get("Homology_class", ""), row.get("Homology_strength", ""), row.get("Homology_conflict", "")),
        ("D", row.get("Domain_class", ""), row.get("Domain_strength", ""), row.get("Domain_conflict", "")),
        ("S", row.get("Structure_class", ""), row.get("Structure_strength", ""), row.get("Structure_conflict", ""))
    ]
    internal_review = any(conflict == "YES" for _, _, _, conflict in dims)
    strong, supportive = [], []
    for name, cls, strength, _ in dims:
        t = tier(cls, strength)
        if t == "STRONG":
            strong.append((name, cls))
        elif t == "SUPPORTIVE":
            supportive.append((name, cls))
    strong_classes = sorted(set(cls for _, cls in strong))
    support_classes = sorted(set(cls for _, cls in supportive))
    if internal_review:
        return "INTERNAL_DIMENSION_REVIEW", ""
    if len(strong_classes) > 1:
        return "STRONG_CLASS_CONFLICT", ""
    if len(strong_classes) == 1:
        candidate = strong_classes[0]
        n_strong_same = sum(1 for _, cls in strong if cls == candidate)
        n_support_same = sum(1 for _, cls in supportive if cls == candidate)
        n_support_opposite = sum(1 for _, cls in supportive if cls != candidate)
        if n_support_opposite > 0:
            return "STRONG_VS_SUPPORT_CONFLICT", candidate
        if n_strong_same >= 2:
            return ("TWO_STRONG_PLUS_SUPPORT" if n_support_same >= 1 else "TWO_STRONG"), candidate
        if n_support_same >= 2:
            return "ONE_STRONG_PLUS_2SUPPORT", candidate
        if n_support_same == 1:
            return "ONE_STRONG_PLUS_SUPPORT", candidate
        return "ONE_STRONG_ONLY", candidate
    if len(support_classes) > 1:
        return "SUPPORT_CLASS_CONFLICT", ""
    if len(support_classes) == 1:
        candidate = support_classes[0]
        n_support_same = sum(1 for _, cls in supportive if cls == candidate)
        return ("MULTI_SUPPORT" if n_support_same >= 2 else "SUPPORT_ONLY"), candidate
    return "NO_EXTERNAL_EVIDENCE", ""

def confidence_from_context(context):
    return {
        "TWO_STRONG_PLUS_SUPPORT": "VERY_HIGH",
        "TWO_STRONG": "HIGH",
        "ONE_STRONG_PLUS_2SUPPORT": "MEDIUM_HIGH",
        "ONE_STRONG_PLUS_SUPPORT": "MEDIUM_HIGH",
        "ONE_STRONG_ONLY": "MEDIUM",
        "MULTI_SUPPORT": "MEDIUM",
        "SUPPORT_ONLY": "LOW",
        "STRONG_CLASS_CONFLICT": "REVIEW",
        "STRONG_VS_SUPPORT_CONFLICT": "REVIEW",
        "SUPPORT_CLASS_CONFLICT": "REVIEW",
        "INTERNAL_DIMENSION_REVIEW": "REVIEW",
        "NO_EXTERNAL_EVIDENCE": "PRIOR_ONLY"
    }[context]

def default_fine(broad):
    return {"DNA": "DNA/Unknown", "LTR": "LTR/Unknown", "LINE": "LINE/Unknown",
            "SINE": "SINE/Unknown", "RC": "RC/Helitron", "nonTE": "nonTE",
            "Unknown": "Unknown"}.get(broad, "Unknown")

def choose_fine_class(row, final_broad):
    if final_broad in ("", "Unknown"):
        return "Unknown", "NONE"
    if final_broad == "nonTE":
        original = row.get("Original_class", "")
        if broad_from_fine(original) == "nonTE":
            return original, "ORIGINAL_NONTE"
        return "nonTE", "EXTERNAL_NONTE"
    original = row.get("Original_class", "")
    if original and broad_from_fine(original) == final_broad and "unknown" not in original.lower():
        return original, "ORIGINAL_COMPATIBLE_FINE"
    ltr_target = row.get("LTR_target_class", "")
    if final_broad == "LTR" and ltr_target and broad_from_fine(ltr_target) == "LTR" and "unknown" not in ltr_target.lower():
        return ltr_target, "LTRRETRIEVER_FINE"
    mch = row.get("MCHelper_class", "")
    if mch and broad_from_fine(mch) == final_broad and "unknown" not in mch.lower():
        return mch, "MCHELPER_COMPATIBLE_FINE"
    tet = row.get("TEtrimmer_class", "")
    if tet and broad_from_fine(tet) == final_broad and "unknown" not in tet.lower():
        return tet, "TETRIMMER_COMPATIBLE_FINE"
    return default_fine(final_broad), "BROAD_ONLY"

def classify(row):
    context, candidate = evidence_context(row)
    confidence = confidence_from_context(context)
    original_broad = row.get("Original_broad", "Unknown") or "Unknown"
    final_broad, action, reason = "Unknown", "", ""
    if original_broad == "nonTE":
        if candidate not in ("", "Unknown", "nonTE") and context in ("TWO_STRONG_PLUS_SUPPORT", "TWO_STRONG"):
            confidence = "REVIEW"
            action = "REVIEW_NONTE_VS_STRONG_TE"
            reason = "Original family is explicit nonTE, but at least two strong external dimensions support a TE class."
        else:
            final_broad = "nonTE"
            confidence = "NONTE"
            action = "KEEP_NONTE"
            reason = "Explicit original nonTE family retained; no two-strong external TE contradiction."
    elif context in CLASSIFIABLE_CONTEXTS:
        if candidate == "nonTE":
            confidence = "REVIEW"
            action = "REVIEW_EXTERNAL_NONTE"
            reason = "External evidence suggests nonTE, but the original family is not explicit nonTE; kept Unknown."
        elif candidate:
            final_broad = candidate
            action = "CLASSIFY_EXTERNAL"
            reason = "Broad class assigned from qualifying independent external evidence."
        else:
            confidence = "REVIEW"
            action = "REVIEW_EMPTY_CANDIDATE"
            reason = "Classifiable evidence context did not provide a broad-class candidate."
    elif context == "SUPPORT_ONLY":
        action = "KEEP_UNKNOWN_SUPPORT_ONLY"
        reason = "Only one supportive external dimension; insufficient for formal broad-class rescue."
    elif context == "NO_EXTERNAL_EVIDENCE":
        action = "KEEP_UNKNOWN_PRIOR_ONLY"
        reason = "No qualifying independent external evidence; prior/composite labels are not used alone for broad-class assignment."
    else:
        confidence = "REVIEW"
        action = "REVIEW_CONFLICT"
        reason = "External evidence dimensions contain class or internal-source conflict."
    final_class, fine_source = choose_fine_class(row, final_broad)
    out = dict(row)
    out.update({
        "Final_evidence_context": context,
        "Final_candidate_broad": candidate,
        "Final_broad_class": final_broad,
        "Final_class": final_class,
        "Final_confidence": confidence,
        "Fine_class_source": fine_source,
        "Final_action": action,
        "Final_reason": reason
    })
    return out

def write_tsv(path, rows, fields=None):
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

def main():
    args = parse_args()
    entries = read_manifest(args.manifest)
    outdir, reviewdir = os.path.abspath(args.output_dir), os.path.abspath(args.review_dir)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(reviewdir, exist_ok=True)
    summaries, class_rows = [], []
    for entry in entries:
        sample = entry["sample_id"]
        with open(entry["evidence_master"], newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            required = {"RepeatModeler_ID","Original_broad","Original_class",
                        "Homology_class","Homology_strength","Homology_conflict",
                        "Domain_class","Domain_strength","Domain_conflict",
                        "Structure_class","Structure_strength","Structure_conflict"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError("%s evidence master missing columns: %s" % (sample, ", ".join(sorted(missing))))
            rows = [classify(row) for row in reader]
        if not rows:
            raise RuntimeError("%s evidence master contains zero rows." % sample)
        outpath = os.path.join(outdir, "%s.v5_final_classification.tsv" % sample)
        reviewpath = os.path.join(reviewdir, "%s.v5_final_review.tsv" % sample)
        write_tsv(outpath, rows)
        write_tsv(reviewpath, [x for x in rows if x["Final_confidence"] == "REVIEW"], list(rows[0].keys()))
        conf = defaultdict(int)
        broad = defaultdict(int)
        for row in rows:
            conf[row["Final_confidence"]] += 1
            broad[row["Final_broad_class"]] += 1
        summary = {"Sample": sample, "Total": len(rows)}
        for key in ("VERY_HIGH","HIGH","MEDIUM_HIGH","MEDIUM","LOW","REVIEW","PRIOR_ONLY","NONTE"):
            summary[key] = conf[key]
        for key in ("DNA","LTR","LINE","SINE","RC","nonTE","Unknown"):
            summary["Final_"+key] = broad[key]
        summaries.append(summary)
        for (confidence, broad_class), count in sorted(
            __import__("collections").Counter((x["Final_confidence"],x["Final_broad_class"]) for x in rows).items()
        ):
            class_rows.append({"Sample": sample, "Confidence": confidence, "Final_broad_class": broad_class, "Count": count})
    write_tsv(os.path.join(outdir, "manifest_samples.v5_final_classification.summary.tsv"), summaries)
    write_tsv(os.path.join(outdir, "manifest_samples.v5_final_classification.by_confidence.tsv"), class_rows,
              ["Sample","Confidence","Final_broad_class","Count"])
    print("Finished strict arbitration for %d sample(s)." % len(entries))

if __name__ == "__main__":
    main()
