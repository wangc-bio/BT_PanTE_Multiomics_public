#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Operational TE classification arbitration.

Generic release derived from operational_arbitration_v5.py. Prior/composite
classifications are retained unless the original evidence rules permit a
rescue or a multi-dimension override. Strict classification is appended as
an independent-evidence confirmation layer.
"""
from __future__ import print_function
import argparse
import csv
import os
from collections import defaultdict

TE_BROAD = {"DNA","LTR","LINE","SINE","RC"}
KNOWN_BROAD = TE_BROAD | {"nonTE"}

def parse_args():
    p = argparse.ArgumentParser(description="Operational arbitration of TE evidence-master tables.")
    p.add_argument("--manifest", required=True, help="TSV with sample_id, evidence_master and strict_classification.")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--review-dir", required=True)
    return p.parse_args()

def resolve_path(value, base):
    value = (value or "").strip()
    if not value:
        raise ValueError("A required manifest path is blank.")
    return os.path.abspath(value if os.path.isabs(value) else os.path.join(base,value))

def read_manifest(path):
    base = os.path.dirname(os.path.abspath(path))
    out, seen = [], set()
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"sample_id","evidence_master","strict_classification"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Manifest missing columns: " + ", ".join(sorted(missing)))
        for line_no,row in enumerate(reader,2):
            sid=(row.get("sample_id") or "").strip()
            if not sid or sid in seen:
                raise ValueError("Blank or duplicate sample_id at line %d" % line_no)
            seen.add(sid)
            item={"sample_id":sid}
            for key in ("evidence_master","strict_classification"):
                p=resolve_path(row.get(key),base)
                if not os.path.isfile(p):
                    raise FileNotFoundError("%s not found: %s" % (key,p))
                item[key]=p
            out.append(item)
    if not out:
        raise ValueError("Manifest contains no sample rows.")
    return out

def broad_from_fine(value):
    s=(value or "").strip().lower()
    if not s:
        return "Unknown"
    if any(x in s for x in ("simple_repeat","simple repeat","satellite","low_complexity","low complexity")) or s in ("trna","rrna","nonte"):
        return "nonTE"
    if "sine" in s:
        return "SINE"
    if "helitron" in s or s.startswith("rc/"):
        return "RC"
    if any(x in s for x in ("ltr","gypsy","copia","bel","pao","erv","dirs")):
        return "LTR"
    if any(x in s for x in ("line","non-ltr","nonltr","cr1","rex","rte","jockey","penelope")):
        return "LINE"
    if any(x in s for x in ("dna","tir","mite","maverick","polinton")):
        return "DNA"
    return "Unknown"

def tier(cls,strength):
    if not cls or cls in ("Unknown","CONFLICT"):
        return "NONE"
    strength=strength or ""
    if strength.startswith("STRONG"):
        return "STRONG"
    if strength.startswith("SUPPORTIVE"):
        return "SUPPORTIVE"
    return "NONE"

def read_table(path,key="RepeatModeler_ID"):
    out={}
    with open(path,newline="") as fh:
        reader=csv.DictReader(fh,delimiter="\t")
        for row in reader:
            fid=(row.get(key) or "").strip()
            if fid:
                out[fid]=row
    return out

def get_external_state(row):
    dims=[
        ("H",row.get("Homology_class",""),row.get("Homology_strength","")),
        ("D",row.get("Domain_class",""),row.get("Domain_strength","")),
        ("S",row.get("Structure_class",""),row.get("Structure_strength",""))
    ]
    strong,supportive=[],[]
    for name,cls,strength in dims:
        t=tier(cls,strength)
        if t=="STRONG":
            strong.append((name,cls))
        elif t=="SUPPORTIVE":
            supportive.append((name,cls))
    strong_classes=sorted(set(cls for _,cls in strong))
    support_classes=sorted(set(cls for _,cls in supportive))
    internal_strong_conflict=False
    for prefix in ("Homology","Domain","Structure"):
        cls=row.get(prefix+"_class","") or ""
        strength=row.get(prefix+"_strength","") or ""
        if cls=="CONFLICT" and strength.startswith("STRONG"):
            internal_strong_conflict=True
    return {"strong":strong,"supportive":supportive,
            "strong_classes":strong_classes,"support_classes":support_classes,
            "internal_strong_conflict":internal_strong_conflict}

def same_count(calls,candidate):
    return sum(1 for _,cls in calls if cls==candidate)

def other_count(calls,candidate):
    return sum(1 for _,cls in calls if cls!=candidate)

def choose_fine_class(row,final_broad,original_class,original_broad):
    if final_broad=="Unknown":
        return "Unknown","NONE"
    if final_broad=="nonTE":
        if original_broad=="nonTE" and original_class:
            return original_class,"ORIGINAL_NONTE"
        return "nonTE","BROAD_ONLY"
    if original_broad==final_broad and original_class and "unknown" not in original_class.lower():
        return original_class,"ORIGINAL_COMPATIBLE_FINE"
    ltr=row.get("LTR_target_class","") or ""
    if final_broad=="LTR" and broad_from_fine(ltr)=="LTR" and "unknown" not in ltr.lower():
        return ltr,"LTRRETRIEVER_FINE"
    mch=row.get("MCHelper_class","") or ""
    if broad_from_fine(mch)==final_broad and "unknown" not in mch.lower():
        return mch,"MCHELPER_COMPATIBLE_FINE"
    tet=row.get("TEtrimmer_class","") or ""
    if broad_from_fine(tet)==final_broad and "unknown" not in tet.lower():
        return tet,"TETRIMMER_COMPATIBLE_FINE"
    defaults={"DNA":"DNA/Unknown","LTR":"LTR/Unknown","LINE":"LINE/Unknown",
              "SINE":"SINE/Unknown","RC":"RC/Helitron","nonTE":"nonTE"}
    return defaults.get(final_broad,"Unknown"),"BROAD_ONLY"

def classify(row,strict_row):
    original_class=row.get("Original_class","") or "Unknown"
    original_broad=row.get("Original_broad","") or "Unknown"
    if original_broad not in KNOWN_BROAD:
        original_broad="Unknown"
    ext=get_external_state(row)
    strong=ext["strong"]
    supportive=ext["supportive"]
    strong_classes=ext["strong_classes"]
    support_classes=ext["support_classes"]
    op_broad="Unknown"
    op_conf="UNRESOLVED"
    op_action="KEEP_UNKNOWN"
    op_reason=""
    weak_conflict="NO"
    manual_review="NO"
    review_reason=""

    if original_broad=="Unknown":
        if ext["internal_strong_conflict"] or len(strong_classes)>1:
            op_action="MANUAL_REVIEW"
            op_conf="REVIEW"
            op_reason="Original Unknown; strong external evidence is internally or cross-dimensionally conflicting."
            manual_review="YES"
            review_reason="STRONG_EXTERNAL_CONFLICT"
        elif len(strong_classes)==1:
            candidate=strong_classes[0]
            nstrong=same_count(strong,candidate)
            nsupport_same=same_count(supportive,candidate)
            nsupport_other=other_count(supportive,candidate)
            if candidate=="nonTE":
                op_action="MANUAL_REVIEW"
                op_conf="REVIEW"
                op_reason="Original Unknown; strong external evidence suggests nonTE. Manual review required."
                manual_review="YES"
                review_reason="EXTERNAL_NONTE_CHALLENGE"
            else:
                op_broad=candidate
                if nstrong>=2 and nsupport_same>=1:
                    op_conf="VERY_HIGH_RESCUE"
                elif nstrong>=2:
                    op_conf="HIGH_RESCUE"
                elif nsupport_same>=1:
                    op_conf="MEDIUM_HIGH_RESCUE"
                else:
                    op_conf="MEDIUM_RESCUE"
                if nsupport_other>0:
                    weak_conflict="YES"
                    op_conf += "_WITH_WEAK_CONFLICT"
                op_action="RESCUE_UNKNOWN"
                op_reason="Original Unknown rescued by at least one STRONG external evidence dimension."
        elif len(support_classes)==1:
            candidate=support_classes[0]
            nsupport=same_count(supportive,candidate)
            if candidate=="nonTE":
                op_broad="Unknown"
                op_conf="UNRESOLVED"
                op_action="KEEP_UNKNOWN"
                op_reason="Original Unknown has only nonTE supportive evidence; not converted to nonTE automatically."
            else:
                op_broad=candidate
                op_action="RESCUE_UNKNOWN"
                op_conf="MEDIUM_RESCUE" if nsupport>=2 else "LOW_RESCUE"
                op_reason=("Original Unknown rescued by concordant SUPPORTIVE evidence from at least two dimensions."
                           if nsupport>=2 else "Original Unknown rescued by one SUPPORTIVE external evidence dimension.")
        else:
            op_broad="Unknown"
            op_conf="UNRESOLVED"
            op_action="KEEP_UNKNOWN"
            if len(support_classes)>1:
                op_reason="Original Unknown has conflicting SUPPORTIVE external classes; no automatic rescue."
                weak_conflict="YES"
            else:
                op_reason="Original Unknown has no usable external evidence."
    else:
        if original_broad=="nonTE":
            te_strong=[x for x in strong_classes if x in TE_BROAD]
            if te_strong:
                op_broad="Unknown"
                op_conf="REVIEW"
                op_action="MANUAL_REVIEW"
                op_reason="Original explicit nonTE is challenged by STRONG TE evidence; automatic override is disabled."
                manual_review="YES"
                review_reason="NONTE_VS_STRONG_TE"
            else:
                op_broad="nonTE"
                op_conf="PRIOR_SUPPORTED" if supportive else "PRIOR_ONLY"
                op_action="PRIOR_RETAINED"
                op_reason="Explicit original nonTE retained; no STRONG TE challenge."
        elif ext["internal_strong_conflict"] or len(strong_classes)>1:
            op_broad="Unknown"
            op_conf="REVIEW"
            op_action="MANUAL_REVIEW"
            op_reason="Known prior encounters strong external class conflict."
            manual_review="YES"
            review_reason="STRONG_EXTERNAL_CONFLICT"
        elif len(strong_classes)==1:
            candidate=strong_classes[0]
            nstrong=same_count(strong,candidate)
            nsupport_same=same_count(supportive,candidate)
            nsupport_other=other_count(supportive,candidate)
            if candidate==original_broad:
                op_broad=original_broad
                op_action="PRIOR_RETAINED"
                if nstrong>=2 and nsupport_same>=1:
                    op_conf="VERY_HIGH_CONFIRMED"
                elif nstrong>=2:
                    op_conf="HIGH_CONFIRMED"
                elif nsupport_same>=1:
                    op_conf="MEDIUM_HIGH_CONFIRMED"
                else:
                    op_conf="MEDIUM_CONFIRMED"
                if nsupport_other>0:
                    weak_conflict="YES"
                    op_conf += "_WITH_WEAK_CONFLICT"
                op_reason="STRONG external evidence is compatible with the original known class."
            elif nstrong==1 and nsupport_same==0:
                op_broad="Unknown"
                op_conf="REVIEW"
                op_action="MANUAL_REVIEW"
                op_reason="A single STRONG external dimension challenges the original known class without a second dimension supporting the new class."
                manual_review="YES"
                review_reason="SINGLE_STRONG_CHALLENGES_PRIOR"
            elif nstrong>=2 or (nstrong==1 and nsupport_same>=1):
                if nsupport_other>0:
                    op_broad="Unknown"
                    op_conf="REVIEW"
                    op_action="MANUAL_REVIEW"
                    op_reason="New class has multi-dimension support but another external dimension supports a different class."
                    manual_review="YES"
                    review_reason="OVERRIDE_WITH_EXTERNAL_CONFLICT"
                else:
                    op_broad=candidate
                    op_action="EXTERNAL_OVERRIDE_PRIOR"
                    if nstrong>=2 and nsupport_same>=1:
                        op_conf="VERY_HIGH_OVERRIDE"
                    elif nstrong>=2:
                        op_conf="HIGH_OVERRIDE"
                    else:
                        op_conf="MEDIUM_HIGH_OVERRIDE"
                    op_reason="Original known class overridden because at least two independent evidence dimensions support the same alternative class."
            else:
                op_broad="Unknown"
                op_conf="REVIEW"
                op_action="MANUAL_REVIEW"
                op_reason="Strong external challenge could not be resolved automatically."
                manual_review="YES"
                review_reason="UNRESOLVED_STRONG_CHALLENGE"
        else:
            op_broad=original_broad
            op_action="PRIOR_RETAINED"
            if not supportive:
                op_conf="PRIOR_ONLY"
                op_reason="Original known class retained because there is no qualifying external evidence and no strong challenge."
            else:
                support_same=same_count(supportive,original_broad)
                support_other=other_count(supportive,original_broad)
                if support_other==0:
                    op_conf="PRIOR_SUPPORTED"
                    op_reason="Original known class retained and supported by compatible SUPPORTIVE external evidence."
                elif support_same>0:
                    op_conf="PRIOR_WITH_WEAK_CONFLICT"
                    weak_conflict="YES"
                    op_reason="Original known class retained; SUPPORTIVE external evidence is mixed but no STRONG contradiction exists."
                else:
                    op_conf=("PRIOR_CHALLENGED_BY_MULTI_SUPPORT"
                             if len(support_classes)==1 and len(supportive)>=2
                             else "PRIOR_WITH_WEAK_CONFLICT")
                    weak_conflict="YES"
                    op_reason="Only SUPPORTIVE evidence challenges the original known class; weak evidence is insufficient to override the prior automatically."

    op_class,fine_source=choose_fine_class(row,op_broad,original_class,original_broad)
    out=dict(row)
    out.update({
        "Evidence_confirmed_class":strict_row.get("Final_class",""),
        "Evidence_confirmed_broad":strict_row.get("Final_broad_class",""),
        "Evidence_confirmed_confidence":strict_row.get("Final_confidence",""),
        "Operational_class":op_class,
        "Operational_broad":op_broad,
        "Operational_confidence":op_conf,
        "Operational_action":op_action,
        "Operational_reason":op_reason,
        "Fine_class_source":fine_source,
        "Weak_conflict_flag":weak_conflict,
        "Manual_review":manual_review,
        "Manual_review_reason":review_reason
    })
    return out

def write_tsv(path,rows,fields=None):
    if fields is None:
        fields=list(rows[0].keys()) if rows else []
    with open(path,"w",newline="") as fh:
        writer=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

def main():
    args=parse_args()
    entries=read_manifest(args.manifest)
    outdir,reviewdir=os.path.abspath(args.output_dir),os.path.abspath(args.review_dir)
    os.makedirs(outdir,exist_ok=True)
    os.makedirs(reviewdir,exist_ok=True)
    summary_rows,action_rows,class_rows=[],[],[]
    for entry in entries:
        sample=entry["sample_id"]
        with open(entry["evidence_master"],newline="") as fh:
            reader=csv.DictReader(fh,delimiter="\t")
            master_fields=set(reader.fieldnames or [])
            required={"RepeatModeler_ID","Original_class","Original_broad",
                      "Homology_class","Homology_strength","Domain_class","Domain_strength",
                      "Structure_class","Structure_strength","LTR_target_class","MCHelper_class","TEtrimmer_class"}
            missing=required-master_fields
            if missing:
                raise ValueError("%s evidence master missing columns: %s"%(sample,", ".join(sorted(missing))))
            master=list(reader)
        strict=read_table(entry["strict_classification"])
        if not master:
            raise RuntimeError("%s evidence master contains zero rows."%sample)
        rows=[classify(row,strict.get(row["RepeatModeler_ID"],{})) for row in master]
        if not rows:
            raise RuntimeError("%s produced zero rows."%sample)
        fields=list(rows[0].keys())
        write_tsv(os.path.join(outdir,"%s.v5_operational_classification.tsv"%sample),rows,fields)
        write_tsv(os.path.join(reviewdir,"%s.v5_operational_manual_review.tsv"%sample),
                  [x for x in rows if x["Manual_review"]=="YES"],fields)
        actions=defaultdict(int); broad=defaultdict(int)
        for row in rows:
            actions[row["Operational_action"]]+=1
            broad[row["Operational_broad"]]+=1
        rescue_n=actions["RESCUE_UNKNOWN"]
        prior_n=actions["PRIOR_RETAINED"]
        override_n=actions["EXTERNAL_OVERRIDE_PRIOR"]
        manual_n=actions["MANUAL_REVIEW"]
        unknown_n=actions["KEEP_UNKNOWN"]
        if rescue_n+prior_n+override_n+manual_n+unknown_n != len(rows):
            raise RuntimeError("%s operational action accounting mismatch."%sample)
        summary_rows.append({
            "Sample":sample,"Total":len(rows),
            "Original_Unknown":sum(1 for x in rows if (x.get("Original_broad") or "Unknown")=="Unknown"),
            "Original_Known":sum(1 for x in rows if (x.get("Original_broad") or "Unknown") in KNOWN_BROAD),
            "Rescue_from_Unknown":rescue_n,"Prior_retained":prior_n,
            "External_override_prior":override_n,"Manual_review":manual_n,
            "Unresolved_Unknown":unknown_n,"Operational_Unknown_total":broad["Unknown"],
            "Operational_DNA":broad["DNA"],"Operational_LTR":broad["LTR"],
            "Operational_LINE":broad["LINE"],"Operational_SINE":broad["SINE"],
            "Operational_RC":broad["RC"],"Operational_nonTE":broad["nonTE"]
        })
        for key,n in sorted(actions.items()):
            action_rows.append({"Sample":sample,"Operational_action":key,"Count":n})
        for key,n in sorted(broad.items()):
            class_rows.append({"Sample":sample,"Operational_broad":key,"Count":n})
    write_tsv(os.path.join(outdir,"manifest_samples.v5_operational.summary.tsv"),summary_rows)
    write_tsv(os.path.join(outdir,"manifest_samples.v5_operational.actions.tsv"),action_rows,
              ["Sample","Operational_action","Count"])
    write_tsv(os.path.join(outdir,"manifest_samples.v5_operational.by_class.tsv"),class_rows,
              ["Sample","Operational_broad","Count"])
    print("Finished operational arbitration for %d sample(s)."%len(entries))

if __name__=="__main__":
    main()
