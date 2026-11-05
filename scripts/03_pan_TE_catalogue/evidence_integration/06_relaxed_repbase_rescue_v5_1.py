#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v5.1 relaxed RepBase rescue for operational-Unknown TE families.

Generic release derived from multi_sample_relaxed_repbase_v5_1.py.
BLAST is run separately; this script prepares the query FASTA or summarizes
the existing 14-column BLAST tabular output.
"""
from __future__ import print_function
import argparse
import csv
import os
import re
from collections import defaultdict

MIN_IDENTITY = 65.0
MIN_QCOV = 30.0
MIN_QUERY_BP = 100
MAX_EVALUE = 1e-3
MIN_WINNER_SUBJECTS = 2
MIN_WINNER_BITSCORE_FRACTION = 0.70
TE_BROAD = {"DNA","LTR","LINE","SINE","RC"}

def parse_args():
    p=argparse.ArgumentParser(description="Prepare or summarize v5.1 relaxed RepBase rescue.")
    p.add_argument("--manifest",required=True)
    p.add_argument("--sample-id",required=True)
    p.add_argument("--mode",choices=("prepare","summarize"),required=True)
    p.add_argument("--output-dir",required=True)
    return p.parse_args()

def resolve_path(value,base):
    value=(value or "").strip()
    if not value:
        raise ValueError("A required manifest path is blank.")
    return os.path.abspath(value if os.path.isabs(value) else os.path.join(base,value))

def read_inputs(path,sample,mode):
    base=os.path.dirname(os.path.abspath(path))
    required={"sample_id","repeatmodeler_fasta","operational_classification","repbase_fasta","repbase_blast_tsv"}
    rows={}
    with open(path,newline="") as fh:
        reader=csv.DictReader(fh,delimiter="\t")
        missing=required-set(reader.fieldnames or [])
        if missing:
            raise ValueError("Manifest missing columns: "+", ".join(sorted(missing)))
        for line_no,row in enumerate(reader,2):
            sid=(row.get("sample_id") or "").strip()
            if not sid or sid in rows:
                raise ValueError("Blank or duplicate sample_id at line %d"%line_no)
            rows[sid]=row
    if sample not in rows:
        raise ValueError("Sample %s is absent from the manifest."%sample)
    row=rows[sample]
    out={}
    for key in required-{"sample_id"}:
        p=resolve_path(row.get(key),base)
        if key=="repbase_blast_tsv":
            if mode=="summarize" and not os.path.isfile(p):
                raise FileNotFoundError("RepBase BLAST table not found: "+p)
        elif not os.path.isfile(p):
            raise FileNotFoundError("%s not found: %s"%(key,p))
        out[key]=p
    return out

def fasta_id(header):
    token=header.strip().lstrip(">").split()[0]
    return token.split("#",1)[0]

def read_fasta(path):
    data={}
    header=None
    seq=[]
    with open(path) as fh:
        for line in fh:
            line=line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    data[fasta_id(header)]="".join(seq)
                header=line
                seq=[]
            else:
                seq.append(line.strip())
    if header is not None:
        data[fasta_id(header)]="".join(seq)
    return data

def wrap(seq,width=60):
    return "\n".join(seq[i:i+width] for i in range(0,len(seq),width))

def broad_from_rbclass(value):
    s=(value or "").strip().lower()
    if not s:
        return "Unknown"
    if (s in ("sat","msat","satellite","trna","rrna","rp5s","multicopy_gene")
        or "satellite" in s or "simple_repeat" in s or "low_complexity" in s):
        return "nonTE"
    if s in ("transposable_element","nonautonomous","unknown") or s.startswith("loki"):
        return "Unknown"
    special={"dada":"DNA","p":"DNA","sola1":"DNA","sola2":"DNA","isl2eu":"DNA",
             "troyka":"LTR","neptune":"LINE","nimb":"LINE","crack":"LINE",
             "daphne":"LINE","hero":"LINE","proto2":"LINE","r2":"LINE"}
    if s in special:
        return special[s]
    if (s.startswith("dna") or any(x in s for x in ("tc1","mariner","hat","harbinger",
        "piggybac","mudr","mule","transib","kolobok","merlin","cacta","ensp","en-spm"))):
        return "DNA"
    if (s.startswith("ltr") or any(x in s for x in ("gypsy","copia","bel","pao","erv","dirs"))):
        return "LTR"
    if (s.startswith("line") or any(x in s for x in ("non-ltr","cr1","rex","jockey","rte","penelope"))):
        return "LINE"
    if s.startswith("sine") or "sine" in s:
        return "SINE"
    if s.startswith("rc") or "helitron" in s:
        return "RC"
    return "Unknown"

def read_operational(path):
    rows={}
    with open(path,newline="") as fh:
        reader=csv.DictReader(fh,delimiter="\t")
        if "RepeatModeler_ID" not in (reader.fieldnames or []) or "Operational_broad" not in (reader.fieldnames or []):
            raise ValueError("Operational classification must contain RepeatModeler_ID and Operational_broad.")
        for row in reader:
            fid=(row.get("RepeatModeler_ID") or "").strip()
            if fid:
                rows[fid]=row
    return rows

def parse_repbase_headers(path):
    mapping={}
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header=line[1:].strip()
            sid=header.split()[0]
            m=re.search(r"\|RBclass=([^|]+)",header)
            rbclass=m.group(1).strip() if m else "Unknown"
            mapping[sid]={"RBclass":rbclass,"Broad":broad_from_rbclass(rbclass)}
    return mapping

def union_bp(intervals):
    if not intervals:
        return 0
    intervals=sorted(intervals)
    start,end=intervals[0]
    total=0
    for s,e in intervals[1:]:
        if s<=end+1:
            end=max(end,e)
        else:
            total+=end-start+1
            start,end=s,e
    return total+end-start+1

def write_tsv(path,rows,fields):
    with open(path,"w",newline="") as fh:
        writer=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

def prepare_unknown(sample,op,seqs,path):
    unknown_ids=[fid for fid,row in op.items() if row.get("Operational_broad","")=="Unknown"]
    missing=[fid for fid in unknown_ids if fid not in seqs]
    if missing:
        raise RuntimeError("%s: %d operational-Unknown IDs missing from original RepeatModeler FASTA"%(sample,len(missing)))
    with open(path,"w") as fo:
        for fid in sorted(unknown_ids):
            fo.write(">%s\n%s\n"%(fid,wrap(seqs[fid])))
    print("%s operational Unknown families = %d"%(sample,len(unknown_ids)))
    print("Written: "+path)

def summarize(sample,op,seqs,rbmap,blast_path,result_dir,lib_dir):
    unknown_ids={fid for fid,row in op.items() if row.get("Operational_broad","")=="Unknown"}
    missing=[fid for fid in unknown_ids if fid not in seqs]
    if missing:
        raise RuntimeError("%s: %d operational-Unknown IDs missing from original RepeatModeler FASTA"%(sample,len(missing)))
    pairs=defaultdict(list)
    with open(blast_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f=line.rstrip("\n").split("\t")
            if len(f)<14:
                continue
            qid,sid=f[0],f[1]
            if qid not in unknown_ids:
                continue
            pairs[(qid,sid)].append({
                "pident":float(f[2]),"alnlen":int(f[3]),
                "qstart":min(int(f[6]),int(f[7])),"qend":max(int(f[6]),int(f[7])),
                "evalue":float(f[10]),"bitscore":float(f[11]),
                "qlen":int(f[12]),"slen":int(f[13])
            })
    qualifying=defaultdict(list)
    pair_rows=[]
    for (qid,sid),hsps in pairs.items():
        qlen=hsps[0]["qlen"]
        qbp=union_bp([(x["qstart"],x["qend"]) for x in hsps])
        qcov=100.0*qbp/qlen if qlen>0 else 0.0
        total_aln=sum(x["alnlen"] for x in hsps)
        weighted_identity=(sum(x["pident"]*x["alnlen"] for x in hsps)/total_aln) if total_aln>0 else 0.0
        best_evalue=min(x["evalue"] for x in hsps)
        total_bitscore=sum(x["bitscore"] for x in hsps)
        meta=rbmap.get(sid,{"RBclass":"Unknown","Broad":"Unknown"})
        passed=(weighted_identity>=MIN_IDENTITY and qcov>=MIN_QCOV and qbp>=MIN_QUERY_BP and best_evalue<=MAX_EVALUE)
        pair_rows.append({
            "RepeatModeler_ID":qid,"RepBase_subject":sid,"RBclass":meta["RBclass"],"Broad":meta["Broad"],
            "Weighted_identity":"%.3f"%weighted_identity,"Query_covered_bp":qbp,"Query_length":qlen,
            "Query_coverage_pct":"%.3f"%qcov,"Best_evalue":"%.3g"%best_evalue,
            "Total_bitscore":"%.3f"%total_bitscore,"HSP_count":len(hsps),
            "Pass_relaxed":"YES" if passed else "NO"
        })
        if passed:
            qualifying[qid].append({"subject":sid,"rbclass":meta["RBclass"],"broad":meta["Broad"],
                                    "bitscore":total_bitscore,"identity":weighted_identity,
                                    "qcov":qcov,"qbp":qbp,"evalue":best_evalue})
    pair_fields=["RepeatModeler_ID","RepBase_subject","RBclass","Broad","Weighted_identity",
                 "Query_covered_bp","Query_length","Query_coverage_pct","Best_evalue",
                 "Total_bitscore","HSP_count","Pass_relaxed"]
    write_tsv(os.path.join(result_dir,"%s.relaxed_repbase.qualifying_pairs.tsv"%sample),
              sorted(pair_rows,key=lambda z:(z["RepeatModeler_ID"],-float(z["Total_bitscore"]))),pair_fields)

    family_rows=[]
    rescue_map={}
    rescue_counts=defaultdict(int)
    fail_counts=defaultdict(int)
    for qid in sorted(unknown_ids):
        hits=qualifying.get(qid,[])
        n_subjects=len(set(x["subject"] for x in hits))
        if not hits:
            decision,reason,winner,confidence="KEEP_UNKNOWN","NO_QUALIFYING_RELAXED_HIT","Unknown","UNRESOLVED"
            winner_n=0;winner_score=0.0;total_score=0.0;fraction=0.0
            fail_counts[reason]+=1
        else:
            informative=[x for x in hits if x["broad"] in TE_BROAD]
            if not informative:
                decision,reason,winner,confidence="KEEP_UNKNOWN","ONLY_UNINFORMATIVE_REPBASE_CLASSES","Unknown","UNRESOLVED"
                winner_n=0;winner_score=0.0;total_score=0.0;fraction=0.0
                fail_counts[reason]+=1
            else:
                score_by_class=defaultdict(float)
                subjects_by_class=defaultdict(set)
                for x in informative:
                    score_by_class[x["broad"]]+=x["bitscore"]
                    subjects_by_class[x["broad"]].add(x["subject"])
                winner=max(score_by_class,key=score_by_class.get)
                winner_score=score_by_class[winner]
                total_score=sum(score_by_class.values())
                fraction=winner_score/total_score if total_score>0 else 0.0
                winner_n=len(subjects_by_class[winner])
                if winner not in TE_BROAD:
                    decision,reason,confidence="KEEP_UNKNOWN","WINNER_IS_NONTE","UNRESOLVED"
                    fail_counts[reason]+=1
                elif winner_n<MIN_WINNER_SUBJECTS:
                    decision,reason,confidence="KEEP_UNKNOWN","WINNER_HAS_LT2_SUBJECTS","UNRESOLVED"
                    fail_counts[reason]+=1
                elif fraction<MIN_WINNER_BITSCORE_FRACTION:
                    decision,reason,confidence="KEEP_UNKNOWN","BROAD_CLASS_CONSENSUS_LT70","UNRESOLVED"
                    fail_counts[reason]+=1
                else:
                    decision,reason,confidence="RELAXED_REPBASE_RESCUE","PASS","LOW_RELAXED_REPBASE"
                    rescue_map[qid]=winner
                    rescue_counts[winner]+=1
        family_rows.append({
            "RepeatModeler_ID":qid,"V5_operational_broad":"Unknown",
            "Qualifying_subjects":n_subjects,"Winner_broad":winner,
            "Winner_subjects":winner_n,"Winner_bitscore":"%.3f"%winner_score,
            "Informative_total_bitscore":"%.3f"%total_score,
            "Winner_bitscore_fraction":"%.4f"%fraction,
            "V5_1_decision":decision,"V5_1_confidence":confidence,"Reason":reason
        })
    family_fields=["RepeatModeler_ID","V5_operational_broad","Qualifying_subjects","Winner_broad",
                   "Winner_subjects","Winner_bitscore","Informative_total_bitscore",
                   "Winner_bitscore_fraction","V5_1_decision","V5_1_confidence","Reason"]
    write_tsv(os.path.join(result_dir,"%s.relaxed_repbase.family_summary.tsv"%sample),family_rows,family_fields)
    write_tsv(os.path.join(result_dir,"%s.relaxed_repbase.rescued.tsv"%sample),
              [x for x in family_rows if x["V5_1_decision"]=="RELAXED_REPBASE_RESCUE"],family_fields)

    v51_rows=[]
    for fid in sorted(op):
        row=dict(op[fid])
        old_class=row.get("Operational_class","") or "Unknown"
        old_broad=row.get("Operational_broad","") or "Unknown"
        if old_broad=="Unknown" and fid in rescue_map:
            new_broad=rescue_map[fid]
            new_class={"DNA":"DNA/Unknown","LTR":"LTR/Unknown","LINE":"LINE/Unknown",
                       "SINE":"SINE/Unknown","RC":"RC/Helitron"}[new_broad]
            action="RELAXED_REPBASE_RESCUE"
            confidence="LOW_RELAXED_REPBASE"
        else:
            new_broad=old_broad;new_class=old_class
            action="UNCHANGED_FROM_V5"
            confidence=row.get("Operational_confidence","")
        row.update({"V5_1_class":new_class,"V5_1_broad":new_broad,
                    "V5_1_action":action,"V5_1_confidence":confidence})
        v51_rows.append(row)
    if not v51_rows:
        raise RuntimeError("%s operational table is empty."%sample)
    if set(seqs)!=set(x["RepeatModeler_ID"] for x in v51_rows):
        raise RuntimeError("%s FASTA IDs and operational-classification IDs differ."%sample)
    class_fields=list(v51_rows[0].keys())
    write_tsv(os.path.join(result_dir,"%s.v5_1_candidate_classification.tsv"%sample),v51_rows,class_fields)
    final_by_id={x["RepeatModeler_ID"]:x for x in v51_rows}
    repeat_lib=os.path.join(lib_dir,"%s.curated_repeat_v5_1_relaxed.fa"%sample)
    te_lib=os.path.join(lib_dir,"%s.curated_TE_v5_1_relaxed.fa"%sample)
    counts=defaultdict(int)
    with open(repeat_lib,"w") as fr,open(te_lib,"w") as ft:
        for fid in sorted(seqs):
            row=final_by_id[fid];cls=row["V5_1_class"];broad=row["V5_1_broad"]
            record=">%s#%s\n%s\n"%(fid,cls,wrap(seqs[fid]))
            fr.write(record)
            if broad!="nonTE":
                ft.write(record)
            counts[broad]+=1
    summary=[("Operational_Unknown",len(unknown_ids)),
             ("Relaxed_RepBase_rescued",len(rescue_map)),
             ("Remaining_Unknown",counts["Unknown"])]
    for broad in ("DNA","LTR","LINE","SINE","RC","nonTE"):
        summary.append(("Final_"+broad,counts[broad]))
    write_tsv(os.path.join(result_dir,"%s.v5_1.summary.tsv"%sample),
              [{"Metric":k,"Count":v} for k,v in summary],["Metric","Count"])
    print("%s v5.1: operational_Unknown=%d rescued=%d remaining_Unknown=%d"%
          (sample,len(unknown_ids),len(rescue_map),counts["Unknown"]))
    print("Repeat library: "+repeat_lib)
    print("TE-only library: "+te_lib)

def main():
    args=parse_args()
    inputs=read_inputs(args.manifest,args.sample_id,args.mode)
    sample=args.sample_id
    root=os.path.join(os.path.abspath(args.output_dir),sample)
    result_dir=os.path.join(root,"result")
    lib_dir=os.path.join(root,"library")
    unknown_fasta=os.path.join(root,"%s.v5_operational.Unknown.fa"%sample)
    blast_file=inputs["repbase_blast_tsv"]
    os.makedirs(result_dir,exist_ok=True)
    os.makedirs(lib_dir,exist_ok=True)
    os.makedirs(os.path.dirname(blast_file),exist_ok=True)
    op=read_operational(inputs["operational_classification"])
    seqs=read_fasta(inputs["repeatmodeler_fasta"])
    if args.mode=="prepare":
        prepare_unknown(sample,op,seqs,unknown_fasta)
    else:
        rbmap=parse_repbase_headers(inputs["repbase_fasta"])
        summarize(sample,op,seqs,rbmap,blast_file,result_dir,lib_dir)

if __name__=="__main__":
    main()
