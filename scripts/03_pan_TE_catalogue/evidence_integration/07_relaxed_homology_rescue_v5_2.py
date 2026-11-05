#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v5.2 relaxed RepBase/Dfam rescue for v5.1 remaining-Unknown families.

Generic release derived from v52_relaxed_homology_rescue.py. BLAST is run
separately; this script prepares the v5.1 remaining-Unknown query FASTA or
summarizes existing RepBase and Dfam 14-column BLAST tables.
"""
from __future__ import print_function
import argparse
import csv
import os
import re
from collections import defaultdict

MIN_IDENTITY=65.0
MIN_QCOV=20.0
MIN_QUERY_BP=80
MAX_EVALUE=1e-3
MIN_SUBJECTS=2
MIN_WINNER_FRACTION=0.70
TOP_N_SUBJECTS=10
TE_BROAD={"DNA","LTR","LINE","SINE","RC"}

def parse_args():
    p=argparse.ArgumentParser(description="Prepare or summarize v5.2 relaxed RepBase/Dfam rescue.")
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
    required={"sample_id","repeatmodeler_fasta","v5_1_classification",
              "repbase_blast_tsv","repbase_fasta","dfam_blast_tsv","dfam_fasta"}
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
    row=rows[sample];out={}
    for key in required-{"sample_id"}:
        p=resolve_path(row.get(key),base)
        if key in ("repbase_blast_tsv","dfam_blast_tsv"):
            if mode=="summarize" and not os.path.isfile(p):
                raise FileNotFoundError("%s not found: %s"%(key,p))
        elif not os.path.isfile(p):
            raise FileNotFoundError("%s not found: %s"%(key,p))
        out[key]=p
    return out

def fasta_id(header):
    token=header.strip().lstrip(">").split()[0]
    return token.split("#",1)[0]

def read_fasta(path):
    data={};header=None;seq=[]
    with open(path) as fh:
        for line in fh:
            line=line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    data[fasta_id(header)]="".join(seq)
                header=line;seq=[]
            else:
                seq.append(line.strip())
    if header is not None:
        data[fasta_id(header)]="".join(seq)
    return data

def wrap(seq,width=60):
    return "\n".join(seq[i:i+width] for i in range(0,len(seq),width))

def union_bp(intervals):
    if not intervals:
        return 0
    intervals=sorted(intervals)
    start,end=intervals[0];total=0
    for s,e in intervals[1:]:
        if s<=end+1:
            end=max(end,e)
        else:
            total+=end-start+1
            start,end=s,e
    return total+end-start+1

def broad_repbase(value):
    s=(value or "").strip().lower()
    if not s:
        return "Unknown"
    if ("satellite" in s or "simple_repeat" in s or "low_complexity" in s
        or s in ("sat","msat","trna","rrna","rp5s","multicopy_gene")):
        return "nonTE"
    if s in ("transposable_element","nonautonomous","unknown") or s.startswith("loki"):
        return "Unknown"
    special={"dada":"DNA","p":"DNA","sola1":"DNA","sola2":"DNA","isl2eu":"DNA",
             "troyka":"LTR","neptune":"LINE","nimb":"LINE","crack":"LINE",
             "daphne":"LINE","hero":"LINE","proto2":"LINE","r2":"LINE"}
    if s in special:
        return special[s]
    if (s.startswith("dna") or any(x in s for x in ("tc1","mariner","hat","harbinger",
        "piggybac","mudr","mule","transib","kolobok","merlin","cacta"))):
        return "DNA"
    if any(x in s for x in ("ltr","gypsy","copia","bel","pao","erv","dirs")):
        return "LTR"
    if any(x in s for x in ("line","non-ltr","cr1","rex","rte","jockey","penelope")):
        return "LINE"
    if "sine" in s:
        return "SINE"
    if s.startswith("rc") or "helitron" in s:
        return "RC"
    return "Unknown"

def parse_repbase_map(path):
    mapping={}
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header=line[1:].strip();sid=header.split()[0]
            m=re.search(r"\|RBclass=([^|]+)",header)
            fine=m.group(1).strip() if m else "Unknown"
            mapping[sid]={"fine":fine,"broad":broad_repbase(fine)}
    return mapping

def parse_dfam_map(path):
    mapping={}
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header=line[1:].strip();sid=header.split()[0]
            m=re.search(r"\|Broad=([^|\s]+)",sid)
            broad=m.group(1) if m else "Unknown"
            fine=sid.split("#",1)[1].split("|Broad=",1)[0] if "#" in sid else "Unknown"
            mapping[sid]={"fine":fine,"broad":broad}
    return mapping

def read_v51(path):
    rows={}
    with open(path,newline="") as fh:
        reader=csv.DictReader(fh,delimiter="\t")
        if not {"RepeatModeler_ID","V5_1_broad","V5_1_class"}.issubset(set(reader.fieldnames or [])):
            raise ValueError("v5.1 table must contain RepeatModeler_ID, V5_1_broad and V5_1_class.")
        for row in reader:
            fid=(row.get("RepeatModeler_ID") or "").strip()
            if fid:
                rows[fid]=row
    return rows

def read_blast(path,allowed_ids,class_map):
    hsps=defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f=line.rstrip("\n").split("\t")
            if len(f)<14:
                continue
            qid,sid=f[0],f[1]
            if qid not in allowed_ids or sid not in class_map:
                continue
            hsps[(qid,sid)].append({
                "pident":float(f[2]),"alnlen":int(f[3]),
                "qstart":min(int(f[6]),int(f[7])),"qend":max(int(f[6]),int(f[7])),
                "evalue":float(f[10]),"bitscore":float(f[11]),"qlen":int(f[12])
            })
    by_query=defaultdict(list)
    for (qid,sid),arr in hsps.items():
        qlen=arr[0]["qlen"]
        qbp=union_bp([(x["qstart"],x["qend"]) for x in arr])
        qcov=100.0*qbp/qlen if qlen>0 else 0.0
        aln_total=sum(x["alnlen"] for x in arr)
        identity=(sum(x["pident"]*x["alnlen"] for x in arr)/aln_total) if aln_total>0 else 0.0
        evalue=min(x["evalue"] for x in arr)
        bitscore=sum(x["bitscore"] for x in arr)
        meta=class_map[sid]
        passed=identity>=MIN_IDENTITY and qcov>=MIN_QCOV and qbp>=MIN_QUERY_BP and evalue<=MAX_EVALUE
        if passed:
            by_query[qid].append({"subject":sid,"broad":meta["broad"],"fine":meta["fine"],
                                  "identity":identity,"qcov":qcov,"qbp":qbp,
                                  "evalue":evalue,"bitscore":bitscore})
    return by_query

def call_database(hits):
    if not hits:
        return {"pass":False,"class":"Unknown","subjects":0,"fraction":0.0,"reason":"NO_QUALIFYING_HIT"}
    informative=[x for x in hits if x["broad"] in TE_BROAD]
    if not informative:
        return {"pass":False,"class":"Unknown","subjects":0,"fraction":0.0,"reason":"NO_CLASSIFIED_TE_HIT"}
    informative=sorted(informative,key=lambda x:-x["bitscore"])[:TOP_N_SUBJECTS]
    scores=defaultdict(float);subjects=defaultdict(set)
    for x in informative:
        scores[x["broad"]]+=x["bitscore"]
        subjects[x["broad"]].add(x["subject"])
    winner=max(scores,key=scores.get)
    total=sum(scores.values())
    fraction=scores[winner]/total if total>0 else 0.0
    nsub=len(subjects[winner])
    if nsub<MIN_SUBJECTS:
        return {"pass":False,"class":winner,"subjects":nsub,"fraction":fraction,"reason":"WINNER_LT2_SUBJECTS"}
    if fraction<MIN_WINNER_FRACTION:
        return {"pass":False,"class":winner,"subjects":nsub,"fraction":fraction,"reason":"WINNER_FRACTION_LT70"}
    return {"pass":True,"class":winner,"subjects":nsub,"fraction":fraction,"reason":"PASS"}

def write_tsv(path,rows,fields):
    with open(path,"w",newline="") as fh:
        writer=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n")
        writer.writeheader();writer.writerows(rows)

def prepare(sample,v51,seqs,query_path):
    ids=[fid for fid,row in v51.items() if row.get("V5_1_broad","")=="Unknown"]
    missing=[fid for fid in ids if fid not in seqs]
    if missing:
        raise RuntimeError("%s: %d v5.1 remaining-Unknown IDs missing from original FASTA"%(sample,len(missing)))
    with open(query_path,"w") as fo:
        for fid in sorted(ids):
            fo.write(">%s\n%s\n"%(fid,wrap(seqs[fid])))
    print("%s v5.1 remaining Unknown = %d"%(sample,len(ids)))
    print("Written: "+query_path)

def summarize(sample,v51,seqs,repmap,dfammap,repbase_blast,dfam_blast,result_dir,lib_dir):
    remaining={fid for fid,row in v51.items() if row.get("V5_1_broad","")=="Unknown"}
    missing=[fid for fid in remaining if fid not in seqs]
    if missing:
        raise RuntimeError("%s: %d v5.1 remaining-Unknown IDs missing from original FASTA"%(sample,len(missing)))
    rep_hits=read_blast(repbase_blast,remaining,repmap)
    dfam_hits=read_blast(dfam_blast,remaining,dfammap)
    rescue={};source_rows=[];counts=defaultdict(int);rescue_class=defaultdict(int)
    for fid in sorted(remaining):
        r=call_database(rep_hits.get(fid,[]))
        d=call_database(dfam_hits.get(fid,[]))
        if r["pass"] and d["pass"]:
            if r["class"]==d["class"]:
                final=r["class"];action="REPBASE_DFAM_CONCORDANT_RESCUE";confidence="MEDIUM_RELAXED_HOMOLOGY"
                rescue[fid]=final;counts["BOTH_CONCORDANT"]+=1;rescue_class[final]+=1
            else:
                final="Unknown";action="HOMOLOGY_DB_CONFLICT";confidence="UNRESOLVED";counts["DB_CONFLICT"]+=1
        elif r["pass"]:
            final=r["class"];action="REPBASE_ONLY_RESCUE";confidence="LOW_RELAXED_REPBASE"
            rescue[fid]=final;counts["REPBASE_ONLY"]+=1;rescue_class[final]+=1
        elif d["pass"]:
            final=d["class"];action="DFAM_ONLY_RESCUE";confidence="LOW_RELAXED_DFAM"
            rescue[fid]=final;counts["DFAM_ONLY"]+=1;rescue_class[final]+=1
        else:
            final="Unknown";action="KEEP_UNKNOWN";confidence="UNRESOLVED";counts["NEITHER"]+=1
        source_rows.append({
            "RepeatModeler_ID":fid,
            "RepBase_pass":"YES" if r["pass"] else "NO","RepBase_class":r["class"],
            "RepBase_subjects":r["subjects"],"RepBase_fraction":"%.4f"%r["fraction"],
            "RepBase_reason":r["reason"],
            "Dfam_pass":"YES" if d["pass"] else "NO","Dfam_class":d["class"],
            "Dfam_subjects":d["subjects"],"Dfam_fraction":"%.4f"%d["fraction"],
            "Dfam_reason":d["reason"],"V5_2_action":action,
            "V5_2_broad":final,"V5_2_confidence":confidence
        })
    source_fields=["RepeatModeler_ID","RepBase_pass","RepBase_class","RepBase_subjects","RepBase_fraction",
                   "RepBase_reason","Dfam_pass","Dfam_class","Dfam_subjects","Dfam_fraction","Dfam_reason",
                   "V5_2_action","V5_2_broad","V5_2_confidence"]
    write_tsv(os.path.join(result_dir,"%s.v52_homology_calls.tsv"%sample),source_rows,source_fields)
    rescue_map={x["RepeatModeler_ID"]:x for x in source_rows}
    defaults={"DNA":"DNA/Unknown","LTR":"LTR/Unknown","LINE":"LINE/Unknown","SINE":"SINE/Unknown","RC":"RC/Helitron"}
    final_rows=[]
    for fid in sorted(v51):
        row=dict(v51[fid])
        old_broad=row.get("V5_1_broad","") or "Unknown"
        old_class=row.get("V5_1_class","") or "Unknown"
        if old_broad=="Unknown" and fid in rescue:
            new_broad=rescue[fid]
            new_class=defaults[new_broad]
            callrow=rescue_map[fid]
            new_action=callrow["V5_2_action"]
            new_confidence=callrow["V5_2_confidence"]
        else:
            new_broad=old_broad;new_class=old_class
            new_action="UNCHANGED_FROM_V5_1"
            new_confidence=row.get("V5_1_confidence","")
        row.update({"V5_2_class":new_class,"V5_2_broad":new_broad,
                    "V5_2_action":new_action,"V5_2_confidence":new_confidence})
        final_rows.append(row)
    if set(seqs)!=set(x["RepeatModeler_ID"] for x in final_rows):
        raise RuntimeError("%s FASTA IDs and v5.1 classification IDs differ."%sample)
    class_fields=list(final_rows[0].keys())
    write_tsv(os.path.join(result_dir,"%s.v5_2_classification.tsv"%sample),final_rows,class_fields)
    final_map={x["RepeatModeler_ID"]:x for x in final_rows}
    repeat_lib=os.path.join(lib_dir,"%s.curated_repeat_v5_2.fa"%sample)
    te_lib=os.path.join(lib_dir,"%s.curated_TE_v5_2.fa"%sample)
    final_counts=defaultdict(int)
    with open(repeat_lib,"w") as fr,open(te_lib,"w") as ft:
        for fid in sorted(seqs):
            row=final_map[fid];cls=row["V5_2_class"];broad=row["V5_2_broad"]
            record=">%s#%s\n%s\n"%(fid,cls,wrap(seqs[fid]))
            fr.write(record)
            if broad!="nonTE":
                ft.write(record)
            final_counts[broad]+=1
    total_rescue=sum(rescue_class.values())
    summary=[("V5_1_Remaining_Unknown",len(remaining)),
             ("RepBase_only_rescue",counts["REPBASE_ONLY"]),
             ("Dfam_only_rescue",counts["DFAM_ONLY"]),
             ("RepBase_Dfam_concordant",counts["BOTH_CONCORDANT"]),
             ("DB_conflict_keep_Unknown",counts["DB_CONFLICT"]),
             ("Neither_keep_Unknown",counts["NEITHER"]),
             ("V5_2_rescued_total",total_rescue)]
    for broad in ("DNA","LTR","LINE","SINE","RC"):
        summary.append(("V5_2_rescued_"+broad,rescue_class[broad]))
    summary.append(("V5_2_remaining_Unknown",final_counts["Unknown"]))
    summary.extend(("Final_"+broad,final_counts[broad]) for broad in ("DNA","LTR","LINE","SINE","RC","Unknown","nonTE"))
    write_tsv(os.path.join(result_dir,"%s.v5_2.summary.tsv"%sample),
              [{"Metric":k,"Count":v} for k,v in summary],["Metric","Count"])
    print("%s v5.2: remaining_unknown=%d rescued=%d final_unknown=%d"%
          (sample,len(remaining),total_rescue,final_counts["Unknown"]))
    print("Repeat library: "+repeat_lib)
    print("TE-only library: "+te_lib)

def main():
    args=parse_args()
    inputs=read_inputs(args.manifest,args.sample_id,args.mode)
    sample=args.sample_id
    work=os.path.join(os.path.abspath(args.output_dir),sample)
    result_dir=os.path.join(work,"result")
    lib_dir=os.path.join(work,"library")
    query_fasta=os.path.join(work,"%s.v5_1_remaining_Unknown.fa"%sample)
    os.makedirs(result_dir,exist_ok=True);os.makedirs(lib_dir,exist_ok=True)
    os.makedirs(os.path.dirname(inputs["dfam_blast_tsv"]),exist_ok=True)
    v51=read_v51(inputs["v5_1_classification"])
    seqs=read_fasta(inputs["repeatmodeler_fasta"])
    if args.mode=="prepare":
        prepare(sample,v51,seqs,query_fasta)
    else:
        repmap=parse_repbase_map(inputs["repbase_fasta"])
        dfammap=parse_dfam_map(inputs["dfam_fasta"])
        summarize(sample,v51,seqs,repmap,dfammap,inputs["repbase_blast_tsv"],
                  inputs["dfam_blast_tsv"],result_dir,lib_dir)

if __name__=="__main__":
    main()
