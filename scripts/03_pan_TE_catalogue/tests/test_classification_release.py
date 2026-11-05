#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Synthetic end-to-end tests for the four Module 03 classification scripts."""
import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EI = ROOT / "scripts" / "03_pan_TE_catalogue" / "evidence_integration"

def write_tsv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

def read_tsv(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def run(script, *args):
    subprocess.run([sys.executable, str(EI / script), *map(str, args)],
                   check=True, capture_output=True, text=True)

def blast_row(qid, sid, pident, alnlen, qlen, bitscore):
    return "\t".join([
        qid, sid, str(pident), str(alnlen), "0", "0",
        "1", str(alnlen), "1", str(alnlen),
        "1e-40", str(bitscore), str(qlen), "100"
    ])

def main():
    with tempfile.TemporaryDirectory(prefix="module03_classification_test_") as tmp:
        work = Path(tmp)
        evidence = work / "input" / "evidence_master.tsv"
        master_fields = [
            "RepeatModeler_ID","Original_class","Original_broad",
            "TEtrimmer_class","MCHelper_class","LTR_target_class",
            "Homology_class","Homology_strength","Homology_conflict",
            "Domain_class","Domain_strength","Domain_conflict",
            "Structure_class","Structure_strength","Structure_conflict"
        ]
        base = {
            "TEtrimmer_class":"Unknown","MCHelper_class":"Unknown","LTR_target_class":"",
            "Homology_class":"","Homology_strength":"NONE","Homology_conflict":"NO",
            "Domain_class":"","Domain_strength":"NONE","Domain_conflict":"NO",
            "Structure_class":"","Structure_strength":"NONE","Structure_conflict":"NO"
        }
        rows=[]
        def add(fid, original_class, original_broad, h=("", "NONE"), d=("", "NONE"), s=("", "NONE")):
            row=dict(base)
            row.update({
                "RepeatModeler_ID":fid,"Original_class":original_class,"Original_broad":original_broad,
                "Homology_class":h[0],"Homology_strength":h[1],
                "Domain_class":d[0],"Domain_strength":d[1],
                "Structure_class":s[0],"Structure_strength":s[1]
            })
            rows.append(row)
        add("A","Unknown","Unknown",("DNA","STRONG"),("DNA","STRONG"))
        add("B","Unknown","Unknown",("LTR","SUPPORTIVE"),("LTR","SUPPORTIVE"))
        add("C","Unknown","Unknown",("DNA","STRONG"),("LTR","STRONG"))
        add("D","LTR/Gypsy","LTR")
        add("E","Unknown","Unknown")
        write_tsv(evidence, rows, master_fields)

        strict_dir=work/"strict"
        strict_review=work/"strict_review"
        op_dir=work/"operational"
        op_review=work/"operational_review"
        strict_file=strict_dir/"TEST.v5_final_classification.tsv"
        op_file=op_dir/"TEST.v5_operational_classification.tsv"

        rm_fasta=work/"input"/"TEST.consensi.fa.classified"
        rm_fasta.write_text("".join(">%s#Unknown\n%s\n"%(x,"A"*100) for x in "ABCDE"))
        repbase_fasta=work/"input"/"RepBase31.07.clean.fa"
        repbase_fasta.write_text(
            ">RB1|RBclass=DNA/CMC-EnSpm\n"+"A"*100+"\n"+
            ">RB2|RBclass=DNA/CMC-EnSpm\n"+"A"*100+"\n"+
            ">RB3|RBclass=LTR/Gypsy\n"+"A"*100+"\n"
        )
        dfam_fasta=work/"input"/"Cyprinidae.TE_only.fa"
        dfam_fasta.write_text(
            ">DF1#LTR/Gypsy|Broad=LTR\n"+"A"*100+"\n"+
            ">DF2#LTR/Gypsy|Broad=LTR\n"+"A"*100+"\n"
        )
        repbase_blast=work/"blast"/"TEST.Unknown_vs_RepBase.tsv"
        dfam_blast=work/"blast"/"TEST.v51Unknown_vs_Dfam.tsv"
        manifest=work/"classification_manifest.tsv"
        manifest_fields=[
            "sample_id","evidence_master","strict_classification","repeatmodeler_fasta",
            "operational_classification","repbase_fasta","repbase_blast_tsv",
            "v5_1_classification","dfam_fasta","dfam_blast_tsv"
        ]
        write_tsv(manifest,[{
            "sample_id":"TEST","evidence_master":str(evidence),"strict_classification":str(strict_file),
            "repeatmodeler_fasta":str(rm_fasta),"operational_classification":str(op_file),
            "repbase_fasta":str(repbase_fasta),"repbase_blast_tsv":str(repbase_blast),
            "v5_1_classification":str(work/"v51"/"TEST"/"result"/"TEST.v5_1_candidate_classification.tsv"),
            "dfam_fasta":str(dfam_fasta),"dfam_blast_tsv":str(dfam_blast)
        }],manifest_fields)

        run("04_final_arbitration.py","--manifest",manifest,"--output-dir",strict_dir,"--review-dir",strict_review)
        strict={x["RepeatModeler_ID"]:x for x in read_tsv(strict_file)}
        assert strict["A"]["Final_broad_class"]=="DNA"
        assert strict["B"]["Final_broad_class"]=="LTR"
        assert strict["C"]["Final_broad_class"]=="Unknown"
        assert strict["C"]["Final_confidence"]=="REVIEW"
        assert strict["D"]["Final_broad_class"]=="Unknown"
        assert strict["D"]["Final_confidence"]=="PRIOR_ONLY"

        run("05_operational_arbitration.py","--manifest",manifest,"--output-dir",op_dir,"--review-dir",op_review)
        operational={x["RepeatModeler_ID"]:x for x in read_tsv(op_file)}
        assert operational["A"]["Operational_broad"]=="DNA"
        assert operational["A"]["Operational_action"]=="RESCUE_UNKNOWN"
        assert operational["B"]["Operational_broad"]=="LTR"
        assert operational["B"]["Operational_action"]=="RESCUE_UNKNOWN"
        assert operational["C"]["Operational_broad"]=="Unknown"
        assert operational["C"]["Operational_action"]=="MANUAL_REVIEW"
        assert operational["D"]["Operational_broad"]=="LTR"
        assert operational["D"]["Operational_action"]=="PRIOR_RETAINED"
        assert operational["E"]["Operational_broad"]=="Unknown"

        repbase_blast.parent.mkdir(parents=True,exist_ok=True)
        repbase_blast.write_text(
            blast_row("C","RB1|RBclass=DNA/CMC-EnSpm",90,100,100,200)+"\n"+
            blast_row("C","RB2|RBclass=DNA/CMC-EnSpm",90,100,100,200)+"\n"
        )
        v51_out=work/"v51"
        run("06_relaxed_repbase_rescue_v5_1.py","--manifest",manifest,"--sample-id","TEST",
            "--mode","prepare","--output-dir",v51_out)
        assert (v51_out/"TEST"/"TEST.v5_operational.Unknown.fa").is_file()
        run("06_relaxed_repbase_rescue_v5_1.py","--manifest",manifest,"--sample-id","TEST",
            "--mode","summarize","--output-dir",v51_out)
        v51_file=v51_out/"TEST"/"result"/"TEST.v5_1_candidate_classification.tsv"
        v51={x["RepeatModeler_ID"]:x for x in read_tsv(v51_file)}
        assert v51["C"]["V5_1_broad"]=="DNA"
        assert v51["C"]["V5_1_action"]=="RELAXED_REPBASE_RESCUE"
        assert v51["E"]["V5_1_broad"]=="Unknown"

        dfam_blast.parent.mkdir(parents=True,exist_ok=True)
        dfam_blast.write_text(
            blast_row("E","DF1#LTR/Gypsy|Broad=LTR",90,80,100,100)+"\n"+
            blast_row("E","DF2#LTR/Gypsy|Broad=LTR",90,80,100,100)+"\n"
        )
        v52_out=work/"v52"
        run("07_relaxed_homology_rescue_v5_2.py","--manifest",manifest,"--sample-id","TEST",
            "--mode","prepare","--output-dir",v52_out)
        assert (v52_out/"TEST"/"TEST.v5_1_remaining_Unknown.fa").is_file()
        run("07_relaxed_homology_rescue_v5_2.py","--manifest",manifest,"--sample-id","TEST",
            "--mode","summarize","--output-dir",v52_out)
        v52_file=v52_out/"TEST"/"result"/"TEST.v5_2_classification.tsv"
        v52={x["RepeatModeler_ID"]:x for x in read_tsv(v52_file)}
        assert v52["C"]["V5_2_broad"]=="DNA"
        assert v52["C"]["V5_2_action"]=="UNCHANGED_FROM_V5_1"
        v52_calls_file=v52_out/"TEST"/"result"/"TEST.v52_homology_calls.tsv"
        v52_calls={x["RepeatModeler_ID"]:x for x in read_tsv(v52_calls_file)}
        assert v52["E"]["V5_2_broad"]=="LTR", "Unexpected v5.2 E classification: final=%r; evidence=%r" % (v52["E"], v52_calls.get("E"))
        assert v52["E"]["V5_2_action"]=="DFAM_ONLY_RESCUE"
        assert (v52_out/"TEST"/"library"/"TEST.curated_repeat_v5_2.fa").is_file()
        assert (v52_out/"TEST"/"library"/"TEST.curated_TE_v5_2.fa").is_file()

    print("Module 03 classification synthetic end-to-end test: PASS")

if __name__=="__main__":
    main()
