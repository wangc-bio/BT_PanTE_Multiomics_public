#!/usr/bin/env bash
#SBATCH --job-name=02_parent_specific_kmer_mapping
#SBATCH --cpus-per-task=36
#SBATCH --mem=160G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
export LC_ALL=C
: "${BSB_DUMP:?}"; : "${TC_DUMP:?}"; : "${BT_FASTA:?}"
: "${BOWTIE:?}"; : "${BOWTIE_BUILD:?}"; : "${BEDTOOLS:?}"; : "${SAMTOOLS:?}"; : "${PYTHON:?}"
: "${OUTDIR:?}"
THREADS=${SLURM_CPUS_PER_TASK:-36}
MIN_COPY=1; MAX_COPY=5; WIN=1000000; MIN_HITS=50; STRONG=0.70
mkdir -p "$OUTDIR" logs; cd "$OUTDIR"
for f in "$BSB_DUMP" "$TC_DUMP" "$BT_FASTA"; do [[ -s "$f" ]] || { echo "Missing $f" >&2; exit 1; }; done
awk -v lo=$MIN_COPY -v hi=$MAX_COPY '$2>=lo && $2<=hi {print $1}' "$BSB_DUMP" | sort -u > BSB.k31.filtered.txt
awk -v lo=$MIN_COPY -v hi=$MAX_COPY '$2>=lo && $2<=hi {print $1}' "$TC_DUMP" | sort -u > TC.k31.filtered.txt
comm -23 BSB.k31.filtered.txt TC.k31.filtered.txt > BSB.specific.k31.txt
comm -13 BSB.k31.filtered.txt TC.k31.filtered.txt > TC.specific.k31.txt
awk '{print ">BSB_kmer_"NR"\n"$1}' BSB.specific.k31.txt > BSB.specific.k31.fa
awk '{print ">TC_kmer_"NR"\n"$1}' TC.specific.k31.txt > TC.specific.k31.fa
"$BOWTIE_BUILD" "$BT_FASTA" BT
"$BOWTIE" -S -f -v 0 -a -m 1 --best --strata -p "$THREADS" -x BT BSB.specific.k31.fa > BSB.k31.unique.sam
"$BOWTIE" -S -f -v 0 -a -m 1 --best --strata -p "$THREADS" -x BT TC.specific.k31.fa > TC.k31.unique.sam
for parent in BSB TC; do
  awk 'BEGIN{OFS="\t"} $1!~/^@/ && $3!="*" {s=$4-1; print $3,s,s+length($10)}' "${parent}.k31.unique.sam" > "${parent}.k31.unique.bed"
done
"$SAMTOOLS" faidx "$BT_FASTA"
cut -f1,2 "${BT_FASTA}.fai" > BT.genome
"$BEDTOOLS" makewindows -g BT.genome -w "$WIN" > win."$WIN".bed
"$BEDTOOLS" coverage -a win."$WIN".bed -b BSB.k31.unique.bed -counts > win.BSB.counts.tsv
"$BEDTOOLS" coverage -a win."$WIN".bed -b TC.k31.unique.bed -counts > win.TC.counts.tsv
# Window-level parental fractions and calls (low-support windows are indistinct).
awk -v min="$MIN_HITS" -v strong="$STRONG" 'BEGIN{OFS="\t"; print "chr","start","end","BSB","TC","total","pBSB","pTC","win_call"}
NR==FNR{b[NR]=$4; c[NR]=$1 FS $2 FS $3; next}
{n=NR-1; split(c[n],x,FS); total=b[n]+$4; pb=total?b[n]/total:0; pt=total?$4/total:0;
call=total<min?"Indistinct":(pb>=strong?"BSB_source":(pt>=strong?"TC_source":"Mixed"));
print x[1],x[2],x[3],b[n],$4,total,pb,pt,call}' win.BSB.counts.tsv win.TC.counts.tsv > BT.k31.win1000000.source.tsv

# Chromosome-level summaries use the same thresholds as the original workflow.
"$PYTHON" - <<'PY'
import math
from collections import Counter, defaultdict
from pathlib import Path

source = Path("BT.k31.win1000000.source.tsv")
out = Path("BT.k31.win1000000.chr_summary.tsv")
chr_stats = defaultdict(lambda: {"sum_pBSB": 0.0, "sum_pTC": 0.0, "n_used": 0,
                                  "BSB_source": 0, "TC_source": 0, "Mixed": 0})
win_counter = Counter()
with source.open() as f:
    header = next(f)
    for line in f:
        chrom, start, end, bsb, tc, total, p_bsb, p_tc, call = line.rstrip("\n").split("\t")
        total = int(total)
        win_counter[call] += 1
        if total >= 50:
            st = chr_stats[chrom]
            st["n_used"] += 1
            st["sum_pBSB"] += float(p_bsb)
            st["sum_pTC"] += float(p_tc)
            st[call] += 1

with out.open("w") as o:
    o.write("\t".join(["chr", "n_windows_used", "mean_pBSB", "mean_pTC",
                       "frac_BSBsource", "frac_TCsource", "frac_Mixed", "chr_call"]) + "\n")
    chr_counter = Counter()
    for chrom in sorted(chr_stats):
        st = chr_stats[chrom]
        n = st["n_used"]
        mean_b = st["sum_pBSB"] / n if n else float("nan")
        mean_t = st["sum_pTC"] / n if n else float("nan")
        frac_b = st["BSB_source"] / n if n else 0.0
        frac_t = st["TC_source"] / n if n else 0.0
        frac_m = st["Mixed"] / n if n else 0.0
        if n == 0:
            call = "IndistinctChr"
        elif mean_b >= 0.70 and frac_t < 0.05:
            call = "BSBChr"
        elif mean_t >= 0.70 and frac_b < 0.05:
            call = "TCChr"
        elif (frac_b > 0.05 and frac_t > 0.05) or frac_m > 0.10:
            call = "RecChr"
        else:
            call = "IndistinctChr"
        chr_counter[call] += 1
        o.write("\t".join([chrom, str(n),
                           "nan" if math.isnan(mean_b) else f"{mean_b:.6f}",
                           "nan" if math.isnan(mean_t) else f"{mean_t:.6f}",
                           f"{frac_b:.6f}", f"{frac_t:.6f}", f"{frac_m:.6f}", call]) + "\n")

print("[INFO] Window calls:", dict(win_counter))
print("[INFO] Chromosome calls:", dict(chr_counter))
print(f"[DONE] {out}")
PY
