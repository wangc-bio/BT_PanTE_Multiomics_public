#!/usr/bin/env bash
#SBATCH --job-name=03_validate_5Mb
#SBATCH --cpus-per-task=4
#SBATCH --mem=160G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
: "${BT_FASTA:?}"; : "${BSB_BED:?}"; : "${TC_BED:?}"; : "${BEDTOOLS:?}"; : "${SAMTOOLS:?}"; : "${OUTDIR:?}"
WIN=5000000; MIN_HITS=50; STRONG=0.70
mkdir -p "$OUTDIR" logs; cd "$OUTDIR"
"$SAMTOOLS" faidx "$BT_FASTA"
cut -f1,2 "${BT_FASTA}.fai" | awk '$1 ~ /^chr[0-9][0-9][AB]$/' > BT.primary.genome
"$BEDTOOLS" makewindows -g BT.primary.genome -w "$WIN" > win."$WIN".bed
"$BEDTOOLS" sort -g BT.primary.genome -i win."$WIN".bed > win."$WIN".sorted.bed
"$BEDTOOLS" coverage -a win."$WIN".sorted.bed -b "$BSB_BED" -counts > win.BSB.counts.tsv
"$BEDTOOLS" coverage -a win."$WIN".sorted.bed -b "$TC_BED" -counts > win.TC.counts.tsv
OUT_WIN="BT.k31.win5000000.onlyChr.source.tsv"
OUT_CHR="BT.k31.win5000000.onlyChr.chr_summary.tsv"
paste win.BSB.counts.tsv win.TC.counts.tsv | awk -v MIN_HITS="$MIN_HITS" -v STRONG="$STRONG" -v genome="BT.primary.genome" -v out_win="$OUT_WIN" -v out_chr="$OUT_CHR" '
BEGIN {
  OFS="\\t"
  while ((getline < genome) > 0) all[$1]=1
  close(genome)
  print "chr","start","end","BSB","TC","total","pBSB","pTC","win_call" > out_win
}
{
  if ($1!=$5 || $2!=$6 || $3!=$7) {
    print "[ERROR] window coordinate mismatch:",$1,$2,$3,"vs",$5,$6,$7 > "/dev/stderr"
    exit 1
  }
  chr=$1; start=$2; end=$3; bsb=$4; tc=$8
  total=bsb+tc
  pBSB=(total>0 ? bsb/total : 0)
  pTC=(total>0 ? tc/total : 0)
  if (total<MIN_HITS) call="Indistinct"
  else if (pBSB>=STRONG) call="BSB_source"
  else if (pTC>=STRONG) call="TC_source"
  else call="Mixed"
  printf "%s\\t%s\\t%s\\t%d\\t%d\\t%d\\t%.6f\\t%.6f\\t%s\\n",chr,start,end,bsb,tc,total,pBSB,pTC,call >> out_win
  win_call[call]++
  if (total>=MIN_HITS) {
    n_used[chr]++; sumB[chr]+=pBSB; sumT[chr]+=pTC
    if (call=="BSB_source") cB[chr]++
    else if (call=="TC_source") cT[chr]++
    else if (call=="Mixed") cM[chr]++
  }
}
END {
  print "chr","n_windows_used","mean_pBSB","mean_pTC","frac_BSBsource","frac_TCsource","frac_Mixed","chr_call" > out_chr
  for (chr in all) {
    nu=(chr in n_used ? n_used[chr] : 0)
    if (nu==0) {
      meanB="nan"; meanT="nan"; fracB=0; fracT=0; fracM=0
      chr_call="IndistinctChr"
    } else {
      meanB=sumB[chr]/nu; meanT=sumT[chr]/nu
      fracB=(chr in cB ? cB[chr]/nu : 0)
      fracT=(chr in cT ? cT[chr]/nu : 0)
      fracM=(chr in cM ? cM[chr]/nu : 0)
      if (meanB>=0.70 && fracT<0.05) chr_call="BSBChr"
      else if (meanT>=0.70 && fracB<0.05) chr_call="TCChr"
      else if ((fracB>0.05 && fracT>0.05) || fracM>0.10) chr_call="RecChr"
      else chr_call="IndistinctChr"
    }
    chr_call_cnt[chr_call]++
    if (meanB=="nan")
      printf "%s\\t%d\\tnan\\tnan\\t%.6f\\t%.6f\\t%.6f\\t%s\\n",chr,nu,fracB,fracT,fracM,chr_call >> out_chr
    else
      printf "%s\\t%d\\t%.6f\\t%.6f\\t%.6f\\t%.6f\\t%.6f\\t%s\\n",chr,nu,meanB,meanT,fracB,fracT,fracM,chr_call >> out_chr
  }
  print "[INFO] 5-Mb window calls:" > "/dev/stderr"
  for (k in win_call) print k,win_call[k] > "/dev/stderr"
  print "[INFO] 5-Mb chromosome calls:" > "/dev/stderr"
  for (k in chr_call_cnt) print k,chr_call_cnt[k] > "/dev/stderr"
}
'
