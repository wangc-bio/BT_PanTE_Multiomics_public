#!/usr/bin/env bash
# Public-release note: source script from the executed analysis. Replace the original local absolute paths and executable locations with paths in your environment before running.
#SBATCH --job-name=02_ATAC_TC_array
#SBATCH --partition=pms
#SBATCH --array=1-9
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

############################
# ĹäÖĂ
############################
WORKDIR="${FB_ROOT:?Set FB_ROOT to the parent FB directory}/ATAC/TC"
MANIFEST="${WORKDIR}/samples.tsv"
REF="${WORKDIR}/01_index/TC.genome.fa"
GSIZE="1052771736"            # ˝¨Ňé»»ÎŞÓĐĐ§»ůŇň×é´óĐˇŁ¨·ÇNĽî»ů×ÜÁżŁ©
PICARD_JAR="${PICARD_JAR:?Set PICARD_JAR to the Picard jar file}"

BLACKLIST=""                  # ŔýŁş/path/to/blacklist.bedŁ¨ÁôżŐÔňĚřąýŁ©
ORG_REGEX=""                  # ŔýŁş'(^chrM$|^MT$)'Ł¨ÁôżŐÔň˛»ČĄ organellarŁ©
TSS_BED="${WORKDIR}/ref/TSS.bed"   # Čô´ćÔÚÔňĽĆËă TSS enrichment

THREADS="${SLURM_CPUS_PER_TASK:-12}"

# żÉ°´ĐčĽÓÔŘ»·ľłŁ¨ČôĂüÁî˛»ÔÚ PATHŁ©
# Activate the required tools in your local job wrapper if they are not already on PATH.
# conda activate macs2
# conda activate deeptools

############################
# ÄżÂĽ
############################
OUT_ALIGN="${WORKDIR}/03_align"
OUT_DEDUP="${WORKDIR}/04_dedup"
OUT_CLEAN="${WORKDIR}/05_clean"
OUT_PEAKS="${WORKDIR}/06_peaks"
OUT_BIGWIG="${WORKDIR}/07_bigwig"
OUT_QC="${WORKDIR}/08_qc"
LOGDIR="${WORKDIR}/logs"
TMPDIR="${WORKDIR}/tmp"
mkdir -p "$OUT_ALIGN" "$OUT_DEDUP" "$OUT_CLEAN" "$OUT_PEAKS" "$OUT_BIGWIG" "$OUT_QC" "$LOGDIR" "$TMPDIR"

############################
# ČˇŃů±ľ
############################
[[ -f "$MANIFEST" ]] || { echo "ERROR: Î´ŐŇµ˝ $MANIFEST"; exit 1; }
# Ë÷ŇýÓ¦ŇŃ˝¨Á˘Ł»ČôÎ´˝¨Á˘Ôň±¨´íŁ¨˛»ÖŘ˝¨Ł©
if [[ ! -f "${REF}.bwt" && ! -f "${REF}.0123" ]]; then
  echo "ERROR: Î´Ľě˛âµ˝ ${REF} µÄ BWA Ë÷ŇýŁ¨*.bwt »ň *.0123Ł©ˇŁÇëĎČ˝¨Á˘Ë÷ŇýˇŁ"; exit 1;
fi

IDX="${SLURM_ARRAY_TASK_ID:?ĐčŇŞ SLURM_ARRAY_TASK_ID}"
read -r _hdr < "$MANIFEST"   # Ěřąý±íÍ·
line=$(awk -v i="$IDX" 'NR==i+1{print; exit}' "$MANIFEST") || { echo "Ńů±ľË÷Ňýł¬˝ç: $IDX"; exit 1; }
sample=$(echo "$line" | cut -f1)
R1=$(echo "$line" | cut -f2)
R2=$(echo "$line" | cut -f3)

echo "==> [$(date)] Task ${IDX} Ńů±ľ: ${sample}"

BAM_SORT="${OUT_ALIGN}/${sample}.sorted.bam"
BAM_MKDUP="${OUT_DEDUP}/${sample}.mkdup.bam"
DUP_METRICS="${OUT_DEDUP}/${sample}.dup.metrics.txt"
BAM_CLEAN="${OUT_CLEAN}/${sample}.clean.bam"
PEAK_PREFIX="${OUT_PEAKS}/${sample}"
BIGWIG_CPM="${OUT_BIGWIG}/${sample}.CPM.bw"

############################
# 1) ¶ÔĆë + ĹĹĐň
############################
if [[ ! -f "${BAM_SORT}" ]]; then
  echo "[`date`] BWA-MEM: ${sample}"
  ( set -x
    bwa mem -t ${THREADS} -M \
      -R "@RG\tID:${sample}\tLB:${sample}\tSM:${sample}\tPL:ILLUMINA" \
      "$REF" "$R1" "$R2" \
    | samtools sort -@ ${THREADS} -o "${BAM_SORT}" -
  ) &> "${LOGDIR}/${sample}.align.log"
  samtools index "${BAM_SORT}"
fi

############################
# 2) ČĄÖŘ¸´
############################
if [[ ! -f "${BAM_MKDUP}" ]]; then
  echo "[`date`] Picard MarkDuplicates: ${sample}"
  ( set -x
    java -Xmx20g -Djava.io.tmpdir="${TMPDIR}" -jar "${PICARD_JAR}" MarkDuplicates \
      I="${BAM_SORT}" O="${BAM_MKDUP}" M="${DUP_METRICS}" \
      REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
  ) &> "${LOGDIR}/${sample}.markdup.log"
  samtools index "${BAM_MKDUP}"
fi

############################
# 3) ŃĎ¸ńąýÂËŁ¨MAPQˇÝ30 + proper pairŁ»ČĄ unmapped/´ÎĽ¶/˛ąłä/QCĘ§°Ü/dupŁ©
############################
if [[ ! -f "${BAM_CLEAN}" ]]; then
  echo "[`date`] ąýÂË: ${sample}"
  ( set -x
    # 1804 = 0x4|0x8|0x100|0x200|0x400
    samtools view -@ ${THREADS} -b -q 30 -f 2 -F 1804 "${BAM_MKDUP}" > "${BAM_CLEAN}"
  ) &> "${LOGDIR}/${sample}.filter.log"

  if [[ -n "${ORG_REGEX}" ]]; then
    keep_list="${OUT_CLEAN}/${sample}.keep.contigs.txt"
    samtools idxstats "${BAM_CLEAN}" | cut -f1 | grep -Ev "${ORG_REGEX}" > "${keep_list}"
    samtools view -@ ${THREADS} -b "${BAM_CLEAN}" $(cat "${keep_list}") > "${OUT_CLEAN}/${sample}.clean.noOrg.bam"
    mv "${OUT_CLEAN}/${sample}.clean.noOrg.bam" "${BAM_CLEAN}"
    rm -f "${keep_list}"
  fi
  samtools index "${BAM_CLEAN}"
fi

############################
# 4) µ÷·ĺŁ¨BAMPEŁ©
############################
if [[ ! -f "${PEAK_PREFIX}_peaks.narrowPeak" ]]; then
  echo "[`date`] MACS2: ${sample}"
  ( set -x
    macs2 callpeak -t "${BAM_CLEAN}" -f BAMPE -g "${GSIZE}" -n "${sample}" \
      --outdir "${OUT_PEAKS}" --call-summits -q 0.01
  ) &> "${LOGDIR}/${sample}.macs2.log"
  if [[ -n "${BLACKLIST}" && -f "${BLACKLIST}" ]]; then
    bedtools intersect -v -a "${PEAK_PREFIX}_peaks.narrowPeak" -b "${BLACKLIST}" \
      > "${PEAK_PREFIX}.peaks.filt.bed"
  fi
fi

############################
# 5) bigWigŁ¨CPMŁ»ČçĐč RPGC ¸ÄÎŞ --normalizeUsing RPGC --effectiveGenomeSize ${GSIZE}Ł©
############################
if [[ ! -f "${BIGWIG_CPM}" ]]; then
  echo "[`date`] bigWig: ${sample}"
  ( set -x
    bamCoverage -b "${BAM_CLEAN}" -o "${BIGWIG_CPM}" --normalizeUsing CPM -p ${THREADS}
  ) &> "${LOGDIR}/${sample}.bigwig.log"
fi

############################
# 6) Ć¬¶Îł¤¶Č QC
############################
if [[ ! -f "${OUT_QC}/${sample}.fraglen.txt" ]]; then
  ( set -x
    bamPEFragmentSize -b "${BAM_CLEAN}" -T "${sample}" -hist "${OUT_QC}/${sample}.fraglen.txt" \
      --samplesLabel "${sample}" -p ${THREADS}
  ) &> "${LOGDIR}/${sample}.fraglen.log"
fi

############################
# 7) TSS ¸»ĽŻŁ¨ČôĚáą© TSS.bedŁ©
############################
if [[ -f "${TSS_BED}" && ! -f "${OUT_QC}/${sample}.TSS_enrich.txt" ]]; then
  MAT="${OUT_QC}/${sample}.TSS.mat.gz"
  TAB="${OUT_QC}/${sample}.TSS_profile.tab"
  PNG="${OUT_QC}/${sample}.TSS_profile.png"
  ( set -x
    computeMatrix reference-point -S "${BIGWIG_CPM}" -R "${TSS_BED}" \
      -a 2000 -b 2000 --binSize 10 --missingDataAsZero \
      -o "${MAT}" -p ${THREADS}
    plotProfile -m "${MAT}" -out "${PNG}" --perGroup --plotTitle "${sample}" \
      --dpi 150 --outFileNameData "${TAB}"
  ) &> "${LOGDIR}/${sample}.tss.log"
  awk -v s="${sample}" '
    BEGIN{maxv=0; n=0}
    /^#/ {next}
    {v=$2; arr[++n]=v; if(v>maxv)maxv=v}
    END{
      bg=0; k=0;
      for(i=1;i<=100 && i<=n;i++){bg+=arr[i];k++}
      for(i=n-99;i<=n;i++){if(i>=1){bg+=arr[i];k++}}
      bgavg=(k>0?bg/k:0);
      tss=(bgavg>0?maxv/bgavg:0);
      printf("%s\t%.4f\n", s, tss);
    }' "${TAB}" > "${OUT_QC}/${sample}.TSS_enrich.txt"
fi

############################
# 8) FRiPŁ¨·ĺÄÚ reads / ×Ü readsŁ©
############################
if [[ ! -f "${OUT_QC}/${sample}.FRiP.txt" ]]; then
  if [[ -f "${PEAK_PREFIX}.peaks.filt.bed" ]]; then
    PEAKBED="${PEAK_PREFIX}.peaks.filt.bed"
  else
    PEAKBED="${PEAK_PREFIX}_peaks.narrowPeak"
  fi
  total=$(samtools view -c "${BAM_CLEAN}")
  inpk=$(bedtools intersect -u -abam "${BAM_CLEAN}" -b "${PEAKBED}" | samtools view -c -)
  awk -v s="${sample}" -v a="${inpk}" -v b="${total}" 'BEGIN{printf "%s\t%.6f\t%d\t%d\n", s, (b>0?a/b:0), a, b}' \
    > "${OUT_QC}/${sample}.FRiP.txt"
fi

############################
# 9) Ěáą© MultiQC µÄÍłĽĆ
############################
[[ -f "${OUT_QC}/${sample}.flagstat.txt" ]] || samtools flagstat "${BAM_CLEAN}" > "${OUT_QC}/${sample}.flagstat.txt"
[[ -f "${OUT_QC}/${sample}.stats.txt"    ]] || samtools stats     "${BAM_CLEAN}" > "${OUT_QC}/${sample}.stats.txt"

echo "==> [$(date)] ÍęłÉ: ${sample}"