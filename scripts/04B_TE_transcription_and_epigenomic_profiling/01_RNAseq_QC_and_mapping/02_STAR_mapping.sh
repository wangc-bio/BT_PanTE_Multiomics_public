#!/usr/bin/env bash
set -euo pipefail

: "${STAR:?Set STAR executable (v2.7.10b)}"
: "${STAR_INDEX:?Set STAR genome index}"
: "${R1:?Set filtered R1 FASTQ.gz}"
: "${R2:?Set filtered R2 FASTQ.gz}"
: "${PREFIX:?Set output prefix}"

# STAR v2.7.10b. The two multimap parameters are the key settings
# retained for downstream Pan-TE quantification.
"$STAR" \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$R1" "$R2" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$PREFIX" \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterMultimapNmax 200 \
  --winAnchorMultimapNmax 200
