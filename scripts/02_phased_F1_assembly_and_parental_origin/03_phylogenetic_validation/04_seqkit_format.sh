#!/usr/bin/env bash
#SBATCH --job-name=04_seqkit_format
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --array=1-24%24
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -Eeuo pipefail

: "${SEQKIT:?Path to SeqKit executable}"
: "${GBLOCK_DIR:?Directory containing chromosome folders chr01–chr24}"
: "${OUT_BASE:?Output base directory}"
CHR=$(printf "chr%02d" "${SLURM_ARRAY_TASK_ID}")
BASE="$OUT_BASE/$CHR"
INDIR="$GBLOCK_DIR/$CHR"
D04="$BASE/04_seq"
D05="$BASE/05_merge"
D06="$BASE/06_out"
mkdir -p "$D04" "$D05" "$D06" logs

echo "[$(date '+%F %T')] $CHR | IN=$INDIR"
cd "$D04"

# Link Gblocks outputs, then sort each alignment and convert it to one sequence per line.
shopt -s nullglob
for f in "$INDIR"/*.muscle-gb; do
  ln -sf "$f" .
done
for i in *.muscle-gb; do
  [[ -e "$i" ]] || continue
  "$SEQKIT" sort "$i" > "$i.3"
  "$SEQKIT" seq "$i.3" -w 0 > "$i.3.4"
done

# Concatenate aligned sequences in matching order to form chromosome-level all.fa.
# After this automated SeqKit/paste step, manually standardize taxon sequence IDs,
# inspect the resulting matrix, and confirm it is suitable for downstream analysis.
if compgen -G "*.muscle-gb.3.4" > /dev/null; then
  paste -d " " *.muscle-gb.3.4 > "$D05/all.fa"
  sed -i 's/ //g' "$D05/all.fa"
  cp -f "$D05/all.fa" "$D06/"
  echo "[$(date '+%F %T')] #04-05_finished $CHR | merged -> $D05/all.fa"
else
  echo "[$(date '+%F %T')] WARN: $CHR has no *.muscle-gb.3.4 files to merge" >&2
  exit 1
fi
