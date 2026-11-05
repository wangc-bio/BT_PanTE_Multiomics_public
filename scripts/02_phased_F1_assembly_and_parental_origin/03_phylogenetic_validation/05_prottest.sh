#!/usr/bin/env bash
#SBATCH --job-name=05_prottest
#SBATCH --cpus-per-task=36
#SBATCH --mem=100G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
: "${PROTTEST_JAR:?Path to ProtTest 3.4 jar}"; : "${ALIGNMENT_DIR:?Directory with chr01–chr24 supermatrices}"
mkdir -p logs
for n in $(seq -w 1 24); do
  chr=chr$n; d="$ALIGNMENT_DIR/$chr"; cd "$d"
  java -jar "$PROTTEST_JAR" -i all.phy -all-distributions -F -AIC -BIC -tc 0.5 -threads 36 -o prottest.out
done
