#!/usr/bin/env bash
set -euo pipefail

: "${R1:?Set R1}"
: "${R2:?Set R2}"
: "${OUT_R1:?Set OUT_R1}"
: "${OUT_R2:?Set OUT_R2}"
: "${HTML:?Set HTML report path}"
: "${JSON:?Set JSON report path}"

# fastp v0.23.2; default filtering settings.
fastp \
  -i "$R1" -I "$R2" \
  -o "$OUT_R1" -O "$OUT_R2" \
  -h "$HTML" -j "$JSON"
