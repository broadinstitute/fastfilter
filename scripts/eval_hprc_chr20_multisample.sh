#!/usr/bin/env bash
# Evaluate merged dipcall truth on chr20 vs joint calls for all samples present in BOTH VCFs.
#
# Confident regions: chr20 rows from HG00097_hap2 dipcall BED in this repo. That BED describes
# HG00097 hap2 confident intervals, not each HPRC sample’s own dipcall mask—substitute a union
# (or per-sample run) when you have sample-matched confident BEDs from the same HPRC release.
#
# Usage (from repo root):
#   ./scripts/eval_hprc_chr20_multisample.sh
# Optional env:
#   TRUTH=... QUERY=... SAMPLES=comma,list JSON_OUT=summary.json
#   OUT_F1=path.png OUT_PR=path.png
#   SCAN_GQ_MIN=0 SCAN_GQ_MAX=40 SCAN_DP_MIN=5 SCAN_DP_MAX=40
#   REGION=chr20:start-end  (required format; default full GRCh38 chr20)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

TRUTH="${TRUTH:-$ROOT/dipcall_merged_chr20.vcf}"
QUERY="${QUERY:-$ROOT/output.vcf.gz}"
SRC_BED="$ROOT/HG00097_hap2_hprc_r2_v1.0.1.dipcall.bed"

# All samples shared by truth and query (23 HPRC / 1kg-style IDs on chr20 merge vs joint callset).
SAMPLES="${SAMPLES:-HG00133,HG00323,HG00597,HG00735,HG01081,HG01433,HG01928,HG01960,HG02074,HG02258,HG02391,HG02809,HG03804,HG03831,HG03942,HG04160,HG06807,NA18522,NA18948,NA18976,NA18982,NA19468,NA20827}"

BED_TMP="${TMPDIR:-/tmp}/fastfilter_chr20_conf_$$.bed"
cleanup() { rm -f "$BED_TMP"; }
trap cleanup EXIT

if [[ ! -f "$TRUTH" ]]; then
  echo "Missing truth VCF: $TRUTH" >&2
  exit 1
fi
if [[ ! -f "$QUERY" ]]; then
  echo "Missing query VCF: $QUERY" >&2
  exit 1
fi
if [[ ! -f "$SRC_BED" ]]; then
  echo "Missing source BED (for chr20 subset): $SRC_BED" >&2
  exit 1
fi

grep -E '^chr20[[:space:]]' "$SRC_BED" > "$BED_TMP"
if [[ ! -s "$BED_TMP" ]]; then
  echo "No chr20 intervals found in $SRC_BED" >&2
  exit 1
fi

OUT_F1="${OUT_F1:-$ROOT/gq_dp_f1_pooled.png}"
OUT_PR="${OUT_PR:-$ROOT/gq_dp_pr_multisample.png}"
JSON_OUT="${JSON_OUT:-}"
SCAN_GQ_MIN="${SCAN_GQ_MIN:-0}"
SCAN_GQ_MAX="${SCAN_GQ_MAX:-40}"
SCAN_DP_MIN="${SCAN_DP_MIN:-5}"
SCAN_DP_MAX="${SCAN_DP_MAX:-40}"
# GRCh38 chr20 length (1-based inclusive end)
REGION="${REGION:-chr20:1-64444167}"

cmd=(
  cargo run -q --
  eval
  --truth "$TRUTH"
  --query "$QUERY"
  --confident-bed "$BED_TMP"
  --region "$REGION"
  --eval-sample "$SAMPLES"
  --scan-gq "$SCAN_GQ_MIN" "$SCAN_GQ_MAX"
  --scan-dp "$SCAN_DP_MIN" "$SCAN_DP_MAX"
  --plot-gq-dp-f1-png "$OUT_F1"
  --plot-gq-dp-pr-png "$OUT_PR"
)
if [[ -n "$JSON_OUT" ]]; then
  cmd+=(--json "$JSON_OUT")
fi

echo "Running: ${cmd[*]}" >&2
"${cmd[@]}"

echo >&2
echo "Wrote: $OUT_F1" >&2
echo "Wrote: $OUT_PR" >&2
[[ -n "$JSON_OUT" ]] && echo "Wrote: $JSON_OUT" >&2
