#!/usr/bin/env bash
# Scan GQ thresholds: filter query VCF with fastfilter (GQ cutoff only), run RTG vcfeval per (sample × GQ),
# append metrics to a TSV (overall + SNV / insertion / deletion from tp/fp/fn VCFs), and run R plots.
#
# Requirements: bcftools, tabix, java, python3; R + ggplot2 for plotting.
#
# Defaults: if --region is omitted, uses chr20:10000000-11000000.
# Use --sample NAME for one cohort, or --all-samples to evaluate every sample present in
# both baseline and query VCFs (intersection of sample names).
#
# Example:
#   ./scripts/gq_dp_vcfeval_scan.sh \
#     --baseline dipcall_merged_chr20.norm.vcf.gz \
#     --query subset_chr20.norm_sorted.chr20_10000000_20000000.g.vcf.gz \
#     --template /path/to/ref.sdf \
#     --evaluation-regions HG00097_hap2_hprc_r2_v1.0.1.dipcall.bed \
#     --all-samples \
#     --gq-min 0 --gq-max 25 \
#     --out-dir /tmp/gq_vcfeval_scan
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PARSE_PY="${ROOT}/scripts/parse_vcfeval_summary.py"
STRAT_PY="${ROOT}/scripts/vcfeval_stratified_pr.py"
PLOT_R="${ROOT}/scripts/plot_gq_dp_vcfeval.R"
PLOT_PR_R="${ROOT}/scripts/plot_gq_dp_vcfeval_pr.R"

BASELINE=""
QUERY=""
TEMPLATE=""
OUT_DIR=""
RTG_JAR="${RTG_JAR:-${ROOT}/scratch/rtg-tools-3.12.1/RTG.jar}"
FASTFILTER_BIN="${FASTFILTER_BIN:-cargo run -q --}"
REGION=""
REGION_USER_SET=0
EVAL_BED=""
BED_REGIONS=""
SAMPLE=""
ALL_SAMPLES=0
SAMPLES=()
GQ_MIN=0
GQ_MAX=25
DRY_RUN=0
SKIP_PLOT=0
EXTRA_VCFEVAL=()

die() { echo "error: $*" >&2; exit 1; }

usage() {
  sed -n '1,35p' "$0" | tail -n +2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --baseline) BASELINE="$2"; shift 2 ;;
    --query) QUERY="$2"; shift 2 ;;
    --template|-t) TEMPLATE="$2"; shift 2 ;;
    --out-dir|-o) OUT_DIR="$2"; shift 2 ;;
    --rtg-jar) RTG_JAR="$2"; shift 2 ;;
    --fastfilter-bin) FASTFILTER_BIN="$2"; shift 2 ;;
    --region) REGION="$2"; REGION_USER_SET=1; shift 2 ;;
    --full-genome) REGION=""; REGION_USER_SET=1; shift ;;
    --evaluation-regions|--confident-bed) EVAL_BED="$2"; shift 2 ;;
    --bed-regions) BED_REGIONS="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --all-samples) ALL_SAMPLES=1; shift ;;
    --gq-min) GQ_MIN="$2"; shift 2 ;;
    --gq-max) GQ_MAX="$2"; shift 2 ;;
    --dp-min|--dp-max)
      echo "warning: DP grid removed; ignoring $1 (GQ-only scan)" >&2
      shift 2
      ;;
    --dry-run) DRY_RUN=1; shift ;;
    --skip-plot) SKIP_PLOT=1; shift ;;
    --vcfeval-arg) EXTRA_VCFEVAL+=("$2"); shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown arg: $1" ;;
  esac
done

[[ -n "$BASELINE" && -n "$QUERY" && -n "$TEMPLATE" && -n "$OUT_DIR" ]] || die "require --baseline --query --template --out-dir"

if [[ $ALL_SAMPLES -eq 1 && -n "$SAMPLE" ]]; then
  die "use either --all-samples or --sample, not both"
fi
if [[ $ALL_SAMPLES -eq 0 && -z "$SAMPLE" ]]; then
  die "require --sample NAME or --all-samples"
fi

if [[ $REGION_USER_SET -eq 0 ]]; then
  REGION="chr20:10000000-11000000"
fi

command -v bcftools >/dev/null || die "bcftools not on PATH"
command -v tabix >/dev/null || die "tabix not on PATH"
command -v java >/dev/null || die "java not on PATH"
command -v python3 >/dev/null || die "python3 not on PATH"
[[ -f "$RTG_JAR" ]] || die "RTG jar not found: $RTG_JAR"
[[ -f "$PARSE_PY" ]] || die "missing $PARSE_PY"
[[ -f "$STRAT_PY" ]] || die "missing $STRAT_PY"
if [[ $DRY_RUN -eq 0 ]]; then
  [[ -e "$TEMPLATE" ]] || die "template SDF not found: $TEMPLATE"
fi

if [[ $ALL_SAMPLES -eq 1 ]]; then
  if [[ $DRY_RUN -eq 1 ]]; then
    SAMPLES=("dry_run_sample")
  else
    SAMPLES=()
    while IFS= read -r line; do
      [[ -n "$line" ]] && SAMPLES+=("$line")
    done < <(comm -12 <(bcftools query -l "$BASELINE" | sort) <(bcftools query -l "$QUERY" | sort))
    [[ ${#SAMPLES[@]} -gt 0 ]] || die "no shared sample names between baseline and query VCFs"
  fi
else
  SAMPLES=("$SAMPLE")
fi

mkdir -p "$OUT_DIR"
TSV="${OUT_DIR}/metrics.tsv"
LOG="${OUT_DIR}/run.log"

{
  echo "# gq_vcfeval_scan $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "# baseline=${BASELINE}"
  echo "# query=${QUERY}"
  echo "# template=${TEMPLATE}"
  echo "# region=${REGION:-<none>}"
  echo "# evaluation_regions=${EVAL_BED:-<none>}"
  echo "# samples=${SAMPLES[*]} (${#SAMPLES[@]} names)"
  echo "# gq ${GQ_MIN}..${GQ_MAX} (GQ-only; no DP threshold)"
} | tee -a "$LOG"

echo -e "sample\tgq_thr\tprecision\trecall\tf1\tprecision_snv\trecall_snv\tf1_snv\tprecision_ins\trecall_ins\tf1_ins\tprecision_del\trecall_del\tf1_del\tvcfeval_dir\tstatus" | tee "$TSV"

run_one_gq() {
  local gq="$1"
  local tag="gq${gq}"
  local fdir="${OUT_DIR}/filtered/${tag}"
  local filtered="${fdir}/calls_filtered.vcf.gz"

  mkdir -p "$fdir"

  if [[ $DRY_RUN -eq 1 ]]; then
    echo "[dry-run] filter -> $filtered ; vcfeval -> ${OUT_DIR}/vcfeval/${tag}/<sample>"
    return 0
  fi

  echo "## fastfilter ${tag} (min_gq=${gq} only)" >>"$LOG"
  if [[ "$FASTFILTER_BIN" == "cargo run -q --" ]]; then
    (cd "$ROOT" && cargo run -q -- filter --input "$QUERY" --output "$filtered" ${REGION:+--region "$REGION"} \
      thresholds --min-gq "$gq") >>"$LOG" 2>&1
  else
    # shellcheck disable=SC2086
    "$FASTFILTER_BIN" filter --input "$QUERY" --output "$filtered" ${REGION:+--region "$REGION"} \
      thresholds --min-gq "$gq" >>"$LOG" 2>&1
  fi

  rm -f "${filtered}.csi" "${filtered}.tbi"
  tabix -fp vcf "$filtered" >>"$LOG" 2>&1

  local s
  for s in "${SAMPLES[@]}"; do
    local edir="${OUT_DIR}/vcfeval/${tag}/${s}"
    rm -rf "$edir"

    local vcfe=( java -jar "$RTG_JAR" vcfeval -b "$BASELINE" -c "$filtered" -o "$edir" -t "$TEMPLATE" )
    if [[ -n "$EVAL_BED" ]]; then
      vcfe+=( --evaluation-regions "$EVAL_BED" )
    fi
    if [[ -n "$BED_REGIONS" ]]; then
      vcfe+=( --bed-regions "$BED_REGIONS" )
    fi
    if [[ -n "$REGION" ]]; then
      vcfe+=( --region "$REGION" )
    fi
    vcfe+=( --sample "$s" )
    if [[ ${#EXTRA_VCFEVAL[@]} -gt 0 ]]; then
      vcfe+=( "${EXTRA_VCFEVAL[@]}" )
    fi

    echo "## vcfeval ${tag} sample=${s}" >>"$LOG"
    if ! "${vcfe[@]}" >>"$LOG" 2>&1; then
      # sample, gq, 12 empty metric cols, edir, status
      echo -e "${s}\t${gq}\t\t\t\t\t\t\t\t\t\t\t\t${edir}\tvcfeval_fail" | tee -a "$TSV"
      continue
    fi

    local summary="${edir}/summary.txt"
    local pr rf f1
    if ! read -r pr rf f1 < <(python3 "$PARSE_PY" "$summary"); then
      pr=""; rf=""; f1=""
    fi
    local st="ok"
    if [[ -z "$pr" ]]; then
      st="parse_fail"
    fi
    local strat
    strat=$(python3 "$STRAT_PY" "$edir" || true)
    echo -e "${s}\t${gq}\t${pr}\t${rf}\t${f1}\t${strat}\t${edir}\t${st}" | tee -a "$TSV"
  done
}

for gq in $(seq "$GQ_MIN" "$GQ_MAX"); do
  run_one_gq "$gq"
done

if [[ $SKIP_PLOT -eq 0 && $DRY_RUN -eq 0 ]]; then
  if command -v Rscript >/dev/null; then
    Rscript "$PLOT_R" "$TSV" "${OUT_DIR}/plots" >>"$LOG" 2>&1 || echo "warning: R metric plot failed (see $LOG)" | tee -a "$LOG"
    Rscript "$PLOT_PR_R" "$TSV" "${OUT_DIR}/plots" >>"$LOG" 2>&1 || echo "warning: R PR plot failed (see $LOG)" | tee -a "$LOG"
  else
    echo "Rscript not found; skip plot. TSV: $TSV" | tee -a "$LOG"
  fi
fi

echo "Done. Metrics: $TSV" | tee -a "$LOG"
