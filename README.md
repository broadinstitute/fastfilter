# fastfilter

`fastfilter` is a Rust console tool to mask low-confidence genotypes in a joint-called DeepVariant + GLNexus callset.

It operates on multi-sample VCF/BCF input and writes a filtered VCF/BCF output.
Only variants with record-level `FILTER=PASS` are modified. For failing genotypes, it sets `GT` to missing (`./.`) while keeping all other FORMAT fields unchanged.

## Build

```bash
cargo build --release
```

## Thresholds mode (GQ/DP)

```bash
fastfilter filter \
  --input in.vcf.gz \
  --output out.vcf.gz \
  thresholds \
  --min-gq 20 \
  --min-dp 10
```

Semantics:
- For each sample at a `FILTER=PASS` site:
  - if `GQ < --min-gq` or `DP < --min-dp` (or the value is missing/unknown), the tool masks that genotype by setting `GT` to `./.`.
- Variants with `FILTER != PASS` are passed through unchanged.

## gbdt-rs mode (model-based per-alt scoring)

This mode uses three pre-converted GBDT models (one per variant class):
- SNV
- insertion (alt_len > ref_len)
- deletion (alt_len < ref_len)

Because `gbdt-rs` loads models in its own format, XGBoost models must be converted first using the `convert_xgboost.py` helper shipped with `gbdt-rs`:
- `convert_xgboost.py` usage is described in the upstream repo README: [mesalock-linux/gbdt-rs](https://github.com/mesalock-linux/gbdt-rs)
- the script logic is here: [examples/convert_xgboost.py](https://raw.githubusercontent.com/mesalock-linux/gbdt-rs/master/examples/convert_xgboost.py)

Your notebook’s `predict_proba[:, 1]` score corresponds to `gbdt-rs`’s binary classification output for `binary:logistic` (loaded via `GBDT::from_xgboost_dump(..., "binary:logistic")`).

CLI:

```bash
fastfilter filter \
  --input in.vcf.gz \
  --output out.vcf.gz \
  gbdt \
  --model-threshold 0.9 \
  --gbdt-objective binary:logistic \
  --gbdt-model-snv snv.gbdt.model \
  --gbdt-model-insertion insertion.gbdt.model \
  --gbdt-model-deletion deletion.gbdt.model
```

Masking semantics:
- Only variants with `FILTER=PASS` are processed.
- For each sample genotype at such a site:
  - the code scores each alt allele instance in the genotype with the model for that variant class,
  - then masks the genotype if `P(correct call)` for **any** alt instance is `< --model-threshold`.
- Masking is done by setting `GT` to missing (`./.`). All other FORMAT fields are kept.

## Optional region processing

You can restrict processing to a genomic interval via `--region`:

```bash
fastfilter filter --input in.vcf.gz --output out.vcf.gz --region "chr20:20000000-21000000" ...
```

This requires an indexed input and uses `tabix` (no streaming fallback).

## Evaluation (truth vs query)

`fastfilter` also provides a lightweight evaluator for SNVs + small indels, intended for fast iteration against assembly-derived truth sets (e.g., dipcall output).

```bash
fastfilter eval \
  --truth truth.vcf.gz \
  --query query.vcf.gz \
  --confident-bed HG00097_hap2_hprc_r2_v1.0.1.dipcall.bed \
  --region "chr20:10000000-11000000"
```

Notes:
- `--confident-bed` is **required**: only variants whose POS lies in these intervals (0-based half-open BED) are evaluated.
- **`--eval-sample NAME`**: repeat or use comma-separated names to evaluate several samples in one run. The same name must exist in both truth and query VCFs. Incompatible with `--truth-sample` / `--query-sample`. With multiple samples you get per-sample TSV sections, a **pooled** summary, pooled `gq_scan` / `gq_dp_scan`, and **JSON** shaped as `{ "pooled": ..., "per_sample": { ... } }`.
- If `--region` is set, evaluation uses the **intersection** of that interval with the confident BED (both must contain the variant POS).
- `--region` is a streaming filter (does not require indexes).
- Matching is **exact on (CHROM, POS, REF, ALT)** after trimming shared REF/ALT prefix/suffix (no reference-based left-alignment).
- Metrics are stratified into **snv / insertion / deletion** based on normalized REF/ALT lengths.
- Insertions and deletions are additionally stratified by **indel length** (\(|len(ALT)-len(REF)|\)).
- By default, `eval` prints **genotype-based** tables only. Add `--emit-allele` to also print allele TP/FP/FN tables.
  - Genotype output is a single table keyed by `(type, len)`, where SNVs have `len=0`.
- You can also grid-scan query GQ thresholds (useful for PR-style plots):

```bash
fastfilter eval \
  --truth truth.vcf.gz \
  --query query.vcf.gz \
  --confident-bed HG00097_hap2_hprc_r2_v1.0.1.dipcall.bed \
  --query-sample HG00097 \
  --scan-gq 0 99 \
  --plot-gq-pr-png gq_pr.png
```

The PR plot shows separate precision-vs-recall curves for `snv`, `insertion`, and `deletion` (plus a faint `total` curve). Curves are computed only on variants inside **`--confident-bed`** (and inside **`--region`**, if set). The PNG title lists the confident BED file name and optional region so figures stay self-describing.

With **multiple `--eval-sample`**, the same `--plot-gq-pr-png` path writes a **2×2 grid** (total / SNV / insertion / deletion). In each panel, **one translucent PR polyline per sample** connects points in **increasing GQ threshold** order; recall and precision axes are auto-padded from the pooled points.

**GQ × DP grid scan:** add `--scan-dp MIN MAX` (requires `--scan-gq`) to evaluate every combination where the query “passes” only if **GQ ≥ gq_thr and DP ≥ dp_thr** (per-sample `FORMAT` tags, default `--gq-tag GQ`, `--dp-tag DP`). This emits a `gq_dp_scan` TSV (all variant classes). **`--plot-gq-pr-png` cannot be used with `--scan-dp`** (that flag is the 1D GQ-only PR curves).

For 2D scans you can plot either or both:

- **`--plot-gq-dp-f1-png`:** F1 heatmap over the **(GQ_thr, DP_thr)** grid (**total** stratum). Large red regions usually mean F1 collapsed (often recall → 0 under strict cutoffs), which is expected.
- **`--plot-gq-dp-pr-png`:** With **one sample**, a **scatter** in recall–precision space (one point per `(GQ≥…, DP≥…)` cell on the **total** stratum), colored by F1. With **multiple `--eval-sample`**, a **2×2 panel** figure (total / SNV / insertion / deletion): **one translucent polyline per sample** in recall–precision space, vertices in **row-major GQ-then-DP** threshold order; the precision axis uses **`--gq-dp-pr-precision-min/max`** (default `0.9–1.0`).

```bash
fastfilter eval \
  --truth truth.vcf.gz \
  --query query.vcf.gz \
  --confident-bed HG00097_hap2_hprc_r2_v1.0.1.dipcall.bed \
  --query-sample HG00097 \
  --scan-gq 0 40 \
  --scan-dp 5 40 \
  --plot-gq-dp-f1-png gq_dp_f1.png \
  --plot-gq-dp-pr-png gq_dp_pr.png
```

Multi-sample example (e.g. merged dipcall truth vs a multi-sample query on chr20), with overlaid GQ×DP PR curves:

```bash
fastfilter eval \
  --truth dipcall_merged_chr20.vcf \
  --query joint_chr20.vcf.gz \
  --confident-bed confident_regions.bed \
  --region "chr20:10000000-20000000" \
  --eval-sample HG00097,HG00099,HG00100 \
  --scan-gq 0 40 \
  --scan-dp 5 40 \
  --plot-gq-dp-f1-png gq_dp_f1_pooled.png \
  --plot-gq-dp-pr-png gq_dp_pr_multisample.png
```

Use a **per-sample confident BED** by running once per sample if your BEDs differ; the same `--confident-bed` is applied to every sample in the multi-sample run.

Example with both a sub-region and a PR figure:

```bash
fastfilter eval \
  --truth truth.vcf.gz \
  --query query.vcf.gz \
  --confident-bed HG00097_hap2_hprc_r2_v1.0.1.dipcall.bed \
  --region "chr20:10000000-11000000" \
  --query-sample HG00097 \
  --scan-gq 0 99 \
  --plot-gq-pr-png chr20_gq_pr.png
```

## Feature expectations (model input columns)

The implementation uses a fixed feature vector per alt-allele instance, derived from FORMAT/INFO fields:
- `DP`, `GQ`, `AB`, `alt_AD`, `ref_AD`, `PL_ref`, `PL_het`, `PL_hom`
- `RNC_bad`
- `site_AF`, `site_AQ`
- `is_snv`, `is_het`, `is_hom_alt`

Before inference, missing/non-finite numeric feature values are filled with `-1.0`, matching the preprocessing described in your notebook.