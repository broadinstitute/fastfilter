use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use anyhow::{anyhow, bail, Context, Result};
use clap::Args;
use plotters::prelude::*;
use plotters::series::LineSeries;
use plotters::style::RGBColor;
use rust_htslib::bcf;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::Read;
use serde::Serialize;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize)]
#[serde(rename_all = "lowercase")]
enum VarClass {
    Snv,
    Insertion,
    Deletion,
}

#[derive(Debug, Args, Clone)]
pub struct EvalArgs {
    /// Truth VCF/BCF (typically .vcf.gz)
    #[arg(long)]
    pub truth: PathBuf,

    /// Query VCF/BCF (typically .vcf.gz)
    #[arg(long)]
    pub query: PathBuf,

    /// Optional region (streaming filter): e.g. chr20:20000000-21000000.
    /// Variants must fall inside this interval **and** inside `--confident-bed`.
    #[arg(long)]
    pub region: Option<String>,

    /// Dipcall (or other) confident/callable regions BED (0-based, half-open).
    /// Required: evaluation is restricted to these intervals; if `--region` is set, only
    /// the intersection of that region with confident intervals is used.
    #[arg(long)]
    pub confident_bed: PathBuf,

    /// Truth sample name (default: the only sample if there is exactly one, else error).
    #[arg(long)]
    pub truth_sample: Option<String>,

    /// Query sample name to evaluate (default: first sample).
    #[arg(long)]
    pub query_sample: Option<String>,

    /// Evaluate multiple samples in one run: same sample name must exist in both truth and query VCFs.
    /// Comma-separated or repeat the flag (`--eval-sample A --eval-sample B`). Incompatible with `--truth-sample` / `--query-sample`.
    #[arg(long = "eval-sample", value_delimiter = ',')]
    pub eval_sample: Vec<String>,

    /// Require that truth and query sample names match exactly.
    #[arg(long, default_value_t = false)]
    pub require_sample_name_match: bool,

    /// If set, only count variants where the chosen sample is non-reference in truth/query.
    #[arg(long, default_value_t = true)]
    pub require_sample_nonref: bool,

    /// Treat missing genotype (./.) as "no call" (excluded from TP/FP/FN).
    #[arg(long, default_value_t = true)]
    pub exclude_nocall: bool,

    /// FORMAT tag for query genotype quality.
    #[arg(long, default_value = "GQ")]
    pub gq_tag: String,

    /// FORMAT tag for per-sample read depth (used with `--scan-dp`).
    #[arg(long, default_value = "DP")]
    pub dp_tag: String,

    /// Grid-scan integer GQ thresholds (inclusive min/max) on the query sample.
    /// Calls with GQ < threshold are treated as not-called for counting TP/FP/FN.
    /// When set, emits a per-threshold table `gq_scan` in addition to other outputs.
    #[arg(long, value_names = ["MIN", "MAX"])]
    pub scan_gq: Option<Vec<i32>>,

    /// Grid-scan integer DP thresholds (inclusive min/max) together with `--scan-gq`.
    /// A call "passes" a cell (gq_thr, dp_thr) only if GQ ≥ gq_thr **and** DP ≥ dp_thr.
    /// Emits `gq_dp_scan` TSV and uses 2D counting (cannot be combined with `--plot-gq-pr-png`).
    #[arg(long, value_names = ["MIN", "MAX"])]
    pub scan_dp: Option<Vec<i32>>,

    /// If set, write a recall-vs-precision **curve** for `--scan-gq` (single sample) or **overlaid translucent PR curves**
    /// (one per `--eval-sample`) in a 2×2 panel layout (total / SNV / ins / del) when multiple samples are used.
    /// GQ-only (`--scan-dp` not set).
    #[arg(long)]
    pub plot_gq_pr_png: Option<PathBuf>,

    /// If set, write an F1 heatmap (GQ × DP) for the **total** stratum to this PNG path.
    /// Requires `--scan-dp` (and `--scan-gq`).
    #[arg(long)]
    pub plot_gq_dp_f1_png: Option<PathBuf>,

    /// If set, write **recall vs precision**: **scatter** (one sample) or **overlaid translucent PR polylines**
    /// (one per sample; path follows GQ×DP row-major threshold order) in 2×2 panels (total / SNV / ins / del).
    /// Precision axis uses `--gq-dp-pr-precision-min/max`. Requires `--scan-dp` and `--scan-gq`.
    #[arg(long)]
    pub plot_gq_dp_pr_png: Option<PathBuf>,

    /// Lower bound of the **precision (y) axis** for `--plot-gq-dp-pr-png` (default zooms to 0.9–1.0).
    #[arg(long, default_value_t = 0.9)]
    pub gq_dp_pr_precision_min: f64,

    /// Upper bound of the **precision (y) axis** for `--plot-gq-dp-pr-png`.
    #[arg(long, default_value_t = 1.0)]
    pub gq_dp_pr_precision_max: f64,

    /// Emit allele presence/absence tables (TP/FP/FN). Off by default.
    #[arg(long, default_value_t = false)]
    pub emit_allele: bool,

    /// Optional path to write JSON summary (also prints TSV to stderr).
    #[arg(long)]
    pub json: Option<PathBuf>,
}

#[derive(Debug, Default, Clone, Serialize)]
struct AlleleCounts {
    tp: u64,
    fp: u64,
    fn_: u64,
    truth_nocall: u64,
    query_nocall: u64,
    compared: u64,
}

#[derive(Debug, Clone, Serialize)]
struct Summary {
    #[serde(skip_serializing_if = "Option::is_none")]
    allele: Option<StratifiedAlleleCounts>,
    genotype: StratifiedGenotypeCounts,
}

/// JSON output when `--eval-sample` lists more than one sample (or multiple pairs).
#[derive(Debug, Serialize)]
struct MultiEvalSummary {
    pooled: Summary,
    per_sample: BTreeMap<String, Summary>,
}

#[derive(Debug, Default, Clone, Serialize)]
struct GenotypeCounts {
    exact: u64,
    compared: u64,
    truth_nocall: u64,
    query_nocall: u64,
}

#[derive(Debug, Default, Clone, Serialize)]
struct StratifiedAlleleCounts {
    total: AlleleCounts,
    snv: AlleleCounts,
    insertion: AlleleCounts,
    deletion: AlleleCounts,
    insertion_by_len: BTreeMap<u64, AlleleCounts>,
    deletion_by_len: BTreeMap<u64, AlleleCounts>,
}

#[derive(Debug, Default, Clone, Serialize)]
struct StratifiedGenotypeCounts {
    total: GenotypeCounts,
    snv: GenotypeCounts,
    insertion: GenotypeCounts,
    deletion: GenotypeCounts,
    insertion_by_len: BTreeMap<u64, GenotypeCounts>,
    deletion_by_len: BTreeMap<u64, GenotypeCounts>,
}

#[derive(Debug, Clone, Eq, PartialEq)]
struct Key {
    chrom: String,
    pos1: i64,
    r: Vec<u8>,
    a: Vec<u8>,
}

impl Ord for Key {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.chrom
            .cmp(&other.chrom)
            .then_with(|| self.pos1.cmp(&other.pos1))
            .then_with(|| self.r.cmp(&other.r))
            .then_with(|| self.a.cmp(&other.a))
    }
}
impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone)]
struct Biallelic {
    key: Key,
    gt: Option<(u8, u8)>, // alt copies: 0,1,2 for this specific ALT; None = nocall
    vc: VarClass,
    indel_len: Option<u64>,
    gq: Option<i32>, // only populated for query; truth uses None
    /// Query FORMAT DP (or missing); only populated when `--scan-dp` is used.
    dp: Option<i32>,
}

pub fn run_eval(args: EvalArgs) -> Result<()> {
    let bed = BedMask::from_path(&args.confident_bed)
        .with_context(|| format!("Failed to read confident BED: {:?}", args.confident_bed))?;

    let region = args
        .region
        .as_deref()
        .map(parse_region)
        .transpose()
        .context("Failed to parse --region")?;

    let truth_reader =
        bcf::Reader::from_path(&args.truth).with_context(|| format!("Failed to open truth: {:?}", args.truth))?;
    let query_reader =
        bcf::Reader::from_path(&args.query).with_context(|| format!("Failed to open query: {:?}", args.query))?;

    let pairs = eval_name_pairs(&args, truth_reader.header(), query_reader.header())?;
    let multi = pairs.len() > 1;

    if args.scan_dp.is_some() && args.scan_gq.is_none() {
        bail!("--scan-dp requires --scan-gq");
    }
    if args.plot_gq_pr_png.is_some() && args.scan_dp.is_some() {
        bail!(
            "--plot-gq-pr-png is only for GQ-only scans; omit --scan-dp or use --plot-gq-dp-f1-png for GQ×DP heatmaps"
        );
    }
    if args.plot_gq_dp_f1_png.is_some() && args.scan_dp.is_none() {
        bail!("--plot-gq-dp-f1-png requires --scan-dp");
    }
    if args.plot_gq_dp_pr_png.is_some() && args.scan_dp.is_none() {
        bail!("--plot-gq-dp-pr-png requires --scan-dp and --scan-gq");
    }

    let gq_scan = args.scan_gq.as_ref().map(|v| parse_gq_scan(v)).transpose()?;
    let dp_scan = args.scan_dp.as_ref().map(|v| parse_dp_scan(v)).transpose()?;

    let scan_template: ScanTemplate = match (&gq_scan, &dp_scan) {
        (Some(gq), None) => ScanTemplate::Gq(gq.clone()),
        (Some(gq), Some(dp)) => ScanTemplate::GqDp(GqDpScan {
            gq: gq.thresholds.clone(),
            dp: dp.thresholds.clone(),
        }),
        (None, None) => ScanTemplate::None,
        (None, Some(_)) => unreachable!("--scan-dp without --scan-gq should have been rejected"),
    };

    let mut pooled_allele = StratifiedAlleleCounts::default();
    let mut pooled_genotype = StratifiedGenotypeCounts::default();
    let mut per_sample_summaries: BTreeMap<String, Summary> = BTreeMap::new();
    let mut gq_per_sample: Vec<StratifiedScanCounts> = Vec::new();
    let mut gqdp_per_sample: Vec<StratifiedScanCounts2D> = Vec::new();

    if multi && args.emit_allele {
        eprintln!("## allele_sample\tsample\tclass\ttp\tfp\tfn\tprecision\trecall\tf1\tcompared\ttruth_nocall\tquery_nocall");
        eprintln!("## allele_indel_len_sample\tsample\ttype\tlen\ttp\tfp\tfn\tprecision\trecall\tf1\tcompared\ttruth_nocall\tquery_nocall");
    }
    if multi {
        eprintln!("## genotype_sample\tsample\ttype\tlen\texact\tcompared\texact_rate\ttruth_nocall\tquery_nocall");
    }
    if multi && matches!(&scan_template, ScanTemplate::Gq(_)) {
        eprintln!("## gq_scan_sample\tsample\tgq_threshold\ttp\tfp\tfn\ttn\tprecision\trecall\tf1");
    }
    if multi && matches!(&scan_template, ScanTemplate::GqDp(_)) {
        eprintln!("## gq_dp_scan_sample\tsample\tgq_thr\tdp_thr\tclass\ttp\tfp\tfn\ttn\tprecision\trecall\tf1");
    }

    for (truth_name, query_name) in &pairs {
        let (allele, genotype, ts) = run_eval_for_pair(
            &args.truth,
            &args.query,
            truth_name,
            query_name,
            &bed,
            &region,
            &args,
            &scan_template,
        )?;

        merge_stratified_allele(&mut pooled_allele, &allele);
        merge_stratified_genotype(&mut pooled_genotype, &genotype);

        let sample_key = if truth_name == query_name {
            truth_name.clone()
        } else {
            format!("{truth_name}__{query_name}")
        };
        if args.emit_allele {
            if multi {
                emit_allele_table_prefixed("allele_sample", Some(&sample_key), &allele);
            } else {
                emit_allele_table(&allele);
            }
        }
        if multi {
            emit_genotype_table_sample(&sample_key, &genotype);
        } else {
            emit_genotype_table(&genotype);
        }

        if let Some(ts) = ts {
            match ts {
                ThresholdScan::Gq { scan, counts } => {
                    if multi {
                        emit_gq_scan_table_sample(&sample_key, &scan, &counts.total);
                    }
                    gq_per_sample.push(counts);
                }
                ThresholdScan::GqDp { scan, counts } => {
                    if multi {
                        emit_gq_dp_scan_table_sample(&sample_key, &scan, &counts);
                    }
                    gqdp_per_sample.push(counts);
                }
            }
        }

        per_sample_summaries.insert(
            sample_key,
            Summary {
                allele: if args.emit_allele { Some(allele) } else { None },
                genotype,
            },
        );
    }

    if multi {
        eprintln!("## genotype_pooled\ttype\tlen\texact\tcompared\texact_rate\ttruth_nocall\tquery_nocall");
        emit_genotype_table_body(&pooled_genotype);
        if args.emit_allele {
            emit_allele_table_pooled(&pooled_allele);
        }
    }

    match &scan_template {
        ScanTemplate::Gq(gq) if !gq_per_sample.is_empty() => {
            let mut pooled = StratifiedScanCounts::new(gq.thresholds.len());
            for c in &gq_per_sample {
                pooled.add_assign(c);
            }
            emit_gq_scan_table(gq, &pooled.total);
            if let Some(png) = &args.plot_gq_pr_png {
                if gq_per_sample.len() > 1 {
                    write_gq_pr_multi_sample_curves_png(gq, &gq_per_sample, png, &args)
                        .with_context(|| format!("Failed to write multi-sample GQ PR curves PNG: {png:?}"))?;
                } else {
                    write_gq_pr_png(gq, &gq_per_sample[0], png, &args)
                        .with_context(|| format!("Failed to write PR PNG: {png:?}"))?;
                }
            }
        }
        ScanTemplate::GqDp(scan_def) if !gqdp_per_sample.is_empty() => {
            let mut pooled = StratifiedScanCounts2D::new(scan_def.gq.len(), scan_def.dp.len());
            for c in &gqdp_per_sample {
                pooled.add_assign(c);
            }
            emit_gq_dp_scan_table(scan_def, &pooled);
            if let Some(png) = &args.plot_gq_dp_f1_png {
                write_gq_dp_f1_heatmap_png(scan_def, &pooled, png, &args)
                    .with_context(|| format!("Failed to write GQ×DP F1 heatmap PNG: {png:?}"))?;
            }
            if let Some(png) = &args.plot_gq_dp_pr_png {
                if gqdp_per_sample.len() > 1 {
                    write_gq_dp_pr_multi_sample_curves_png(scan_def, &gqdp_per_sample, png, &args)
                        .with_context(|| format!("Failed to write GQ×DP multi-sample PR curves PNG: {png:?}"))?;
                } else {
                    write_gq_dp_pr_scatter_png(scan_def, &gqdp_per_sample[0], png, &args)
                        .with_context(|| format!("Failed to write GQ×DP PR scatter PNG: {png:?}"))?;
                }
            }
        }
        _ => {}
    }

    if args.plot_gq_pr_png.is_some() && gq_scan.is_none() {
        bail!("--plot-gq-pr-png requires --scan-gq");
    }
    if args.plot_gq_dp_f1_png.is_some() && args.scan_dp.is_none() {
        bail!("--plot-gq-dp-f1-png requires --scan-dp and --scan-gq");
    }
    if args.plot_gq_dp_pr_png.is_some() && args.scan_dp.is_none() {
        bail!("--plot-gq-dp-pr-png requires --scan-dp and --scan-gq");
    }

    let pooled_summary = Summary {
        allele: if args.emit_allele {
            Some(pooled_allele)
        } else {
            None
        },
        genotype: pooled_genotype,
    };

    if let Some(out) = &args.json {
        let f = File::create(out).with_context(|| format!("Failed to create json: {out:?}"))?;
        if multi {
            serde_json::to_writer_pretty(
                f,
                &MultiEvalSummary {
                    pooled: pooled_summary,
                    per_sample: per_sample_summaries,
                },
            )
            .context("Failed to write json")?;
        } else {
            serde_json::to_writer_pretty(f, &pooled_summary).context("Failed to write json")?;
        }
    }
    Ok(())
}

#[derive(Debug, Clone)]
struct GqScan {
    thresholds: Vec<i32>,
}

#[derive(Debug, Clone)]
struct GqDpScan {
    gq: Vec<i32>,
    dp: Vec<i32>,
}

/// Template to build a fresh per-sample threshold accumulator.
#[derive(Debug, Clone)]
enum ScanTemplate {
    None,
    Gq(GqScan),
    GqDp(GqDpScan),
}

fn parse_int_range(v: &[i32], flag: &str) -> Result<Vec<i32>> {
    if v.len() != 2 {
        bail!("{flag} expects exactly two integers: MIN MAX");
    }
    let min = v[0];
    let max = v[1];
    if max < min {
        bail!("{flag} MAX must be >= MIN");
    }
    Ok((min..=max).collect())
}

fn parse_gq_scan(v: &Vec<i32>) -> Result<GqScan> {
    Ok(GqScan {
        thresholds: parse_int_range(v, "--scan-gq")?,
    })
}

fn parse_dp_scan(v: &Vec<i32>) -> Result<GqScan> {
    Ok(GqScan {
        thresholds: parse_int_range(v, "--scan-dp")?,
    })
}

/// GQ-only grid scan vs combined GQ×DP grid scan.
enum ThresholdScan {
    Gq {
        scan: GqScan,
        counts: StratifiedScanCounts,
    },
    GqDp {
        scan: GqDpScan,
        counts: StratifiedScanCounts2D,
    },
}

fn monotonic_hi(thresholds: &[i32], value: i32) -> Option<usize> {
    let mut hi = None;
    for (i, thr) in thresholds.iter().enumerate() {
        if *thr <= value {
            hi = Some(i);
        } else {
            break;
        }
    }
    hi
}

#[derive(Debug, Default, Clone)]
struct ScanCounts {
    tp: u64,
    fp: u64,
    fn_: u64,
    tn: u64,
}

#[derive(Debug, Default, Clone)]
struct StratifiedScanCounts {
    total: Vec<ScanCounts>,
    snv: Vec<ScanCounts>,
    insertion: Vec<ScanCounts>,
    deletion: Vec<ScanCounts>,
}

impl StratifiedScanCounts {
    fn new(len: usize) -> Self {
        Self {
            total: vec![ScanCounts::default(); len],
            snv: vec![ScanCounts::default(); len],
            insertion: vec![ScanCounts::default(); len],
            deletion: vec![ScanCounts::default(); len],
        }
    }

    fn slice_mut(&mut self, vc: VarClass) -> &mut [ScanCounts] {
        match vc {
            VarClass::Snv => &mut self.snv,
            VarClass::Insertion => &mut self.insertion,
            VarClass::Deletion => &mut self.deletion,
        }
    }

    fn add_assign(&mut self, o: &Self) {
        assert_eq!(self.total.len(), o.total.len());
        for i in 0..self.total.len() {
            self.total[i].tp += o.total[i].tp;
            self.total[i].fp += o.total[i].fp;
            self.total[i].fn_ += o.total[i].fn_;
            self.total[i].tn += o.total[i].tn;
            self.snv[i].tp += o.snv[i].tp;
            self.snv[i].fp += o.snv[i].fp;
            self.snv[i].fn_ += o.snv[i].fn_;
            self.snv[i].tn += o.snv[i].tn;
            self.insertion[i].tp += o.insertion[i].tp;
            self.insertion[i].fp += o.insertion[i].fp;
            self.insertion[i].fn_ += o.insertion[i].fn_;
            self.insertion[i].tn += o.insertion[i].tn;
            self.deletion[i].tp += o.deletion[i].tp;
            self.deletion[i].fp += o.deletion[i].fp;
            self.deletion[i].fn_ += o.deletion[i].fn_;
            self.deletion[i].tn += o.deletion[i].tn;
        }
    }
}

#[derive(Debug, Clone)]
struct StratifiedScanCounts2D {
    n_dp: usize,
    total: Vec<ScanCounts>,
    snv: Vec<ScanCounts>,
    insertion: Vec<ScanCounts>,
    deletion: Vec<ScanCounts>,
}

impl StratifiedScanCounts2D {
    fn new(n_gq: usize, n_dp: usize) -> Self {
        let len = n_gq * n_dp;
        Self {
            n_dp,
            total: vec![ScanCounts::default(); len],
            snv: vec![ScanCounts::default(); len],
            insertion: vec![ScanCounts::default(); len],
            deletion: vec![ScanCounts::default(); len],
        }
    }

    fn idx(&self, i_gq: usize, i_dp: usize) -> usize {
        i_gq * self.n_dp + i_dp
    }

    fn slice_mut(&mut self, vc: VarClass) -> &mut [ScanCounts] {
        match vc {
            VarClass::Snv => &mut self.snv,
            VarClass::Insertion => &mut self.insertion,
            VarClass::Deletion => &mut self.deletion,
        }
    }

    fn add_assign(&mut self, o: &Self) {
        assert_eq!(self.total.len(), o.total.len());
        for i in 0..self.total.len() {
            self.total[i].tp += o.total[i].tp;
            self.total[i].fp += o.total[i].fp;
            self.total[i].fn_ += o.total[i].fn_;
            self.total[i].tn += o.total[i].tn;
            self.snv[i].tp += o.snv[i].tp;
            self.snv[i].fp += o.snv[i].fp;
            self.snv[i].fn_ += o.snv[i].fn_;
            self.snv[i].tn += o.snv[i].tn;
            self.insertion[i].tp += o.insertion[i].tp;
            self.insertion[i].fp += o.insertion[i].fp;
            self.insertion[i].fn_ += o.insertion[i].fn_;
            self.insertion[i].tn += o.insertion[i].tn;
            self.deletion[i].tp += o.deletion[i].tp;
            self.deletion[i].fp += o.deletion[i].fp;
            self.deletion[i].fn_ += o.deletion[i].fn_;
            self.deletion[i].tn += o.deletion[i].tn;
        }
    }
}

fn merge_stratified_allele(dst: &mut StratifiedAlleleCounts, src: &StratifiedAlleleCounts) {
    merge_allele(&mut dst.total, &src.total);
    merge_allele(&mut dst.snv, &src.snv);
    merge_allele(&mut dst.insertion, &src.insertion);
    merge_allele(&mut dst.deletion, &src.deletion);
    for (k, v) in &src.insertion_by_len {
        merge_allele(dst.insertion_by_len.entry(*k).or_default(), v);
    }
    for (k, v) in &src.deletion_by_len {
        merge_allele(dst.deletion_by_len.entry(*k).or_default(), v);
    }
}

fn merge_allele(dst: &mut AlleleCounts, src: &AlleleCounts) {
    dst.tp += src.tp;
    dst.fp += src.fp;
    dst.fn_ += src.fn_;
    dst.truth_nocall += src.truth_nocall;
    dst.query_nocall += src.query_nocall;
    dst.compared += src.compared;
}

fn merge_stratified_genotype(dst: &mut StratifiedGenotypeCounts, src: &StratifiedGenotypeCounts) {
    merge_genotype(&mut dst.total, &src.total);
    merge_genotype(&mut dst.snv, &src.snv);
    merge_genotype(&mut dst.insertion, &src.insertion);
    merge_genotype(&mut dst.deletion, &src.deletion);
    for (k, v) in &src.insertion_by_len {
        merge_genotype(dst.insertion_by_len.entry(*k).or_default(), v);
    }
    for (k, v) in &src.deletion_by_len {
        merge_genotype(dst.deletion_by_len.entry(*k).or_default(), v);
    }
}

fn merge_genotype(dst: &mut GenotypeCounts, src: &GenotypeCounts) {
    dst.exact += src.exact;
    dst.compared += src.compared;
    dst.truth_nocall += src.truth_nocall;
    dst.query_nocall += src.query_nocall;
}

fn eval_name_pairs(
    args: &EvalArgs,
    truth_h: &rust_htslib::bcf::header::HeaderView,
    query_h: &rust_htslib::bcf::header::HeaderView,
) -> Result<Vec<(String, String)>> {
    if !args.eval_sample.is_empty() {
        if args.truth_sample.is_some() || args.query_sample.is_some() {
            bail!("Do not combine --eval-sample with --truth-sample / --query-sample");
        }
        let mut seen = HashSet::new();
        let mut pairs_out: Vec<(String, String)> = Vec::new();
        for s in &args.eval_sample {
            if !seen.insert(s.clone()) {
                continue;
            }
            sample_index_required(truth_h, s)?;
            sample_index_required(query_h, s)?;
            pairs_out.push((s.clone(), s.clone()));
        }
        if pairs_out.is_empty() {
            bail!("--eval-sample listed no unique samples");
        }
        return Ok(pairs_out);
    }
    let t = truth_sample_name(truth_h, args.truth_sample.as_deref())?;
    let q = query_sample_name(query_h, args.query_sample.as_deref())?;
    if args.require_sample_name_match && t != q {
        bail!(
            "Sample name mismatch (use --require-sample-name-match=false to override): truth={:?} query={:?}",
            t,
            q
        );
    }
    Ok(vec![(t, q)])
}

fn threshold_scan_from_template(tpl: &ScanTemplate) -> Option<ThresholdScan> {
    match tpl {
        ScanTemplate::None => None,
        ScanTemplate::Gq(gq) => Some(ThresholdScan::Gq {
            scan: gq.clone(),
            counts: StratifiedScanCounts::new(gq.thresholds.len()),
        }),
        ScanTemplate::GqDp(s) => Some(ThresholdScan::GqDp {
            scan: GqDpScan {
                gq: s.gq.clone(),
                dp: s.dp.clone(),
            },
            counts: StratifiedScanCounts2D::new(s.gq.len(), s.dp.len()),
        }),
    }
}

fn run_eval_for_pair(
    truth_path: &PathBuf,
    query_path: &PathBuf,
    truth_name: &str,
    query_name: &str,
    bed: &BedMask,
    region: &Option<Region>,
    args: &EvalArgs,
    scan_template: &ScanTemplate,
) -> Result<(StratifiedAlleleCounts, StratifiedGenotypeCounts, Option<ThresholdScan>)> {
    let mut truth_reader =
        bcf::Reader::from_path(truth_path).with_context(|| format!("Failed to open truth: {truth_path:?}"))?;
    let mut query_reader =
        bcf::Reader::from_path(query_path).with_context(|| format!("Failed to open query: {query_path:?}"))?;

    let truth_sample_i = sample_index_required(truth_reader.header(), truth_name)?;
    let query_sample_i = sample_index_required(query_reader.header(), query_name)?;

    let mut threshold_scan = threshold_scan_from_template(scan_template);

    let need_gq = threshold_scan.is_some();
    let need_dp = matches!(&threshold_scan, Some(ThresholdScan::GqDp { .. }));

    let mut truth_it = BiallelicIter::new(
        &mut truth_reader,
        truth_sample_i,
        bed.clone(),
        region.clone(),
        None,
        None,
    )?;
    let mut query_it = BiallelicIter::new(
        &mut query_reader,
        query_sample_i,
        bed.clone(),
        region.clone(),
        need_gq.then(|| args.gq_tag.as_bytes().to_vec()),
        need_dp.then(|| args.dp_tag.as_bytes().to_vec()),
    )?;

    let mut allele = StratifiedAlleleCounts::default();
    let mut genotype = StratifiedGenotypeCounts::default();

    let mut t = truth_it.next()?;
    let mut q = query_it.next()?;

    while t.is_some() || q.is_some() {
        match (&t, &q) {
            (Some(tv), Some(qv)) => match tv.key.cmp(&qv.key) {
                std::cmp::Ordering::Less => {
                    score_side_only(&mut allele, &mut genotype, tv, true, args);
                    if let Some(ts) = &mut threshold_scan {
                        match ts {
                            ThresholdScan::Gq { scan, counts } => {
                                update_gq_scan_truth_only(counts, scan, tv, args);
                            }
                            ThresholdScan::GqDp { scan, counts } => {
                                update_gq_dp_scan_truth_only(counts, scan, tv, args);
                            }
                        }
                    }
                    t = truth_it.next()?;
                }
                std::cmp::Ordering::Greater => {
                    score_side_only(&mut allele, &mut genotype, qv, false, args);
                    if let Some(ts) = &mut threshold_scan {
                        match ts {
                            ThresholdScan::Gq { scan, counts } => {
                                update_gq_scan_query_only(counts, scan, qv, args);
                            }
                            ThresholdScan::GqDp { scan, counts } => {
                                update_gq_dp_scan_query_only(counts, scan, qv, args);
                            }
                        }
                    }
                    q = query_it.next()?;
                }
                std::cmp::Ordering::Equal => {
                    score_both(&mut allele, &mut genotype, tv, qv, args);
                    if let Some(ts) = &mut threshold_scan {
                        match ts {
                            ThresholdScan::Gq { scan, counts } => {
                                update_gq_scan(counts, scan, tv, qv, args);
                            }
                            ThresholdScan::GqDp { scan, counts } => {
                                update_gq_dp_scan(counts, scan, tv, qv, args);
                            }
                        }
                    }
                    t = truth_it.next()?;
                    q = query_it.next()?;
                }
            },
            (Some(tv), None) => {
                score_side_only(&mut allele, &mut genotype, tv, true, args);
                if let Some(ts) = &mut threshold_scan {
                    match ts {
                        ThresholdScan::Gq { scan, counts } => {
                            update_gq_scan_truth_only(counts, scan, tv, args);
                        }
                        ThresholdScan::GqDp { scan, counts } => {
                            update_gq_dp_scan_truth_only(counts, scan, tv, args);
                        }
                    }
                }
                t = truth_it.next()?;
            }
            (None, Some(qv)) => {
                score_side_only(&mut allele, &mut genotype, qv, false, args);
                if let Some(ts) = &mut threshold_scan {
                    match ts {
                        ThresholdScan::Gq { scan, counts } => {
                            update_gq_scan_query_only(counts, scan, qv, args);
                        }
                        ThresholdScan::GqDp { scan, counts } => {
                            update_gq_dp_scan_query_only(counts, scan, qv, args);
                        }
                    }
                }
                q = query_it.next()?;
            }
            (None, None) => break,
        }
    }

    Ok((allele, genotype, threshold_scan))
}

fn update_gq_scan(
    out: &mut StratifiedScanCounts,
    scan: &GqScan,
    t: &Biallelic,
    q: &Biallelic,
    args: &EvalArgs,
) {
    // For thresholding we treat "filtered by GQ" as "not called" for the query side.
    let t_present = t.gt.map(|(a, b)| (a as u64 + b as u64) > 0).unwrap_or(false);
    let q_present_base = q.gt.map(|(a, b)| (a as u64 + b as u64) > 0).unwrap_or(false);

    if args.require_sample_nonref && !t_present && !q_present_base {
        return;
    }

    // Missing GQ is treated as -1 (fails any threshold >=0).
    let gq = q.gq.unwrap_or(-1);

    // Find largest index i such that thresholds[i] <= gq.
    let mut hi = None;
    for (i, thr) in scan.thresholds.iter().enumerate() {
        if *thr <= gq {
            hi = Some(i);
        } else {
            break;
        }
    }

    match (t_present, q_present_base) {
        (true, true) => {
            // TP for thresholds <= gq, FN for thresholds > gq.
            if let Some(i) = hi {
                for c in out.total.iter_mut().take(i + 1) {
                    c.tp += 1;
                }
                for c in out.total.iter_mut().skip(i + 1) {
                    c.fn_ += 1;
                }
                let by = out.slice_mut(t.vc);
                for c in by.iter_mut().take(i + 1) {
                    c.tp += 1;
                }
                for c in by.iter_mut().skip(i + 1) {
                    c.fn_ += 1;
                }
            } else {
                // gq below min => FN for all thresholds
                for c in out.total.iter_mut() {
                    c.fn_ += 1;
                }
                let by = out.slice_mut(t.vc);
                for c in by.iter_mut() {
                    c.fn_ += 1;
                }
            }
        }
        (true, false) => {
            // FN for all thresholds
            for c in out.total.iter_mut() {
                c.fn_ += 1;
            }
            let by = out.slice_mut(t.vc);
            for c in by.iter_mut() {
                c.fn_ += 1;
            }
        }
        (false, true) => {
            // FP for thresholds <= gq, otherwise TN.
            if let Some(i) = hi {
                for c in out.total.iter_mut().take(i + 1) {
                    c.fp += 1;
                }
                for c in out.total.iter_mut().skip(i + 1) {
                    c.tn += 1;
                }
                let by = out.slice_mut(t.vc);
                for c in by.iter_mut().take(i + 1) {
                    c.fp += 1;
                }
                for c in by.iter_mut().skip(i + 1) {
                    c.tn += 1;
                }
            } else {
                // gq below min => TN for all thresholds
                for c in out.total.iter_mut() {
                    c.tn += 1;
                }
                let by = out.slice_mut(t.vc);
                for c in by.iter_mut() {
                    c.tn += 1;
                }
            }
        }
        (false, false) => {
            // TN
            for c in out.total.iter_mut() {
                c.tn += 1;
            }
            let by = out.slice_mut(t.vc);
            for c in by.iter_mut() {
                c.tn += 1;
            }
        }
    }
}

fn update_gq_scan_truth_only(out: &mut StratifiedScanCounts, _scan: &GqScan, t: &Biallelic, args: &EvalArgs) {
    let t_present = t.gt.map(|(a, b)| (a as u64 + b as u64) > 0).unwrap_or(false);
    if args.require_sample_nonref && !t_present {
        return;
    }
    if t_present {
        for c in out.total.iter_mut() {
            c.fn_ += 1;
        }
        let by = out.slice_mut(t.vc);
        for c in by.iter_mut() {
            c.fn_ += 1;
        }
    } else {
        for c in out.total.iter_mut() {
            c.tn += 1;
        }
        let by = out.slice_mut(t.vc);
        for c in by.iter_mut() {
            c.tn += 1;
        }
    }
}

fn update_gq_scan_query_only(out: &mut StratifiedScanCounts, scan: &GqScan, q: &Biallelic, args: &EvalArgs) {
    let q_present_base = q.gt.map(|(a, b)| (a as u64 + b as u64) > 0).unwrap_or(false);
    if args.require_sample_nonref && !q_present_base {
        return;
    }
    let gq = q.gq.unwrap_or(-1);
    let mut hi = None;
    for (i, thr) in scan.thresholds.iter().enumerate() {
        if *thr <= gq {
            hi = Some(i);
        } else {
            break;
        }
    }
    if q_present_base {
        // FP for thresholds <= gq, else TN
        if let Some(i) = hi {
            for c in out.total.iter_mut().take(i + 1) {
                c.fp += 1;
            }
            for c in out.total.iter_mut().skip(i + 1) {
                c.tn += 1;
            }
            let by = out.slice_mut(q.vc);
            for c in by.iter_mut().take(i + 1) {
                c.fp += 1;
            }
            for c in by.iter_mut().skip(i + 1) {
                c.tn += 1;
            }
        } else {
            for c in out.total.iter_mut() {
                c.tn += 1;
            }
            let by = out.slice_mut(q.vc);
            for c in by.iter_mut() {
                c.tn += 1;
            }
        }
    } else {
        for c in out.total.iter_mut() {
            c.tn += 1;
        }
        let by = out.slice_mut(q.vc);
        for c in by.iter_mut() {
            c.tn += 1;
        }
    }
}

fn update_gq_dp_scan(
    out: &mut StratifiedScanCounts2D,
    scan: &GqDpScan,
    t: &Biallelic,
    q: &Biallelic,
    args: &EvalArgs,
) {
    let t_present = t.gt.map(|(a, b)| (a as u64 + b as u64) > 0).unwrap_or(false);
    let q_present_base = q.gt.map(|(a, b)| (a as u64 + b as u64) > 0).unwrap_or(false);

    if args.require_sample_nonref && !t_present && !q_present_base {
        return;
    }

    let gq = q.gq.unwrap_or(-1);
    let dp = q.dp.unwrap_or(-1);
    let i_max = monotonic_hi(&scan.gq, gq);
    let j_max = monotonic_hi(&scan.dp, dp);

    let pass_ij = |i: usize, j: usize| -> bool {
        i_max.map_or(false, |im| i <= im) && j_max.map_or(false, |jm| j <= jm)
    };

    match (t_present, q_present_base) {
        (true, true) => {
            for i in 0..scan.gq.len() {
                for j in 0..scan.dp.len() {
                    let idx = out.idx(i, j);
                    if pass_ij(i, j) {
                        out.total[idx].tp += 1;
                        out.slice_mut(t.vc)[idx].tp += 1;
                    } else {
                        out.total[idx].fn_ += 1;
                        out.slice_mut(t.vc)[idx].fn_ += 1;
                    }
                }
            }
        }
        (true, false) => {
            for i in 0..scan.gq.len() {
                for j in 0..scan.dp.len() {
                    let idx = out.idx(i, j);
                    out.total[idx].fn_ += 1;
                    out.slice_mut(t.vc)[idx].fn_ += 1;
                }
            }
        }
        (false, true) => {
            for i in 0..scan.gq.len() {
                for j in 0..scan.dp.len() {
                    let idx = out.idx(i, j);
                    if pass_ij(i, j) {
                        out.total[idx].fp += 1;
                        out.slice_mut(t.vc)[idx].fp += 1;
                    } else {
                        out.total[idx].tn += 1;
                        out.slice_mut(t.vc)[idx].tn += 1;
                    }
                }
            }
        }
        (false, false) => {
            for i in 0..scan.gq.len() {
                for j in 0..scan.dp.len() {
                    let idx = out.idx(i, j);
                    out.total[idx].tn += 1;
                    out.slice_mut(t.vc)[idx].tn += 1;
                }
            }
        }
    }
}

fn update_gq_dp_scan_truth_only(out: &mut StratifiedScanCounts2D, _scan: &GqDpScan, t: &Biallelic, args: &EvalArgs) {
    let t_present = t.gt.map(|(a, b)| (a as u64 + b as u64) > 0).unwrap_or(false);
    if args.require_sample_nonref && !t_present {
        return;
    }
    if t_present {
        for c in out.total.iter_mut() {
            c.fn_ += 1;
        }
        let by = out.slice_mut(t.vc);
        for c in by.iter_mut() {
            c.fn_ += 1;
        }
    } else {
        for c in out.total.iter_mut() {
            c.tn += 1;
        }
        let by = out.slice_mut(t.vc);
        for c in by.iter_mut() {
            c.tn += 1;
        }
    }
}

fn update_gq_dp_scan_query_only(out: &mut StratifiedScanCounts2D, scan: &GqDpScan, q: &Biallelic, args: &EvalArgs) {
    let q_present_base = q.gt.map(|(a, b)| (a as u64 + b as u64) > 0).unwrap_or(false);
    if args.require_sample_nonref && !q_present_base {
        return;
    }
    let gq = q.gq.unwrap_or(-1);
    let dp = q.dp.unwrap_or(-1);
    let i_max = monotonic_hi(&scan.gq, gq);
    let j_max = monotonic_hi(&scan.dp, dp);
    let pass_ij = |i: usize, j: usize| -> bool {
        i_max.map_or(false, |im| i <= im) && j_max.map_or(false, |jm| j <= jm)
    };

    if q_present_base {
        for i in 0..scan.gq.len() {
            for j in 0..scan.dp.len() {
                let idx = out.idx(i, j);
                if pass_ij(i, j) {
                    out.total[idx].fp += 1;
                    out.slice_mut(q.vc)[idx].fp += 1;
                } else {
                    out.total[idx].tn += 1;
                    out.slice_mut(q.vc)[idx].tn += 1;
                }
            }
        }
    } else {
        for c in out.total.iter_mut() {
            c.tn += 1;
        }
        let by = out.slice_mut(q.vc);
        for c in by.iter_mut() {
            c.tn += 1;
        }
    }
}

fn emit_gq_scan_table(scan: &GqScan, counts: &[ScanCounts]) {
    eprintln!("## gq_scan\tgq_threshold\ttp\tfp\tfn\ttn\tprecision\trecall\tf1");
    for (thr, c) in scan.thresholds.iter().zip(counts.iter()) {
        let p = precision(c.tp, c.fp);
        let r = recall(c.tp, c.fn_);
        let f = f1(p, r);
        eprintln!(
            "gq_scan\t{thr}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}",
            c.tp, c.fp, c.fn_, c.tn, p, r, f
        );
    }
}

fn emit_gq_dp_scan_table(scan: &GqDpScan, strat: &StratifiedScanCounts2D) {
    eprintln!("## gq_dp_scan\tgq_thr\tdp_thr\tclass\ttp\tfp\tfn\ttn\tprecision\trecall\tf1");
    for (i, gq_thr) in scan.gq.iter().enumerate() {
        for (j, dp_thr) in scan.dp.iter().enumerate() {
            let idx = strat.idx(i, j);
            let rows = [
                ("total", &strat.total[idx]),
                ("snv", &strat.snv[idx]),
                ("insertion", &strat.insertion[idx]),
                ("deletion", &strat.deletion[idx]),
            ];
            for (class, c) in rows {
                let p = precision(c.tp, c.fp);
                let r = recall(c.tp, c.fn_);
                let f = f1(p, r);
                eprintln!(
                    "gq_dp_scan\t{gq_thr}\t{dp_thr}\t{class}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}",
                    c.tp, c.fp, c.fn_, c.tn, p, r, f
                );
            }
        }
    }
}

fn emit_allele_table_prefixed(prefix: &str, sample: Option<&str>, s: &StratifiedAlleleCounts) {
    let sm = sample.expect("emit_allele_table_prefixed requires a sample");
    emit_allele_row_sample(prefix, sm, "total", &s.total);
    emit_allele_row_sample(prefix, sm, "snv", &s.snv);
    emit_allele_row_sample(prefix, sm, "insertion", &s.insertion);
    emit_allele_row_sample(prefix, sm, "deletion", &s.deletion);
    if !s.insertion_by_len.is_empty() || !s.deletion_by_len.is_empty() {
        for (l, c) in &s.insertion_by_len {
            emit_allele_len_row_sample(sm, "insertion", *l, c);
        }
        for (l, c) in &s.deletion_by_len {
            emit_allele_len_row_sample(sm, "deletion", *l, c);
        }
    }
}

fn emit_allele_row_sample(prefix: &str, sample: &str, class: &str, c: &AlleleCounts) {
    let prec = precision(c.tp, c.fp);
    let rec = recall(c.tp, c.fn_);
    let f = f1(prec, rec);
    eprintln!(
        "{prefix}\t{sample}\t{class}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}",
        c.tp, c.fp, c.fn_, prec, rec, f, c.compared, c.truth_nocall, c.query_nocall
    );
}

fn emit_allele_len_row_sample(sample: &str, t: &str, len: u64, c: &AlleleCounts) {
    let prec = precision(c.tp, c.fp);
    let rec = recall(c.tp, c.fn_);
    let f = f1(prec, rec);
    eprintln!(
        "allele_indel_len_sample\t{sample}\t{t}\t{len}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}",
        c.tp, c.fp, c.fn_, prec, rec, f, c.compared, c.truth_nocall, c.query_nocall
    );
}

fn emit_allele_table_pooled(s: &StratifiedAlleleCounts) {
    eprintln!("## allele_pooled\tclass\ttp\tfp\tfn\tprecision\trecall\tf1\tcompared\ttruth_nocall\tquery_nocall");
    emit_allele_row_pooled("total", &s.total);
    emit_allele_row_pooled("snv", &s.snv);
    emit_allele_row_pooled("insertion", &s.insertion);
    emit_allele_row_pooled("deletion", &s.deletion);
    if !s.insertion_by_len.is_empty() || !s.deletion_by_len.is_empty() {
        eprintln!("## allele_indel_len_pooled\ttype\tlen\ttp\tfp\tfn\tprecision\trecall\tf1\tcompared\ttruth_nocall\tquery_nocall");
        for (l, c) in &s.insertion_by_len {
            emit_allele_len_row_pooled("insertion", *l, c);
        }
        for (l, c) in &s.deletion_by_len {
            emit_allele_len_row_pooled("deletion", *l, c);
        }
    }
}

fn emit_allele_row_pooled(class: &str, c: &AlleleCounts) {
    let prec = precision(c.tp, c.fp);
    let rec = recall(c.tp, c.fn_);
    let f = f1(prec, rec);
    eprintln!(
        "allele_pooled\t{class}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}",
        c.tp, c.fp, c.fn_, prec, rec, f, c.compared, c.truth_nocall, c.query_nocall
    );
}

fn emit_allele_len_row_pooled(t: &str, len: u64, c: &AlleleCounts) {
    let prec = precision(c.tp, c.fp);
    let rec = recall(c.tp, c.fn_);
    let f = f1(prec, rec);
    eprintln!(
        "allele_indel_len_pooled\t{t}\t{len}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}",
        c.tp, c.fp, c.fn_, prec, rec, f, c.compared, c.truth_nocall, c.query_nocall
    );
}

fn emit_genotype_table_sample(sample: &str, s: &StratifiedGenotypeCounts) {
    emit_genotype_row_sample(sample, "snv", 0, &s.snv);
    for (l, c) in &s.insertion_by_len {
        emit_genotype_row_sample(sample, "insertion", *l, c);
    }
    for (l, c) in &s.deletion_by_len {
        emit_genotype_row_sample(sample, "deletion", *l, c);
    }
}

fn emit_genotype_row_sample(sample: &str, t: &str, len: u64, c: &GenotypeCounts) {
    let rate = if c.compared == 0 {
        0.0
    } else {
        c.exact as f64 / c.compared as f64
    };
    eprintln!(
        "genotype_sample\t{sample}\t{t}\t{len}\t{}\t{}\t{:.6}\t{}\t{}",
        c.exact, c.compared, rate, c.truth_nocall, c.query_nocall
    );
}

fn emit_genotype_table_body(s: &StratifiedGenotypeCounts) {
    emit_genotype_row_pooled("snv", 0, &s.snv);
    for (l, c) in &s.insertion_by_len {
        emit_genotype_row_pooled("insertion", *l, c);
    }
    for (l, c) in &s.deletion_by_len {
        emit_genotype_row_pooled("deletion", *l, c);
    }
}

fn emit_genotype_row_pooled(t: &str, len: u64, c: &GenotypeCounts) {
    let rate = if c.compared == 0 {
        0.0
    } else {
        c.exact as f64 / c.compared as f64
    };
    eprintln!(
        "genotype_pooled\t{t}\t{len}\t{}\t{}\t{:.6}\t{}\t{}",
        c.exact, c.compared, rate, c.truth_nocall, c.query_nocall
    );
}

fn emit_gq_scan_table_sample(sample: &str, scan: &GqScan, counts: &[ScanCounts]) {
    for (thr, c) in scan.thresholds.iter().zip(counts.iter()) {
        let p = precision(c.tp, c.fp);
        let r = recall(c.tp, c.fn_);
        let f = f1(p, r);
        eprintln!(
            "gq_scan_sample\t{sample}\t{thr}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}",
            c.tp, c.fp, c.fn_, c.tn, p, r, f
        );
    }
}

fn emit_gq_dp_scan_table_sample(sample: &str, scan: &GqDpScan, strat: &StratifiedScanCounts2D) {
    for (i, gq_thr) in scan.gq.iter().enumerate() {
        for (j, dp_thr) in scan.dp.iter().enumerate() {
            let idx = strat.idx(i, j);
            let rows = [
                ("total", &strat.total[idx]),
                ("snv", &strat.snv[idx]),
                ("insertion", &strat.insertion[idx]),
                ("deletion", &strat.deletion[idx]),
            ];
            for (class, c) in rows {
                let p = precision(c.tp, c.fp);
                let r = recall(c.tp, c.fn_);
                let f = f1(p, r);
                eprintln!(
                    "gq_dp_scan_sample\t{sample}\t{gq_thr}\t{dp_thr}\t{class}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}",
                    c.tp, c.fp, c.fn_, c.tn, p, r, f
                );
            }
        }
    }
}

fn scan_at_2d(strat: &StratifiedScanCounts2D, vc: VarClass, idx: usize) -> &ScanCounts {
    match vc {
        VarClass::Snv => &strat.snv[idx],
        VarClass::Insertion => &strat.insertion[idx],
        VarClass::Deletion => &strat.deletion[idx],
    }
}

fn scan_at_1d(strat: &StratifiedScanCounts, vc: VarClass, idx: usize) -> &ScanCounts {
    match vc {
        VarClass::Snv => &strat.snv[idx],
        VarClass::Insertion => &strat.insertion[idx],
        VarClass::Deletion => &strat.deletion[idx],
    }
}

/// Tab10-like colors; `.mix(alpha)` keeps lines readable when many samples overlap.
fn sample_line_color(si: usize) -> RGBColor {
    const PALETTE: [RGBColor; 10] = [
        RGBColor(31, 119, 180),
        RGBColor(255, 127, 14),
        RGBColor(44, 160, 44),
        RGBColor(214, 39, 40),
        RGBColor(148, 103, 189),
        RGBColor(140, 86, 75),
        RGBColor(227, 119, 194),
        RGBColor(127, 127, 127),
        RGBColor(188, 189, 34),
        RGBColor(23, 190, 207),
    ];
    PALETTE[si % PALETTE.len()]
}

fn pr_polyline_1d(strat: &StratifiedScanCounts, pick: usize) -> Vec<(f64, f64)> {
    // pick: 0 total, 1 snv, 2 insertion, 3 deletion
    let sl = match pick {
        0 => &strat.total[..],
        1 => &strat.snv[..],
        2 => &strat.insertion[..],
        3 => &strat.deletion[..],
        _ => &strat.total[..],
    };
    sl.iter().map(pr_point).collect()
}

fn pr_polyline_2d_rowmajor(strat: &StratifiedScanCounts2D, scan: &GqDpScan, pick: usize) -> Vec<(f64, f64)> {
    let mut out = Vec::with_capacity(scan.gq.len() * scan.dp.len());
    for i in 0..scan.gq.len() {
        for j in 0..scan.dp.len() {
            let idx = strat.idx(i, j);
            let c = match pick {
                0 => &strat.total[idx],
                1 => &strat.snv[idx],
                2 => &strat.insertion[idx],
                3 => &strat.deletion[idx],
                _ => &strat.total[idx],
            };
            out.push(pr_point(c));
        }
    }
    out
}

fn write_gq_pr_multi_sample_curves_png(
    _scan: &GqScan,
    per_sample: &[StratifiedScanCounts],
    out_path: &PathBuf,
    args: &EvalArgs,
) -> Result<()> {
    if per_sample.is_empty() {
        bail!("No samples for PR curves");
    }
    let n = per_sample[0].total.len();
    for s in per_sample {
        if s.total.len() != n {
            bail!("Internal error: GQ scan length mismatch across samples");
        }
    }

    let panel_titles = ["total", "snv", "insertion", "deletion"];
    let sub_w = 520u32;
    let sub_h = 420u32;
    let root = BitMapBackend::new(out_path, (sub_w * 2 + 60, sub_h * 2 + 90)).into_drawing_area();
    root.fill(&WHITE)?;

    let areas = root.split_evenly((2, 2));
    let alpha_line = 0.42f64;
    let header = pr_plot_caption(args);

    for (pi, title) in panel_titles.iter().enumerate() {
        let panel = &areas[pi];
        let pick = pi;

        let series_per: Vec<Vec<(f64, f64)>> =
            per_sample.iter().map(|st| pr_polyline_1d(st, pick)).collect();
        let (xmin, xmax) = bounds01(series_per.iter().flat_map(|v| v.iter().map(|(r, _)| *r)));
        let (ymin, ymax) = bounds01(series_per.iter().flat_map(|v| v.iter().map(|(_, p)| *p)));

        let cap = if pi == 0 {
            format!("{header}\n{title} ({} samples)", per_sample.len())
        } else {
            format!("{title} ({} samples)", per_sample.len())
        };

        let mut chart = ChartBuilder::on(panel)
            .margin(10)
            .caption(cap, ("sans-serif", 10))
            .x_label_area_size(42)
            .y_label_area_size(46)
            .build_cartesian_2d(xmin..xmax, ymin..ymax)?;

        chart
            .configure_mesh()
            .x_desc("Recall = TP / (TP + FN)")
            .y_desc("Precision = TP / (TP + FP)")
            .draw()?;

        for (si, pts) in series_per.iter().enumerate() {
            let style = sample_line_color(si).mix(alpha_line).stroke_width(2);
            chart.draw_series(LineSeries::new(pts.clone(), style))?;
        }
    }

    root.present()?;
    Ok(())
}

fn write_gq_dp_pr_multi_sample_curves_png(
    scan: &GqDpScan,
    per_sample: &[StratifiedScanCounts2D],
    out_path: &PathBuf,
    args: &EvalArgs,
) -> Result<()> {
    let n_cells = scan.gq.len() * scan.dp.len();
    if n_cells == 0 || per_sample.is_empty() {
        bail!("No data for GQ×DP PR curves");
    }
    for s in per_sample {
        if s.total.len() != n_cells {
            bail!("Internal error: gq×dp sample length mismatch");
        }
    }

    let pmin = args.gq_dp_pr_precision_min;
    let pmax = args.gq_dp_pr_precision_max;
    if !(pmin.is_finite() && pmax.is_finite() && pmin < pmax && pmin >= 0.0 && pmax <= 1.0) {
        bail!("--gq-dp-pr-precision-min / --gq-dp-pr-precision-max must satisfy 0 <= min < max <= 1");
    }

    let panel_titles = ["total", "snv", "insertion", "deletion"];
    let sub_w = 520u32;
    let sub_h = 420u32;
    let root = BitMapBackend::new(out_path, (sub_w * 2 + 60, sub_h * 2 + 90)).into_drawing_area();
    root.fill(&WHITE)?;

    let areas = root.split_evenly((2, 2));
    let alpha_line = 0.42f64;
    let header = pr_plot_caption(args);

    for (pi, title) in panel_titles.iter().enumerate() {
        let panel = &areas[pi];
        let pick = pi;

        let series_per: Vec<Vec<(f64, f64)>> = per_sample
            .iter()
            .map(|st| pr_polyline_2d_rowmajor(st, scan, pick))
            .collect();
        let (xmin, xmax) = bounds01(series_per.iter().flat_map(|v| v.iter().map(|(r, _)| *r)));
        let ymin = pmin;
        let ymax = pmax;

        let cap = if pi == 0 {
            format!("{header}\n{title} ({} samples)", per_sample.len())
        } else {
            format!("{title} ({} samples)", per_sample.len())
        };

        let mut chart = ChartBuilder::on(panel)
            .margin(10)
            .caption(cap, ("sans-serif", 10))
            .x_label_area_size(42)
            .y_label_area_size(46)
            .build_cartesian_2d(xmin..xmax, ymin..ymax)?;

        chart
            .configure_mesh()
            .x_desc("Recall = TP / (TP + FN)")
            .y_desc("Precision (zoomed)")
            .draw()?;

        for (si, pts) in series_per.iter().enumerate() {
            let style = sample_line_color(si).mix(alpha_line).stroke_width(2);
            chart.draw_series(LineSeries::new(pts.clone(), style))?;
        }
    }

    root.present()?;
    Ok(())
}

fn f1_to_rgb(f1v: f64) -> RGBColor {
    let t = f1v.clamp(0.0, 1.0);
    let r = ((1.0 - t) * 220.0) as u8;
    let g = (t * 200.0) as u8;
    let b = 40u8;
    RGBColor(r, g, b)
}

fn write_gq_dp_f1_heatmap_png(
    scan: &GqDpScan,
    strat: &StratifiedScanCounts2D,
    out_path: &PathBuf,
    args: &EvalArgs,
) -> Result<()> {
    let n_gq = scan.gq.len();
    let n_dp = scan.dp.len();
    if strat.total.len() != n_gq * n_dp {
        bail!("Internal error: gq×dp scan length mismatch");
    }

    let root = BitMapBackend::new(out_path, (900, 720)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut caption = pr_plot_caption(args);
    caption.push_str("\nF1 heatmap (total; GQ ≥ x, DP ≥ y)");
    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .caption(caption, ("sans-serif", 18))
        .x_label_area_size(55)
        .y_label_area_size(55)
        .build_cartesian_2d(0f64..n_gq as f64, 0f64..n_dp as f64)?;

    chart
        .configure_mesh()
        .x_desc("GQ threshold (inclusive)")
        .y_desc("DP threshold (inclusive)")
        .x_label_formatter(&|x| {
            let i = *x as usize;
            if i < n_gq {
                format!("{}", scan.gq[i])
            } else {
                String::new()
            }
        })
        .y_label_formatter(&|y| {
            let j = *y as usize;
            if j < n_dp {
                format!("{}", scan.dp[j])
            } else {
                String::new()
            }
        })
        .disable_mesh()
        .draw()?;

    for i in 0..n_gq {
        for j in 0..n_dp {
            let c = &strat.total[strat.idx(i, j)];
            let p = precision(c.tp, c.fp);
            let r = recall(c.tp, c.fn_);
            let f1v = f1(p, r);
            let color = f1_to_rgb(f1v);
            chart.draw_series(std::iter::once(Rectangle::new(
                [(i as f64, j as f64), ((i + 1) as f64, (j + 1) as f64)],
                color.filled(),
            )))?;
            if n_gq <= 25 && n_dp <= 25 {
                let label = format!("{:.2}", f1v);
                chart.draw_series(std::iter::once(Text::new(
                    label,
                    (i as f64 + 0.05, j as f64 + 0.12),
                    ("sans-serif", 8).into_font().color(&BLACK),
                )))?;
            }
        }
    }

    root.present()?;
    Ok(())
}

/// One point per `(GQ≥thr, DP≥thr)` cell: **x = recall**, **y = precision** (same convention as `--plot-gq-pr-png`).
/// Color encodes F1 so overlapping points are easier to read than a single hue.
fn write_gq_dp_pr_scatter_png(
    scan: &GqDpScan,
    strat: &StratifiedScanCounts2D,
    out_path: &PathBuf,
    args: &EvalArgs,
) -> Result<()> {
    let n_gq = scan.gq.len();
    let n_dp = scan.dp.len();
    if strat.total.len() != n_gq * n_dp {
        bail!("Internal error: gq×dp scan length mismatch");
    }

    let mut pts: Vec<(f64, f64, f64, i32, i32)> = Vec::with_capacity(n_gq * n_dp);
    for i in 0..n_gq {
        for j in 0..n_dp {
            let c = &strat.total[strat.idx(i, j)];
            let p = precision(c.tp, c.fp);
            let r = recall(c.tp, c.fn_);
            let f1v = f1(p, r);
            pts.push((r, p, f1v, scan.gq[i], scan.dp[j]));
        }
    }

    let all_x = pts.iter().map(|(r, _, _, _, _)| *r);
    let (xmin, xmax) = bounds01(all_x);
    let ymin = args.gq_dp_pr_precision_min;
    let ymax = args.gq_dp_pr_precision_max;
    if !(ymin.is_finite()
        && ymax.is_finite()
        && ymin < ymax
        && ymin >= 0.0
        && ymax <= 1.0)
    {
        bail!(
            "--gq-dp-pr-precision-min / --gq-dp-pr-precision-max must satisfy 0 <= min < max <= 1"
        );
    }

    let root = BitMapBackend::new(out_path, (900, 720)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut caption = pr_plot_caption(args);
    caption.push_str("\nGQ×DP PR scatter (total): each point = one threshold pair");
    caption.push_str("\nx=Recall, y=Precision; color=F1");
    caption.push_str(&format!(
        "\nPrecision y-axis: {:.3}..{:.3} (use --gq-dp-pr-precision-min/max to adjust)",
        ymin, ymax
    ));
    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .caption(caption, ("sans-serif", 15))
        .x_label_area_size(55)
        .y_label_area_size(60)
        .build_cartesian_2d(xmin..xmax, ymin..ymax)?;

    chart
        .configure_mesh()
        .x_desc("Recall = TP / (TP + FN)")
        .y_desc("Precision = TP / (TP + FP)")
        .draw()?;

    let r_px: i32 = if pts.len() > 400 {
        2
    } else if pts.len() > 150 {
        3
    } else {
        4
    };

    for (r, p, f1v, _, _) in &pts {
        let color = f1_to_rgb(*f1v).mix(0.78);
        chart.draw_series(std::iter::once(Circle::new(
            (*r, *p),
            r_px,
            color.filled(),
        )))?;
    }

    // Label top F1 threshold pairs (same semantics as `gq_dp_scan` total row).
    let mut order: Vec<usize> = (0..pts.len()).collect();
    order.sort_by(|&a, &b| pts[b].2.partial_cmp(&pts[a].2).unwrap_or(std::cmp::Ordering::Equal));

    let x_leg = xmin + (xmax - xmin) * 0.58;
    let mut y_leg = ymax - (ymax - ymin) * 0.06;
    let row_h = (ymax - ymin) * 0.042;
    chart.draw_series(std::iter::once(Text::new(
        "Top F1 (total):",
        (x_leg, y_leg),
        ("sans-serif", 12).into_font().color(&BLACK),
    )))?;
    y_leg -= row_h;
    for idx in order.iter().take(4) {
        let (_, _, f1v, gq_thr, dp_thr) = pts[*idx];
        let label = format!("GQ≥{gq_thr} DP≥{dp_thr}  F1={f1v:.3}");
        chart.draw_series(std::iter::once(Text::new(
            label,
            (x_leg, y_leg),
            ("sans-serif", 11).into_font().color(&BLACK),
        )))?;
        y_leg -= row_h;
    }

    root.present()?;
    Ok(())
}

fn pr_plot_caption(args: &EvalArgs) -> String {
    let bed_name = args
        .confident_bed
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("(confident bed path)");
    match &args.region {
        Some(r) => format!(
            "GQ Precision vs Recall (SNV/InDel)\nConfident BED: {bed_name}\nAlso restricted to: {r}"
        ),
        None => format!("GQ Precision vs Recall (SNV/InDel)\nConfident BED: {bed_name}"),
    }
}

fn write_gq_pr_png(
    scan: &GqScan,
    strat: &StratifiedScanCounts,
    out_path: &PathBuf,
    args: &EvalArgs,
) -> Result<()> {
    if strat.total.len() != scan.thresholds.len()
        || strat.snv.len() != scan.thresholds.len()
        || strat.insertion.len() != scan.thresholds.len()
        || strat.deletion.len() != scan.thresholds.len()
    {
        bail!("Internal error: scan thresholds and counts length mismatch");
    }

    let pts_total: Vec<(f64, f64)> = strat.total.iter().map(|c| pr_point(c)).collect();
    let pts_snv: Vec<(f64, f64)> = strat.snv.iter().map(|c| pr_point(c)).collect();
    let pts_ins: Vec<(f64, f64)> = strat.insertion.iter().map(|c| pr_point(c)).collect();
    let pts_del: Vec<(f64, f64)> = strat.deletion.iter().map(|c| pr_point(c)).collect();

    let all_x = pts_total
        .iter()
        .chain(pts_snv.iter())
        .chain(pts_ins.iter())
        .chain(pts_del.iter())
        .map(|(x, _)| *x);
    let all_y = pts_total
        .iter()
        .chain(pts_snv.iter())
        .chain(pts_ins.iter())
        .chain(pts_del.iter())
        .map(|(_, y)| *y);

    let (xmin, xmax) = bounds01(all_x);
    let (ymin, ymax) = bounds01(all_y);

    let root = BitMapBackend::new(out_path, (900, 700)).into_drawing_area();
    root.fill(&WHITE)?;

    let caption = pr_plot_caption(args);
    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .caption(caption, ("sans-serif", 22))
        .x_label_area_size(55)
        .y_label_area_size(60)
        .build_cartesian_2d(xmin..xmax, ymin..ymax)?;

    chart
        .configure_mesh()
        .x_desc("Recall = TP / (TP + FN)")
        .y_desc("Precision = TP / (TP + FP)")
        .draw()?;

    // Curves + small markers.
    chart.draw_series(LineSeries::new(pts_total.clone(), &BLACK.mix(0.55)))?;
    chart.draw_series(pts_total.iter().map(|&(x, y)| Circle::new((x, y), 2, BLACK.mix(0.55).filled())))?;

    chart.draw_series(LineSeries::new(pts_snv.clone(), &BLUE))?;
    chart.draw_series(pts_snv.iter().map(|&(x, y)| Circle::new((x, y), 2, BLUE.filled())))?;

    chart.draw_series(LineSeries::new(pts_ins.clone(), &GREEN))?;
    chart.draw_series(pts_ins.iter().map(|&(x, y)| Circle::new((x, y), 2, GREEN.filled())))?;

    chart.draw_series(LineSeries::new(pts_del.clone(), &MAGENTA))?;
    chart.draw_series(pts_del.iter().map(|&(x, y)| Circle::new((x, y), 2, MAGENTA.filled())))?;

    // Simple legend in the top-right corner.
    let x0 = xmin + (xmax - xmin) * 0.62;
    let mut y = ymax - (ymax - ymin) * 0.10;
    let row_h = (ymax - ymin) * 0.06;
    let dx = (xmax - xmin) * 0.02;

    chart.draw_series(std::iter::once(Circle::new((x0, y), 3, BLACK.mix(0.55).filled())))?;
    chart.draw_series(std::iter::once(Text::new(
        "total",
        (x0 + dx, y),
        ("sans-serif", 14).into_font().color(&BLACK),
    )))?;
    y -= row_h;
    chart.draw_series(std::iter::once(Circle::new((x0, y), 3, BLUE.filled())))?;
    chart.draw_series(std::iter::once(Text::new(
        "snv",
        (x0 + dx, y),
        ("sans-serif", 14).into_font().color(&BLACK),
    )))?;
    y -= row_h;
    chart.draw_series(std::iter::once(Circle::new((x0, y), 3, GREEN.filled())))?;
    chart.draw_series(std::iter::once(Text::new(
        "insertion",
        (x0 + dx, y),
        ("sans-serif", 14).into_font().color(&BLACK),
    )))?;
    y -= row_h;
    chart.draw_series(std::iter::once(Circle::new((x0, y), 3, MAGENTA.filled())))?;
    chart.draw_series(std::iter::once(Text::new(
        "deletion",
        (x0 + dx, y),
        ("sans-serif", 14).into_font().color(&BLACK),
    )))?;

    // Annotate a few points along the TOTAL curve so you can map them to thresholds.
    let n = scan.thresholds.len();
    let picks = [0usize, n / 2, n.saturating_sub(1)];
    for &i in &picks {
        if i >= n {
            continue;
        }
        let (x, y) = pts_total[i];
        let label = format!("GQ≥{}", scan.thresholds[i]);
        chart.draw_series(std::iter::once(Text::new(
            label,
            (x, y),
            ("sans-serif", 14).into_font().color(&BLACK),
        )))?;
    }

    root.present()?;
    Ok(())
}

fn pr_point(c: &ScanCounts) -> (f64, f64) {
    let r = if (c.tp + c.fn_) == 0 {
        0.0
    } else {
        c.tp as f64 / (c.tp + c.fn_) as f64
    };
    let p = if (c.tp + c.fp) == 0 {
        0.0
    } else {
        c.tp as f64 / (c.tp + c.fp) as f64
    };
    (r, p)
}

fn bounds01<I: IntoIterator<Item = f64>>(vals: I) -> (f64, f64) {
    let mut minv: f64 = 1.0;
    let mut maxv: f64 = 0.0;
    let mut any = false;
    for v in vals {
        let v = v.clamp(0.0, 1.0);
        minv = minv.min(v);
        maxv = maxv.max(v);
        any = true;
    }
    if !any {
        return (0.0, 1.0);
    }
    // Add a small pad so curves aren't glued to borders.
    let pad: f64 = 0.02;
    let lo = (minv - pad).clamp(0.0, 1.0);
    let hi = (maxv + pad).clamp(0.0, 1.0);
    if (hi - lo) < 1e-6 {
        (0.0, 1.0)
    } else {
        (lo, hi)
    }
}

fn emit_allele_table(s: &StratifiedAlleleCounts) {
    eprintln!("## allele\tclass\ttp\tfp\tfn\tprecision\trecall\tf1\tcompared\ttruth_nocall\tquery_nocall");
    emit_allele_row("total", &s.total);
    emit_allele_row("snv", &s.snv);
    emit_allele_row("insertion", &s.insertion);
    emit_allele_row("deletion", &s.deletion);

    if !s.insertion_by_len.is_empty() || !s.deletion_by_len.is_empty() {
        eprintln!("## allele_indel_len\ttype\tlen\ttp\tfp\tfn\tprecision\trecall\tf1\tcompared\ttruth_nocall\tquery_nocall");
        for (l, c) in &s.insertion_by_len {
            emit_allele_len_row("insertion", *l, c);
        }
        for (l, c) in &s.deletion_by_len {
            emit_allele_len_row("deletion", *l, c);
        }
    }
}

fn emit_allele_row(label: &str, c: &AlleleCounts) {
    let prec = precision(c.tp, c.fp);
    let rec = recall(c.tp, c.fn_);
    let f1 = f1(prec, rec);
    eprintln!(
        "allele\t{label}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}",
        c.tp,
        c.fp,
        c.fn_,
        prec,
        rec,
        f1,
        c.compared,
        c.truth_nocall,
        c.query_nocall
    );
}

fn emit_allele_len_row(t: &str, len: u64, c: &AlleleCounts) {
    let prec = precision(c.tp, c.fp);
    let rec = recall(c.tp, c.fn_);
    let f1 = f1(prec, rec);
    eprintln!(
        "allele_indel_len\t{t}\t{len}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}",
        c.tp, c.fp, c.fn_, prec, rec, f1, c.compared, c.truth_nocall, c.query_nocall
    );
}

fn emit_genotype_table(s: &StratifiedGenotypeCounts) {
    eprintln!("## genotype\ttype\tlen\texact\tcompared\texact_rate\ttruth_nocall\tquery_nocall");
    emit_genotype_len_row("snv", 0, &s.snv);
    for (l, c) in &s.insertion_by_len {
        emit_genotype_len_row("insertion", *l, c);
    }
    for (l, c) in &s.deletion_by_len {
        emit_genotype_len_row("deletion", *l, c);
    }
}

fn emit_genotype_len_row(t: &str, len: u64, c: &GenotypeCounts) {
    let rate = if c.compared == 0 { 0.0 } else { c.exact as f64 / c.compared as f64 };
    eprintln!(
        "genotype\t{t}\t{len}\t{}\t{}\t{:.6}\t{}\t{}",
        c.exact, c.compared, rate, c.truth_nocall, c.query_nocall
    );
}

fn precision(tp: u64, fp: u64) -> f64 {
    let denom = tp + fp;
    if denom == 0 { 0.0 } else { tp as f64 / denom as f64 }
}
fn recall(tp: u64, fn_: u64) -> f64 {
    let denom = tp + fn_;
    if denom == 0 { 0.0 } else { tp as f64 / denom as f64 }
}
fn f1(p: f64, r: f64) -> f64 {
    if p == 0.0 && r == 0.0 { 0.0 } else { 2.0 * p * r / (p + r) }
}

fn score_side_only(
    allele: &mut StratifiedAlleleCounts,
    gt: &mut StratifiedGenotypeCounts,
    v: &Biallelic,
    is_truth: bool,
    args: &EvalArgs,
) {
    // v.gt encodes ALT copies (0..2) for this ALT. None means nocall.
    if v.gt.is_none() {
        if is_truth {
            inc_allele_truth_nocall(allele, v.vc, v.indel_len);
            inc_gt_truth_nocall(gt, v.vc, v.indel_len);
        } else {
            inc_allele_query_nocall(allele, v.vc, v.indel_len);
            inc_gt_query_nocall(gt, v.vc, v.indel_len);
        }
        return;
    }
    let (c1, c2) = v.gt.unwrap();
    let alt_copies = c1 as u64 + c2 as u64;
    if args.require_sample_nonref && alt_copies == 0 {
        return;
    }
    inc_allele_compared(allele, v.vc, v.indel_len);
    if is_truth {
        inc_allele_fn(allele, v.vc, v.indel_len);
    } else {
        inc_allele_fp(allele, v.vc, v.indel_len);
    }
}

fn score_both(
    allele: &mut StratifiedAlleleCounts,
    gt: &mut StratifiedGenotypeCounts,
    t: &Biallelic,
    q: &Biallelic,
    args: &EvalArgs,
) {
    let vc = t.vc;

    let t_gt = t.gt;
    let q_gt = q.gt;

    if t_gt.is_none() {
        inc_allele_truth_nocall(allele, vc, t.indel_len);
        inc_gt_truth_nocall(gt, vc, t.indel_len);
    }
    if q_gt.is_none() {
        inc_allele_query_nocall(allele, vc, t.indel_len);
        inc_gt_query_nocall(gt, vc, t.indel_len);
    }

    if args.exclude_nocall && (t_gt.is_none() || q_gt.is_none()) {
        return;
    }

    let t_alt = t_gt.map(|(a, b)| a as u64 + b as u64).unwrap_or(0);
    let q_alt = q_gt.map(|(a, b)| a as u64 + b as u64).unwrap_or(0);

    if args.require_sample_nonref && t_alt == 0 && q_alt == 0 {
        return;
    }

    inc_allele_compared(allele, vc, t.indel_len);
    match (t_alt > 0, q_alt > 0) {
        (true, true) => {
            inc_allele_tp(allele, vc, t.indel_len);
        }
        (false, true) => {
            inc_allele_fp(allele, vc, t.indel_len);
        }
        (true, false) => {
            inc_allele_fn(allele, vc, t.indel_len);
        }
        (false, false) => {}
    }

    inc_gt_compared(gt, vc, t.indel_len);
    if t_alt == q_alt {
        inc_gt_exact(gt, vc, t.indel_len);
    }
}

fn allele_class_bucket_mut<'a>(s: &'a mut StratifiedAlleleCounts, vc: VarClass) -> &'a mut AlleleCounts {
    match vc {
        VarClass::Snv => &mut s.snv,
        VarClass::Insertion => &mut s.insertion,
        VarClass::Deletion => &mut s.deletion,
    }
}

fn allele_len_bucket_mut<'a>(
    s: &'a mut StratifiedAlleleCounts,
    vc: VarClass,
    indel_len: Option<u64>,
) -> Option<&'a mut AlleleCounts> {
    match (vc, indel_len) {
        (VarClass::Insertion, Some(l)) => Some(s.insertion_by_len.entry(l).or_default()),
        (VarClass::Deletion, Some(l)) => Some(s.deletion_by_len.entry(l).or_default()),
        _ => None,
    }
}

fn genotype_class_bucket_mut<'a>(s: &'a mut StratifiedGenotypeCounts, vc: VarClass) -> &'a mut GenotypeCounts {
    match vc {
        VarClass::Snv => &mut s.snv,
        VarClass::Insertion => &mut s.insertion,
        VarClass::Deletion => &mut s.deletion,
    }
}

fn genotype_len_bucket_mut<'a>(
    s: &'a mut StratifiedGenotypeCounts,
    vc: VarClass,
    indel_len: Option<u64>,
) -> Option<&'a mut GenotypeCounts> {
    match (vc, indel_len) {
        (VarClass::Insertion, Some(l)) => Some(s.insertion_by_len.entry(l).or_default()),
        (VarClass::Deletion, Some(l)) => Some(s.deletion_by_len.entry(l).or_default()),
        _ => None,
    }
}

fn inc_allele_compared(s: &mut StratifiedAlleleCounts, vc: VarClass, indel_len: Option<u64>) {
    s.total.compared += 1;
    allele_class_bucket_mut(s, vc).compared += 1;
    if let Some(c) = allele_len_bucket_mut(s, vc, indel_len) {
        c.compared += 1;
    }
}
fn inc_allele_tp(s: &mut StratifiedAlleleCounts, vc: VarClass, indel_len: Option<u64>) {
    s.total.tp += 1;
    allele_class_bucket_mut(s, vc).tp += 1;
    if let Some(c) = allele_len_bucket_mut(s, vc, indel_len) {
        c.tp += 1;
    }
}
fn inc_allele_fp(s: &mut StratifiedAlleleCounts, vc: VarClass, indel_len: Option<u64>) {
    s.total.fp += 1;
    allele_class_bucket_mut(s, vc).fp += 1;
    if let Some(c) = allele_len_bucket_mut(s, vc, indel_len) {
        c.fp += 1;
    }
}
fn inc_allele_fn(s: &mut StratifiedAlleleCounts, vc: VarClass, indel_len: Option<u64>) {
    s.total.fn_ += 1;
    allele_class_bucket_mut(s, vc).fn_ += 1;
    if let Some(c) = allele_len_bucket_mut(s, vc, indel_len) {
        c.fn_ += 1;
    }
}
fn inc_allele_truth_nocall(s: &mut StratifiedAlleleCounts, vc: VarClass, indel_len: Option<u64>) {
    s.total.truth_nocall += 1;
    allele_class_bucket_mut(s, vc).truth_nocall += 1;
    if let Some(c) = allele_len_bucket_mut(s, vc, indel_len) {
        c.truth_nocall += 1;
    }
}
fn inc_allele_query_nocall(s: &mut StratifiedAlleleCounts, vc: VarClass, indel_len: Option<u64>) {
    s.total.query_nocall += 1;
    allele_class_bucket_mut(s, vc).query_nocall += 1;
    if let Some(c) = allele_len_bucket_mut(s, vc, indel_len) {
        c.query_nocall += 1;
    }
}

fn inc_gt_compared(s: &mut StratifiedGenotypeCounts, vc: VarClass, indel_len: Option<u64>) {
    s.total.compared += 1;
    genotype_class_bucket_mut(s, vc).compared += 1;
    if let Some(c) = genotype_len_bucket_mut(s, vc, indel_len) {
        c.compared += 1;
    }
}
fn inc_gt_exact(s: &mut StratifiedGenotypeCounts, vc: VarClass, indel_len: Option<u64>) {
    s.total.exact += 1;
    genotype_class_bucket_mut(s, vc).exact += 1;
    if let Some(c) = genotype_len_bucket_mut(s, vc, indel_len) {
        c.exact += 1;
    }
}
fn inc_gt_truth_nocall(s: &mut StratifiedGenotypeCounts, vc: VarClass, indel_len: Option<u64>) {
    s.total.truth_nocall += 1;
    genotype_class_bucket_mut(s, vc).truth_nocall += 1;
    if let Some(c) = genotype_len_bucket_mut(s, vc, indel_len) {
        c.truth_nocall += 1;
    }
}
fn inc_gt_query_nocall(s: &mut StratifiedGenotypeCounts, vc: VarClass, indel_len: Option<u64>) {
    s.total.query_nocall += 1;
    genotype_class_bucket_mut(s, vc).query_nocall += 1;
    if let Some(c) = genotype_len_bucket_mut(s, vc, indel_len) {
        c.query_nocall += 1;
    }
}

fn truth_sample_name(header: &rust_htslib::bcf::header::HeaderView, name: Option<&str>) -> Result<String> {
    let samples = header.samples();
    if samples.is_empty() {
        bail!("Truth VCF has no samples");
    }
    if let Some(name) = name {
        return Ok(name.to_string());
    }
    if samples.len() == 1 {
        return Ok(String::from_utf8_lossy(samples[0]).to_string());
    }
    bail!(
        "Truth VCF has {} samples; pass --truth-sample. First few: {:?}",
        samples.len(),
        samples
            .iter()
            .take(5)
            .map(|x| String::from_utf8_lossy(x).to_string())
            .collect::<Vec<_>>()
    );
}

fn query_sample_name(header: &rust_htslib::bcf::header::HeaderView, name: Option<&str>) -> Result<String> {
    let samples = header.samples();
    if samples.is_empty() {
        bail!("Query VCF has no samples");
    }
    if let Some(name) = name {
        return Ok(name.to_string());
    }
    Ok(String::from_utf8_lossy(samples[0]).to_string())
}

fn sample_index_required(header: &rust_htslib::bcf::header::HeaderView, name: &str) -> Result<usize> {
    let samples = header.samples();
    if samples.is_empty() {
        bail!("VCF has no samples");
    }
    for (i, s) in samples.iter().enumerate() {
        if *s == name.as_bytes() {
            return Ok(i);
        }
    }
    Err(anyhow!(
        "Sample {name:?} not found. First few: {:?}",
        samples
            .iter()
            .take(5)
            .map(|x| String::from_utf8_lossy(x).to_string())
            .collect::<Vec<_>>()
    ))
}

#[derive(Clone)]
struct BedMask {
    // chrom -> sorted, merged intervals [start0, end0)
    map: HashMap<String, Vec<(u64, u64)>>,
}

impl BedMask {
    fn from_path(path: &PathBuf) -> Result<Self> {
        let f = File::open(path)?;
        let r = BufReader::new(f);
        let mut map: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        for line in r.lines() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 3 {
                continue;
            }
            let chrom = parts[0].to_string();
            let start: u64 = parts[1].parse()?;
            let end: u64 = parts[2].parse()?;
            if end <= start {
                continue;
            }
            map.entry(chrom).or_default().push((start, end));
        }
        for v in map.values_mut() {
            v.sort_by_key(|x| x.0);
            // merge overlaps
            let mut merged: Vec<(u64, u64)> = Vec::with_capacity(v.len());
            for (s, e) in v.drain(..) {
                if let Some(last) = merged.last_mut() {
                    if s <= last.1 {
                        last.1 = last.1.max(e);
                        continue;
                    }
                }
                merged.push((s, e));
            }
            *v = merged;
        }
        Ok(Self { map })
    }

    fn contains_pos1(&self, chrom: &str, pos1: i64) -> bool {
        let pos0 = (pos1 - 1) as u64;
        let Some(v) = self.map.get(chrom) else { return false };
        // binary search interval containing pos0
        let mut lo = 0usize;
        let mut hi = v.len();
        while lo < hi {
            let mid = (lo + hi) / 2;
            let (s, e) = v[mid];
            if pos0 < s {
                hi = mid;
            } else if pos0 >= e {
                lo = mid + 1;
            } else {
                return true;
            }
        }
        false
    }
}

struct BiallelicIter<'a> {
    reader: &'a mut bcf::Reader,
    sample_i: usize,
    bed: BedMask,
    region: Option<Region>,
    gq_tag: Option<Vec<u8>>,
    dp_tag: Option<Vec<u8>>,
    // buffered record and next alt index to emit
    rec: Option<bcf::Record>,
    next_alt_i: usize,
}

impl<'a> BiallelicIter<'a> {
    fn new(
        reader: &'a mut bcf::Reader,
        sample_i: usize,
        bed: BedMask,
        region: Option<Region>,
        gq_tag: Option<Vec<u8>>,
        dp_tag: Option<Vec<u8>>,
    ) -> Result<Self> {
        Ok(Self {
            reader,
            sample_i,
            bed,
            region,
            gq_tag,
            dp_tag,
            rec: None,
            next_alt_i: 1,
        })
    }

    fn next(&mut self) -> Result<Option<Biallelic>> {
        loop {
            if self.rec.is_none() {
                match self.reader.records().next() {
                    None => return Ok(None),
                    Some(r) => {
                        self.rec = Some(r?);
                        self.next_alt_i = 1;
                    }
                }
            }
            let rec = self.rec.as_ref().unwrap();
            let alleles = rec.alleles();
            if alleles.len() <= 1 {
                self.rec = None;
                continue;
            }
            if self.next_alt_i >= alleles.len() {
                self.rec = None;
                continue;
            }

            let rid = rec.rid().ok_or_else(|| anyhow!("Record missing RID/CHROM"))?;
            let chrom = String::from_utf8_lossy(self.reader.header().rid2name(rid)?).to_string();
            let pos1 = (rec.pos() as i64) + 1;

            if let Some(region) = &self.region {
                if chrom != region.chrom {
                    // If input is sorted, we can early-skip by chrom string ordering.
                    // But we keep it simple: just skip records outside desired chrom.
                    self.rec = None;
                    continue;
                }
                if pos1 < region.start1 || pos1 > region.end1 {
                    self.rec = None;
                    continue;
                }
            }

            if !self.bed.contains_pos1(&chrom, pos1) {
                // skip whole record if its POS isn't inside confident BED.
                self.rec = None;
                continue;
            }

            let r = alleles[0].to_vec();
            let a = alleles[self.next_alt_i].to_vec();
            self.next_alt_i += 1;

            // Skip symbolic/non-ref alts.
            if a.starts_with(b"<") || a == b"*" || a.contains(&b'[') || a.contains(&b']') {
                continue;
            }

            let (pos1n, rn, an) = normalize_no_ref(pos1, &r, &a);
            let gt = gt_alt_copies(rec, self.sample_i, self.next_alt_i - 1)?;
            let vc = variant_class(&rn, &an);
            let gq = if let Some(tag) = &self.gq_tag {
                extract_sample_format_i32(rec, self.sample_i, tag)?
            } else {
                None
            };
            let dp = if let Some(tag) = &self.dp_tag {
                extract_sample_format_i32(rec, self.sample_i, tag)?
            } else {
                None
            };

            return Ok(Some(Biallelic {
                key: Key { chrom, pos1: pos1n, r: rn, a: an },
                gt,
                vc,
                indel_len: indel_len(vc, &r, &a),
                gq,
                dp,
            }));
        }
    }
}

fn extract_sample_format_i32(rec: &bcf::Record, sample_i: usize, tag: &[u8]) -> Result<Option<i32>> {
    let buf = match rec.format(tag).integer() {
        Ok(b) => b,
        Err(_) => return Ok(None),
    };
    // FORMAT integer buffers are [sample][value_index].
    let Some(sample_vals) = buf.get(sample_i) else {
        return Ok(None);
    };
    let Some(v) = sample_vals.get(0) else {
        return Ok(None);
    };
    if *v == i32::MIN { Ok(None) } else { Ok(Some(*v)) }
}

fn gt_alt_copies(rec: &bcf::Record, sample_i: usize, alt_allele_index: usize) -> Result<Option<(u8, u8)>> {
    let gts = rec.genotypes()?;
    let gt = gts.get(sample_i);
    if gt.len() < 2 {
        return Ok(None);
    }
    let a1 = allele_index(gt[0]);
    let a2 = allele_index(gt[1]);
    if a1.is_none() || a2.is_none() {
        return Ok(None);
    }
    let a1 = a1.unwrap();
    let a2 = a2.unwrap();
    let c1 = (a1 == alt_allele_index as i32) as u8;
    let c2 = (a2 == alt_allele_index as i32) as u8;
    Ok(Some((c1, c2)))
}

fn allele_index(a: GenotypeAllele) -> Option<i32> {
    match a {
        GenotypeAllele::Unphased(i) | GenotypeAllele::Phased(i) => Some(i),
        GenotypeAllele::UnphasedMissing | GenotypeAllele::PhasedMissing => None,
    }
}

// Minimal normalization that does not require a reference sequence:
// - trim common suffix
// - trim common prefix (but keep at least 1 base in each allele), incrementing POS accordingly
fn normalize_no_ref(pos1: i64, r: &[u8], a: &[u8]) -> (i64, Vec<u8>, Vec<u8>) {
    let mut r = r.to_vec();
    let mut a = a.to_vec();
    let mut pos1 = pos1;

    // trim common suffix
    while r.len() > 1 && a.len() > 1 {
        if r.last() == a.last() {
            r.pop();
            a.pop();
        } else {
            break;
        }
    }

    // trim common prefix
    while r.len() > 1 && a.len() > 1 {
        if r[0] == a[0] {
            r.remove(0);
            a.remove(0);
            pos1 += 1;
        } else {
            break;
        }
    }

    (pos1, r, a)
}

fn variant_class(r: &[u8], a: &[u8]) -> VarClass {
    if r.len() == 1 && a.len() == 1 {
        VarClass::Snv
    } else if a.len() > r.len() {
        VarClass::Insertion
    } else {
        // Includes deletions and same-length MNVs; treat same-length non-SNVs as deletion? No.
        // We keep the three requested buckets by mapping all non-insertions to deletion unless it's a single-base SNV.
        // In typical dipcall-normalized truth/query for small variants, same-length non-SNVs should be rare.
        VarClass::Deletion
    }
}

fn indel_len(vc: VarClass, r: &[u8], a: &[u8]) -> Option<u64> {
    match vc {
        VarClass::Insertion => Some((a.len().saturating_sub(r.len())) as u64),
        VarClass::Deletion => Some((r.len().saturating_sub(a.len())) as u64),
        VarClass::Snv => None,
    }
}

#[derive(Debug, Clone)]
struct Region {
    chrom: String,
    start1: i64,
    end1: i64,
}

fn parse_region(s: &str) -> Result<Region> {
    // Expected: chrom:start-end (1-based, inclusive). This matches common CLI usage.
    let (chrom, rest) = s
        .split_once(':')
        .ok_or_else(|| anyhow!("Region must look like chr20:100-200, got {s:?}"))?;
    let (start_s, end_s) = rest
        .split_once('-')
        .ok_or_else(|| anyhow!("Region must look like chr20:100-200, got {s:?}"))?;
    let start1: i64 = start_s.parse()?;
    let end1: i64 = end_s.parse()?;
    if start1 <= 0 || end1 <= 0 || end1 < start1 {
        bail!("Invalid region range: {s:?}");
    }
    Ok(Region { chrom: chrom.to_string(), start1, end1 })
}

