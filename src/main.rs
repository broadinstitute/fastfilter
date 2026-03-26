use std::path::PathBuf;

use anyhow::{anyhow, Context, Result};
use clap::{Args, Parser, Subcommand};
use rayon::ThreadPoolBuilder;
use rust_htslib::bcf::{self, Format};
use rust_htslib::bcf::Read;

mod eval;
mod filter;
mod vcfeval;

#[derive(Debug, Parser)]
#[command(
    name = "fastfilter",
    version,
    about = "Fast per-genotype masking for DeepVariant+GLNexus joint-callsets (thresholds or gbdt-rs)."
)]
struct Cli {
    #[command(subcommand)]
    command: Command,

    /// Number of worker threads for record-batch parallelism.
    #[arg(long, default_value_t = num_cpus::get())]
    threads: usize,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    /// Filter/mask per-sample genotypes in a multi-sample VCF/BCF.
    Filter(FilterArgs),
    /// Evaluate a query VCF against a truth VCF (SNVs + small indels).
    Eval(eval::EvalArgs),
    /// Run RTG Tools `vcfeval` (representation-aware; requires reference SDF).
    Vcfeval(vcfeval::VcfevalArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .ok();

    match cli.command {
        Command::Filter(args) => run_filter(args),
        Command::Eval(args) => eval::run_eval(args),
        Command::Vcfeval(args) => vcfeval::run_vcfeval(args),
    }
}

#[derive(Debug, Args)]
pub struct FilterArgs {
    /// Input multi-sample VCF/BCF (typically .vcf.gz)
    #[arg(long)]
    input: PathBuf,

    /// Output filtered VCF (or .vcf.gz if you pass a path ending in .gz)
    #[arg(long)]
    output: PathBuf,

    /// Optional region: e.g. chr20:20000000-21000000
    /// Region extraction requires an indexed input and uses `tabix` (no streaming fallback).
    #[arg(long)]
    region: Option<String>,

    /// How to decide failing genotypes.
    #[command(subcommand)]
    mode: Mode,

    /// Batch size in number of VCF records (increases throughput).
    #[arg(long, default_value_t = 1000)]
    batch_size: usize,

    /// Maximum ploidy to write back when masking genotypes.
    #[arg(long, default_value_t = 2)]
    max_ploidy: usize,
}

#[derive(Debug, Subcommand, Clone)]
pub enum Mode {
    /// Threshold mode: mask genotype if GQ/DP do not meet cutoffs.
    Thresholds {
        /// Minimum genotype quality (per sample GQ in FORMAT).
        #[arg(long)]
        min_gq: Option<i32>,

        /// Minimum read depth (per sample DP in FORMAT).
        #[arg(long)]
        min_dp: Option<i32>,

        /// FORMAT tag for GQ.
        #[arg(long, default_value = "GQ")]
        gq_tag: String,

        /// FORMAT tag for DP.
        #[arg(long, default_value = "DP")]
        dp_tag: String,
    },

    /// gbdt-rs mode: mask genotype if P(correct call) for ANY alt allele instance is < threshold.
    ///
    /// Your training notebook uses separate models per variant_class:
    /// - SNV
    /// - insertion (alt_len > ref_len)
    /// - deletion (alt_len < ref_len)
    Gbdt {
        /// Mask when predicted P(correct call) < model-threshold.
        #[arg(long)]
        model_threshold: f64,

        /// gbdt-rs model for SNV variant_class (converted by examples/convert_xgboost.py).
        #[arg(long)]
        gbdt_model_snv: PathBuf,

        /// gbdt-rs model for insertion variant_class (converted by examples/convert_xgboost.py).
        #[arg(long)]
        gbdt_model_insertion: PathBuf,

        /// gbdt-rs model for deletion variant_class (converted by examples/convert_xgboost.py).
        #[arg(long)]
        gbdt_model_deletion: PathBuf,

        /// Objective used when loading converted xgboost models.
        /// For binary classification trained with predict_proba[:,1], use `binary:logistic`.
        #[arg(long, default_value = "binary:logistic")]
        gbdt_objective: String,
    },
}

fn run_filter(args: FilterArgs) -> Result<()> {
    let (tmp_guard, input_path_for_reader) = prepare_input(&args.input, args.region.as_deref())?;

    let mut reader = bcf::Reader::from_path(&input_path_for_reader)
        .with_context(|| format!("Failed to open input: {input_path_for_reader:?}"))?;

    let sample_names: Vec<String> = reader
        .header()
        .samples()
        .iter()
        .map(|s| String::from_utf8_lossy(s).to_string())
        .collect();
    let stats = std::sync::Arc::new(filter::Stats::new(sample_names));

    // Writer requires a full Header (not just HeaderView).
    let header = bcf::header::Header::from_template(reader.header());
    let mut writer = bcf::Writer::from_path(&args.output, &header, false, Format::Vcf)
        .with_context(|| format!("Failed to create writer: {:?}", args.output))?;

    let filter_ctx = filter::FilterCtx::new(args.mode, args.max_ploidy)?;

    let batch_size = args.batch_size;
    let mut records: Vec<bcf::Record> = Vec::with_capacity(batch_size);

    for result in reader.records() {
        let record = result?;
        records.push(record);
        if records.len() >= batch_size {
            let processed = filter::process_batch(records, &filter_ctx, Some(stats.clone()))?;
            for rec in processed {
                writer.write(&rec)?;
            }
            records = Vec::with_capacity(batch_size);
        }
    }

    if !records.is_empty() {
        let processed = filter::process_batch(records, &filter_ctx, Some(stats.clone()))?;
        for rec in processed {
            writer.write(&rec)?;
        }
    }

    eprintln!("## Per-sample genotype stats");
    eprintln!("sample\tseen_genotypes\tmasked_genotypes\tmasked_rate");
    for (name, seen, masked) in stats.snapshot() {
        let rate = if seen == 0 { 0.0 } else { (masked as f64) / (seen as f64) };
        eprintln!("{name}\t{seen}\t{masked}\t{rate:.6}");
    }

    drop(tmp_guard);
    Ok(())
}

pub fn prepare_input(input: &PathBuf, region: Option<&str>) -> Result<(Option<tempfile::NamedTempFile>, PathBuf)> {
    if let Some(region) = region {
        which::which("tabix").context("`--region` requires `tabix` on PATH")?;

        let mut tmp = tempfile::NamedTempFile::new().context("Failed to create temp file")?;
        let tmp_path = tmp.path().to_path_buf();

        // `tabix -h` includes the VCF header so rust-htslib can parse the result standalone.
        let output = std::process::Command::new("tabix")
            .arg("-h")
            .arg(input)
            .arg(region)
            .output()
            .with_context(|| format!("Failed to run `tabix -h {:?} {region}`", input))?;

        if !output.status.success() {
            return Err(anyhow!(
                "tabix failed for region {region}. stderr: {}",
                String::from_utf8_lossy(&output.stderr)
            ));
        }

        use std::io::Write;
        tmp.write_all(&output.stdout)?;
        tmp.flush()?;

        Ok((Some(tmp), tmp_path))
    } else {
        Ok((None, input.clone()))
    }
}
