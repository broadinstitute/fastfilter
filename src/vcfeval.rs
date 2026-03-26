use std::path::PathBuf;
use std::process::{Command, Stdio};

use anyhow::{bail, Context, Result};
use clap::Args;

/// Run [RTG vcfeval](https://github.com/RealTimeGenomics/rtg-tools) for baseline-vs-calls comparison
/// (representation-aware matching). Requires a reference genome in RTG **SDF** form (`rtg format` from FASTA).
#[derive(Debug, Args, Clone)]
pub struct VcfevalArgs {
    /// Path to `RTG.jar`
    #[arg(long, default_value = "scratch/rtg-tools-3.12.1/RTG.jar")]
    pub rtg_jar: PathBuf,

    /// Baseline / truth VCF (`vcfeval -b`)
    #[arg(long, alias = "truth")]
    pub baseline: PathBuf,

    /// Calls / query VCF (`vcfeval -c`)
    #[arg(long, alias = "query")]
    pub calls: PathBuf,

    /// Reference template SDF (`vcfeval -t`). Build with e.g. `rtg format -o ref.sdf ref.fa`.
    #[arg(long, short = 't')]
    pub template: PathBuf,

    /// Output directory for vcfeval outputs (`vcfeval -o`)
    #[arg(long, short = 'o')]
    pub output: PathBuf,

    /// Only read VCF records overlapping this BED (`vcfeval --bed-regions`)
    #[arg(long)]
    pub bed_regions: Option<PathBuf>,

    /// Evaluate within these regions (`vcfeval --evaluation-regions`). Typical for truth high-confidence BEDs.
    #[arg(long = "evaluation-regions", visible_alias = "confident-bed")]
    pub evaluation_regions: Option<PathBuf>,

    /// Restrict input loci (`vcfeval --region`), e.g. `chr20:10000000-11000000`
    #[arg(long)]
    pub region: Option<String>,

    /// Sample name(s) for multi-sample VCFs (`vcfeval --sample`). Use `baseline_sample,calls_sample` if names differ.
    #[arg(long)]
    pub sample: Option<String>,

    /// Pass all variants including non-PASS (`vcfeval --all-records`)
    #[arg(long, default_value_t = false)]
    pub all_records: bool,

    /// Thread count for RTG (`vcfeval -T`)
    #[arg(long, short = 'T')]
    pub rtg_threads: Option<usize>,

    /// Extra arguments forwarded to `vcfeval` (repeat flag), e.g. `--vcfeval-arg --squash-ploidy`
    #[arg(long = "vcfeval-arg", value_name = "ARG")]
    pub vcfeval_arg: Vec<String>,
}

pub fn run_vcfeval(args: VcfevalArgs) -> Result<()> {
    if !args.rtg_jar.is_file() {
        bail!("RTG jar not found (expected a file): {:?}", args.rtg_jar);
    }
    if !args.template.exists() {
        bail!("Reference template SDF not found: {:?}", args.template);
    }

    let mut cmd = Command::new("java");
    cmd.arg("-jar")
        .arg(&args.rtg_jar)
        .arg("vcfeval")
        .arg("-b")
        .arg(&args.baseline)
        .arg("-c")
        .arg(&args.calls)
        .arg("-o")
        .arg(&args.output)
        .arg("-t")
        .arg(&args.template)
        .stdin(Stdio::null());

    if let Some(p) = &args.bed_regions {
        cmd.arg("--bed-regions").arg(p);
    }
    if let Some(p) = &args.evaluation_regions {
        cmd.arg("--evaluation-regions").arg(p);
    }
    if let Some(r) = &args.region {
        cmd.arg("--region").arg(r);
    }
    if let Some(s) = &args.sample {
        cmd.arg("--sample").arg(s);
    }
    if args.all_records {
        cmd.arg("--all-records");
    }
    if let Some(n) = args.rtg_threads {
        cmd.arg("-T").arg(n.to_string());
    }
    for a in &args.vcfeval_arg {
        cmd.arg(a);
    }

    eprintln!("## Running: {:?}", cmd);
    let status = cmd
        .status()
        .with_context(|| format!("Failed to spawn `java` (is it on PATH?) for {:?}", cmd))?;
    if !status.success() {
        bail!("vcfeval exited with {:?}", status.code());
    }
    Ok(())
}
