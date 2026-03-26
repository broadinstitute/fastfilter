use std::collections::HashSet;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use anyhow::{anyhow, bail, Context, Result};
use rayon::prelude::*;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf;

use crate::Mode;

use gbdt::decision_tree::Data;
use gbdt::gradient_boost::GBDT;

#[derive(Clone)]
pub struct ThresholdCfg {
    pub min_gq: Option<i32>,
    pub min_dp: Option<i32>,
    pub gq_tag: Vec<u8>,
    pub dp_tag: Vec<u8>,
}

#[derive(Clone)]
pub struct GbdtCfg {
    pub model_threshold: f64,
    pub snv: Arc<GBDT>,
    pub insertion: Arc<GBDT>,
    pub deletion: Arc<GBDT>,
}

#[derive(Clone)]
pub enum FilterMode {
    Thresholds(ThresholdCfg),
    Gbdt(GbdtCfg),
}

#[derive(Clone)]
pub struct FilterCtx {
    pub mode: FilterMode,
    pub max_ploidy: usize,
}

pub struct Stats {
    pub sample_names: Vec<String>,
    seen_genotypes: Vec<AtomicU64>,
    masked_genotypes: Vec<AtomicU64>,
}

impl Stats {
    pub fn new(sample_names: Vec<String>) -> Self {
        let n = sample_names.len();
        Self {
            sample_names,
            seen_genotypes: (0..n).map(|_| AtomicU64::new(0)).collect(),
            masked_genotypes: (0..n).map(|_| AtomicU64::new(0)).collect(),
        }
    }

    pub fn seen(&self, sample_i: usize) {
        self.seen_genotypes[sample_i].fetch_add(1, Ordering::Relaxed);
    }

    pub fn masked(&self, sample_i: usize) {
        self.masked_genotypes[sample_i].fetch_add(1, Ordering::Relaxed);
    }

    pub fn snapshot(&self) -> Vec<(String, u64, u64)> {
        let mut out = Vec::with_capacity(self.sample_names.len());
        for i in 0..self.sample_names.len() {
            out.push((
                self.sample_names[i].clone(),
                self.seen_genotypes[i].load(Ordering::Relaxed),
                self.masked_genotypes[i].load(Ordering::Relaxed),
            ));
        }
        out
    }
}

impl FilterCtx {
    pub fn new(cli_mode: Mode, max_ploidy: usize) -> Result<Self> {
        if max_ploidy == 0 {
            bail!("--max-ploidy must be >= 1");
        }

        let mode = match cli_mode {
            Mode::Thresholds {
                min_gq,
                min_dp,
                gq_tag,
                dp_tag,
            } => {
                if min_gq.is_none() && min_dp.is_none() {
                    bail!("Thresholds mode requires at least one of --min-gq or --min-dp");
                }
                FilterMode::Thresholds(ThresholdCfg {
                    min_gq,
                    min_dp,
                    gq_tag: gq_tag.into_bytes(),
                    dp_tag: dp_tag.into_bytes(),
                })
            }
            Mode::Gbdt {
                model_threshold,
                gbdt_model_snv,
                gbdt_model_insertion,
                gbdt_model_deletion,
                gbdt_objective,
            } => {
                let objective = gbdt_objective;
                let snv = Arc::new(
                    GBDT::from_xgboost_dump(
                        gbdt_model_snv
                            .to_str()
                            .ok_or_else(|| anyhow!("Invalid UTF-8 in --gbdt-model-snv path"))?,
                        &objective,
                    )
                    .context("Failed to load SNV gbdt model")?,
                );
                let insertion = Arc::new(
                    GBDT::from_xgboost_dump(
                        gbdt_model_insertion
                            .to_str()
                            .ok_or_else(|| anyhow!("Invalid UTF-8 in --gbdt-model-insertion path"))?,
                        &objective,
                    )
                    .context("Failed to load insertion gbdt model")?,
                );
                let deletion = Arc::new(
                    GBDT::from_xgboost_dump(
                        gbdt_model_deletion
                            .to_str()
                            .ok_or_else(|| anyhow!("Invalid UTF-8 in --gbdt-model-deletion path"))?,
                        &objective,
                    )
                    .context("Failed to load deletion gbdt model")?,
                );
                FilterMode::Gbdt(GbdtCfg {
                    model_threshold,
                    snv,
                    insertion,
                    deletion,
                })
            }
        };

        Ok(Self { mode, max_ploidy })
    }
}

pub fn process_batch(
    records: Vec<bcf::Record>,
    ctx: &FilterCtx,
    stats: Option<Arc<Stats>>,
) -> Result<Vec<bcf::Record>> {
    let ctx = ctx.clone();
    let stats = stats.clone();
    let out = records
        .into_par_iter()
        .map(|mut rec| {
            filter_record(&mut rec, &ctx, stats.as_deref())?;
            Ok(rec)
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(out)
}

fn filter_record(rec: &mut bcf::Record, ctx: &FilterCtx, stats: Option<&Stats>) -> Result<()> {
    if let Some(stats) = stats {
        let sample_count = rec.sample_count() as usize;
        for sample_i in 0..sample_count {
            stats.seen(sample_i);
        }
    }

    // Process records that are unfiltered (FILTER='.') or explicitly PASS.
    // In VCF, '.' means "no filters have been applied" (i.e., effectively PASS).
    let has_any_filter = {
        let mut filters = rec.filters();
        filters.next().is_some()
    };
    if has_any_filter && !rec.has_filter("PASS".as_bytes()) {
        return Ok(());
    }

    match &ctx.mode {
        FilterMode::Thresholds(thr) => filter_thresholds(rec, thr, ctx.max_ploidy, stats),
        FilterMode::Gbdt(xgb) => filter_gbdt(rec, xgb, ctx.max_ploidy, stats),
    }
}

fn filter_thresholds(
    rec: &mut bcf::Record,
    thr: &ThresholdCfg,
    max_ploidy: usize,
    stats: Option<&Stats>,
) -> Result<()> {
    let sample_count = rec.sample_count() as usize;

    // Extract FORMAT/GQ and FORMAT/DP arrays only if requested.
    let gq_opt = thr.min_gq.map(|_| rec.format(&thr.gq_tag).integer());
    let dp_opt = thr.min_dp.map(|_| rec.format(&thr.dp_tag).integer());

    let gq = match gq_opt {
        Some(Ok(buf)) => Some(buf),
        Some(Err(e)) => return Err(anyhow!("Missing/invalid GQ tag: {e}")),
        None => None,
    };
    let dp = match dp_opt {
        Some(Ok(buf)) => Some(buf),
        Some(Err(e)) => return Err(anyhow!("Missing/invalid DP tag: {e}")),
        None => None,
    };

    let mut fail_mask = vec![false; sample_count];

    for i in 0..sample_count {
        let mut failing = false;

        if let (Some(min_gq), Some(gq_buf)) = (thr.min_gq, &gq) {
            let v = gq_buf[i][0];
            if v < 0 || v < min_gq {
                failing = true;
            }
        }
        if let (Some(min_dp), Some(dp_buf)) = (thr.min_dp, &dp) {
            let v = dp_buf[i][0];
            if v < 0 || v < min_dp {
                failing = true;
            }
        }

        // Missing values in FORMAT integer fields should be treated as failing per requirements.
        if (thr.min_gq.is_some() && gq.is_none()) || (thr.min_dp.is_some() && dp.is_none()) {
            failing = true;
        }

        fail_mask[i] = failing;
    }

    mask_gt_by_sample(rec, &fail_mask, max_ploidy, stats)
}

#[derive(Copy, Clone, Debug)]
enum VariantClass {
    Snv,
    Insertion,
    Deletion,
}

fn variant_class_from_lengths(ref_len: usize, alt_len: usize) -> VariantClass {
    if ref_len == 1 && alt_len == 1 {
        VariantClass::Snv
    } else if alt_len > ref_len {
        VariantClass::Insertion
    } else {
        VariantClass::Deletion
    }
}

fn symbol_is_non_ref(allele: &[u8]) -> bool {
    // Minimal symbolic checks needed for is_snv definition.
    matches!(
        std::str::from_utf8(allele).unwrap_or(""),
        "<NON_REF>" | "<*>"
    )
}

fn compute_is_snv_site(alleles: &[&[u8]]) -> f32 {
    // Notebook definition:
    // is_snv = (len(ref) == 1 AND all non-symbolic alt alleles have len == 1)
    if alleles.is_empty() {
        return f32::NAN;
    }
    let ref_len = alleles[0].len();
    let mut all_alt_len1 = true;
    for alt_seq in alleles.iter().skip(1) {
        if symbol_is_non_ref(alt_seq) {
            continue;
        }
        if alt_seq.len() != 1 {
            all_alt_len1 = false;
            break;
        }
    }
    if ref_len == 1 && all_alt_len1 {
        1.0_f32
    } else {
        0.0_f32
    }
}

fn rnc_bad_from_raw(raw: &[u8]) -> f32 {
    // Notebook definition:
    // RNC_bad = any(c in rnc_str for c in ("I","U"))
    // After preprocessing, missing values are filled with -1.
    let s = std::str::from_utf8(raw).unwrap_or("");
    if s.is_empty() {
        f32::NAN
    } else if s.contains('I') || s.contains('U') {
        1.0_f32
    } else {
        0.0_f32
    }
}

fn filter_gbdt(
    rec: &mut bcf::Record,
    xgb: &GbdtCfg,
    max_ploidy: usize,
    stats: Option<&Stats>,
) -> Result<()> {
    let sample_count = rec.sample_count() as usize;

    // Site-level alleles for deriving is_snv and per-alt variant_class.
    let alleles = rec.alleles();
    if alleles.is_empty() {
        bail!("Record has no alleles");
    }
    let ref_len = alleles[0].len();
    let is_snv_site = compute_is_snv_site(&alleles);

    // INFO/AF and INFO/AQ (first element).
    let site_af: Option<f32> = rec
        .info(b"AF")
        .float()
        .ok()
        .and_then(|opt| opt.map(|v| v[0]));
    let site_aq: Option<f32> = rec
        .info(b"AQ")
        .float()
        .ok()
        .and_then(|opt| opt.map(|v| v[0]));

    // Extract per-sample GT and required FORMAT arrays.
    // NOTE: We compute per-alt feature rows in a per-record loop and then run one predict() per variant_class.
    let mut buf = bcf::record::Buffer::new();
    let gts = rec.genotypes_shared_buffer(&mut buf)?;

    let dp_tag = b"DP";
    let gq_tag = b"GQ";
    let ad_tag = b"AD";
    let pl_tag = b"PL";
    let rnc_tag = b"RNC";

    let dp = rec.format(dp_tag).integer().ok();
    let gq = rec.format(gq_tag).integer().ok();
    let ad = rec.format(ad_tag).integer().ok();
    let pl = rec.format(pl_tag).integer().ok();
    let rnc = rec.format(rnc_tag).string().ok();

    let mut fail_mask = vec![false; sample_count];

    // Build feature rows per variant_class.
    let mut rows_flat: [Vec<f32>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    let mut row_samples: [Vec<usize>; 3] = [Vec::new(), Vec::new(), Vec::new()];

    for sample_i in 0..sample_count {
        let gt = gts.get(sample_i);

        // Extract allele indices; missing alleles force masking because we cannot reliably compute is_het/is_hom_alt.
        let mut allele_indices: Vec<i32> = Vec::with_capacity(max_ploidy);
        for a in gt.iter() {
            let idx = a.index().map(|u| u as i32);
            if let Some(idx) = idx {
                allele_indices.push(idx);
            } else {
                allele_indices.push(-1);
            }
            if allele_indices.len() >= max_ploidy {
                break;
            }
        }

        if allele_indices.len() < 2 {
            // Non-diploid or missing: be conservative.
            fail_mask[sample_i] = true;
            continue;
        }

        let a1 = allele_indices[0];
        let a2 = allele_indices[1];
        if a1 < 0 || a2 < 0 {
            fail_mask[sample_i] = true;
            continue;
        }

        let is_het = if a1 != a2 { 1.0_f32 } else { 0.0_f32 };
        let is_hom_alt = if a1 == a2 && a1 > 0 { 1.0_f32 } else { 0.0_f32 };

        // For training: one row per alt allele present in the genotype.
        let mut alt_indices: HashSet<i32> = HashSet::new();
        for &ai in allele_indices.iter() {
            if ai > 0 {
                alt_indices.insert(ai);
            }
        }

        if alt_indices.is_empty() {
            // Homozygous reference genotype; no alt allele rows => never masked by model in this notebook logic.
            continue;
        }

        // Pull per-genotype FORMAT values once per sample.
        let dp_v = dp
            .as_ref()
            .map(|x| x[sample_i][0])
            .and_then(|v| if v < 0 { None } else { Some(v as f32) })
            .unwrap_or(f32::NAN);
        let gq_v = gq
            .as_ref()
            .map(|x| x[sample_i][0])
            .and_then(|v| if v < 0 { None } else { Some(v as f32) })
            .unwrap_or(f32::NAN);

        // PL: first 3 entries [0]=ref/ref, [1]=ref/het, [2]=hom-alt.
        let (pl_ref, pl_het, pl_hom) = match pl.as_ref() {
            Some(pl_buf) => {
                // PL is encoded as i32; we treat missing/unavailable as NaN.
                let row = pl_buf[sample_i];
                let get = |j: usize| {
                    row.get(j)
                        .copied()
                        .and_then(|v| if v < 0 { None } else { Some(v as f32) })
                        .unwrap_or(f32::NAN)
                };
                (get(0), get(1), get(2))
            }
            None => (f32::NAN, f32::NAN, f32::NAN),
        };

        // AD per allele: AD[0]=ref, AD[alt_idx]=the alt allele copy.
        let (ref_ad, total_ad) = match ad.as_ref() {
            Some(ad_buf) => {
                let row = ad_buf[sample_i];
                let ref_ad_v = row
                    .get(0)
                    .copied()
                    .and_then(|v| if v < 0 { None } else { Some(v as f32) })
                    .unwrap_or(f32::NAN);
                // If any AD entry is missing/invalid, treat totals as missing => AB becomes NaN.
                let mut total: f32 = 0.0;
                for &v in row.iter() {
                    if v < 0 {
                        total = f32::NAN;
                        break;
                    }
                    total += v as f32;
                }
                (ref_ad_v, total)
            }
            None => (f32::NAN, f32::NAN),
        };

        // RNC_bad: derived from FORMAT/RNC string, any of ("I","U") => 1 else 0.
        let rnc_bad = match rnc.as_ref() {
            Some(rnc_buf) => rnc_bad_from_raw(rnc_buf[sample_i]),
            None => f32::NAN,
        };

        for &alt_idx_i32 in alt_indices.iter() {
            let alt_idx = alt_idx_i32 as usize;
            let alt_seq = alleles.get(alt_idx).ok_or_else(|| anyhow!("alt_idx {alt_idx} out of range"))?;

            // alt_AD, AB:
            let alt_ad = match ad.as_ref() {
                Some(ad_buf) => ad_buf[sample_i]
                    .get(alt_idx)
                    .copied()
                    .and_then(|v| if v < 0 { None } else { Some(v as f32) })
                    .unwrap_or(f32::NAN),
                None => f32::NAN,
            };
            let ab = if !ad.is_none() {
                // AB = alt_AD / sum(AD); notebook says: if AD missing or totals are 0 => NaN.
                if !total_ad.is_finite() || total_ad == 0.0 {
                    f32::NAN
                } else {
                    alt_ad / total_ad
                }
            } else {
                f32::NAN
            };

            // Variant_class based on ref/alt lengths for this alt instance.
            let alt_len = alt_seq.len();
            let vc = variant_class_from_lengths(ref_len, alt_len);
            let class_idx = match vc {
                VariantClass::Snv => 0,
                VariantClass::Insertion => 1,
                VariantClass::Deletion => 2,
            };

            let site_af_v = site_af.filter(|&v| v >= 0.0).unwrap_or(f32::NAN);
            let site_aq_v = site_aq.filter(|&v| v >= 0.0).unwrap_or(f32::NAN);

            // Feature order must match the notebook training columns.
            // Columns (order matters): DP, GQ, AB, alt_AD, ref_AD, PL_ref, PL_het, PL_hom,
            // RNC_bad, site_AF, site_AQ, is_snv, is_het, is_hom_alt
            let mut row: [f32; 14] = [
                dp_v,
                gq_v,
                ab,
                alt_ad,
                ref_ad,
                pl_ref,
                pl_het,
                pl_hom,
                rnc_bad,
                site_af_v,
                site_aq_v,
                is_snv_site,
                is_het,
                is_hom_alt,
            ];

            // Notebook preprocessing: fill missing values with -1 before predict_proba.
            for v in &mut row {
                if !(*v).is_finite() {
                    *v = -1.0;
                }
            }

            rows_flat[class_idx].extend_from_slice(&row);
            row_samples[class_idx].push(sample_i);
        }
    }

    let model_threshold = xgb.model_threshold as f32;

    // Predict per variant_class and mark failing genotypes if ANY alt-instance score is below threshold.
    for class_idx in 0..3 {
        if row_samples[class_idx].is_empty() {
            continue;
        }
        let nrows = row_samples[class_idx].len();

        // Build one DataVec row per alt-instance feature vector.
        // `rows_flat[class_idx]` is a flat Vec<f32> with 14 columns per row.
        let mut data_vec: Vec<Data> = Vec::with_capacity(nrows);
        let flat = &rows_flat[class_idx];
        if flat.len() != nrows * 14 {
            bail!(
                "Unexpected flattened feature vector length for class {class_idx}: got {}, expected {}",
                flat.len(),
                nrows * 14
            );
        }
        for j in 0..nrows {
            let start = j * 14;
            let end = start + 14;
            data_vec.push(Data::new_test_data(flat[start..end].to_vec(), None));
        }

        let booster = match class_idx {
            0 => &xgb.snv,
            1 => &xgb.insertion,
            _ => &xgb.deletion,
        };
        let preds = booster.predict(&data_vec);
        if preds.len() != nrows {
            bail!("Unexpected predict size: got {}, expected {}", preds.len(), nrows);
        }
        for (j, &score) in preds.iter().enumerate() {
            if score < model_threshold {
                fail_mask[row_samples[class_idx][j]] = true;
            }
        }
    }

    mask_gt_by_sample(rec, &fail_mask, max_ploidy, stats)
}

fn genotype_is_fully_missing(gt: &bcf::record::Genotype, max_ploidy: usize) -> bool {
    let mut saw_any = false;
    for a in gt.iter().take(max_ploidy) {
        saw_any = true;
        if a.index().is_some() {
            return false;
        }
    }
    // If GT has fewer alleles than max_ploidy, we treat it as missing for "did we change it?" purposes.
    if !saw_any {
        return true;
    }
    true
}

fn mask_gt_by_sample(
    rec: &mut bcf::Record,
    fail_mask: &[bool],
    max_ploidy: usize,
    stats: Option<&Stats>,
) -> Result<()> {
    let sample_count = rec.sample_count() as usize;
    if fail_mask.len() != sample_count {
        bail!("fail_mask length mismatch");
    }

    let mut buf = bcf::record::Buffer::new();
    let gts = rec.genotypes_shared_buffer(&mut buf)?;

    // rust-htslib 0.44.x only exposes `push_genotypes` (flattened GT array).
    // We write exactly `max_ploidy` alleles per sample, padding/truncating as needed.
    let mut gt_flat: Vec<GenotypeAllele> = Vec::with_capacity(sample_count * max_ploidy);

    for sample_i in 0..sample_count {
        let gt = gts.get(sample_i);
        let was_fully_missing = genotype_is_fully_missing(&gt, max_ploidy);
        let mut alleles_out: Vec<GenotypeAllele> = Vec::with_capacity(max_ploidy);
        for a in gt.iter().take(max_ploidy) {
            match a.index() {
                Some(idx) => alleles_out.push(GenotypeAllele::Unphased(idx as i32)),
                None => alleles_out.push(GenotypeAllele::UnphasedMissing),
            }
        }

        while alleles_out.len() < max_ploidy {
            alleles_out.push(GenotypeAllele::UnphasedMissing);
        }

        if fail_mask[sample_i] {
            if let Some(stats) = stats {
                if !was_fully_missing {
                    stats.masked(sample_i);
                }
            }
            gt_flat.extend(std::iter::repeat(GenotypeAllele::UnphasedMissing).take(max_ploidy));
        } else {
            gt_flat.extend(alleles_out);
        }
    }

    rec.push_genotypes(&gt_flat)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use gbdt::config::Config;
    use rust_htslib::bcf::{Reader, Read};

    fn train_constant_gbdt(feature_size: usize, label: f32) -> GBDT {
        let mut cfg = Config::new();
        cfg.set_feature_size(feature_size);
        cfg.set_max_depth(2);
        cfg.set_min_leaf_size(1);
        cfg.set_iterations(3);
        cfg.set_shrinkage(0.1);
        cfg.set_loss("LogLikelyhood");

        let mut gbdt = GBDT::new(&cfg);
        let feature = vec![0.0; feature_size];
        let mut dv = Vec::new();
        for _ in 0..8 {
            dv.push(Data::new_training_data(feature.clone(), 1.0, label, None));
        }
        gbdt.fit(&mut dv);
        gbdt
    }

    #[test]
    fn variant_class_from_lengths_basic() {
        assert!(matches!(
            variant_class_from_lengths(1, 1),
            VariantClass::Snv
        ));
        assert!(matches!(
            variant_class_from_lengths(1, 2),
            VariantClass::Insertion
        ));
        assert!(matches!(
            variant_class_from_lengths(2, 1),
            VariantClass::Deletion
        ));
    }

    #[test]
    fn preprocess_missing_to_minus_one() {
        let mut row: [f32; 14] = [0.0; 14];
        row[3] = f32::NAN;
        row[4] = f32::INFINITY;
        for v in &mut row {
            if !(*v).is_finite() {
                *v = -1.0;
            }
        }
        assert_eq!(row[3], -1.0);
        assert_eq!(row[4], -1.0);
    }

    #[test]
    fn compute_is_snv_site_respects_symbolic() {
        let ref_seq = b"A";
        let alt_real = b"T";
        let alt_symbol = b"<NON_REF>";
        let alleles: Vec<&[u8]> = vec![ref_seq.as_ref(), alt_real.as_ref(), alt_symbol.as_ref()];
        assert_eq!(compute_is_snv_site(&alleles), 1.0);

        let alt_real_long = b"TT";
        let alleles2: Vec<&[u8]> = vec![ref_seq.as_ref(), alt_real_long.as_ref()];
        assert_eq!(compute_is_snv_site(&alleles2), 0.0);
    }

    #[test]
    fn test_rnc_bad_from_raw() {
        assert_eq!(rnc_bad_from_raw(b"I,foo"), 1.0);
        assert_eq!(rnc_bad_from_raw(b"U"), 1.0);
        assert_eq!(rnc_bad_from_raw(b"X"), 0.0);
        assert!(rnc_bad_from_raw(b"").is_nan());
    }

    #[test]
    fn gbdt_mode_masks_alt_sample_only() {
        // Two samples:
        // - s1: 0/0 (no alt allele) => should NOT be masked
        // - s2: 0/1 (alt allele present) => should be masked
        let vcf = r#"##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Genotype likelihoods">
##FORMAT=<ID=RNC,Number=1,Type=String,Description="RNC">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AQ,Number=A,Type=Float,Description="Allele Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	1	.	A	T	.	PASS	AF=0.1;AQ=0.9	GT:GQ:DP:AD:PL:RNC	0/0:50:30:30,0:0,100,100:A	0/1:10:5:1,10:100,0,100:I
"#;

        let tmp_dir = tempfile::tempdir().expect("tmpdir");
        let input_path = tmp_dir.path().join("in.vcf");
        std::fs::write(&input_path, vcf).expect("write vcf");

        let mut reader = Reader::from_path(&input_path).expect("reader");
        let mut record = reader.records().next().unwrap().expect("record");

        let snv = Arc::new(train_constant_gbdt(14, -1.0));
        let insertion = Arc::new(train_constant_gbdt(14, -1.0));
        let deletion = Arc::new(train_constant_gbdt(14, -1.0));

        let cfg = GbdtCfg {
            model_threshold: 0.99,
            snv,
            insertion,
            deletion,
        };

        filter_gbdt(&mut record, &cfg, 2, None).expect("filter_gbdt");

        let gts = record.genotypes().expect("genotypes");
        assert_eq!(gts.get(0).to_string(), "0/0");
        assert_eq!(gts.get(1).to_string(), "./.");
    }
}

