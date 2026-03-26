#!/usr/bin/env python3
"""
Per-stratum precision/recall/F1 from RTG vcfeval output VCFs (SNV / insertion / deletion).

Uses baseline representation for recall (tp-baseline + fn) and call representation for
precision (tp + fp), matching RTG's overall definitions per stratum.

Same-length multi-base REF/ALT (MNPs) are omitted from these three buckets (rare in matched tp/fp/fn).
"""
from __future__ import annotations

import gzip
import sys
from collections import defaultdict
from pathlib import Path


def first_alt(alt_field: str) -> str:
    return alt_field.split(",")[0].strip()


def variant_class(ref: str, alt1: str) -> str | None:
    if not ref or not alt1:
        return None
    if alt1.startswith("<") or "[" in alt1 or "]" in alt1:
        return None
    if len(ref) == 1 and len(alt1) == 1:
        return "snv"
    if len(alt1) > len(ref):
        return "ins"
    if len(alt1) < len(ref):
        return "del"
    if len(ref) == len(alt1) and len(ref) > 1:
        return "mnp"
    return None


def count_by_class(vcf_gz: Path) -> dict[str, int]:
    out: dict[str, int] = defaultdict(int)
    with gzip.open(vcf_gz, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            ref, alt_field = parts[3], parts[4]
            vc = variant_class(ref, first_alt(alt_field))
            if vc is None:
                continue
            out[vc] += 1
    return dict(out)


def f1(p: float, r: float) -> float:
    if p + r <= 0:
        return 0.0
    return 2.0 * p * r / (p + r)


def pr_for_stratum(
    tp_b: dict[str, int],
    fn: dict[str, int],
    tp_c: dict[str, int],
    fp: dict[str, int],
    key: str,
) -> tuple[str, str, str]:
    tpb = tp_b.get(key, 0)
    fnn = fn.get(key, 0)
    tpc = tp_c.get(key, 0)
    fpp = fp.get(key, 0)
    r_d = tpb + fnn
    p_d = tpc + fpp
    recall = tpb / r_d if r_d > 0 else float("nan")
    prec = tpc / p_d if p_d > 0 else float("nan")
    f1v = f1(prec, recall) if prec == prec and recall == recall else float("nan")

    def fmt(x: float) -> str:
        return "NA" if x != x else f"{x:.6g}"

    return fmt(prec), fmt(recall), fmt(f1v)


def main() -> None:
    if len(sys.argv) != 2:
        print("usage: vcfeval_stratified_pr.py <vcfeval_output_dir>", file=sys.stderr)
        sys.exit(2)
    d = Path(sys.argv[1])
    req = ["tp-baseline.vcf.gz", "tp.vcf.gz", "fn.vcf.gz", "fp.vcf.gz"]
    for name in req:
        if not (d / name).is_file():
            print("", file=sys.stdout)
            sys.exit(1)

    tp_b = count_by_class(d / "tp-baseline.vcf.gz")
    tp_c = count_by_class(d / "tp.vcf.gz")
    fn = count_by_class(d / "fn.vcf.gz")
    fp = count_by_class(d / "fp.vcf.gz")

    cols: list[str] = []
    for key in ("snv", "ins", "del"):
        cols.extend(pr_for_stratum(tp_b, fn, tp_c, fp, key))

    print("\t".join(cols))


if __name__ == "__main__":
    main()
