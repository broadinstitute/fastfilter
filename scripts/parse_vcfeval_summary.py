#!/usr/bin/env python3
"""
Parse RTG Tools vcfeval summary.txt and print one TSV line of metrics.

Tries several layouts used across RTG versions. Prefers the "no threshold" / all-variants
block when present; otherwise takes the first Precision/Sensitivity/F-measure triple.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path


def _find_floats_after_labels(text: str, labels: tuple[str, ...]) -> dict[str, float | None]:
    out: dict[str, float | None] = {k: None for k in labels}
    for label in labels:
        # "Precision   0.99" or "    Precision: 0.99" or "Precision\t0.99"
        m = re.search(
            rf"(?im)^\s*{re.escape(label)}\s*[:\s]\s*([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\s*$",
            text,
            re.MULTILINE,
        )
        if m:
            out[label] = float(m.group(1))
            continue
        # Same line anywhere: " ... Precision 0.99 ..."
        m = re.search(
            rf"(?i)\b{re.escape(label)}\s*[:\s]\s*([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\b",
            text,
        )
        if m:
            out[label] = float(m.group(1))
    return out


def _slice_no_threshold(text: str) -> str:
    """Prefer text after 'No threshold' / 'no score' if present."""
    lower = text.lower()
    for needle in (
        "no score threshold",
        "no threshold",
        "without threshold",
    ):
        i = lower.find(needle)
        if i >= 0:
            return text[i : i + 8000]
    return text


def _parse_rtg_table(text: str) -> dict[str, float | None] | None:
    """RTG 3.x summary.txt: whitespace table with Threshold / Precision / Sensitivity / F-measure columns."""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip() and not ln.strip().startswith("-")]
    header_idx = None
    for i, ln in enumerate(lines):
        if ln.startswith("Threshold") and "Precision" in ln and "F-measure" in ln:
            header_idx = i
            break
    if header_idx is None:
        return None
    cols = re.split(r"\s+", lines[header_idx])
    try:
        ip = cols.index("Precision")
        is_ = cols.index("Sensitivity")
        im = cols.index("F-measure")
    except ValueError:
        return None
    # Prefer the "None" threshold row (no score cutoff on the ROC); else last data row.
    chosen: list[str] | None = None
    for ln in lines[header_idx + 1 :]:
        parts = re.split(r"\s+", ln)
        if len(parts) <= max(ip, is_, im):
            continue
        if parts[0].lower() == "none":
            chosen = parts
            break
    if chosen is None:
        for ln in reversed(lines[header_idx + 1 :]):
            parts = re.split(r"\s+", ln)
            if len(parts) > max(ip, is_, im) and parts[0] not in ("Threshold",):
                chosen = parts
                break
    if chosen is None:
        return None
    try:
        return {
            "precision": float(chosen[ip]),
            "recall": float(chosen[is_]),
            "f1": float(chosen[im]),
        }
    except (ValueError, IndexError):
        return None


def parse_summary(path: Path) -> dict[str, float | None]:
    text = path.read_text(errors="replace")
    tab = _parse_rtg_table(text)
    if tab is not None:
        return tab
    focus = _slice_no_threshold(text)
    labels = ("Precision", "Sensitivity", "F-measure")
    m = _find_floats_after_labels(focus, labels)
    if m["Precision"] is None:
        m = _find_floats_after_labels(text, labels)
    # Some docs use Recall instead of Sensitivity
    if m.get("Sensitivity") is None:
        alt = _find_floats_after_labels(focus, ("Precision", "Recall", "F-measure"))
        if alt["Recall"] is not None:
            m["Sensitivity"] = alt["Recall"]
            m["Precision"] = m["Precision"] or alt["Precision"]
            m["F-measure"] = m["F-measure"] or alt["F-measure"]
    return {
        "precision": m["Precision"],
        "recall": m["Sensitivity"],
        "f1": m["F-measure"],
    }


def main() -> None:
    if len(sys.argv) != 2:
        print("usage: parse_vcfeval_summary.py <summary.txt>", file=sys.stderr)
        sys.exit(2)
    p = Path(sys.argv[1])
    if not p.is_file():
        print("", end="")
        sys.exit(1)
    d = parse_summary(p)
    pr, rc, f1 = d["precision"], d["recall"], d["f1"]

    def fmt(x: float | None) -> str:
        return "" if x is None else f"{x:.9g}"

    print("\t".join((fmt(pr), fmt(rc), fmt(f1))))


if __name__ == "__main__":
    main()
