#!/usr/bin/env python3
"""
Post-filter script for dgra-prefilter output.

Applies two additional hard filters WITHOUT modifying dgra-prefilter code:
1. QC hard filter (GATK metrics)
2. Deep intron removal — keep only coding+UTR and splice-site ±100bp
   within gene loci; discard deep intronic variants.

Pure coordinate-based filtering (no DGRA tag dependency). This allows the
pipeline order: prefilter → postfilter → annotate, so annotation only
processes the already-trimmed variant set (~1.4M instead of ~3M).
"""
from __future__ import annotations

import gzip
import json
import logging
import os
import time
from collections import Counter, defaultdict
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# QC thresholds (moderate / relaxed)
# ---------------------------------------------------------------------------
QC_THRESHOLDS = {
    "QD": {"min": 1.5},
    "FS": {"max": 70.0},
    "SOR": {"max": 4.0},
    "MQ": {"min": 35.0},
    "ReadPosRankSum": {"min": -7.0},
    "MQRankSum": {"min": -13.0},
    "BaseQRankSum": {"min": -13.0},
}


def parse_info(info_str: str) -> dict[str, str]:
    d: dict[str, str] = {}
    for part in info_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d


def fails_qc(info_dict: dict[str, str]) -> bool:
    """Return True if variant fails QC thresholds."""
    for key, limits in QC_THRESHOLDS.items():
        val_str = info_dict.get(key)
        if val_str is None:
            continue
        try:
            val = float(val_str)
        except ValueError:
            continue
        if "min" in limits and val < limits["min"]:
            return True
        if "max" in limits and val > limits["max"]:
            return True
    return False


# ---------------------------------------------------------------------------
# Interval utilities
# ---------------------------------------------------------------------------
def load_bed(path: str | Path) -> dict[str, list[tuple[int, int]]]:
    """Load BED intervals by chromosome. Returns 0-based [start, end)."""
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            intervals[chrom].append((start, end))
    for chrom in intervals:
        intervals[chrom].sort()
    return dict(intervals)


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge sorted overlapping intervals."""
    if not intervals:
        return []
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(m[0], m[1]) for m in merged]


def overlaps(pos: int, intervals: list[tuple[int, int]]) -> bool:
    """Binary search for overlap with sorted intervals."""
    lo, hi = 0, len(intervals)
    while lo < hi:
        mid = (lo + hi) // 2
        s, e = intervals[mid]
        if pos < s:
            hi = mid
        elif pos >= e:
            lo = mid + 1
        else:
            return True
    return False


def build_splice_regions(
    coding_utr: dict[str, list[tuple[int, int]]], flank: int = 100
) -> dict[str, list[tuple[int, int]]]:
    """Build splice-site regions: ±flank bp around each exon boundary."""
    splice: dict[str, list[tuple[int, int]]] = {}
    for chrom, intervals in coding_utr.items():
        raw: list[tuple[int, int]] = []
        for start, end in intervals:
            raw.append((max(0, start - flank), start + flank))
            raw.append((max(0, end - flank), end + flank))
        raw.sort()
        splice[chrom] = merge_intervals(raw)
    return splice


# ---------------------------------------------------------------------------
# Main post-filter (coordinate-only, no DGRA dependency)
# ---------------------------------------------------------------------------
def postfilter(
    input_vcf: str | Path,
    output_vcf: str | Path,
    ref_dir: str | Path,
    preset: str = "comprehensive",
    qc_filter: bool = True,
    remove_deep_intron: bool = True,
    splice_flank: int = 100,
) -> dict:
    """Apply post-filtering to dgra-prefilter output using pure coordinates.

    Args:
        input_vcf: Path to filtered VCF from dgra-prefilter.
        output_vcf: Output path.
        ref_dir: Directory containing reference BED files.
        preset: Preset name determining which regulatory BED to use.
        qc_filter: Apply QC hard filter.
        remove_deep_intron: Remove deep intronic variants (keep only
            coding+UTR and splice-site ±flank bp).
        splice_flank: Flank size around exon boundaries (bp).

    Returns:
        Dictionary with filtering statistics.
    """
    start_time = time.time()
    input_vcf = Path(input_vcf)
    output_vcf = Path(output_vcf)
    ref_dir = Path(ref_dir).expanduser()

    if not input_vcf.exists():
        raise FileNotFoundError(f"Input VCF not found: {input_vcf}")

    # -----------------------------------------------------------------------
    # Load all retention-region BEDs (coordinate-only, no DGRA tags needed)
    # -----------------------------------------------------------------------
    logger.info("Loading reference BEDs from %s", ref_dir)
    coding_utr = load_bed(ref_dir / "gencode_v44_coding_exon_utr.bed")
    splice_regions = build_splice_regions(coding_utr, flank=splice_flank)
    ncrna = load_bed(ref_dir / "gencode_v44_ncrna_loci.bed")
    clinvar = load_bed(ref_dir / "clinvar_pathogenic_GRCh38.bed")

    if preset == "regulatory-minimal":
        regulatory = load_bed(ref_dir / "encode_screen_v3_pls_pels.bed")
    else:  # comprehensive or coding-only
        regulatory = load_bed(ref_dir / "encode_screen_v3_ccres.bed")

    # Count stats
    stats = {
        "input": 0,
        "qc_fail": 0,
        "deep_intron_removed": 0,
        "retained": 0,
        "by_region": Counter(),
        "by_reason": Counter(),
    }

    opener = gzip.open if str(input_vcf).endswith(".gz") else open
    out_opener = gzip.open if str(output_vcf).endswith(".gz") else open

    with opener(input_vcf, "rt") as fin, out_opener(output_vcf, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            stats["input"] += 1
            parts = line.strip().split("\t")
            chrom = parts[0]
            pos = int(parts[1]) - 1  # 0-based for BED comparison
            info_dict = parse_info(parts[7])

            # Step 1: QC filter
            if qc_filter and fails_qc(info_dict):
                stats["qc_fail"] += 1
                stats["by_reason"]["qc_fail"] += 1
                continue

            # Step 2: Coordinate-based retention
            # A variant is retained if it falls in ANY of the following:
            #   - coding exon + UTR
            #   - splice-site ±flank
            #   - ncRNA loci
            #   - regulatory elements
            #   - ClinVar safety net
            # Otherwise it's a deep intron (or other non-interesting region) → discard
            retained = False
            region_tags: list[str] = []

            if overlaps(pos, coding_utr.get(chrom, [])):
                retained = True
                region_tags.append("coding")
            if overlaps(pos, splice_regions.get(chrom, [])):
                retained = True
                region_tags.append("splice")
            if overlaps(pos, ncrna.get(chrom, [])):
                retained = True
                region_tags.append("ncrna")
            if overlaps(pos, regulatory.get(chrom, [])):
                retained = True
                region_tags.append("regulatory")
            if overlaps(pos, clinvar.get(chrom, [])):
                retained = True
                region_tags.append("clinvar")

            if remove_deep_intron and not retained:
                stats["deep_intron_removed"] += 1
                stats["by_reason"]["deep_intron"] += 1
                continue

            stats["retained"] += 1
            region_key = ",".join(sorted(region_tags)) if region_tags else "other"
            stats["by_region"][region_key] += 1
            fout.write(line)

    elapsed = time.time() - start_time

    report = {
        "input_variants": stats["input"],
        "qc_filtered": stats["qc_fail"],
        "deep_intron_removed": stats["deep_intron_removed"],
        "retained_variants": stats["retained"],
        "retention_rate": stats["retained"] / stats["input"] if stats["input"] > 0 else 0,
        "elapsed_seconds": round(elapsed, 1),
        "qc_thresholds": QC_THRESHOLDS,
        "splice_flank_bp": splice_flank,
        "by_region": dict(stats["by_region"]),
        "by_reason": dict(stats["by_reason"]),
    }

    report_path = output_vcf.with_suffix(output_vcf.suffix + ".postfilter.json")
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)

    logger.info("Post-filter complete: %d → %d (%.1f%%) in %.1fs",
                stats["input"], stats["retained"],
                report["retention_rate"] * 100, elapsed)
    logger.info("Report: %s", report_path)

    return report


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Post-filter dgra-prefilter output (coordinate-only, no DGRA dependency)"
    )
    parser.add_argument("-i", "--input", required=True, help="Input filtered VCF")
    parser.add_argument("-o", "--output", required=True, help="Output VCF")
    parser.add_argument("--ref-dir", default="~/.dgra-prefilter/refs",
                        help="Reference BED directory")
    parser.add_argument("--preset", default="comprehensive",
                        choices=["comprehensive", "coding-only", "regulatory-minimal"],
                        help="Preset determining regulatory BED to use")
    parser.add_argument("--no-qc", action="store_true", help="Skip QC filter")
    parser.add_argument("--keep-deep-intron", action="store_true",
                        help="Keep deep intronic variants")
    parser.add_argument("--splice-flank", type=int, default=100,
                        help="Splice-site flank size (bp)")

    args = parser.parse_args()

    report = postfilter(
        input_vcf=args.input,
        output_vcf=args.output,
        ref_dir=args.ref_dir,
        preset=args.preset,
        qc_filter=not args.no_qc,
        remove_deep_intron=not args.keep_deep_intron,
        splice_flank=args.splice_flank,
    )

    print(f"\n{'='*60}")
    print(f"Post-filter report")
    print(f"{'='*60}")
    print(f"Input variants:     {report['input_variants']:,}")
    print(f"QC filtered:        {report['qc_filtered']:,}")
    print(f"Deep intron removed:{report['deep_intron_removed']:,}")
    print(f"Retained:           {report['retained_variants']:,} ({report['retention_rate']:.1%})")
    print(f"Elapsed:            {report['elapsed_seconds']:.1f}s")
    print(f"\nRetained by region:")
    for region, count in sorted(report["by_region"].items(), key=lambda x: -x[1]):
        print(f"  {region}: {count:,}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
