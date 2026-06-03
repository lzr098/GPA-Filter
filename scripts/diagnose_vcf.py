#!/usr/bin/env python3
"""Diagnose a VCF file and output a JSON quality assessment report."""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import sys
from collections import Counter
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# Chromosome lengths for GRCh38 (used for outside-region placement)
GRCH38_CHROM_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
    "chrM": 16569,
}

GRCH37_CHROM_LENGTHS = {
    "1": 249250621, "2": 243199373, "3": 198022430,
    "4": 191154276, "5": 180915260, "6": 171115067,
    "7": 159138663, "8": 146364022, "9": 141213431,
    "10": 135534747, "11": 135006516, "12": 133851895,
    "13": 115169878, "14": 107349540, "15": 102531392,
    "16": 90354753, "17": 81195210, "18": 78077248,
    "19": 59128983, "20": 63025520, "21": 48129895,
    "22": 51304566, "X": 155270560, "Y": 59373566,
    "MT": 16569,
}


def open_vcf(path: str | Path):
    """Open a VCF file (plain or gzipped) for reading text."""
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_gt(sample_col: str) -> str:
    """Extract GT value from a VCF sample column."""
    if not sample_col or sample_col == ".":
        return "./."
    parts = sample_col.split(":")
    gt = parts[0]
    return gt if gt else "./."


def infer_genome(
    has_chr_prefix: bool,
    chroms_seen: set[str],
    mitochondria_id: str | None,
) -> str:
    """Infer genome build from chromosome naming patterns."""
    grch38_chroms = set(GRCH38_CHROM_LENGTHS.keys())
    grch37_chroms = set(GRCH37_CHROM_LENGTHS.keys())

    if has_chr_prefix:
        # chr prefix strongly suggests GRCh38 (or hg19 with chr)
        if mitochondria_id == "chrM":
            return "GRCh38"
        if "chrM" in chroms_seen:
            return "GRCh38"
        # hg19 also uses chr prefix but usually chrMT
        if "chrMT" in chroms_seen:
            return "GRCh37"
        # Default with chr prefix
        return "GRCh38"
    else:
        # No chr prefix
        if mitochondria_id == "MT":
            return "GRCh37"
        if "MT" in chroms_seen:
            return "GRCh37"
        if mitochondria_id == "chrM":
            return "GRCh38"
        # Check overlap with known chromosomes
        normalized = chroms_seen & (grch37_chroms | grch38_chroms)
        if normalized.issubset(grch37_chroms) and not normalized & grch38_chroms:
            return "GRCh37"
        return "unknown"


def diagnose_vcf(
    input_path: str | Path,
    max_variants: int = 100_000,
) -> dict:
    """Diagnose a VCF file and return a JSON-compatible report dictionary.

    Args:
        input_path: Path to input VCF (.vcf or .vcf.gz).
        max_variants: Maximum number of variants to scan (default 100,000).

    Returns:
        Dictionary containing diagnosis report.
    """
    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input VCF not found: {input_path}")

    total_variants = 0
    filter_counts: Counter = Counter()
    gt_counts: Counter = Counter()
    ref_eq_alt_count = 0
    no_alt_count = 0

    chroms_seen: set[str] = set()
    has_chr_prefix = False
    mitochondria_id: str | None = None
    filter_uniformly_dot = True

    with open_vcf(input_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            if total_variants >= max_variants:
                break

            total_variants += 1
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info = parts[:8]
            chroms_seen.add(chrom)

            if chrom.startswith("chr"):
                has_chr_prefix = True

            if chrom in ("chrM", "chrMT", "M", "MT"):
                mitochondria_id = chrom

            filter_counts[filt] += 1
            if filt != ".":
                filter_uniformly_dot = False

            if ref == alt:
                ref_eq_alt_count += 1

            if not alt or alt == ".":
                no_alt_count += 1

            # Parse genotype if sample columns exist
            if len(parts) >= 10:
                sample_col = parts[9]
                gt = parse_gt(sample_col)
                gt_counts[gt] += 1

    inferred_genome = infer_genome(has_chr_prefix, chroms_seen, mitochondria_id)

    # Determine quality assessment
    has_ref_eq_alt = ref_eq_alt_count > 0
    has_uncalled = gt_counts.get("./.", 0) > 0 or gt_counts.get(".", 0) > 0
    has_hom_ref = gt_counts.get("0/0", 0) > 0 or gt_counts.get("0|0", 0) > 0
    has_low_qual = any(
        k != "PASS" and k != "." for k in filter_counts.keys()
    )

    if inferred_genome == "GRCh37":
        quality_assessment = "needs_liftover"
    elif has_ref_eq_alt or has_uncalled or has_hom_ref or has_low_qual:
        quality_assessment = "raw_caller_output"
    elif filter_uniformly_dot and not has_hom_ref and not has_uncalled:
        # No FILTER field, all variants are called
        quality_assessment = "clean_genotyped"
    elif not has_ref_eq_alt and not has_uncalled and not has_hom_ref:
        quality_assessment = "clean_genotyped"
    else:
        quality_assessment = "ambiguous"

    # Build recommendations
    recommendations: list[str] = []
    if quality_assessment == "needs_liftover":
        recommendations.append(
            "VCF appears to be GRCh37. Liftover to GRCh38 is required before processing."
        )
    elif quality_assessment == "raw_caller_output":
        recommendations.append(
            "VCF contains uncalled genotypes, homozygous reference calls, or non-PASS filters. "
            "Applying conservative preprocessing filters."
        )
    elif quality_assessment == "ambiguous":
        recommendations.append(
            "VCF quality is ambiguous. Applying conservative preprocessing filters."
        )
    else:
        recommendations.append("VCF is clean and genotyped. No preprocessing needed.")

    if has_ref_eq_alt:
        recommendations.append(
            f"Found {ref_eq_alt_count} variants where REF == ALT. These should be removed."
        )

    report = {
        "total_variants": total_variants,
        "filter_distribution": dict(filter_counts),
        "gt_distribution": dict(gt_counts),
        "ref_eq_alt_count": ref_eq_alt_count,
        "no_alt_count": no_alt_count,
        "has_chr_prefix": has_chr_prefix,
        "mitochondria_id": mitochondria_id,
        "inferred_genome": inferred_genome,
        "quality_assessment": quality_assessment,
        "recommendations": recommendations,
    }

    return report


def main() -> None:
    """CLI entry point for VCF diagnosis."""
    parser = argparse.ArgumentParser(
        description="Diagnose a VCF file and output a JSON quality assessment report."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input VCF file (.vcf or .vcf.gz)"
    )
    parser.add_argument(
        "-n",
        "--max-variants",
        type=int,
        default=100_000,
        help="Maximum number of variants to scan (default: 100000)",
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output JSON report path"
    )

    args = parser.parse_args()

    report = diagnose_vcf(args.input, max_variants=args.max_variants)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(report, f, indent=2)

    logger.info("Diagnosis complete. Report written to %s", output_path)
    logger.info("Quality assessment: %s", report["quality_assessment"])
    logger.info("Inferred genome: %s", report["inferred_genome"])


if __name__ == "__main__":
    main()
