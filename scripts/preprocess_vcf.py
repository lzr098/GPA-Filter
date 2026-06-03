#!/usr/bin/env python3
"""Preprocess a VCF file based on a diagnose report."""

from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Path to bcftools binary
BCFTOOLS = os.environ.get("BCFTOOLS", "/Users/zhaorongli/.micromamba/envs/bcftools/bin/bcftools")


def run_bcftools(args: list[str], stdin=None) -> subprocess.CompletedProcess:
    """Run bcftools with the given arguments."""
    cmd = [BCFTOOLS] + args
    logger.debug("Running: %s", " ".join(cmd))
    result = subprocess.run(
        cmd,
        stdin=stdin,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        logger.error("bcftools failed: %s", result.stderr)
        raise subprocess.CalledProcessError(
            result.returncode, cmd, output=result.stdout, stderr=result.stderr
        )
    return result


def count_variants(vcf_path: str | Path) -> int:
    """Count variants in a VCF file using bcftools."""
    try:
        result = run_bcftools(["index", "-n", str(vcf_path)])
        return int(result.stdout.strip())
    except Exception:
        # Fallback: count non-header lines
        import gzip

        path = Path(vcf_path)
        opener = gzip.open if str(path).endswith(".gz") else open
        count = 0
        with opener(path, "rt") as f:
            for line in f:
                if not line.startswith("#"):
                    count += 1
        return count


def preprocess_vcf(
    input_path: str | Path,
    diagnose_report_path: str | Path,
    output_path: str | Path,
) -> dict:
    """Preprocess a VCF based on its diagnose report.

    Args:
        input_path: Input VCF file (.vcf or .vcf.gz).
        diagnose_report_path: Path to JSON diagnose report.
        output_path: Output VCF file path.

    Returns:
        Dictionary with preprocessing results.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    diagnose_report_path = Path(diagnose_report_path)

    if not input_path.exists():
        raise FileNotFoundError(f"Input VCF not found: {input_path}")
    if not diagnose_report_path.exists():
        raise FileNotFoundError(f"Diagnose report not found: {diagnose_report_path}")

    with open(diagnose_report_path) as f:
        report = json.load(f)

    assessment = report.get("quality_assessment", "ambiguous")
    inferred_genome = report.get("inferred_genome", "unknown")
    filter_distribution = report.get("filter_distribution", {})
    has_non_dot_filter = any(k != "." for k in filter_distribution.keys())

    input_count = count_variants(input_path)
    result = {
        "input_variants": input_count,
        "output_variants": input_count,
        "action": assessment,
        "output_path": str(output_path),
    }

    if assessment == "needs_liftover":
        logger.error(
            "VCF appears to be %s. Liftover to GRCh38 is required before processing.",
            inferred_genome,
        )
        sys.exit(1)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    if assessment == "clean_genotyped":
        # Copy or symlink through
        logger.info("VCF is clean_genotyped. Copying through to output.")
        if str(input_path).endswith(".gz") and str(output_path).endswith(".gz"):
            shutil.copy2(input_path, output_path)
        else:
            # If compression mismatches, use bcftools view
            run_bcftools(
                ["view", "-Oz", "-o", str(output_path), str(input_path)]
            )
        result["output_variants"] = count_variants(output_path)
        return result

    if assessment == "ambiguous":
        logger.warning(
            "VCF quality is ambiguous. Applying conservative preprocessing filters."
        )

    # assessment in ("raw_caller_output", "ambiguous")
    logger.info("Applying preprocessing filters to VCF.")

    # Build bcftools filter pipeline
    # Step 1: Remove REF == ALT
    cmd1 = [BCFTOOLS, "view", "-e", "REF=ALT", str(input_path)]

    # Step 2: Remove uncalled and homozygous reference genotypes
    gt_expr = 'FMT/GT!="0/0" && FMT/GT!="0|0" && FMT/GT!="." && FMT/GT!="./." && FMT/GT!=".|."'

    # Step 3: Filter to PASS only if FILTER is not uniformly "."
    if has_non_dot_filter:
        # We need PASS filter too
        cmd = [
            BCFTOOLS,
            "view",
            "-e",
            "REF=ALT",
            str(input_path),
        ]
        # Pipe to second bcftools for GT and FILTER
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        cmd2 = [
            BCFTOOLS,
            "view",
            "-i",
            gt_expr,
            "-f",
            "PASS,.",
            "-Oz",
            "-o",
            str(output_path),
        ]
        result_proc = subprocess.run(cmd2, stdin=p1.stdout, capture_output=True, text=True)
        p1.wait()
        if p1.returncode != 0:
            logger.error("bcftools step 1 failed")
            raise subprocess.CalledProcessError(p1.returncode, cmd)
        if result_proc.returncode != 0:
            logger.error("bcftools step 2 failed: %s", result_proc.stderr)
            raise subprocess.CalledProcessError(
                result_proc.returncode, cmd2, output=result_proc.stdout, stderr=result_proc.stderr
            )
    else:
        # No PASS filter needed, just GT and REF!=ALT
        cmd = [
            BCFTOOLS,
            "view",
            "-e",
            "REF=ALT",
            str(input_path),
        ]
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        cmd2 = [
            BCFTOOLS,
            "view",
            "-i",
            gt_expr,
            "-Oz",
            "-o",
            str(output_path),
        ]
        result_proc = subprocess.run(cmd2, stdin=p1.stdout, capture_output=True, text=True)
        p1.wait()
        if p1.returncode != 0:
            logger.error("bcftools step 1 failed")
            raise subprocess.CalledProcessError(p1.returncode, cmd)
        if result_proc.returncode != 0:
            logger.error("bcftools step 2 failed: %s", result_proc.stderr)
            raise subprocess.CalledProcessError(
                result_proc.returncode, cmd2, output=result_proc.stdout, stderr=result_proc.stderr
            )

    result["output_variants"] = count_variants(output_path)
    removed = result["input_variants"] - result["output_variants"]
    logger.info(
        "Preprocessing complete: %d → %d (%d removed, %.1f%%)",
        result["input_variants"],
        result["output_variants"],
        removed,
        (removed / result["input_variants"] * 100) if result["input_variants"] > 0 else 0,
    )

    return result


def main() -> None:
    """CLI entry point for VCF preprocessing."""
    parser = argparse.ArgumentParser(
        description="Preprocess a VCF file based on a diagnose report."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input VCF file (.vcf or .vcf.gz)"
    )
    parser.add_argument(
        "--diagnose-report", required=True, help="Path to JSON diagnose report"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output VCF file path"
    )

    args = parser.parse_args()

    preprocess_vcf(args.input, args.diagnose_report, args.output)


if __name__ == "__main__":
    main()
