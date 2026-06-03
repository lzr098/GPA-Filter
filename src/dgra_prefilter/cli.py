"""CLI entry point for dgra-prefilter."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from dgra_prefilter.constants import DEFAULT_REF_DIR, LOG_FORMAT


def _build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the dgra-prefilter CLI.

    Returns:
        Configured ArgumentParser instance.
    """
    parser = argparse.ArgumentParser(
        prog="dgra-prefilter",
        description=(
            "Genomic region prefilter for whole-genome VCF files. "
            "Retains variants in gene, ncRNA, and regulatory regions, "
            "plus ClinVar pathogenic sites as a safety net."
        ),
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        type=Path,
        help="Input VCF/VCF.gz/BCF file path",
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        type=Path,
        help="Output VCF path (.vcf.gz for compressed output)",
    )
    parser.add_argument(
        "-g", "--genome",
        default="GRCh38",
        choices=["GRCh38"],
        help="Genome version (default: GRCh38)",
    )
    parser.add_argument(
        "-p", "--preset",
        default="comprehensive",
        choices=["comprehensive", "coding-only", "regulatory-minimal"],
        help="Filter preset (default: comprehensive)",
    )
    parser.add_argument(
        "--ref-dir",
        default=DEFAULT_REF_DIR,
        type=Path,
        help=f"Reference BED files directory (default: {DEFAULT_REF_DIR})",
    )
    parser.add_argument(
        "--report",
        default=None,
        type=Path,
        help="JSON report output path (default: same dir as output)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Skip genome version validation",
    )
    parser.add_argument(
        "--update-refs",
        action="store_true",
        default=False,
        help="Update reference data before filtering",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="Enable verbose (DEBUG) logging",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="dgra-prefilter 1.0.0",
    )

    return parser


def _args_to_config(args: argparse.Namespace) -> dict:
    """Convert parsed CLI arguments to a dictionary for PrefilterConfig.

    Args:
        args: Parsed command-line arguments.

    Returns:
        Dictionary of PrefilterConfig constructor arguments.
    """
    return {
        "input_path": args.input,
        "output_path": args.output,
        "genome": args.genome,
        "preset": args.preset,
        "ref_dir": args.ref_dir,
        "report_path": args.report,
        "force": args.force,
        "update_refs": args.update_refs,
    }


def main() -> None:
    """Main CLI entry point.

    Parses command-line arguments, configures logging, and delegates
    to the prefilter_vcf() function.
    """
    parser = _build_parser()
    args = parser.parse_args()

    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format=LOG_FORMAT,
        stream=sys.stderr,
    )
    logger = logging.getLogger("dgra_prefilter")
    logger.setLevel(log_level)

    # Import here to avoid circular imports at module level
    from dgra_prefilter.core import (
        BcftoolsNotFoundError,
        GenomeMismatchError,
        PrefilterConfig,
        RefDataMissingError,
        VCFProcessingError,
        prefilter_vcf,
    )

    try:
        config_kwargs = _args_to_config(args)
        result = prefilter_vcf(**config_kwargs)

        # Print summary to stdout
        stats = result.stats
        print(
            f"Filtering complete: {stats.retained_variants:,}/{stats.input_variants:,} "
            f"variants retained ({stats.retention_rate:.1%})"
        )
        print(f"  Region-only: {stats.region_only_variants:,}")
        print(f"  Safety net-only: {stats.safetynet_only_variants:,}")
        print(f"  Both region and safety net: {stats.region_and_safetynet_variants:,}")
        print(f"  ClinVar safety net hits: {stats.clinvar_count:,}")
        print(f"  Elapsed: {stats.elapsed_seconds:.1f}s")
        print(f"  Output: {result.output_path}")
        print(f"  Report: {result.report_path}")

    except FileNotFoundError as exc:
        logger.error("File not found: %s", exc)
        sys.exit(1)
    except GenomeMismatchError as exc:
        logger.error("Genome mismatch: %s", exc)
        sys.exit(2)
    except BcftoolsNotFoundError as exc:
        logger.error("Dependency missing: %s", exc)
        sys.exit(3)
    except RefDataMissingError as exc:
        logger.error("Reference data missing: %s", exc)
        sys.exit(4)
    except VCFProcessingError as exc:
        logger.error("VCF processing error: %s", exc)
        sys.exit(5)
    except Exception as exc:
        logger.error("Unexpected error: %s", exc, exc_info=True)
        sys.exit(99)


if __name__ == "__main__":
    main()
