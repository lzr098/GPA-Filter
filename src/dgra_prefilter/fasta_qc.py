#!/usr/bin/env python3
"""
GRCh38 FASTA QC for dgra-prefilter
v0.1.0 - 2026-06-07

Optional REF allele verification using local reference genome.
Safe to import even if FASTA is unavailable (no-op).
"""

import logging
from pathlib import Path
from typing import Optional, Tuple

logger = logging.getLogger(__name__)

_FASTA_MODULE = Path("/Users/zhaorongli/.workbuddy/scripts/grch38_fasta_local.py")


def _import_fasta():
    if not _FASTA_MODULE.exists():
        return None
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location("grch38_fasta_local", _FASTA_MODULE)
        if spec is None or spec.loader is None:
            return None
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod
    except Exception:
        return None


def is_fasta_available() -> bool:
    """Check if local GRCh38 FASTA is available."""
    fasta = _import_fasta()
    if fasta is None:
        return False
    try:
        return fasta.FASTA_PATH.exists()
    except Exception:
        return False


def log_fasta_status() -> None:
    """Log whether local GRCh38 FASTA is available."""
    if is_fasta_available():
        logger.info("Local GRCh38 FASTA available for optional REF verification")
    else:
        logger.debug("Local GRCh38 FASTA not available (skipping REF verification)")


def verify_vcf_ref(chrom: str, pos: int, ref: str) -> Tuple[bool, Optional[str], Optional[str]]:
    """Verify a single variant's REF allele against GRCh38.

    Returns:
        (is_correct, actual_ref, error_msg)
    """
    fasta = _import_fasta()
    if fasta is None:
        return True, None, "FASTA module unavailable"
    try:
        ok, actual = fasta.verify_ref(chrom, pos, ref)
        if ok:
            return True, actual, None
        return False, actual, f"REF mismatch: VCF='{ref}', genome='{actual}'"
    except Exception as e:
        return True, None, str(e)
