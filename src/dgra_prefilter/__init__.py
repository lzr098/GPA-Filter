"""dgra-prefilter: Genomic region prefilter for whole-genome VCF files.

This package provides tools for filtering whole-genome VCF files to retain
only variants in biologically relevant regions (genes, ncRNA, regulatory
elements) and known pathogenic sites from ClinVar.
"""

from __future__ import annotations

from dgra_prefilter.core import (
    BcftoolsNotFoundError,
    FilterResult,
    FilterStats,
    GenomeMismatchError,
    PrefilterConfig,
    RefDataMissingError,
    VCFProcessingError,
    annotate_vcf_file,
    prefilter_vcf,
)

__all__ = [
    "prefilter_vcf",
    "annotate_vcf_file",
    "FilterResult",
    "FilterStats",
    "PrefilterConfig",
    "GenomeMismatchError",
    "BcftoolsNotFoundError",
    "RefDataMissingError",
    "VCFProcessingError",
]

__version__ = "1.0.0"
