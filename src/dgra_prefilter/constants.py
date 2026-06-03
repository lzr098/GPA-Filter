"""Shared constants for dgra-prefilter."""

from __future__ import annotations

from enum import IntEnum

# Standard chromosomes (GRCh38)
CHROMOSOMES: list[str] = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
    "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
    "chrX", "chrY", "chrM",
]

# BED file name constants
GENCODE_GENE_BED = "gencode_v44_gene_loci.bed"
GENCODE_NCRNA_BED = "gencode_v44_ncrna_loci.bed"
GENCODE_CODING_EXON_UTR_BED = "gencode_v44_coding_exon_utr.bed"
ENCODE_CCRE_BED = "encode_screen_v3_ccres.bed"
ENCODE_PLS_PELS_BED = "encode_screen_v3_pls_pels.bed"
FANTOM5_BED = "fantom5_enhancers_promoters.bed"
VISTA_BED = "vista_enhancers.bed"
CLINVAR_BED = "clinvar_pathogenic_GRCh38.bed"
OMIM_BED = "omim_pathogenic_GRCh38.bed"
MERGED_RETAINED_BED = "merged_retained_regions.bed"

# All BED filenames that should exist in refs/
ALL_BED_FILES: list[str] = [
    GENCODE_GENE_BED,
    GENCODE_NCRNA_BED,
    GENCODE_CODING_EXON_UTR_BED,
    ENCODE_CCRE_BED,
    ENCODE_PLS_PELS_BED,
    FANTOM5_BED,
    VISTA_BED,
    CLINVAR_BED,
    OMIM_BED,
    MERGED_RETAINED_BED,
]

# Default reference data directory
DEFAULT_REF_DIR = "~/.dgra-prefilter/refs"

# Log format template
LOG_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"

# GitHub release URL for reference data
GITHUB_RELEASE_URL = (
    "https://github.com/lzr098/GPA-Filter/releases/latest/download/refs.tar.gz"
)


class ErrorCode(IntEnum):
    """Error codes for dgra-prefilter."""

    # Input errors (1xxx)
    INPUT_FILE_NOT_FOUND = 1001
    INPUT_FORMAT_ERROR = 1002
    GENOME_MISMATCH = 1003
    EMPTY_VCF = 1004

    # Reference data errors (2xxx)
    REF_DIR_NOT_FOUND = 2001
    REF_FILE_MISSING = 2002
    REF_FILE_CORRUPT = 2003
    REF_UPDATE_FAILED = 2004

    # System dependency errors (3xxx)
    BCFTOOLS_NOT_FOUND = 3001
    BCFTOOLS_VERSION_LOW = 3002

    # Runtime errors (4xxx)
    BCFTOOLS_EXECUTION_ERROR = 4001
    TEMP_FILE_ERROR = 4002
    ANNOTATION_ERROR = 4003
