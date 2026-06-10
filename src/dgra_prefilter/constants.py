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
GENCODE_5UTR_BED = "gencode_5utr.bed"
GENCODE_CDS_BED = "gencode_cds.bed"
GENCODE_3UTR_BED = "gencode_3utr.bed"
GENCODE_SPLICE_BED = "gencode_splice_sites.bed"
ENCODE_CCRE_BED = "encode_screen_v3_ccres.bed"
ENCODE_PLS_PELS_BED = "encode_screen_v3_pls_pels.bed"
ENCODE_BALANCED_BED = "encode_screen_v3_balanced.bed"
FANTOM5_BED = "fantom5_enhancers_promoters.bed"
VISTA_BED = "vista_enhancers.bed"
CLINVAR_BED = "clinvar_pathogenic_GRCh38.bed"
OMIM_BED = "omim_pathogenic_GRCh38.bed"
MERGED_RETAINED_BED = "merged_retained_regions.bed"
ENSEMBL_REGULATORY_BED = "ensembl_regulatory_features.bed"

# Region sub-tag constants
REGION_TAG_GENE_5UTR = "gene_5utr"
REGION_TAG_GENE_CDS = "gene_cds"
REGION_TAG_GENE_3UTR = "gene_3utr"
REGION_TAG_GENE_SPLICE = "gene_splice"
REGION_TAG_REGULATORY_PLS = "regulatory_pls"
REGION_TAG_REGULATORY_PELS = "regulatory_pels"
REGION_TAG_REGULATORY_DELS = "regulatory_dels"
REGION_TAG_REGULATORY_CTCF = "regulatory_ctcf"
REGION_TAG_REGULATORY_PROMOTER = "regulatory_promoter"
REGION_TAG_REGULATORY_ENHANCER = "regulatory_enhancer"
REGION_TAG_REGULATORY_OPEN_CHROMATIN = "regulatory_open_chromatin"
REGION_TAG_REGULATORY_TF_BINDING = "regulatory_tf_binding"

# All BED filenames that should exist in refs/
ALL_BED_FILES: list[str] = [
    GENCODE_GENE_BED,
    GENCODE_NCRNA_BED,
    GENCODE_CODING_EXON_UTR_BED,
    GENCODE_5UTR_BED,
    GENCODE_CDS_BED,
    GENCODE_3UTR_BED,
    GENCODE_SPLICE_BED,
    ENCODE_CCRE_BED,
    ENCODE_PLS_PELS_BED,
    ENCODE_BALANCED_BED,
    FANTOM5_BED,
    VISTA_BED,
    CLINVAR_BED,
    OMIM_BED,
    MERGED_RETAINED_BED,
    ENSEMBL_REGULATORY_BED,
]

# Default reference data directory
DEFAULT_REF_DIR = "~/.dgra-prefilter/refs"

# Log format template
LOG_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"

# GitHub release URL for reference data
GITHUB_RELEASE_URL = (
    "https://github.com/lzr098/GPA-Filter/releases/latest/download/refs.tar.gz"
)

# Ensembl Regulatory Build URL
ENSEMBL_REGULATORY_URL = (
    "ftp://ftp.ensembl.org/pub/release-115/regulation/homo_sapiens/GRCh38/"
    "annotation/Homo_sapiens.GRCh38.regulatory_features.v115.gff3.gz"
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
