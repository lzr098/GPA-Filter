"""Shared fixtures for dgra-prefilter tests."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Minimal VCF header (GRCh38)
# ---------------------------------------------------------------------------
VCF_HEADER = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##assembly=GRCh38
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

# ---------------------------------------------------------------------------
# Helper: write a VCF from lines
# ---------------------------------------------------------------------------

def write_vcf(path: Path, header: str, records: list[str]) -> Path:
    """Write a minimal VCF file with the given header and record lines."""
    with open(path, "w") as f:
        f.write(header)
        for rec in records:
            f.write(rec.rstrip("\n") + "\n")
    return path


def write_bed(path: Path, intervals: list[tuple[str, int, int]]) -> Path:
    """Write a BED file from (chrom, start, end) tuples."""
    with open(path, "w") as f:
        for chrom, start, end in intervals:
            f.write(f"{chrom}\t{start}\t{end}\n")
    return path


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def tmp_dir(tmp_path: Path) -> Path:
    """Provide a temporary directory (alias for tmp_path)."""
    return tmp_path


@pytest.fixture()
def simple_bed(tmp_path: Path) -> Path:
    """A small BED file with 4 intervals on chr1."""
    bed_path = tmp_path / "test.bed"
    write_bed(bed_path, [
        ("chr1", 100, 200),
        ("chr1", 300, 400),
        ("chr2", 500, 600),
        ("chr2", 700, 800),
    ])
    return bed_path


@pytest.fixture()
def overlapping_bed(tmp_path: Path) -> Path:
    """A BED file with overlapping intervals on chr1."""
    bed_path = tmp_path / "overlapping.bed"
    write_bed(bed_path, [
        ("chr1", 100, 250),
        ("chr1", 200, 350),
        ("chr1", 340, 400),
    ])
    return bed_path


@pytest.fixture()
def tmp_vcf(tmp_path: Path) -> Path:
    """A small VCF with variants spanning gene/ncrna/regulatory/safetynet/intergenic regions.

    Uses GRCh38 header. Positions chosen to be testable against the
    mini_ref_dir BED files (see mini_ref_dir fixture).
    """
    records = [
        # Gene region hit (chr1:150, within BED chr1 100-200)
        "chr1\t150\t.\tA\tG\t30\tPASS\tDP=10",
        # ncRNA region hit (chr1:550, within BED chr1 500-600)
        "chr1\t550\t.\tC\tT\t30\tPASS\tDP=12",
        # Regulatory region hit (chr1:1100, within BED chr1 1000-1200)
        "chr1\t1100\t.\tG\tA\t30\tPASS\tDP=8",
        # ClinVar safetynet hit (chr1:2050, within BED chr1 2000-2100)
        "chr1\t2050\t.\tT\tC\t30\tPASS\tDP=15",
        # Intergenic (chr1:9999, not in any BED)
        "chr1\t9999\t.\tA\tT\t30\tPASS\tDP=5",
        # Gene + ClinVar overlap (chr1:108, within gene 100-200 AND clinvar 100-150)
        "chr1\t108\t.\tG\tC\t30\tPASS\tDP=20",
        # Another intergenic on chr2
        "chr2\t50\t.\tA\tG\t30\tPASS\tDP=3",
        # Gene on chr2 (within chr2 500-600)
        "chr2\t550\t.\tC\tA\t30\tPASS\tDP=11",
    ]
    vcf_path = tmp_path / "test.vcf"
    return write_vcf(vcf_path, VCF_HEADER, records)


@pytest.fixture()
def empty_vcf(tmp_path: Path) -> Path:
    """A VCF with header only, no variant records."""
    vcf_path = tmp_path / "empty.vcf"
    return write_vcf(vcf_path, VCF_HEADER, [])


def write_bed_named(
    path: Path,
    intervals: list[tuple[str, int, int, str]],
) -> Path:
    """Write a 4-column BED file from (chrom, start, end, name) tuples."""
    with open(path, "w") as f:
        for chrom, start, end, name in intervals:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\n")
    return path


@pytest.fixture()
def mini_ref_dir(tmp_path: Path) -> Path:
    """Create a minimal reference directory with all BED files required by presets.

    BED file layout (0-based half-open):
    - gencode_v44_gene_loci.bed: chr1 100-200, chr2 500-600
    - gencode_v44_coding_exon_utr.bed: chr1 100-180, chr2 500-580
    - gencode_5utr.bed: chr1 100-120
    - gencode_cds.bed: chr1 120-160
    - gencode_3utr.bed: chr1 160-180
    - gencode_splice_sites.bed: chr1 180-185
    - gencode_v44_ncrna_loci.bed: chr1 500-600
    - encode_screen_v3_ccres.bed: chr1 1000-1200
    - encode_screen_v3_pls_pels.bed: chr1 1000-1150
    - encode_screen_v3_balanced.bed: chr1 1000-1200 with types
    - fantom5_enhancers_promoters.bed: chr1 1050-1100
    - vista_enhancers.bed: chr1 1100-1150
    - clinvar_pathogenic_GRCh38.bed: chr1 100-150 (3star), chr1 2000-2100 (1star)
    - omim_pathogenic_GRCh38.bed: (empty)
    """
    ref_dir = tmp_path / "refs"
    ref_dir.mkdir()

    write_bed(ref_dir / "gencode_v44_gene_loci.bed", [
        ("chr1", 100, 200),
        ("chr2", 500, 600),
    ])
    write_bed(ref_dir / "gencode_v44_coding_exon_utr.bed", [
        ("chr1", 100, 180),
        ("chr2", 500, 580),
    ])
    write_bed(ref_dir / "gencode_5utr.bed", [
        ("chr1", 100, 120),
    ])
    write_bed(ref_dir / "gencode_cds.bed", [
        ("chr1", 120, 160),
    ])
    write_bed(ref_dir / "gencode_3utr.bed", [
        ("chr1", 160, 180),
    ])
    write_bed(ref_dir / "gencode_splice_sites.bed", [
        ("chr1", 180, 185),
    ])
    write_bed(ref_dir / "gencode_v44_ncrna_loci.bed", [
        ("chr1", 500, 600),
    ])
    write_bed(ref_dir / "encode_screen_v3_ccres.bed", [
        ("chr1", 1000, 1200),
    ])
    write_bed(ref_dir / "encode_screen_v3_pls_pels.bed", [
        ("chr1", 1000, 1150),
    ])
    write_bed_named(ref_dir / "encode_screen_v3_balanced.bed", [
        ("chr1", 1000, 1050, "PLS"),
        ("chr1", 1050, 1100, "pELS"),
        ("chr1", 1100, 1150, "dELS"),
        ("chr1", 1150, 1200, "CTCF"),
    ])
    write_bed_named(ref_dir / "ensembl_regulatory_features.bed", [
        ("chr1", 1200, 1300, "promoter"),
        ("chr1", 1300, 1400, "enhancer"),
        ("chr1", 1400, 1450, "ctcf"),
        ("chr1", 1450, 1500, "open_chromatin"),
        ("chr1", 1500, 1550, "tf_binding"),
    ])
    write_bed(ref_dir / "fantom5_enhancers_promoters.bed", [
        ("chr1", 1050, 1100),
    ])
    write_bed(ref_dir / "vista_enhancers.bed", [
        ("chr1", 1100, 1150),
    ])
    write_bed_named(ref_dir / "clinvar_pathogenic_GRCh38.bed", [
        ("chr1", 100, 150, "3"),
        ("chr1", 2000, 2100, "1"),
    ])
    # omim - empty file (just a header comment)
    (ref_dir / "omim_pathogenic_GRCh38.bed").write_text("# placeholder\n")

    return ref_dir


@pytest.fixture()
def clinvar_bed(tmp_path: Path) -> Path:
    """A ClinVar P/LP BED file with specific test positions."""
    bed_path = tmp_path / "clinvar_pathogenic_GRCh38.bed"
    write_bed(bed_path, [
        ("chr1", 100, 150),
        ("chr1", 2000, 2100),
        ("chr3", 300, 400),
    ])
    return bed_path
