#!/usr/bin/env python3
"""Build reference BED files for dgra-prefilter.

Downloads raw data from GENCODE, ENCODE, FANTOM5, Vista, and ClinVar,
then converts each source into BED format for use by the prefilter.

Usage:
    python scripts/build_refs.py [--output-dir refs/] [--skip-download]
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

# Add src to path so we can import dgra_prefilter
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from dgra_prefilter.bed_utils import BedUtils
from dgra_prefilter.constants import (
    CLINVAR_BED,
    ENCODE_BALANCED_BED,
    ENCODE_CCRE_BED,
    ENCODE_PLS_PELS_BED,
    ENSEMBL_REGULATORY_BED,
    ENSEMBL_REGULATORY_URL,
    FANTOM5_BED,
    GENCODE_3UTR_BED,
    GENCODE_5UTR_BED,
    GENCODE_CDS_BED,
    GENCODE_CODING_EXON_UTR_BED,
    GENCODE_GENE_BED,
    GENCODE_NCRNA_BED,
    GENCODE_SPLICE_BED,
    MERGED_RETAINED_BED,
    OMIM_BED,
    VISTA_BED,
)

logger = logging.getLogger(__name__)

# Shared local data assets (preferred over download)
SHARED_DATA_DIR = Path.home() / ".workbuddy" / "data"
SHARED_GENCODE_GTF = SHARED_DATA_DIR / "gencode" / "gencode.v44.annotation.gtf.gz"

# Data source URLs
GENCODE_GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_44/gencode.v44.annotation.gtf.gz"
)
ENCODE_SCREEN_URL = (
    "https://www.encodeproject.org/files/ENCFF420VPZ/"
    "@@download/ENCFF420VPZ.bed.gz"
)
FANTOM5_ENHANCER_URL = (
    "https://dbarchive.biosciencedbc.jp/data/fantom5/datafiles/"
    "reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz"
)
FANTOM5_CAGE_URL = (
    "https://dbarchive.biosciencedbc.jp/data/fantom5/datafiles/"
    "reprocessed/hg38_latest/extra/CAGE_peaks/"
    "hg38_fair+new_CAGE_peaks_phase1and2.bed.gz"
)
VISTA_URL = (
    "https://hgdownload.soe.ucsc.edu/gbdb/hg38/vistaEnhancers/"
    "vistaEnhancers.bb"
)
CLINVAR_VCF_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/"
    "clinvar.vcf.gz"
)

# GENCODE biotype categories
PROTEIN_CODING_BIOTYPE = "protein_coding"
NCRNA_BIOTYPES = {
    "3prime_overlapping_ncRNA",
    "antisense",
    "bidirectional_promoter_lncRNA",
    "lincRNA",
    "lncRNA",
    "macrolncRNA",
    "miRNA",
    "misc_RNA",
    "non_coding",
    "polymorphic_pseudogene",
    "processed_transcript",
    "rRNA",
    "ribozyme",
    "sRNA",
    "scRNA",
    "scaRNA",
    "snRNA",
    "snoRNA",
    "transcribed_processed_pseudogene",
    "transcribed_unitary_pseudogene",
    "transcribed_unprocessed_pseudogene",
    "translated_processed_pseudogene",
    "translated_unprocessed_pseudogene",
    "unitary_pseudogene",
    "unprocessed_pseudogene",
    "vaultRNA",
}

# Exon/UTR feature types for coding-only preset
CODING_FEATURE_TYPES = {"exon", "five_prime_utr", "three_prime_utr", "UTR", "CDS"}


def _download_file(url: str, output_path: Path) -> Path:
    """Download a file using curl or wget.

    Args:
        url: URL to download from.
        output_path: Local path to save the file.

    Returns:
        Path to the downloaded file.

    Raises:
        RuntimeError: If download fails.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Cache check: skip download if file exists and non-empty
    if output_path.exists() and output_path.stat().st_size > 0:
        logger.info("Using cached file: %s (%d bytes)", output_path, output_path.stat().st_size)
        return output_path

    logger.info("Downloading %s to %s", url, output_path)

    for cmd_name in ["curl", "wget"]:
        try:
            if cmd_name == "curl":
                subprocess.run(
                    ["curl", "-L", "-o", str(output_path), url],
                    check=True,
                    capture_output=True,
                )
            else:
                subprocess.run(
                    ["wget", "-O", str(output_path), url],
                    check=True,
                    capture_output=True,
                )
            return output_path
        except (subprocess.CalledProcessError, FileNotFoundError):
            continue

    raise RuntimeError(
        f"Failed to download {url}. Please install curl or wget."
    )


def write_version_file(
    output_dir: Path,
    bed_name: str,
    version_info: str,
) -> None:
    """Write a .version file for a BED file.

    Args:
        output_dir: Directory to write the version file.
        bed_name: BED filename.
        version_info: Version information string.
    """
    version_path = output_dir / f"{bed_name}.version"
    with open(version_path, "w") as f:
        f.write(version_info + "\n")
    logger.info("Version file: %s -> %s", version_path, version_info)


def build_gencode_beds(gtf_path: Path, output_dir: Path) -> None:
    """Build GENCODE gene, ncRNA, coding exon+UTR, and sub-region BED files from GTF.

    Args:
        gtf_path: Path to the GENCODE GTF file (may be gzipped).
        output_dir: Directory to write output BED files.
    """
    logger.info("Building GENCODE BEDs from %s", gtf_path)

    # Gene loci: protein-coding genes (full transcript span)
    gene_loci: dict[str, list[tuple[int, int]]] = {}
    # ncRNA loci
    ncrna_loci: dict[str, list[tuple[int, int]]] = {}
    # Coding exon+UTR intervals per gene
    coding_exon_utr: dict[str, dict[str, list[tuple[int, int]]]] = {}
    # Sub-region intervals
    five_utr: dict[str, list[tuple[int, int]]] = {}
    cds_regions: dict[str, list[tuple[int, int]]] = {}
    three_utr: dict[str, list[tuple[int, int]]] = {}
    # Exon coordinates per transcript for splice-site derivation
    transcript_exons: dict[str, dict[str, list[tuple[int, int]]]] = {}
    # Transcript CDS ranges for UTR classification (tid -> (strand, cds_min, cds_max))
    transcript_cds: dict[str, tuple[str, int, int]] = {}
    # Deferred UTR entries: (chrom, bed_start, bed_end, transcript_id, strand)
    utr_entries: list[tuple[str, int, int, str, str]] = []

    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom_raw = parts[0]
            feature = parts[2]
            start = int(parts[3])  # GTF is 1-based
            end = int(parts[4])    # GTF is 1-based, closed
            strand = parts[6]
            attributes_str = parts[8]

            # Skip non-standard chromosomes
            chrom = BedUtils.normalize_chrom(chrom_raw)
            if chrom not in {
                "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
                "chr21", "chr22", "chrX", "chrY", "chrM",
            }:
                continue

            # Parse gene_type and gene_name from attributes
            gene_type = _parse_gtf_attribute(attributes_str, "gene_type")
            if not gene_type:
                gene_type = _parse_gtf_attribute(attributes_str, "gene_biotype")
            gene_name = _parse_gtf_attribute(attributes_str, "gene_name")
            transcript_id = _parse_gtf_attribute(attributes_str, "transcript_id")

            # Convert GTF 1-based closed to BED 0-based half-open
            bed_start = start - 1
            bed_end = end

            # Gene loci (full transcript span for gene-level records)
            if feature == "gene":
                if gene_type == PROTEIN_CODING_BIOTYPE:
                    if chrom not in gene_loci:
                        gene_loci[chrom] = []
                    gene_loci[chrom].append((bed_start, bed_end))
                elif gene_type in NCRNA_BIOTYPES:
                    if chrom not in ncrna_loci:
                        ncrna_loci[chrom] = []
                    ncrna_loci[chrom].append((bed_start, bed_end))

            # Collect exon/UTR intervals for protein-coding genes
            if feature in CODING_FEATURE_TYPES and gene_type == PROTEIN_CODING_BIOTYPE:
                if gene_name:
                    if gene_name not in coding_exon_utr:
                        coding_exon_utr[gene_name] = {}
                    if chrom not in coding_exon_utr[gene_name]:
                        coding_exon_utr[gene_name][chrom] = []
                    coding_exon_utr[gene_name][chrom].append((bed_start, bed_end))

            # Sub-region collection for protein-coding genes
            if gene_type == PROTEIN_CODING_BIOTYPE:
                if feature == "five_prime_utr":
                    if chrom not in five_utr:
                        five_utr[chrom] = []
                    five_utr[chrom].append((bed_start, bed_end))
                elif feature == "three_prime_utr":
                    if chrom not in three_utr:
                        three_utr[chrom] = []
                    three_utr[chrom].append((bed_start, bed_end))
                elif feature == "UTR":
                    # Defer classification until we know CDS ranges
                    if transcript_id:
                        utr_entries.append((chrom, bed_start, bed_end, transcript_id, strand))
                elif feature == "CDS":
                    if chrom not in cds_regions:
                        cds_regions[chrom] = []
                    cds_regions[chrom].append((bed_start, bed_end))
                    if transcript_id:
                        if transcript_id not in transcript_cds:
                            transcript_cds[transcript_id] = (strand, bed_start, bed_end)
                        else:
                            prev_strand, prev_min, prev_max = transcript_cds[transcript_id]
                            transcript_cds[transcript_id] = (
                                prev_strand,
                                min(prev_min, bed_start),
                                max(prev_max, bed_end),
                            )
                elif feature == "exon" and transcript_id:
                    if transcript_id not in transcript_exons:
                        transcript_exons[transcript_id] = {}
                    if chrom not in transcript_exons[transcript_id]:
                        transcript_exons[transcript_id][chrom] = []
                    transcript_exons[transcript_id][chrom].append(
                        (bed_start, bed_end)
                    )

    # Classify deferred UTR entries
    for chrom, bed_start, bed_end, transcript_id, strand in utr_entries:
        if transcript_id not in transcript_cds:
            continue
        _, cds_min, cds_max = transcript_cds[transcript_id]
        if strand == "+":
            if bed_end <= cds_min:
                if chrom not in five_utr:
                    five_utr[chrom] = []
                five_utr[chrom].append((bed_start, bed_end))
            elif bed_start >= cds_max:
                if chrom not in three_utr:
                    three_utr[chrom] = []
                three_utr[chrom].append((bed_start, bed_end))
        else:  # strand == "-"
            if bed_end <= cds_min:
                if chrom not in three_utr:
                    three_utr[chrom] = []
                three_utr[chrom].append((bed_start, bed_end))
            elif bed_start >= cds_max:
                if chrom not in five_utr:
                    five_utr[chrom] = []
                five_utr[chrom].append((bed_start, bed_end))

    # Sort and merge gene loci
    for chrom in gene_loci:
        gene_loci[chrom].sort(key=lambda x: x[0])
        gene_loci[chrom] = BedUtils.merge_intervals(gene_loci[chrom])

    BedUtils.write_bed(gene_loci, output_dir / GENCODE_GENE_BED)
    write_version_file(output_dir, GENCODE_GENE_BED, "GENCODE v44")

    # Sort and merge ncRNA loci
    for chrom in ncrna_loci:
        ncrna_loci[chrom].sort(key=lambda x: x[0])
        ncrna_loci[chrom] = BedUtils.merge_intervals(ncrna_loci[chrom])

    BedUtils.write_bed(ncrna_loci, output_dir / GENCODE_NCRNA_BED)
    write_version_file(output_dir, GENCODE_NCRNA_BED, "GENCODE v44")

    # Merge exon+UTR intervals per gene, then across all genes
    all_coding_intervals: dict[str, list[tuple[int, int]]] = {}
    for gene_name, chroms in coding_exon_utr.items():
        for chrom, intervals in chroms.items():
            intervals.sort(key=lambda x: x[0])
            merged = BedUtils.merge_intervals(intervals)
            if chrom not in all_coding_intervals:
                all_coding_intervals[chrom] = []
            all_coding_intervals[chrom].extend(merged)

    # Sort and merge across genes
    for chrom in all_coding_intervals:
        all_coding_intervals[chrom].sort(key=lambda x: x[0])
        all_coding_intervals[chrom] = BedUtils.merge_intervals(
            all_coding_intervals[chrom]
        )

    BedUtils.write_bed(all_coding_intervals, output_dir / GENCODE_CODING_EXON_UTR_BED)
    write_version_file(output_dir, GENCODE_CODING_EXON_UTR_BED, "GENCODE v44")

    # Sort and merge sub-region BEDs
    for chrom in five_utr:
        five_utr[chrom].sort(key=lambda x: x[0])
        five_utr[chrom] = BedUtils.merge_intervals(five_utr[chrom])
    BedUtils.write_bed(five_utr, output_dir / GENCODE_5UTR_BED)
    write_version_file(output_dir, GENCODE_5UTR_BED, "GENCODE v44")

    for chrom in cds_regions:
        cds_regions[chrom].sort(key=lambda x: x[0])
        cds_regions[chrom] = BedUtils.merge_intervals(cds_regions[chrom])
    BedUtils.write_bed(cds_regions, output_dir / GENCODE_CDS_BED)
    write_version_file(output_dir, GENCODE_CDS_BED, "GENCODE v44")

    for chrom in three_utr:
        three_utr[chrom].sort(key=lambda x: x[0])
        three_utr[chrom] = BedUtils.merge_intervals(three_utr[chrom])
    BedUtils.write_bed(three_utr, output_dir / GENCODE_3UTR_BED)
    write_version_file(output_dir, GENCODE_3UTR_BED, "GENCODE v44")

    # Derive splice sites from exon boundaries (±3 bp)
    splice_sites: dict[str, list[tuple[int, int]]] = {}
    for transcript_id, chroms in transcript_exons.items():
        for chrom, exons in chroms.items():
            if len(exons) < 2:
                continue
            exons_sorted = sorted(exons, key=lambda x: x[0])
            for i in range(len(exons_sorted) - 1):
                _, exon_end = exons_sorted[i]
                next_start, _ = exons_sorted[i + 1]
                # Splice donor: last 3 bp of exon
                donor_start = max(0, exon_end - 3)
                donor_end = exon_end + 3
                # Splice acceptor: first 3 bp of next exon
                acceptor_start = max(0, next_start - 3)
                acceptor_end = next_start + 3
                if chrom not in splice_sites:
                    splice_sites[chrom] = []
                splice_sites[chrom].append((donor_start, donor_end))
                splice_sites[chrom].append((acceptor_start, acceptor_end))

    for chrom in splice_sites:
        splice_sites[chrom].sort(key=lambda x: x[0])
        splice_sites[chrom] = BedUtils.merge_intervals(splice_sites[chrom])
    BedUtils.write_bed(splice_sites, output_dir / GENCODE_SPLICE_BED)
    write_version_file(output_dir, GENCODE_SPLICE_BED, "GENCODE v44")

    logger.info(
        "GENCODE BEDs built: %d gene intervals, %d ncRNA intervals, "
        "%d 5'UTR, %d CDS, %d 3'UTR, %d splice",
        sum(len(v) for v in gene_loci.values()),
        sum(len(v) for v in ncrna_loci.values()),
        sum(len(v) for v in five_utr.values()),
        sum(len(v) for v in cds_regions.values()),
        sum(len(v) for v in three_utr.values()),
        sum(len(v) for v in splice_sites.values()),
    )


def _parse_gtf_attribute(attributes_str: str, key: str) -> str:
    """Parse a key from a GTF attribute string.

    GTF attributes are formatted as: key "value"; key "value"; ...

    Args:
        attributes_str: The attributes column from a GTF line.
        key: Attribute key to find.

    Returns:
        The attribute value, or empty string if not found.
    """
    # Handle both quote formats: key "value" and key "value";
    for part in attributes_str.split(";"):
        part = part.strip()
        if part.startswith(key + " ") or part.startswith(key + '"'):
            # Remove key and quotes
            value = part[len(key):].strip().strip('"')
            return value
    return ""


def build_encode_beds(tsv_path: Path, output_dir: Path) -> None:
    """Build ENCODE cCRE BED files from SCREEN TSV/BED data.

    Supports two input formats:
    - Legacy SCREEN TSV: type in 6th column (0-based index 5)
    - ENCODE portal BED (ENCFF420VPZ): type in 10th column (0-based index 9)

    Generates three files:
    - encode_screen_v3_ccres.bed: All cCRE types (3 columns)
    - encode_screen_v3_pls_pels.bed: Only PLS and pELS types (3 columns)
    - encode_screen_v3_balanced.bed: PLS + pELS + dELS + CTCF (4 columns)

    Args:
        tsv_path: Path to ENCODE SCREEN TSV/BED file (may be gzipped).
        output_dir: Directory to write output BED files.
    """
    logger.info("Building ENCODE BEDs from %s", tsv_path)

    ccre_intervals: dict[str, list[tuple[int, int]]] = {}
    pls_pels_intervals: dict[str, list[tuple[int, int]]] = {}
    balanced_intervals: dict[str, list[tuple[int, int, str]]] = {}

    # Auto-detect type column by peeking at the first data line
    type_col = 5  # default legacy SCREEN TSV
    opener = gzip.open if str(tsv_path).endswith(".gz") else open
    with opener(tsv_path, "rt") as f:
        header_skipped = False
        for line in f:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            if not header_skipped:
                # Heuristic: skip header only if it doesn't look like a chrom line
                if not parts[0].startswith("chr"):
                    header_skipped = True
                    continue
                header_skipped = True
                # Detect format: ENCODE portal BED has cCRE type in column 10
                if len(parts) >= 10 and parts[9] in {
                    "PLS", "pELS", "dELS", "CTCF", "CTCF-only",
                    "CA-CTCF", "CA", "CA-H3K4me3", "CA-TF", "TF",
                }:
                    type_col = 9
                    logger.info("Detected ENCODE portal BED format (type column 10)")
                elif len(parts) >= 6 and parts[5].split(",")[0].strip() in {
                    "PLS", "pELS", "dELS", "CTCF",
                }:
                    logger.info("Detected legacy SCREEN TSV format (type column 6)")

            if len(parts) < 3:
                continue

            chrom_raw = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue

            # cCRE type
            type_field = parts[type_col] if len(parts) > type_col else ""
            ccre_type = type_field.split(",")[0].strip()
            # Normalize CTCF-only -> CTCF for balanced preset
            ccre_type_balanced = "CTCF" if ccre_type in ("CTCF", "CTCF-only", "CA-CTCF") else ccre_type

            chrom = BedUtils.normalize_chrom(chrom_raw)

            # All cCREs
            if chrom not in ccre_intervals:
                ccre_intervals[chrom] = []
            ccre_intervals[chrom].append((start, end))

            # PLS/pELS only (per U2 decision: filter by cCRE type)
            if ccre_type in ("PLS", "pELS"):
                if chrom not in pls_pels_intervals:
                    pls_pels_intervals[chrom] = []
                pls_pels_intervals[chrom].append((start, end))

            # Balanced preset: PLS + pELS + dELS + CTCF with type in 4th column
            if ccre_type in ("PLS", "pELS", "dELS") or ccre_type_balanced == "CTCF":
                if chrom not in balanced_intervals:
                    balanced_intervals[chrom] = []
                balanced_intervals[chrom].append((start, end, ccre_type_balanced))

    # Sort and merge
    for chrom in ccre_intervals:
        ccre_intervals[chrom].sort(key=lambda x: x[0])
        ccre_intervals[chrom] = BedUtils.merge_intervals(ccre_intervals[chrom])

    for chrom in pls_pels_intervals:
        pls_pels_intervals[chrom].sort(key=lambda x: x[0])
        pls_pels_intervals[chrom] = BedUtils.merge_intervals(
            pls_pels_intervals[chrom]
        )

    for chrom in balanced_intervals:
        balanced_intervals[chrom].sort(key=lambda x: x[0])

    BedUtils.write_bed(ccre_intervals, output_dir / ENCODE_CCRE_BED)
    write_version_file(output_dir, ENCODE_CCRE_BED, "ENCODE SCREEN v3")

    BedUtils.write_bed(pls_pels_intervals, output_dir / ENCODE_PLS_PELS_BED)
    write_version_file(output_dir, ENCODE_PLS_PELS_BED, "ENCODE SCREEN v3 (PLS/pELS)")

    _write_named_bed(balanced_intervals, output_dir / ENCODE_BALANCED_BED)
    write_version_file(
        output_dir, ENCODE_BALANCED_BED, "ENCODE SCREEN v3 (PLS/pELS/dELS/CTCF)"
    )

    logger.info(
        "ENCODE BEDs built: %d cCRE intervals, %d PLS/pELS intervals, "
        "%d balanced intervals",
        sum(len(v) for v in ccre_intervals.values()),
        sum(len(v) for v in pls_pels_intervals.values()),
        sum(len(v) for v in balanced_intervals.values()),
    )


def build_fantom5_bed(
    enhancer_path: Path,
    cage_path: Path,
    output_dir: Path,
) -> None:
    """Build FANTOM5 enhancers/promoters BED file.

    Merges FANTOM5 enhancers and CAGE peak (promoter/TSS) intervals into a
    single BED for the regulatory-balanced and comprehensive presets.

    Args:
        enhancer_path: Path to FANTOM5 enhancer BED/BED.GZ file.
        cage_path: Path to FANTOM5 CAGE peak BED/BED.GZ file.
        output_dir: Directory to write output BED file.
    """
    logger.info("Building FANTOM5 BED from %s and %s", enhancer_path, cage_path)

    intervals: dict[str, list[tuple[int, int]]] = {}

    for src_path in (enhancer_path, cage_path):
        if not src_path or not src_path.exists():
            continue
        opener = gzip.open if str(src_path).endswith(".gz") else open
        with opener(src_path, "rt") as f:
            for line in f:
                if line.startswith("#") or line.startswith("track"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    parts = line.rstrip("\n").split()
                if len(parts) < 3:
                    continue
                chrom = BedUtils.normalize_chrom(parts[0])
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue

                if chrom not in intervals:
                    intervals[chrom] = []
                intervals[chrom].append((start, end))

    for chrom in intervals:
        intervals[chrom].sort(key=lambda x: x[0])
        intervals[chrom] = BedUtils.merge_intervals(intervals[chrom])

    BedUtils.write_bed(intervals, output_dir / FANTOM5_BED)
    write_version_file(output_dir, FANTOM5_BED, "FANTOM5 hg38 enhancers + CAGE peaks")

    logger.info(
        "FANTOM5 BED built: %d intervals",
        sum(len(v) for v in intervals.values()),
    )


def build_vista_bed(tsv_path: Path, output_dir: Path) -> None:
    """Build Vista Enhancer Browser BED file.

    Args:
        tsv_path: Path to Vista data file.
        output_dir: Directory to write output BED file.
    """
    logger.info("Building Vista BED from %s", tsv_path)

    intervals: dict[str, list[tuple[int, int]]] = {}
    opener = gzip.open if str(tsv_path).endswith(".gz") else open
    with opener(tsv_path, "rt") as f:
        for line in f:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                parts = line.rstrip("\n").split()
            if len(parts) < 3:
                continue
            chrom = BedUtils.normalize_chrom(parts[0])
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue

            if chrom not in intervals:
                intervals[chrom] = []
            intervals[chrom].append((start, end))

    for chrom in intervals:
        intervals[chrom].sort(key=lambda x: x[0])
        intervals[chrom] = BedUtils.merge_intervals(intervals[chrom])

    BedUtils.write_bed(intervals, output_dir / VISTA_BED)
    write_version_file(output_dir, VISTA_BED, "Vista Enhancer Browser")

    logger.info(
        "Vista BED built: %d intervals",
        sum(len(v) for v in intervals.values()),
    )


def _write_named_bed(
    intervals: dict[str, list[tuple[int, int, str]]],
    output: Path,
) -> None:
    """Write named intervals to a 4-column BED file.

    Args:
        intervals: Dictionary {chrom: [(start, end, name), ...]}.
        output: Output file path.
    """
    output.parent.mkdir(parents=True, exist_ok=True)

    def chrom_sort_key(c: str) -> tuple[int, str]:
        name = c.replace("chr", "")
        if name.isdigit():
            return (0, name.zfill(2))
        elif name == "X":
            return (1, "X")
        elif name == "Y":
            return (2, "Y")
        elif name == "M":
            return (3, "M")
        else:
            return (4, name)

    sorted_chroms = sorted(intervals.keys(), key=chrom_sort_key)
    with open(output, "w") as f:
        for chrom in sorted_chroms:
            for start, end, name in intervals[chrom]:
                f.write(f"{chrom}\t{start}\t{end}\t{name}\n")


def _clnrevstat_to_star(clnrevstat: str) -> int:
    """Map ClinVar CLNREVSTAT to star level.

    Star levels:
    - 0: no_assertion, no_criteria, conflicting (filtered out)
    - 1: criteria_provided,_single_submitter
    - 2: criteria_provided,_multiple_submitters,_no_conflicts
    - 3: reviewed_by_expert_panel
    - 4: practice_guideline

    Args:
        clnrevstat: CLNREVSTAT field value.

    Returns:
        Integer star level (0-4).
    """
    if not clnrevstat:
        return 0
    rev_lower = clnrevstat.lower()
    if "practice_guideline" in rev_lower:
        return 4
    if "reviewed_by_expert_panel" in rev_lower:
        return 3
    if "criteria_provided,_multiple_submitters,_no_conflicts" in rev_lower:
        return 2
    if "criteria_provided,_conflicting_interpretations" in rev_lower:
        return 0
    if "criteria_provided" in rev_lower:
        return 1
    return 0


def build_clinvar_bed(vcf_path: Path, output_dir: Path) -> None:
    """Build ClinVar P/LP BED file from ClinVar VCF.

    Filters for variants where CLNSIG contains 'Pathogenic' or
    'Likely_pathogenic' (including combined states like
    'Pathogenic/Likely_pathogenic'). Writes a 4-column BED where the
    4th column is the ClinVar review-status star level (1-4). Zero-star
    records are filtered out.

    Args:
        vcf_path: Path to ClinVar VCF (may be gzipped).
        output_dir: Directory to write output BED file.
    """
    logger.info("Building ClinVar BED from %s", vcf_path)

    intervals: dict[str, list[tuple[int, int, str]]] = {}
    opener = gzip.open if str(vcf_path).endswith(".gz") else open

    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom = BedUtils.normalize_chrom(parts[0])
            try:
                pos = int(parts[1])  # 1-based VCF position
            except ValueError:
                continue

            info = parts[7]

            # Check CLNSIG for Pathogenic or Likely_pathogenic
            clnsig = ""
            clnrevstat = ""
            for field_str in info.split(";"):
                if field_str.startswith("CLNSIG="):
                    clnsig = field_str[7:]
                elif field_str.startswith("CLNREVSTAT="):
                    clnrevstat = field_str[11:]

            # Per U3: include any CLNSIG containing Pathogenic or Likely_pathogenic
            is_pathogenic = (
                "Pathogenic" in clnsig or "Likely_pathogenic" in clnsig
            )
            if not is_pathogenic:
                continue

            # Filter out zero-star records
            star_level = _clnrevstat_to_star(clnrevstat)
            if star_level == 0:
                continue

            # Convert VCF 1-based position to BED 0-based half-open
            ref = parts[3]
            bed_start = pos - 1
            bed_end = pos  # At minimum, cover the variant position

            # For longer REF alleles, extend to cover the full reference
            if len(ref) > 1:
                bed_end = pos - 1 + len(ref)

            if chrom not in intervals:
                intervals[chrom] = []
            intervals[chrom].append((bed_start, bed_end, str(star_level)))

    # Sort (no merge: preserve star-level distinctions at overlapping positions)
    for chrom in intervals:
        intervals[chrom].sort(key=lambda x: x[0])

    _write_named_bed(intervals, output_dir / CLINVAR_BED)
    write_version_file(
        output_dir, CLINVAR_BED, f"ClinVar {datetime.now(timezone.utc).strftime('%Y-%m-%d')}"
    )

    logger.info(
        "ClinVar BED built: %d P/LP intervals",
        sum(len(v) for v in intervals.values()),
    )


def build_omim_bed(output_dir: Path) -> None:
    """Build empty OMIM BED file (v2.0 placeholder).

    Args:
        output_dir: Directory to write the empty BED file.
    """
    logger.info("Building OMIM BED placeholder (empty, v2.0)")
    omim_path = output_dir / OMIM_BED
    with open(omim_path, "w") as f:
        f.write('track name="OMIM_pathogenic"\n')
    write_version_file(output_dir, OMIM_BED, "OMIM placeholder (v2.0)")


def build_ensembl_regulatory_bed(gff3_path: Path, output_dir: Path) -> None:
    """Build Ensembl Regulatory Build BED file from GFF3.

    Parses Ensembl regulatory features GFF3, keeps selected feature types,
    converts to 0-based half-open BED with type in the 4th column.

    Args:
        gff3_path: Path to the Ensembl GFF3.gz file.
        output_dir: Directory to write output BED file.
    """
    logger.info("Building Ensembl Regulatory BED from %s", gff3_path)

    feature_types = {
        "promoter",
        "enhancer",
        "CTCF_binding_site",
        "open_chromatin_region",
        "TF_binding_site",
    }

    # Map Ensembl type names to shorter annotation tag names
    type_map = {
        "promoter": "promoter",
        "enhancer": "enhancer",
        "CTCF_binding_site": "ctcf",
        "open_chromatin_region": "open_chromatin",
        "TF_binding_site": "tf_binding",
    }

    intervals: dict[str, list[tuple[int, int, str]]] = {}
    opener = gzip.open if str(gff3_path).endswith(".gz") else open

    with opener(gff3_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom_raw = parts[0]
            try:
                start = int(parts[3])  # GFF3 is 1-based inclusive
                end = int(parts[4])    # GFF3 is 1-based inclusive
            except ValueError:
                continue

            feature_type = parts[2].strip()

            # Also check attributes for type if column 3 is generic
            if feature_type not in feature_types:
                attributes = parts[8]
                for attr in attributes.split(";"):
                    attr = attr.strip()
                    if attr.startswith("type=") or attr.startswith("Type="):
                        attr_type = attr.split("=", 1)[1]
                        if attr_type in feature_types:
                            feature_type = attr_type
                        break

            if feature_type not in feature_types:
                continue

            chrom = BedUtils.normalize_chrom(chrom_raw)
            bed_start = start - 1  # Convert to 0-based half-open
            bed_end = end
            bed_name = type_map.get(feature_type, feature_type)

            if chrom not in intervals:
                intervals[chrom] = []
            intervals[chrom].append((bed_start, bed_end, bed_name))

    # Sort by start coordinate (no merge: preserve type distinctions)
    for chrom in intervals:
        intervals[chrom].sort(key=lambda x: x[0])

    _write_named_bed(intervals, output_dir / ENSEMBL_REGULATORY_BED)
    write_version_file(output_dir, ENSEMBL_REGULATORY_BED, "Ensembl Regulatory Build v115")

    logger.info(
        "Ensembl Regulatory BED built: %d intervals",
        sum(len(v) for v in intervals.values()),
    )


def build_merged_retained(output_dir: Path) -> None:
    """Build merged_retained_regions.bed from comprehensive preset BEDs.

    This is the pre-merged comprehensive region set for quick filtering.

    Args:
        output_dir: Directory containing individual BED files and output.
    """
    logger.info("Building merged retained regions BED")

    bed_paths = [
        output_dir / GENCODE_GENE_BED,
        output_dir / GENCODE_NCRNA_BED,
        output_dir / ENCODE_CCRE_BED,
        output_dir / FANTOM5_BED,
        output_dir / VISTA_BED,
        output_dir / ENSEMBL_REGULATORY_BED,
    ]

    # Filter to existing files
    existing_paths = [p for p in bed_paths if p.exists()]
    if not existing_paths:
        logger.warning("No region BED files found for merging")
        return

    BedUtils.merge_bed_files(existing_paths, output_dir / MERGED_RETAINED_BED)
    write_version_file(
        output_dir,
        MERGED_RETAINED_BED,
        f"Merged from comprehensive preset {datetime.now(timezone.utc).strftime('%Y-%m-%d')}",
    )

    logger.info("Merged retained regions BED built")


def main() -> None:
    """Main entry point for the build script."""
    parser = argparse.ArgumentParser(
        description="Build reference BED files for dgra-prefilter"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("refs"),
        help="Output directory for BED files (default: refs/)",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        default=False,
        help="Skip downloading raw data, use existing files",
    )
    parser.add_argument(
        "--gencode-gtf",
        type=Path,
        default=None,
        help="Path to GENCODE GTF file (skip download if provided)",
    )
    parser.add_argument(
        "--encode-tsv",
        type=Path,
        default=None,
        help="Path to ENCODE SCREEN TSV file (skip download if provided)",
    )
    parser.add_argument(
        "--fantom5-file",
        type=Path,
        default=None,
        help="Path to FANTOM5 data file (skip download if provided)",
    )
    parser.add_argument(
        "--vista-file",
        type=Path,
        default=None,
        help="Path to Vista data file (skip download if provided)",
    )
    parser.add_argument(
        "--clinvar-vcf",
        type=Path,
        default=None,
        help="Path to ClinVar VCF file (skip download if provided)",
    )
    parser.add_argument(
        "--ensembl-file",
        type=Path,
        default=None,
        help="Path to Ensembl Regulatory Build GFF3 file (skip download if provided)",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="Enable verbose logging",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    )

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    cache_dir = output_dir / ".cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    tmp_path = cache_dir

    # Build GENCODE BEDs: prefer shared asset, then CLI arg, then download
    if args.gencode_gtf:
        gencode_path = args.gencode_gtf
    elif SHARED_GENCODE_GTF.exists():
        logger.info("Using shared GENCODE GTF: %s", SHARED_GENCODE_GTF)
        gencode_path = SHARED_GENCODE_GTF
    elif not args.skip_download:
        gencode_path = tmp_path / "gencode.v44.annotation.gtf.gz"
        _download_file(GENCODE_GTF_URL, gencode_path)
    else:
        logger.warning("No GENCODE GTF provided and --skip-download set, skipping")
        gencode_path = None

    if gencode_path and gencode_path.exists():
        build_gencode_beds(gencode_path, output_dir)
    else:
        logger.warning("GENCODE GTF not available, creating empty BED files")
        for bed_name in [
            GENCODE_GENE_BED,
            GENCODE_NCRNA_BED,
            GENCODE_CODING_EXON_UTR_BED,
            GENCODE_5UTR_BED,
            GENCODE_CDS_BED,
            GENCODE_3UTR_BED,
            GENCODE_SPLICE_BED,
        ]:
            (output_dir / bed_name).touch()
            write_version_file(output_dir, bed_name, "missing")

    # Build ENCODE BEDs
    if args.encode_tsv:
        encode_path = args.encode_tsv
    elif not args.skip_download:
        encode_path = tmp_path / "encode_screen_v3.tsv.gz"
        try:
            _download_file(ENCODE_SCREEN_URL, encode_path)
        except RuntimeError:
            logger.warning("ENCODE SCREEN download failed")
            encode_path = None
    else:
        encode_path = None

    if encode_path and encode_path.exists():
        try:
            build_encode_beds(encode_path, output_dir)
        except (OSError, gzip.BadGzipFile, ValueError) as exc:
            logger.warning("ENCODE SCREEN data invalid (%s), creating empty BED files", exc)
            for bed_name in [ENCODE_CCRE_BED, ENCODE_PLS_PELS_BED, ENCODE_BALANCED_BED]:
                bed_path = output_dir / bed_name
                if bed_path.exists() and bed_path.stat().st_size > 0:
                    logger.info("Preserving existing BED: %s", bed_name)
                else:
                    bed_path.touch()
                    write_version_file(output_dir, bed_name, "missing")
    else:
        logger.warning("ENCODE SCREEN not available, creating empty BED files")
        for bed_name in [ENCODE_CCRE_BED, ENCODE_PLS_PELS_BED, ENCODE_BALANCED_BED]:
            bed_path = output_dir / bed_name
            if bed_path.exists() and bed_path.stat().st_size > 0:
                logger.info("Preserving existing BED: %s", bed_name)
            else:
                bed_path.touch()
                write_version_file(output_dir, bed_name, "missing")

    # Build FANTOM5 BED (enhancers + CAGE peaks/promoters)
    if args.fantom5_file:
        fantom5_enhancer_path = args.fantom5_file
    elif not args.skip_download:
        fantom5_enhancer_path = tmp_path / "fantom5_enhancers.bed.gz"
        try:
            _download_file(FANTOM5_ENHANCER_URL, fantom5_enhancer_path)
        except RuntimeError:
            logger.warning("FANTOM5 enhancer download failed")
            fantom5_enhancer_path = None
    else:
        fantom5_enhancer_path = None

    if args.skip_download:
        fantom5_cage_path = tmp_path / "fantom5_cage_peaks.bed.gz"
        if not fantom5_cage_path.exists():
            fantom5_cage_path = None
    elif fantom5_enhancer_path:
        fantom5_cage_path = tmp_path / "fantom5_cage_peaks.bed.gz"
        try:
            _download_file(FANTOM5_CAGE_URL, fantom5_cage_path)
        except RuntimeError:
            logger.warning("FANTOM5 CAGE peak download failed")
            fantom5_cage_path = None
    else:
        fantom5_cage_path = None

    if fantom5_enhancer_path and fantom5_enhancer_path.exists():
        build_fantom5_bed(fantom5_enhancer_path, fantom5_cage_path, output_dir)
    else:
        logger.warning("FANTOM5 not available, creating empty BED file")
        (output_dir / FANTOM5_BED).touch()
        write_version_file(output_dir, FANTOM5_BED, "missing")

    # Build Vista BED (download bigBed and convert with bigBedToBed if available)
    if args.vista_file:
        vista_path = args.vista_file
    elif not args.skip_download:
        vista_bb_path = tmp_path / "vista_enhancers.bb"
        try:
            _download_file(VISTA_URL, vista_bb_path)
        except RuntimeError:
            logger.warning("Vista bigBed download failed")
            vista_bb_path = None
        vista_path = tmp_path / "vista_enhancers.bed"
        if vista_bb_path and vista_bb_path.exists():
            bigbedtobed = shutil.which("bigBedToBed")
            if not bigbedtobed:
                # Try common UCSC binary install paths
                for candidate in [
                    Path.home() / "bin" / "bigBedToBed",
                    Path("/usr/local/bin/bigBedToBed"),
                ]:
                    if candidate.exists():
                        bigbedtobed = str(candidate)
                        break
            if bigbedtobed:
                try:
                    subprocess.run(
                        [bigbedtobed, str(vista_bb_path), str(vista_path)],
                        check=True,
                        capture_output=True,
                    )
                except (subprocess.CalledProcessError, FileNotFoundError):
                    logger.warning("bigBedToBed conversion failed")
                    vista_path = None
            else:
                logger.warning("bigBedToBed not found; cannot convert Vista bigBed")
                vista_path = None
    else:
        vista_path = None

    if vista_path and vista_path.exists():
        build_vista_bed(vista_path, output_dir)
    else:
        logger.warning("Vista not available, creating empty BED file")
        (output_dir / VISTA_BED).touch()
        write_version_file(output_dir, VISTA_BED, "missing")

    # Build ClinVar BED
    if args.clinvar_vcf:
        clinvar_path = args.clinvar_vcf
    elif not args.skip_download:
        clinvar_path = tmp_path / "clinvar.vcf.gz"
        try:
            _download_file(CLINVAR_VCF_URL, clinvar_path)
        except RuntimeError:
            logger.warning("ClinVar download failed")
            clinvar_path = None
    else:
        clinvar_path = None

    if clinvar_path and clinvar_path.exists():
        build_clinvar_bed(clinvar_path, output_dir)
    else:
        logger.warning("ClinVar not available, creating empty BED file")
        (output_dir / CLINVAR_BED).touch()
        write_version_file(output_dir, CLINVAR_BED, "missing")

    # Build Ensembl Regulatory BED
    if args.ensembl_file:
        ensembl_path = args.ensembl_file
    elif not args.skip_download:
        ensembl_path = tmp_path / "ensembl_regulatory.gff3.gz"
        try:
            _download_file(ENSEMBL_REGULATORY_URL, ensembl_path)
        except RuntimeError:
            logger.warning("Ensembl Regulatory Build download failed")
            ensembl_path = None
    else:
        ensembl_path = None

    if ensembl_path and ensembl_path.exists():
        build_ensembl_regulatory_bed(ensembl_path, output_dir)
    else:
        logger.warning("Ensembl Regulatory Build not available, creating empty BED file")
        (output_dir / ENSEMBL_REGULATORY_BED).touch()
        write_version_file(output_dir, ENSEMBL_REGULATORY_BED, "missing")

    # Build OMIM placeholder
    build_omim_bed(output_dir)

    # Build merged retained regions
    build_merged_retained(output_dir)

    logger.info("All reference BED files built in %s", output_dir)


if __name__ == "__main__":
    main()
