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
    ENCODE_CCRE_BED,
    ENCODE_PLS_PELS_BED,
    FANTOM5_BED,
    GENCODE_CODING_EXON_UTR_BED,
    GENCODE_GENE_BED,
    GENCODE_NCRNA_BED,
    MERGED_RETAINED_BED,
    OMIM_BED,
    VISTA_BED,
)

logger = logging.getLogger(__name__)

# Data source URLs
GENCODE_GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_44/gencode.v44.annotation.gtf.gz"
)
ENCODE_SCREEN_URL = (
    "https://screen.encodeproject.org/download/tsv/?"
    "format=gz&assembly=GRCh38&accession=ENCSR000AIZ&fileType=bed"
)
FANTOM5_URL = (
    "http://fantom.gsc.riken.jp/5/data/hg38/robust/"
    "hg38_enhancers.bed"
)
VISTA_URL = "https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?form=download"
CLINVAR_VCF_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/"
    "clinvar_20260530.vcf.gz"
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
    logger.info("Downloading %s to %s", url, output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

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
    """Build GENCODE gene, ncRNA, and coding exon+UTR BED files from GTF.

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

    logger.info(
        "GENCODE BEDs built: %d gene intervals, %d ncRNA intervals",
        sum(len(v) for v in gene_loci.values()),
        sum(len(v) for v in ncrna_loci.values()),
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
    """Build ENCODE cCRE BED files from SCREEN TSV data.

    Generates two files:
    - encode_screen_v3_ccres.bed: All cCRE types
    - encode_screen_v3_pls_pels.bed: Only PLS and pELS types

    Args:
        tsv_path: Path to ENCODE SCREEN TSV/TSV.GZ file.
        output_dir: Directory to write output BED files.
    """
    logger.info("Building ENCODE BEDs from %s", tsv_path)

    ccre_intervals: dict[str, list[tuple[int, int]]] = {}
    pls_pels_intervals: dict[str, list[tuple[int, int]]] = {}

    opener = gzip.open if str(tsv_path).endswith(".gz") else open
    with opener(tsv_path, "rt") as f:
        header_skipped = False
        for line in f:
            if line.startswith("#"):
                continue
            if not header_skipped:
                # Heuristic: skip header only if it doesn't look like a chrom line
                if not line.startswith("chr"):
                    header_skipped = True
                    continue
                header_skipped = True

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue

            chrom_raw = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue

            # cCRE type is in the 6th column, format "type,CTCF-state" or just "type"
            type_field = parts[5] if len(parts) > 5 else ""
            ccre_type = type_field.split(",")[0].strip()

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

    # Sort and merge
    for chrom in ccre_intervals:
        ccre_intervals[chrom].sort(key=lambda x: x[0])
        ccre_intervals[chrom] = BedUtils.merge_intervals(ccre_intervals[chrom])

    for chrom in pls_pels_intervals:
        pls_pels_intervals[chrom].sort(key=lambda x: x[0])
        pls_pels_intervals[chrom] = BedUtils.merge_intervals(
            pls_pels_intervals[chrom]
        )

    BedUtils.write_bed(ccre_intervals, output_dir / ENCODE_CCRE_BED)
    write_version_file(output_dir, ENCODE_CCRE_BED, "ENCODE SCREEN v3")

    BedUtils.write_bed(pls_pels_intervals, output_dir / ENCODE_PLS_PELS_BED)
    write_version_file(output_dir, ENCODE_PLS_PELS_BED, "ENCODE SCREEN v3 (PLS/pELS)")

    logger.info(
        "ENCODE BEDs built: %d cCRE intervals, %d PLS/pELS intervals",
        sum(len(v) for v in ccre_intervals.values()),
        sum(len(v) for v in pls_pels_intervals.values()),
    )


def build_fantom5_bed(tsv_path: Path, output_dir: Path) -> None:
    """Build FANTOM5 enhancers/promoters BED file.

    Args:
        tsv_path: Path to FANTOM5 data file (BED or TSV).
        output_dir: Directory to write output BED file.
    """
    logger.info("Building FANTOM5 BED from %s", tsv_path)

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

    BedUtils.write_bed(intervals, output_dir / FANTOM5_BED)
    write_version_file(output_dir, FANTOM5_BED, "FANTOM5")

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


def build_clinvar_bed(vcf_path: Path, output_dir: Path) -> None:
    """Build ClinVar P/LP BED file from ClinVar VCF.

    Filters for variants where CLNSIG contains 'Pathogenic' or
    'Likely_pathogenic' (including combined states like
    'Pathogenic/Likely_pathogenic'). Handles multi-allelic sites
    by first splitting with bcftools norm.

    Args:
        vcf_path: Path to ClinVar VCF (may be gzipped).
        output_dir: Directory to write output BED file.
    """
    logger.info("Building ClinVar BED from %s", vcf_path)

    intervals: dict[str, list[tuple[int, int]]] = {}
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
            for field_str in info.split(";"):
                if field_str.startswith("CLNSIG="):
                    clnsig = field_str[7:]
                    break

            # Per U3: include any CLNSIG containing Pathogenic or Likely_pathogenic
            is_pathogenic = (
                "Pathogenic" in clnsig or "Likely_pathogenic" in clnsig
            )
            if not is_pathogenic:
                continue

            # Convert VCF 1-based position to BED 0-based half-open
            # For SNPs: [pos-1, pos)
            # For indels: we'd need REF/ALT parsing, but for BED we use
            # a minimal interval [pos-1, pos) which covers the variant
            ref = parts[3]
            bed_start = pos - 1
            bed_end = pos  # At minimum, cover the variant position

            # For longer REF alleles, extend to cover the full reference
            if len(ref) > 1:
                bed_end = pos - 1 + len(ref)

            if chrom not in intervals:
                intervals[chrom] = []
            intervals[chrom].append((bed_start, bed_end))

    # Sort and merge
    for chrom in intervals:
        intervals[chrom].sort(key=lambda x: x[0])
        intervals[chrom] = BedUtils.merge_intervals(intervals[chrom])

    BedUtils.write_bed(intervals, output_dir / CLINVAR_BED)
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

    with tempfile.TemporaryDirectory(prefix="dgra_build_refs_") as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Build GENCODE BEDs
        if args.gencode_gtf:
            gencode_path = args.gencode_gtf
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
            for bed_name in [GENCODE_GENE_BED, GENCODE_NCRNA_BED, GENCODE_CODING_EXON_UTR_BED]:
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
                for bed_name in [ENCODE_CCRE_BED, ENCODE_PLS_PELS_BED]:
                    (output_dir / bed_name).touch()
                    write_version_file(output_dir, bed_name, "missing")
        else:
            logger.warning("ENCODE SCREEN not available, creating empty BED files")
            for bed_name in [ENCODE_CCRE_BED, ENCODE_PLS_PELS_BED]:
                (output_dir / bed_name).touch()
                write_version_file(output_dir, bed_name, "missing")

        # Build FANTOM5 BED
        if args.fantom5_file:
            fantom5_path = args.fantom5_file
        elif not args.skip_download:
            fantom5_path = tmp_path / "fantom5_enhancers.bed"
            try:
                _download_file(FANTOM5_URL, fantom5_path)
            except RuntimeError:
                logger.warning("FANTOM5 download failed")
                fantom5_path = None
        else:
            fantom5_path = None

        if fantom5_path and fantom5_path.exists():
            build_fantom5_bed(fantom5_path, output_dir)
        else:
            logger.warning("FANTOM5 not available, creating empty BED file")
            (output_dir / FANTOM5_BED).touch()
            write_version_file(output_dir, FANTOM5_BED, "missing")

        # Build Vista BED
        if args.vista_file:
            vista_path = args.vista_file
        elif not args.skip_download:
            vista_path = tmp_path / "vista_enhancers.tsv"
            try:
                _download_file(VISTA_URL, vista_path)
            except RuntimeError:
                logger.warning("Vista download failed")
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

        # Build OMIM placeholder
        build_omim_bed(output_dir)

        # Build merged retained regions
        build_merged_retained(output_dir)

    logger.info("All reference BED files built in %s", output_dir)


if __name__ == "__main__":
    main()
