#!/usr/bin/env python3
"""Generate synthetic demo VCFs for end-to-end pipeline testing."""

from __future__ import annotations

import argparse
import gzip
import os
import random
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

# GRCh38 chromosome lengths
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

BCFTOOLS = os.environ.get("BCFTOOLS", "/Users/zhaorongli/.micromamba/envs/bcftools/bin/bcftools")
BGZIP = os.environ.get("BGZIP", "/Users/zhaorongli/.micromamba/envs/bcftools/bin/bgzip")

NUCS = ["A", "C", "G", "T"]


def load_bed_intervals(bed_path: str | Path) -> dict[str, list[tuple[int, int]]]:
    """Load BED intervals by chromosome as 0-based [start, end)."""
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            intervals[chrom].append((start, end))
    for chrom in intervals:
        intervals[chrom].sort()
    return dict(intervals)


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge sorted overlapping intervals."""
    if not intervals:
        return []
    merged: list[tuple[int, int]] = [intervals[0]]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def overlaps(pos: int, intervals: list[tuple[int, int]]) -> bool:
    """Binary search for overlap with sorted intervals (0-based)."""
    lo, hi = 0, len(intervals)
    while lo < hi:
        mid = (lo + hi) // 2
        s, e = intervals[mid]
        if pos < s:
            hi = mid
        elif pos >= e:
            lo = mid + 1
        else:
            return True
    return False


def build_combined_intervals(
    gene_intervals: dict[str, list[tuple[int, int]]],
    reg_intervals: dict[str, list[tuple[int, int]]],
) -> dict[str, list[tuple[int, int]]]:
    """Build merged combined intervals from gene and regulatory BEDs."""
    combined: dict[str, list[tuple[int, int]]] = {}
    all_chroms = set(gene_intervals.keys()) | set(reg_intervals.keys())
    for chrom in all_chroms:
        raw: list[tuple[int, int]] = []
        raw.extend(gene_intervals.get(chrom, []))
        raw.extend(reg_intervals.get(chrom, []))
        raw.sort()
        combined[chrom] = merge_intervals(raw)
    return combined


def sample_position_inside(
    intervals: dict[str, list[tuple[int, int]]],
    chrom_lengths: dict[str, int],
    rng: random.Random,
) -> tuple[str, int]:
    """Sample a random position inside a BED interval."""
    chroms = list(intervals.keys())
    while True:
        chrom = rng.choice(chroms)
        chrom_ivs = intervals[chrom]
        if not chrom_ivs:
            continue
        start, end = rng.choice(chrom_ivs)
        if end - start <= 0:
            continue
        pos = rng.randint(start, end - 1)  # 0-based
        # Ensure position is within chromosome bounds
        if pos < 0 or pos >= chrom_lengths.get(chrom, 0):
            continue
        return chrom, pos


def sample_position_outside(
    intervals: dict[str, list[tuple[int, int]]],
    chrom_lengths: dict[str, int],
    rng: random.Random,
    min_distance: int = 10000,
) -> tuple[str, int]:
    """Sample a random position outside all BED intervals."""
    chroms = list(chrom_lengths.keys())
    max_attempts = 10000
    for _ in range(max_attempts):
        chrom = rng.choice(chroms)
        length = chrom_lengths[chrom]
        if length <= 0:
            continue
        pos = rng.randint(0, length - 1)
        chrom_ivs = intervals.get(chrom, [])
        # Check if inside any interval
        if overlaps(pos, chrom_ivs):
            continue
        # Check if too close to any interval boundary
        too_close = False
        for s, e in chrom_ivs:
            if abs(pos - s) < min_distance or abs(pos - e) < min_distance:
                too_close = True
                break
        if too_close:
            continue
        return chrom, pos
    # Fallback: pick a position near the end of a chromosome
    chrom = rng.choice(chroms)
    length = chrom_lengths[chrom]
    pos = max(0, length - rng.randint(100000, 500000))
    return chrom, pos


def generate_random_nuc(rng: random.Random) -> str:
    """Generate a random nucleotide."""
    return rng.choice(NUCS)


def generate_alt(ref: str, rng: random.Random) -> str:
    """Generate a random ALT different from REF."""
    alt = rng.choice(NUCS)
    while alt == ref:
        alt = rng.choice(NUCS)
    return alt


def generate_info(rng: random.Random) -> dict[str, str]:
    """Generate random but realistic INFO values."""
    return {
        "QD": f"{rng.uniform(1.5, 35.0):.3f}",
        "FS": f"{rng.uniform(0.0, 60.0):.3f}",
        "SOR": f"{rng.uniform(0.5, 3.5):.3f}",
        "MQ": f"{rng.uniform(25.0, 60.0):.1f}",
        "ReadPosRankSum": f"{rng.uniform(-5.0, 5.0):.3f}",
        "MQRankSum": f"{rng.uniform(-10.0, 10.0):.3f}",
        "BaseQRankSum": f"{rng.uniform(-10.0, 10.0):.3f}",
    }


def build_vcf_header(
    chrom_lengths: dict[str, int],
    sample_name: str = "DEMO",
    reference: str = "GRCh38",
) -> list[str]:
    """Build VCF 4.2 header lines."""
    lines: list[str] = []
    lines.append("##fileformat=VCFv4.2")
    lines.append(f"##reference={reference}")

    for chrom, length in sorted(chrom_lengths.items(), key=lambda x: (x[0][3:] if x[0].startswith("chr") else x[0])):
        lines.append(f"##contig=<ID={chrom},length={length}>")

    lines.append('##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">')
    lines.append('##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher\'s exact test to detect strand bias">')
    lines.append('##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">')
    lines.append('##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">')
    lines.append('##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">')
    lines.append('##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping quality">')
    lines.append('##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref base qualities">')
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    lines.append('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">')
    lines.append('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">')
    lines.append('##FILTER=<ID=PASS,Description="All filters passed">')
    lines.append('##FILTER=<ID=LowQual,Description="Low quality">')
    lines.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}")
    return lines


def generate_vcf_records(
    n_variants: int,
    chrom_lengths: dict[str, int],
    gene_intervals: dict[str, list[tuple[int, int]]],
    reg_intervals: dict[str, list[tuple[int, int]]],
    placement: str,
    gt_distribution: dict[str, float],
    ref_eq_alt_rate: float,
    low_qual_rate: float,
    rng: random.Random,
) -> list[str]:
    """Generate VCF record lines.

    Args:
        n_variants: Number of variants to generate.
        chrom_lengths: Chromosome lengths.
        gene_intervals: Gene loci intervals.
        reg_intervals: Regulatory intervals.
        placement: "gene", "regulatory", "gene_or_regulatory", "outside", or "mixed".
        gt_distribution: Probability distribution for genotypes.
        ref_eq_alt_rate: Fraction of variants with REF == ALT.
        low_qual_rate: Fraction of variants with FILTER=LowQual.
        rng: Random generator.

    Returns:
        List of VCF record lines.
    """
    combined = build_combined_intervals(gene_intervals, reg_intervals)
    records: list[str] = []
    gt_choices = list(gt_distribution.keys())
    gt_weights = [gt_distribution[k] for k in gt_choices]

    for i in range(n_variants):
        # Determine placement
        if placement == "gene_or_regulatory":
            chrom, pos_0based = sample_position_inside(combined, chrom_lengths, rng)
        elif placement == "mixed":
            if rng.random() < 0.5:
                chrom, pos_0based = sample_position_inside(combined, chrom_lengths, rng)
            else:
                chrom, pos_0based = sample_position_outside(combined, chrom_lengths, rng)
        elif placement == "gene":
            chrom, pos_0based = sample_position_inside(gene_intervals, chrom_lengths, rng)
        elif placement == "regulatory":
            chrom, pos_0based = sample_position_inside(reg_intervals, chrom_lengths, rng)
        else:
            chrom, pos_0based = sample_position_outside(combined, chrom_lengths, rng)

        pos_1based = pos_0based + 1
        ref = generate_random_nuc(rng)

        if rng.random() < ref_eq_alt_rate:
            alt = ref
        else:
            alt = generate_alt(ref, rng)

        qual = f"{rng.uniform(10.0, 999.0):.2f}"

        if rng.random() < low_qual_rate:
            filt = "LowQual"
        else:
            filt = "PASS"

        info_vals = generate_info(rng)
        info_str = ";".join(f"{k}={v}" for k, v in info_vals.items())

        gt = rng.choices(gt_choices, weights=gt_weights, k=1)[0]
        dp = rng.randint(5, 100)
        gq = rng.randint(10, 99)
        sample = f"{gt}:{dp}:{gq}"

        record = (
            f"{chrom}\t{pos_1based}\t.\t{ref}\t{alt}\t{qual}\t{filt}\t{info_str}\t"
            f"GT:DP:GQ\t{sample}"
        )
        records.append(record)

    return records


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    """Return a sort key for chromosome names."""
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom.isdigit():
        return (0, int(chrom))
    if chrom == "X":
        return (1, "")
    if chrom == "Y":
        return (2, "")
    if chrom in ("M", "MT"):
        return (3, "")
    return (4, chrom)


def write_vcf(records: list[str], header: list[str], output_path: Path) -> None:
    """Write VCF records to a bgzip-compressed file.

    Records are sorted by chromosome and position to ensure valid VCF.
    """
    # Parse and sort records
    parsed: list[tuple[str, int, str]] = []
    for record in records:
        parts = record.split("\t")
        chrom = parts[0]
        pos = int(parts[1])
        parsed.append((chrom, pos, record))

    parsed.sort(key=lambda x: (chrom_sort_key(x[0]), x[1]))

    tmp_path = output_path.with_suffix(".tmp.vcf")
    with open(tmp_path, "w") as f:
        for line in header:
            f.write(line + "\n")
        for _, _, record in parsed:
            f.write(record + "\n")

    # Compress with bgzip
    subprocess.run([BGZIP, "-f", str(tmp_path)], check=True)
    compressed = Path(str(tmp_path) + ".gz")
    compressed.rename(output_path)

    # Create tabix index
    subprocess.run([BCFTOOLS, "index", "-t", str(output_path)], check=True)


def generate_full_grch38(
    output_dir: Path,
    gene_intervals: dict[str, list[tuple[int, int]]],
    reg_intervals: dict[str, list[tuple[int, int]]],
    rng: random.Random,
) -> Path:
    """Generate full_grch38.vcf.gz — raw caller output style."""
    path = output_dir / "full_grch38.vcf.gz"
    header = build_vcf_header(GRCH38_CHROM_LENGTHS, reference="GRCh38")
    records = generate_vcf_records(
        n_variants=5000,
        chrom_lengths=GRCH38_CHROM_LENGTHS,
        gene_intervals=gene_intervals,
        reg_intervals=reg_intervals,
        placement="mixed",
        gt_distribution={"0/0": 0.30, "0/1": 0.40, "1/1": 0.25, "./.": 0.05},
        ref_eq_alt_rate=0.01,
        low_qual_rate=0.10,
        rng=rng,
    )
    write_vcf(records, header, path)
    return path


def generate_full_grch37(
    output_dir: Path,
    gene_intervals_grch37: dict[str, list[tuple[int, int]]],
    reg_intervals_grch37: dict[str, list[tuple[int, int]]],
    rng: random.Random,
) -> Path:
    """Generate full_grch37.vcf.gz — raw caller output style, GRCh37 coords."""
    path = output_dir / "full_grch37.vcf.gz"
    header = build_vcf_header(GRCH37_CHROM_LENGTHS, reference="GRCh37")
    records = generate_vcf_records(
        n_variants=5000,
        chrom_lengths=GRCH37_CHROM_LENGTHS,
        gene_intervals=gene_intervals_grch37,
        reg_intervals=reg_intervals_grch37,
        placement="mixed",
        gt_distribution={"0/0": 0.30, "0/1": 0.40, "1/1": 0.25, "./.": 0.05},
        ref_eq_alt_rate=0.01,
        low_qual_rate=0.10,
        rng=rng,
    )
    write_vcf(records, header, path)
    return path


def generate_genotyped_grch38(
    output_dir: Path,
    gene_intervals: dict[str, list[tuple[int, int]]],
    reg_intervals: dict[str, list[tuple[int, int]]],
    rng: random.Random,
) -> Path:
    """Generate genotyped_grch38.vcf.gz — clean genotyped."""
    path = output_dir / "genotyped_grch38.vcf.gz"
    header = build_vcf_header(GRCH38_CHROM_LENGTHS, reference="GRCh38")
    records = generate_vcf_records(
        n_variants=2000,
        chrom_lengths=GRCH38_CHROM_LENGTHS,
        gene_intervals=gene_intervals,
        reg_intervals=reg_intervals,
        placement="gene_or_regulatory",
        gt_distribution={"0/1": 0.60, "1/1": 0.40},
        ref_eq_alt_rate=0.0,
        low_qual_rate=0.0,
        rng=rng,
    )
    write_vcf(records, header, path)
    return path


def generate_genotyped_grch37(
    output_dir: Path,
    gene_intervals_grch37: dict[str, list[tuple[int, int]]],
    reg_intervals_grch37: dict[str, list[tuple[int, int]]],
    rng: random.Random,
) -> Path:
    """Generate genotyped_grch37.vcf.gz — clean genotyped, GRCh37 coords."""
    path = output_dir / "genotyped_grch37.vcf.gz"
    header = build_vcf_header(GRCH37_CHROM_LENGTHS, reference="GRCh37")
    records = generate_vcf_records(
        n_variants=2000,
        chrom_lengths=GRCH37_CHROM_LENGTHS,
        gene_intervals=gene_intervals_grch37,
        reg_intervals=reg_intervals_grch37,
        placement="gene_or_regulatory",
        gt_distribution={"0/1": 0.60, "1/1": 0.40},
        ref_eq_alt_rate=0.0,
        low_qual_rate=0.0,
        rng=rng,
    )
    write_vcf(records, header, path)
    return path


def strip_chr_prefix(intervals: dict[str, list[tuple[int, int]]]) -> dict[str, list[tuple[int, int]]]:
    """Remove 'chr' prefix from chromosome names in intervals dict."""
    result: dict[str, list[tuple[int, int]]] = {}
    for chrom, ivs in intervals.items():
        new_chrom = chrom[3:] if chrom.startswith("chr") else chrom
        result[new_chrom] = ivs
    return result


def generate_all_demo_vcfs(output_dir: Path, ref_dir: Path | None = None, seed: int = 42) -> dict[str, Path]:
    """Generate all 4 demo VCFs and return a mapping of names to paths.

    This function is used by the e2e test suite.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if ref_dir is None:
        ref_dir = Path("~/.dgra-prefilter/refs").expanduser()
    else:
        ref_dir = Path(ref_dir).expanduser()

    rng = random.Random(seed)

    gene_intervals = load_bed_intervals(ref_dir / "gencode_v44_gene_loci.bed")
    reg_intervals = load_bed_intervals(ref_dir / "encode_screen_v3_ccres.bed")

    gene_intervals_grch37 = strip_chr_prefix(gene_intervals)
    reg_intervals_grch37 = strip_chr_prefix(reg_intervals)

    paths: dict[str, Path] = {
        "full_grch38": generate_full_grch38(output_dir, gene_intervals, reg_intervals, rng),
        "full_grch37": generate_full_grch37(output_dir, gene_intervals_grch37, reg_intervals_grch37, rng),
        "genotyped_grch38": generate_genotyped_grch38(output_dir, gene_intervals, reg_intervals, rng),
        "genotyped_grch37": generate_genotyped_grch37(output_dir, gene_intervals_grch37, reg_intervals_grch37, rng),
    }
    return paths


def main() -> None:
    """CLI entry point for demo VCF generation."""
    parser = argparse.ArgumentParser(description="Generate synthetic demo VCFs.")
    parser.add_argument(
        "-o",
        "--output-dir",
        default="tests/demo_vcfs",
        help="Output directory for demo VCFs",
    )
    parser.add_argument(
        "--ref-dir",
        default="~/.dgra-prefilter/refs",
        help="Reference BED files directory",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility",
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    ref_dir = Path(args.ref_dir).expanduser()

    paths = generate_all_demo_vcfs(output_dir, ref_dir=ref_dir, seed=args.seed)

    for name, path in paths.items():
        print(f"Generated {name}: {path}")

    print(f"All demo VCFs generated in {output_dir}")


if __name__ == "__main__":
    main()
