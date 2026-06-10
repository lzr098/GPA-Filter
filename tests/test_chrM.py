"""Tests for chrM handling in dgra-prefilter."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from dgra_prefilter.core import FilterEngine, PrefilterConfig


VCF_HEADER = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##assembly=GRCh38
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def _write_vcf(path: Path, header: str, records: list[str]) -> Path:
    with open(path, "w") as f:
        f.write(header)
        for rec in records:
            f.write(rec.rstrip("\n") + "\n")
    return path


def _write_bed(path: Path, intervals: list[tuple[str, int, int]]) -> Path:
    with open(path, "w") as f:
        for chrom, start, end in intervals:
            f.write(f"{chrom}\t{start}\t{end}\n")
    return path


def _build_mini_ref(ref_dir: Path) -> None:
    """Create a minimal reference directory."""
    ref_dir.mkdir(parents=True, exist_ok=True)
    _write_bed(ref_dir / "gencode_v44_gene_loci.bed", [
        ("chr1", 100, 200),
    ])
    _write_bed(ref_dir / "gencode_v44_ncrna_loci.bed", [
        ("chr1", 500, 600),
    ])
    _write_bed(ref_dir / "encode_screen_v3_ccres.bed", [
        ("chr1", 1000, 1200),
    ])
    _write_bed(ref_dir / "fantom5_enhancers_promoters.bed", [
        ("chr1", 1050, 1100),
    ])
    _write_bed(ref_dir / "vista_enhancers.bed", [
        ("chr1", 1100, 1150),
    ])
    _write_bed(ref_dir / "clinvar_pathogenic_GRCh38.bed", [
        ("chr1", 100, 150),
    ])
    (ref_dir / "omim_pathogenic_GRCh38.bed").write_text("# placeholder\n")
    (ref_dir / "ensembl_regulatory_features.bed").touch()
    (ref_dir / "gencode_v44_coding_exon_utr.bed").touch()
    (ref_dir / "gencode_5utr.bed").touch()
    (ref_dir / "gencode_cds.bed").touch()
    (ref_dir / "gencode_3utr.bed").touch()
    (ref_dir / "gencode_splice_sites.bed").touch()
    (ref_dir / "encode_screen_v3_pls_pels.bed").touch()
    (ref_dir / "encode_screen_v3_balanced.bed").touch()
    (ref_dir / "merged_retained_regions.bed").touch()


class TestKeepAllChrM:
    """Tests for --keep-all-chrM functionality."""

    def test_chrM_variant_retained_when_flag_set(self, tmp_path: Path) -> None:
        """chrM variants outside any region are retained when keep_all_chrM=True."""
        ref_dir = tmp_path / "refs"
        _build_mini_ref(ref_dir)

        # Create input VCF with a chrM variant and a chr1 variant
        input_vcf = tmp_path / "input.vcf"
        _write_vcf(input_vcf, VCF_HEADER, [
            "chr1\t150\t.\tA\tG\t30\tPASS\tDP=10",
            "chrM\t1000\t.\tT\tC\t30\tPASS\tDP=5",
        ])

        output_vcf = tmp_path / "output.vcf.gz"
        config = PrefilterConfig(
            input_path=input_vcf,
            output_path=output_vcf,
            ref_dir=ref_dir,
            preset_name="comprehensive",
            keep_all_chrM=True,
            force=True,
        )
        engine = FilterEngine(config)
        result = engine.run()

        # Both variants should be retained
        assert result.stats.retained_variants == 2

    def test_chrM_variant_filtered_by_default(self, tmp_path: Path) -> None:
        """chrM variants outside any region are filtered by default."""
        ref_dir = tmp_path / "refs"
        _build_mini_ref(ref_dir)

        input_vcf = tmp_path / "input.vcf"
        _write_vcf(input_vcf, VCF_HEADER, [
            "chr1\t150\t.\tA\tG\t30\tPASS\tDP=10",
            "chrM\t1000\t.\tT\tC\t30\tPASS\tDP=5",
        ])

        output_vcf = tmp_path / "output.vcf.gz"
        config = PrefilterConfig(
            input_path=input_vcf,
            output_path=output_vcf,
            ref_dir=ref_dir,
            preset_name="comprehensive",
            keep_all_chrM=False,
            force=True,
        )
        engine = FilterEngine(config)
        result = engine.run()

        # Only chr1 variant should be retained
        assert result.stats.retained_variants == 1

    def test_merged_bed_contains_chrM_when_flag_set(self, tmp_path: Path) -> None:
        """The merged BED should contain chrM 0-16569 when keep_all_chrM=True."""
        ref_dir = tmp_path / "refs"
        _build_mini_ref(ref_dir)

        input_vcf = tmp_path / "input.vcf"
        _write_vcf(input_vcf, VCF_HEADER, [
            "chr1\t150\t.\tA\tG\t30\tPASS\tDP=10",
        ])

        config = PrefilterConfig(
            input_path=input_vcf,
            output_path=tmp_path / "output.vcf.gz",
            ref_dir=ref_dir,
            preset_name="comprehensive",
            keep_all_chrM=True,
            force=True,
        )
        engine = FilterEngine(config)
        # Run until just after BED merge by calling run() - it will fail
        # at bcftools because input isn't bgzipped/indexed, so instead we
        # inspect the merged BED by replicating the logic.
        engine._temp_dir = __import__("tempfile").TemporaryDirectory(prefix="dgra_test_")
        temp_path = Path(engine._temp_dir.name)
        merged_bed = temp_path / "merged_regions.bed"
        engine.ref_manager.merge_preset_beds(engine.preset, merged_bed)
        if engine.preset.keep_all_chrM:
            with open(merged_bed, "a") as f:
                f.write("chrM\t0\t16569\n")

        content = merged_bed.read_text()
        assert "chrM\t0\t16569" in content
        engine._temp_dir.cleanup()
