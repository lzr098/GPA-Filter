"""Benchmark/comparison tests for preset retention rates."""

from __future__ import annotations

from pathlib import Path

import pytest

from dgra_prefilter.bed_utils import BedUtils
from dgra_prefilter.presets import get_preset


# ======================================================================
# regulatory-balanced vs regulatory-minimal retention benchmark
# ======================================================================

class TestRegulatoryBalancedVsMinimal:
    """Benchmark verifying retention rate difference between
    regulatory-balanced and regulatory-minimal presets.

    regulatory-balanced includes dELS and CTCF in addition to PLS/pELS,
    so it should retain strictly more variants than regulatory-minimal
    when those additional categories are present.
    """

    @pytest.fixture()
    def benchmark_ref_dir(self, tmp_path: Path) -> Path:
        """Reference dir where regulatory-balanced has extra intervals."""
        ref_dir = tmp_path / "refs"
        ref_dir.mkdir()

        def _write_bed(path: Path, intervals: list[tuple[str, int, int]]) -> None:
            with open(path, "w") as f:
                for chrom, start, end in intervals:
                    f.write(f"{chrom}\t{start}\t{end}\n")

        def _write_named(path: Path, intervals: list[tuple[str, int, int, str]]) -> None:
            with open(path, "w") as f:
                for chrom, start, end, name in intervals:
                    f.write(f"{chrom}\t{start}\t{end}\t{name}\n")

        _write_bed(ref_dir / "gencode_v44_gene_loci.bed", [
            ("chr1", 0, 10),
        ])
        _write_bed(ref_dir / "gencode_v44_ncrna_loci.bed", [
            ("chr1", 20, 30),
        ])
        # regulatory-minimal uses pls_pels (narrower)
        _write_bed(ref_dir / "encode_screen_v3_pls_pels.bed", [
            ("chr1", 40, 50),
            ("chr1", 50, 60),
        ])
        # regulatory-balanced adds dELS and CTCF regions
        _write_named(ref_dir / "encode_screen_v3_balanced.bed", [
            ("chr1", 40, 50, "PLS"),
            ("chr1", 50, 60, "pELS"),
            ("chr1", 70, 80, "dELS"),
            ("chr1", 90, 100, "CTCF"),
        ])
        _write_bed(ref_dir / "clinvar_pathogenic_GRCh38.bed", [])
        _write_bed(ref_dir / "omim_pathogenic_GRCh38.bed", [])

        return ref_dir

    def test_balanced_covers_more_regulatory_regions(
        self, benchmark_ref_dir: Path
    ) -> None:
        """Balanced preset merged BED must be a strict superset of minimal."""
        minimal_beds = get_preset("regulatory-minimal").get_region_bed_names()
        balanced_beds = get_preset("regulatory-balanced").get_region_bed_names()

        # Both presets use gene and ncRNA BEDs
        assert "gencode_v44_gene_loci.bed" in minimal_beds
        assert "gencode_v44_gene_loci.bed" in balanced_beds
        assert "gencode_v44_ncrna_loci.bed" in minimal_beds
        assert "gencode_v44_ncrna_loci.bed" in balanced_beds

        # regulatory-minimal uses pls_pels; balanced uses balanced
        assert "encode_screen_v3_pls_pels.bed" in minimal_beds
        assert "encode_screen_v3_balanced.bed" in balanced_beds

        # Load and merge both regulatory BEDs
        minimal_reg = BedUtils.load_bed(
            benchmark_ref_dir / "encode_screen_v3_pls_pels.bed"
        )
        balanced_reg = BedUtils.load_bed_named(
            benchmark_ref_dir / "encode_screen_v3_balanced.bed"
        )

        # balanced_reg must contain all minimal intervals (as raw coverage)
        for chrom, intervals in minimal_reg.items():
            for start, end in intervals:
                # At least one balanced interval must cover this point
                balanced_intervals = balanced_reg.get(chrom, [])
                covered = any(
                    b_start <= start < b_end and b_start < end <= b_end
                    for b_start, b_end, _ in balanced_intervals
                )
                assert covered, (
                    f"regulatory-balanced does not cover {chrom}:{start}-{end} "
                    "from regulatory-minimal"
                )

        # Total span of balanced must be greater than minimal
        minimal_span = sum(
            end - start
            for chroms in minimal_reg.values()
            for start, end in chroms
        )
        balanced_span = sum(
            end - start
            for chroms in balanced_reg.values()
            for start, end, _ in chroms
        )
        assert balanced_span > minimal_span, (
            "regulatory-balanced should cover more bases than regulatory-minimal"
        )
