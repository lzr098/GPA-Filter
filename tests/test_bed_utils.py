"""Tests for dgra_prefilter.bed_utils — highest priority module.

Key areas:
  - Coordinate conversion (BED 0-based half-open vs VCF 1-based)
  - Chromosome normalization (chr prefix)
  - Interval merging
  - Overlap query via bisect
  - BED file I/O
"""

from __future__ import annotations

from pathlib import Path

import pytest

from dgra_prefilter.bed_utils import BedUtils


# ======================================================================
# normalize_chrom
# ======================================================================

class TestNormalizeChrom:
    """Tests for BedUtils.normalize_chrom()."""

    def test_adds_chr_prefix(self) -> None:
        assert BedUtils.normalize_chrom("1") == "chr1"

    def test_already_has_chr_prefix(self) -> None:
        assert BedUtils.normalize_chrom("chr1") == "chr1"

    def test_x_chromosome(self) -> None:
        assert BedUtils.normalize_chrom("X") == "chrX"

    def test_chr_x_chromosome(self) -> None:
        assert BedUtils.normalize_chrom("chrX") == "chrX"

    def test_mitochondrial(self) -> None:
        assert BedUtils.normalize_chrom("M") == "chrM"

    def test_chr_mitochondrial(self) -> None:
        assert BedUtils.normalize_chrom("chrM") == "chrM"

    def test_y_chromosome(self) -> None:
        assert BedUtils.normalize_chrom("Y") == "chrY"

    def test_22(self) -> None:
        assert BedUtils.normalize_chrom("22") == "chr22"


# ======================================================================
# load_bed
# ======================================================================

class TestLoadBed:
    """Tests for BedUtils.load_bed()."""

    def test_load_simple_bed(self, simple_bed: Path) -> None:
        data = BedUtils.load_bed(simple_bed)
        assert "chr1" in data
        assert "chr2" in data
        assert len(data["chr1"]) == 2
        assert len(data["chr2"]) == 2

    def test_intervals_sorted_by_start(self, simple_bed: Path) -> None:
        data = BedUtils.load_bed(simple_bed)
        for chrom, intervals in data.items():
            starts = [iv[0] for iv in intervals]
            assert starts == sorted(starts)

    def test_normalizes_chr_prefix(self, tmp_path: Path) -> None:
        """BED file without 'chr' prefix should be normalized."""
        bed = tmp_path / "no_chr.bed"
        bed.write_text("1\t100\t200\n2\t300\t400\n")
        data = BedUtils.load_bed(bed)
        assert "chr1" in data
        assert "chr2" in data

    def test_skips_comment_lines(self, tmp_path: Path) -> None:
        bed = tmp_path / "comments.bed"
        bed.write_text("# comment line\nchr1\t100\t200\n")
        data = BedUtils.load_bed(bed)
        assert len(data) == 1
        assert data["chr1"] == [(100, 200)]

    def test_skips_track_header(self, tmp_path: Path) -> None:
        bed = tmp_path / "track.bed"
        bed.write_text('track name="test"\nchr1\t100\t200\n')
        data = BedUtils.load_bed(bed)
        assert data["chr1"] == [(100, 200)]

    def test_skips_blank_lines(self, tmp_path: Path) -> None:
        bed = tmp_path / "blanks.bed"
        bed.write_text("\n\nchr1\t100\t200\n\n")
        data = BedUtils.load_bed(bed)
        assert data["chr1"] == [(100, 200)]

    def test_skips_malformed_lines(self, tmp_path: Path) -> None:
        bed = tmp_path / "malformed.bed"
        bed.write_text("chr1\t100\t200\nbad_line\nchr2\t300\t400\n")
        data = BedUtils.load_bed(bed)
        assert "chr1" in data
        assert "chr2" in data
        assert len(data) == 2

    def test_file_not_found(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="BED file not found"):
            BedUtils.load_bed(tmp_path / "nonexistent.bed")

    def test_space_separated_fallback(self, tmp_path: Path) -> None:
        bed = tmp_path / "spaces.bed"
        bed.write_text("chr1 100 200\n")
        data = BedUtils.load_bed(bed)
        assert data["chr1"] == [(100, 200)]

    def test_invalid_coordinates_skipped(self, tmp_path: Path) -> None:
        bed = tmp_path / "bad_coords.bed"
        bed.write_text("chr1\tabc\tdef\nchr2\t300\t400\n")
        data = BedUtils.load_bed(bed)
        assert "chr1" not in data
        assert data["chr2"] == [(300, 400)]


# ======================================================================
# merge_intervals
# ======================================================================

class TestMergeIntervals:
    """Tests for BedUtils.merge_intervals()."""

    def test_empty_list(self) -> None:
        assert BedUtils.merge_intervals([]) == []

    def test_single_interval(self) -> None:
        assert BedUtils.merge_intervals([(100, 200)]) == [(100, 200)]

    def test_non_overlapping(self) -> None:
        result = BedUtils.merge_intervals([(100, 200), (300, 400)])
        assert result == [(100, 200), (300, 400)]

    def test_overlapping(self) -> None:
        result = BedUtils.merge_intervals([(100, 250), (200, 350)])
        assert result == [(100, 350)]

    def test_adjacent_merged(self) -> None:
        """Adjacent intervals (gap=0) should be merged."""
        result = BedUtils.merge_intervals([(100, 200), (200, 300)])
        assert result == [(100, 300)]

    def test_contained_interval(self) -> None:
        result = BedUtils.merge_intervals([(100, 400), (150, 250)])
        assert result == [(100, 400)]

    def test_multiple_merges(self) -> None:
        result = BedUtils.merge_intervals([
            (100, 250), (200, 350), (340, 400),
        ])
        assert result == [(100, 400)]

    def test_three_way_merge(self) -> None:
        result = BedUtils.merge_intervals([(10, 20), (15, 25), (20, 30)])
        assert result == [(10, 30)]


# ======================================================================
# query_overlap — CRITICAL: coordinate boundary tests
# ======================================================================

class TestQueryOverlap:
    """Tests for BedUtils.query_overlap() — VCF 1-based → BED 0-based conversion.

    BED uses 0-based half-open [start, end).
    VCF uses 1-based inclusive [pos, pos].

    For BED interval [100, 200):
      - VCF pos 101 → BED pos 100 → 100 < 200 ✓ IN
      - VCF pos 200 → BED pos 199 → 199 < 200 ✓ IN
      - VCF pos 201 → BED pos 200 → 200 < 200 ✗ OUT
      - VCF pos 100 → BED pos 99  → 99 < 100 ✗ OUT
    """

    @pytest.fixture()
    def bed_data(self) -> dict[str, list[tuple[int, int]]]:
        """Simple bed_data with chr1 interval [100, 200)."""
        return {"chr1": [(100, 200)]}

    def test_inside_interval(self, bed_data: dict) -> None:
        """VCF pos 101 (BED pos 100) is inside [100, 200)."""
        result = BedUtils.query_overlap(bed_data, "chr1", 101)
        assert len(result) == 1
        assert result[0] == (100, 200)

    def test_at_end_boundary(self, bed_data: dict) -> None:
        """VCF pos 200 (BED pos 199) is inside [100, 200) — last included position."""
        result = BedUtils.query_overlap(bed_data, "chr1", 200)
        assert len(result) == 1

    def test_just_past_end(self, bed_data: dict) -> None:
        """VCF pos 201 (BED pos 200) is OUTSIDE [100, 200)."""
        result = BedUtils.query_overlap(bed_data, "chr1", 201)
        assert result == []

    def test_just_before_start(self, bed_data: dict) -> None:
        """VCF pos 100 (BED pos 99) is OUTSIDE [100, 200)."""
        result = BedUtils.query_overlap(bed_data, "chr1", 100)
        assert result == []

    def test_at_start_boundary(self, bed_data: dict) -> None:
        """VCF pos 101 (BED pos 100) equals start — inside."""
        result = BedUtils.query_overlap(bed_data, "chr1", 101)
        assert len(result) == 1

    def test_chrom_not_in_data(self, bed_data: dict) -> None:
        result = BedUtils.query_overlap(bed_data, "chr5", 150)
        assert result == []

    def test_chrom_normalization_in_query(self, bed_data: dict) -> None:
        """Query with '1' should match 'chr1' after normalization."""
        result = BedUtils.query_overlap(bed_data, "1", 150)
        assert len(result) == 1

    def test_multiple_intervals(self) -> None:
        """Query should find the correct interval among multiple."""
        bed_data = {"chr1": [(100, 200), (300, 400), (500, 600)]}
        # Hit middle interval
        result = BedUtils.query_overlap(bed_data, "chr1", 350)
        assert result == [(300, 400)]
        # Hit first interval
        result2 = BedUtils.query_overlap(bed_data, "chr1", 150)
        assert result2 == [(100, 200)]

    def test_between_intervals(self) -> None:
        """Position between two intervals should return empty."""
        bed_data = {"chr1": [(100, 200), (300, 400)]}
        result = BedUtils.query_overlap(bed_data, "chr1", 250)
        assert result == []

    def test_empty_bed_data(self) -> None:
        result = BedUtils.query_overlap({}, "chr1", 150)
        assert result == []

    def test_empty_intervals_list(self) -> None:
        result = BedUtils.query_overlap({"chr1": []}, "chr1", 150)
        assert result == []

    def test_zero_position(self) -> None:
        """VCF position 0 is invalid but should not crash."""
        bed_data = {"chr1": [(0, 10)]}
        # VCF pos 1 → BED pos 0, inside [0, 10)
        result = BedUtils.query_overlap(bed_data, "chr1", 1)
        assert len(result) == 1


# ======================================================================
# merge_bed_files
# ======================================================================

class TestMergeBedFiles:
    """Tests for BedUtils.merge_bed_files()."""

    def test_merge_two_files(self, tmp_path: Path) -> None:
        bed1 = tmp_path / "a.bed"
        bed1.write_text("chr1\t100\t200\n")
        bed2 = tmp_path / "b.bed"
        bed2.write_text("chr1\t150\t250\n")

        output = tmp_path / "merged.bed"
        BedUtils.merge_bed_files([bed1, bed2], output)

        data = BedUtils.load_bed(output)
        assert data["chr1"] == [(100, 250)]

    def test_merge_non_overlapping(self, tmp_path: Path) -> None:
        bed1 = tmp_path / "a.bed"
        bed1.write_text("chr1\t100\t200\n")
        bed2 = tmp_path / "b.bed"
        bed2.write_text("chr2\t300\t400\n")

        output = tmp_path / "merged.bed"
        BedUtils.merge_bed_files([bed1, bed2], output)

        data = BedUtils.load_bed(output)
        assert data["chr1"] == [(100, 200)]
        assert data["chr2"] == [(300, 400)]

    def test_merge_missing_file(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            BedUtils.merge_bed_files(
                [tmp_path / "nonexistent.bed"],
                tmp_path / "out.bed",
            )


# ======================================================================
# write_bed
# ======================================================================

class TestWriteBed:
    """Tests for BedUtils.write_bed()."""

    def test_write_and_read_roundtrip(self, tmp_path: Path) -> None:
        intervals = {"chr1": [(100, 200)], "chr2": [(300, 400)]}
        output = tmp_path / "output.bed"
        BedUtils.write_bed(intervals, output)

        data = BedUtils.load_bed(output)
        assert data["chr1"] == [(100, 200)]
        assert data["chr2"] == [(300, 400)]

    def test_chromosome_sorting(self, tmp_path: Path) -> None:
        """Written chromosomes should be in natural order."""
        intervals = {
            "chr2": [(300, 400)],
            "chr10": [(100, 200)],
            "chr1": [(500, 600)],
            "chrX": [(700, 800)],
        }
        output = tmp_path / "sorted.bed"
        BedUtils.write_bed(intervals, output)

        lines = output.read_text().strip().split("\n")
        chroms = [l.split("\t")[0] for l in lines]
        assert chroms == ["chr1", "chr2", "chr10", "chrX"]

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        output = tmp_path / "deep" / "nested" / "out.bed"
        BedUtils.write_bed({"chr1": [(100, 200)]}, output)
        assert output.exists()

    def test_empty_intervals(self, tmp_path: Path) -> None:
        output = tmp_path / "empty.bed"
        BedUtils.write_bed({}, output)
        assert output.exists()
        assert output.read_text() == ""
