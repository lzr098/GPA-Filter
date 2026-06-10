"""Tests for dgra_prefilter.ref_manager — Reference file management."""

from __future__ import annotations

from pathlib import Path

import pytest

from dgra_prefilter.presets import get_preset
from dgra_prefilter.ref_manager import RefManager


class TestRefManagerInit:
    """Tests for RefManager initialization."""

    def test_ref_dir_resolved(self, tmp_path: Path) -> None:
        rm = RefManager(tmp_path / "refs")
        assert rm.ref_dir.is_absolute()

    def test_tilde_expansion(self) -> None:
        rm = RefManager(Path("~/.dgra-prefilter/refs"))
        assert "~" not in str(rm.ref_dir)


class TestGetBedPath:
    """Tests for RefManager.get_bed_path()."""

    def test_returns_full_path(self, tmp_path: Path) -> None:
        rm = RefManager(tmp_path)
        result = rm.get_bed_path("gencode_v44_gene_loci.bed")
        assert result == tmp_path / "gencode_v44_gene_loci.bed"


class TestValidateRefs:
    """Tests for RefManager.validate_refs()."""

    def test_valid_refs_comprehensive(self, mini_ref_dir: Path) -> None:
        rm = RefManager(mini_ref_dir)
        preset = get_preset("comprehensive")
        assert rm.validate_refs(preset) is True

    def test_valid_refs_coding_only(self, mini_ref_dir: Path) -> None:
        rm = RefManager(mini_ref_dir)
        preset = get_preset("coding-only")
        assert rm.validate_refs(preset) is True

    def test_valid_refs_regulatory_minimal(self, mini_ref_dir: Path) -> None:
        rm = RefManager(mini_ref_dir)
        preset = get_preset("regulatory-minimal")
        assert rm.validate_refs(preset) is True

    def test_missing_region_bed(self, tmp_path: Path) -> None:
        rm = RefManager(tmp_path)
        preset = get_preset("comprehensive")
        with pytest.raises(FileNotFoundError, match="Required reference file missing"):
            rm.validate_refs(preset)

    def test_missing_safetynet_bed_warns(self, tmp_path: Path) -> None:
        """Missing safety net BED logs a warning but does not raise."""
        rm = RefManager(tmp_path)
        preset = get_preset("comprehensive")
        # Create region BEDs but not ClinVar/OMIM
        for name in preset.get_region_bed_names():
            (tmp_path / name).write_text("chr1\t100\t200\n")
        # Should not raise; warnings are logged
        assert rm.validate_refs(preset) is True


class TestMergePresetBeds:
    """Tests for RefManager.merge_preset_beds()."""

    def test_merge_comprehensive(self, mini_ref_dir: Path, tmp_path: Path) -> None:
        rm = RefManager(mini_ref_dir)
        preset = get_preset("comprehensive")
        output = tmp_path / "merged.bed"
        rm.merge_preset_beds(preset, output)
        assert output.exists()
        # Verify merged data has at least chr1
        from dgra_prefilter.bed_utils import BedUtils
        data = BedUtils.load_bed(output)
        assert "chr1" in data

    def test_merge_coding_only(self, mini_ref_dir: Path, tmp_path: Path) -> None:
        rm = RefManager(mini_ref_dir)
        preset = get_preset("coding-only")
        output = tmp_path / "merged.bed"
        rm.merge_preset_beds(preset, output)
        assert output.exists()


class TestGetVersions:
    """Tests for RefManager.get_versions()."""

    def test_reads_version_files(self, mini_ref_dir: Path) -> None:
        # Write a version file
        (mini_ref_dir / "gencode_v44_gene_loci.bed.version").write_text("v44\n")
        rm = RefManager(mini_ref_dir)
        versions = rm.get_versions()
        assert "gencode_v44_gene_loci.bed" in versions
        assert versions["gencode_v44_gene_loci.bed"] == "v44"

    def test_empty_version_file_skipped(self, mini_ref_dir: Path) -> None:
        (mini_ref_dir / "empty.bed.version").write_text("")
        rm = RefManager(mini_ref_dir)
        versions = rm.get_versions()
        assert "empty.bed" not in versions

    def test_nonexistent_dir(self, tmp_path: Path) -> None:
        rm = RefManager(tmp_path / "nonexistent")
        versions = rm.get_versions()
        assert versions == {}


class TestCheckFileNonempty:
    """Tests for RefManager._check_file_nonempty()."""

    def test_nonempty_file(self, tmp_path: Path) -> None:
        f = tmp_path / "test.txt"
        f.write_text("data\n")
        rm = RefManager(tmp_path)
        assert rm._check_file_nonempty(f) is True

    def test_empty_file(self, tmp_path: Path) -> None:
        f = tmp_path / "empty.txt"
        f.write_text("")
        rm = RefManager(tmp_path)
        assert rm._check_file_nonempty(f) is False

    def test_comments_only(self, tmp_path: Path) -> None:
        f = tmp_path / "comments.txt"
        f.write_text("# comment\ntrack name=x\n")
        rm = RefManager(tmp_path)
        assert rm._check_file_nonempty(f) is False

    def test_nonexistent_file(self, tmp_path: Path) -> None:
        rm = RefManager(tmp_path)
        assert rm._check_file_nonempty(tmp_path / "nope.txt") is False
