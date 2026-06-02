"""Tests for dgra_prefilter.core — Core filtering engine (pure Python logic only).

Note: bcftools-dependent methods are tested via mocking where needed.
We focus on testing pure Python logic: PrefilterConfig, FilterStats,
exception classes, and the top-level API parameter mapping.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from dgra_prefilter.core import (
    BcftoolsNotFoundError,
    FilterEngine,
    FilterResult,
    FilterStats,
    GenomeMismatchError,
    PrefilterConfig,
    RefDataMissingError,
    VCFProcessingError,
)


# ======================================================================
# PrefilterConfig
# ======================================================================

class TestPrefilterConfig:
    """Tests for PrefilterConfig dataclass."""

    def test_default_values(self, tmp_path: Path) -> None:
        config = PrefilterConfig(
            input_path=tmp_path / "input.vcf",
            output_path=tmp_path / "output.vcf",
        )
        assert config.genome == "GRCh38"
        assert config.preset_name == "comprehensive"
        assert config.force is False
        assert config.update_refs is False

    def test_path_normalization(self, tmp_path: Path) -> None:
        config = PrefilterConfig(
            input_path="input.vcf",
            output_path="output.vcf",
            ref_dir="~/refs",
        )
        assert config.input_path.is_absolute()
        assert config.output_path.is_absolute()
        assert config.ref_dir.is_absolute()
        assert "~" not in str(config.ref_dir)

    def test_report_path_none_by_default(self, tmp_path: Path) -> None:
        config = PrefilterConfig(
            input_path=tmp_path / "input.vcf",
            output_path=tmp_path / "output.vcf",
        )
        assert config.report_path is None

    def test_report_path_normalized(self, tmp_path: Path) -> None:
        config = PrefilterConfig(
            input_path=tmp_path / "input.vcf",
            output_path=tmp_path / "output.vcf",
            report_path="report.json",
        )
        assert config.report_path is not None
        assert config.report_path.is_absolute()


# ======================================================================
# FilterStats
# ======================================================================

class TestFilterStats:
    """Tests for FilterStats dataclass."""

    def test_default_values(self) -> None:
        stats = FilterStats()
        assert stats.input_variants == 0
        assert stats.retained_variants == 0
        assert stats.region_counts == {}
        assert stats.clinvar_count == 0
        assert stats.omim_count == 0

    def test_retention_rate_zero_input(self) -> None:
        stats = FilterStats(input_variants=0, retained_variants=0)
        assert stats.retention_rate == 0.0

    def test_retention_rate_calculation(self) -> None:
        stats = FilterStats(input_variants=100, retained_variants=25)
        assert stats.retention_rate == 0.25

    def test_retention_rate_full(self) -> None:
        stats = FilterStats(input_variants=100, retained_variants=100)
        assert stats.retention_rate == 1.0

    def test_retention_rate_partial(self) -> None:
        stats = FilterStats(input_variants=1000, retained_variants=333)
        assert abs(stats.retention_rate - 0.333) < 0.01


# ======================================================================
# Exception classes
# ======================================================================

class TestExceptions:
    """Tests for custom exception classes."""

    def test_genome_mismatch_error(self) -> None:
        err = GenomeMismatchError("GRCh38", "GRCh37")
        assert err.expected == "GRCh38"
        assert err.found == "GRCh37"
        assert "GRCh38" in str(err)
        assert "GRCh37" in str(err)

    def test_bcftools_not_found_error(self) -> None:
        err = BcftoolsNotFoundError()
        assert "bcftools" in str(err).lower()

    def test_ref_data_missing_error(self) -> None:
        err = RefDataMissingError("Missing file")
        assert "Missing file" in str(err)

    def test_vcf_processing_error(self) -> None:
        err = VCFProcessingError(["bcftools", "view"], 1, "some error")
        assert err.cmd == ["bcftools", "view"]
        assert err.returncode == 1
        assert err.stderr == "some error"
        assert "exit code 1" in str(err)


# ======================================================================
# FilterResult
# ======================================================================

class TestFilterResult:
    """Tests for FilterResult dataclass."""

    def test_fields(self, tmp_path: Path) -> None:
        stats = FilterStats(input_variants=100)
        result = FilterResult(
            output_path=tmp_path / "output.vcf",
            report_path=tmp_path / "report.json",
            stats=stats,
        )
        assert result.stats.input_variants == 100
        assert result.output_path == tmp_path / "output.vcf"
        assert result.report_path == tmp_path / "report.json"


# ======================================================================
# FilterEngine — input validation (pure Python checks)
# ======================================================================

class TestFilterEngineValidation:
    """Tests for FilterEngine input validation that doesn't need bcftools."""

    def test_missing_input_file_raises(self, tmp_path: Path) -> None:
        config = PrefilterConfig(
            input_path=tmp_path / "nonexistent.vcf",
            output_path=tmp_path / "output.vcf",
        )
        engine = FilterEngine(config)
        with pytest.raises(FileNotFoundError, match="Input VCF not found"):
            engine._validate_input()
