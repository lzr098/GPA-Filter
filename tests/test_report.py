"""Tests for dgra_prefilter.report — JSON report generation."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dgra_prefilter.core import FilterStats, PrefilterConfig
from dgra_prefilter.report import ReportGenerator


class TestReportGenerator:
    """Tests for ReportGenerator.generate()."""

    @pytest.fixture()
    def stats(self) -> FilterStats:
        """Sample FilterStats for testing."""
        return FilterStats(
            input_variants=1000,
            retained_variants=250,
            region_only_variants=200,
            safetynet_only_variants=30,
            region_and_safetynet_variants=20,
            region_counts={"gene": 150, "ncrna": 30, "regulatory": 40},
            clinvar_count=50,
            omim_count=0,
            elapsed_seconds=5.123,
            ref_data_versions={"gencode_v44_gene_loci.bed": "v44"},
            preset_config={
                "name": "comprehensive",
                "gene": "True",
                "ncrna": "True",
            },
        )

    @pytest.fixture()
    def config(self, tmp_path: Path) -> PrefilterConfig:
        return PrefilterConfig(
            input_path=tmp_path / "input.vcf.gz",
            output_path=tmp_path / "output.vcf.gz",
            genome="GRCh38",
            preset_name="comprehensive",
            ref_dir=tmp_path / "refs",
        )

    def test_generates_valid_json(self, stats: FilterStats, config: PrefilterConfig, tmp_path: Path) -> None:
        output = tmp_path / "report.json"
        ReportGenerator.generate(stats, config, output)
        data = json.loads(output.read_text())
        assert isinstance(data, dict)

    def test_report_structure(self, stats: FilterStats, config: PrefilterConfig, tmp_path: Path) -> None:
        output = tmp_path / "report.json"
        ReportGenerator.generate(stats, config, output)
        data = json.loads(output.read_text())

        # Top-level keys
        assert "tool" in data
        assert "version" in data
        assert "timestamp" in data
        assert "input" in data
        assert "output" in data
        assert "regions" in data
        assert "safety_net" in data
        assert "preset" in data
        assert "ref_data_versions" in data
        assert "elapsed_seconds" in data

    def test_tool_name_and_version(self, stats: FilterStats, config: PrefilterConfig, tmp_path: Path) -> None:
        output = tmp_path / "report.json"
        ReportGenerator.generate(stats, config, output)
        data = json.loads(output.read_text())

        assert data["tool"] == "dgra-prefilter"
        assert data["version"] == "1.0.0"

    def test_input_section(self, stats: FilterStats, config: PrefilterConfig, tmp_path: Path) -> None:
        output = tmp_path / "report.json"
        ReportGenerator.generate(stats, config, output)
        data = json.loads(output.read_text())

        assert data["input"]["variants"] == 1000
        assert data["input"]["genome"] == "GRCh38"
        assert "path" in data["input"]

    def test_output_section(self, stats: FilterStats, config: PrefilterConfig, tmp_path: Path) -> None:
        output = tmp_path / "report.json"
        ReportGenerator.generate(stats, config, output)
        data = json.loads(output.read_text())

        assert data["output"]["variants"] == 250
        assert data["output"]["retention_rate"] == 0.25

    def test_regions_section(self, stats: FilterStats, config: PrefilterConfig, tmp_path: Path) -> None:
        output = tmp_path / "report.json"
        ReportGenerator.generate(stats, config, output)
        data = json.loads(output.read_text())

        assert data["regions"]["gene"] == 150
        assert data["regions"]["ncrna"] == 30
        assert data["regions"]["regulatory"] == 40
        assert data["regions"]["region_only"] == 200
        assert data["regions"]["safetynet_only"] == 30
        assert data["regions"]["region_and_safetynet"] == 20

    def test_safety_net_section(self, stats: FilterStats, config: PrefilterConfig, tmp_path: Path) -> None:
        output = tmp_path / "report.json"
        ReportGenerator.generate(stats, config, output)
        data = json.loads(output.read_text())

        assert data["safety_net"]["ClinVar"] == 50
        assert data["safety_net"]["OMIM"] == 0

    def test_creates_parent_dirs(self, stats: FilterStats, config: PrefilterConfig, tmp_path: Path) -> None:
        output = tmp_path / "deep" / "nested" / "report.json"
        ReportGenerator.generate(stats, config, output)
        assert output.exists()

    def test_timestamp_is_iso_format(self, stats: FilterStats, config: PrefilterConfig, tmp_path: Path) -> None:
        output = tmp_path / "report.json"
        ReportGenerator.generate(stats, config, output)
        data = json.loads(output.read_text())

        ts = data["timestamp"]
        # Should be parseable as ISO format
        assert "T" in ts

    def test_ref_data_versions(self, stats: FilterStats, config: PrefilterConfig, tmp_path: Path) -> None:
        output = tmp_path / "report.json"
        ReportGenerator.generate(stats, config, output)
        data = json.loads(output.read_text())

        assert "gencode_v44_gene_loci.bed" in data["ref_data_versions"]

    def test_empty_stats(self, tmp_path: Path) -> None:
        stats = FilterStats()
        config = PrefilterConfig(
            input_path=tmp_path / "input.vcf.gz",
            output_path=tmp_path / "output.vcf.gz",
        )
        output = tmp_path / "empty_report.json"
        ReportGenerator.generate(stats, config, output)
        data = json.loads(output.read_text())

        assert data["input"]["variants"] == 0
        assert data["output"]["retention_rate"] == 0.0
