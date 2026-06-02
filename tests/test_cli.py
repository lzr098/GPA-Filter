"""Tests for dgra_prefilter.cli — CLI argument parsing and execution."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from dgra_prefilter.cli import _args_to_config, _build_parser


# ======================================================================
# Argument parser
# ======================================================================

class TestBuildParser:
    """Tests for _build_parser()."""

    def test_parser_returns_argument_parser(self) -> None:
        parser = _build_parser()
        assert parser.prog == "dgra-prefilter"

    def test_required_args(self) -> None:
        parser = _build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])  # Missing -i and -o

    def test_minimal_args(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["-i", "input.vcf", "-o", "output.vcf"])
        assert args.input == Path("input.vcf")
        assert args.output == Path("output.vcf")

    def test_all_args(self) -> None:
        parser = _build_parser()
        args = parser.parse_args([
            "-i", "input.vcf",
            "-o", "output.vcf",
            "-g", "GRCh38",
            "-p", "coding-only",
            "--ref-dir", "/data/refs",
            "--report", "report.json",
            "--force",
            "--update-refs",
            "-v",
        ])
        assert args.genome == "GRCh38"
        assert args.preset == "coding-only"
        assert args.ref_dir == Path("/data/refs")
        assert args.report == Path("report.json")
        assert args.force is True
        assert args.update_refs is True
        assert args.verbose is True

    def test_default_preset(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["-i", "in.vcf", "-o", "out.vcf"])
        assert args.preset == "comprehensive"

    def test_default_genome(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["-i", "in.vcf", "-o", "out.vcf"])
        assert args.genome == "GRCh38"

    def test_default_force_false(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["-i", "in.vcf", "-o", "out.vcf"])
        assert args.force is False

    def test_default_verbose_false(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["-i", "in.vcf", "-o", "out.vcf"])
        assert args.verbose is False

    def test_invalid_preset_exits(self) -> None:
        parser = _build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([
                "-i", "in.vcf", "-o", "out.vcf",
                "-p", "nonexistent",
            ])

    def test_invalid_genome_exits(self) -> None:
        parser = _build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([
                "-i", "in.vcf", "-o", "out.vcf",
                "-g", "GRCh37",
            ])


# ======================================================================
# _args_to_config
# ======================================================================

class TestArgsToConfig:
    """Tests for _args_to_config()."""

    def test_maps_all_fields(self) -> None:
        parser = _build_parser()
        args = parser.parse_args([
            "-i", "input.vcf",
            "-o", "output.vcf",
            "-p", "coding-only",
            "--force",
        ])
        config_dict = _args_to_config(args)

        assert config_dict["input_path"] == Path("input.vcf")
        assert config_dict["output_path"] == Path("output.vcf")
        assert config_dict["preset_name"] == "coding-only"
        assert config_dict["force"] is True
        assert config_dict["update_refs"] is False
        assert config_dict["genome"] == "GRCh38"

    def test_report_none_by_default(self) -> None:
        parser = _build_parser()
        args = parser.parse_args(["-i", "in.vcf", "-o", "out.vcf"])
        config_dict = _args_to_config(args)
        assert config_dict["report_path"] is None
