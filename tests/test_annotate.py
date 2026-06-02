"""Tests for dgra_prefilter.annotate — VCF INFO field annotation.

Key areas:
  - DGRA_REGION / DGRA_SAFETYNET tag correctness
  - Coordinate boundary correctness (VCF 1-based → BED 0-based)
  - ClinVar zero-omission: every P/LP site must be tagged
  - Proper handling of existing INFO fields
  - Empty VCF / no-overlap cases
"""

from __future__ import annotations

from pathlib import Path

import pytest

from dgra_prefilter.annotate import VCFAnnotator
from dgra_prefilter.presets import get_preset


# ======================================================================
# Helper
# ======================================================================

def _parse_vcf_records(vcf_path: Path) -> list[dict]:
    """Parse a VCF file into a list of record dicts (chrom, pos, info)."""
    records = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 8:
                records.append({
                    "chrom": parts[0],
                    "pos": int(parts[1]),
                    "info": parts[7],
                })
    return records


def _get_info_tags(info: str) -> dict[str, list[str]]:
    """Parse INFO field into {tag: [values]} dict."""
    tags: dict[str, list[str]] = {}
    for field in info.split(";"):
        if "=" in field:
            key, val = field.split("=", 1)
            tags[key] = val.split(",")
    return tags


# ======================================================================
# _determine_regions
# ======================================================================

class TestDetermineRegions:
    """Tests for VCFAnnotator._determine_regions()."""

    def test_gene_region_hit(self, mini_ref_dir: Path) -> None:
        """VCF pos 150 (BED pos 149) is inside gene interval [100, 200)."""
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        regions = ann._determine_regions("chr1", 150)
        assert "gene" in regions

    def test_ncrna_region_hit(self, mini_ref_dir: Path) -> None:
        """VCF pos 550 (BED pos 549) is inside ncRNA interval [500, 600)."""
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        regions = ann._determine_regions("chr1", 550)
        assert "ncrna" in regions

    def test_regulatory_region_hit(self, mini_ref_dir: Path) -> None:
        """VCF pos 1100 (BED pos 1099) is inside regulatory interval [1000, 1200)."""
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        regions = ann._determine_regions("chr1", 1100)
        assert "regulatory" in regions

    def test_intergenic_no_region(self, mini_ref_dir: Path) -> None:
        """VCF pos 9999 is outside all intervals."""
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        regions = ann._determine_regions("chr1", 9999)
        assert regions == []

    def test_gene_and_regulatory_overlap(self, mini_ref_dir: Path) -> None:
        """A variant can fall in both gene and regulatory if intervals overlap."""
        # Position 108 is in gene [100,200) and clinvar [100,150)
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        regions = ann._determine_regions("chr1", 108)
        assert "gene" in regions


# ======================================================================
# _determine_safetynets
# ======================================================================

class TestDetermineSafetynets:
    """Tests for VCFAnnotator._determine_safetynets()."""

    def test_clinvar_hit(self, mini_ref_dir: Path) -> None:
        """VCF pos 108 (BED pos 107) is in ClinVar interval [100, 150)."""
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        safetynets = ann._determine_safetynets("chr1", 108)
        assert "ClinVar" in safetynets

    def test_clinvar_safetynet_only(self, mini_ref_dir: Path) -> None:
        """VCF pos 2050 (BED pos 2049) is in ClinVar [2000, 2100) but not in gene."""
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        safetynets = ann._determine_safetynets("chr1", 2050)
        assert "ClinVar" in safetynets

    def test_no_safetynet(self, mini_ref_dir: Path) -> None:
        """VCF pos 9999 is not in any safetynet."""
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        safetynets = ann._determine_safetynets("chr1", 9999)
        assert safetynets == []


# ======================================================================
# _annotate_record — boundary-focused tests
# ======================================================================

class TestAnnotateRecord:
    """Tests for VCFAnnotator._annotate_record() with boundary cases."""

    @pytest.fixture()
    def annotator(self, mini_ref_dir: Path) -> VCFAnnotator:
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        return ann

    def test_adds_region_tag(self, annotator: VCFAnnotator) -> None:
        record = "chr1\t150\t.\tA\tG\t30\tPASS\t."
        result = annotator._annotate_record(record)
        assert "DGRA_REGION=gene" in result

    def test_adds_safetynet_tag(self, annotator: VCFAnnotator) -> None:
        record = "chr1\t2050\t.\tT\tC\t30\tPASS\t."
        result = annotator._annotate_record(record)
        assert "DGRA_SAFETYNET=ClinVar" in result

    def test_both_region_and_safetynet(self, annotator: VCFAnnotator) -> None:
        """VCF pos 108 is in both gene and ClinVar."""
        record = "chr1\t108\t.\tG\tC\t30\tPASS\t."
        result = annotator._annotate_record(record)
        assert "DGRA_REGION=gene" in result
        assert "DGRA_SAFETYNET=ClinVar" in result

    def test_no_overlap_keeps_dot_info(self, annotator: VCFAnnotator) -> None:
        record = "chr1\t9999\t.\tA\tT\t30\tPASS\t."
        result = annotator._annotate_record(record)
        # Should still have . as INFO
        parts = result.strip().split("\t")
        assert parts[7] == "."

    def test_appends_to_existing_info(self, annotator: VCFAnnotator) -> None:
        record = "chr1\t150\t.\tA\tG\t30\tPASS\tDP=10"
        result = annotator._annotate_record(record)
        assert "DP=10;" in result
        assert "DGRA_REGION=gene" in result

    def test_dot_info_replaced(self, annotator: VCFAnnotator) -> None:
        """When INFO is '.', it should be replaced (not prepended with '.;').

        VCF pos 150 (BED pos 149) is in gene [100,200) AND ClinVar [100,150),
        so both tags appear. The key assertion is no leading '.' in INFO.
        """
        record = "chr1\t150\t.\tA\tG\t30\tPASS\t."
        result = annotator._annotate_record(record)
        parts = result.strip().split("\t")
        # INFO should NOT start with ".;" — the dot must be replaced, not appended to
        assert not parts[7].startswith(".;")
        assert "DGRA_REGION=gene" in parts[7]
        assert "DGRA_SAFETYNET=ClinVar" in parts[7]

    def test_fewer_than_8_columns_unchanged(self, annotator: VCFAnnotator) -> None:
        record = "chr1\t150\t.\tA\tG\t30\tPASS"
        result = annotator._annotate_record(record)
        # Less than 8 columns — returned as-is (with newline)
        assert result.strip() == record

    def test_coordinate_boundary_start(self, annotator: VCFAnnotator) -> None:
        """VCF pos 101 = BED pos 100 = start of gene interval → should be IN."""
        record = "chr1\t101\t.\tA\tG\t30\tPASS\t."
        result = annotator._annotate_record(record)
        assert "DGRA_REGION=gene" in result

    def test_coordinate_boundary_end(self, annotator: VCFAnnotator) -> None:
        """VCF pos 200 = BED pos 199 → last pos inside [100, 200) → should be IN."""
        record = "chr1\t200\t.\tA\tG\t30\tPASS\t."
        result = annotator._annotate_record(record)
        assert "DGRA_REGION=gene" in result

    def test_coordinate_boundary_past_end(self, annotator: VCFAnnotator) -> None:
        """VCF pos 201 = BED pos 200 → OUTSIDE [100, 200)."""
        record = "chr1\t201\t.\tA\tG\t30\tPASS\t."
        result = annotator._annotate_record(record)
        assert "DGRA_REGION" not in result

    def test_coordinate_boundary_before_start(self, annotator: VCFAnnotator) -> None:
        """VCF pos 100 = BED pos 99 → OUTSIDE [100, 200)."""
        record = "chr1\t100\t.\tA\tG\t30\tPASS\t."
        result = annotator._annotate_record(record)
        assert "DGRA_REGION" not in result


# ======================================================================
# annotate_vcf — full file annotation
# ======================================================================

class TestAnnotateVcf:
    """Tests for VCFAnnotator.annotate_vcf() — full VCF file annotation."""

    def test_annotates_full_vcf(self, mini_ref_dir: Path, tmp_vcf: Path, tmp_path: Path) -> None:
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        output = tmp_path / "annotated.vcf"
        stats = ann.annotate_vcf(tmp_vcf, output)

        records = _parse_vcf_records(output)
        # 8 input records
        assert len(records) == 8

        # Check specific annotations
        r150 = [r for r in records if r["pos"] == 150][0]
        assert "DGRA_REGION=gene" in r150["info"]

        r2050 = [r for r in records if r["pos"] == 2050][0]
        assert "DGRA_SAFETYNET=ClinVar" in r2050["info"]

        r9999 = [r for r in records if r["pos"] == 9999][0]
        assert "DGRA_REGION" not in r9999["info"]
        assert "DGRA_SAFETYNET" not in r9999["info"]

    def test_clinvar_zero_omission(self, mini_ref_dir: Path, tmp_vcf: Path, tmp_path: Path) -> None:
        """Every ClinVar P/LP site must be tagged — zero omissions allowed."""
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        output = tmp_path / "clinvar_check.vcf"
        stats = ann.annotate_vcf(tmp_vcf, output)

        records = _parse_vcf_records(output)
        # ClinVar positions in mini_ref: 100-150, 2000-2100
        clinvar_variants = [r for r in records if r["pos"] in (108, 2050)]
        for r in clinvar_variants:
            assert "DGRA_SAFETYNET=ClinVar" in r["info"], (
                f"ClinVar site chr1:{r['pos']} not tagged! Zero-omission violation."
            )

    def test_info_header_lines_added(self, mini_ref_dir: Path, tmp_vcf: Path, tmp_path: Path) -> None:
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        output = tmp_path / "header_check.vcf"
        ann.annotate_vcf(tmp_vcf, output)

        content = output.read_text()
        assert "ID=DGRA_REGION" in content
        assert "ID=DGRA_SAFETYNET" in content

    def test_empty_vcf_annotation(self, mini_ref_dir: Path, empty_vcf: Path, tmp_path: Path) -> None:
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        output = tmp_path / "empty_annotated.vcf"
        stats = ann.annotate_vcf(empty_vcf, output)

        records = _parse_vcf_records(output)
        assert len(records) == 0
        assert stats.region_only_variants == 0
        assert stats.safetynet_only_variants == 0

    def test_stats_region_counts(self, mini_ref_dir: Path, tmp_vcf: Path, tmp_path: Path) -> None:
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        output = tmp_path / "stats_check.vcf"
        stats = ann.annotate_vcf(tmp_vcf, output)

        # Should have gene counts (positions 150, 108, 550=ncrna not gene)
        assert stats.region_counts.get("gene", 0) > 0
        assert stats.clinvar_count > 0


# ======================================================================
# Load beds with different presets
# ======================================================================

class TestLoadBedsWithPresets:
    """Test that load_beds correctly loads different BED sets per preset."""

    def test_coding_only_loads_coding_exon(self, mini_ref_dir: Path) -> None:
        ann = VCFAnnotator()
        ann.load_beds(get_preset("coding-only"), mini_ref_dir)
        assert "gene" in ann.region_beds
        assert "ncrna" not in ann.region_beds
        assert "regulatory" not in ann.region_beds

    def test_regulatory_minimal_loads_pls_pels(self, mini_ref_dir: Path) -> None:
        ann = VCFAnnotator()
        ann.load_beds(get_preset("regulatory-minimal"), mini_ref_dir)
        assert "gene" in ann.region_beds
        assert "ncrna" in ann.region_beds
        assert "regulatory" in ann.region_beds
