"""Tests for dgra_prefilter.annotate — VCF INFO field annotation.

Key areas:
  - DGRA_REGION / DGRA_SAFETYNET tag correctness
  - Coordinate boundary correctness (VCF 1-based → BED 0-based)
  - ClinVar zero-omission: every P/LP site must be tagged
  - ClinVar star-level tags (ClinVar_1star, ClinVar_3star)
  - UTR/CDS/splice sub-region tags
  - regulatory-balanced sub-tags
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

    def test_gene_sub_regions(self, mini_ref_dir: Path) -> None:
        """Sub-region BEDs produce gene_5utr, gene_cds, gene_3utr, gene_splice."""
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        assert "gene_5utr" in ann._determine_regions("chr1", 110)
        assert "gene_cds" in ann._determine_regions("chr1", 130)
        assert "gene_3utr" in ann._determine_regions("chr1", 170)
        assert "gene_splice" in ann._determine_regions("chr1", 182)

    def test_regulatory_balanced_sub_tags(self, mini_ref_dir: Path) -> None:
        """regulatory-balanced preset emits regulatory_pls/pels/dels/ctcf."""
        ann = VCFAnnotator()
        ann.load_beds(get_preset("regulatory-balanced"), mini_ref_dir)
        assert "regulatory_pls" in ann._determine_regions("chr1", 1020)
        assert "regulatory_pels" in ann._determine_regions("chr1", 1070)
        assert "regulatory_dels" in ann._determine_regions("chr1", 1120)
        assert "regulatory_ctcf" in ann._determine_regions("chr1", 1170)
        # Summary tag preserved
        assert "regulatory" in ann._determine_regions("chr1", 1020)


# ======================================================================
# _determine_safetynets
# ======================================================================

class TestDetermineSafetynets:
    """Tests for VCFAnnotator._determine_safetynets()."""

    def test_clinvar_hit(self, mini_ref_dir: Path) -> None:
        """VCF pos 108 (BED pos 107) is in ClinVar interval [100, 150) 3-star."""
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        safetynets = ann._determine_safetynets("chr1", 108)
        assert "ClinVar_3star" in safetynets

    def test_clinvar_safetynet_only(self, mini_ref_dir: Path) -> None:
        """VCF pos 2050 (BED pos 2049) is in ClinVar [2000, 2100) 1-star."""
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        safetynets = ann._determine_safetynets("chr1", 2050)
        assert "ClinVar_1star" in safetynets

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

    def test_adds_safetynet_star_tag(self, annotator: VCFAnnotator) -> None:
        record = "chr1\t2050\t.\tT\tC\t30\tPASS\t."
        result = annotator._annotate_record(record)
        assert "DGRA_SAFETYNET=ClinVar_1star" in result

    def test_both_region_and_safetynet(self, annotator: VCFAnnotator) -> None:
        """VCF pos 108 is in both gene and ClinVar 3-star."""
        record = "chr1\t108\t.\tG\tC\t30\tPASS\t."
        result = annotator._annotate_record(record)
        assert "DGRA_REGION=gene" in result
        assert "DGRA_SAFETYNET=ClinVar_3star" in result

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
        assert "DGRA_SAFETYNET=ClinVar_3star" in parts[7]

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
        assert "DGRA_SAFETYNET=ClinVar_1star" in r2050["info"]

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
        # ClinVar positions in mini_ref: 100-150 (3-star), 2000-2100 (1-star)
        clinvar_variants = [r for r in records if r["pos"] in (108, 2050)]
        for r in clinvar_variants:
            tags = _get_info_tags(r["info"]).get("DGRA_SAFETYNET", [])
            assert any(t.startswith("ClinVar_") for t in tags), (
                f"ClinVar site chr1:{r['pos']} not tagged! Zero-omission violation."
            )

    def test_clinvar_star_tags(self, mini_ref_dir: Path, tmp_vcf: Path, tmp_path: Path) -> None:
        """ClinVar annotations include star level (ClinVar_3star, ClinVar_1star)."""
        ann = VCFAnnotator()
        ann.load_beds(get_preset("comprehensive"), mini_ref_dir)
        output = tmp_path / "star_check.vcf"
        ann.annotate_vcf(tmp_vcf, output)

        records = _parse_vcf_records(output)
        r108 = [r for r in records if r["pos"] == 108][0]
        r2050 = [r for r in records if r["pos"] == 2050][0]

        assert "ClinVar_3star" in r108["info"]
        assert "ClinVar_1star" in r2050["info"]
        # Base ClinVar tag should not appear when star info is available
        assert "DGRA_SAFETYNET=ClinVar;" not in r108["info"]
        assert "DGRA_SAFETYNET=ClinVar\t" not in r108["info"]

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

    def test_sub_region_counts(self, mini_ref_dir: Path, tmp_path: Path) -> None:
        """Sub-region counts are tracked separately in annotation stats."""
        from tests.conftest import VCF_HEADER, write_vcf

        records = [
            "chr1\t110\t.\tA\tG\t30\tPASS\tDP=10",
            "chr1\t130\t.\tA\tG\t30\tPASS\tDP=10",
            "chr1\t170\t.\tA\tG\t30\tPASS\tDP=10",
            "chr1\t182\t.\tA\tG\t30\tPASS\tDP=10",
        ]
        vcf_path = tmp_path / "sub_regions.vcf"
        write_vcf(vcf_path, VCF_HEADER, records)

        ann = VCFAnnotator()
        ann.load_beds(get_preset("coding-only"), mini_ref_dir)
        output = tmp_path / "sub_annotated.vcf"
        stats = ann.annotate_vcf(vcf_path, output)

        assert stats.region_counts.get("gene_5utr", 0) == 1
        assert stats.region_counts.get("gene_cds", 0) == 1
        assert stats.region_counts.get("gene_3utr", 0) == 1
        assert stats.region_counts.get("gene_splice", 0) == 1
        # Summary gene count equals sum of sub-categories
        assert stats.region_counts.get("gene", 0) == 4

    def test_regulatory_balanced_annotation(self, mini_ref_dir: Path, tmp_path: Path) -> None:
        """regulatory-balanced outputs regulatory_pls/pels/dels/ctcf tags."""
        from tests.conftest import VCF_HEADER, write_vcf

        records = [
            "chr1\t1020\t.\tA\tG\t30\tPASS\tDP=10",
            "chr1\t1070\t.\tA\tG\t30\tPASS\tDP=10",
            "chr1\t1120\t.\tA\tG\t30\tPASS\tDP=10",
            "chr1\t1170\t.\tA\tG\t30\tPASS\tDP=10",
        ]
        vcf_path = tmp_path / "balanced.vcf"
        write_vcf(vcf_path, VCF_HEADER, records)

        ann = VCFAnnotator()
        ann.load_beds(get_preset("regulatory-balanced"), mini_ref_dir)
        output = tmp_path / "balanced_annotated.vcf"
        stats = ann.annotate_vcf(vcf_path, output)

        assert stats.region_counts.get("regulatory_pls", 0) == 1
        assert stats.region_counts.get("regulatory_pels", 0) == 1
        assert stats.region_counts.get("regulatory_dels", 0) == 1
        assert stats.region_counts.get("regulatory_ctcf", 0) == 1
        assert stats.region_counts.get("regulatory", 0) == 4


# ======================================================================
# Load beds with different presets
# ======================================================================

class TestLoadBedsWithPresets:
    """Test that load_beds correctly loads different BED sets per preset."""

    def test_coding_only_loads_sub_region_beds(self, mini_ref_dir: Path) -> None:
        ann = VCFAnnotator()
        ann.load_beds(get_preset("coding-only"), mini_ref_dir)
        # Summary gene bed loaded for backward-compatible annotation
        assert "gene" in ann.region_beds
        # Sub-region beds loaded
        assert "gene_5utr" in ann.region_simple_sub_beds
        assert "gene_cds" in ann.region_simple_sub_beds
        assert "gene_3utr" in ann.region_simple_sub_beds
        assert "gene_splice" in ann.region_simple_sub_beds
        assert "ncrna" not in ann.region_beds
        assert "regulatory" not in ann.region_beds

    def test_regulatory_minimal_loads_pls_pels(self, mini_ref_dir: Path) -> None:
        ann = VCFAnnotator()
        ann.load_beds(get_preset("regulatory-minimal"), mini_ref_dir)
        assert "gene" in ann.region_beds
        assert "ncrna" in ann.region_beds
        assert "regulatory" in ann.region_beds

    def test_regulatory_balanced_loads_balanced_bed(self, mini_ref_dir: Path) -> None:
        ann = VCFAnnotator()
        ann.load_beds(get_preset("regulatory-balanced"), mini_ref_dir)
        assert "gene" in ann.region_beds
        assert "ncrna" in ann.region_beds
        assert "regulatory" in ann.region_beds
        assert "regulatory" in ann.region_named_sub_beds


# ======================================================================
# Ensembl regulatory annotations
# ======================================================================

class TestEnsemblRegulatory:
    """Tests for Ensembl Regulatory Build annotation."""

    def test_ensembl_sub_tags(self, mini_ref_dir: Path) -> None:
        """Ensembl regulatory_source emits promoter/enhancer/ctcf/open_chromatin/tf_binding."""
        from dataclasses import replace

        preset = replace(get_preset("comprehensive"), regulatory_source="ensembl")
        ann = VCFAnnotator()
        ann.load_beds(preset, mini_ref_dir)
        assert "regulatory_promoter" in ann._determine_regions("chr1", 1250)
        assert "regulatory_enhancer" in ann._determine_regions("chr1", 1350)
        assert "regulatory_ctcf" in ann._determine_regions("chr1", 1425)
        assert "regulatory_open_chromatin" in ann._determine_regions("chr1", 1475)
        assert "regulatory_tf_binding" in ann._determine_regions("chr1", 1525)
        # Summary tag preserved
        assert "regulatory" in ann._determine_regions("chr1", 1250)

    def test_both_regulatory_sources(self, mini_ref_dir: Path) -> None:
        """regulatory_source='both' loads both FANTOM5 and Ensembl."""
        from dataclasses import replace

        preset = replace(get_preset("comprehensive"), regulatory_source="both")
        ann = VCFAnnotator()
        ann.load_beds(preset, mini_ref_dir)
        # FANTOM5 interval at 1050-1100
        assert "regulatory" in ann._determine_regions("chr1", 1075)
        # Ensembl interval at 1200-1300
        assert "regulatory_promoter" in ann._determine_regions("chr1", 1250)
