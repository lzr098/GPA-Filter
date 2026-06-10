"""Tests for dgra_prefilter.presets — preset configuration correctness."""

from __future__ import annotations

import pytest

from dgra_prefilter.presets import PRESETS, PresetConfig, get_preset
from dgra_prefilter.safetynet import ClinVarSafetyNet, OMIMSafetyNet


# ======================================================================
# PresetConfig dataclass
# ======================================================================

class TestPresetConfig:
    """Tests for PresetConfig."""

    def test_frozen(self) -> None:
        p = PRESETS["comprehensive"]
        with pytest.raises(AttributeError):
            p.gene = False  # type: ignore[misc]

    def test_get_region_bed_names_comprehensive(self) -> None:
        p = PRESETS["comprehensive"]
        beds = p.get_region_bed_names()
        assert "gencode_v44_gene_loci.bed" in beds
        assert "gencode_v44_ncrna_loci.bed" in beds
        assert "encode_screen_v3_ccres.bed" in beds
        assert "fantom5_enhancers_promoters.bed" in beds
        assert "vista_enhancers.bed" in beds
        # comprehensive uses full cCREs, not pls_pels_only or balanced
        assert "encode_screen_v3_pls_pels.bed" not in beds
        assert "encode_screen_v3_balanced.bed" not in beds

    def test_get_region_bed_names_coding_only(self) -> None:
        p = PRESETS["coding-only"]
        beds = p.get_region_bed_names()
        # coding-only uses the sub-region BEDs
        assert "gencode_5utr.bed" in beds
        assert "gencode_cds.bed" in beds
        assert "gencode_3utr.bed" in beds
        assert "gencode_splice_sites.bed" in beds
        assert "gencode_v44_gene_loci.bed" not in beds
        assert "gencode_v44_ncrna_loci.bed" not in beds
        # No regulatory
        assert "encode_screen_v3_ccres.bed" not in beds
        assert "fantom5_enhancers_promoters.bed" not in beds
        assert "vista_enhancers.bed" not in beds

    def test_get_region_bed_names_regulatory_minimal(self) -> None:
        p = PRESETS["regulatory-minimal"]
        beds = p.get_region_bed_names()
        assert "gencode_v44_gene_loci.bed" in beds
        assert "gencode_v44_ncrna_loci.bed" in beds
        # regulatory-minimal uses pls_pels_only
        assert "encode_screen_v3_pls_pels.bed" in beds
        assert "encode_screen_v3_ccres.bed" not in beds
        assert "encode_screen_v3_balanced.bed" not in beds
        # No FANTOM5 or VISTA
        assert "fantom5_enhancers_promoters.bed" not in beds
        assert "vista_enhancers.bed" not in beds

    def test_get_region_bed_names_regulatory_balanced(self) -> None:
        p = PRESETS["regulatory-balanced"]
        beds = p.get_region_bed_names()
        assert "gencode_v44_gene_loci.bed" in beds
        assert "gencode_v44_ncrna_loci.bed" in beds
        assert "encode_screen_v3_balanced.bed" in beds
        assert "encode_screen_v3_ccres.bed" not in beds
        assert "encode_screen_v3_pls_pels.bed" not in beds
        assert "fantom5_enhancers_promoters.bed" not in beds
        assert "vista_enhancers.bed" not in beds

    def test_get_safetynet_providers_comprehensive(self) -> None:
        p = PRESETS["comprehensive"]
        providers = p.get_safetynet_providers()
        assert len(providers) == 2
        assert any(isinstance(p, ClinVarSafetyNet) for p in providers)
        assert any(isinstance(p, OMIMSafetyNet) for p in providers)

    def test_get_safetynet_providers_coding_only(self) -> None:
        p = PRESETS["coding-only"]
        providers = p.get_safetynet_providers()
        assert len(providers) == 2
        assert any(isinstance(p, ClinVarSafetyNet) for p in providers)
        assert any(isinstance(p, OMIMSafetyNet) for p in providers)

    def test_get_safetynet_providers_regulatory_minimal(self) -> None:
        p = PRESETS["regulatory-minimal"]
        providers = p.get_safetynet_providers()
        assert len(providers) == 2
        assert any(isinstance(p, ClinVarSafetyNet) for p in providers)
        assert any(isinstance(p, OMIMSafetyNet) for p in providers)

    def test_get_safetynet_providers_regulatory_balanced(self) -> None:
        p = PRESETS["regulatory-balanced"]
        providers = p.get_safetynet_providers()
        assert len(providers) == 2
        assert any(isinstance(p, ClinVarSafetyNet) for p in providers)
        assert any(isinstance(p, OMIMSafetyNet) for p in providers)


# ======================================================================
# get_preset
# ======================================================================

class TestGetPreset:
    """Tests for get_preset() function."""

    def test_comprehensive(self) -> None:
        p = get_preset("comprehensive")
        assert p.name == "comprehensive"

    def test_coding_only(self) -> None:
        p = get_preset("coding-only")
        assert p.name == "coding-only"

    def test_regulatory_minimal(self) -> None:
        p = get_preset("regulatory-minimal")
        assert p.name == "regulatory-minimal"

    def test_regulatory_balanced(self) -> None:
        p = get_preset("regulatory-balanced")
        assert p.name == "regulatory-balanced"

    def test_unknown_preset_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown preset"):
            get_preset("nonexistent")


# ======================================================================
# Preset field correctness
# ======================================================================

class TestPresetFieldCorrectness:
    """Verify each preset has the expected field values."""

    def test_comprehensive_flags(self) -> None:
        p = PRESETS["comprehensive"]
        assert p.gene is True
        assert p.ncrna is True
        assert p.regulatory_encode is True
        assert p.regulatory_encode_pls_pels_only is False
        assert p.regulatory_encode_balanced is False
        assert p.regulatory_fantom5 is True
        assert p.regulatory_vista is True
        assert p.safetynet_clinvar is True
        assert p.safetynet_omim is True
        assert p.regulatory_source == "fantom5"
        assert p.keep_all_chrM is False

    def test_coding_only_flags(self) -> None:
        p = PRESETS["coding-only"]
        assert p.gene is True
        assert p.ncrna is False
        assert p.regulatory_encode is False
        assert p.regulatory_encode_pls_pels_only is False
        assert p.regulatory_encode_balanced is False
        assert p.regulatory_fantom5 is False
        assert p.regulatory_vista is False
        assert p.safetynet_clinvar is True
        assert p.safetynet_omim is True
        assert p.regulatory_source == "fantom5"
        assert p.keep_all_chrM is False

    def test_regulatory_minimal_flags(self) -> None:
        p = PRESETS["regulatory-minimal"]
        assert p.gene is True
        assert p.ncrna is True
        assert p.regulatory_encode is True
        assert p.regulatory_encode_pls_pels_only is True
        assert p.regulatory_encode_balanced is False
        assert p.regulatory_fantom5 is False
        assert p.regulatory_vista is False
        assert p.safetynet_clinvar is True
        assert p.safetynet_omim is True
        assert p.regulatory_source == "fantom5"
        assert p.keep_all_chrM is False

    def test_regulatory_balanced_flags(self) -> None:
        p = PRESETS["regulatory-balanced"]
        assert p.gene is True
        assert p.ncrna is True
        assert p.regulatory_encode is True
        assert p.regulatory_encode_pls_pels_only is False
        assert p.regulatory_encode_balanced is True
        assert p.regulatory_fantom5 is False
        assert p.regulatory_vista is False
        assert p.safetynet_clinvar is True
        assert p.safetynet_omim is True
        assert p.regulatory_source == "fantom5"
        assert p.keep_all_chrM is False

    def test_four_presets_exist(self) -> None:
        assert len(PRESETS) == 4
        assert set(PRESETS.keys()) == {
            "comprehensive",
            "coding-only",
            "regulatory-minimal",
            "regulatory-balanced",
        }

    def test_regulatory_source_fantom5(self) -> None:
        p = PRESETS["comprehensive"]
        beds = p.get_region_bed_names()
        assert "fantom5_enhancers_promoters.bed" in beds
        assert "ensembl_regulatory_features.bed" not in beds

    def test_regulatory_source_ensembl(self) -> None:
        from dataclasses import replace
        p = replace(PRESETS["comprehensive"], regulatory_source="ensembl")
        beds = p.get_region_bed_names()
        assert "fantom5_enhancers_promoters.bed" not in beds
        assert "ensembl_regulatory_features.bed" in beds

    def test_regulatory_source_both(self) -> None:
        from dataclasses import replace
        p = replace(PRESETS["comprehensive"], regulatory_source="both")
        beds = p.get_region_bed_names()
        assert "fantom5_enhancers_promoters.bed" in beds
        assert "ensembl_regulatory_features.bed" in beds
