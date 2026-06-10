"""Tests for dgra_prefilter.safetynet — SafetyNetProvider ABC + ClinVarSafetyNet."""

from __future__ import annotations

from pathlib import Path

import pytest

from dgra_prefilter.safetynet import ClinVarSafetyNet, OMIMSafetyNet, SafetyNetProvider


# ======================================================================
# SafetyNetProvider ABC
# ======================================================================

class TestSafetyNetProviderABC:
    """Tests for the SafetyNetProvider abstract base class."""

    def test_cannot_instantiate_abc(self) -> None:
        with pytest.raises(TypeError):
            SafetyNetProvider()  # type: ignore[abstract]

    def test_subclass_must_implement_methods(self) -> None:
        class IncompleteProvider(SafetyNetProvider):
            pass

        with pytest.raises(TypeError):
            IncompleteProvider()  # type: ignore[abstract]

    def test_subclass_with_all_methods(self) -> None:
        class FullProvider(SafetyNetProvider):
            def get_bed_path(self, ref_dir: Path) -> Path:
                return ref_dir / "test.bed"

            def get_tag(self) -> str:
                return "TestTag"

            def is_available(self, ref_dir: Path) -> bool:
                return True

        p = FullProvider()
        assert p.get_tag() == "TestTag"
        assert p.get_bed_path(Path("/data")) == Path("/data/test.bed")
        assert p.is_available(Path("/data")) is True


# ======================================================================
# ClinVarSafetyNet
# ======================================================================

class TestClinVarSafetyNet:
    """Tests for ClinVarSafetyNet."""

    def test_get_tag(self) -> None:
        sn = ClinVarSafetyNet()
        assert sn.get_tag() == "ClinVar"

    def test_format_tag(self) -> None:
        assert ClinVarSafetyNet.format_tag("3") == "ClinVar_3star"
        assert ClinVarSafetyNet.format_tag("1") == "ClinVar_1star"

    def test_get_bed_path(self, tmp_path: Path) -> None:
        sn = ClinVarSafetyNet()
        assert sn.get_bed_path(tmp_path) == tmp_path / "clinvar_pathogenic_GRCh38.bed"

    def test_bed_filename_constant(self) -> None:
        assert ClinVarSafetyNet.BED_FILENAME == "clinvar_pathogenic_GRCh38.bed"

    def test_is_available_with_data(self, clinvar_bed: Path) -> None:
        sn = ClinVarSafetyNet()
        ref_dir = clinvar_bed.parent
        assert sn.is_available(ref_dir) is True

    def test_is_available_missing_file(self, tmp_path: Path) -> None:
        sn = ClinVarSafetyNet()
        assert sn.is_available(tmp_path) is False

    def test_is_available_empty_file(self, tmp_path: Path) -> None:
        sn = ClinVarSafetyNet()
        bed_path = tmp_path / "clinvar_pathogenic_GRCh38.bed"
        bed_path.write_text("")
        assert sn.is_available(tmp_path) is False

    def test_is_available_header_only(self, tmp_path: Path) -> None:
        sn = ClinVarSafetyNet()
        bed_path = tmp_path / "clinvar_pathogenic_GRCh38.bed"
        bed_path.write_text("# header\ntrack name=test\n")
        assert sn.is_available(tmp_path) is False

    def test_is_available_comment_and_data(self, tmp_path: Path) -> None:
        sn = ClinVarSafetyNet()
        bed_path = tmp_path / "clinvar_pathogenic_GRCh38.bed"
        bed_path.write_text("# header\nchr1\t100\t200\n")
        assert sn.is_available(tmp_path) is True


# ======================================================================
# OMIMSafetyNet
# ======================================================================

class TestOMIMSafetyNet:
    """Tests for OMIMSafetyNet."""

    def test_get_tag(self) -> None:
        sn = OMIMSafetyNet()
        assert sn.get_tag() == "OMIM"

    def test_get_bed_path(self, tmp_path: Path) -> None:
        sn = OMIMSafetyNet()
        assert sn.get_bed_path(tmp_path) == tmp_path / "omim_pathogenic_GRCh38.bed"

    def test_bed_filename_constant(self) -> None:
        assert OMIMSafetyNet.BED_FILENAME == "omim_pathogenic_GRCh38.bed"

    def test_is_available_missing_file(self, tmp_path: Path) -> None:
        sn = OMIMSafetyNet()
        assert sn.is_available(tmp_path) is False

    def test_is_available_with_data(self, tmp_path: Path) -> None:
        sn = OMIMSafetyNet()
        bed_path = tmp_path / "omim_pathogenic_GRCh38.bed"
        bed_path.write_text("chr1\t100\t200\n")
        assert sn.is_available(tmp_path) is True

    def test_is_available_empty_file(self, tmp_path: Path) -> None:
        sn = OMIMSafetyNet()
        bed_path = tmp_path / "omim_pathogenic_GRCh38.bed"
        bed_path.write_text("")
        assert sn.is_available(tmp_path) is False
