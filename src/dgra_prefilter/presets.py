"""Preset configurations for dgra-prefilter."""

from __future__ import annotations

from dataclasses import dataclass

from dgra_prefilter.constants import (
    ENCODE_BALANCED_BED,
    ENCODE_CCRE_BED,
    ENCODE_PLS_PELS_BED,
    ENSEMBL_REGULATORY_BED,
    FANTOM5_BED,
    GENCODE_3UTR_BED,
    GENCODE_5UTR_BED,
    GENCODE_CDS_BED,
    GENCODE_CODING_EXON_UTR_BED,
    GENCODE_GENE_BED,
    GENCODE_NCRNA_BED,
    GENCODE_SPLICE_BED,
    VISTA_BED,
)
from dgra_prefilter.safetynet import ClinVarSafetyNet, OMIMSafetyNet, SafetyNetProvider


@dataclass(frozen=True)
class PresetConfig:
    """Immutable preset configuration defining which regions and safety nets to use."""

    name: str
    gene: bool
    ncrna: bool
    regulatory_encode: bool
    regulatory_encode_pls_pels_only: bool
    regulatory_encode_balanced: bool
    regulatory_fantom5: bool
    regulatory_vista: bool
    safetynet_clinvar: bool
    safetynet_omim: bool
    regulatory_source: str = "fantom5"
    keep_all_chrM: bool = False
    splice_window: int = 0

    def get_region_bed_names(self) -> list[str]:
        """Return the list of region BED filenames needed by this preset.

        Returns:
            List of BED filenames (e.g., 'gencode_v44_gene_loci.bed').
        """
        beds: list[str] = []
        if self.gene:
            if self.name == "coding-only":
                beds.extend([
                    GENCODE_5UTR_BED,
                    GENCODE_CDS_BED,
                    GENCODE_3UTR_BED,
                    GENCODE_SPLICE_BED,
                ])
            elif self.splice_window > 0:
                beds.append(GENCODE_CODING_EXON_UTR_BED)
            else:
                beds.append(GENCODE_GENE_BED)
        if self.ncrna:
            beds.append(GENCODE_NCRNA_BED)
        if self.regulatory_encode:
            if self.regulatory_encode_balanced:
                beds.append(ENCODE_BALANCED_BED)
            elif self.regulatory_encode_pls_pels_only:
                beds.append(ENCODE_PLS_PELS_BED)
            else:
                beds.append(ENCODE_CCRE_BED)
        if self.regulatory_fantom5:
            if self.regulatory_source == "fantom5":
                beds.append(FANTOM5_BED)
            elif self.regulatory_source == "ensembl":
                beds.append(ENSEMBL_REGULATORY_BED)
            elif self.regulatory_source == "both":
                beds.append(FANTOM5_BED)
                beds.append(ENSEMBL_REGULATORY_BED)
        if self.regulatory_vista:
            beds.append(VISTA_BED)
        return beds

    def get_safetynet_providers(self) -> list[SafetyNetProvider]:
        """Return the list of safety net providers enabled by this preset.

        Returns:
            List of SafetyNetProvider instances.
        """
        providers: list[SafetyNetProvider] = []
        if self.safetynet_clinvar:
            providers.append(ClinVarSafetyNet())
        if self.safetynet_omim:
            providers.append(OMIMSafetyNet())
        return providers


PRESETS: dict[str, PresetConfig] = {
    "comprehensive": PresetConfig(
        name="comprehensive",
        gene=True,
        ncrna=True,
        regulatory_encode=True,
        regulatory_encode_pls_pels_only=False,
        regulatory_encode_balanced=False,
        regulatory_fantom5=True,
        regulatory_vista=True,
        safetynet_clinvar=True,
        safetynet_omim=False,
    ),
    "coding-only": PresetConfig(
        name="coding-only",
        gene=True,
        ncrna=False,
        regulatory_encode=False,
        regulatory_encode_pls_pels_only=False,
        regulatory_encode_balanced=False,
        regulatory_fantom5=False,
        regulatory_vista=False,
        safetynet_clinvar=True,
        safetynet_omim=False,
    ),
    "regulatory-minimal": PresetConfig(
        name="regulatory-minimal",
        gene=True,
        ncrna=True,
        regulatory_encode=True,
        regulatory_encode_pls_pels_only=True,
        regulatory_encode_balanced=False,
        regulatory_fantom5=False,
        regulatory_vista=False,
        safetynet_clinvar=True,
        safetynet_omim=False,
    ),
    "regulatory-balanced": PresetConfig(
        name="regulatory-balanced",
        gene=True,
        ncrna=True,
        regulatory_encode=True,
        regulatory_encode_pls_pels_only=False,
        regulatory_encode_balanced=True,
        regulatory_fantom5=False,
        regulatory_vista=False,
        safetynet_clinvar=True,
        safetynet_omim=False,
    ),
    "comprehensive-splice100": PresetConfig(
        name="comprehensive-splice100",
        gene=True,
        ncrna=True,
        regulatory_encode=True,
        regulatory_encode_pls_pels_only=False,
        regulatory_encode_balanced=False,
        regulatory_fantom5=True,
        regulatory_vista=True,
        safetynet_clinvar=True,
        safetynet_omim=False,
        splice_window=100,
    ),
}


def get_preset(name: str) -> PresetConfig:
    """Look up a preset by name.

    Args:
        name: Preset name (comprehensive, comprehensive-splice100, coding-only, regulatory-minimal, regulatory-balanced).

    Returns:
        The corresponding PresetConfig.

    Raises:
        ValueError: If the preset name is not recognized.
    """
    if name not in PRESETS:
        valid = ", ".join(sorted(PRESETS.keys()))
        raise ValueError(
            f"Unknown preset '{name}'. Valid presets: {valid}"
        )
    return PRESETS[name]
