"""Preset configurations for dgra-prefilter."""

from __future__ import annotations

from dataclasses import dataclass

from dgra_prefilter.safetynet import ClinVarSafetyNet, SafetyNetProvider


@dataclass(frozen=True)
class PresetConfig:
    """Immutable preset configuration defining which regions and safety nets to use."""

    name: str
    gene: bool
    ncrna: bool
    regulatory_encode: bool
    regulatory_encode_pls_pels_only: bool
    regulatory_fantom5: bool
    regulatory_vista: bool
    safetynet_clinvar: bool
    safetynet_omim: bool

    def get_region_bed_names(self) -> list[str]:
        """Return the list of region BED filenames needed by this preset.

        Returns:
            List of BED filenames (e.g., 'gencode_v44_gene_loci.bed').
        """
        beds: list[str] = []
        if self.gene:
            if self.name == "coding-only":
                beds.append("gencode_v44_coding_exon_utr.bed")
            else:
                beds.append("gencode_v44_gene_loci.bed")
        if self.ncrna:
            beds.append("gencode_v44_ncrna_loci.bed")
        if self.regulatory_encode:
            if self.regulatory_encode_pls_pels_only:
                beds.append("encode_screen_v3_pls_pels.bed")
            else:
                beds.append("encode_screen_v3_ccres.bed")
        if self.regulatory_fantom5:
            beds.append("fantom5_enhancers_promoters.bed")
        if self.regulatory_vista:
            beds.append("vista_enhancers.bed")
        return beds

    def get_safetynet_providers(self) -> list[SafetyNetProvider]:
        """Return the list of safety net providers enabled by this preset.

        Returns:
            List of SafetyNetProvider instances.
        """
        providers: list[SafetyNetProvider] = []
        if self.safetynet_clinvar:
            providers.append(ClinVarSafetyNet())
        # v2.0: if self.safetynet_omim: providers.append(OMIMSafetyNet())
        return providers


PRESETS: dict[str, PresetConfig] = {
    "comprehensive": PresetConfig(
        name="comprehensive",
        gene=True,
        ncrna=True,
        regulatory_encode=True,
        regulatory_encode_pls_pels_only=False,
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
        regulatory_fantom5=False,
        regulatory_vista=False,
        safetynet_clinvar=True,
        safetynet_omim=False,
    ),
}


def get_preset(name: str) -> PresetConfig:
    """Look up a preset by name.

    Args:
        name: Preset name (comprehensive, coding-only, or regulatory-minimal).

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
