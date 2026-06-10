"""Safety net providers for dgra-prefilter."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path


class SafetyNetProvider(ABC):
    """Abstract base class for safety net providers.

    A safety net ensures that known pathogenic variants are never filtered out,
    even if they fall outside the retained genomic regions.
    """

    @abstractmethod
    def get_bed_path(self, ref_dir: Path) -> Path:
        """Return the path to the pathogenic variants BED file.

        Args:
            ref_dir: Directory containing reference BED files.

        Returns:
            Path to the safety net BED file.
        """
        ...

    @abstractmethod
    def get_tag(self) -> str:
        """Return the DGRA_SAFETYNET INFO tag value.

        Returns:
            Tag string (e.g., 'ClinVar', 'OMIM').
        """
        ...

    @abstractmethod
    def is_available(self, ref_dir: Path) -> bool:
        """Check whether safety net data is available.

        A safety net is available when its BED file exists and is non-empty.

        Args:
            ref_dir: Directory containing reference BED files.

        Returns:
            True if the data file exists and contains data.
        """
        ...


class ClinVarSafetyNet(SafetyNetProvider):
    """ClinVar Pathogenic/Likely_pathogenic variants safety net.

    Extracts variant coordinates from ClinVar where CLNSIG contains
    'Pathogenic' or 'Likely_pathogenic' (including combined states
    like 'Pathogenic/Likely_pathogenic'). The BED file is expected to
    have a 4th column containing the ClinVar review-status star level
    (1-4) so annotations can emit star-specific tags (ClinVar_1star,
    ClinVar_2star, etc.). Zero-star records are filtered at build time.
    """

    BED_FILENAME = "clinvar_pathogenic_GRCh38.bed"

    def get_bed_path(self, ref_dir: Path) -> Path:
        """Return the ClinVar pathogenic variants BED file path.

        Args:
            ref_dir: Directory containing reference BED files.

        Returns:
            Path to clinvar_pathogenic_GRCh38.bed.
        """
        return ref_dir / self.BED_FILENAME

    def get_tag(self) -> str:
        """Return the base safety net tag.

        The per-variant annotation uses the star level from the BED's
        4th column via format_tag(). This method returns the provider
        identity tag used for logging and filenames.

        Returns:
            The string 'ClinVar'.
        """
        return "ClinVar"

    @staticmethod
    def format_tag(star_level: str) -> str:
        """Format a ClinVar tag including the star level.

        Args:
            star_level: Star level as a string (e.g., '1', '2', '3', '4').

        Returns:
            Tag string such as 'ClinVar_3star'.
        """
        return f"ClinVar_{star_level}star"

    def is_available(self, ref_dir: Path) -> bool:
        """Check whether ClinVar safety net data is available.

        Args:
            ref_dir: Directory containing reference BED files.

        Returns:
            True if clinvar_pathogenic_GRCh38.bed exists and is non-empty.
        """
        path = self.get_bed_path(ref_dir)
        if not path.exists():
            return False
        # Check that the file has data beyond a possible header line
        try:
            with open(path, "r") as f:
                for line in f:
                    stripped = line.strip()
                    if stripped and not stripped.startswith("track") and not stripped.startswith("#"):
                        return True
            return False
        except OSError:
            return False


# v2.0 implementation — activated 2026-06-07 with local OMIM SQLite
class OMIMSafetyNet(SafetyNetProvider):
    """OMIM Mendelian disease gene safety net (v2.0).

    Protects all known Mendelian disease genes from being filtered out.
    BED file contains gene coordinates for all OMIM entries with disease
    associations (prefix: #, %, or + with phenotypeMap entries).

    Generated via: bcftools query to extract gene coordinates from OMIM SQLite.
    """

    BED_FILENAME = "omim_pathogenic_GRCh38.bed"

    def get_bed_path(self, ref_dir: Path) -> Path:
        """Return the OMIM pathogenic variants BED file path."""
        return ref_dir / self.BED_FILENAME

    def get_tag(self) -> str:
        """Return the OMIM safety net tag."""
        return "OMIM"

    def is_available(self, ref_dir: Path) -> bool:
        """Check whether OMIM safety net data is available.

        Args:
            ref_dir: Directory containing reference BED files.

        Returns:
            True if omim_pathogenic_GRCh38.bed exists and is non-empty.
        """

        path = self.get_bed_path(ref_dir)
        if not path.exists():
            return False
        try:
            with open(path, "r") as f:
                for line in f:
                    stripped = line.strip()
                    if stripped and not stripped.startswith("track") and not stripped.startswith("#"):
                        return True
            return False
        except OSError:
            return False
