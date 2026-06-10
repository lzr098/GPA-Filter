"""Reference data manager for dgra-prefilter."""

from __future__ import annotations

import gzip
import logging
import subprocess
import tarfile
import tempfile
from pathlib import Path

from dgra_prefilter.bed_utils import BedUtils
from dgra_prefilter.constants import GITHUB_RELEASE_URL
from dgra_prefilter.presets import PresetConfig

logger = logging.getLogger(__name__)


class RefManager:
    """Manager for reference BED files used by the prefilter.

    Handles validation, path resolution, BED merging, version tracking,
    and reference data updates.
    """

    def __init__(self, ref_dir: Path) -> None:
        """Initialize the reference manager.

        Args:
            ref_dir: Path to the directory containing reference BED files.
        """
        self.ref_dir = ref_dir.expanduser().resolve()

    def validate_refs(self, preset: PresetConfig) -> bool:
        """Validate that all reference files required by a preset exist.

        Args:
            preset: The preset configuration to validate.

        Returns:
            True if all required files exist.

        Raises:
            FileNotFoundError: If any required reference file is missing.
        """
        # Check region BED files
        bed_names = preset.get_region_bed_names()
        for name in bed_names:
            path = self.get_bed_path(name)
            if not self._check_file_exists(path):
                raise FileNotFoundError(
                    f"Required reference file missing: {path}. "
                    f"Run 'dgra-prefilter --update-refs' to download."
                )

        # Check safety net BED files
        for provider in preset.get_safetynet_providers():
            bed_path = provider.get_bed_path(self.ref_dir)
            # Safety net files may legitimately be empty (e.g., OMIM),
            # but they must at least exist. Missing safety nets are warned
            # about rather than failing hard so the pipeline can still run.
            if not bed_path.exists():
                logger.warning(
                    "Safety net file missing: %s. Skipping %s safety net.",
                    bed_path,
                    provider.get_tag(),
                )

        return True

    def get_bed_path(self, bed_name: str) -> Path:
        """Get the full path for a named BED file.

        Args:
            bed_name: BED filename (e.g., 'gencode_v44_gene_loci.bed').

        Returns:
            Full path to the BED file in the reference directory.
        """
        return self.ref_dir / bed_name

    def merge_preset_beds(self, preset: PresetConfig, output: Path) -> Path:
        """Merge all region BED files required by a preset into one file.

        Args:
            preset: The preset configuration.
            output: Path for the merged output BED file.

        Returns:
            Path to the merged BED file.
        """
        bed_names = preset.get_region_bed_names()
        paths = [self.get_bed_path(name) for name in bed_names]
        return BedUtils.merge_bed_files(paths, output)

    def update_refs(self) -> None:
        """Update reference data by downloading from GitHub Release.

        MVP: Full update only (no incremental). Downloads the latest
        pre-built reference data package and extracts it to ref_dir.

        Raises:
            RuntimeError: If the update fails.
        """
        logger.info("Updating reference data from %s", GITHUB_RELEASE_URL)

        self.ref_dir.mkdir(parents=True, exist_ok=True)

        try:
            # Download to a temporary file
            with tempfile.TemporaryDirectory() as tmp_dir:
                archive_path = Path(tmp_dir) / "refs.tar.gz"

                # Try downloading with curl
                try:
                    subprocess.run(
                        ["curl", "-L", "-o", str(archive_path), GITHUB_RELEASE_URL],
                        check=True,
                        capture_output=True,
                        text=True,
                    )
                except (subprocess.CalledProcessError, FileNotFoundError):
                    # Fallback to wget
                    try:
                        subprocess.run(
                            ["wget", "-O", str(archive_path), GITHUB_RELEASE_URL],
                            check=True,
                            capture_output=True,
                            text=True,
                        )
                    except (subprocess.CalledProcessError, FileNotFoundError):
                        raise RuntimeError(
                            "Failed to download reference data. "
                            "Please install curl or wget, or download manually "
                            f"from {GITHUB_RELEASE_URL}"
                        )

                # Extract the archive
                with tarfile.open(archive_path, "r:gz") as tar:
                    tar.extractall(path=self.ref_dir, filter="data")

                logger.info("Reference data updated successfully in %s", self.ref_dir)

        except RuntimeError:
            raise
        except Exception as exc:
            raise RuntimeError(
                f"Failed to update reference data: {exc}"
            ) from exc

    def get_versions(self) -> dict[str, str]:
        """Read version information from .version files in the reference directory.

        Each BED file may have a corresponding .version file (e.g.,
        gencode_v44_gene_loci.bed.version) containing version info.

        Returns:
            Dictionary mapping data source name to version string.
        """
        versions: dict[str, str] = {}
        if not self.ref_dir.exists():
            return versions

        for version_file in sorted(self.ref_dir.glob("*.version")):
            # Derive source name from the .version file
            # e.g., gencode_v44_gene_loci.bed.version -> gencode_v44_gene_loci.bed
            bed_name = version_file.name[: -len(".version")]
            try:
                version_text = version_file.read_text().strip()
                if version_text:
                    versions[bed_name] = version_text
            except OSError:
                logger.warning("Could not read version file: %s", version_file)

        return versions

    def _check_file_exists(self, path: Path) -> bool:
        """Check that a file exists.

        Args:
            path: Path to check.

        Returns:
            True if the file exists.
        """
        return path.exists()

    def _check_file_nonempty(self, path: Path) -> bool:
        """Check that a file exists and contains data beyond header lines.

        Args:
            path: Path to check.

        Returns:
            True if the file has at least one data line.
        """
        if not path.exists():
            return False
        try:
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt") as f:
                for line in f:
                    stripped = line.strip()
                    if stripped and not stripped.startswith("#") and not stripped.startswith("track"):
                        return True
            return False
        except OSError:
            return False
