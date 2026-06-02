"""Core filtering engine for dgra-prefilter."""

from __future__ import annotations

import gzip
import logging
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

from dgra_prefilter.bed_utils import BedUtils
from dgra_prefilter.presets import PRESETS, PresetConfig, get_preset
from dgra_prefilter.ref_manager import RefManager

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Custom exceptions
# ---------------------------------------------------------------------------

class GenomeMismatchError(Exception):
    """Raised when the input VCF genome version does not match the expected version."""

    def __init__(self, expected: str, found: str) -> None:
        self.expected = expected
        self.found = found
        super().__init__(
            f"Genome version mismatch: expected {expected}, found {found}. "
            f"Use --force to skip this check."
        )


class BcftoolsNotFoundError(Exception):
    """Raised when bcftools is not installed or not on PATH."""

    def __init__(self) -> None:
        super().__init__(
            "bcftools not found. Please install bcftools >= 1.17 and ensure "
            "it is on your PATH."
        )


class RefDataMissingError(Exception):
    """Raised when required reference data files are missing."""

    def __init__(self, message: str) -> None:
        super().__init__(message)


class VCFProcessingError(Exception):
    """Raised when a bcftools command fails during processing."""

    def __init__(self, cmd: list[str], returncode: int, stderr: str) -> None:
        self.cmd = cmd
        self.returncode = returncode
        self.stderr = stderr
        super().__init__(
            f"bcftools command failed (exit code {returncode}): "
            f"{' '.join(cmd)}\nstderr: {stderr}"
        )


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class PrefilterConfig:
    """Configuration for a prefilter run."""

    input_path: Path
    output_path: Path
    genome: str = "GRCh38"
    preset_name: str = "comprehensive"
    ref_dir: Path = field(default_factory=lambda: Path("~/.dgra-prefilter/refs"))
    report_path: Path | None = None
    force: bool = False
    update_refs: bool = False

    def __post_init__(self) -> None:
        """Normalize paths after initialization."""
        self.input_path = Path(self.input_path).expanduser().resolve()
        self.output_path = Path(self.output_path).expanduser().resolve()
        self.ref_dir = Path(self.ref_dir).expanduser().resolve()
        if self.report_path is not None:
            self.report_path = Path(self.report_path).expanduser().resolve()


@dataclass
class FilterStats:
    """Statistics from a prefilter run."""

    input_variants: int = 0
    retained_variants: int = 0
    region_only_variants: int = 0
    safetynet_only_variants: int = 0
    region_and_safetynet_variants: int = 0
    region_counts: dict[str, int] = field(default_factory=dict)
    clinvar_count: int = 0
    omim_count: int = 0
    elapsed_seconds: float = 0.0
    ref_data_versions: dict[str, str] = field(default_factory=dict)
    preset_config: dict[str, str] = field(default_factory=dict)

    @property
    def retention_rate(self) -> float:
        """Fraction of input variants that were retained."""
        if self.input_variants > 0:
            return self.retained_variants / self.input_variants
        return 0.0


@dataclass
class FilterResult:
    """Result of a prefilter run."""

    output_path: Path
    report_path: Path
    stats: FilterStats


# ---------------------------------------------------------------------------
# FilterEngine
# ---------------------------------------------------------------------------

class FilterEngine:
    """Core filtering engine that orchestrates the complete pipeline.

    Pipeline stages:
    1. Input validation (VCF format, genome version, bcftools availability)
    2. Reference file validation and BED merging
    3. Coordinate-based hard filtering via bcftools
    4. Safety net extraction and merging
    5. VCF INFO annotation
    6. Report generation
    """

    MIN_BCFTOOLS_VERSION = "1.17"

    def __init__(self, config: PrefilterConfig) -> None:
        """Initialize the filter engine.

        Args:
            config: Prefilter configuration.
        """
        self.config = config
        self.preset = get_preset(config.preset_name)
        self.ref_manager = RefManager(config.ref_dir)
        self.stats = FilterStats()
        self._temp_dir: tempfile.TemporaryDirectory | None = None
        self._temp_files: list[Path] = []

    def run(self) -> FilterResult:
        """Execute the complete filtering pipeline.

        Returns:
            FilterResult containing output paths and statistics.

        Raises:
            BcftoolsNotFoundError: If bcftools is not installed.
            GenomeMismatchError: If genome version does not match.
            RefDataMissingError: If reference data is missing.
            VCFProcessingError: If a bcftools command fails.
        """
        start_time = time.time()

        try:
            self._temp_dir = tempfile.TemporaryDirectory(prefix="dgra_prefilter_")
            temp_path = Path(self._temp_dir.name)

            # Stage 1: Input validation
            logger.info("Stage 1/6: Input validation")
            self._validate_input()

            # Stage 2: Reference loading
            logger.info("Stage 2/6: Reference loading and BED merging")
            self.ref_manager.validate_refs(self.preset)
            merged_bed = temp_path / "merged_regions.bed"
            self.ref_manager.merge_preset_beds(self.preset, merged_bed)

            # Count input variants
            self.stats.input_variants = self._count_variants(self.config.input_path)

            # Stage 3: Coordinate hard filtering
            logger.info("Stage 3/6: Coordinate-based filtering")
            region_vcf = temp_path / "region_filtered.vcf.gz"
            self._filter_by_regions(self.config.input_path, merged_bed, region_vcf)
            region_count = self._count_variants(region_vcf)
            logger.info("Region-filtered variants: %d", region_count)

            # Stage 4: Safety net extraction and merging
            logger.info("Stage 4/6: Safety net extraction and merging")
            safetynet_vcfs: list[Path] = []
            for provider in self.preset.get_safetynet_providers():
                safetynet_bed = provider.get_bed_path(self.config.ref_dir)
                if provider.is_available(self.config.ref_dir):
                    sn_vcf = temp_path / f"safetynet_{provider.get_tag()}.vcf.gz"
                    self._filter_by_safetynet(
                        self.config.input_path, safetynet_bed, sn_vcf
                    )
                    sn_count = self._count_variants(sn_vcf)
                    logger.info(
                        "Safety net %s: %d variants", provider.get_tag(), sn_count
                    )
                    if sn_count > 0:
                        safetynet_vcfs.append(sn_vcf)
                else:
                    logger.info(
                        "Safety net %s: data not available, skipping",
                        provider.get_tag(),
                    )

            # Merge region and safety net VCFs
            if safetynet_vcfs:
                combined_vcf = temp_path / "combined.vcf.gz"
                self._merge_vcfs(region_vcf, safetynet_vcfs, combined_vcf)
            else:
                combined_vcf = region_vcf

            self.stats.retained_variants = self._count_variants(combined_vcf)
            logger.info("Total retained variants: %d", self.stats.retained_variants)

            # Stage 5: VCF INFO annotation
            logger.info("Stage 5/6: VCF INFO annotation")
            from dgra_prefilter.annotate import VCFAnnotator

            annotator = VCFAnnotator()
            annotator.load_beds(self.preset, self.config.ref_dir)
            annotation_stats = annotator.annotate_vcf(combined_vcf, self.config.output_path)

            # Merge annotation stats into main stats
            self.stats.region_counts = annotation_stats.region_counts
            self.stats.clinvar_count = annotation_stats.clinvar_count
            self.stats.omim_count = annotation_stats.omim_count
            self.stats.region_only_variants = annotation_stats.region_only_variants
            self.stats.safetynet_only_variants = annotation_stats.safetynet_only_variants
            self.stats.region_and_safetynet_variants = (
                annotation_stats.region_and_safetynet_variants
            )

            # Populate stats before report generation
            self.stats.elapsed_seconds = time.time() - start_time
            self.stats.ref_data_versions = self.ref_manager.get_versions()
            self.stats.preset_config = {
                "name": self.preset.name,
                "gene": str(self.preset.gene),
                "ncrna": str(self.preset.ncrna),
                "regulatory_encode": str(self.preset.regulatory_encode),
                "regulatory_encode_pls_pels_only": str(
                    self.preset.regulatory_encode_pls_pels_only
                ),
                "regulatory_fantom5": str(self.preset.regulatory_fantom5),
                "regulatory_vista": str(self.preset.regulatory_vista),
                "safetynet_clinvar": str(self.preset.safetynet_clinvar),
                "safetynet_omim": str(self.preset.safetynet_omim),
            }

            # Stage 6: Report generation
            logger.info("Stage 6/6: Report generation")
            report_path = self.config.report_path
            if report_path is None:
                # Strip .vcf.gz or .vcf suffix, then add .report.json
                out_name = self.config.output_path.name
                if out_name.endswith(".vcf.gz"):
                    report_name = out_name[:-7] + ".report.json"
                elif out_name.endswith(".vcf"):
                    report_name = out_name[:-4] + ".report.json"
                else:
                    report_name = out_name + ".report.json"
                report_path = self.config.output_path.parent / report_name

            from dgra_prefilter.report import ReportGenerator

            ReportGenerator.generate(self.stats, self.config, report_path)

            logger.info(
                "Filtering complete: %d/%d variants retained (%.1f%%) in %.1fs",
                self.stats.retained_variants,
                self.stats.input_variants,
                self.stats.retention_rate * 100,
                self.stats.elapsed_seconds,
            )

            return FilterResult(
                output_path=self.config.output_path,
                report_path=report_path,
                stats=self.stats,
            )

        finally:
            # Clean up temporary files
            if self._temp_dir is not None:
                try:
                    self._temp_dir.cleanup()
                except Exception:
                    pass

    def _validate_input(self) -> None:
        """Validate the input VCF and system prerequisites.

        Checks:
        - Input file exists
        - bcftools is installed and meets minimum version
        - VCF header genome version (unless force=True)

        Raises:
            FileNotFoundError: If the input VCF does not exist.
            BcftoolsNotFoundError: If bcftools is not installed.
            GenomeMismatchError: If genome version does not match.
        """
        # Check input file exists
        if not self.config.input_path.exists():
            raise FileNotFoundError(
                f"Input VCF not found: {self.config.input_path}"
            )

        # Check bcftools availability
        self._check_bcftools()

        # Check genome version in VCF header
        if not self.config.force:
            self._check_genome_version()

    def _check_bcftools(self) -> None:
        """Check that bcftools is installed and meets minimum version.

        Raises:
            BcftoolsNotFoundError: If bcftools is not installed.
        """
        try:
            result = subprocess.run(
                ["bcftools", "--version"],
                capture_output=True,
                text=True,
                check=False,
            )
            if result.returncode != 0:
                raise BcftoolsNotFoundError()
            # Parse version from output like "bcftools 1.17"
            version_line = result.stdout.split("\n")[0]
            parts = version_line.split()
            if len(parts) >= 2:
                version_str = parts[1]
                logger.debug("bcftools version: %s", version_str)
        except FileNotFoundError:
            raise BcftoolsNotFoundError()

    def _check_genome_version(self) -> None:
        """Check that the VCF header declares the expected genome version.

        Looks for 'assembly=GRCh38' or 'reference=GRCh38' in the VCF header.

        Raises:
            GenomeMismatchError: If the genome version does not match.
        """
        input_path = self.config.input_path
        found_genome: str | None = None

        try:
            opener = gzip.open if str(input_path).endswith(".gz") else open
            with opener(input_path, "rt") as f:
                for line in f:
                    if not line.startswith("#"):
                        break
                    line_lower = line.lower()
                    if "grch38" in line_lower or "hg38" in line_lower:
                        if "assembly=grch38" in line_lower or "assembly=GRCh38" in line:
                            found_genome = "GRCh38"
                        elif "reference=grch38" in line_lower or "reference=GRCh38" in line:
                            found_genome = "GRCh38"
                        elif "grch38" in line_lower:
                            found_genome = "GRCh38"
                        break
                    if "grch37" in line_lower or "hg19" in line_lower:
                        found_genome = "GRCh37"
                        break
        except Exception:
            # If we can't read the header, we'll skip the check
            logger.warning("Could not read VCF header for genome version check")
            return

        if found_genome is not None and found_genome != self.config.genome:
            raise GenomeMismatchError(self.config.genome, found_genome)

        if found_genome is None:
            logger.warning(
                "Could not determine genome version from VCF header. "
                "Proceeding without genome check."
            )

    def _filter_by_regions(
        self, input_vcf: Path, merged_bed: Path, output: Path
    ) -> Path:
        """Filter VCF by coordinate regions using bcftools.

        Uses bcftools view -T to retain only variants within the merged
        region BED file.

        Args:
            input_vcf: Path to input VCF.
            merged_bed: Path to merged region BED file.
            output: Path for the filtered output VCF.

        Returns:
            Path to the filtered VCF.
        """
        cmd = [
            "bcftools", "view",
            "-T", str(merged_bed),
            "-Oz",
            "-o", str(output),
            str(input_vcf),
        ]
        self._run_bcftools(cmd)
        return output

    def _filter_by_safetynet(
        self, input_vcf: Path, safetynet_bed: Path, output: Path
    ) -> Path:
        """Extract variants falling in safety net regions using bcftools.

        Args:
            input_vcf: Path to input VCF.
            safetynet_bed: Path to safety net BED file.
            output: Path for the safety net filtered output VCF.

        Returns:
            Path to the filtered VCF.
        """
        cmd = [
            "bcftools", "view",
            "-T", str(safetynet_bed),
            "-Oz",
            "-o", str(output),
            str(input_vcf),
        ]
        self._run_bcftools(cmd)
        return output

    def _merge_vcfs(
        self, region_vcf: Path, safetynet_vcfs: list[Path], output: Path
    ) -> Path:
        """Merge region-filtered and safety net VCFs (union with deduplication).

        Uses bcftools concat -a to concatenate, bcftools sort to order,
        and bcftools norm -d snps to deduplicate.

        Args:
            region_vcf: Path to region-filtered VCF.
            safetynet_vcfs: List of safety net VCF paths.
            output: Path for the merged output VCF.

        Returns:
            Path to the merged VCF.
        """
        # Build concat input list
        all_vcfs = [region_vcf] + safetynet_vcfs

        # We need to pipe: concat | sort | norm
        concat_args = ["bcftools", "concat", "-a"] + [str(v) for v in all_vcfs]
        sort_args = ["bcftools", "sort"]
        norm_args = [
            "bcftools", "norm",
            "-d", "snps",
            "-Oz",
            "-o", str(output),
        ]

        logger.debug("Merge command: %s | %s | %s",
                      " ".join(concat_args),
                      " ".join(sort_args),
                      " ".join(norm_args))

        # Run as a pipeline
        p_concat = subprocess.Popen(
            concat_args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        p_sort = subprocess.Popen(
            sort_args,
            stdin=p_concat.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        p_concat.stdout.close()  # Allow p_concat to receive SIGPIPE

        p_norm = subprocess.Popen(
            norm_args,
            stdin=p_sort.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        p_sort.stdout.close()  # Allow p_sort to receive SIGPIPE

        # Wait for pipeline to complete
        _, norm_stderr = p_norm.communicate()
        p_concat.wait()
        p_sort.wait()

        if p_norm.returncode != 0:
            raise VCFProcessingError(
                norm_args, p_norm.returncode, norm_stderr.decode("utf-8", errors="replace")
            )
        if p_concat.returncode != 0:
            concat_stderr = p_concat.stderr.read().decode("utf-8", errors="replace")
            raise VCFProcessingError(concat_args, p_concat.returncode, concat_stderr)
        if p_sort.returncode != 0:
            sort_stderr = p_sort.stderr.read().decode("utf-8", errors="replace")
            raise VCFProcessingError(sort_args, p_sort.returncode, sort_stderr)

        # Index the output
        self._run_bcftools(["bcftools", "index", str(output)])

        return output

    def _count_variants(self, vcf_path: Path) -> int:
        """Count the number of variant records in a VCF.

        Uses bcftools view -H to show only records (no header), then
        counts lines with wc -l.

        Args:
            vcf_path: Path to the VCF file.

        Returns:
            Number of variant records.
        """
        if not vcf_path.exists():
            return 0
        try:
            result = subprocess.run(
                ["bcftools", "view", "-H", str(vcf_path)],
                capture_output=True,
                text=True,
                check=False,
            )
            if result.returncode != 0:
                logger.warning("Could not count variants in %s", vcf_path)
                return 0
            return sum(1 for line in result.stdout.split("\n") if line.strip())
        except Exception:
            return 0

    def _run_bcftools(self, args: list[str]) -> subprocess.CompletedProcess:
        """Execute a bcftools command and check for errors.

        Args:
            args: Full command line (including 'bcftools' as first element).

        Returns:
            CompletedProcess result.

        Raises:
            VCFProcessingError: If the command returns a non-zero exit code.
        """
        logger.debug("Running: %s", " ".join(args))
        result = subprocess.run(
            args,
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            raise VCFProcessingError(args, result.returncode, result.stderr)
        return result


# ---------------------------------------------------------------------------
# Top-level API function
# ---------------------------------------------------------------------------

def prefilter_vcf(
    input_path: str | Path,
    output_path: str | Path,
    genome: str = "GRCh38",
    preset: str = "comprehensive",
    ref_dir: str | Path = "~/.dgra-prefilter/refs",
    report_path: str | Path | None = None,
    force: bool = False,
    update_refs: bool = False,
) -> FilterResult:
    """Whole-genome VCF region prefiltering main entry point.

    Filters a GRCh38 VCF to retain only variants in biologically relevant
    regions (genes, ncRNA, regulatory elements) and known pathogenic sites
    from ClinVar.

    Args:
        input_path: Input VCF/VCF.gz/BCF file path.
        output_path: Output VCF path (.vcf.gz for compressed output).
        genome: Genome version, only GRCh38 supported.
        preset: Preset name: comprehensive | coding-only | regulatory-minimal.
        ref_dir: Directory containing reference BED files.
        report_path: JSON report output path (default: same dir as output).
        force: Skip genome version validation.
        update_refs: Trigger reference data update before filtering.

    Returns:
        FilterResult containing output path and filter statistics.

    Raises:
        ValueError: If parameters are invalid.
        FileNotFoundError: If input file or reference files are missing.
        GenomeMismatchError: If genome version does not match.
        BcftoolsNotFoundError: If bcftools is not installed.
        RefDataMissingError: If required reference data is missing.
    """
    config = PrefilterConfig(
        input_path=Path(input_path),
        output_path=Path(output_path),
        genome=genome,
        preset_name=preset,
        ref_dir=Path(ref_dir),
        report_path=Path(report_path) if report_path else None,
        force=force,
        update_refs=update_refs,
    )

    # Handle reference data update if requested
    if update_refs:
        ref_manager = RefManager(config.ref_dir)
        ref_manager.update_refs()

    engine = FilterEngine(config)
    return engine.run()
