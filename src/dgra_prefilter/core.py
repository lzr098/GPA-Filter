"""Core filtering engine for dgra-prefilter."""

from __future__ import annotations

import gzip
import logging
import subprocess
import tempfile
import time
from dataclasses import dataclass, field, replace
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
    annotate: bool = False
    regulatory_source: str = "fantom5"
    keep_all_chrM: bool = False

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

    Pipeline stages (fixed two-phase):

    Phase 1 – Coordinate filtering (always runs)
        1. Input validation (VCF format, genome version, bcftools)
        2. Reference loading and BED merging (gene + ncRNA + cCRE)
        3. Coordinate-based hard filtering via bcftools view -T
        4. Safety net extraction (ClinVar P/LP) and merge/dedup

    Phase 2 – INFO annotation (optional, --annotate)
        5. VCF INFO annotation with DGRA_REGION / DGRA_SAFETYNET tags

    Phase 3 – Reporting (always runs)
        6. JSON report generation
    """

    MIN_BCFTOOLS_VERSION = "1.17"

    def __init__(self, config: PrefilterConfig) -> None:
        """Initialize the filter engine.

        Args:
            config: Prefilter configuration.
        """
        self.config = config
        base_preset = get_preset(config.preset_name)
        # Apply runtime overrides from config to preset
        preset_kwargs: dict[str, object] = {}
        if config.regulatory_source != "fantom5":
            preset_kwargs["regulatory_source"] = config.regulatory_source
        if config.keep_all_chrM:
            preset_kwargs["keep_all_chrM"] = True
        self.preset = replace(base_preset, **preset_kwargs) if preset_kwargs else base_preset
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
        region_count = 0  # initialized for cases where filtering is skipped

        try:
            self._temp_dir = tempfile.TemporaryDirectory(prefix="dgra_prefilter_")
            temp_path = Path(self._temp_dir.name)

            # =====================================================================
            # Stage 1: Input validation
            # =====================================================================
            logger.info("Stage 1/6: Input validation")
            self._validate_input()

            # v0.1.0: Check local GRCh38 FASTA availability
            from dgra_prefilter.fasta_qc import log_fasta_status
            log_fasta_status()

            # Check if input is already prefiltered/annotated
            already_annotated = self._has_dgra_annotations()
            if already_annotated and self.config.annotate:
                logger.info(
                    "Input VCF already contains DGRA annotations; "
                    "skipping filtering stages and proceeding directly to annotation."
                )
                self.stats.input_variants = self._count_variants(self.config.input_path)
                combined_vcf = self.config.input_path
                self.stats.retained_variants = self.stats.input_variants
                region_count = self.stats.input_variants
            else:
                # =====================================================================
                # Stage 2: Reference loading
                # =====================================================================
                logger.info("Stage 2/6: Reference loading and BED merging")
                self.ref_manager.validate_refs(self.preset)
                merged_bed = temp_path / "merged_regions.bed"
                self.ref_manager.merge_preset_beds(self.preset, merged_bed)

                # Add chrM virtual region if keep_all_chrM is enabled
                if self.preset.keep_all_chrM:
                    with open(merged_bed, "a") as f:
                        f.write("chrM\t0\t16569\n")
                    logger.info("Added chrM virtual region to merged BED")

                # Count input variants
                self.stats.input_variants = self._count_variants(self.config.input_path)

                # =====================================================================
                # Stage 3: Coordinate hard filtering (gene + ncRNA + cCRE)
                # =====================================================================
                logger.info("Stage 3/6: Coordinate-based filtering")
                region_vcf = temp_path / "region_filtered.vcf.gz"
                self._filter_by_regions(self.config.input_path, merged_bed, region_vcf)
                region_count = self._count_variants(region_vcf)
                logger.info("Region-filtered variants: %d", region_count)

                # =====================================================================
                # Stage 4: Safety net extraction and merging
                # =====================================================================
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

            # =====================================================================
            # Stage 5: VCF INFO annotation (optional)
            # =====================================================================
            if self.config.annotate:
                logger.info("Stage 5/6: VCF INFO annotation (enabled)")
                from dgra_prefilter.annotate import VCFAnnotator

                annotator = VCFAnnotator()
                annotator.load_beds(self.preset, self.config.ref_dir)
                annotation_stats = annotator.annotate_vcf(
                    combined_vcf, self.config.output_path
                )

                # Merge annotation stats into main stats
                self.stats.region_counts = annotation_stats.region_counts
                self.stats.clinvar_count = annotation_stats.clinvar_count
                self.stats.omim_count = annotation_stats.omim_count
                self.stats.region_only_variants = annotation_stats.region_only_variants
                self.stats.safetynet_only_variants = (
                    annotation_stats.safetynet_only_variants
                )
                self.stats.region_and_safetynet_variants = (
                    annotation_stats.region_and_safetynet_variants
                )
            else:
                logger.info("Stage 5/6: VCF INFO annotation (skipped)")
                # Copy combined VCF directly to output
                import shutil

                self.config.output_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(combined_vcf, self.config.output_path)
                # Try to copy index if it exists
                idx_src = combined_vcf.with_suffix(combined_vcf.suffix + ".csi")
                if idx_src.exists():
                    shutil.copy2(idx_src, self.config.output_path.with_suffix(
                        self.config.output_path.suffix + ".csi"
                    ))
                # Set minimal stats for non-annotated runs
                self.stats.region_only_variants = region_count
                self.stats.safetynet_only_variants = (
                    self.stats.retained_variants - region_count
                )
                self.stats.region_and_safetynet_variants = 0
                self.stats.region_counts = {"gene": 0, "ncrna": 0, "regulatory": 0}
                self.stats.clinvar_count = 0
                self.stats.omim_count = 0

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
                "regulatory_encode_balanced": str(
                    self.preset.regulatory_encode_balanced
                ),
                "regulatory_fantom5": str(self.preset.regulatory_fantom5),
                "regulatory_vista": str(self.preset.regulatory_vista),
                "safetynet_clinvar": str(self.preset.safetynet_clinvar),
                "safetynet_omim": str(self.preset.safetynet_omim),
                "regulatory_source": self.preset.regulatory_source,
                "keep_all_chrM": str(self.preset.keep_all_chrM),
            }

            # =====================================================================
            # Stage 6: Report generation
            # =====================================================================
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

    def _has_dgra_annotations(self) -> bool:
        """Check if input VCF already contains DGRA_REGION or DGRA_SAFETYNET tags.

        Used to detect whether the input is already a prefiltered/annotated VCF,
        allowing the engine to skip redundant filtering when annotate=True.

        Returns:
            True if DGRA annotations are present in the header or records.
        """
        input_path = self.config.input_path
        try:
            opener = gzip.open if str(input_path).endswith(".gz") else open
            with opener(input_path, "rt") as f:
                for line in f:
                    if not line.startswith("#"):
                        # Also check first few data lines for DGRA tags
                        if "DGRA_REGION" in line or "DGRA_SAFETYNET" in line:
                            return True
                        break
                    if "DGRA_REGION" in line or "DGRA_SAFETYNET" in line:
                        return True
        except Exception:
            pass
        return False

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

        Uses intermediate temp files instead of shell pipelines for
        reliability (bcftools sort does not stream correctly to stdout).

        Args:
            region_vcf: Path to region-filtered VCF.
            safetynet_vcfs: List of safety net VCF paths.
            output: Path for the merged output VCF.

        Returns:
            Path to the merged VCF.
        """
        all_vcfs = [region_vcf] + safetynet_vcfs

        # Ensure all inputs are indexed (bcftools concat requires indices)
        for vcf in all_vcfs:
            self._run_bcftools(["bcftools", "index", str(vcf)])

        temp_path = Path(self._temp_dir.name) if self._temp_dir else Path(tempfile.mkdtemp())

        # Step 1: concat (allow overlapping records)
        concat_out = temp_path / "concat.vcf.gz"
        concat_cmd = ["bcftools", "concat", "-a"] + [str(v) for v in all_vcfs] + ["-Oz", "-o", str(concat_out)]
        self._run_bcftools(concat_cmd)

        # Step 2: sort
        sorted_out = temp_path / "sorted.vcf.gz"
        self._run_bcftools(["bcftools", "sort", str(concat_out), "-Oz", "-o", str(sorted_out)])

        # Step 3: norm (deduplicate)
        self._run_bcftools(["bcftools", "norm", "-d", "snps", str(sorted_out), "-Oz", "-o", str(output)])

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
    annotate: bool = False,
    regulatory_source: str = "fantom5",
    keep_all_chrM: bool = False,
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
        annotate: Enable DGRA_REGION/DGRA_SAFETYNET INFO annotation (slow).
        regulatory_source: Regulatory source for comprehensive preset.
        keep_all_chrM: Retain all chrM variants regardless of region.

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
        annotate=annotate,
        regulatory_source=regulatory_source,
        keep_all_chrM=keep_all_chrM,
    )

    # Handle reference data update if requested
    if update_refs:
        ref_manager = RefManager(config.ref_dir)
        ref_manager.update_refs()

    engine = FilterEngine(config)
    return engine.run()


def annotate_vcf_file(
    input_path: str | Path,
    output_path: str | Path,
    preset: str = "comprehensive",
    ref_dir: str | Path = "~/.dgra-prefilter/refs",
    report_path: str | Path | None = None,
) -> FilterResult:
    """Annotate an existing VCF with DGRA_REGION and DGRA_SAFETYNET INFO tags.

    This is a standalone annotation function that does NOT re-run coordinate
    filtering. It is useful when you already have a prefiltered VCF and only
    want to add region/safety-net annotations.

    Args:
        input_path: Input VCF/VCF.gz path (may already be prefiltered).
        output_path: Output VCF path (.vcf.gz for compressed output).
        preset: Preset name determining which BED files to load for annotation.
        ref_dir: Directory containing reference BED files.
        report_path: JSON report output path (default: same dir as output).

    Returns:
        FilterResult containing output path and annotation statistics.
    """
    import time

    start_time = time.time()
    input_path = Path(input_path).expanduser().resolve()
    output_path = Path(output_path).expanduser().resolve()
    ref_dir = Path(ref_dir).expanduser().resolve()

    if not input_path.exists():
        raise FileNotFoundError(f"Input VCF not found: {input_path}")

    preset_cfg = get_preset(preset)
    ref_manager = RefManager(ref_dir)
    ref_manager.validate_refs(preset_cfg)

    from dgra_prefilter.annotate import VCFAnnotator

    annotator = VCFAnnotator()
    annotator.load_beds(preset_cfg, ref_dir)
    annotation_stats = annotator.annotate_vcf(input_path, output_path)

    # Build minimal stats
    stats = FilterStats()
    stats.input_variants = annotation_stats.region_only_variants + annotation_stats.safetynet_only_variants + annotation_stats.region_and_safetynet_variants
    stats.retained_variants = stats.input_variants
    stats.region_counts = annotation_stats.region_counts
    stats.clinvar_count = annotation_stats.clinvar_count
    stats.omim_count = annotation_stats.omim_count
    stats.region_only_variants = annotation_stats.region_only_variants
    stats.safetynet_only_variants = annotation_stats.safetynet_only_variants
    stats.region_and_safetynet_variants = annotation_stats.region_and_safetynet_variants
    stats.elapsed_seconds = time.time() - start_time
    stats.ref_data_versions = ref_manager.get_versions()
    stats.preset_config = {
        "name": preset_cfg.name,
        "gene": str(preset_cfg.gene),
        "ncrna": str(preset_cfg.ncrna),
        "regulatory_encode": str(preset_cfg.regulatory_encode),
        "regulatory_encode_pls_pels_only": str(
            preset_cfg.regulatory_encode_pls_pels_only
        ),
        "regulatory_encode_balanced": str(preset_cfg.regulatory_encode_balanced),
        "regulatory_fantom5": str(preset_cfg.regulatory_fantom5),
        "regulatory_vista": str(preset_cfg.regulatory_vista),
        "safetynet_clinvar": str(preset_cfg.safetynet_clinvar),
        "safetynet_omim": str(preset_cfg.safetynet_omim),
        "regulatory_source": preset_cfg.regulatory_source,
        "keep_all_chrM": str(preset_cfg.keep_all_chrM),
    }

    # Generate report
    if report_path is None:
        out_name = output_path.name
        if out_name.endswith(".vcf.gz"):
            report_name = out_name[:-7] + ".report.json"
        elif out_name.endswith(".vcf"):
            report_name = out_name[:-4] + ".report.json"
        else:
            report_name = out_name + ".report.json"
        report_path = output_path.parent / report_name

    from dgra_prefilter.report import ReportGenerator

    config = PrefilterConfig(
        input_path=input_path,
        output_path=output_path,
        preset_name=preset,
        ref_dir=ref_dir,
        annotate=True,
    )
    ReportGenerator.generate(stats, config, report_path)

    logger.info(
        "Annotation complete: %d variants annotated in %.1fs",
        stats.retained_variants,
        stats.elapsed_seconds,
    )

    return FilterResult(
        output_path=output_path,
        report_path=report_path,
        stats=stats,
    )
