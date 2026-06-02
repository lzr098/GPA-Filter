"""VCF INFO field annotator for dgra-prefilter."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

from dgra_prefilter.bed_utils import BedUtils
from dgra_prefilter.constants import (
    CLINVAR_BED,
    ENCODE_CCRE_BED,
    ENCODE_PLS_PELS_BED,
    FANTOM5_BED,
    GENCODE_CODING_EXON_UTR_BED,
    GENCODE_GENE_BED,
    GENCODE_NCRNA_BED,
    OMIM_BED,
    VISTA_BED,
)
from dgra_prefilter.core import FilterStats
from dgra_prefilter.presets import PresetConfig

logger = logging.getLogger(__name__)


class VCFAnnotator:
    """VCF INFO field annotator.

    Adds DGRA_REGION and DGRA_SAFETYNET INFO tags to each variant record
    based on overlap with loaded BED files.
    """

    # VCF header INFO field definitions
    REGION_INFO_DEF = (
        '##INFO=<ID=DGRA_REGION,Number=.,Type=String,'
        'Description="Retained region type: gene|ncrna|regulatory">'
    )
    SAFETYNET_INFO_DEF = (
        '##INFO=<ID=DGRA_SAFETYNET,Number=.,Type=String,'
        'Description="Safety net source: ClinVar|OMIM">'
    )

    def __init__(self) -> None:
        """Initialize the annotator with empty BED data."""
        self.region_beds: dict[str, dict[str, list[tuple[int, int]]]] = {}
        self.safetynet_beds: dict[str, dict[str, list[tuple[int, int]]]] = {}

    def load_beds(self, preset: PresetConfig, ref_dir: Path) -> None:
        """Load all BED files needed for annotation based on the preset.

        Region BEDs are mapped to category labels (gene, ncrna, regulatory).
        Safety net BEDs are mapped by their tag names (ClinVar, OMIM).

        Args:
            preset: The preset configuration.
            ref_dir: Directory containing reference BED files.
        """
        ref_dir = ref_dir.expanduser().resolve()

        # Load region BEDs
        if preset.gene:
            if preset.name == "coding-only":
                gene_bed_path = ref_dir / GENCODE_CODING_EXON_UTR_BED
                self.region_beds["gene"] = BedUtils.load_bed(gene_bed_path)
            else:
                gene_bed_path = ref_dir / GENCODE_GENE_BED
                self.region_beds["gene"] = BedUtils.load_bed(gene_bed_path)

        if preset.ncrna:
            ncrna_bed_path = ref_dir / GENCODE_NCRNA_BED
            self.region_beds["ncrna"] = BedUtils.load_bed(ncrna_bed_path)

        # Regulatory regions
        if preset.regulatory_encode:
            if preset.regulatory_encode_pls_pels_only:
                encode_bed_path = ref_dir / ENCODE_PLS_PELS_BED
            else:
                encode_bed_path = ref_dir / ENCODE_CCRE_BED
            # Merge with other regulatory sources under one label
            reg_data: dict[str, list[tuple[int, int]]] = {}
            encode_data = BedUtils.load_bed(encode_bed_path)
            for chrom, intervals in encode_data.items():
                reg_data[chrom] = list(intervals)

            if preset.regulatory_fantom5:
                fantom5_data = BedUtils.load_bed(ref_dir / FANTOM5_BED)
                for chrom, intervals in fantom5_data.items():
                    if chrom not in reg_data:
                        reg_data[chrom] = []
                    reg_data[chrom].extend(intervals)

            if preset.regulatory_vista:
                vista_data = BedUtils.load_bed(ref_dir / VISTA_BED)
                for chrom, intervals in vista_data.items():
                    if chrom not in reg_data:
                        reg_data[chrom] = []
                    reg_data[chrom].extend(intervals)

            # Sort and merge all regulatory intervals
            for chrom in reg_data:
                reg_data[chrom].sort(key=lambda x: x[0])
                reg_data[chrom] = BedUtils.merge_intervals(reg_data[chrom])

            self.region_beds["regulatory"] = reg_data

        # Load safety net BEDs
        if preset.safetynet_clinvar:
            clinvar_path = ref_dir / CLINVAR_BED
            if clinvar_path.exists():
                self.safetynet_beds["ClinVar"] = BedUtils.load_bed(clinvar_path)
            else:
                self.safetynet_beds["ClinVar"] = {}

        if preset.safetynet_omim:
            omim_path = ref_dir / OMIM_BED
            if omim_path.exists():
                self.safetynet_beds["OMIM"] = BedUtils.load_bed(omim_path)
            else:
                self.safetynet_beds["OMIM"] = {}

    def annotate_vcf(self, input_vcf: Path, output_vcf: Path) -> FilterStats:
        """Annotate a VCF file with DGRA_REGION and DGRA_SAFETYNET INFO tags.

        Reads the input VCF line by line, adds INFO annotations for each
        variant, and writes the result to the output VCF.

        Args:
            input_vcf: Path to input VCF (may be gzipped).
            output_vcf: Path for annotated output VCF (gzipped if .vcf.gz).

        Returns:
            FilterStats with region_counts, clinvar_count, omim_count,
            region_only_variants, safetynet_only_variants, and
            region_and_safetynet_variants populated.
        """
        stats = FilterStats()
        region_counts: dict[str, int] = {"gene": 0, "ncrna": 0, "regulatory": 0}
        clinvar_count = 0
        omim_count = 0
        region_only = 0
        safetynet_only = 0
        region_and_safetynet = 0

        # Determine if we need to compress the output
        compress_output = str(output_vcf).endswith(".gz")

        # Determine input opener
        input_gz = str(input_vcf).endswith(".gz")

        output_vcf.parent.mkdir(parents=True, exist_ok=True)

        opener_in = gzip.open if input_gz else open

        if compress_output:
            out_f = gzip.open(output_vcf, "wt")
        else:
            out_f = open(output_vcf, "w")  # noqa: SIM115

        try:
            with opener_in(input_vcf, "rt") as in_f:
                for line in in_f:
                    if line.startswith("#"):
                        if line.startswith("##"):
                            # Insert our INFO headers before the first #CHROM line
                            # We'll add them when we see the last ## line before #CHROM
                            pass
                        if line.startswith("#CHROM"):
                            # Write the INFO headers before #CHROM
                            out_f.write(self.REGION_INFO_DEF + "\n")
                            out_f.write(self.SAFETYNET_INFO_DEF + "\n")
                            out_f.write(line)
                        else:
                            out_f.write(line)
                    else:
                        # Parse and annotate the variant record
                        annotated = self._annotate_record(line)

                        # Count regions and safety nets for statistics
                        parts = annotated.rstrip("\n").split("\t")
                        if len(parts) >= 8:
                            info = parts[7]
                            regions: list[str] = []
                            safetynets: list[str] = []

                            for field_str in info.split(";"):
                                if field_str.startswith("DGRA_REGION="):
                                    regions = field_str.split("=", 1)[1].split(",")
                                elif field_str.startswith("DGRA_SAFETYNET="):
                                    safetynets = field_str.split("=", 1)[1].split(",")

                            has_region = len(regions) > 0
                            has_safetynet = len(safetynets) > 0

                            if has_region and has_safetynet:
                                region_and_safetynet += 1
                            elif has_region:
                                region_only += 1
                            elif has_safetynet:
                                safetynet_only += 1

                            for r in regions:
                                if r in region_counts:
                                    region_counts[r] += 1

                            if "ClinVar" in safetynets:
                                clinvar_count += 1
                            if "OMIM" in safetynets:
                                omim_count += 1

                        out_f.write(annotated)
        finally:
            out_f.close()

        # Index the output if compressed
        if compress_output:
            try:
                import subprocess
                subprocess.run(
                    ["bcftools", "index", str(output_vcf)],
                    capture_output=True,
                    check=False,
                )
            except Exception:
                logger.debug("Could not index output VCF (non-critical)")

        stats.region_counts = region_counts
        stats.clinvar_count = clinvar_count
        stats.omim_count = omim_count
        stats.region_only_variants = region_only
        stats.safetynet_only_variants = safetynet_only
        stats.region_and_safetynet_variants = region_and_safetynet

        return stats

    def _determine_regions(self, chrom: str, pos: int) -> list[str]:
        """Determine which region categories a variant falls within.

        Args:
            chrom: Chromosome name (may or may not have 'chr' prefix).
            pos: 1-based VCF position.

        Returns:
            List of region labels (e.g., ['gene', 'regulatory']).
        """
        regions: list[str] = []
        for label, bed_data in self.region_beds.items():
            if BedUtils.query_overlap(bed_data, chrom, pos):
                regions.append(label)
        return regions

    def _determine_safetynets(self, chrom: str, pos: int) -> list[str]:
        """Determine which safety nets a variant falls within.

        Args:
            chrom: Chromosome name (may or may not have 'chr' prefix).
            pos: 1-based VCF position.

        Returns:
            List of safety net tags (e.g., ['ClinVar']).
        """
        safetynets: list[str] = []
        for tag, bed_data in self.safetynet_beds.items():
            if bed_data and BedUtils.query_overlap(bed_data, chrom, pos):
                safetynets.append(tag)
        return safetynets

    def _annotate_record(self, record: str) -> str:
        """Add DGRA_REGION and DGRA_SAFETYNET INFO tags to a VCF record.

        Args:
            record: A single VCF record line.

        Returns:
            The annotated VCF record line.
        """
        parts = record.rstrip("\n").split("\t")
        if len(parts) < 8:
            return record

        chrom = parts[0]
        try:
            pos = int(parts[1])
        except ValueError:
            return record + "\n"

        info = parts[7]

        # Determine region and safety net tags
        regions = self._determine_regions(chrom, pos)
        safetynets = self._determine_safetynets(chrom, pos)

        # Build new INFO fields
        new_info_parts: list[str] = []

        if regions:
            new_info_parts.append(f"DGRA_REGION={','.join(regions)}")
        if safetynets:
            new_info_parts.append(f"DGRA_SAFETYNET={','.join(safetynets)}")

        # Append to existing INFO
        if info == ".":
            parts[7] = ";".join(new_info_parts) if new_info_parts else "."
        else:
            if new_info_parts:
                parts[7] = info + ";" + ";".join(new_info_parts)

        return "\t".join(parts) + "\n"
