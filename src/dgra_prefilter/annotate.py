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

    @staticmethod
    def _scan_chromosome(
        positions: list[int],
        bed_data: dict[str, list[tuple[int, int]]],
        chrom: str,
    ) -> list[bool]:
        """Two-pointer sweep to find which positions overlap BED intervals.

        Both positions and intervals must be sorted. Intervals are assumed
        non-overlapping (merged). Runs in O(n + m) where n = len(positions)
        and m = len(intervals).

        Args:
            positions: Sorted list of 1-based VCF positions.
            bed_data: Loaded BED dictionary {chrom: [(start, end), ...]}.
            chrom: Chromosome name (will be normalized).

        Returns:
            List of booleans, same length as positions, True if overlapping.
        """
        chrom = BedUtils.normalize_chrom(chrom)
        intervals = bed_data.get(chrom, [])
        if not intervals or not positions:
            return [False] * len(positions)

        n = len(positions)
        results = [False] * n
        vi = 0
        ii = 0
        m = len(intervals)

        while vi < n and ii < m:
            # Convert VCF 1-based to BED 0-based
            bed_pos = positions[vi] - 1
            i_start, i_end = intervals[ii]
            if bed_pos < i_start:
                vi += 1
            elif bed_pos >= i_end:
                ii += 1
            else:
                results[vi] = True
                vi += 1

        return results

    def _annotate_chromosome(
        self,
        lines: list[str],
        chrom: str,
    ) -> tuple[list[str], dict[str, int], dict[str, int], int, int, int]:
        """Annotate all variants for a single chromosome in one batch.

        Uses two-pointer sweep for O(n + m) performance per BED.
        Positions are sorted internally before scanning; results mapped
        back to original line order.

        Args:
            lines: List of raw VCF record lines for this chromosome.
            chrom: Chromosome name.

        Returns:
            Tuple of (annotated_lines, region_counts, clinvar_count,
            omim_count, region_only, safetynet_only, region_and_safetynet).
        """
        n = len(lines)
        # Extract positions and parsed parts; keep original index for mapping back
        indexed: list[tuple[int, int, list[str]]] = []
        for orig_idx, line in enumerate(lines):
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 8:
                try:
                    pos = int(parts[1])
                except ValueError:
                    pos = 0
            else:
                pos = 0
            indexed.append((pos, orig_idx, parts))

        # Sort by position for the sweep
        indexed.sort(key=lambda x: x[0])
        sorted_positions = [x[0] for x in indexed]

        # Determine region tags for each variant (on sorted order)
        sorted_region_tags: list[set[str]] = [set() for _ in range(n)]
        for label, bed in self.region_beds.items():
            hits = self._scan_chromosome(sorted_positions, bed, chrom)
            for sorted_idx, hit in enumerate(hits):
                if hit:
                    sorted_region_tags[sorted_idx].add(label)

        # Determine safety net tags for each variant (on sorted order)
        sorted_safetynet_tags: list[set[str]] = [set() for _ in range(n)]
        for label, bed in self.safetynet_beds.items():
            hits = self._scan_chromosome(sorted_positions, bed, chrom)
            for sorted_idx, hit in enumerate(hits):
                if hit:
                    sorted_safetynet_tags[sorted_idx].add(label)

        # Map tags back to original order
        region_tags: list[set[str]] = [set() for _ in range(n)]
        safetynet_tags: list[set[str]] = [set() for _ in range(n)]
        parsed_parts: list[list[str]] = [[] for _ in range(n)]
        for sorted_idx, (_, orig_idx, parts) in enumerate(indexed):
            region_tags[orig_idx] = sorted_region_tags[sorted_idx]
            safetynet_tags[orig_idx] = sorted_safetynet_tags[sorted_idx]
            parsed_parts[orig_idx] = parts

        # Build annotated lines and collect stats
        annotated: list[str] = []
        region_counts: dict[str, int] = {"gene": 0, "ncrna": 0, "regulatory": 0}
        clinvar_count = 0
        omim_count = 0
        region_only = 0
        safetynet_only = 0
        region_and_safetynet = 0

        for i in range(n):
            parts = parsed_parts[i]
            if not parts or len(parts) < 8:
                annotated.append(lines[i])
                continue

            info = parts[7]
            regions = region_tags[i]
            safetynets = safetynet_tags[i]

            # Build new INFO fields
            new_info_parts: list[str] = []
            if regions:
                new_info_parts.append(f"DGRA_REGION={','.join(sorted(regions))}")
            if safetynets:
                new_info_parts.append(f"DGRA_SAFETYNET={','.join(sorted(safetynets))}")

            if info == ".":
                parts[7] = ";".join(new_info_parts) if new_info_parts else "."
            else:
                if new_info_parts:
                    parts[7] = info + ";" + ";".join(new_info_parts)

            annotated.append("\t".join(parts) + "\n")

            # Stats
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

        return (
            annotated,
            region_counts,
            clinvar_count,
            omim_count,
            region_only,
            safetynet_only,
            region_and_safetynet,
        )

    def annotate_vcf(self, input_vcf: Path, output_vcf: Path) -> FilterStats:
        """Annotate a VCF file with DGRA_REGION and DGRA_SAFETYNET INFO tags.

        Processes the VCF chromosome-by-chromosome using two-pointer sweep
        for O(n + m) performance instead of per-variant bisect.

        Args:
            input_vcf: Path to input VCF (may be gzipped).
            output_vcf: Path for annotated output VCF (gzipped if .vcf.gz).

        Returns:
            FilterStats with region_counts, clinvar_count, omim_count,
            region_only_variants, safetynet_only_variants, and
            region_and_safetynet_variants populated.
        """
        stats = FilterStats()
        total_region_counts: dict[str, int] = {"gene": 0, "ncrna": 0, "regulatory": 0}
        total_clinvar = 0
        total_omim = 0
        total_region_only = 0
        total_safetynet_only = 0
        total_region_and_safetynet = 0

        compress_output = str(output_vcf).endswith(".gz")
        input_gz = str(input_vcf).endswith(".gz")
        output_vcf.parent.mkdir(parents=True, exist_ok=True)

        opener_in = gzip.open if input_gz else open

        if compress_output:
            out_f = gzip.open(output_vcf, "wt")
        else:
            out_f = open(output_vcf, "w")  # noqa: SIM115

        try:
            with opener_in(input_vcf, "rt") as in_f:
                chrom_buffer: list[str] = []
                current_chrom: str | None = None

                for line in in_f:
                    if line.startswith("#"):
                        if line.startswith("#CHROM"):
                            out_f.write(self.REGION_INFO_DEF + "\n")
                            out_f.write(self.SAFETYNET_INFO_DEF + "\n")
                        out_f.write(line)
                        continue

                    parts = line.rstrip("\n").split("\t")
                    chrom = parts[0] if parts else ""

                    if chrom != current_chrom and current_chrom is not None and chrom_buffer:
                        # Flush previous chromosome
                        (
                            annotated_lines,
                            rc,
                            cc,
                            oc,
                            ro,
                            so,
                            ras,
                        ) = self._annotate_chromosome(chrom_buffer, current_chrom)
                        for al in annotated_lines:
                            out_f.write(al)
                        for k, v in rc.items():
                            total_region_counts[k] += v
                        total_clinvar += cc
                        total_omim += oc
                        total_region_only += ro
                        total_safetynet_only += so
                        total_region_and_safetynet += ras
                        chrom_buffer = []

                    current_chrom = chrom
                    chrom_buffer.append(line)

                # Flush final chromosome
                if current_chrom is not None and chrom_buffer:
                    (
                        annotated_lines,
                        rc,
                        cc,
                        oc,
                        ro,
                        so,
                        ras,
                    ) = self._annotate_chromosome(chrom_buffer, current_chrom)
                    for al in annotated_lines:
                        out_f.write(al)
                    for k, v in rc.items():
                        total_region_counts[k] += v
                    total_clinvar += cc
                    total_omim += oc
                    total_region_only += ro
                    total_safetynet_only += so
                    total_region_and_safetynet += ras
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

        stats.region_counts = total_region_counts
        stats.clinvar_count = total_clinvar
        stats.omim_count = total_omim
        stats.region_only_variants = total_region_only
        stats.safetynet_only_variants = total_safetynet_only
        stats.region_and_safetynet_variants = total_region_and_safetynet

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
