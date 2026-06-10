"""VCF INFO field annotator for dgra-prefilter."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

from dgra_prefilter.bed_utils import BedUtils
from dgra_prefilter.constants import (
    CLINVAR_BED,
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
    OMIM_BED,
    REGION_TAG_GENE_3UTR,
    REGION_TAG_GENE_5UTR,
    REGION_TAG_GENE_CDS,
    REGION_TAG_GENE_SPLICE,
    REGION_TAG_REGULATORY_CTCF,
    REGION_TAG_REGULATORY_DELS,
    REGION_TAG_REGULATORY_ENHANCER,
    REGION_TAG_REGULATORY_OPEN_CHROMATIN,
    REGION_TAG_REGULATORY_PELS,
    REGION_TAG_REGULATORY_PLS,
    REGION_TAG_REGULATORY_PROMOTER,
    REGION_TAG_REGULATORY_TF_BINDING,
    VISTA_BED,
)
from dgra_prefilter.core import FilterStats
from dgra_prefilter.presets import PresetConfig
from dgra_prefilter.safetynet import ClinVarSafetyNet

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
        # Summary region BEDs (gene, ncrna, regulatory)
        self.region_beds: dict[str, dict[str, list[tuple[int, int]]]] = {}
        # Simple sub-region BEDs where the dict key is the full tag
        # (e.g., gene_5utr, gene_cds). These BEDs have only 3 columns.
        self.region_simple_sub_beds: dict[
            str, dict[str, list[tuple[int, int]]]
        ] = {}
        # Named sub-region BEDs where the 4th column carries the sub-tag
        # (e.g., regulatory types PLS/pELS/dELS/CTCF).
        self.region_named_sub_beds: dict[
            str, dict[str, list[tuple[int, int, str]]]
        ] = {}
        self.safetynet_beds: dict[str, dict[str, list[tuple[int, int]]]] = {}
        self.safetynet_sub_beds: dict[
            str, dict[str, list[tuple[int, int, str]]]
        ] = {}

    def load_beds(self, preset: PresetConfig, ref_dir: Path) -> None:
        """Load all BED files needed for annotation based on the preset.

        Region BEDs are mapped to category labels (gene, ncrna, regulatory).
        Sub-region BEDs with a 4th name column are loaded for detailed
        annotation (e.g., gene_5utr, regulatory_pls). Safety net BEDs are
        mapped by their tag names (ClinVar, OMIM).

        Args:
            preset: The preset configuration.
            ref_dir: Directory containing reference BED files.
        """
        ref_dir = ref_dir.expanduser().resolve()

        # Load region BEDs
        if preset.gene:
            # Coding-only filters by sub-regions; for the summary
            # "gene" tag we still load the full gene loci BED when
            # available so backward-compatible output is preserved.
            if preset.name == "coding-only":
                gene_bed_path = ref_dir / GENCODE_GENE_BED
                if gene_bed_path.exists():
                    self.region_beds["gene"] = BedUtils.load_bed(gene_bed_path)

            # Always load the full gene loci BED for summary annotation
            # when gene regions are enabled.
            if "gene" not in self.region_beds:
                self.region_beds["gene"] = BedUtils.load_bed(
                    ref_dir / GENCODE_GENE_BED
                )

            # Load sub-region BEDs if they exist (used by coding-only and
            # for richer annotation in other presets).
            sub_beds = {
                REGION_TAG_GENE_5UTR: GENCODE_5UTR_BED,
                REGION_TAG_GENE_CDS: GENCODE_CDS_BED,
                REGION_TAG_GENE_3UTR: GENCODE_3UTR_BED,
                REGION_TAG_GENE_SPLICE: GENCODE_SPLICE_BED,
            }
            for tag, bed_name in sub_beds.items():
                bed_path = ref_dir / bed_name
                if bed_path.exists():
                    self.region_simple_sub_beds[tag] = BedUtils.load_bed(
                        bed_path
                    )

        if preset.ncrna:
            ncrna_bed_path = ref_dir / GENCODE_NCRNA_BED
            self.region_beds["ncrna"] = BedUtils.load_bed(ncrna_bed_path)

        # Regulatory regions
        if preset.regulatory_encode:
            if preset.regulatory_encode_balanced:
                encode_bed_path = ref_dir / ENCODE_BALANCED_BED
                encode_data = BedUtils.load_bed_named(encode_bed_path)
                # Sub-tags come directly from the 4th column
                self.region_named_sub_beds["regulatory"] = encode_data
                # Summary regulatory BED merged from all intervals
                reg_data: dict[str, list[tuple[int, int]]] = {}
                for chrom, intervals in encode_data.items():
                    reg_data[chrom] = [(s, e) for s, e, _ in intervals]
                for chrom in reg_data:
                    reg_data[chrom].sort(key=lambda x: x[0])
                    reg_data[chrom] = BedUtils.merge_intervals(reg_data[chrom])
                self.region_beds["regulatory"] = reg_data
            elif preset.regulatory_encode_pls_pels_only:
                encode_bed_path = ref_dir / ENCODE_PLS_PELS_BED
                self.region_beds["regulatory"] = BedUtils.load_bed(
                    encode_bed_path
                )
            else:
                encode_bed_path = ref_dir / ENCODE_CCRE_BED
                # Merge with other regulatory sources under one label
                reg_data = {}
                encode_data = BedUtils.load_bed(encode_bed_path)
                for chrom, intervals in encode_data.items():
                    reg_data[chrom] = list(intervals)

                if preset.regulatory_fantom5:
                    if preset.regulatory_source == "fantom5":
                        fantom5_data = BedUtils.load_bed(ref_dir / FANTOM5_BED)
                        for chrom, intervals in fantom5_data.items():
                            if chrom not in reg_data:
                                reg_data[chrom] = []
                            reg_data[chrom].extend(intervals)
                    elif preset.regulatory_source == "ensembl":
                        ensembl_path = ref_dir / ENSEMBL_REGULATORY_BED
                        if ensembl_path.exists():
                            ensembl_data = BedUtils.load_bed_named(ensembl_path)
                            self.region_named_sub_beds["regulatory"] = ensembl_data
                            for chrom, intervals in ensembl_data.items():
                                if chrom not in reg_data:
                                    reg_data[chrom] = []
                                reg_data[chrom].extend([(s, e) for s, e, _ in intervals])
                    elif preset.regulatory_source == "both":
                        fantom5_data = BedUtils.load_bed(ref_dir / FANTOM5_BED)
                        for chrom, intervals in fantom5_data.items():
                            if chrom not in reg_data:
                                reg_data[chrom] = []
                            reg_data[chrom].extend(intervals)
                        ensembl_path = ref_dir / ENSEMBL_REGULATORY_BED
                        if ensembl_path.exists():
                            ensembl_data = BedUtils.load_bed_named(ensembl_path)
                            # Merge Ensembl named sub-beds with any existing
                            if "regulatory" not in self.region_named_sub_beds:
                                self.region_named_sub_beds["regulatory"] = {}
                            for chrom, intervals in ensembl_data.items():
                                if chrom not in self.region_named_sub_beds["regulatory"]:
                                    self.region_named_sub_beds["regulatory"][chrom] = []
                                self.region_named_sub_beds["regulatory"][chrom].extend(intervals)
                            for chrom in self.region_named_sub_beds["regulatory"]:
                                self.region_named_sub_beds["regulatory"][chrom].sort(key=lambda x: x[0])
                            for chrom, intervals in ensembl_data.items():
                                if chrom not in reg_data:
                                    reg_data[chrom] = []
                                reg_data[chrom].extend([(s, e) for s, e, _ in intervals])

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
                try:
                    self.safetynet_sub_beds["ClinVar"] = (
                        BedUtils.load_bed_named(clinvar_path)
                    )
                except (FileNotFoundError, ValueError):
                    # Fall back to 3-column BED for backward compatibility
                    self.safetynet_beds["ClinVar"] = BedUtils.load_bed(
                        clinvar_path
                    )
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

    @staticmethod
    def _scan_chromosome_named(
        positions: list[int],
        bed_data: dict[str, list[tuple[int, int, str]]],
        chrom: str,
    ) -> list[set[str]]:
        """Two-pointer sweep returning the set of interval names per position.

        Args:
            positions: Sorted list of 1-based VCF positions.
            bed_data: Loaded named BED dictionary {chrom: [(start, end, name), ...]}.
            chrom: Chromosome name (will be normalized).

        Returns:
            List of sets of interval names, same length as positions.
        """
        chrom = BedUtils.normalize_chrom(chrom)
        intervals = bed_data.get(chrom, [])
        n = len(positions)
        results: list[set[str]] = [set() for _ in range(n)]
        if not intervals or not positions:
            return results

        vi = 0
        ii = 0
        m = len(intervals)

        while vi < n and ii < m:
            bed_pos = positions[vi] - 1
            i_start, i_end, i_name = intervals[ii]
            if bed_pos < i_start:
                vi += 1
            elif bed_pos >= i_end:
                ii += 1
            else:
                results[vi].add(i_name)
                # A position can overlap multiple named intervals;
                # advance variant index but keep interval index.
                vi += 1

        return results

    def _map_sub_tags_to_summary(self, sub_tag: str) -> str | None:
        """Map a sub-region tag to its summary tag.

        Args:
            sub_tag: Sub-region tag like 'gene_5utr' or 'regulatory_pls'.

        Returns:
            Summary tag ('gene', 'regulatory') or None.
        """
        if sub_tag.startswith("gene_"):
            return "gene"
        if sub_tag.startswith("regulatory_"):
            return "regulatory"
        return None

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

        # Determine simple sub-region tags (e.g., gene_5utr) and ensure
        # summary tags are present.
        for label, bed in self.region_simple_sub_beds.items():
            hits = self._scan_chromosome(sorted_positions, bed, chrom)
            for sorted_idx, hit in enumerate(hits):
                if hit:
                    sorted_region_tags[sorted_idx].add(label)
                    summary = self._map_sub_tags_to_summary(label)
                    if summary is not None:
                        sorted_region_tags[sorted_idx].add(summary)

        # Determine named sub-region tags (e.g., regulatory_pls). The 4th
        # BED column carries the sub-type.
        for label, bed in self.region_named_sub_beds.items():
            hits = self._scan_chromosome_named(sorted_positions, bed, chrom)
            for sorted_idx, names in enumerate(hits):
                if not names:
                    continue
                for name in names:
                    sorted_region_tags[sorted_idx].add(
                        f"{label}_{name.lower()}"
                    )
                summary = self._map_sub_tags_to_summary(label)
                if summary is not None:
                    sorted_region_tags[sorted_idx].add(summary)

        # Determine safety net tags for each variant (on sorted order)
        sorted_safetynet_tags: list[set[str]] = [set() for _ in range(n)]
        for label, bed in self.safetynet_beds.items():
            hits = self._scan_chromosome(sorted_positions, bed, chrom)
            for sorted_idx, hit in enumerate(hits):
                if hit:
                    sorted_safetynet_tags[sorted_idx].add(label)

        # Determine star-level safety net tags from named BEDs
        for label, bed in self.safetynet_sub_beds.items():
            hits = self._scan_chromosome_named(sorted_positions, bed, chrom)
            for sorted_idx, names in enumerate(hits):
                for name in names:
                    if label == "ClinVar":
                        sorted_safetynet_tags[sorted_idx].add(
                            ClinVarSafetyNet.format_tag(name)
                        )
                    else:
                        sorted_safetynet_tags[sorted_idx].add(name)

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
        region_counts: dict[str, int] = {
            "gene": 0,
            "ncrna": 0,
            "regulatory": 0,
            REGION_TAG_GENE_5UTR: 0,
            REGION_TAG_GENE_CDS: 0,
            REGION_TAG_GENE_3UTR: 0,
            REGION_TAG_GENE_SPLICE: 0,
            REGION_TAG_REGULATORY_PLS: 0,
            REGION_TAG_REGULATORY_PELS: 0,
            REGION_TAG_REGULATORY_DELS: 0,
            REGION_TAG_REGULATORY_CTCF: 0,
            REGION_TAG_REGULATORY_PROMOTER: 0,
            REGION_TAG_REGULATORY_ENHANCER: 0,
            REGION_TAG_REGULATORY_OPEN_CHROMATIN: 0,
            REGION_TAG_REGULATORY_TF_BINDING: 0,
        }
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
                new_info_parts.append(
                    f"DGRA_SAFETYNET={','.join(sorted(safetynets))}"
                )

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

            # Count ClinVar hits including star-level tags
            has_clinvar = any(
                t == "ClinVar" or t.startswith("ClinVar_")
                for t in safetynets
            )
            if has_clinvar:
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
        total_region_counts: dict[str, int] = {
            "gene": 0,
            "ncrna": 0,
            "regulatory": 0,
            REGION_TAG_GENE_5UTR: 0,
            REGION_TAG_GENE_CDS: 0,
            REGION_TAG_GENE_3UTR: 0,
            REGION_TAG_GENE_SPLICE: 0,
            REGION_TAG_REGULATORY_PLS: 0,
            REGION_TAG_REGULATORY_PELS: 0,
            REGION_TAG_REGULATORY_DELS: 0,
            REGION_TAG_REGULATORY_CTCF: 0,
            REGION_TAG_REGULATORY_PROMOTER: 0,
            REGION_TAG_REGULATORY_ENHANCER: 0,
            REGION_TAG_REGULATORY_OPEN_CHROMATIN: 0,
            REGION_TAG_REGULATORY_TF_BINDING: 0,
        }
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
        for label, bed_data in self.region_simple_sub_beds.items():
            if BedUtils.query_overlap(bed_data, chrom, pos):
                if label not in regions:
                    regions.append(label)
                summary = self._map_sub_tags_to_summary(label)
                if summary is not None and summary not in regions:
                    regions.append(summary)
        for label, bed_data in self.region_named_sub_beds.items():
            overlaps = BedUtils.query_overlap_named(bed_data, chrom, pos)
            if not overlaps:
                continue
            for _, _, name in overlaps:
                sub_tag = f"{label}_{name.lower()}"
                if sub_tag not in regions:
                    regions.append(sub_tag)
            summary = self._map_sub_tags_to_summary(label)
            if summary is not None and summary not in regions:
                regions.append(summary)
        return regions

    def _determine_safetynets(self, chrom: str, pos: int) -> list[str]:
        """Determine which safety nets a variant falls within.

        Args:
            chrom: Chromosome name (may or may not have 'chr' prefix).
            pos: 1-based VCF position.

        Returns:
            List of safety net tags (e.g., ['ClinVar_3star']).
        """
        safetynets: list[str] = []
        for tag, bed_data in self.safetynet_beds.items():
            if bed_data and BedUtils.query_overlap(bed_data, chrom, pos):
                safetynets.append(tag)
        for tag, bed_data in self.safetynet_sub_beds.items():
            if bed_data:
                overlaps = BedUtils.query_overlap_named(bed_data, chrom, pos)
                for _, _, name in overlaps:
                    if tag == "ClinVar":
                        safetynets.append(ClinVarSafetyNet.format_tag(name))
                    else:
                        safetynets.append(name)
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
