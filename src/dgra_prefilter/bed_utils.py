"""BED file utility functions for dgra-prefilter."""

from __future__ import annotations

import bisect
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


class BedUtils:
    """Static utility class for BED file I/O, interval merging, and overlap queries."""

    @staticmethod
    def normalize_chrom(chrom: str) -> str:
        """Normalize chromosome name to include 'chr' prefix.

        Args:
            chrom: Chromosome name (e.g., '1', 'chr1', 'X', 'chrX').

        Returns:
            Chromosome name with 'chr' prefix (e.g., 'chr1', 'chrX').
        """
        if not chrom.startswith("chr"):
            return f"chr{chrom}"
        return chrom

    @staticmethod
    def load_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
        """Load a BED file into a dictionary keyed by chromosome.

        Each chromosome maps to a sorted list of (start, end) tuples.
        Automatically normalizes chromosome names to include 'chr' prefix.
        Skips comment lines ('#') and track header lines.

        Args:
            path: Path to the BED file.

        Returns:
            Dictionary {chrom: [(start, end), ...]} with intervals sorted
            by start coordinate.

        Raises:
            FileNotFoundError: If the BED file does not exist.
        """
        if not path.exists():
            raise FileNotFoundError(f"BED file not found: {path}")

        data: dict[str, list[tuple[int, int]]] = {}
        with open(path, "r") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or stripped.startswith("track"):
                    continue
                parts = stripped.split("\t")
                if len(parts) < 3:
                    # Try space-separated as fallback
                    parts = stripped.split()
                if len(parts) < 3:
                    continue
                chrom = BedUtils.normalize_chrom(parts[0])
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    logger.warning("Skipping invalid BED line: %s", stripped)
                    continue
                if chrom not in data:
                    data[chrom] = []
                data[chrom].append((start, end))

        # Sort intervals by start coordinate within each chromosome
        for chrom in data:
            data[chrom].sort(key=lambda x: x[0])

        return data

    @staticmethod
    def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
        """Merge overlapping and adjacent intervals using a greedy algorithm.

        Input intervals must be sorted by start coordinate. Two intervals
        are merged if they overlap or are adjacent (gap = 0).

        Args:
            intervals: Sorted list of (start, end) tuples.

        Returns:
            List of merged (start, end) tuples with no overlaps.
        """
        if not intervals:
            return []

        merged: list[tuple[int, int]] = [intervals[0]]
        for start, end in intervals[1:]:
            last_start, last_end = merged[-1]
            if start <= last_end:
                # Overlapping or adjacent; extend the current interval
                merged[-1] = (last_start, max(last_end, end))
            else:
                merged.append((start, end))

        return merged

    @staticmethod
    def merge_bed_files(paths: list[Path], output: Path) -> Path:
        """Merge multiple BED files into a single file with overlapping intervals merged.

        Reads all files, collects intervals by chromosome, merges overlapping
        intervals within each chromosome, and writes the sorted result.

        Args:
            paths: List of BED file paths to merge.
            output: Output path for the merged BED file.

        Returns:
            Path to the merged BED file.

        Raises:
            FileNotFoundError: If any input BED file does not exist.
        """
        all_intervals: dict[str, list[tuple[int, int]]] = {}

        for bed_path in paths:
            data = BedUtils.load_bed(bed_path)
            for chrom, intervals in data.items():
                if chrom not in all_intervals:
                    all_intervals[chrom] = []
                all_intervals[chrom].extend(intervals)

        # Sort and merge per chromosome
        for chrom in all_intervals:
            all_intervals[chrom].sort(key=lambda x: x[0])
            all_intervals[chrom] = BedUtils.merge_intervals(all_intervals[chrom])

        BedUtils.write_bed(all_intervals, output)
        return output

    @staticmethod
    def query_overlap(
        bed_data: dict[str, list[tuple[int, int]]],
        chrom: str,
        pos: int,
    ) -> list[tuple[int, int]]:
        """Query whether a genomic position overlaps with any interval in bed_data.

        Uses bisect for O(log n) lookup. The position is treated as a VCF
        1-based coordinate, which is converted to BED 0-based half-open
        for comparison: bed_start <= (pos - 1) < bed_end.

        Args:
            bed_data: Dictionary {chrom: [(start, end), ...]} from load_bed().
            chrom: Chromosome name (will be normalized).
            pos: 1-based VCF position.

        Returns:
            List of (start, end) tuples that contain the position.
        """
        chrom = BedUtils.normalize_chrom(chrom)
        if chrom not in bed_data:
            return []

        intervals = bed_data[chrom]
        if not intervals:
            return []

        # Convert VCF 1-based to BED 0-based
        bed_pos = pos - 1

        # Use bisect to find the insertion point
        # We search for intervals where start <= bed_pos
        starts = [iv[0] for iv in intervals]
        idx = bisect.bisect_right(starts, bed_pos)

        # Check intervals from idx-1 backwards for overlap
        results: list[tuple[int, int]] = []
        for i in range(idx - 1, -1, -1):
            start, end = intervals[i]
            if start <= bed_pos < end:
                results.append((start, end))
            elif end <= bed_pos:
                # Since intervals are sorted by start and non-overlapping after merge,
                # once we find an interval that ends before our position, we can stop
                break

        return results

    @staticmethod
    def write_bed(
        intervals: dict[str, list[tuple[int, int]]],
        output: Path,
    ) -> None:
        """Write intervals to a BED file.

        Intervals are sorted by chromosome (natural order) and then by
        start coordinate within each chromosome.

        Args:
            intervals: Dictionary {chrom: [(start, end), ...]}.
            output: Output file path.
        """
        output.parent.mkdir(parents=True, exist_ok=True)

        # Sort chromosomes in natural order
        def chrom_sort_key(c: str) -> tuple[int, str]:
            """Sort key for chromosome names."""
            name = c.replace("chr", "")
            if name.isdigit():
                return (0, name.zfill(2))
            elif name == "X":
                return (1, "X")
            elif name == "Y":
                return (2, "Y")
            elif name == "M":
                return (3, "M")
            else:
                return (4, name)

        sorted_chroms = sorted(intervals.keys(), key=chrom_sort_key)

        with open(output, "w") as f:
            for chrom in sorted_chroms:
                for start, end in intervals[chrom]:
                    f.write(f"{chrom}\t{start}\t{end}\n")
