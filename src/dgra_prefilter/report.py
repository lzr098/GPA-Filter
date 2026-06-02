"""JSON report generator for dgra-prefilter."""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path

from dgra_prefilter.core import FilterStats, PrefilterConfig

logger = logging.getLogger(__name__)

_TOOL_NAME = "dgra-prefilter"
_TOOL_VERSION = "1.0.0"


class ReportGenerator:
    """Generates JSON reports for prefilter runs."""

    @staticmethod
    def generate(
        stats: FilterStats,
        config: PrefilterConfig,
        output_path: Path,
    ) -> Path:
        """Generate a structured JSON report for a prefilter run.

        The report includes input/output statistics, region counts,
        safety net hits, reference data versions, and timing information.

        Args:
            stats: Filter statistics from the run.
            config: Prefilter configuration used for the run.
            output_path: Path for the JSON report file.

        Returns:
            Path to the generated report file.
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)

        report = {
            "tool": _TOOL_NAME,
            "version": _TOOL_VERSION,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "input": {
                "path": str(config.input_path),
                "genome": config.genome,
                "variants": stats.input_variants,
            },
            "output": {
                "path": str(config.output_path),
                "variants": stats.retained_variants,
                "retention_rate": round(stats.retention_rate, 6),
            },
            "regions": {
                "gene": stats.region_counts.get("gene", 0),
                "ncrna": stats.region_counts.get("ncrna", 0),
                "regulatory": stats.region_counts.get("regulatory", 0),
                "region_only": stats.region_only_variants,
                "safetynet_only": stats.safetynet_only_variants,
                "region_and_safetynet": stats.region_and_safetynet_variants,
            },
            "safety_net": {
                "ClinVar": stats.clinvar_count,
                "OMIM": stats.omim_count,
            },
            "preset": stats.preset_config or {
                "name": config.preset_name,
            },
            "ref_data_versions": stats.ref_data_versions,
            "elapsed_seconds": round(stats.elapsed_seconds, 3),
        }

        with open(output_path, "w") as f:
            json.dump(report, f, indent=2, ensure_ascii=False)

        logger.info("Report written to %s", output_path)
        return output_path
