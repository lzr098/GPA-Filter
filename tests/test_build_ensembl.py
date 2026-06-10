"""Tests for Ensembl regulatory build in scripts/build_refs.py."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from dgra_prefilter.bed_utils import BedUtils


class TestBuildEnsemblRegulatoryBed:
    """Tests for build_ensembl_regulatory_bed()."""

    def test_parses_gff3_and_writes_bed(self, tmp_path: Path) -> None:
        """A minimal GFF3 is parsed into the expected 4-column BED."""
        from scripts.build_refs import build_ensembl_regulatory_bed

        gff3_path = tmp_path / "test.gff3.gz"
        lines = [
            "##gff-version 3",
            '1\tEnsembl\tpromoter\t101\t200\t.\t+\t.\tID=regulatory_1;type=promoter',
            '1\tEnsembl\tenhancer\t201\t300\t.\t+\t.\tID=regulatory_2;type=enhancer',
            '1\tEnsembl\tCTCF_binding_site\t301\t400\t.\t+\t.\tID=regulatory_3',
            '1\tEnsembl\topen_chromatin_region\t401\t500\t.\t+\t.\tID=regulatory_4',
            '1\tEnsembl\tTF_binding_site\t501\t600\t.\t+\t.\tID=regulatory_5',
            # Feature type not in the allowed set should be skipped
            '1\tEnsembl\tinsulator\t601\t700\t.\t+\t.\tID=regulatory_6',
        ]
        with gzip.open(gff3_path, "wt") as f:
            for line in lines:
                f.write(line + "\n")

        output_dir = tmp_path / "refs"
        output_dir.mkdir()
        build_ensembl_regulatory_bed(gff3_path, output_dir)

        bed_path = output_dir / "ensembl_regulatory_features.bed"
        assert bed_path.exists()

        data = BedUtils.load_bed_named(bed_path)
        assert "chr1" in data
        assert len(data["chr1"]) == 5

        # Verify 0-based conversion and type mapping
        assert data["chr1"][0] == (100, 200, "promoter")
        assert data["chr1"][1] == (200, 300, "enhancer")
        assert data["chr1"][2] == (300, 400, "ctcf")
        assert data["chr1"][3] == (400, 500, "open_chromatin")
        assert data["chr1"][4] == (500, 600, "tf_binding")

    def test_skips_unwanted_feature_types(self, tmp_path: Path) -> None:
        """Only allowed feature types are retained."""
        from scripts.build_refs import build_ensembl_regulatory_bed

        gff3_path = tmp_path / "test.gff3"
        with open(gff3_path, "w") as f:
            f.write("1\tEnsembl\tpromoter\t101\t200\t.\t+\t.\t.\n")
            f.write("1\tEnsembl\tinsulator\t201\t300\t.\t+\t.\t.\n")

        output_dir = tmp_path / "refs"
        output_dir.mkdir()
        build_ensembl_regulatory_bed(gff3_path, output_dir)

        bed_path = output_dir / "ensembl_regulatory_features.bed"
        data = BedUtils.load_bed_named(bed_path)
        assert len(data["chr1"]) == 1
        assert data["chr1"][0] == (100, 200, "promoter")

    def test_reads_type_from_attributes(self, tmp_path: Path) -> None:
        """Type can be extracted from GFF3 attributes when column 3 is generic."""
        from scripts.build_refs import build_ensembl_regulatory_bed

        gff3_path = tmp_path / "test.gff3"
        with open(gff3_path, "w") as f:
            f.write('1\tEnsembl\tregulatory_feature\t101\t200\t.\t+\t.\ttype=promoter\n')
            f.write('1\tEnsembl\tregulatory_feature\t201\t300\t.\t+\t.\ttype=enhancer\n')

        output_dir = tmp_path / "refs"
        output_dir.mkdir()
        build_ensembl_regulatory_bed(gff3_path, output_dir)

        bed_path = output_dir / "ensembl_regulatory_features.bed"
        data = BedUtils.load_bed_named(bed_path)
        assert len(data["chr1"]) == 2
        assert data["chr1"][0] == (100, 200, "promoter")
        assert data["chr1"][1] == (200, 300, "enhancer")
