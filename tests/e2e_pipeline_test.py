#!/usr/bin/env python3
"""End-to-end pipeline test for dgra-prefilter."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

# Ensure dgra_prefilter package is importable
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / "src"))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).parent.parent
TESTS_DIR = Path(__file__).parent

import dgra_prefilter  # noqa: E402

sys.path.insert(0, str(PROJECT_ROOT))
import importlib.util
_spec = importlib.util.spec_from_file_location("generate_demo_vcfs", str(TESTS_DIR / "generate_demo_vcfs.py"))
_generate_demo_vcfs_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_generate_demo_vcfs_mod)
generate_all_demo_vcfs = _generate_demo_vcfs_mod.generate_all_demo_vcfs

BCFTOOLS = os.environ.get(
    "BCFTOOLS", "/Users/zhaorongli/.micromamba/envs/bcftools/bin/bcftools"
)
DGRA_PREFILTER = os.environ.get(
    "DGRA_PREFILTER",
    shutil.which("dgra-prefilter")
    or "/Users/zhaorongli/.workbuddy/binaries/python/envs/default/bin/dgra-prefilter",
)
POSTFILTER = PROJECT_ROOT / "scripts" / "postfilter.py"
DIAGNOSE = PROJECT_ROOT / "scripts" / "diagnose_vcf.py"
PREPROCESS = PROJECT_ROOT / "scripts" / "preprocess_vcf.py"
REF_DIR = Path("~/.dgra-prefilter/refs").expanduser()

# Temp directories are created per-test-class to avoid sandbox restrictions.
DEMO_VCFS_DIR: Path | None = None
TMP_DIR: Path | None = None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def run_cmd(cmd: list[str], check: bool = True, cwd: Path | None = None) -> subprocess.CompletedProcess:
    """Run a shell command and return the result."""
    env = os.environ.copy()
    env["PATH"] = f"{os.path.dirname(BCFTOOLS)}:{env.get('PATH', '')}"
    env["BCFTOOLS_PATH"] = BCFTOOLS
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=cwd,
        env=env,
    )
    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode, cmd, output=result.stdout, stderr=result.stderr
        )
    return result


def count_vcf_variants(vcf_path: Path) -> int:
    """Count variants in a VCF using bcftools."""
    try:
        result = run_cmd([BCFTOOLS, "index", "-n", str(vcf_path)])
        return int(result.stdout.strip())
    except Exception:
        # Fallback
        import gzip

        opener = gzip.open if str(vcf_path).endswith(".gz") else open
        count = 0
        with opener(vcf_path, "rt") as f:
            for line in f:
                if not line.startswith("#"):
                    count += 1
        return count


def bcftools_index(vcf_path: Path) -> None:
    """Index a VCF with bcftools (force overwrite if index exists)."""
    run_cmd([BCFTOOLS, "index", "-f", str(vcf_path)])


def bgzip_recompress(vcf_path: Path) -> Path:
    """Recompress a gzip VCF to bgzip format for bcftools indexing."""
    tmp_path = vcf_path.with_suffix(".tmp.vcf.gz")
    result = subprocess.run(
        f"gunzip -c {vcf_path} | {BCFTOOLS} view -Oz -o {tmp_path} -",
        shell=True, capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"bgzip recompress failed: {result.stderr}")
    vcf_path.unlink()
    tmp_path.rename(vcf_path)
    bcftools_index(vcf_path)
    return vcf_path


def bcftools_view_check(vcf_path: Path) -> None:
    """Verify a VCF can be read by bcftools (works for plain gzip too)."""
    run_cmd([BCFTOOLS, "view", "-H", str(vcf_path)])


# ---------------------------------------------------------------------------
# Test cases
# ---------------------------------------------------------------------------
class TestE2EPipeline(unittest.TestCase):
    """End-to-end pipeline tests."""

    @classmethod
    def setUpClass(cls) -> None:
        """Generate demo VCFs once before all tests."""
        global DEMO_VCFS_DIR, TMP_DIR
        cls._demo_tmp = tempfile.TemporaryDirectory(prefix="dgra_e2e_demo_")
        cls._run_tmp = tempfile.TemporaryDirectory(prefix="dgra_e2e_run_")
        DEMO_VCFS_DIR = Path(cls._demo_tmp.name)
        TMP_DIR = Path(cls._run_tmp.name)
        DEMO_VCFS_DIR.mkdir(parents=True, exist_ok=True)
        TMP_DIR.mkdir(parents=True, exist_ok=True)
        cls.demo_paths = generate_all_demo_vcfs(DEMO_VCFS_DIR, ref_dir=REF_DIR)

    @classmethod
    def tearDownClass(cls) -> None:
        """Clean up temporary directories."""
        cls._demo_tmp.cleanup()
        cls._run_tmp.cleanup()

    def setUp(self) -> None:
        """Clean temp files for each test."""
        for f in TMP_DIR.iterdir():
            if f.is_file():
                f.unlink()

    # -----------------------------------------------------------------------
    # full_grch38
    # -----------------------------------------------------------------------
    def test_full_grch38_pipeline(self) -> None:
        """Test full pipeline with full_grch38 (raw_caller_output)."""
        vcf = self.demo_paths["full_grch38"]
        diagnose_report = TMP_DIR / "full_grch38_diagnose.json"
        preprocessed = TMP_DIR / "full_grch38_preprocessed.vcf.gz"
        filtered = TMP_DIR / "full_grch38_filtered.vcf.gz"
        postfiltered = TMP_DIR / "full_grch38_postfiltered.vcf.gz"
        annotated = TMP_DIR / "full_grch38_annotated.vcf.gz"

        # Step 1: diagnose
        run_cmd([
            sys.executable, str(DIAGNOSE),
            "-i", str(vcf),
            "-o", str(diagnose_report),
        ])
        with open(diagnose_report) as f:
            report = json.load(f)
        self.assertEqual(report["quality_assessment"], "raw_caller_output")
        self.assertEqual(report["inferred_genome"], "GRCh38")

        input_count = count_vcf_variants(vcf)

        # Step 2: preprocess
        run_cmd([
            sys.executable, str(PREPROCESS),
            "-i", str(vcf),
            "--diagnose-report", str(diagnose_report),
            "-o", str(preprocessed),
        ])
        preproc_count = count_vcf_variants(preprocessed)
        removed_pct = (input_count - preproc_count) / input_count * 100
        # Assertion 1: preprocess removes ~35% of variants
        self.assertGreater(removed_pct, 25.0, f"Expected >25% removed, got {removed_pct:.1f}%")
        self.assertLess(removed_pct, 50.0, f"Expected <50% removed, got {removed_pct:.1f}%")

        # Step 3: dgra-prefilter
        run_cmd([
            DGRA_PREFILTER,
            "-i", str(preprocessed),
            "-o", str(filtered),
            "-g", "GRCh38",
            "-p", "comprehensive",
            "--ref-dir", str(REF_DIR),
        ])
        filtered_count = count_vcf_variants(filtered)
        # Assertion 2: prefilter retains >0 variants
        self.assertGreater(filtered_count, 0, "Prefilter should retain >0 variants")

        # Step 4: postfilter
        run_cmd([
            sys.executable, str(POSTFILTER),
            "-i", str(filtered),
            "-o", str(postfiltered),
            "--ref-dir", str(REF_DIR),
            "--preset", "comprehensive",
        ])
        postfiltered_count = count_vcf_variants(postfiltered)
        # Assertion 3: postfilter further reduces
        self.assertLess(
            postfiltered_count,
            filtered_count,
            "Postfilter should further reduce variant count",
        )

        # Step 5: annotate
        dgra_prefilter.annotate_vcf_file(
            input_path=postfiltered,
            output_path=annotated,
            preset="comprehensive",
            ref_dir=REF_DIR,
        )

        # Assertion 8: output VCFs valid (bcftools can read/index them)
        bcftools_index(preprocessed)
        bcftools_index(filtered)
        bcftools_view_check(postfiltered)
        bcftools_view_check(annotated)

    # -----------------------------------------------------------------------
    # full_grch37
    # -----------------------------------------------------------------------
    def test_full_grch37_preprocess_exits(self) -> None:
        """Test that full_grch37 preprocess exits with code 1 (needs_liftover)."""
        vcf = self.demo_paths["full_grch37"]
        diagnose_report = TMP_DIR / "full_grch37_diagnose.json"
        preprocessed = TMP_DIR / "full_grch37_preprocessed.vcf.gz"

        # Step 1: diagnose
        run_cmd([
            sys.executable, str(DIAGNOSE),
            "-i", str(vcf),
            "-o", str(diagnose_report),
        ])
        with open(diagnose_report) as f:
            report = json.load(f)
        self.assertEqual(report["quality_assessment"], "needs_liftover")
        self.assertEqual(report["inferred_genome"], "GRCh37")

        # Step 2: preprocess should exit with code 1
        # Assertion 4: preprocess exits code 1
        with self.assertRaises(subprocess.CalledProcessError) as cm:
            run_cmd([
                sys.executable, str(PREPROCESS),
                "-i", str(vcf),
                "--diagnose-report", str(diagnose_report),
                "-o", str(preprocessed),
            ])
        self.assertEqual(cm.exception.returncode, 1)

    # -----------------------------------------------------------------------
    # genotyped_grch38
    # -----------------------------------------------------------------------
    def test_genotyped_grch38_pipeline(self) -> None:
        """Test full pipeline with genotyped_grch38 (clean_genotyped)."""
        vcf = self.demo_paths["genotyped_grch38"]
        diagnose_report = TMP_DIR / "genotyped_grch38_diagnose.json"
        preprocessed = TMP_DIR / "genotyped_grch38_preprocessed.vcf.gz"
        filtered = TMP_DIR / "genotyped_grch38_filtered.vcf.gz"
        postfiltered = TMP_DIR / "genotyped_grch38_postfiltered.vcf.gz"
        annotated = TMP_DIR / "genotyped_grch38_annotated.vcf.gz"

        # Step 1: diagnose
        run_cmd([
            sys.executable, str(DIAGNOSE),
            "-i", str(vcf),
            "-o", str(diagnose_report),
        ])
        with open(diagnose_report) as f:
            report = json.load(f)
        self.assertEqual(report["quality_assessment"], "clean_genotyped")
        self.assertEqual(report["inferred_genome"], "GRCh38")

        input_count = count_vcf_variants(vcf)

        # Step 2: preprocess (copy through)
        run_cmd([
            sys.executable, str(PREPROCESS),
            "-i", str(vcf),
            "--diagnose-report", str(diagnose_report),
            "-o", str(preprocessed),
        ])
        preproc_count = count_vcf_variants(preprocessed)
        # Assertion 5: preprocess copies with no loss
        self.assertEqual(
            preproc_count,
            input_count,
            "Clean genotyped VCF should not lose variants during preprocess",
        )

        # Step 3: dgra-prefilter
        run_cmd([
            DGRA_PREFILTER,
            "-i", str(preprocessed),
            "-o", str(filtered),
            "-g", "GRCh38",
            "-p", "comprehensive",
            "--ref-dir", str(REF_DIR),
        ])
        filtered_count = count_vcf_variants(filtered)
        self.assertGreater(filtered_count, 0, "Prefilter should retain >0 variants")

        # Step 4: postfilter
        run_cmd([
            sys.executable, str(POSTFILTER),
            "-i", str(filtered),
            "-o", str(postfiltered),
            "--ref-dir", str(REF_DIR),
            "--preset", "comprehensive",
        ])
        postfiltered_count = count_vcf_variants(postfiltered)

        # Step 5: annotate
        dgra_prefilter.annotate_vcf_file(
            input_path=postfiltered,
            output_path=annotated,
            preset="comprehensive",
            ref_dir=REF_DIR,
        )

        # Assertion 6: full pipeline completes (implicitly true if we get here)
        # Assertion 8: output VCFs valid
        bcftools_index(preprocessed)
        bcftools_index(filtered)
        bcftools_view_check(postfiltered)
        bcftools_view_check(annotated)

    # -----------------------------------------------------------------------
    # genotyped_grch37
    # -----------------------------------------------------------------------
    def test_genotyped_grch37_preprocess_exits(self) -> None:
        """Test that genotyped_grch37 preprocess exits with code 1."""
        vcf = self.demo_paths["genotyped_grch37"]
        diagnose_report = TMP_DIR / "genotyped_grch37_diagnose.json"
        preprocessed = TMP_DIR / "genotyped_grch37_preprocessed.vcf.gz"

        # Step 1: diagnose
        run_cmd([
            sys.executable, str(DIAGNOSE),
            "-i", str(vcf),
            "-o", str(diagnose_report),
        ])
        with open(diagnose_report) as f:
            report = json.load(f)
        self.assertEqual(report["quality_assessment"], "needs_liftover")
        self.assertEqual(report["inferred_genome"], "GRCh37")

        # Assertion 7: preprocess exits code 1
        with self.assertRaises(subprocess.CalledProcessError) as cm:
            run_cmd([
                sys.executable, str(PREPROCESS),
                "-i", str(vcf),
                "--diagnose-report", str(diagnose_report),
                "-o", str(preprocessed),
            ])
        self.assertEqual(cm.exception.returncode, 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
