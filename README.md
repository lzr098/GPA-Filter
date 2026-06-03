# dgra-prefilter

Genomic region prefilter for whole-genome VCF files based on GENCODE, ENCODE, and ClinVar.

## What it does

dgra-prefilter is a **static position mapping table** — it pre-defines which genomic regions are "interesting" so you can rapidly filter WGS VCF files to only those variants falling within biologically significant positions. It does NOT perform variant annotation; it simply keeps or discards variants based on their genomic coordinates.

The tool builds BED files from three data sources, then uses `bcftools -T` to filter VCF variants by position:

![GPA Filter — genomic region coverage](docs/genomic-region-coverage.svg)

### Three data layers

| Layer | Source | What it covers |
|-------|--------|---------------|
| **Gene loci** | GENCODE v44 GTF | All transcript extents (min→max coordinates) for protein-coding genes and ncRNA |
| **Regulatory elements** | ENCODE SCREEN v3 | cCREs — promoter-like sequences (PLS), proximal/pDistal enhancer-like sequences (pELS/dELS) |
| **Safety net** | ClinVar | All pathogenic/likely pathogenic (P/LP) variants — guarantees zero omission of known pathogenic variants, even if they fall outside gene loci or regulatory regions |

### How it works

1. **Build phase** (`build_refs.py`): Downloads GENCODE GTF, ENCODE SCREEN BED, and ClinVar VCF → parses and converts each into sorted BED files → merges into preset-specific composite BED files
2. **Filter phase** (CLI/API): Loads the appropriate BED for your chosen preset → runs `bcftools view -T` for coordinate-based filtering → merges in ClinVar safety net hits → annotates retained variants with `PREFILTER` INFO tag

The BED reference files are **static** — once built, they don't change until you manually re-run `build_refs.py` with updated source data. No live API calls during filtering.

## Installation

```bash
pip install dgra-prefilter
```

## Usage

### CLI

```bash
dgra-prefilter \
  --input sample.vcf.gz \
  --output filtered.vcf.gz \
  --preset comprehensive
```

### Python API

```python
from dgra_prefilter import prefilter_vcf

result = prefilter_vcf(
    input_path="sample.vcf.gz",
    output_path="filtered.vcf.gz",
    preset="comprehensive",
)
print(f"Retained {result.stats.retained_variants}/{result.stats.input_variants} variants")
```

## Presets

| Preset | Gene loci | Regulatory | Safety net |
|--------|-----------|------------|------------|
| `comprehensive` | Full transcript (exon+intron+UTR) | ALL cCREs (PLS+pELS+dELS) | ClinVar P/LP |
| `coding-only` | Exon + UTR only | None | ClinVar P/LP |
| `regulatory-minimal` | Full transcript | PLS/pELS only | ClinVar P/LP |

- **comprehensive**: Maximum sensitivity. Retains all variants in genes, ncRNA, and regulatory elements. Best for discovery or when sensitivity is paramount.
- **coding-only**: Minimal region set. Only protein-coding exon and UTR variants, plus ClinVar P/LP safety net. Best for focused clinical pipelines.
- **regulatory-minimal**: Balances gene coverage with key regulatory elements (promoters and proximal enhancers). Skips distal enhancers.

## Updating reference data

Reference BED files are static. To update with new GENCODE/ENCODE/ClinVar releases:

```bash
python scripts/build_refs.py --output-dir refs/
```

Then re-deploy the skill or reinstall the package.

## Requirements

- Python >= 3.9
- bcftools >= 1.17
