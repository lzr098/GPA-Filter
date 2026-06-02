# dgra-prefilter

Genomic region prefilter for whole-genome VCF files based on GENCODE, ENCODE, and ClinVar.

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

| Preset | Gene | ncRNA | Regulatory | Safety Net |
|--------|------|-------|------------|------------|
| comprehensive | Full transcript | All | ENCODE + FANTOM5 + Vista | ClinVar |
| coding-only | Exon+UTR only | None | None | ClinVar |
| regulatory-minimal | Full transcript | All | ENCODE PLS/pELS only | ClinVar |

## Requirements

- Python >= 3.9
- bcftools >= 1.17
