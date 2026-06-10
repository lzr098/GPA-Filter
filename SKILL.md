---
name: dgra-prefilter
description: |
  全基因组 VCF 基因组区域预过滤模块（GPA Filter）v1.0.0。基于 GENCODE 基因座、ncRNA 区域、ENCODE 调控元件和 ClinVar 致病变异安全网，对 GRCh38 VCF 进行硬过滤。零外部 Python 运行时依赖，核心过滤委托 bcftools。三种预设（comprehensive/coding-only/regulatory-minimal），支持 CLI 和 Python API。

  **当以下情况时使用此 Skill**：
  (1) 用户提到"VCF 预过滤"、"基因组区域过滤"、"变异筛选"
  (2) 全基因组 VCF 文件过大需要精简后再做致病性分析
  (3) 需要仅保留基因区/ncRNA/调控元件/ClinVar 致病位点
  (4) 需要更新预过滤参考数据（GENCODE/ENCODE/ClinVar）
  (5) 任何涉及"prefilter"、"基因组预过滤"、"过滤VCF"的场景

  **需要系统安装 bcftools >= 1.17。**
---

# dgra-prefilter: Genomic Region Prefilter

## 概述

输入单样本 GRCh38 VCF/VCF.gz，基于预构建的保留区域 BED 文件（GENCODE 全转录本基因座、ncRNA 基因座、ENCODE 实验验证调控元件），以及 ClinVar 已知致病位点安全网，输出精简后的 VCF 和过滤统计报告。

## Trigger

当用户提到以下关键词时触发：
- VCF 预过滤 / 基因组区域过滤 / 变异筛选
- prefilter / 基因组预过滤 / 过滤 VCF
- 保留致病位点 / ClinVar 安全网
- 编码区过滤 / 调控元件过滤
- GPA Filter

## Parameters

| 参数 | 类型 | 必需 | 默认值 | 说明 |
|------|------|------|--------|------|
| input | string | 是 | — | 输入 VCF/VCF.gz/BCF 文件路径 |
| output | string | 否 | 自动生成 | 输出文件路径（.vcf.gz 自动压缩） |
| preset | string | 否 | comprehensive | 预设配置：comprehensive / coding-only / regulatory-minimal / regulatory-balanced |
| genome | string | 否 | GRCh38 | 基因组版本（仅支持 GRCh38） |
| ref_dir | string | 否 | ~/.dgra-prefilter/refs | 参考 BED 文件目录 |
| report | string | 否 | 与 output 同目录 | JSON 报告输出路径 |
| force | boolean | 否 | false | 跳过基因组版本校验 |
| update_refs | boolean | 否 | false | 触发参考数据更新 |
| annotate | boolean | 否 | false | 启用 DGRA_REGION/DGRA_SAFETYNET INFO 标注（较慢） |
| regulatory_source | string | 否 | fantom5 | 调控数据来源：fantom5 / ensembl / both |
| keep_all_chrM | boolean | 否 | false | 保留所有 chrM 变异 |

### Preset 说明

| Preset | 基因区 | ncRNA | 调控元件 | 安全网 |
|--------|--------|-------|----------|--------|
| comprehensive | 全转录本 | 全部 | ENCODE + FANTOM5 + Vista | ClinVar + OMIM |
| coding-only | 仅外显子+UTR | 无 | 无 | ClinVar + OMIM |
| regulatory-minimal | 全转录本 | 全部 | 仅 ENCODE PLS/pELS | ClinVar + OMIM |
| regulatory-balanced | 全转录本 | 全部 | ENCODE PLS/pELS/dELS/CTCF | ClinVar + OMIM |

## Output

### 文件输出
- 过滤后的 VCF/VCF.gz 文件（默认无 INFO 标注；加 `--annotate` 后含 DGRA_REGION / DGRA_SAFETYNET）
- JSON 统计报告

### 固定 Pipeline（两阶段）

**Phase 1 – 坐标过滤（始终运行）**
1. 输入校验（VCF 格式、基因组版本、bcftools）
2. 加载参考 BED（gene + ncRNA + cCRE）并合并
3. `bcftools view -T` 坐标硬过滤
4. ClinVar 安全网提取 + 合并去重

**Phase 2 – INFO 标注（可选，`--annotate`）**
5. 逐行添加 DGRA_REGION / DGRA_SAFETYNET 标签（Python bisect，较慢）

**Phase 3 – 报告（始终运行）**
6. 生成 JSON 统计报告

### Skill 返回格式
```
过滤完成！
- 输入变异：{input_variants:,}
- 保留变异：{retained_variants:,}（{retention_rate:.1%}）
- 安全网命中：ClinVar={clinvar_count}
- 耗时：{elapsed_seconds:.1f}s
- 输出文件：{output_path}
- 报告文件：{report_path}
```

## Installation

```bash
pip install dgra-prefilter
```

Requires: Python >= 3.9, bcftools >= 1.17

## Example

用户：「帮我过滤这个 VCF，只保留基因区和 ClinVar 致病变异」
→ dgra-prefilter --input /path/to/sample.vcf.gz --preset comprehensive

用户：「这个 VCF 只需要编码区，帮我预过滤」
→ dgra-prefilter --input /path/to/sample.vcf.gz --preset coding-only

用户：「用最小调控区域过滤这个 VCF」
→ dgra-prefilter --input /path/to/sample.vcf.gz --preset regulatory-minimal

用户：「用平衡调控区域过滤这个 VCF」
→ dgra-prefilter --input /path/to/sample.vcf.gz --preset regulatory-balanced

用户：「过滤后还要标注每个变异落在哪个区域」
→ dgra-prefilter --input /path/to/sample.vcf.gz --preset comprehensive --annotate

用户：「使用 Ensembl 调控数据过滤」
→ dgra-prefilter --input /path/to/sample.vcf.gz --preset comprehensive --regulatory-source ensembl

用户：「保留所有 chrM 变异」
→ dgra-prefilter --input /path/to/sample.vcf.gz --preset comprehensive --keep-all-chrM

用户：「帮我更新参考数据」
→ dgra-prefilter --update-refs

## Python API

```python
from dgra_prefilter import prefilter_vcf

# 默认：只过滤，不标注（最快）
result = prefilter_vcf(
    input_path="sample.vcf.gz",
    output_path="filtered.vcf.gz",
    preset="comprehensive",
)

# 启用标注（较慢）
result = prefilter_vcf(
    input_path="sample.vcf.gz",
    output_path="filtered.vcf.gz",
    preset="comprehensive",
    annotate=True,
)

## GitHub

https://github.com/lzr098/GPA-Filter
