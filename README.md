<h1 align="center">
  <code>dgra-prefilter</code>
</h1>

<p align="center">
  <a href="https://github.com/lzr098/GPA-Filter"><img src="https://img.shields.io/badge/version-1.1.0-blue" alt="version"></a>
  <a href="#"><img src="https://img.shields.io/badge/requires-bcftools-orange" alt="bcftools"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-lightgrey" alt="license"></a>
</p>

<p align="center">
  <strong>Genomic region prefilter for whole-genome VCF files</strong><br>
  GENCODE · ENCODE · ClinVar · Zero runtime Python dependencies
</p>

<p align="center">
  <a href="#中文文档">中文</a> · <a href="#english-documentation">English</a>
</p>

---

<h2 id="中文文档">📖 中文文档</h2>

### 一句话

**dgra-prefilter**（GPA Filter）是一个全基因组 VCF 预过滤工具。基于预构建的保留区域 BED 文件（GENCODE 基因座、ncRNA 区域、ENCODE 调控元件）和 ClinVar 致病变异安全网，通过 `bcftools` 坐标硬过滤，将 WGS VCF 快速精简到生物学上"值得关注"的变异子集。

> ⚡ **核心特点**：零外部 Python 运行时依赖，过滤速度完全取决于 bcftools 的 BED 交集性能。百万级变异通常可在数秒内完成。

---

### 📋 目录

- [功能特性](#功能特性)
- [系统要求](#系统要求)
- [部署方式](#部署方式)
  - [WorkBuddy 部署](#workbuddy-部署)
  - [CodeX / OpenClaw 部署](#codex--openclaw-部署)
  - [通用部署（pip install / git clone）](#通用部署)
- [依赖安装](#依赖安装)
- [快速开始](#快速开始)
- [CLI 使用指南](#cli-使用指南)
- [Python API 使用](#python-api-使用)
- [Preset 说明](#preset-说明)
- [参考数据管理](#参考数据管理)
- [输出说明](#输出说明)
- [常见问题](#常见问题)
- [测试](#测试)

---

### 功能特性

| 特性 | 说明 |
|------|------|
| 🧬 **三层数据过滤** | GENCODE 基因座 + ncRNA + ENCODE 调控元件 |
| 🛡️ **ClinVar 安全网** | 保留所有 ClinVar 致病/可能致病变异（P/LP），即使落在基因区外 |
| ⚡ **bcftools 硬过滤** | 基于坐标交集，无实时 API 调用，速度极快 |
| 📦 **五种预设策略** | comprehensive / comprehensive-splice100 / coding-only / regulatory-minimal / regulatory-balanced |
| 🎛️ **交互式预设选择** | 运行时通过菜单选择过滤策略 (`-I` / `--interactive`) |
| 🔧 **macOS 引号兼容** | 自动处理拖入文件路径中的 Unicode/ASCII 引号差异 |
| 🏷️ **可选 INFO 标注** | DGRA_REGION / DGRA_SAFETYNET VCF 标签 |
| 🔄 **参考数据可更新** | 支持手动更新 GENCODE / ENCODE / ClinVar 数据 |
| 🎯 **零 Python 依赖** | 仅依赖系统 bcftools，无外部 Python 包 |

---

### 系统要求

| 项目 | 要求 |
|------|------|
| Python | ≥ 3.9 |
| bcftools | ≥ 1.17（系统级依赖，**必须预先安装**） |
| 操作系统 | macOS / Linux / Windows (WSL) |
| 基因组版本 | 仅支持 GRCh38 |

---

### 部署方式

#### WorkBuddy 部署

WorkBuddy 是 macOS 上的 AI 助手桌面应用，支持 Skill 扩展。

**步骤 1：打开 Skill 目录**

```bash
open ~/.workbuddy/skills/
```

**步骤 2：克隆仓库到 Skill 目录**

```bash
cd ~/.workbuddy/skills/
git clone https://github.com/lzr098/GPA-Filter.git dgra-prefilter
```

**步骤 3：安装依赖并构建参考数据**

```bash
cd dgra-prefilter
pip install -e .                # 安装 Python 包
python scripts/build_refs.py    # 构建参考 BED 文件（首次必需）
```

**步骤 4：重启 WorkBuddy**

- WorkBuddy 自动扫描 `~/.workbuddy/skills/` 目录
- 重启应用后，dgra-prefilter Skill 即可使用

> 💡 **提示**：确保系统已安装 bcftools ≥ 1.17。WorkBuddy 也可能在 `~/.workbuddy/binaries/bcftools/` 管理 bcftools 安装。

#### CodeX / OpenClaw 部署

CodeX（或 OpenClaw）是命令行/IDE 集成的 AI 编程助手。

**步骤 1：找到 Skill 目录**

```bash
mkdir -p ~/.codex/skills   # 或 ~/.openclaw/skills
cd ~/.codex/skills/
```

**步骤 2：克隆仓库**

```bash
git clone https://github.com/lzr098/GPA-Filter.git dgra-prefilter
```

**步骤 3：安装**

```bash
cd dgra-prefilter
pip install -e .
python scripts/build_refs.py
```

**步骤 4：重启 CodeX**

重启 IDE 或 CodeX 扩展，dgra-prefilter 即可在对话中使用。

#### 通用部署

**方式 A：pip 安装（推荐最终用户）**

```bash
pip install dgra-prefilter

# 构建参考数据
python -m dgra_prefilter.build_refs  # 或通过源码运行 scripts/build_refs.py
```

**方式 B：源码安装（推荐开发者）**

```bash
git clone https://github.com/lzr098/GPA-Filter.git
cd GPA-Filter
pip install -e ".[dev]"           # 包含测试依赖
python scripts/build_refs.py       # 构建参考 BED

# 验证安装
dgra-prefilter --version
```

---

### 依赖安装

#### bcftools（系统级必需）

bcftools 是 dgra-prefilter 的唯一外部依赖，**必须预先安装**。

**macOS（Homebrew）**：

```bash
brew install bcftools
bcftools --version    # 确认 >= 1.17
```

**Linux（Ubuntu/Debian）**：

```bash
sudo apt-get update
sudo apt-get install bcftools
bcftools --version
```

**Linux（CentOS/RHEL）**：

```bash
sudo yum install bcftools
# 或从源码编译
```

**Conda**：

```bash
conda install -c bioconda bcftools
```

#### Python 包

```bash
# 生产环境（零额外 Python 依赖）
pip install dgra-prefilter

# 开发环境
pip install -e ".[dev]"

# 构建参考数据（需要 requests）
pip install -e ".[build]"
```

---

### 快速开始

```bash
# 基础过滤：comprehensive preset（最大敏感度）
dgra-prefilter \
  --input sample.vcf.gz \
  --output filtered.vcf.gz \
  --preset comprehensive

# 交互式选择预设（推荐）
dgra-prefilter \
  --input sample.vcf.gz \
  --output filtered.vcf.gz \
  --interactive

# Exon/UTR + 100bp 剪接区（去除深度内含子）
dgra-prefilter \
  --input sample.vcf.gz \
  --output filtered.vcf.gz \
  --preset comprehensive-splice100

# 仅保留编码区（最精简）
dgra-prefilter \
  --input sample.vcf.gz \
  --output coding_only.vcf.gz \
  --preset coding-only

# 过滤 + 区域标注（较慢，约 1.5~2 倍时间）
dgra-prefilter \
  --input sample.vcf.gz \
  --output annotated.vcf.gz \
  --preset comprehensive \
  --annotate

# 过滤前更新参考数据
dgra-prefilter \
  --input sample.vcf.gz \
  --output filtered.vcf.gz \
  --update-refs
```

---

### CLI 使用指南

#### 必选参数

| 参数 | 说明 | 示例 |
|------|------|------|
| `-i, --input PATH` | 输入 VCF/VCF.gz/BCF | `sample.vcf.gz` |
| `-o, --output PATH` | 输出 VCF 路径 | `filtered.vcf.gz` |

#### 可选参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-p, --preset` | 预设策略：`comprehensive` / `comprehensive-splice100` / `coding-only` / `regulatory-minimal` / `regulatory-balanced` | `comprehensive` |
| `-I, --interactive` | 交互式选择预设（运行前弹出菜单） | 否 |
| `-g, --genome` | 基因组版本（仅 GRCh38） | `GRCh38` |
| `--ref-dir PATH` | 参考 BED 文件目录 | `~/.dgra-prefilter/refs` |
| `--report PATH` | JSON 报告输出路径 | 与 output 同目录 |
| `--annotate` | 启用 DGRA_REGION / DGRA_SAFETYNET 标注 | 否 |
| `--update-refs` | 过滤前更新参考数据 | 否 |
| `--regulatory-source` | 调控数据来源：`fantom5` / `ensembl` / `both` | `fantom5` |
| `--keep-all-chrM` | 保留所有 chrM 变异（无视区域） | 否 |
| `--force` | 跳过基因组版本校验 | 否 |
| `-v, --verbose` | 启用 DEBUG 级别日志 | 否 |
| `--version` | 显示版本号 | — |

#### 完整示例

```bash
dgra-prefilter \
  --input /data/wgs_sample.vcf.gz \
  --output /data/filtered.vcf.gz \
  --preset regulatory-minimal \
  --annotate \
  --report /data/filter_report.json \
  --verbose
```

#### 退出码

| 退出码 | 含义 |
|--------|------|
| 0 | 成功 |
| 1 | 文件未找到 |
| 2 | 基因组版本不匹配 |
| 3 | bcftools 未安装 |
| 4 | 参考数据缺失 |
| 5 | VCF 处理错误 |
| 99 | 未知错误 |

---

### Python API 使用

```python
from dgra_prefilter import prefilter_vcf

# 基础过滤
result = prefilter_vcf(
    input_path="sample.vcf.gz",
    output_path="filtered.vcf.gz",
    preset="comprehensive",
)

print(f"保留: {result.stats.retained_variants:,} / {result.stats.input_variants:,}")
print(f"保留率: {result.stats.retention_rate:.1%}")
print(f"输出: {result.output_path}")
print(f"报告: {result.report_path}")

# 交互式选择预设
result = prefilter_vcf(
    input_path="sample.vcf.gz",
    output_path="filtered.vcf.gz",
    interactive=True,
)

# 过滤 + 标注
result = prefilter_vcf(
    input_path="sample.vcf.gz",
    output_path="annotated.vcf.gz",
    preset="comprehensive",
    annotate=True,
)

# 使用 splice100 预设（去除深度内含子）
result = prefilter_vcf(
    input_path="sample.vcf.gz",
    output_path="filtered.vcf.gz",
    preset="comprehensive-splice100",
)

# 使用 PrefilterConfig（类型安全）
from dgra_prefilter import PrefilterConfig

config = PrefilterConfig(
    input_path="sample.vcf.gz",
    output_path="filtered.vcf.gz",
    preset_name="coding-only",
    annotate=False,
)
# 注意：PrefilterConfig 为 dataclass，需配合 core 内部逻辑使用
```

#### 异常处理

```python
from dgra_prefilter import (
    prefilter_vcf,
    BcftoolsNotFoundError,
    GenomeMismatchError,
    RefDataMissingError,
    VCFProcessingError,
)

try:
    result = prefilter_vcf(input_path="sample.vcf.gz", output_path="out.vcf.gz")
except BcftoolsNotFoundError:
    print("错误：未找到 bcftools，请先安装")
except GenomeMismatchError as e:
    print(f"错误：基因组版本不匹配 - {e}")
except RefDataMissingError:
    print("错误：参考数据缺失，请运行 build_refs.py")
except VCFProcessingError as e:
    print(f"错误：VCF 处理失败 - {e}")
```

---

### Preset 说明

| Preset | 基因区 | ncRNA | 调控元件 | ClinVar 安全网 | 保留率* |
|--------|--------|-------|----------|----------------|--------|
| **comprehensive** | 全转录本（外显子+内含子+UTR） | 全部 | ENCODE cCRE + FANTOM5 + Vista | P/LP | ~50-60% |
| **comprehensive-splice100** | Exon/UTR + 100bp 剪接区 | 全部 | ENCODE cCRE + FANTOM5 + Vista | P/LP | ~25-30% |
| **coding-only** | 仅外显子+UTR | 无 | 无 | P/LP | ~1-3% |
| **regulatory-minimal** | 全转录本 | 全部 | 仅 ENCODE PLS/pELS | P/LP | ~10-20% |
| **regulatory-balanced** | 全转录本 | 全部 | ENCODE 平衡子集 | P/LP | ~15-25% |

- **comprehensive**：最大敏感度，适合发现研究或敏感度优先的场景
- **comprehensive-splice100**：去除深度内含子，保留剪接位点附近区域，平衡敏感度与精简度
- **coding-only**：最小区域集，仅蛋白编码外显子和 UTR，适合聚焦临床流程
- **regulatory-minimal**：平衡基因覆盖与关键调控元件，跳过远端增强子
- **regulatory-balanced**：更均衡的调控元件覆盖策略

> *保留率基于全基因组 VCF（~500 万变异）估算，实际值因样本而异。

---

### 参考数据管理

#### 首次构建

```bash
# 方式 1：通过 pip 安装的模块
python -m dgra_prefilter.build_refs

# 方式 2：通过源码
python scripts/build_refs.py --output-dir ~/.dgra-prefilter/refs
```

build_refs.py 会下载并解析：
- GENCODE v44 GTF → 基因座 BED
- ENCODE SCREEN v3 → cCRE BED
- ClinVar VCF → 致病位点 BED

#### 更新参考数据

```bash
# 过滤时自动更新
dgra-prefilter --input sample.vcf.gz --output out.vcf.gz --update-refs

# 或单独更新
python scripts/build_refs.py --output-dir ~/.dgra-prefilter/refs
```

#### 参考数据目录结构

```
~/.dgra-prefilter/refs/
├── gencode_gene.bed
├── gencode_ncrna.bed
├── gencode_coding_exon_utr.bed
├── encode_ccre.bed
├── encode_pls_pels.bed
├── fantom5.bed
├── vista.bed
├── clinvar_pathogenic.bed
├── omim.bed
└── manifest.json
```

---

### 输出说明

#### 过滤后的 VCF

- 保留落在保留区域内的变异 + ClinVar 安全网命中变异
- 默认无新增 INFO 标签（最快）
- `--annotate` 时添加：
  - `DGRA_REGION`：区域类型（`gene` / `ncrna` / `regulatory`）
  - `DGRA_SAFETYNET`：安全网来源（`ClinVar` / `OMIM`）

#### JSON 报告示例

```json
{
  "input_variants": 3204567,
  "retained_variants": 98765,
  "retention_rate": 0.0308,
  "region_only_variants": 95000,
  "safetynet_only_variants": 1200,
  "region_and_safetynet_variants": 2565,
  "clinvar_count": 3765,
  "omim_count": 0,
  "elapsed_seconds": 4.2,
  "preset": "comprehensive",
  "ref_data_versions": {
    "gencode": "v44",
    "encode": "v3",
    "clinvar": "20240603"
  }
}
```

---

### 常见问题

| 问题 | 原因 | 解决 |
|------|------|------|
| `BcftoolsNotFoundError` | bcftools 未安装或不在 PATH | 安装 bcftools ≥ 1.17 |
| `RefDataMissingError` | 参考 BED 文件缺失 | 运行 `python scripts/build_refs.py` |
| `GenomeMismatchError` | VCF 基因组版本非 GRCh38 | 使用 `--force` 跳过校验，或确认输入为 GRCh38 |
| `VCFProcessingError` | bcftools 命令执行失败 | 检查 VCF 格式是否合法，查看 stderr |
| 保留变异为 0 | BED 文件为空或路径错误 | 检查 `--ref-dir` 指向的目录 |
| `--annotate` 非常慢 | 标注需逐变异查 BED 交集 | 这是预期行为，约 1.5~2 倍过滤时间 |

---

### 测试

```bash
# 安装开发依赖
pip install -e ".[dev]"

# 运行测试
pytest

# 指定测试文件
pytest tests/test_core.py
pytest tests/test_cli.py

# 详细输出
pytest -v
```

---

<hr>

<h2 id="english-documentation">📖 English Documentation</h2>

### One-Liner

**dgra-prefilter** (GPA Filter) is a whole-genome VCF prefiltering tool. Based on pre-built retention region BED files (GENCODE gene loci, ncRNA regions, ENCODE regulatory elements) and a ClinVar pathogenic variant safety net, it rapidly filters WGS VCF files to a biologically "interesting" variant subset through coordinate-based hard filtering via `bcftools`.

> ⚡ **Core feature**: Zero external Python runtime dependencies. Filtering speed depends entirely on bcftools BED intersection performance. Millions of variants typically complete in seconds.

---

### Table of Contents

- [Features](#features)
- [System Requirements](#system-requirements)
- [Deployment](#deployment)
  - [WorkBuddy](#workbuddy-deployment)
  - [CodeX / OpenClaw](#codex--openclaw-deployment)
  - [Generic Deployment](#generic-deployment)
- [Dependency Installation](#dependency-installation)
- [Quick Start](#quick-start)
- [CLI Usage Guide](#cli-usage-guide)
- [Python API Usage](#python-api-usage)
- [Presets](#presets)
- [Reference Data Management](#reference-data-management)
- [Output](#output)
- [FAQ](#faq)
- [Testing](#testing)

---

### Features

| Feature | Description |
|---------|-------------|
| 🧬 **Three-layer filtering** | GENCODE loci + ncRNA + ENCODE regulatory elements |
| 🛡️ **ClinVar safety net** | Retains all ClinVar pathogenic/likely pathogenic (P/LP) variants, even outside gene loci |
| ⚡ **bcftools hard filter** | Coordinate-based intersection, no live API calls, extremely fast |
| 📦 **Five preset strategies** | comprehensive / comprehensive-splice100 / coding-only / regulatory-minimal / regulatory-balanced |
| 🎛️ **Interactive preset selection** | Choose filter strategy via runtime menu (`-I` / `--interactive`) |
| 🔧 **macOS quote compatibility** | Auto-handles Unicode/ASCII quotation mark differences in dragged file paths |
| 🏷️ **Optional INFO annotation** | DGRA_REGION / DGRA_SAFETYNET VCF tags |
| 🔄 **Updatable references** | Manual update of GENCODE / ENCODE / ClinVar data supported |
| 🎯 **Zero Python deps** | Only requires system bcftools, no external Python packages |

---

### System Requirements

| Item | Requirement |
|------|-------------|
| Python | ≥ 3.9 |
| bcftools | ≥ 1.17 (system dependency, **must be pre-installed**) |
| OS | macOS / Linux / Windows (WSL) |
| Genome | GRCh38 only |

---

### Deployment

#### WorkBuddy Deployment

WorkBuddy is an AI assistant desktop app for macOS that supports Skill extensions.

**Step 1: Open the Skill directory**

```bash
open ~/.workbuddy/skills/
```

**Step 2: Clone into the Skill directory**

```bash
cd ~/.workbuddy/skills/
git clone https://github.com/lzr098/GPA-Filter.git dgra-prefilter
```

**Step 3: Install dependencies and build reference data**

```bash
cd dgra-prefilter
pip install -e .                # Install Python package
python scripts/build_refs.py    # Build reference BED files (required first time)
```

**Step 4: Restart WorkBuddy**

- WorkBuddy auto-scans `~/.workbuddy/skills/`
- Restart the app, then dgra-prefilter Skill is ready

> 💡 **Tip**: Ensure bcftools ≥ 1.17 is installed. WorkBuddy may also manage bcftools under `~/.workbuddy/binaries/bcftools/`.

#### CodeX / OpenClaw Deployment

CodeX (or OpenClaw) is a command-line/IDE-integrated AI programming assistant.

**Step 1: Locate the Skill directory**

```bash
mkdir -p ~/.codex/skills   # or ~/.openclaw/skills
cd ~/.codex/skills/
```

**Step 2: Clone the repository**

```bash
git clone https://github.com/lzr098/GPA-Filter.git dgra-prefilter
```

**Step 3: Install**

```bash
cd dgra-prefilter
pip install -e .
python scripts/build_refs.py
```

**Step 4: Restart CodeX**

Restart the IDE or CodeX extension, dgra-prefilter is available in chat.

#### Generic Deployment

**Option A: pip install (recommended for end users)**

```bash
pip install dgra-prefilter

# Build reference data
python -m dgra_prefilter.build_refs  # or run scripts/build_refs.py from source
```

**Option B: Source install (recommended for developers)**

```bash
git clone https://github.com/lzr098/GPA-Filter.git
cd GPA-Filter
pip install -e ".[dev]"           # Includes test dependencies
python scripts/build_refs.py       # Build reference BEDs

# Verify installation
dgra-prefilter --version
```

---

### Dependency Installation

#### bcftools (System-level Required)

bcftools is the only external dependency of dgra-prefilter and **must be pre-installed**.

**macOS (Homebrew)**:

```bash
brew install bcftools
bcftools --version    # Verify >= 1.17
```

**Linux (Ubuntu/Debian)**:

```bash
sudo apt-get update
sudo apt-get install bcftools
bcftools --version
```

**Linux (CentOS/RHEL)**:

```bash
sudo yum install bcftools
# Or compile from source
```

**Conda**:

```bash
conda install -c bioconda bcftools
```

#### Python Packages

```bash
# Production (zero extra Python dependencies)
pip install dgra-prefilter

# Development
pip install -e ".[dev]"

# Building reference data (requires requests)
pip install -e ".[build]"
```

---

### Quick Start

```bash
# Basic filtering: comprehensive preset (maximum sensitivity)
dgra-prefilter \
  --input sample.vcf.gz \
  --output filtered.vcf.gz \
  --preset comprehensive

# Interactive preset selection (recommended)
dgra-prefilter \
  --input sample.vcf.gz \
  --output filtered.vcf.gz \
  --interactive

# Exon/UTR + 100bp splice window (removes deep intronic)
dgra-prefilter \
  --input sample.vcf.gz \
  --output filtered.vcf.gz \
  --preset comprehensive-splice100

# Keep coding regions only (most compact)
dgra-prefilter \
  --input sample.vcf.gz \
  --output coding_only.vcf.gz \
  --preset coding-only

# Filter + annotate (slower, ~1.5-2x time)
dgra-prefilter \
  --input sample.vcf.gz \
  --output annotated.vcf.gz \
  --preset comprehensive \
  --annotate

# Update references before filtering
dgra-prefilter \
  --input sample.vcf.gz \
  --output filtered.vcf.gz \
  --update-refs
```

---

### CLI Usage Guide

#### Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `-i, --input PATH` | Input VCF/VCF.gz/BCF | `sample.vcf.gz` |
| `-o, --output PATH` | Output VCF path | `filtered.vcf.gz` |

#### Optional Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `-p, --preset` | Preset: `comprehensive` / `comprehensive-splice100` / `coding-only` / `regulatory-minimal` / `regulatory-balanced` | `comprehensive` |
| `-I, --interactive` | Interactive preset selection (menu before filtering) | No |
| `-g, --genome` | Genome version (GRCh38 only) | `GRCh38` |
| `--ref-dir PATH` | Reference BED files directory | `~/.dgra-prefilter/refs` |
| `--report PATH` | JSON report output path | Same dir as output |
| `--annotate` | Enable DGRA_REGION / DGRA_SAFETYNET tags | No |
| `--update-refs` | Update reference data before filtering | No |
| `--regulatory-source` | Regulatory source: `fantom5` / `ensembl` / `both` | `fantom5` |
| `--keep-all-chrM` | Retain all chrM variants regardless of region | No |
| `--force` | Skip genome version validation | No |
| `-v, --verbose` | Enable DEBUG logging | No |
| `--version` | Show version | — |

#### Complete Example

```bash
dgra-prefilter \
  --input /data/wgs_sample.vcf.gz \
  --output /data/filtered.vcf.gz \
  --preset regulatory-minimal \
  --annotate \
  --report /data/filter_report.json \
  --verbose
```

#### Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | File not found |
| 2 | Genome version mismatch |
| 3 | bcftools not installed |
| 4 | Reference data missing |
| 5 | VCF processing error |
| 99 | Unexpected error |

---

### Python API Usage

```python
from dgra_prefilter import prefilter_vcf

# Basic filtering
result = prefilter_vcf(
    input_path="sample.vcf.gz",
    output_path="filtered.vcf.gz",
    preset="comprehensive",
)

print(f"Retained: {result.stats.retained_variants:,} / {result.stats.input_variants:,}")
print(f"Retention: {result.stats.retention_rate:.1%}")
print(f"Output: {result.output_path}")
print(f"Report: {result.report_path}")

# Interactive preset selection
result = prefilter_vcf(
    input_path="sample.vcf.gz",
    output_path="filtered.vcf.gz",
    interactive=True,
)

# Filter + annotate
result = prefilter_vcf(
    input_path="sample.vcf.gz",
    output_path="annotated.vcf.gz",
    preset="comprehensive",
    annotate=True,
)

# Use splice100 preset (removes deep intronic)
result = prefilter_vcf(
    input_path="sample.vcf.gz",
    output_path="filtered.vcf.gz",
    preset="comprehensive-splice100",
)

# Using PrefilterConfig (type-safe)
from dgra_prefilter import PrefilterConfig

config = PrefilterConfig(
    input_path="sample.vcf.gz",
    output_path="filtered.vcf.gz",
    preset_name="coding-only",
    annotate=False,
)
# Note: PrefilterConfig is a dataclass, use with core internal logic
```

#### Exception Handling

```python
from dgra_prefilter import (
    prefilter_vcf,
    BcftoolsNotFoundError,
    GenomeMismatchError,
    RefDataMissingError,
    VCFProcessingError,
)

try:
    result = prefilter_vcf(input_path="sample.vcf.gz", output_path="out.vcf.gz")
except BcftoolsNotFoundError:
    print("Error: bcftools not found, please install first")
except GenomeMismatchError as e:
    print(f"Error: Genome version mismatch - {e}")
except RefDataMissingError:
    print("Error: Reference data missing, run build_refs.py")
except VCFProcessingError as e:
    print(f"Error: VCF processing failed - {e}")
```

---

### Presets

| Preset | Gene Loci | ncRNA | Regulatory | ClinVar Safety Net | Retention* |
|--------|-----------|-------|------------|---------------------|------------|
| **comprehensive** | Full transcript (exon+intron+UTR) | All | ENCODE cCRE + FANTOM5 + Vista | P/LP | ~50-60% |
| **comprehensive-splice100** | Exon/UTR + 100bp splice window | All | ENCODE cCRE + FANTOM5 + Vista | P/LP | ~25-30% |
| **coding-only** | Exon + UTR only | None | None | P/LP | ~1-3% |
| **regulatory-minimal** | Full transcript | All | ENCODE PLS/pELS only | P/LP | ~10-20% |
| **regulatory-balanced** | Full transcript | All | ENCODE balanced subset | P/LP | ~15-25% |

- **comprehensive**: Maximum sensitivity, best for discovery or sensitivity-first scenarios
- **comprehensive-splice100**: Removes deep intronic regions while keeping splice junctions; balances sensitivity and compactness
- **coding-only**: Minimal region set, only protein-coding exons and UTRs, best for focused clinical pipelines
- **regulatory-minimal**: Balances gene coverage with key regulatory elements, skips distal enhancers
- **regulatory-balanced**: More balanced regulatory element coverage strategy

> *Retention rates are estimates based on whole-genome VCF (~5M variants); actual values vary by sample.

---

### Reference Data Management

#### First-time Build

```bash
# Via pip-installed module
python -m dgra_prefilter.build_refs

# Via source
python scripts/build_refs.py --output-dir ~/.dgra-prefilter/refs
```

build_refs.py downloads and parses:
- GENCODE v44 GTF → gene loci BED
- ENCODE SCREEN v3 → cCRE BED
- ClinVar VCF → pathogenic site BED

#### Update References

```bash
# Auto-update during filtering
dgra-prefilter --input sample.vcf.gz --output out.vcf.gz --update-refs

# Or standalone update
python scripts/build_refs.py --output-dir ~/.dgra-prefilter/refs
```

#### Reference Data Directory

```
~/.dgra-prefilter/refs/
├── gencode_gene.bed
├── gencode_ncrna.bed
├── gencode_coding_exon_utr.bed
├── encode_ccre.bed
├── encode_pls_pels.bed
├── fantom5.bed
├── vista.bed
├── clinvar_pathogenic.bed
├── omim.bed
└── manifest.json
```

---

### Output

#### Filtered VCF

- Retains variants in retention regions + ClinVar safety net hits
- Default: no new INFO tags (fastest)
- With `--annotate`, adds:
  - `DGRA_REGION`: region type (`gene` / `ncrna` / `regulatory`)
  - `DGRA_SAFETYNET`: safety net source (`ClinVar` / `OMIM`)

#### JSON Report Example

```json
{
  "input_variants": 3204567,
  "retained_variants": 98765,
  "retention_rate": 0.0308,
  "region_only_variants": 95000,
  "safetynet_only_variants": 1200,
  "region_and_safetynet_variants": 2565,
  "clinvar_count": 3765,
  "omim_count": 0,
  "elapsed_seconds": 4.2,
  "preset": "comprehensive",
  "ref_data_versions": {
    "gencode": "v44",
    "encode": "v3",
    "clinvar": "20240603"
  }
}
```

---

### FAQ

| Issue | Cause | Solution |
|-------|-------|----------|
| `BcftoolsNotFoundError` | bcftools not installed or not in PATH | Install bcftools ≥ 1.17 |
| `RefDataMissingError` | Reference BED files missing | Run `python scripts/build_refs.py` |
| `GenomeMismatchError` | VCF genome not GRCh38 | Use `--force` to skip, or confirm GRCh38 input |
| `VCFProcessingError` | bcftools command failed | Check VCF validity, inspect stderr |
| 0 retained variants | BED files empty or wrong path | Check `--ref-dir` directory |
| `--annotate` very slow | Annotation requires per-variant BED lookup | Expected, ~1.5-2x filter time |

---

### Testing

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Specific test file
pytest tests/test_core.py
pytest tests/test_cli.py

# Verbose output
pytest -v
```

---

## Architecture

```
Input VCF
    │
    ▼
┌─────────────────────┐
│  Input Validation   │  VCF format, genome version (GRCh38), bcftools check
└─────────────────────┘
    │
    ▼
┌─────────────────────┐
│  Load Reference BED │  Gene + ncRNA + Regulatory (per preset)
│  + ClinVar Safety   │  Pathogenic/Likely pathogenic variants
└─────────────────────┘
    │
    ▼
┌─────────────────────┐
│  bcftools view -T   │  Coordinate intersection (hard filter)
│  + Safety Net Merge │  Union of region hits + ClinVar hits, deduplicated
└─────────────────────┘
    │
    ├─ annotate ─► DGRA_REGION / DGRA_SAFETYNET INFO tags (optional)
    │
    ▼
┌─────────────────────┐
│  JSON Report        │  Statistics, retention rate, data versions
└─────────────────────┘
    │
    ▼
Output VCF + Report
```

---

## Project Structure

```
dgra-prefilter/
├── src/dgra_prefilter/           # Python package
│   ├── __init__.py               # Public API exports
│   ├── cli.py                    # CLI entry point
│   ├── core.py                   # Core filtering engine
│   ├── annotate.py               # VCF INFO annotation
│   ├── bed_utils.py              # BED file utilities
│   ├── presets.py                # Preset configurations
│   ├── ref_manager.py            # Reference data management
│   ├── report.py                 # JSON report generation
│   ├── safetynet.py              # ClinVar/OMIM safety net
│   └── constants.py              # Default paths and constants
├── scripts/
│   └── build_refs.py             # Reference BED builder
├── tests/                        # Test suite
│   ├── test_core.py
│   ├── test_cli.py
│   ├── test_annotate.py
│   └── ...
├── refs/                         # Reference BED files (generated)
├── docs/
│   └── genomic-region-coverage.svg
├── pyproject.toml                # Package configuration
├── README.md                     # This file
└── .gitignore
```

---

## Data Sources

| Layer | Source | Version |
|-------|--------|---------|
| Gene loci | GENCODE | v44 |
| Regulatory elements | ENCODE SCREEN | v3 |
| ncRNA | GENCODE | v44 |
| Safety net | ClinVar | Latest |
| Additional | OMIM, FANTOM5, Vista | — |

---

## License

MIT

---

## Related Skills · 相关技能

| Skill · 技能 | Repo · 仓库 | Purpose · 用途 |
|---|---|---|
| **GPA** | [lzr098/dgra-genomic-risk](https://github.com/lzr098/dgra-genomic-risk) | Whole-genome phenotype association |
| **variant-impact** | [lzr098/variant-impact](https://github.com/lzr098/variant-impact) | Single variant ACMG classification |
| **disease-risk-query** | [lzr098/Disease-Risk-Query](https://github.com/lzr098/Disease-Risk-Query) | Disease-specific genetic risk |

---

**Maintainer**: [@lzr098](https://github.com/lzr098)  
**Current Version**: 1.1.0  
**Last Updated**: 2026-06-10
