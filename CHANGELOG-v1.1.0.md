# dgra-prefilter v1.1.0 变更摘要（Review 用）

> 20 个文件修改，3 个新测试文件，190 项测试全部通过（17.08s）

---

## 一、稳定性修复（Sprint 1 Group A — 7 项）

| # | 变更 | 文件 | 说明 |
|---|------|------|------|
| A1 | bcftools 版本检查 | `core.py`, `constants.py` | 解析版本号并与 `MIN_BCFTOOLS_VERSION = "1.17"` 比较，版本不足时抛出 `BcftoolsVersionError`（exit code 3002） |
| A2 | bcftools 超时 | `core.py`, `cli.py` | 新增 `--bcftools-timeout` 参数（默认 300s），`_run_bcftools()` 捕获 `TimeoutExpired` |
| A3 | header-only DGRA 检测 | `core.py` | `_has_dgra_annotations()` 仅扫描 `##` 开头行，遇到数据行立即停止，避免误匹配 |
| A4 | PID 前缀临时目录 | `core.py`, `build_refs.py` | `tempfile.TemporaryDirectory(prefix=f"dgra_prefilter_{os.getpid()}_")` 避免并发冲突 |
| A5 | 文件句柄安全 | `annotate.py` | 嵌套 `with` 语句包裹 `gzip.open()`，避免 `NameError` |
| A6 | 退出码统一 | `cli.py`, `constants.py` | 所有 `sys.exit()` 使用 `ErrorCode` 枚举值，不再硬编码 |
| A7 | `--strict-genome` | `cli.py`, `core.py` | 新增 flag，无法推断 genome 时报错而非警告 |

## 二、参考数据构建（Sprint 1 Group B — 4 项）

| # | 变更 | 文件 | 说明 |
|---|------|------|------|
| B1 | ClinVar Conflicting + 星级 | `build_refs.py` | 先 `bcftools norm -m -any` 拆分 multi-allelic；解析 `CLNSIGCONF`（Conflicting 但含 Pathogenic 也保留）；解析 `CLNREVSTAT` 映射为 1-4 star；仅保留 >=1-star；BED 第 4 列格式 `ClinVar_{N}star` |
| B2 | ClinVar URL + FANTOM5 HTTPS | `build_refs.py` | URL 改为无日期最新版 `clinvar.vcf.gz`；新增 `--clinvar-date` 参数；下载前 HTTP HEAD 检查；FANTOM5 改为 HTTPS |
| B3 | GENCODE UTR/CDS/splice BED | `build_refs.py`, `constants.py` | 新增 `gencode_5utr.bed`、`gencode_cds.bed`、`gencode_3utr.bed`、`gencode_splice_sites.bed`（外显子上下游各延伸 50bp） |
| B4 | ENCODE 类型标签 + balanced | `build_refs.py`, `constants.py`, `bed_utils.py` | ENCODE BED 增加第 4 列类型标签（PLS/pELS/dELS/CTCF）；新增 `encode_screen_v3_balanced.bed`（仅 PLS+pELS+dELS+CTCF）；`write_bed()` 支持可选第 4 列 |

## 三、生物学集成（Sprint 2 Group C — 4 项）

| # | 变更 | 文件 | 说明 |
|---|------|------|------|
| C1 | OMIM 安全网默认启用 | `presets.py`, `safetynet.py`, `cli.py` | `safetynet_omim=True`；解注释 `OMIMSafetyNet`；缺失时 warning 跳过不报错 |
| C2 | ClinVar star tag | `safetynet.py`, `annotate.py` | `DGRA_SAFETYNET` 输出星级标签如 `ClinVar_3star`、`ClinVar_1star`（0-star 已过滤） |
| C3 | UTR/剪接子区域注释 | `annotate.py`, `presets.py`, `report.py`, `constants.py` | `DGRA_REGION` 新增子标签：`gene_5utr`、`gene_cds`、`gene_3utr`、`gene_splice`；保留向后兼容汇总标签 `gene` |
| C4 | regulatory-balanced preset | `presets.py`, `cli.py`, `annotate.py`, `constants.py` | 新增 `--preset regulatory-balanced`：保留 PLS+pELS+dELS+CTCF；输出 `regulatory_pls/pels/dels/ctcf` + 汇总 `regulatory` |

## 四、可选数据源 + chrM（Sprint 3 — 3 项）

| # | 变更 | 文件 | 说明 |
|---|------|------|------|
| F1 | Ensembl Regulatory Build | `build_refs.py`, `constants.py`, `presets.py`, `cli.py`, `core.py`, `annotate.py` | 新增可选 `build_ensembl_regulatory_bed()` 从 GFF3 解析（promoter/enhancer/CTCF/open_chromatin/TF_binding）；新增 `--regulatory-source` 参数（`fantom5` 默认 / `ensembl` / `both`）；**不替换 FANTOM5 默认** |
| F2 | chrM 处理 | `presets.py`, `core.py`, `cli.py` | 新增 `--keep-all-chrM` flag（默认关闭）；开启时向 merged BED 追加 `chrM\t0\t16569`，保留全部 chrM 变异 |
| G1 | 全局一致性 | 所有文件 | 移除所有"产前"限定描述（无匹配）；ErrorCode 枚举统一使用；SKILL.md 更新 |

## 五、向后兼容保证

- 所有现有 CLI 参数默认行为不变
- 所有现有 preset（comprehensive / coding-only / regulatory-minimal）默认输出不变
- `DGRA_REGION` 汇总标签（`gene`、`regulatory`、`ncrna`）始终保留
- `DGRA_SAFETYNET` 格式兼容（旧格式 `ClinVar` → 新格式 `ClinVar_3star` 等）

## 六、测试覆盖

- 190 项测试全部通过（新增 15 项）
- 新增测试文件：`tests/test_benchmark.py`、`tests/test_build_ensembl.py`、`tests/test_chrM.py`

---

## Commit 建议

```bash
cd /Users/zhaorongli/WorkBuddy/2026-06-02-23-14-32/dgra-prefilter
git add -A
git commit -m "feat: v1.1.0 — stability, sub-regions, star tags, balanced preset, Ensembl optional, chrM

Sprint 1 (A1-A7, B1-B4):
- bcftools version check, timeout, strict-genome
- ClinVar Conflicting + star + multi-allelic norm
- GENCODE UTR/CDS/splice BED, ENCODE balanced BED

Sprint 2 (C1-C4):
- OMIM safety net default enabled
- ClinVar star tag annotation
- UTR 5/3 + splice sub-region tags (backward compatible)
- regulatory-balanced preset

Sprint 3 (F1-F2, G1):
- Ensembl Regulatory Build optional support (fantom5/ensembl/both)
- --keep-all-chrM flag
- Global consistency cleanup"
```
