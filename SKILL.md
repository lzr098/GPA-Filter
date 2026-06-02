# dgra-prefilter

## Description
全基因组 VCF 基因组区域预过滤模块。基于 GENCODE 基因座、ncRNA 区域、ENCODE 调控元件和 ClinVar 致病变异安全网，对 GRCh38 VCF 进行硬过滤。过滤后仅保留具有生物学功能意义的区域变异和已知致病位点，大幅减少下游分析的计算负担。

## Trigger
当用户提到以下关键词时触发：
- VCF 预过滤
- 基因组区域过滤
- 变异筛选
- prefilter
- 基因组预过滤
- 过滤 VCF
- 保留致病位点
- ClinVar 安全网
- 编码区过滤
- 调控元件过滤

## Parameters

| 参数 | 类型 | 必需 | 默认值 | 说明 |
|------|------|------|--------|------|
| input | string | 是 | — | 输入 VCF/VCF.gz/BCF 文件路径 |
| output | string | 否 | 自动生成 | 输出文件路径（.vcf.gz 自动压缩） |
| preset | string | 否 | comprehensive | 预设配置：comprehensive / coding-only / regulatory-minimal |
| genome | string | 否 | GRCh38 | 基因组版本（仅支持 GRCh38） |
| ref_dir | string | 否 | ~/.dgra-prefilter/refs | 参考 BED 文件目录 |
| report | string | 否 | 与 output 同目录 | JSON 报告输出路径 |
| force | boolean | 否 | false | 跳过基因组版本校验 |
| update_refs | boolean | 否 | false | 触发参考数据更新 |

### Preset 说明

| Preset | 基因区 | ncRNA | 调控元件 | 安全网 |
|--------|--------|-------|----------|--------|
| comprehensive | 全转录本 | 全部 | ENCODE + FANTOM5 + Vista | ClinVar |
| coding-only | 仅外显子+UTR | 无 | 无 | ClinVar |
| regulatory-minimal | 全转录本 | 全部 | 仅 ENCODE PLS/pELS | ClinVar |

## Output

### 文件输出
- 过滤后的 VCF/VCF.gz 文件（含 DGRA_REGION 和 DGRA_SAFETYNET INFO 标注）
- JSON 统计报告

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

## Example

用户：「帮我过滤这个 VCF，只保留基因区和 ClinVar 致病变异」
→ dgra-prefilter --input /path/to/sample.vcf.gz --preset comprehensive

用户：「这个 VCF 只需要编码区，帮我预过滤」
→ dgra-prefilter --input /path/to/sample.vcf.gz --preset coding-only

用户：「用最小调控区域过滤这个 VCF」
→ dgra-prefilter --input /path/to/sample.vcf.gz --preset regulatory-minimal

用户：「帮我更新参考数据」
→ dgra-prefilter --update-refs
