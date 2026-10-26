# 序列清单与本地工具

## 环境与入口

在独立环境安装 `requirements.txt`，不要为此升级 HPC 的 Evo2 环境。脚本不下载数据库，不在云端加载 `.pt`，也不执行自动采购。

```bash
python -m pip install -r requirements.txt
python scripts/vector_tools.py --help
python scripts/evo2_job.py --help
```

## 构建清单

以 `templates/construct.example.json` 为结构样板，按 `templates/construct.schema.json` 校验；其中空数组/待填值不构成可用设计。真实原始序列文件单独保留，路径相对 manifest。

- `sources`：id、文件、format、可选 record_id/预期序列哈希，以及包含 URL、版本、边界依据、reviewed 的 evidence。多记录文件必须指定 record_id。实验室直接提供的原文件使用 user_provided，不强制伪造公共 URL/accession。
- `assembly`：按最终顺序列出来源片段，坐标均为 1-based inclusive，strand 为 1/-1。不得用缺失序列占位。
- `features`：显式最终坐标；`segments` 按生物学转录/翻译顺序给出，支持跨环状原点和负链。CDS 必须有 expected_protein 与 complete_cds；内含子不能混进连续 CDS。
- `new_expression_target_ids`：显式列出本次新设计的目标蛋白 CDS，不包括骨架筛选/维持基因或标签组件本身。de_novo 必填；部分重建的新产物同样登记，避免误向用户询问是否给抗性基因加标签。
- `tag_decisions`：按 `templates/tag-decision.example.json` 为上述目标逐个记录用户标签决定，pending/缺少用户决定不能最终导出；添加时须有已确认方案、文献证据 ID、精确序列来源及位于目标 CDS 内的标签注释。新增融合还须填写 `templates/fusion-plan.example.json` 所示的 fusion_plan。原有复刻 manifest 可不新增该字段；改变原标签需另记。
- `protected_regions`：同样的坐标结构，用于候选保护；可提供应匹配的 DNA 字符串。
- `restriction_sites`：名称、ACGT 识别序列及 must_be_absent；检查包含 DNA 双链和环状跨原点，但不模拟完整酶切反应。
- `replication_reference_source_id`：声明 exact_plasmid 时必填，要求构建与基准全长等价；软件无法代替文献身份核验。
- `unresolved_issues`：任何非空项都阻断正式导出。

完整片段会继承 GenBank 中包含在片段内的注释；边界截断的注释不自动补造，需重新检查。实际剪接、细胞功能和构建成功不由本工具验证。

## 构建与导出

```bash
python scripts/vector_tools.py build \
  --manifest project/construct.json --out project/design-v1
```

输出目录必须不存在或为空。失败只写 QC，不导出伪完整载体；通过后生成 `construct.gb`、`construct.fasta`、`features.tsv`、`sources.tsv`、`plasmid_map.svg`、`expression_map.svg`、`qc_report.json` 与说明。图为数据驱动的线性全长图和表达元件图，不是原生 SnapGene 文件或环状美术图。

完整序列仅接受 ACGT；N、RNA U、空序列和坐标错误被拒绝。哈希按规范化大写 DNA 计算，原文件不修改。CDS 检查包括翻译、完整起止要求及预期蛋白；不支持把计算预测剪接称为已验证。

标签判断检查登记的决定、来源字段和注释归属。新融合的 component 必须按翻译方向覆盖整条 CDS 的编码区（末端 stop 除外）、保持密码子边界，并与实际来源按声明的 stop 处理后完全相符；目标/标签之间必须有已注释柔性连接肽，连接肽 DNA 翻译须等于申报肽段。最终 stop 与来源记录核对，整体 CDS 再做完整翻译。

计算一致性不证明用户决定/论文结论真实，不识别全部未声明标签，也不证明连接肽在该蛋白中的柔性或功能效果。它不自动删改碱基；失败应先修正设计、来源记录或注释，不能隐去错误。复杂非连续 ORF/特殊翻译机制需专门规则，不当作普通融合处理。

测试 fixture 只在 `--allow-test-fixtures` 明确指定时允许，输出标为测试用途，不能作真实构建。

## 等长候选检查

```bash
python scripts/vector_tools.py compare \
  --manifest project/construct.json \
  --baseline project/design-v1/construct.fasta \
  --candidate project/candidate.fa --out project/candidate-qc.json

python scripts/vector_tools.py apply-edits \
  --manifest project/construct.json \
  --baseline project/design-v1/construct.fasta \
  --edits project/edits.json --out project/candidate.fa
```

`edits.json` 是列表，例如 `[{"pos": 12, "ref": "A", "alt": "G"}]`，位置为 1-based；这只是格式例子，不是推荐改动。仅允许等长、不重叠、REF 完全匹配的 DNA 替换，且必须通过受保护区域、CDS 和 must_be_absent 位点检查。

不同长度候选返回不可比/需新坐标，不能沿用原版坐标检查 intron 或其他插入删除。手工新增 CDS 的缺失 expected_protein 不能绕过。

## SVG 与跨文件一致性

全长图使用编号、方向及图例避免标签挤压，图例给真实坐标；表达图从同一 feature 数据生成。最终写出后重新读入 GenBank/FASTA，检查序列、拓扑、feature 坐标、方向和提取序列一致，未通过不保留最终导出。

## 可选的只读序列审查

按 `references/evidence-driven-review.md`，在完整GenBank上运行 `scripts/sequence_review.py`。其ATG-only leader扫描和精确调控片段比较只输出可复核观察；不补猜TSS/leader，不覆盖非AUG、近似重复或实际加工/表达，不自动改写序列。缺少所需注释标未评估，不可报告“无风险”。
