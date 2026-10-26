---
name: plant-expression-vector-design
description: "Design and reproduce plant expression vectors from papers, supplements, gene IDs, NCBI and plasmid repositories. Use for 植物表达载体、文献复刻、多基因表达框、标签选择、intron、密码子优化、载体图谱 and HPC EasyEvo2 evaluation. Prioritize reported designs and laboratory backbones. Not for general biology explanations or unconstrained sequence generation."
compatibility: "Sequence utilities require Python 3.10+, Biopython and jsonschema. LAB_VECTOR_DIR is checked on the user's Mac each invocation; the check script uses only Python's standard library. Evo2 reuses the existing HPC environment with a runtime .pt path. Missing device/network/HPC access must be reported."
metadata:
  version: "1.5.0"
  language: "zh-CN"
---

# 植物表达载体设计

## 工作原则

默认复刻文献实际使用的设计，不默认追求最高表达。优先级：原报道精确复刻、最小必要适配、检索后仍缺失部分的新设计；用户要求的优化另存候选，不覆盖复刻版。

主动读取全文和附件，按基因号到 NCBI 获取真实序列与注释；不把查序列的常规工作交回用户。名称、示意图、蛋白序列或模型评分都不能替代原始 DNA。找不到、访问失败、编号歧义和版本冲突分开记录。

保留作者使用的转录本、启动子边界、CDS、UTR、intron、终止区域、标签、接头、方向和多基因机制。作者已有密码子优化时查找其实际版本，不以天然 CDS 或重新优化版本冒充复刻。涉及 RNA 结构/剪接时保护原始核酸；RNA 输入只处理提供方向。

实验室已有且兼容的同版本骨架优先；更换骨架必须标为适配，不称完整质粒复刻。未知项写 NA。计算设计、实物测序和功能验证独立分层，不捏造数据或实验结果。

每次调用先读取 `references/lab-vector-directory.md`，解析 `LAB_VECTOR_DIR` 并在用户 Mac 上重新检查存在性、目录类型和读取权限。默认路径只存于 `templates/lab-library.config.json`；不复用旧检查，不把云端查不到当作 Mac 目录不存在。检查失败不自动建目录、换路径或下载骨架绕过；可以继续文献工作，但库存可用性保持未确认。

## 执行清单

1. **检查库并建立项目**：先按目录规则运行 `scripts/check_lab_vectors.py` 或等价的本机只读探测，记录本次路径、时间、设备和状态；无权访问本机标未执行。通过后才在该库做有界文件检索，不默认扫描全部内容。填写 `templates/project-brief.md`，确认物种/组织、瞬时或稳定、核或细胞器、目标文献和骨架清单；关键歧义列为阻断项。
2. **追溯证据**：读取 `references/evidence-and-retrieval.md`。逐篇重建构建；用 `templates/evidence.tsv`、`templates/sequence-registry.tsv`、`templates/retrieval-log.tsv` 记录全文/附件位置、支持实验、来源、版本和缺口。完成条件：每个元件有来源或明确缺失状态。
3. **审查架构与标签**：读取 `references/design-decisions.md` 和 `references/evidence-driven-review.md`，逐项填写 `templates/design-review.tsv`，不只比较启动子强度。新设计按 `references/tag-selection.md` 先问标签需求；已回答不重复问，未答不默认。需要时查原始研究并填 `templates/tag-decision.example.json`。复刻保留原标签/intron/密码子版本；用本次有效的 `LAB_VECTOR_DIR` 查骨架并填 `templates/lab-backbones.tsv`，不把目录可读当序列可用。
4. **序列级构建**：按 `references/sequence-tools.md` 建清单，新增标签执行 `references/fusion-junctions.md`，登记 fusion_plan 的来源、STOP、柔性连接肽并审查完整ORF。未确认/接头失败阻断最终导出。运行 `scripts/vector_tools.py build` 输出序列和图谱；有真实GenBank后运行 `scripts/sequence_review.py` 标记已注释leader的ATG及精确调控重复，不能据此自动删除/改写。缺序列只交草案。
5. **Evo2 评估/优化**：需要模型时读取 `references/evo2-hpc.md`。用户已确认其 EasyEvo2 可识别独立 Evo2-7B `.pt`；复用该环境，不要求下载或迁移缓存。到实际准备 Evo2 作业时再询问权重绝对路径，以 `weights_path` 参数传入；创建技能时不要求路径，也不写死路径。按 `templates/hpc-profile.example.json` 配置已核验调用，使用 `scripts/evo2_job.py` 准备、执行和收集。未配置/未运行标 NA。
6. **候选与验收**：先冻结原版。仅在许可区域改动；等长替换可用 `vector_tools.py apply-edits` 和 `compare` 检查蛋白/保护区域/位点。插入删除须另外建立候选注释与坐标，不沿用旧坐标。高模型分不能抵消 QC 失败。
7. **图文交付**：按 `references/delivery-and-validation.md`、`references/report-boundaries.md`、`references/report-visual-style.md` 填写报告。图示输入用 `templates/report-figures.example.json`，由 `scripts/render_report_figures.py` 生成概念图，不冒充序列图谱。报告必须包含适用边界、已知缺点/权衡、无法回答的问题及优化验证计划，给出证据、对照与阶段判据；关键序列缺口同步写入 unresolved_issues，不交伪完整载体。

## 工具入口与状态

本包含已实现的只读目录检查、本地导出/QC、等长候选检查、显式窗口准备及命令适配器；序列工具依赖见 `requirements.txt`，目录检查脚本只需标准库。使用命令前查看脚本 `--help`，格式和示例按对应参考文件执行。

`scripts/evo2_job.py` 的 `--weights-path` 是本技能自己的参数，不代表 EasyEvo2 有同名参数。根据用户环境的实际 help/脚本绑定到真实参数或环境变量，不猜 `--weights`。配置未验证时不能执行；远程准备不检查云端是否存在 HPC 权重，实际 HPC 执行时检查。

Evo2 仅提供模型偏好，不宣称表达倍数、剪接成功率或载体稳定性；无数据时不输出校准成功概率。缺少网络/HPC 权限则保留诊断并完成可做部分，不伪称已检索或已评分。不要自动重装 HPC 环境、提交无界搜索、购买材料或上传私有序列至外部服务。

## 示例

输入：“按这两篇论文复刻两个基因，优先实验室骨架，再用 EasyEvo2 比较。”

应先检查本次 `LAB_VECTOR_DIR`；可读后按名称有界检索骨架，再获得原构建/附件并按基因号核查 NCBI。列出原骨架与现有骨架差别，保持文献原 intron/优化 CDS，分别交付复刻与适配版。到 Evo2 阶段询问 `.pt` 路径和必要配置；实际调用未完成则报告未执行，不给示例分数。
