# 标签融合：终止密码子与柔性连接肽

新增或改变标签融合时读取本文件。用户要求目标与标签之间加入柔性连接肽，并确定各片段的终止密码子怎样处理；不能只拼接两个名义 CDS，也不统一删除所有出现的 TAA/TAG/TGA。

## 用明确的 N 到 C 架构消除歧义

对普通连续 ORF，先画出从 N 到 C 的蛋白顺序，再处理来源片段：

- **目标–柔性 linker–标签**：目标片段末端若有 stop，必须从融合插入片段移除；来源本就无 stop 时记为原本没有。完整融合产物的 stop 在最后的标签/尾部模块之后。
- **标签–柔性 linker–目标**：检查上游标签来源片段的 stop，不能让它提前终止；目标是最后模块时，其 stop 可作为完整融合的最终 stop。若目标后还有其他模块，重新按整体顺序判断，不能机械保留。
- **目标 N 部分–linker–内部标签–linker–目标 C 部分**：两侧都按接头审查，内部无终止、保持读框，终止放在完整产物末端。

Addgene 的供稿说明明确提醒 C 端融合要去掉上游蛋白 stop，并核对下游荧光蛋白的读框；该具体载体是线虫体系，不能据此推断其植物适用性（[KG#45 原始质粒说明](https://www.addgene.org/110926/)）。

每个片段分别记录来源 DNA、方向/坐标、原始末端三联体、是否真的存在 stop，以及 `remove_for_fusion`、`none_present` 或 `retained_as_final`。只针对已识别的末端 stop 改动，不做全序列搜索替换。最终 stop 可来自最后模块、骨架或明确的新设计，必须给实际来源。

以上指普通连续 ORF，不把 2A、IRES、双表达框、特殊读穿或特殊编码规则当作同一种融合；不支持的架构应另行设计，不绕过本检查。

## 柔性连接肽是独立设计元件

对新增目标–标签接头，加入独立的柔性连接肽注释，记录肽段、DNA、氨基酸长度、来源/编码方案、柔性依据及相关文献。富 Gly/Ser 的候选可以检索比较，但不固定某个序列或重复次数，也不默认越长越好。

植物原始研究使用了 `(GGGGS)2` 并验证了相应构建的表达、纯化/切割及功能，但没有直接比较无连接肽或不同长度，不能把整个构建成功归因于该 linker（[Islam et al., Plant Biotechnology Journal, DOI:10.1111/pbi.13040](https://pmc.ncbi.nlm.nih.gov/articles/PMC6523591/)）。

另一植物研究的免疫印迹显示，带 `(GGGGS)2` 的二聚 Fc 融合仍有明显降解，说明该实例中加入柔性 linker 并不保证稳定；不要因此泛化到任意蛋白（[Gattinger et al., Frontiers in Plant Science, 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8299721/)）。

实际项目仍按标签规则检索相近蛋白/植物体系，区分直接比较、构建可行性和间接外推。对照和功能读数应能判断具体目标；不能以“序列富 Gly/Ser”或 Evo2 高分证明该融合功能不受影响。

忠实复刻时不擅自改写作者原始接头。原文没有柔性 linker 而用户希望加入时，保留原版并另列适配候选、差异和确认项，不把改造版称为完整复刻。

## fusion_plan 数据约定

以 `templates/fusion-plan.example.json` 为模板，嵌入相应 tag_decisions 条目的 fusion_plan。模板只是结构，不能把占位或空 components 当作完成设计。

- `assembly_mode`：当前为 continuous_orf。
- `components`：按翻译方向排列。每项引用独立最终 feature、角色、原始 source_id/source_segments/source_strand、stop_action 和 expected_peptide。CDS/片段坐标沿用 1-based inclusive；joined/负链均按生物学顺序。
- `role`：target、tag、flexible_linker、signal 或 other。标签组件须与 tag_feature_ids 一致，目标/标签两类之间每次切换都要有柔性 linker。
- `stop_action`：来源末端有 stop 时，内部模块必须 remove_for_fusion；只有最后模块可 retained_as_final。原本没有 stop 用 none_present，不能伪称删除。
- `expected_peptide`：该组件编码体，不含末端 stop；必须与来源按上述处理后的 DNA 及最终注释一致。
- `linker_flexibility_basis`、`linker_evidence_ids`：记录为什么把该组件作为柔性连接肽及文献依据，不是通用功能保证。
- `final_stop`：准确 codon、provider（last_component/vector/designed）和对应 source_id/source_segments/source_strand，不能只写“已加 stop”。

整体 CDS 的 expected_protein 包含全部目标、linker、标签、信号和克隆留下的实际氨基酸；信号或位点等也必须显式成为组件。不要无依据移除内部目标起始 Met，任何必要改变应有来源/确认。

## 计算验收与失败处理

脚本要求组件按翻译方向、密码子边界连续覆盖完整 CDS 的编码体，末端只有所声明的最终 stop；检查源片段 stop_action、来源/最终组件 DNA 一致、连接肽翻译、N/C/内部位置与方向，随后再验证完整融合翻译。

缺少组件来源、错误的 stop 处理、缺失柔性 linker、接头移码、内部 stop 或与预期产物不符时阻断最终导出。脚本不自动删改碱基；先修正清单/序列并重新检查，必要时询问用户。

原样复刻的既有标签条目可不强制补录新 fusion_plan，但不能据此声称已完成逐来源接头审计。新增或改造标签必须执行；同时保留原始输入、改动、来源证据和验证边界。
