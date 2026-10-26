# 文献驱动的完整设计审查

所有载体任务在架构确定前读取本文件，交付前重查；逐项填写 `templates/design-review.tsv`。详细证据见 `references/plant-vector-evidence.md`，记录检索截止日；这些文献不是对新项目的性能保证。

## 项目级判定

先声明宿主/基因型、组织/时期、瞬时或稳定、核表达或专用体系、目标是功能保真/定位/诱导/机制还是产量。保留用户的原文复刻版，任何改变的依据和代价另记。没有兼容证据时不直接转用另一植物体系的最佳部件。

强度不是唯一目标：GreenGate 中瞬时检测到荧光的 WUS 构建未获得稳定转化植株，文中将其与组成型 WUS 的已知效应联系；这提示具体功能目标和表达时空背景必须单独审查，不能把瞬时信号当稳定应用成功（[Lampropoulos et al., 2013](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0083043)）。

## 必查项目

| ID | 要检查什么 | 必须形成的决策或验证 |
|---|---|---|
| R01 | 宿主、组织、瞬时/稳定、表达目的 | 对比原文条件；不明项写缺口，不默认最大表达 |
| R02 | 精确启动子边界、TSS、5′UTR/leader、主起始环境 | 区分元件与编码接头；未知边界不靠任意上下游长度补齐 |
| R03 | leader 中的 uAUG/uORF 与非 AUG 可能性 | 候选仅为筛查，不自动删除；必要时用实际 leader 变体与参考报告系统验证 |
| R04 | 3′UTR、终止/加工区域的精确版本及边界 | 不把翻译 STOP 等同于转录终止；保留原始终止片段及间距，必要时检测 read-through/3′端 |
| R05 | 启动子–leader–CDS–终止区域组合 | 在配对背景测试候选，不照抄其他背景中的强度排名 |
| R06 | intron、密码子、宿主剪接与 RNA 保护 | DNA 翻译通过不证明 RNA 加工正确；异源/改写 CDS 可有隐性剪接风险，按问题加入 RNA 层读数 |
| R07 | 多表达框方向、邻近增强影响、重复调控序列 | 同时检查独立框和整合组合；精确重复只作事实标记，不宣称已沉默或已重组 |
| R08 | 独立表达框、2A 或融合机制及检测内参 | 2A 需检查各产物和残留/未分离产物；不保证等量；内参稳定性与串扰单独验证 |
| R09 | 组装体系/级别、Type IIS 位点、overhang 和留下的序列 | 使用该工具包实际语法与序列，不把“同属 Golden Gate”当兼容；不擅改调控/RNA 区去位点 |
| R10 | 骨架复制/维持、菌株/辅助质粒、细菌与植物筛选 | 检查真实版本和兼容性；骨架不是天然中性，未知不推断 |
| R11 | 标签、定位、接头、STOP 与柔性连接肽 | 执行现有标签询问和融合检查；保持功能对照，不以信号证明功能不受影响 |
| R12 | 质粒结构、RNA、蛋白、功能、稳定插入/遗传 | 分阶段读数；稳定体系另查边界/拷贝/骨架/分离等适用指标，多事件比较而非单株外推 |

每一行给 `review_status`（reviewed/gap/not_applicable）、当前证据、具体适用性、动作、是否阻断哪项交付及待解决问题。not_applicable 需原因，不能整表默认通过。若设计仅瞬时表达，稳定插入/遗传子项可以不适用，但不能因此声称适用于稳定转化。

## 关键证据对应规则

- **组合而非单元件排名**：在 Nicotiana 瞬时系统，选定启动子背景下更换终止子可产生很大表达差异，强弱随组合而变；记录完整组合和检测条件，不永久指定一个最优终止子（[Shakhova et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9740895/)）。
- **3′端不是装饰**：tHSP 删除/转移片段的 GFP、read-through、small-RNA 与 RT-PCR 结果说明该研究中的位置和序列重要；不能随意截短终止片段，也不把其增强元件复制次数当普遍剂量规则（[de Felippes et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9149433/)）。
- **相邻表达框要一起测**：特定 head-to-head 配置下，2x35S/CmYLCV 激活 GmGy1，但不是所有候选或弱启动子都有此效应；不把它泛化为所有 35S 组合必然串扰（[Chamness et al., 2023](https://onlinelibrary.wiley.com/doi/10.1002/tpg2.20312)）。
- **uORF 需要功能证据**：2026 年 leader 变体在烟草报告系统中表现不一，部分非 AUG 候选也被测试；不能把一次 ATG 扫描当全 uORF 检测，更不能把所有上游起始位点一律删掉（[Haas et al., 2026](https://ripe.illinois.edu/sites/ripe.illinois.edu/files/2026-04/Plant%20Direct%20-%202026%20-%20Haas%20-%20A%20Fluorescence%E2%80%90Based%20Transient%20Expression%20Assay%20for%20the%20Analysis%20of%20Upstream%20Open%20Reading.pdf)）。
- **RNA 读数独立**：水稻/拟南芥的 intron 效果及苔藓异源 cDNA 的异常剪接都依赖背景；保留原版并按问题检验成熟 RNA，不能承诺加 intron 或提高 GC 即能解决（[Morello et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3003800/)，[Top et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8361020/)）。
- **2A 不等于无条件完全分离**：原始研究中胞内定位与上游编码内容影响分离率，需检查每个具体构建；2A 架构不得冒充单蛋白柔性 linker 融合（[Burén et al., 2012](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0051973)）。
- **语法/骨架需版本化**：共用语法是标准提案，不是所有试剂互换证明；B3 的实验也表明所测 ORI 配对与读数相关，而非任意骨架可替换（[Patron et al., 2015](https://repository.lsu.edu/cgi/viewcontent.cgi?article=4216&context=biosci_pubs)，[De Saeger et al., 2026](https://academic.oup.com/plphys/article/201/4/kiag610/8772012)）。
- **验证转化后的真实结构**：少量 Arabidopsis 插入株系的长读长/光学图谱显示复杂 T-DNA 与骨架结构；不能将质粒测序正确、筛选阳性或一段 PCR 外推为完整单拷贝事件，也不将该样本频率外推给用户项目（[Jupe et al., 2019](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007819)）。

## 可执行但有限的序列观察

完整 GenBank 输出后，可运行：

```bash
python scripts/sequence_review.py --input project/design-v1/construct.gb --out project/review-v1
```

输出 `sequence_review.json` 与 `review_flags.tsv`：只筛查已注释 5′UTR 的 ATG 候选和已注释调控元件的精确重复。未注释/方向不明即未评估；不覆盖非 AUG、近似重复、全基因组同源、RNA 二级结构或功能效应，也不生成“风险概率”。它不修改任何序列，不能替代 R01–R12 的人工/代理证据审查。

## 交付门槛与验证层级

1. 序列身份、边界、读框或组装结构的关键缺口同步写 manifest 的 unresolved_issues，阻断相应最终序列导出。
2. 功能、强度、剪接、分离、稳定性等未验证项保留在边界/验证计划中；仍可交付清楚标注的计算设计，不能宣称已经改善。
3. 最小验证按目的设置 DNA→RNA→蛋白/加工→功能→适用的稳定遗传层级；每层读数、阴阳性/无标签等关键对照、失败后回退分别写明。
4. 保留无改动复刻基线。新建议不足时优先补证，不强迫优化；Evo2 结果只作模型偏好，不能替代上述层级。
5. 报告按图文风格区分“构建流程箭头”和“生物学因果证据”，并保留用户要求的局限、未解问题和可否决优化计划。
