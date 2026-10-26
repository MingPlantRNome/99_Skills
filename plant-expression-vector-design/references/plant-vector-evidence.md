# 植物表达载体设计：定向文献证据集

核查截止日：2026-10-26。此为围绕本技能缺口的定向检索，不是穷尽式系统综述；下列为 11 篇原始研究、1 篇标准观点和1篇近期综述。逐项保留实验系统与外推限制，具体项目需重新核查物种、组织、构建版本和来源。

## E01：组装语法是标准，不是通用性能保证

**Patron NJ, Orzaez D, Marillonnet S, et al.** Standards for plant synthetic biology: A common syntax for exchange of DNA parts. New Phytologist, 2015;208:13–19. DOI:10.1111/nph.13532。类别为 Viewpoint/标准倡议，不是新的性能比较实验（[全文](https://repository.lsu.edu/cgi/viewcontent.cgi?article=4216&context=biosci_pubs)）。

- **内容与依据**：提出12个融合接口及Type IIS兼容准则，讨论BsaI和多级组装所涉其他酶、载体选择标记与非回文接口（[全文](https://repository.lsu.edu/cgi/viewcontent.cgi?article=4216&context=biosci_pubs)）。
- **强项/限制**：适合解释语法与元件边界，不证明所有工具包互换，也不证明去内部酶切位点的改写功能中性（[全文](https://repository.lsu.edu/cgi/viewcontent.cgi?article=4216&context=biosci_pubs)）。
- **技能动作**：登记工具包/级别/实际位点与接口，禁止未经确认改写调控或RNA保护区。

## E02：GreenGate 的实际构建、标签和表达背景

**Lampropoulos A, Sutikovic Z, Wenzl C, et al.** GreenGate - A Novel, Versatile, and Efficient Cloning System for Plant Transgenesis. PLOS ONE, 2013;8:e83043. DOI:10.1371/journal.pone.0083043（[全文](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0083043)）。

- **系统/方法**：Arabidopsis 稳定转化、N. benthamiana 瞬时表达；组装测序、共聚焦、荧光和WUS相关表型，见Figures 2–4（[全文](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0083043)）。
- **关键结果**：作者分别验证模块化构建和多表达框定位；含Ω元件的特定WUS方案未显示预期增强，组成型WUS方案可瞬时检出信号却未获得稳定植株，作者讨论了其功能背景（[全文](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0083043)）。
- **强项/限制**：有结构及植物读数；比较涉及具体模块/接口，不能得出所有GreenGate优于其他平台或某leader普遍更强。pGreen相关维持需求依赖所用具体骨架（[全文](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0083043)）。
- **技能动作**：检查 dummy/标签/stop 提供者、辅助维持要求、克隆留下的序列，并区分瞬时和稳定用途。

## E03：启动子与终止子的组合依赖

**Shakhova ES, Markina NM, Mitiouchkina T, et al.** Systematic Comparison of Plant Promoters in Nicotiana spp. Expression Systems. International Journal of Molecular Sciences, 2022;23:15441. DOI:10.3390/ijms232315441（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9740895/)）。

- **系统/方法**：N. benthamiana 叶片与N. tabacum BY-2 cell packs的瞬时EGFP，含p19；比较19个启动子、部分启动子–终止子组合及输入剂量（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9740895/)）。
- **关键结果**：在pdel5_MtHP背景，更换终止区域带来超过50倍的表达差异；终止子影响程度随启动子而变，输入剂量也非在所有组合中线性对应信号（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9740895/)）。
- **强项/限制**：组合和两种检测体系有直接比较；不覆盖任意植物、稳定转化、所有目标蛋白，不能抽成固定终止子排名（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9740895/)）。
- **技能动作**：比较完整表达组合，记录p19/输入/归一化背景，不默认最强组成型表达。

## E04：终止区域边界、read-through 与小RNA

**Fenselau de Felippes F, Shand K, Waterhouse PM.** Identification of a Transferrable Terminator Element That Inhibits Small RNA Production and Improves Transgene Expression Levels. Frontiers in Plant Science, 2022;13:877793. DOI:10.3389/fpls.2022.877793（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9149433/)）。

- **系统/方法**：N. benthamiana GFP瞬时检测及Arabidopsis稳定群体；tHSP删除/转移片段、RT-qPCR read-through、small RNA、RT-PCR和cRT-PCR（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9149433/)）。
- **关键结果**：删除tHSP前32nt减少GFP并增加read-through；向一些终止子转移片段可改善表达，但位置和拷贝效应依赖终止子，两个拷贝一般并非额外增强。多数比较的poly(A)位点分布没有明显改变，不能把全部效果归为同一个poly(A)机制（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9149433/)）。
- **强项/限制**：连接RNA与蛋白层读数，也有稳定群体结果；不证明任意终止片段、间距或重复次数都有效（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9149433/)）。
- **技能动作**：保留精确3′区域及间距；必要时列read-through、加工和沉默相关验证，不擅自截短或叠加终止片段。

## E05：相邻表达框的增强影响

**Chamness JC, Kumar J, Cruz AJ, et al.** An extensible vector toolkit and parts library for advanced engineering of plant genomes. The Plant Genome, 2023;16:e20312. DOI:10.1002/tpg2.20312（[全文](https://onlinelibrary.wiley.com/doi/10.1002/tpg2.20312)）。

- **系统/方法**：N. benthamiana叶片，双荧光素酶与反向相邻表达框；§3.2.3/Figure 4测试选择标记启动子对GmGy1/GmHSP的作用（[全文](https://onlinelibrary.wiley.com/doi/10.1002/tpg2.20312)）。
- **关键结果**：2x35S与CmYLCV在该配置下显著激活GmGy1；GmHSP未出现相同显著结果，所测管家启动子也未显示相同激活。on-target强度并不能简单预测off-target激活（[全文](https://onlinelibrary.wiley.com/doi/10.1002/tpg2.20312)）。
- **强项/限制**：直接考察跨框效应；样本是特定方向、距离和启动子组合，不支持所有病毒启动子必然串扰（[全文](https://onlinelibrary.wiley.com/doi/10.1002/tpg2.20312)）。
- **技能动作**：多表达框要做整合组合/背景控制，参考报告基因不预设恒定，不能把元件各自可用当组合已验证。

## E06：intron 增强和宿主加工并不通用

**Morello L, Gianì S, Troina F, Breviario D.** Testing the IMEter on rice introns and other aspects of intron-mediated enhancement of gene expression. Journal of Experimental Botany, 2011;62:533–544；在线2010. DOI:10.1093/jxb/erq273（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3003800/)）。

- **系统/方法**：水稻愈伤组织、Arabidopsis细胞及烟草；同位置替换intron，GUS/GFP与RT-PCR检测，包含剪接位点变体（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3003800/)）。
- **关键结果**：水稻匹配背景的增强约4–17倍；部分intron转到双子叶体系后加工和表达不同，rice-trained IMEter评分未解释实际增强差异（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3003800/)）。
- **强项/限制**：同背景比较并检测RNA；特定leader和侧翼外显子不能代表任意intron位置或所有植物（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3003800/)）。
- **技能动作**：复刻保留原intron，新增时明确位置/侧翼/宿主并验证加工，预测分数不代替实验。

## E07：无intron cDNA也可能发生宿主异常剪接

**Top O, Milferstaedt SWL, van Gessel N, et al.** Expression of a human cDNA in moss results in spliced mRNAs and fragmentary protein isoforms. Communications Biology, 2021;4:964. DOI:10.1038/s42003-021-02486-3（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8361020/)）。

- **系统/方法**：Physcomitrium patens，外源FIX/FH构建；RT-PCR和转录本测序、免疫印迹、质谱及Citrine成像（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8361020/)）。
- **关键结果**：原FIX cDNA出现异常剪接和蛋白片段；局部位点改写尚有异常，组合改写/优化版本仅检出全长产物。该方案同时改变剪接邻域和同义编码，不能把因果单独归于“密码子最优”（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8361020/)）。
- **强项/限制**：RNA、蛋白及质谱证据互补；苔藓外源人基因不代表所有植物天然CDS，也不确定统一GC阈值（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8361020/)）。
- **技能动作**：异源CDS按宿主评估加工，DNA翻译正确不等于真实RNA正确；不全面改写RNA机制实验序列。

## E08：2A 分离效率依赖具体融合背景

**Burén S, Ortega-Villasante C, Ötvös K, et al.** Use of the Foot-and-Mouth Disease Virus 2A Peptide Co-Expression System to Study Intracellular Protein Trafficking in Arabidopsis. PLOS ONE, 2012;7:e51973. DOI:10.1371/journal.pone.0051973（[全文](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0051973)）。

- **系统/方法**：Arabidopsis原生质体与N. benthamiana叶片，定位信号/上游蛋白内容变体，免疫印迹比较分离和未分离产物（[全文](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0051973)）。
- **关键结果**：所测CFP-2A背景由胞质约20%分离变为ER靶向约60%；改变上游内容又出现不同效率，其他构建也非一致结果（[全文](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0051973)）。
- **强项/限制**：直接量化产物；不是所有2A、顺序或定位的通用效率，不证明各最终蛋白等量/等功能（[全文](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0051973)）。
- **技能动作**：独立模块、2A与单蛋白融合分支隔离，按构建检测各产物、残留和功能。

## E09：稳定插入结构须独立验证

**Jupe F, Rivkin AC, Michael TP, et al.** The complex architecture and epigenomic impact of plant T-DNA insertions. PLOS Genetics, 2019;15:e1007819. DOI:10.1371/journal.pgen.1007819（[全文](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007819)）。

- **系统/方法**：Arabidopsis T-DNA株系的光学图谱、ONT长读长、RNA/small-RNA、bisulfite和染色质分析（[全文](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007819)）。
- **关键结果**：所测株系出现完整/部分T-DNA及骨架串联与重排，转基因和局部表观状态也不一致；基因组插入不能仅由原载体图谱推定（[全文](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007819)）。
- **强项/限制**：结构和组学互证；所测少量具体株系不能给用户任意项目提供发生概率或统一插入结构（[全文](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007819)）。
- **技能动作**：稳定任务另列完整插入/边界、拷贝、骨架和遗传验证；不要求每项目必做同一平台，但必须选择能回答目标问题的证据。

## E10：5′leader 与 uORF 需要配对功能测试

**Haas BE, Saif F, Bolger E, et al.** A Fluorescence-Based Transient Expression Assay for the Analysis of Upstream Open Reading Frames in Plants. Plant Direct, 2026;10:e70163. DOI:10.1002/pld3.70163（[全文PDF](https://ripe.illinois.edu/sites/ripe.illinois.edu/files/2026-04/Plant%20Direct%20-%202026%20-%20Haas%20-%20A%20Fluorescence%E2%80%90Based%20Transient%20Expression%20Assay%20for%20the%20Analysis%20of%20Upstream%20Open%20Reading.pdf)）。

- **系统/方法**：N. benthamiana双荧光瞬时报告，测试Arabidopsis、莴苣、大豆和豇豆leader，结合部分RT-qPCR；参考框采用OCS而非复制同一35S背景（[全文PDF](https://ripe.illinois.edu/sites/ripe.illinois.edu/files/2026-04/Plant%20Direct%20-%202026%20-%20Haas%20-%20A%20Fluorescence%E2%80%90Based%20Transient%20Expression%20Assay%20for%20the%20Analysis%20of%20Upstream%20Open%20Reading.pdf)）。
- **关键结果**：某些uORF起始位点变体提高报告比值，另一些没有显著提高；研究也测试ACG/TTG候选，不能把ATG-only扫描当完整uORF分析（[全文PDF](https://ripe.illinois.edu/sites/ripe.illinois.edu/files/2026-04/Plant%20Direct%20-%202026%20-%20Haas%20-%20A%20Fluorescence%E2%80%90Based%20Transient%20Expression%20Assay%20for%20the%20Analysis%20of%20Upstream%20Open%20Reading.pdf)）。
- **强项/限制**：leader配对变体和独立参照；在烟草瞬时体系而非各来源作物稳定编辑，报告变动不等于天然基因功能已改变（[全文PDF](https://ripe.illinois.edu/sites/ripe.illinois.edu/files/2026-04/Plant%20Direct%20-%202026%20-%20Haas%20-%20A%20Fluorescence%E2%80%90Based%20Transient%20Expression%20Assay%20for%20the%20Analysis%20of%20Upstream%20Open%20Reading.pdf)）。
- **技能动作**：leader是显式部件；序列候选只标待验证，禁止自动“清除所有上游ATG”。

## E11：骨架 ORI 配对不是默认中性

**De Saeger J, Vermeersch M, Aesaert S, et al.** BackBone Builder (B3): a modular Golden Gate standard with a compatible Agrobacterium parts library. Plant Physiology, 2026;201(4). DOI:10.1093/plphys/kiag610（[全文](https://academic.oup.com/plphys/article/201/4/kiag610/8772012)）。

- **系统/方法**：16种E. coli/Agrobacterium ORI组合的结构验证与细菌报告，N. benthamiana RUBY检测，另有B104玉米稳定转化及digital PCR事件质量检测（[全文](https://academic.oup.com/plphys/article/201/4/kiag610/8772012)）。
- **关键结果**：配对ORI影响所测细菌与叶片读数，部分配对出现明显下降；特定B3骨架可完成玉米转化，但未建立任意骨架的普遍最佳排名（[全文](https://academic.oup.com/plphys/article/201/4/kiag610/8772012)）。
- **强项/限制**：组合骨架、瞬时与稳定用途均有验证；有限背景不决定其他目标的优劣，本文不据此向用户推荐通用ORI或直接套用转化率（[全文](https://academic.oup.com/plphys/article/201/4/kiag610/8772012)）。
- **技能动作**：库存优先需附版本、复制/维持、菌株和用途兼容检查，不因同是二元载体就认为等价。

## E12：柔性 linker 不是稳定性的保证

**Gattinger P, Izadi S, Grünwald-Gruber C, et al.** The Instability of Dimeric Fc-Fusions Expressed in Plants Can Be Solved by Monomeric Fc Technology. Frontiers in Plant Science, 2021;12:671728. DOI:10.3389/fpls.2021.671728（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8299721/)）。

- **系统/方法**：N. benthamiana分泌型融合，免疫印迹、纯化、native PAGE和质谱断点分析，含不同域/顺序/连接版本（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8299721/)）。
- **关键结果**：该研究的(G4S)2融合仍出现明显片段，改变结构背景后结果不同；不能将linker、标签位置或单一序列视为全部降解原因（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8299721/)）。
- **强项/限制**：有全长与碎片等多层读数；特定分泌蛋白研究不能预测任意植物蛋白，机制解释仍含假说（[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8299721/)）。
- **技能动作**：保留用户要求的柔性连接肽，但用全长/片段和功能对照验证，不承诺功能中性或稳定提高。

## E13：近期互操作性综述，仅作导航

**Gélinas Bélanger J, Hanaishi S, Péret B.** The plant Golden Gate toolkit: compatibility, interoperability, and harmonized assembly in plant synthetic biology. Plant Methods, 2026. DOI:10.1186/s13007-026-01563-0（[全文](https://link.springer.com/article/10.1186/s13007-026-01563-0)）。

- **性质/用途**：综述比较工具包和标准，并非新的跨所有平台实测；兼容性需回到具体版本、接口与所引原研究核验（[全文](https://link.springer.com/article/10.1186/s13007-026-01563-0)）。
- **技能动作**：作为查找新工具包的入口，不替代原始质粒序列、实际语法或功能实验。

## 检索限制

另查阅 Engler et al. 2014 的 MoClo Plant Parts Kit 页面，其说明了96个标准化部件和补充文件入口，但本次页面未提供完整方法/性能数据，故未把它作为新的定量比较证据（[原始页面，DOI:10.1021/sb4001504](https://pubs.acs.org/doi/full/10.1021/sb4001504)）。

不同论文系统、目标、p19/输入、检测归一化及转化方式不同，不做跨文献原始倍数排行榜。未取得的序列、未读到的补充实验和不可确认参数写NA，不用本证据集填补具体项目的数据。
