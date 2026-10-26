# 植物表达载体设计技能

版本：v1.5.0。默认以文献复刻为先，只有必要适配、明确优化或经检索仍缺失的部分才提出新设计，并保留原版与完整来源记录。

## 安装与文件入口

- [v1.5 可安装 ZIP](dist/plant-expression-vector-design-v1.5.zip)：根目录直接包含 `SKILL.md`，可作为完整技能包导入。
- [主技能](SKILL.md)：执行流程和按需读取入口。
- [文献证据集](references/plant-vector-evidence.md)：11篇原始研究、1篇标准观点和1篇综述，附实验、限制和原始链接。
- [设计审查规则](references/evidence-driven-review.md)：R01–R12审查。
- [图文报告模板](templates/design-report.md)：包含标签/接头、边界、未解问题与优化计划；附图为格式示例，真实项目须重绘。

## 核心行为

1. 从论文、附件、NCBI和原始质粒记录追溯精确序列与版本，优先适用的实验室骨架。
2. 每次调用重新检查 `LAB_VECTOR_DIR`，源库只读；未连接Mac不能当作目录不存在。
3. 从头设计先问是否添加标签，需要时据文献比较候选，审核STOP、柔性连接肽、接头和完整融合翻译。
4. 输出GenBank、FASTA、注释、图谱和QC，缺少关键序列或检查失败不导出伪完整载体。
5. 复用用户现有HPC EasyEvo2环境，`.pt`路径在实际任务中再询问；配置未验证或未运行时不生成模型分数。
6. 审查leader/终止区域、元件组合、多表达框、组装/骨架兼容及分阶段验证，区分序列事实、模型偏好和功能实验。

## 本地软件检查

在独立Python环境安装序列工具依赖并运行测试，不为此改变HPC的Evo2部署：

```bash
python -m pip install -r requirements.txt pytest
python -m pytest tests -q
```

可选PNG生成依赖见 `requirements-report.txt`，还需要系统Cairo和中文字体。SVG输出与序列工具不要求安装PNG依赖。

v1.5在发布前通过148项本地软件测试和技能格式校验；测试不等于实际HPC接入、植物表达或功能验证。新序列审查仅做已注释leader的ATG候选和精确调控片段重复观察，不自动删除序列，也不预测功能或风险概率。

## 私有配置边界

本包包含用户明确指定的Mac载体目录字符串，见 `templates/lab-library.config.json`，但不包含该目录中的任何载体文件。没有打包模型权重、真实HPC运行配置或凭据；后续私有项目输出和权重不应提交到本技能目录。

源码按 `references/`、`scripts/`、`templates/` 和 `tests/` 组织。`dist/` 内的ZIP是已验证的v1.5分发快照；若日后修改源码，应重新测试和打包，避免源码与安装包版本不一致。
