# EasyEvo2 / Evo2-7B：运行时权重与受约束评估

## 用户已确认的配置

用户在自己的 HPC 上部署 EasyEvo2，使用独立 Evo2-7B `.pt` 文件，现有环境能够识别该文件。接受此事实，不以公开上游默认接口否定用户环境，不要求重新下载或迁移 Hugging Face 缓存。

用户明确要求：`.pt` 地址作为后续运行参数，开始实际项目时再询问。使用字段 `weights_path`，初始值为 null；不在本技能包保存真实服务器地址、路径、凭据或模型权重。

## 首次运行配置

1. 到 Evo2 阶段询问：“本次使用的 Evo2-7B `.pt` 文件在 HPC 上的绝对路径是什么？”若用户已提供本次路径则不再询问。
2. 同轮补齐真正必要的现有运行命令/help、模型变体/版本、调度规则、资源上限与数据授权。保留调度器 GPU 可见性，不默认 Slurm。
3. 用 `templates/hpc-profile.example.json` 建立私有配置，填写已核验的 argv 和环境变量映射。`{weights_path}` 必须在 argv 或已确认的加载环境变量中被使用；不要猜 EasyEvo2 的权重 flag。
4. `verified` 只有实际检查部署后才置 true；`score_type` 也需核验，不凭列名 probability 推断。先单序列/no-change 基线和少量候选试跑。

公开上游提供 `score`/`vep`，但其参数和加载封装不一定与用户部署相同；参考上游代码只用于定位核验点（[EasyEvo2 命令](https://github.com/ylab-hi/EasyEvo2)，[score](https://raw.githubusercontent.com/ylab-hi/EasyEvo2/main/easyevo2/score.py)，[加载封装](https://raw.githubusercontent.com/ylab-hi/EasyEvo2/main/easyevo2/model.py)）。

## 本包的准备与执行

`--weights-path` 是本包脚本自己的参数，不声明为 EasyEvo2 原生参数。运行模板不使用 shell 字符串，不用 eval；参数按 argv 列表传递，支持含空格的实际路径。

```bash
python scripts/evo2_job.py prepare \
  --input project/candidates.fa \
  --profile project/hpc-profile.json \
  --weights-path "/用户届时提供的绝对路径/model.pt" \
  --window-size 512 --step-size 256 --out project/evo2-job
```

准备阶段只检查路径格式，不在云端检查远端文件。生成显式窗口 FASTA、输入哈希、参数和作业清单；不给缺失权重虚构默认路径。窗口覆盖尾端，但不是逐位点效应预测；同长度、同背景且相同窗口位置才能配对。

将任务目录与本包脚本放到已授权 HPC 环境，在获分配 GPU 的环境执行：

```bash
# 默认只展示经校验的执行计划，不运行模型
python scripts/evo2_job.py run --job project/evo2-job/job.json

# 实际运行需在 HPC、有资源和执行授权时明确指定
python scripts/evo2_job.py run --job project/evo2-job/job.json --execute

python scripts/evo2_job.py collect \
  --job project/evo2-job/job.json --baseline-id reported_baseline
```

模板 placeholders 仅允许 `{input}`、`{output}`、`{weights_path}`、`{model_type}`、`{window_size}`、`{step_size}`、`{device}`。不支持任意 shell 展开。可以在经核验 argv 前添加现有容器/启动器，但本包不自动 SSH、排队或重装部署；无访问权限则交付任务包待运行。

默认窗口由本脚本预先切出，避免不透明整质粒输入。公开 EasyEvo2 的窗口函数还会首先输出整条输入并且不补尾窗，收集时只接收任务登记的窗口 ID，不把额外行混进排名（[窗口实现](https://raw.githubusercontent.com/ylab-hi/EasyEvo2/main/easyevo2/utils.py)）。

## 评分规则

Evo2 官方评分为条件对数概率，支持 mean/sum 和可选反向互补平均；与任何表达量或剪接概率不是同一指标（[官方接口](https://raw.githubusercontent.com/ArcInstitute/evo2/main/evo2/models.py)，[评分实现](https://raw.githubusercontent.com/ArcInstitute/evo2/main/evo2/scoring.py)）。

- 固定检查点、环境、方向、BOS、padding、窗口背景与计分规则。当前适配器只收集明确核验为 mean_log_likelihood 的结果。
- DNA 默认以提供方向评分；若需要反向互补必须单独登记和对应。RNA 不自动反向扫描。
- 输出差值为候选减基线；不同长度、位置、上下文或模型不直接排名，标不可比。
- 检验同一序列不同 batch 的一致性，确认实际 token/计分位置与 padding 无误。默认包装器元数据不能证明这些正确，检查失败时不能推荐候选。
- 收集器会核对相同序列哈希的窗口评分在声明的 score_tolerance 内一致；这只是数值一致性检查，不是生物学效应阈值。
- 原始分数、可比性、QC 和功能校准状态分列；无校准时不写表达倍数或成功率。

原论文可控染色质设计结合了其他功能预测器且主要验证人/小鼠细胞，不能作为植物表达载体优化的直接证明（[Brixi et al., Nature, 2026](https://www.nature.com/articles/s41586-026-10176-5)）。

## 优化与降级

冻结原复刻版，只比较许可区域的有限候选。先核对翻译、RNA 保护区域、intron/边界、接口与位点，再按文献证据、最小改动和模型辅助比较；不能用模型分数弥补约束失败。

HPC 不可达、路径未给、环境未核验、排队、超时、OOM、缺行/重复行/NaN 都明确记为未执行或失败。`run_receipt.json` 不存在或执行失败时 `collect` 拒绝把文件当实际运行；外部执行的结果需另外登记真实运行证据，不能伪造 receipt。

候选最好也可能仍是原版；模型生成不是自动“优化成功”。任何插入/删除、天然调控元件替换或 RNA 功能改写需要新注释及独立实验验证。
