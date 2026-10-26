# 文献、附件与序列追溯

## 检索顺序

1. 从全文、Methods、图注、补充表/文件和数据可用性声明列出构建与序列；保留原始文件、URL、文件名、表号/页码和时间。
2. 优先原构建的完整质粒/表达框文件，交叉核验原文明确引用的质粒库编号、GenBank accession 和作者仓库。来源冲突先比较版本，不按优先顺序盲目覆盖。
3. 提取附件中的序列、基因号、引物及边界。区分完整构建、CDS、引物、接头和示意图；OCR 结果核对原图，字符歧义阻断精确序列输出。
4. 依据 Gene ID/基因座号/accession 主动检索 NCBI；结合物种/TaxID、转录本和文献材料消歧，不静默选最长转录本或同名基因。
5. 追查元件原始论文和附件；仍未找到时记录已尝试来源。无法访问不等于没有；作者邮件可起草但不自动发送。
6. 只对有记录的缺口提出替代，保留已复刻部分。仅有蛋白号不能反向翻译后称原始 DNA；天然 CDS 不能覆盖作者未公开的优化 CDS。

NCBI Datasets 支持 Gene ID、带物种的 symbol 和 RefSeq accession，并可获取 gene/RNA/CDS/protein 与元数据（[官方说明](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genes/download-gene-data-package/)）。

有官方 `datasets` CLI 时，可用以下已文档化的形式获取数据包；占位变量由项目的真实编号填写，不把普通 GenBank 编号误作受支持的 RefSeq 输入（[官方命令示例](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genes/download-gene-data-package/)）。

```bash
datasets download gene gene-id "$GENE_ID" \
  --include gene,rna,cds,protein --filename ncbi_gene.zip
datasets download gene accession "$REFSEQ_ACCESSION" \
  --include gene,rna,cds,protein --filename ncbi_accession.zip
```

其他记录用 NCBI Nucleotide/Protein 原始记录和可用的官方下载途径。工具不可用时使用获准的检索/下载能力，不用伪 accession 或未实际下载的文件路径代替。

NCBI 的序列与注释修订并不总是同时改变 accession.version；保存原文件、注释、CDS 坐标与获取日期，不能假定当前版本等于作者当时版本（[NCBI RefSeq FAQ](https://www.ncbi.nlm.nih.gov/books/NBK50679/)）。

## 证据与来源登记

每个关键论文记录标题、作者、期刊、年、DOI/URL、物种/系统、实际方法/支持实验、主结果、优势、限制与检索截止日期。只读摘要时明确标记，不能宣称完整重建。

序列记录至少有原标识、标识类型、TaxID、材料/转录本、accession.version、原文件与哈希、坐标/方向、边界依据、对应文献构建、差异和核验方式。直接提供、按模板/引物重建、提出的新设计分开标记。

Addgene 的 full/verified 标签还需查看 NGS、全质粒测序或参考/Sanger 组装等来源；不能把所有完整记录都视为实物全质粒实测（[Addgene 序列说明](https://www.addgene.org/111834/sequences/)）。

## 阻断和复刻等级

- **完整质粒复刻**：有基准全长文件，核对完整序列；环状起点和整体反向表示的等价性须明确。
- **表达框复刻/骨架适配**：表达框一致但骨架变更，不称完整复刻。
- **文献约束重建**：模板、引物、改动与边界可追溯，但无原始全长直接比对。
- **部分复刻＋新设计**：逐元件说明缺口、已检索范围和替代依据。
- **从头设计**：注明未找到适用报道的范围/日期，不宣称全世界不存在。

未知边界、未解决字符、版本冲突、缺少必需序列均阻断最终完整序列交付；可先交付文献矩阵与明确标注的非比例草图。
