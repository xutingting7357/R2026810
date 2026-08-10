# R2026810
## 项目简介
本仓库为基于TCGA、GEO公共数据集的肿瘤预后模型完整复现代码，包含数据预处理、差异基因筛选、功能富集、预后模型构建、免疫浸润、临床相关性、药物分子网络全套分析流程，全部基于R语言实现。

## 目录结构说明
1-data                     # 输入数据目录（表达矩阵、临床基线信息、样本分组注释）
2-DRGs                     # 差异表达基因筛选代码
3-Diffanalysis             # 差异分析可视化（火山图、热图）
4-GOKEGG                   # GO/KEGG 功能富集分析
5-PrognosticModel          # 预后风险评分模型构建
6-PrognosticAnalysis-TCGA  # TCGA 队列预后生存、ROC、单 / 多因素 Cox 分析
7-PrognosticAnalysis-GEO   # GEO 外部验证队列预后验证
9-clinical correlation     # 风险评分与临床指标关联分析
10-GSEA-risk               # 基于风险分组的 GSEA 富集分析
11-CIBERSORT               # 肿瘤微环境免疫浸润分析
12-drug-molecular network  # 风险基因药物预测与分子网络构建

## 数据说明
仓库仅上传分析脚本，原始表达矩阵（count/FPKM CSV）、完整临床表格因文件体积过大未存入仓库。
数据集压缩包下载地址：（后续填入Zenodo/云盘链接）
使用说明：下载后解压放入 `1-data` 文件夹即可匹配代码路径。

## 运行环境
- R ≥ 4.2.0
- 核心依赖包：TCGAbiolinks, limma, DESeq2, clusterProfiler, survival, glmnet, pheatmap, ggplot2, CIBERSORT, igraph

## 使用步骤
1. 下载本仓库全部代码
2. 获取配套数据集并放置至 `1-data`
3. 按文件夹数字顺序依次运行R脚本完成全套分析

## 许可协议
MIT License
