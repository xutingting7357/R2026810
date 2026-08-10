# 1 环境变量准备----
## 1.1 清空环境----
rm(list=ls())
gc()

## 1.2 R包安装----
# 定义一个函数来检查并安装所需的R包
install_if_missing <- function(packages) {
  installed_packages <- rownames(installed.packages())
  for (pkg in packages) {
    if (!(pkg %in% installed_packages)) {
      BiocManager::install(pkg, dependencies = TRUE, update = FALSE)
    }
  }
}

# 检查安装的包列表
required_packages <- c("magrittr","DESeq2","dplyr","ggplot2","ggrepel","ggvenn",
                       "pheatmap","RCircos")

# 调用函数进行检查和安装
install_if_missing(required_packages)

## 1.3 加载R包----
library(magrittr)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggvenn)
library(pheatmap)
library(RCircos)

## 1.4 创建输出文件夹----
dir.create("input/")
dir.create("output/") 

## 1.5 定义变量----
pheno_abbr <- "D"    # 表型缩写
treat_name = "STAD"    # 疾病组命名
con_name = "Control"  # 对照组命名

# 2 读取数据----
RGs <- read.csv("input/DRGs.csv")

## 排序样本，tumor在前
TCGA_group <- data.table::fread("input/PD.csv", data.table = F)[, c(1, 2)]%>% 
  dplyr::arrange(Group)
colnames(TCGA_group)[1] <- "RNAseq样本编号"
head(TCGA_group)
tail(TCGA_group)

## 矩阵样本与TCGA_group样本顺序一致
TCGA_mat <- data.table::fread("input/Matrix_Counts.csv", data.table = F) %>% 
  tibble::column_to_rownames("V1") %>% 
  dplyr::select(TCGA_group$RNAseq样本编号) 

# TCGA_mat <- TCGA_mat[,TCGA_group$RNAseq样本编号]

identical(TCGA_group$RNAseq样本编号,colnames(TCGA_mat)) ## TRUE顺序一致

# 3 差异分析----
group_list <- factor(TCGA_group$Group, levels = c(treat_name,con_name))
colData <- data.frame(row.names = colnames(TCGA_mat), group_list = group_list)
dds <- DESeqDataSetFromMatrix(countData = round(TCGA_mat),
                              colData = colData,
                              design = ~ group_list)
dds2 <- DESeq(dds)
res <-  results(dds2, contrast = c("group_list",treat_name,con_name))

# 转为数据框并以padj排序
res1 <- res %>% data.frame() %>% dplyr::arrange(padj) 

length(which((abs(res1$log2FoldChange) > 4) & (res1$padj < 0.05)))
length(which((abs(res1$log2FoldChange) > 2) & (res1$padj < 0.05)))
length(which((abs(res1$log2FoldChange) > 1) & (res1$padj < 0.05)))
length(which((abs(res1$log2FoldChange) > 0.5) & (res1$padj < 0.05)))
length(which((abs(res1$log2FoldChange) > 0) & (res1$padj < 0.05)))

res1 <- na.omit(res1)
res1$gene_symbol <- row.names(res1)
# 修改：log2FoldChange阈值，padj
logFC_cutoff <- 0
res.up <- res1 %>% dplyr::filter(log2FoldChange > logFC_cutoff & padj < 0.05)
res.down <- res1 %>% dplyr::filter(log2FoldChange < -logFC_cutoff & padj < 0.05)

# 4. 输出结果----
res_final <- res1 %>% mutate(group = case_when(
  gene_symbol %in% res.up$gene_symbol ~ "Up",
  gene_symbol %in% res.down$gene_symbol ~ "Down",
  TRUE ~ "Not"))

DEGs <- res_final[res_final$group == "Up"|res_final$group == "Down",]
RDEGs <- intersect(rownames(DEGs),RGs$x)

write.csv(res_final, file = paste0("output/1-diffAnalysis_logFC=",logFC_cutoff,
                                   "_padj=0.05.csv"))
write.csv(DEGs$gene_symbol,"output/2-DEGs.csv",row.names = F)
write.csv(RDEGs,"output/3-RDEGs.csv",row.names = F)


# 5 结果可视化----
## 5.1 火山图----
# 修改：P.Value或者adj.p与上面保持一致
dat_volcanoplot <- res_final[, c('log2FoldChange', 'padj', 'group')]
dat_volcanoplot$logP <- -log10(dat_volcanoplot$padj)
dat_volcanoplot$change <- factor(dat_volcanoplot$group,
                                 levels = c('Down', 'Not', 'Up'))


data_label <- dat_volcanoplot[rownames(dat_volcanoplot) %in% RDEGs,]
data_label$gene_symbol <- rownames(data_label)

# 火山图添加基因标签
# p_volcano <- ggplot(dat_volcanoplot, aes(log2FoldChange, logP, color = change)) +
#   geom_point(alpha = 0.6) +
#   theme_bw() +
#   # 修改：P.Value或者adj.p
#   labs(x = 'log2FoldChange', y = '-Log10(padj)', color = 'Significance') +
#   geom_hline(yintercept = -log10(0.05), lty = 2) +
#   geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 2) +
#   scale_x_continuous(limits = c(-5, 5)) +
#   scale_y_continuous(limits = c(0, 240)) +
#   scale_color_manual(values = c('#4DBBD5', 'grey', '#E64B35'))+
#   ggtitle("TCGA-BRCA")+
#   theme(plot.title = element_text(hjust = 0.5, size = 16))+
#   geom_point(size= 1,shape=1,color="black",data=data_label)+
#   geom_text_repel(data = data_label, aes(label = gene_symbol),
#                   color="black",size = 3,
#                   cex = 2,max.overlaps = 30,min.segment.length = 0)
# p_volcano

p_volcano <- ggplot(dat_volcanoplot, aes(log2FoldChange, logP, color = change)) +
  geom_point(alpha = 0.6) +
  theme_bw() +
  # 修改：P.Value或者adj.p
  labs(x = 'log2FoldChange', y = '-Log10(padj)', color = 'Significance') +
  geom_hline(yintercept = -log10(0.05), lty = 2) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 2) +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(0, 75)) +
  scale_color_manual(values = c('#4DBBD5', 'grey', '#E64B35'))+
  ggtitle("TCGA-STAD")+
  theme(plot.title = element_text(hjust = 0.5, size = 16)) +
  theme(panel.grid = element_blank())  # 添加这行代码去除所有网格线
p_volcano
ggsave(file = 'output/4-vocano.pdf', p_volcano, width = 10, height = 7)


## 5.2 韦恩图----
p_venn <- ggvenn(
  setNames(list(DEGs$gene_symbol, RGs$x), c("DEGs", paste0(pheno_abbr, "RGs"))),
  # 下面要与列表中的命名一致
  c('DEGs', paste0(pheno_abbr, "RGs")),
  # 不展示比例
  show_percentage = F,
  fill_alpha = 0.5,
  stroke_color = NA,
  fill_color = c('#a5dff9', '#FDD692'))
p_venn
ggsave(file = 'output/5_venn.pdf', p_venn, width = 5, height = 5)


## 5.3 热图----
# 利用FPKM中的表达谱绘制热图
## 画图Normal组在前
TCGA_group <- data.table::fread("input/PD.csv", data.table = F)[, c(1, 2)]%>% 
  dplyr::arrange(desc(Group))
head(TCGA_group)
colnames(TCGA_group)[1] <- "RNAseq样本编号"
TCGA_mat_FPKM <- data.table::fread("input/Matrix_FPKM.csv", data.table = F) %>% 
  tibble::column_to_rownames("V1") %>% 
  dplyr::select(TCGA_group$RNAseq样本编号) 

identical(TCGA_group$RNAseq样本编号,colnames(TCGA_mat_FPKM)) ## TRUE顺序一致


dat_RDEGs <- DEGs[RDEGs,]

# 当表型相关差异基因较多时选择logFC_top20展示热图
if (length(RDEGs)>20){
  RDEGs_sorted <- dat_RDEGs[order(dat_RDEGs$log2FoldChange), ]
  top_20 <- RDEGs_sorted[c(1:10,(nrow(RDEGs_sorted) - 9):nrow(RDEGs_sorted)), ]
  DEG_heatmap <- top_20$gene_symbol
}else{
  DEG_heatmap <- RDEGs
}
dat_heatmap <- TCGA_mat_FPKM[DEG_heatmap,]

p_heatmap <- pheatmap(dat_heatmap,
                      # 对行归一化（每个基因在样本中的表达量）
                      scale = 'row',
                      # 列注释
                      annotation_col = data.frame(Group = TCGA_group$Group, 
                                          row.names = TCGA_group$RNAseq样本编号),
                      # 列注释颜色
                      annotation_colors = list(Group = rlang::set_names(
                        c('#a5dff9', '#FDD692'), c(con_name, treat_name))),
                      # 热图颜色
                      color = colorRampPalette(
                        c('#4DBBD5', '#4DBBD5', 'white', '#E64B35', '#E64B35'))(50),
                      # 颜色范围（注意这里的50要和上面括号中的50保持一致）
                      breaks = c(seq(-3, 3, length = 50)),
                      # 行聚类和列聚类
                      cluster_cols = F,
                      cluster_rows = T,
                      labels_col = '',
                      border_color = NA,
                      main = "TCGA-STAD")
p_heatmap
ggsave(file = 'output/6-HeatMap.pdf', p_heatmap, width = 8, height = 6)

# ## 5.4 染色体定位图----
# chr <- read.csv("input/chromosome_locate_input.csv") # 输入文件，不用改
# chr_gene <- chr[which(chr$Gene %in% RDEGs),]
# 
# pdf(width = 8.1,height = 8,file = "output/7-Chromosome_Localization.pdf")
# data(UCSC.HG38.Human.CytoBandIdeogram)
# cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
# RCircos.Set.Core.Components(cyto.info)
# RCircos.Set.Plot.Area()
# RCircos.Chromosome.Ideogram.Plot()
# RCircos.Gene.Connector.Plot(chr_gene, track.num = 1, side = "in")
# RCircos.Gene.Name.Plot(chr_gene, name.col=4,track.num=2, side="in")
# dev.off()

## 5.5 输出GSEA富集分析输入文件----
write.csv(res_final[,2,drop = F],"output/8-GSEA_input.csv")

## 5.6 保存工作环境----
save.image("3-DiffAnalysis.RData.gz",compress = "gzip")




