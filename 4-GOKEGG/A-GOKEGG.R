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
required_packages <- c("clusterProfiler","org.Hs.eg.db","ggplot2","enrichplot",
                       "RColorBrewer","ggrepel","GOplot",'rlang',"stringr")

# 调用函数进行检查和安装
install_if_missing(required_packages)

## 1.3 加载R包----
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(RColorBrewer)
library(ggrepel)
library(GOplot)
library(rlang)
library(stringr)

## 1.4 创建文件夹----
dir.create("input/")
dir.create("output/")

# 2 GO,KEGG富集分析 ----
DEGs <- read.csv("input/3-RDEGs.csv",header = T)
dat_DEGs <- read.csv("input/1-diffAnalysis_logFC=0_padj=0.05.csv")
rownames(dat_DEGs) <- dat_DEGs$X
dat_DEGs <- dat_DEGs[DEGs$x,-1]
colnames(dat_DEGs)[2] <- "logFC"

## 2.1 基因名-ENTREZID转换----
genelist <- bitr(DEGs$x, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

## 2.2 GO富集分析----
raw_RGs_GO <- enrichGO(gene = genelist$ENTREZID, OrgDb ='org.Hs.eg.db',
                       pAdjustMethod="BH",
                       # BP, CC, MF
                       ont = "ALL",
                       pvalueCutoff = 0.05,qvalueCutoff = 0.25,
                       # 基因ID转换为基因名
                       readable = T)

# 提取GO分析结果，按p.adj升序排序，BP、CC、MF各取top5的条目
RGs_GO <- as.data.frame(raw_RGs_GO)
RGs_GO <- RGs_GO[order(RGs_GO$p.adjust,decreasing = F),]
RGs_GO_5 <- rbind(RGs_GO %>% dplyr::filter(ONTOLOGY == 'BP') %>% .[1:5, ],
                   RGs_GO %>% dplyr::filter(ONTOLOGY == 'CC') %>% .[1:5, ],
                   RGs_GO %>% dplyr::filter(ONTOLOGY == 'MF') %>% .[1:5, ]
                   ) %>% na.omit()

## 2.2 KEGG分析----
# 设置clusterProfiler包内部的一些标准
R.utils::setOption( "clusterProfiler.download.method",'wininet')

raw_RGs_KEGG <- enrichKEGG(gene = genelist$ENTREZID, organism ="hsa", 
                 keyType = "kegg", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff =0.25) %>%
                 setReadable(OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")

# 提取KEGG分析结果
RGs_KEGG <- as.data.frame(raw_RGs_KEGG)
RGs_KEGG <- RGs_KEGG[order(RGs_KEGG$p.adjust,decreasing = F),]
RGs_KEGG_5 <- RGs_KEGG[1:5, ] %>% na.omit() 
RGs_KEGG_5 <- RGs_KEGG_5[,-c(1:2)]
RGs_KEGG_5$ONTOLOGY <- "KEGG"
RGs_KEGG_5 <- RGs_KEGG_5[, c("ONTOLOGY", setdiff(names(RGs_KEGG_5), "ONTOLOGY"))]
identical(colnames(RGs_KEGG_5),colnames(RGs_GO_5))

GO_KEGG <- rbind(RGs_GO_5,RGs_KEGG_5)

## 2.3 保存GO,KEGG分析结果----
write.csv(GO_KEGG,file="output/1-GO_KEGG_result.csv",row.names = T,quote = F)


# 3 GO,KEGG可视化----
# 合并GO,KEGG的富集分析结果
# 将GeneRatio列转为小数，便于绘图
GO_KEGG$GeneRatio <- sapply(GO_KEGG$GeneRatio, 
                             function(x) as.numeric(gsub("/.*", "", x)) / 
                               as.numeric(gsub(".*/", "", x)))

# 保留top5的富集结果
raw_RGs_GO@result <- raw_RGs_GO@result[raw_RGs_GO@result$ID %in% GO_KEGG$ID,]

# 确保 ONTOLOGY 列是因子类型，并指定水平顺序
GO_KEGG$ONTOLOGY <- factor(GO_KEGG$ONTOLOGY, levels = c("BP", "CC", "MF", "KEGG"))

# 柱状图
ggplot(data = GO_KEGG,
       aes(x = Description,y = -log10(p.adjust),fill = ONTOLOGY))+
  geom_bar(stat = "identity")+
  #根据ONTOLOGY进行分面
  facet_wrap(~ ONTOLOGY, scales = "free_x",ncol = 4)+
  scale_fill_manual(values = c("#FFC640","#9BC7B2","#FFFFB3","#BEBADA"))+
  #翻转坐标轴
  #coord_flip()+
  xlab("Term")+
  theme_bw()+
  theme(axis.text.x=element_text(face = "bold",
                                 angle = 40,vjust = 1, hjust = 1 ))
ggsave("output/2-GO_KEGG_bar.pdf",width = 16,height = 7)


# 气泡图
# 设置换行宽度
GO_KEGG$Description_wrapped <- str_wrap(GO_KEGG$Description, width = 30)

ggplot(data = GO_KEGG, aes(x = Description_wrapped, 
                           y = GeneRatio, colour = p.adjust)) +
  geom_point(aes(size = Count), shape = 16, stat = "identity") +
  scale_color_gradient(high = "#4DBBD5", low = "#E64B35") +
  facet_wrap(~ ONTOLOGY, scales = "free_x", ncol = 4) +
  scale_fill_manual(values = c("#FFC640", "#9BC7B2", "#FFFFB3", "#BEBADA")) +
  #翻转坐标轴
  #coord_flip()+
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold",angle = 30,vjust = 0.5, 
                                   hjust = 1, size = 10,lineheight = 0.9),
        panel.spacing = unit(1, "lines"))

ggsave("output/3-GO_KEGG_bubble.pdf",width = 16,height = 5)

# GO-BP网络图
pdf(file="output/4-GO_net-BP.pdf", width=30, height=24)
cnetplot(raw_RGs_GO,
         categorySize="pvalue", 
         showCategory = GO_KEGG$Description[which(GO_KEGG$ONTOLOGY == "BP")],
         # 将基因名与logFC对应起来
         foldChange = set_names(dat_DEGs[genelist$SYMBOL, 'logFC'], 
                                genelist$SYMBOL),
         color.params = list(category = '#3C5488'),
         # 调整标签字体大小
         cex_label_category = 2,    # 调整GO术语标签大小
         cex_label_gene = 2,        # 调整基因标签大小
         # 画成环状
         # circular = TRUE,
         # 展示图中全部信息
         node_label = "all",
         colorEdge = F)+
  scale_color_gradient(high = "#E64B35",low = "#4DBBD5", 
                       limits = c(-1, 1)) +
  labs(color = 'LogFC', size = 'Count') 
dev.off()

# GO-CC网络图
pdf(file="output/5-GO_net-CC.pdf", width=30, height=24)
cnetplot(raw_RGs_GO, 
         categorySize="pvalue", 
         showCategory = GO_KEGG$Description[which(GO_KEGG$ONTOLOGY == "CC")],
         foldChange = set_names(dat_DEGs[genelist$SYMBOL, 'logFC'], 
                                genelist$SYMBOL),
         color.params = list(category = '#3C5488'),
         # 调整标签字体大小
         cex_label_category = 2,    # 调整GO术语标签大小
         cex_label_gene = 2,        # 调整基因标签大小
         #circular = TRUE,
         node_label = "all",
         colorEdge = F)+
  scale_color_gradient(high = "#E64B35",low = "#4DBBD5", 
                       limits = c(-1, 1)) +
  labs(color = 'LogFC', size = 'Count')
dev.off()

# MF网络图
pdf(file="output/6-GO_net-MF.pdf", width=30, height=24)
cnetplot(raw_RGs_GO,
         categorySize="pvalue", 
         showCategory = GO_KEGG$Description[which(GO_KEGG$ONTOLOGY == "MF")],
         foldChange = set_names(dat_DEGs[genelist$SYMBOL, 'logFC'], 
                                genelist$SYMBOL),
         color.params = list(category = '#3C5488'),
         # 调整标签字体大小
         cex_label_category = 2,    # 调整GO术语标签大小
         cex_label_gene = 2,        # 调整基因标签大小
         #circular = TRUE,
         node_label = "all",
         colorEdge = F)+
  scale_color_gradient(high = "#E64B35",low = "#4DBBD5", 
                       limits = c(-2, 2)) +
  labs(color = 'LogFC', size = 'Count')

dev.off()

# KEGG网络图
pdf(file="output/7-GO_net-KEGG.pdf", width=30, height=24)
cnetplot(raw_RGs_KEGG,
         categorySize="pvalue", 
         showCategory = GO_KEGG$Description[which(GO_KEGG$ONTOLOGY == "KEGG")],
         foldChange = set_names(dat_DEGs[genelist$SYMBOL, 'logFC'], 
                                genelist$SYMBOL),
         color.params = list(category = '#3C5488'),
         # 调整标签字体大小
         cex_label_category = 2,    # 调整GO术语标签大小
         cex_label_gene = 2,        # 调整基因标签大小
         #circular = TRUE ,
         node_label = "all",
         colorEdge = F)+
  scale_color_gradient(high = "#E64B35",low = "#4DBBD5", 
                       limits = c(-1.5, 1.5)) +
  labs(color = 'LogFC', size = 'Count')

dev.off()

# 准备bubble图输入数据
dat_bubble <- data.frame(Category=GO_KEGG$ONTOLOGY, ID=GO_KEGG$ID,
                 Term=GO_KEGG$Description, 
                 Genes = gsub("/", ", ", GO_KEGG$geneID), 
                 adj_pval = GO_KEGG$p.adjust)

symbol_logFC <- data.frame(ID = dat_DEGs$gene_symbol, logFC = dat_DEGs$logFC)
row.names(symbol_logFC)=symbol_logFC[,1]
circ <- circle_dat(dat_bubble, symbol_logFC)

pdf(file="output/8-GO_bubble.pdf", width=8, height=6)
# 去除重复的通路名
circ1 <- circ[!duplicated(circ$ID),]
circ1levels <- as.character(circ1$category[which(!duplicated(circ1$category))])

circ1$category <- factor(circ1$category,levels= circ1levels)
ggplot(circ1,aes(x=zscore,y=-log10(adj_pval)))+
  geom_point(aes(size=count,color=category),alpha=0.8)+
  scale_color_brewer(palette = "Accent")+
  theme_bw()+
  theme(
    legend.position = c("none"),
    # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )+
  geom_text_repel(
    data = circ1[-log10(circ1$adj_pval)>-2,],
    # 添加符合标准的ID
    aes(label = circ1[-log10(circ1$adj_pval)>-2,]$ID),
    size = 3, max.overlaps = 20,
    segment.color = "black", show.legend = FALSE )+
  # 根据通路类别进行分面
  facet_grid(.~category)
dev.off()

# EMAP图
# 获取相似性评分矩阵
dat_emap <- pairwise_termsim(raw_RGs_GO)

set.seed(2024)
p <- emapplot(dat_emap,showCategory = 20,
  # 标签不重叠
  repel = F,
  # 布局方式
  # 可选：'star','circle','gem','dh','graphopt','grid','mds',
  # 'randomly','fr','kk','drl','lgl'
  layout.params = list(layout = 'fr')) +
  labs(fill = 'pvalue') +
  scale_fill_gradient(high = "#4DBBD5",low = "#E64B35")

ggsave(file = 'output/9_emap_network.pdf', p, width = 9, height = 8)

