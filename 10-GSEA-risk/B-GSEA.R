# 1 环境变量准备----
## 1.1 清空环境----
rm(list=ls())
options(stringsAsFactors = F)
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
required_packages <- c("clusterProfiler", "org.Hs.eg.db","tidyverse","ggplot2",
                       'enrichplot','ggridges','GseaVis')

#调用函数进行检查和安装
install_if_missing(required_packages)

## 1.3 加载R包----
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ggplot2)
library(enrichplot)
library(ggridges)
library(GseaVis)

## 1.4 创建文件夹----
dir.create("input/")
dir.create("output/")

# 2 基因ID转换与排序----
GSEA_input <- read.csv("output/5-GSEA_input.csv")
colnames(GSEA_input)[1] = "gene"
GSEA_input <- GSEA_input[!duplicated(GSEA_input$gene),]

# 添加ENTREZID列
log2FoldChange_list <- bitr(GSEA_input$gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
dat_GSEA <- inner_join(GSEA_input,log2FoldChange_list,by=c("gene"="SYMBOL")) 

# 按照log2FoldChange值对基因进行排序
dat_GSEA =dat_GSEA[order(dat_GSEA$log2FoldChange,decreasing = T),]
log2FoldChange_list <-  dat_GSEA$log2FoldChange
names(log2FoldChange_list) = as.character(dat_GSEA[,'gene'])
head(log2FoldChange_list)

# 3 MSigDb的GSEA富集分析----
geneset_c2 <- read.gmt('input/c2.all.v2024.1.Hs.symbols.gmt')

# 进行GSEA富集分析并设置种子，重复结果
GSEA<-GSEA(log2FoldChange_list,seed=2022,
           exponent = 1,
           minGSSize = 10,
           maxGSSize = 500,
           eps = 1e-10, # 设置计算 p 值的边界
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           TERM2GENE = geneset_c2) 

GSEA_result <- as.data.frame(GSEA)
colnames(GSEA_result)
max(GSEA_result$qvalue)
GSEA_result <- GSEA_result[GSEA_result$qvalue < 0.05,]

# 根据NES从高到底对结果进行排序
GSEA_result <- GSEA_result[order(GSEA_result$NES,decreasing = T),]
write.csv(GSEA_result,file="output/GSEA_MSigDb_result.csv",
          row.names = T,quote = T)

# 筛选所需通路后重新保存并读取，可利用GSEA_筛选明星通路脚本进行筛选并保存
GSEA_result_filter <- read.csv("output/GSEA_term_selected.csv")
GSEA_result_filter <- as.data.frame(GSEA_result_filter)
GSEA_result_filter <- na.omit(GSEA_result_filter)
GSEA@result <- GSEA@result[GSEA@result$ID %in% GSEA_result_filter$ID,]

# 4 可视化----
## 4.1 山峦图 ----
y_labels <- GSEA_result_filter$ID
y_labels <- gsub("_"," ",y_labels)
y_labels <- str_to_title(y_labels)
y_labels <- str_wrap(y_labels, width = 40)

p_ridge <- ridgeplot(GSEA,showCategory=length(GSEA_result_filter$ID),
                     fill = "p.adjust")+
  scale_fill_gradient(
    high = '#4DBBD5',
    low = '#E64B35',
    #limits = c(0, 0.05)
  ) +
  scale_y_discrete(labels = y_labels) 
ggsave(file = "output/GSEA_ridge.pdf", p_ridge, width=10, height =5)


## 4.2 经典通路图 ----
for(i in 1:length(GSEA_result_filter$ID)) { 
  title  <- GSEA@result$Description[i]
  p<-gseaNb(object = GSEA,
            geneSetID = title,
            addPval = T,
            htCol=c('#4DBBD5', '#E64B35'),
            geneSize=4,
            # curveCol=c('#E64B35','#4DBBD5'),
            rankCol = c('#4DBBD5',"white", '#E64B35'),
            rmSegment = T,
            termWidth = 100,
            pvalX = 0.75,pvalY = 0.75,
            pvalSize=4,
            pDigit=3,
            nesDigit = 3,
            pCol = 'black',
            rankSeq =10000,
            #Selectadd = val,
            pHjust = 0,
            rmPrefix = FALSE)
  pdf(file= paste0("output/",i,"-",title,".pdf"),width=5, height=6)
  print(p)
  dev.off()
}

