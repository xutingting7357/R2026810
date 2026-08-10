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
required_packages <- c("dplyr")

# 调用函数进行检查和安装
install_if_missing(required_packages)

## 1.3 加载R包----
library(dplyr)

## 1.4 创建文件夹----
dir.create("input/")
dir.create("output/")

## 1.5 定义变量----
tumor_name <- "STAD"

# 2 数据清洗----
PD <- data.table::fread("output/PD.csv",data.table = F)

## 2.1 提取蛋白编码基因----
genetype <- data.table::fread("input/id_trans.csv", data.table = F) %>%
  dplyr::select(1, 2, 4)
gene_mrna <- genetype %>%
  dplyr::filter(gene_type == "protein_coding")


## 2.2 FPKM数据清洗----
dat_fpkm <- data.table::fread("input/TCGA-STAD.txt",
                                     check.names = F,data.table = F)

# 只保留蛋白编码基因
dat_fpkm <- dat_fpkm[dat_fpkm$gene_id %in% gene_mrna$gene_name,] 

# 如果有多个基因取表达平均值
dat_fpkm <- aggregate(dat_fpkm, by = dat_fpkm$gene_id %>% list(),FUN = mean)
rownames(dat_fpkm) <- dat_fpkm$Group.1
dat_fpkm <- dat_fpkm[,-(1:2)]


# 范围过大时进行log处理
range(dat_fpkm)
dat_fpkm <- log2(dat_fpkm+1)
dat_fpkm <- dat_fpkm[,PD$ID]

## 2.3 COUNTS数据清洗----
dat_counts <- data.table::fread("input/TCGA-STAD-COUNTS.txt",
                                check.names = F,data.table = F)

#只保留蛋白编码基因
dat_counts <- dat_counts[dat_counts$gene_id %in% gene_mrna$gene_name,]
# 如果有多个基因取表达平均值
dat_counts <- aggregate(dat_counts, by = dat_counts$gene_id %>% list(),FUN = mean)
rownames(dat_counts) <- dat_counts$Group.1
dat_counts <- dat_counts[,-(1:2)]

range(dat_counts)
dat_counts <- dat_counts[,PD$ID]


#去除0方差的基因
dat_counts <- dat_counts[apply(dat_counts,1,var) != 0,] 
#去除低表达量的基因
dat_counts <- dat_counts[rowSums(dat_counts) > 1,] 

# fpkm基因与count基因对齐
dat_fpkm <- dat_fpkm[rownames(dat_counts),] 
write.csv(dat_fpkm,"output/Matrix_FPKM.csv")
write.csv(dat_counts,"output/Matrix_Counts.csv")
write.csv(dat_counts[,PD$Group == "STAD"],"output/Disease_Matrix_Counts.csv")

#从fpkm基因中提取疾病组
group <- data.table::fread("output/sample_group.csv",data.table = F)
group <- group[match(colnames(dat_fpkm), group$ID),]
all(group$ID == colnames(dat_fpkm))
write.csv(dat_fpkm[,group$status == tumor_name],
          "output/Disease_Matrix_FPKM.csv")

