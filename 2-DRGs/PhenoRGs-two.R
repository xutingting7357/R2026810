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
required_packages <- c("dplyr", "openxlsx", "rvest","msig","stringr")

# 调用函数进行检查和安装
install_if_missing(required_packages)

## 1.3 加载R包----
library(dplyr)
library(openxlsx)
library(rvest)
library(msig)
library(stringr)

## 1.4 创建输出文件夹----
dir.create("input/")
dir.create("output/")

## 1.5 定义变量----
# pheno1 <- "Ferroptosis"
# pheno1_file1 <- "input/ferroptosis_driver.csv"
# pheno1_file2 <- "input/Ferroptosis-铁死亡.xlsx"
# pheno1_file3 <- "input/ferroptosis_marker.csv"
# pheno1_file4 <- "input/ferroptosis_suppressor.csv"
# pheno1_file5 <- "input/ferroptosis_unclassified.csv"

pheno2 <- "deubiquitination"
pheno2_file1 <- "input/GeneCards-SearchResults.csv"
pheno2_file2 <- "input/去泛素化.xlsx"

pheno_abbr <- "D"

# 2 表型1合并、去重----

## 2.1 GeneCards来源----
# pheno1_RGs_genecards <- data.table::fread(pheno1_file1,data.table = F)
# # 先读取第一个文件
# pheno1_RGs_genecards <- data.table::fread(pheno1_file1, data.table = FALSE)
# 
# # 逐个合并其他文件
# pheno1_RGs_genecards <- rbind(pheno1_RGs_genecards, data.table::fread(pheno1_file3, data.table = FALSE))
# pheno1_RGs_genecards <- rbind(pheno1_RGs_genecards, data.table::fread(pheno1_file4, data.table = FALSE))
# pheno1_RGs_genecards <- rbind(pheno1_RGs_genecards, data.table::fread(pheno1_file5, data.table = FALSE))
# ## 选择蛋白编码基因并依据相关性分数筛选基因（可修改）
# # pheno1_RGs_genecards_filter <- pheno1_RGs_genecards %>%
# #   dplyr::filter(uniformgenetype == 'gene with protein product' )
# pheno1_RGs_genecards_filter <- pheno1_RGs_genecards#基因太少保留所有基因
# ## 统一命名基因列
# pheno1_RGs_genecards_filter <- pheno1_RGs_genecards_filter %>%
#   dplyr::select("symbol") %>%
#   dplyr::rename("gene" = "symbol")
# 
# ###铁死亡数据库来源基因去重输出564个基因
# pheno1_data_unique1 <- pheno1_RGs_genecards %>%
#   distinct(symbol, .keep_all = TRUE)

## 2.2 文献来源----
# 定义函数：读取 Excel 文件中指定 sheet 的第一列(视情况修改列数)数据并合并
Merge_sheets_onecol_F <- function(file_path, selected_sheets) {
  all_data <- lapply(selected_sheets, function(sheet) {
    data <- readxl::read_excel(file_path, sheet = sheet,col_names = T)
    data[, 1]
  })
  merged_data <- unlist(all_data)
  return(merged_data)
}

# 指定要合并的 sheet 名称（需改）
# selected_sheets <- c("T01","T02","T03","T04","T05","T06","T07","T08","T09")
# pheno1_RGs_paper <- data.frame(Merge_sheets_onecol_F(pheno1_file2,
#                                                      selected_sheets))
# colnames(pheno1_RGs_paper)[1] <- "gene"
# 
# ###铁死亡文献来源基因去重输出601个基因
# pheno1_data_unique2 <- pheno1_RGs_paper %>%
#   distinct(gene, .keep_all = TRUE)
# 
# ## 2.3 表型基因的合并、去重（取并集）----
# pheno1_RGs <- union(pheno1_RGs_genecards_filter$gene, pheno1_RGs_paper$gene) %>%
#   as.data.frame()
# colnames(pheno1_RGs) <- "symbol"


# 3 表型2合并、去重----

## 3.1 GeneCards来源----
pheno2_RGs_genecards <- data.table::fread(pheno2_file1,data.table = F)

## 选择蛋白编码基因并依据相关性分数筛选基因（可修改）
pheno2_RGs_genecards_filter <- pheno2_RGs_genecards %>% 
  dplyr::filter(Category == 'Protein Coding' &`Relevance score` > 1)

## 重新命名基因列
pheno2_RGs_genecards_filter <- pheno2_RGs_genecards_filter %>% 
  dplyr::select("Gene Symbol") %>% 
  dplyr::rename("gene" = "Gene Symbol")


## 3.2 文献来源----
# 读取 Excel 文件中指定 sheet 的第一列数据并合并（需改）
selected_sheets <- c("T01","T02","T03")
pheno2_RGs_paper <- data.frame(Merge_sheets_onecol_F(pheno2_file2, 
                                                     selected_sheets))
colnames(pheno2_RGs_paper)[1] <- "gene"


## 3.3 表型基因的合并、去重（取并集）----
pheno2_RGs <- union(pheno2_RGs_genecards_filter$gene, pheno2_RGs_paper$gene) %>% 
  as.data.frame()
colnames(pheno2_RGs) <- "symbol"


# 4 两个表型基因取交集----
# RGs <- intersect(pheno1_RGs, pheno2_RGs)
RGs <- pheno2_RGs
# 表型基因和所有GEO和TCGA数据集基因再取交集
GEO_genes_1 <- read.csv("input/GSE84437_Datasets_Matrix.csv")[,1]
GEO_genes_2 <- read.csv("input/Matrix_FPKM.csv")[,1]
GEO_genes <- intersect(GEO_genes_1, GEO_genes_2)

pheno_RGs <- intersect(RGs$symbol, GEO_genes)
length(pheno_RGs)

# 5 输出结果保存为csv文件----
write.csv(pheno_RGs, file = paste0("output/",pheno_abbr,"RGs.csv"), row.names = F)

