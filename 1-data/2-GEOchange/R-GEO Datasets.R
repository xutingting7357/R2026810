# 1 环境变量准备----
## 1.1 清空环境----
rm(list=ls())
gc()
##增加网络延迟时间
options(timeout = 60000)
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
required_packages <- c("GEOquery","dplyr","magrittr","ggplot2","sva","limma",
                       "FactoMineR","factoextra")

# 调用函数进行检查和安装
install_if_missing(required_packages)

## 1.3 加载R包----
library(GEOquery)
library(dplyr)
library(magrittr)
library(ggplot2)
library(sva)
library(limma)
library(FactoMineR)
library(factoextra)

## 1.4 创建文件夹----
dir.create("input/")
dir.create("output/")

## 1.5 定义变量----
GSE_1 <- "GSE84437"   # 修改：GEO数据集
#GSE_2 <- "GSE98566"
# GSE_3 <- "GSE29221"
treat_name = "STAD"    # 修改：疾病组命名
con_name = "Control"  # 对照组命名


# 2 获取数据----
## 2.1 下载原始矩阵数据----
# 本地获取数据
aset_1 <- getGEO(GSE_1,destdir = "input", AnnotGPL = F, getGPL = F)##下载得到原始数据

#GSE98564_raw <-data.table::fread("input/GSE98564_series_matrix.txt",header = T) ##读取原始数据


#aset_2 <- getGEO(GSE_2,destdir = "input", AnnotGPL = F, getGPL = F)##下载得到原始数据

#GSE98566_raw <-data.table::fread("input/GSE98566_series_matrix.txt",header = T) ##读取原始数据
# aset_3 <- getGEO(GSE_3,destdir = "input", AnnotGPL = F, getGPL = F)

## 2.2 定义所需函数----
# 判断是否需要log2处理
log2_if_need <- function(data)
{
  print(c("before:",range(data)))
  ex <- as.data.frame(data)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  ## 如果LogC=T,就进行循环，如果存在小于0的数，就用非数填入
  if (LogC) {
    for(i in 1:ncol(ex)){
      ex[which(ex[,i] < 0),i] <- NaN
    }
    data <- log2(data + 1)
    print("log2 transform finished")
    print(c("after:",range(data)))
  }else{
    print("log2 transform not needed")
  }
  return(data)
}

# 整理分组信息
geo_group <- function(pdata, sample_id, group_id, treat_name, con_name, 
                      treat_key, con_key, desc){
  geogroup <- pdata %>% 
    dplyr::select(ID = sample_id, group0 = group_id) %>%
    dplyr::mutate(group = dplyr::case_when(
      grepl(treat_key, group0) ~ treat_name,
      grepl(con_key, group0) ~ con_name
    )) %>% 
    dplyr::select(1,3) %>% 
    dplyr::filter(group == treat_name | group == con_name) %>% 
    dplyr::arrange(group)
  if(desc){
    geogroup <- geogroup %>% 
      dplyr::arrange(desc(group))
  }
  message("Please check the order of samples!")
  return(geogroup)
}

# 整理芯片平台
gpl_pro <- function(gpl_name, probe_id, gene_id)
{
  gpl_name$gene_id <- strsplit(gpl_name[[gene_id]], split = " /// ") %>%
    lapply(function(x){
      x1 <- x
      if(length(x1) > 1) x1 <- x1[!grepl("^OTTHU", x1)]
      if(length(x1) == 0) x1 <- x[1]
      vec_gene = strsplit(x, split = " // ") %>% lapply(
        function(x) if(length(x) > 1) 
        return(x[2]) else return(x)) %>%
        unlist() %>% unique()
      vec_gene <- vec_gene[!grepl("^OTTHU", vec_gene)]
      vec_gene <- paste0(vec_gene, collapse = " /// ")
      return(vec_gene)
    }) %>% unlist()
  gpl_name <- gpl_name %>% 
    dplyr::select(probe_id, gene_id) %>% # 选择探针名和基因名列
    dplyr::filter(!grepl("///", .[,2])) %>% # 去除基因名中有///的
    dplyr::filter(.[,2] != "") %>% # 去除基因名为空的
    dplyr::filter(!grepl("---", .[,2])) %>% # 去除基因名中有---的
    dplyr::mutate(ID = as.character(ID)) %>% # 探针列转字符
    dplyr::rename(GENE_SYMBOL = gene_id) %>%  # 重命名基因列为GENE_SYMBOL
    dplyr::select(probe_id, GENE_SYMBOL)
  return(gpl_name)
}

# 探针ID转换，同一基因取均值
geo_matrix <- function(data_matrix, data_group, gpl_name){
  data_matrix <- data_matrix[, data_group$ID] %>% # 按group筛选并排序样本
    data.frame() %>% 
    dplyr::mutate(across(!where(is.numeric), as.numeric)) %>% 
    dplyr::mutate(ID = row.names(.) %>% as.character()) %>% 
    dplyr::left_join(., gpl_name, by = "ID") %>% 
    na.omit()
  gene_row <- data_matrix$GENE_SYMBOL
  data_matrix <- data_matrix %>% 
    dplyr::select(-ID, -GENE_SYMBOL) %>% 
    as.matrix()
  rownames(data_matrix) <- gene_row
  data_matrix <- data_matrix %>%  
    limma::avereps() %>%
    data.frame()
  print(dim(data_matrix))
  return(data_matrix)
}

## 2.3 数据清洗----
# 需修改*_matrix.txt.gz文件
GSE_matrix_1 <- as.data.frame(exprs(aset_1$GSE84437_series_matrix.txt.gz))
range(GSE_matrix_1)
#去除NA值
GSE_matrix_1 <- na.omit(GSE_matrix_1)
# #负值处理
# GSE_matrix_1 <- GSE_matrix_1 + abs(apply(GSE_matrix_1,1,min))
range(GSE_matrix_1)
GSE_matrix_1 <- log2_if_need(GSE_matrix_1)
range(GSE_matrix_1)

# GSE_matrix_2 <- as.data.frame(exprs(aset_2$GSE98566_series_matrix.txt.gz))
# range(GSE_matrix_2)
# GSE_matrix_2 <- na.omit(GSE_matrix_2)
# GSE_matrix_2 <- log2_if_need(GSE_matrix_2)
# range(GSE_matrix_2)

# GSE_matrix_3 <- as.data.frame(exprs(aset_3$GSE29221_series_matrix.txt.gz))
# range(GSE_matrix_3)
# GSE_matrix_3 = GSE_matrix_3 - apply(GSE_matrix_3,1,min)
# GSE_matrix_3 <- na.omit(GSE_matrix_3)
# GSE_matrix_3 <- log2_if_need(GSE_matrix_3)
# range(GSE_matrix_3)

## 2.4 整理分组信息----
PD_1 <- pData(aset_1$GSE84437_series_matrix.txt.gz)
# # 如果有研究疾病之外的其他疾病样本，运行下面一行代码去除
# PD_1 <- PD_1[PD_1$characteristics_ch1!="subject id: s6234",]  # 修改：去除抗性个体（3行代码）
# PD_1 <- PD_1[PD_1$characteristics_ch1!="subject id: s6207",]
# PD_1 <- PD_1[PD_1$characteristics_ch1!="subject id: s6341",]
# 
# 
# PD_1 <- PD_1[PD_1$characteristics_ch1.3=="study phase: experimental",]  # 修改：只保留处理组行

group_1 <- geo_group(pdata = PD_1,
                    # 样本编号列名；分组信息列名（修改）
                    sample_id = "geo_accession", group_id = "characteristics_ch1", 
                    # 疾病组命名；对照组命名；是否降序可修改）
                    treat_name = treat_name, con_name = con_name, desc = F, 
                    # 疾病组关键字符串；对照组关键字符串（修改）
                    treat_key = "cancer", con_key = "non2") ##只保留疾病组

# PD_2 <- pData(aset_2$GSE98566_series_matrix.txt.gz)
# 
# # 如果有研究疾病之外的其他疾病样本，运行下面一行代码去除
# PD_2 <- PD_2[PD_2$characteristics_ch1!="subject id: s6234",]  # 修改：去除抗性个体（3行代码）
# PD_2 <- PD_2[PD_2$characteristics_ch1!="subject id: s6207",]
# PD_2 <- PD_2[PD_2$characteristics_ch1!="subject id: s6341",]
# 
# 
# PD_2 <- PD_2[PD_2$characteristics_ch1.3=="study phase: experimental",]  # 修改：只保留处理组行
# 
# 
# #PD_2 <- PD_2[PD_2$`disease status:ch1`!="SSD patient",]  # 修改：去除SSD样本
# group_2 <- geo_group(pdata = PD_2,
#                      # 样本编号列名；分组信息列名（修改）
#                      sample_id = "geo_accession", group_id = "title", 
#                      # 疾病组命名；对照组命名；是否降序（修改）
#                      treat_name = treat_name, con_name = con_name, desc = F, 
#                      # 疾病组关键字符串；对照组关键字符串（修改）
#                      treat_key = "Sleep", con_key = "Control") 
# 
# table(group_2$group)
# PD_3 <- pData(aset_3$GSE29221_series_matrix.txt.gz)
# group_3 <- geo_group(pdata = PD_3,
#                      # 样本编号列名；分组信息列名（修改）
#                      sample_id = "geo_accession", group_id = "title", 
#                      # 疾病组命名；对照组命名；是否降序（修改）
#                      treat_name = treat_name, con_name = con_name, desc = F, 
#                      # 疾病组关键字符串；对照组关键字符串（修改）
#                      treat_key = "Diabetic", con_key = "Non-diabetic") 

##整理临床信息
clinal_1 <- PD_1[,c(10:13,38,36,37)]
colnames(clinal_1) <- c("Group","Age","Gender","Stage_T","Stage_N","OS","OS.time")

# clinal_1$OS <- gsub("1", "dead", clinal_1$OS)
# clinal_1$OS <- gsub("0", "alive", clinal_1$OS)
# clinal_1$OS.time <- sub(".*?:\\s*", "", clinal_1$OS.time)
clinal_1$Age <- sub(".*?:\\s*", "", clinal_1$Age)
clinal_1$Gender <- sub(".*?:\\s*", "", clinal_1$Gender)
clinal_1$Stage_T <- sub(".*?:\\s*", "", clinal_1$Stage_T)

clinal_1$Group <- "STAD"

clinal_1$OS.time <- as.numeric(clinal_1$OS.time)  # 转换为数值
clinal_1$OS.time <- clinal_1$OS.time * 30  # 月 → 天
clinal_1 <- subset(clinal_1, OS.time != 0)        # 删除 OS.time 为 0 的行

library(tibble)
clinal_1 <- rownames_to_column(clinal_1, var = "geo_accession")  # 行名 → 新列"ID"

# 删除NA
clinal_1 <- na.omit(clinal_1)
# 删除空格（假设检查"Status"和"OS"列）
# clinal_1 <- clinal_1[clinal_1$Status != " " & clinal_1$OS != " ", ]

# 输出结果----
write.csv(clinal_1,"output/GSE_PD.csv",row.names = F)

##只保留group文件中有临床信息的样本
group_1 <- group_1[group_1$ID %in% clinal_1$geo_accession, ]

#根据修改后的数据集分组修改信息对数据集矩阵进行校正修改
GSE_matrix_1 <- GSE_matrix_1[,group_1$ID]
# GSE_matrix_2 <- GSE_matrix_2[,group_2$ID]

#再次确认
identical(colnames(GSE_matrix_1),group_1$ID)
# identical(colnames(GSE_matrix_2),group_2$ID)




## 2.5 探针转换,同一基因取均值----
GPL_1 <- Table(getGEO(unique(PD_1$"platform_id"), 
                      destdir = "input", AnnotGPL = F))##
# GPL_1 <- Table(getGEO(filename = "input/GPL10558.annot.gz", AnnotGPL = F))
# 探针和基因ID视情况修改为相应名称
gpl_1 <- gpl_pro(gpl_name = GPL_1, probe_id = "ID", gene_id = "Symbol")
GSE_matrix_processed_1 <- geo_matrix(data_matrix = GSE_matrix_1, 
                                     data_group = group_1, gpl_name = gpl_1)



#GPL_2 <- Table(getGEO(unique(PD_2$"platform_id"), 
#                      destdir = "input", AnnotGPL = F))
# GPL_2 <- Table(getGEO(filename = "input/GPL6244.annot.gz", AnnotGPL = F))
# 
# gpl_2 <- gpl_pro(gpl_name = GPL_2, probe_id = "ID", gene_id = "Gene symbol")
# GSE_matrix_processed_2 <- geo_matrix(data_matrix = GSE_matrix_2, 
#                                      data_group = group_2, gpl_name = gpl_2)


# GPL_3 <- Table(getGEO(unique(PD_3$"platform_id"), 
#                       destdir = "input", AnnotGPL = F))
# ##手动下载GPL_3 <- Table(getGEO(filename = "input/GPL6947.annot.gz", AnnotGPL = F))
# gpl_3 <- gpl_pro(gpl_name = GPL_3, probe_id = "ID", gene_id = "Gene symbol")
# GSE_matrix_processed_3 <- geo_matrix(data_matrix = GSE_matrix_3, 
#                      data_group = group_3, gpl_name = gpl_3)

# 3 数据集合并（两个）----
## 3.1 保留数据集中共有基因的数据部分----
# same_gene <- intersect(rownames(GSE_matrix_processed_1), 
#                        rownames(GSE_matrix_processed_2))
# comGSE_1 <- GSE_matrix_processed_1[same_gene,]
# comGSE_2 <- GSE_matrix_processed_2[same_gene,]
# 
# identical(rownames(comGSE_1), rownames(comGSE_2))
# count_matrix <- cbind(comGSE_1, comGSE_2) %>% na.omit()
# boxplot(count_matrix, outline = F, las = 2) # 查看数据集间是否有批次
# 
# 
# ## 3.2 去除批次效应----
# # "1"代表 data1,"2"代表 data2
# batch <- c(rep("1", length(group_1$group)),
#            rep("2", length(group_2$group)))
# 
# adjusted_counts <- ComBat(count_matrix, batch = batch) 
# boxplot(adjusted_counts, outline = F, las = 2) # 再次查看数据集间是否有批次
# 
# # 标准化
# adjusted_counts <- normalizeBetweenArrays(adjusted_counts)
# # 去除方差为0的行
# adjusted_counts <- adjusted_counts[apply(adjusted_counts, 1, var) != 0,]
# adjusted_counts <- as.data.frame(adjusted_counts)
# boxplot(adjusted_counts, outline = F, las = 2)
# 
# ## 3.3 整理合并后的分组文件----
# combined_pd <- rbind(group_1,group_2)
# adjusted_pd <- combined_pd[order(combined_pd$group,decreasing = F),]
# 
# # 根据合并数据集的分组信息对合并数据集样本进行重新排列
# adjusted_counts<-adjusted_counts[,adjusted_pd$ID]
# # 再次确认
# identical(colnames(adjusted_counts),adjusted_pd$ID)

# 4 绘制校正前后boxplot,PCA图----
## 4.1 定义绘图函数----
# mytheme
mytheme <- 
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        plot.title = element_text(size = 7)) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, size = 0.75 * 0.47)) +
  theme(axis.line = element_line(size = 0.75 * 0.47),
        axis.text = element_text(size = 6, color = "black"), 
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.75 * 0.47)) +
  theme(legend.key = element_rect(color = NA,fill = NA),
        legend.key.size = unit(c(0.3, 0.3), "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(),
        legend.box.margin = margin(),
        legend.box.spacing = unit(0, "cm"),
        legend.background = element_blank(), legend.spacing = unit(0, "cm"),
        legend.box.background = element_blank())


# 转长数据
dat_to_long <- function(raw_dat,combined_pd,cla){
  raw_dat <- cbind(rownames(raw_dat),raw_dat)
  rownames(raw_dat) <- NULL
  colnames(raw_dat)[1] <- "Gene"
  rownames(combined_pd) <- NULL
  raw_dat<-raw_dat[,-1]

  dat_long <- raw_dat %>%
    t() %>%
    data.frame() %>%
    dplyr::mutate(sample = rownames(.)) %>%
    dplyr::mutate(group = cla) %>%
    dplyr::mutate(sample = factor(.$sample, levels = .$sample)) %>%
    tidyr::gather(key = geneid, value, - c(sample, group))
  dat_long$sample <- factor(dat_long$sample,
                            levels = as.vector(unique(dat_long$sample)))
  return(dat_long)
}

# Boxplot绘图
dat_boxplot <- function(dat_long,GSEID,desc,group_col,boxplot_width,boxplot_high)
{
boxplot <- 
  ggplot(dat_long, aes(sample, value)) + 
  geom_boxplot(aes(fill = group), lwd = 0.1 *0.47, outlier.shape = NA) + 
  scale_fill_manual(values = group_col) +
  labs(title = paste0(GSEID," ",desc," Normalization")) + 
  mytheme + # 自己的主题
  theme(legend.position = 'top', 
        legend.direction = "horizontal", 
        legend.title = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  scale_y_continuous(limits = round(boxplot.stats(dat_long$value)$stats[c(1, 5)]*1.05, 2)) +
  coord_cartesian(clip = "off")
boxplot
ggsave(paste0("output/", GSEID, " Boxplot ", desc, ".pdf"), 
       plot = boxplot, units = "cm", 
       width = boxplot_width, height = boxplot_high)
}     

# # PCA绘图
# pca_plot <- function(data,GSEID,desc,group_col,pca_width,pca_height,cla){
#   data <- as.matrix(t(data))
#   iris.pca<- PCA(data,graph =F,scale.unit = TRUE)
#   pdf(paste0("output/", GSEID," PCA ", desc, ".pdf"), 
#       width = pca_width, height = pca_height)
#   ind.p <- fviz_pca_ind(iris.pca, geom.ind = "point", col.ind = cla,
#                               palette =group_col,
#                               addEllipses = TRUE,
#                               legend.title = "Group")
#   ind.p <- ggpubr::ggpar(ind.p, title = paste(GSEID,desc, "Normalization"), 
#                 ggtheme = theme_bw())
#   print(ind.p)
#   dev.off()
# }


## 4.2 画校正前的boxplot图----
# group_col_1 = c("#FDD692","#A593E0")  # 修改：自定义绘图分组颜色
# cla <- c(rep(GSE_1, length(group_1$group)),
#          rep(GSE_2, length(group_2$group)))
# 
# dat_before_long <- dat_to_long(count_matrix,combined_pd,cla = cla)
# boxplot_before <- dat_boxplot(dat_before_long,GSEID = "Combined Datasets",
#                               desc = "Before", group_col_1,
#                               boxplot_width=8.1, boxplot_high=5)
# 
# ## 4.3 画校正后的boxplot图----
# dat_after_long <- dat_to_long(adjusted_counts,adjusted_pd,cla)
# boxplot_after <- dat_boxplot(dat_after_long,GSEID = "Combined Datasets",
#                              desc = "After",group_col_1,
#                              boxplot_width=8.1, boxplot_high=5)
# 
# ## 4.4 绘制校正前后的PCA图----
# pca_before <- pca_plot(count_matrix,GSEID = "Combined Datasets",
#                        desc = "Before",group_col=group_col_1,
#                        pca_width = 6, pca_height = 4.5,cla=cla)
# pca_after <- pca_plot(adjusted_counts,GSEID = "Combined Datasets",
#                       desc = "After",group_col=group_col_1,
#                       pca_width = 6, pca_height = 4.5,cla=cla)

# 5 绘制训练集和验证集数据boxplot----
## 5.1训练集--
# 画校正前的boxplot图
group_col_1 = c("#FDD692")  # 修改：只有疾病
gs <- factor(group_1$group, levels = c(con_name,treat_name))

test_before_long <- dat_to_long(GSE_matrix_processed_1,group_1,gs)
boxplot_before <- dat_boxplot(test_before_long,GSEID = GSE_1,
                              desc = "Before", group_col_1,
                              boxplot_width=16, boxplot_high=5)

# 标准化
GSE_matrix_adjusted_1 <- normalizeBetweenArrays(GSE_matrix_processed_1)
GSE_matrix_adjusted_1 <- GSE_matrix_adjusted_1[apply(GSE_matrix_adjusted_1,
                                                     1, var) != 0,]
GSE_matrix_adjusted_1 <- as.data.frame(GSE_matrix_adjusted_1)

# 画校正后的boxplot图
test_after_long <- dat_to_long(GSE_matrix_adjusted_1,group_1,gs)
boxplot_after <- dat_boxplot(test_after_long,GSEID = GSE_1,
                             desc = "After",group_col_1,
                             boxplot_width=16, boxplot_high=5)

# 
# ## 绘制校正前后的PCA图----
# pca_before <- pca_plot(GSE_matrix_processed_1,GSEID = "GSE76427",
#                        desc = "Before",group_col=group_col_1,
#                        pca_width = 6, pca_height = 4.5,cla=gs)
# pca_after <- pca_plot(GSE_matrix_adjusted_1,GSEID = "GSE76427",
#                       desc = "After",group_col=group_col_1,
#                       pca_width = 6, pca_height = 4.5,cla=gs)


# ## 5.2验证集--
# # 画校正前的boxplot图
# group_col_1 = c("#8EC0E4","#ED9282")  # 修改：自定义绘图分组颜色
# gs <- factor(group_2$group, levels = c(con_name,treat_name))
# 
# test_before_long <- dat_to_long(GSE_matrix_processed_2,group_2,gs)
# boxplot_before <- dat_boxplot(test_before_long,GSEID = GSE_2,
#                               desc = "Before", group_col_1,
#                               boxplot_width=8.1, boxplot_high=5)
# 
# # 标准化
# GSE_matrix_adjusted_2 <- normalizeBetweenArrays(GSE_matrix_processed_2)
# GSE_matrix_adjusted_2 <- GSE_matrix_adjusted_2[apply(GSE_matrix_adjusted_2,
#                                                      1, var) != 0,]
# GSE_matrix_adjusted_2 <- as.data.frame(GSE_matrix_adjusted_2)
# 
# # 画校正后的boxplot图
# test_after_long <- dat_to_long(GSE_matrix_adjusted_2,group_2,gs)
# boxplot_after <- dat_boxplot(test_after_long,GSEID = GSE_2,
#                              desc = "After",group_col_1,
#                              boxplot_width=8.1, boxplot_high=5)
# 
# 
# ## 绘制校正前后的PCA图----
# pca_before <- pca_plot(GSE_matrix_processed_2,GSEID = "GSE98566",
#                        desc = "Before",group_col=group_col_1,
#                        pca_width = 6, pca_height = 4.5,cla=gs)
# pca_after <- pca_plot(GSE_matrix_adjusted_2,GSEID = "GSE98566",
#                       desc = "After",group_col=group_col_1,
#                       pca_width = 6, pca_height = 4.5,cla=gs)
# 6 输出结果----
## 6.1 保存为csv文件----
# write.csv(adjusted_pd,"output/Combined_Datasets_Group.csv",row.names = F)
# write.csv(adjusted_counts,"output/Combined_Datasets_Matrix.csv")

write.csv(group_1, 
          file = paste0("output/",GSE_1,"_Datasets_Group.csv"), row.names = F)
write.csv(GSE_matrix_adjusted_1, 
          file = paste0("output/",GSE_1,"_Datasets_Matrix.csv"))

# write.csv(group_2, 
#           file = paste0("output/",GSE_2,"_Datasets_Group.csv"), row.names = F)
# write.csv(GSE_matrix_adjusted_2, 
#           file = paste0("output/",GSE_2,"_Datasets_Matrix.csv"))



## 6.2 保存为工作环境----
save.image("my_workspace.RData.gz",compress = "gzip")

