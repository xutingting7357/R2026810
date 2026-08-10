# 1 环境变量准备----
## 1.1 清空环境----
##需要在每个代码块之前清空所有文件及缓存
rm(list=ls())
gc()
##所有用到的R包都在最开始确认下载和加载

## 1.2 R包安装----
#定义一个函数来检查并安装所需的R包
install_if_missing <- function(packages) {
  installed_packages <- rownames(installed.packages())
  for (pkg in packages) {
    if (!(pkg %in% installed_packages)) {
      BiocManager::install(pkg, dependencies = TRUE, update = FALSE)
    }
  }
}

#检查安装的包列表
##当输入代码很长超过页面一半的时候需要换行
required_packages <- c("dplyr", "stringr","ggplot2","ggpubr","FSA")

#调用函数进行检查和安装
install_if_missing(required_packages)

## 1.3 加载R包----
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(FSA)

dir.create("input")
dir.create("output")

# 2 数据处理----
## 2.1 载入数据----
lassogene <- read.csv("input/8-LASSO_hubgenes.csv")
data <- data.table::fread("input/Matrix_FPKM.csv")
data <- as.data.frame(data)
rownames(data) <- data$V1
data <- data[,-1]
pd <- read.csv("input/PD.csv")

# 空格赋值NA
pd[pd == ""] <- NA
# 删除 "Age", "Gender", "Stage_M", "Stage_N", "Stage_T"列中包含 NA 的行
pd <- pd[complete.cases(pd[, c("Age", "Gender", "Stage_M", "Stage_N", "Stage_T")]), ]


## 2.2 规范临床分期名称----
pd[pd$Age <= 60,]$Age <- "<= 60"
pd[pd$Age > 60,]$Age <- "> 60"
table(pd$Age)
table(pd$Gender)
table(pd$Stage_M)
table(pd$Stage_N)
table(pd$Stage_T)
# 提取肿瘤样本临床信息
pd <- pd[pd$Group == "STAD",]
data <- data[,pd$ID]

# 3 Age分组箱线图----
## 3.1 整理画Age分组箱线图的数据----
identical(pd$ID,colnames(data))
dat_Age <- as.data.frame(cbind(pd$Age,t(data[lassogene$x,])))

dat_boxplot_Age <- dat_Age[pd$ID,2:length(dat_Age)]
longdat_boxplot_Age <- dat_boxplot_Age %>% 
  dplyr::mutate(group = dat_Age$V1) %>% #新建一列group存放风险分组信息
  tidyr::gather(key = gene, value = value, - "group") %>% #宽格式的数据转换为长格式的数据
  dplyr::mutate(gene = factor(.$gene, levels = colnames(dat_boxplot_Age)))
mylist <- list()
longdat_boxplot_Age$group <- factor(longdat_boxplot_Age$group,#调整箱线顺序
                                levels = c("<= 60", "> 60"))
longdat_boxplot_Age$value <- as.numeric(longdat_boxplot_Age$value)

## 3.2 设置主题----
mytheme <-
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), # 图像
        plot.title = element_text(size = 7)) +
  theme(panel.background = element_blank(), # 面板
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = 0.75 * 0.47)) + # 添加外框
  theme(axis.line = element_line(linewidth = 0.75 * 0.47), # 坐标轴
        axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(linewidth = 0.75 * 0.47)) +
  theme(legend.key = element_rect(fill = NA, color = NA,linewidth = 0.75 * 0.47), # 图注
        legend.key.size = unit(c(0.3, 0.3), "cm"),
        legend.title = element_text(size = 7, vjust = 1),
        legend.text = element_text(size = 6),
        legend.margin = margin(r = 0.08,  t = 0.03, b = 0.01,
                               l = 0.08, unit = "cm"), # 每个图例周围的边距
        legend.justification = c(0.5, 0.5),
        legend.box.margin = margin(-5, -7.5, -5, -5), # 完整图例区域周围的边距
        legend.spacing = unit(0, "cm"), # 图例之间的间距
        legend.box.spacing = unit(0, "cm"), # 打印区域和图例框之间的间距
        legend.background = element_blank(),
        legend.box.background = element_blank())

## 3.3 绘制Age分组箱线图----
p1 <- ggplot(longdat_boxplot_Age, aes(x = gene, y = value, fill = group)) +
  geom_boxplot(width=0.7,size=0.1,outlier.color = NA) +
  scale_fill_manual(values = c("<= 60" = "#a5dff9","> 60" = "#FDD692")) +
  stat_compare_means(method= "wilcox.test", size = 6 * 0.35, vjust = 0.5,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")),
                     label = "p.signif") +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_text(hjust = 1, vjust = 1), 
        legend.position = 'top', 
  ) +
  xlab('') + ylab('score') + labs(fill='Group') + #分组标签的标题
  coord_cartesian(clip = "off") # 去除蒙版
p1
ggsave("output/1-boxplot_Age.pdf", plot = p1, units = "cm",
       width = 16.6, height = 7)

# 4 Gender分组箱线图----
## 4.1 整理画Gender分组箱线图的数据----
dat_Gender <- as.data.frame(cbind(pd$Gender,t(data[lassogene$x,])))
dat_boxplot_Gender <- dat_Gender[pd$ID,2:length(dat_Gender)]
longdat_boxplot_Gender <- dat_boxplot_Gender %>% 
  dplyr::mutate(group = dat_Gender$V1) %>% #新建一列group存放风险分组信息
  tidyr::gather(key = gene, value = value, - "group") %>% #宽格式的数据转换为长格式的数据
  dplyr::mutate(gene = factor(.$gene, levels = colnames(dat_boxplot_Gender)))
mylist <- list()
longdat_boxplot_Gender$group <- factor(longdat_boxplot_Gender$group,#调整箱线顺序
                                    levels = c("MALE", "FEMALE"))
longdat_boxplot_Gender$value <- as.numeric(longdat_boxplot_Gender$value)

## 4.2 绘制Gender分组箱线图----
p2 <- ggplot(longdat_boxplot_Gender, aes(x = gene, y = value, fill = group)) +
  geom_boxplot(width=0.7,size=0.1,outlier.color = NA) +
  scale_fill_manual(values = c("MALE" = "#a5dff9","FEMALE" = "#FDD692")) +
  stat_compare_means(method        = "wilcox.test", size = 6 * 0.35, vjust = 0.5,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")),
                     label = "p.signif") +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_text(hjust = 1, vjust = 1), 
        legend.position = 'top', 
  ) +
  xlab('') + ylab('score') + labs(fill='Group') + #分组标签的标题
  coord_cartesian(clip = "off") # 去除蒙版
p2
ggsave("output/2-boxplot_Gender.pdf", plot = p2, units = "cm",
       width = 16.6, height = 7)

# 5 Stage_M分组箱线图----
## 5.1 整理Stage_M分组画箱线图的数据----
dat_Stage_M <- as.data.frame(cbind(pd$Stage_M,t(data[lassogene$x,])))
dat_boxplot_Stage_M <- dat_Stage_M[pd$ID,2:length(dat_Stage_M)]
longdat_boxplot_Stage_M <- dat_boxplot_Stage_M %>% 
  dplyr::mutate(group = dat_Stage_M$V1) %>% #新建一列group存放风险分组信息
  tidyr::gather(key = gene, value = value, - "group") %>% #宽格式的数据转换为长格式的数据
  dplyr::mutate(gene = factor(.$gene, levels = colnames(dat_boxplot_Stage_M)))
mylist <- list()
longdat_boxplot_Stage_M$group <- factor(longdat_boxplot_Stage_M$group,#调整箱线顺序
                                       levels = c("M0", "M1"))
longdat_boxplot_Stage_M$value <- as.numeric(longdat_boxplot_Stage_M$value)

## 5.2 绘制Stage_M分组箱线图----
p3 <- ggplot(longdat_boxplot_Stage_M, aes(x = gene, y = value, fill = group)) +
  geom_boxplot(width=0.7,size=0.1,outlier.color = NA) +
  scale_fill_manual(values = c("M0" = "#a5dff9","M1" = "#FDD692")) +
  stat_compare_means(method = "wilcox.test", size = 6 * 0.35, vjust = 0.5,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")),
                     label = "p.signif") +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_text(hjust = 1, vjust = 1), 
        legend.position = 'top', 
  ) +
  xlab('') + ylab('score') + labs(fill='Group') + #分组标签的标题
  coord_cartesian(clip = "off") # 去除蒙版
p3
ggsave("output/3-boxplot_Stage_M.pdf", plot = p3, units = "cm",
       width = 16.6, height = 7)

# 6 Stage_N分组箱线图----
## 6.1 整理画Stage_N分组箱线图的数据----
dat_Stage_N <- as.data.frame(cbind(pd$Stage_N,t(data[lassogene$x,])))
dat_boxplot_Stage_N <- dat_Stage_N[pd$ID,2:length(dat_Stage_N)]
longdat_boxplot_Stage_N <- dat_boxplot_Stage_N %>% 
  dplyr::mutate(group = dat_Stage_N$V1) %>% #新建一列group存放风险分组信息
  tidyr::gather(key = gene, value = value, - "group") %>% #宽格式的数据转换为长格式的数据
  dplyr::mutate(gene = factor(.$gene, levels = colnames(dat_boxplot_Stage_N)))
mylist <- list()
longdat_boxplot_Stage_N$group <- factor(longdat_boxplot_Stage_N$group,#调整箱线顺序
                                       levels = c("N0","N1","N2","N3"))
longdat_boxplot_Stage_N$value <- as.numeric(longdat_boxplot_Stage_N$value)

## 6.2 绘制Stage_N分组箱线图----
# 分组列表
genes <- lassogene$x
ofset_x <- c(-0.3,-0.1,0.1,0.3)#每个分组相对于基因x轴坐标偏差，根据分组数量增删
groups <- sort(unique(longdat_boxplot_Stage_N$group))
# 创建一个没有显著性标记的箱线图
p4 <- ggplot(longdat_boxplot_Stage_N, aes(x = gene, y = value, fill = group)) +
  geom_boxplot(width=0.7,size=0.1,outlier.color = NA) +
  scale_fill_manual(values = c("N0" = "#a5dff9","N1"="#FDD692",
                               "N2" = "#ED9282","N3"="#6AAFE6")) +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_text(hjust = 1, vjust = 1),
        legend.position = 'top',
  ) +
  xlab('') + ylab('score') + labs(fill='Group') +
  coord_cartesian(clip = "off")
kk <- 0
# 对每个基因的每个分组进行显著性测试，并添加显著性标记
for (x in genes) {
  # 初始化偏移量
  kk <- 1+kk
  offset <- 0
  
  # 对所有分组进行 Kruskal-Wallis 测试
  kruskal_result <- kruskal.test(
    value ~ group,
    data = subset(longdat_boxplot_Stage_N, gene == x)
  )
  
  # 如果 Kruskal-Wallis 测试的结果显著，那么进行 Dunn 测试
  if (kruskal_result$p.value < 0.05) {
    dunn_result <- dunnTest(
      value ~ group,
      data = subset(longdat_boxplot_Stage_N, gene == x),
      method = "bonferroni"
    )
    
    # 对每一对分组进行检查
    for (i in 1:(length(groups) - 1)) {
      for (j in (i + 1):length(groups)) {
        # 找到这一对分组的 p 值
        p_value <- dunn_result$res[which(
          dunn_result$res$Comparison == paste(groups[i],
                                              groups[j], sep = " - ")), "P.adj"]
        
        # 计算标记的位置
        label_y <- max(subset(longdat_boxplot_Stage_N, 
                              gene == x)$value)-0.2 + offset
        x_p<- groups[i]
        ofset_x_p <- ofset_x+kk
        left_bracket_x <- ofset_x_p[i]
        right_bracket_x <- ofset_x_p[j]
        label_x <- (left_bracket_x + right_bracket_x) / 2
        
        # 根据 p 值决定标记
        if (p_value < 0.001) {
          label <- '***'
        } else if (p_value < 0.01) {
          label <- '**'
        } else if (p_value < 0.05) {
          label <- '*'
        } else {
          label <- 'ns'
        }
        
        # 添加显著性标记
        if(label!="ns"){
          p4 <- p4 + annotate(
            "text",
            x = label_x,
            y = label_y,
            label = label,
            vjust = -0.5,
            size = 1.75
          )
          
          p4 <- p4 + annotate(
            "segment",
            x = left_bracket_x, xend = right_bracket_x,
            y = label_y+0.5, yend = label_y+0.5,
            arrow = arrow(ends = "first", angle = 90, 
                          length = unit(0.005, "npc")),size=0.1
          ) + annotate(
            "segment",
            x = left_bracket_x, xend = right_bracket_x,
            y = label_y+0.5, yend = label_y+0.5,
            arrow = arrow(ends = "last", angle = 90, 
                          length = unit(0.005, "npc"),),size=0.1
          )}
        
        # 更新偏移量
        offset <- offset + 0.3
      }
    }
  }
}
ggsave("output/4-boxplot_Stage_N.pdf", plot = p4, units = "cm",
       width = 16.6, height = 7)


# 7 Stage_T分组箱线图----
## 7.1 整理画Stage_T分组箱线图的数据----
dat_Stage_T <- as.data.frame(cbind(pd$Stage_T,t(data[lassogene$x,])))
dat_boxplot_Stage_T <- dat_Stage_T[pd$ID,2:length(dat_Stage_T)]
longdat_boxplot_Stage_T <- dat_boxplot_Stage_T %>% 
  dplyr::mutate(group = dat_Stage_T$V1) %>% #新建一列group存放风险分组信息
  tidyr::gather(key = gene, value = value, - "group") %>% #宽格式的数据转换为长格式的数据
  dplyr::mutate(gene = factor(.$gene, levels = colnames(dat_boxplot_Stage_T)))
mylist <- list()
longdat_boxplot_Stage_T$group <- factor(longdat_boxplot_Stage_T$group,#调整箱线顺序
                                       levels = c("T1","T2","T3","T4"))
longdat_boxplot_Stage_T$value <- as.numeric(longdat_boxplot_Stage_T$value)

## 7.2 绘制Stage_T分组箱线图----
ofset_x <- c(-0.3,-0.1,0.1,0.3)#每个分组相对于基因x轴坐标偏差，根据分组数量增删
groups <- sort(unique(longdat_boxplot_Stage_T$group))

# 创建一个没有显著性标记的箱线图
p5 <- ggplot(longdat_boxplot_Stage_T, aes(x = gene, y = value, fill = group)) +
  geom_boxplot(width=0.7,size=0.1,outlier.color = NA) +
  scale_fill_manual(values = c("T1" = "#a5dff9","T2"="#FDD692",
                               "T3"="#6AAFE6","T4" = "#ED9282")) +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_text(hjust = 1, vjust = 1),
        legend.position = 'top',
  ) +
  xlab('') + ylab('score') + labs(fill='Group') +
  coord_cartesian(clip = "off")
kk <- 0
# 对每个基因的每个分组进行显著性测试，并添加显著性标记
for (x in genes) {
  # 初始化偏移量
  kk <- 1+kk
  offset <- 0
  
  # 对所有分组进行 Kruskal-Wallis 测试
  kruskal_result <- kruskal.test(
    value ~ group,
    data = subset(longdat_boxplot_Stage_T, gene == x)
  )
  
  # 如果 Kruskal-Wallis 测试的结果显著，那么进行 Dunn 测试
  if (kruskal_result$p.value < 0.05) {
    dunn_result <- dunnTest(
      value ~ group,
      data = subset(longdat_boxplot_Stage_T, gene == x),
      method = "bonferroni"
    )
    
    # 对每一对分组进行检查
    for (i in 1:(length(groups) - 1)) {
      for (j in (i + 1):length(groups)) {
        # 找到这一对分组的 p 值
        p_value <- dunn_result$res[which(
          dunn_result$res$Comparison == paste(groups[i], groups[j], 
                                              sep = " - ")), "P.adj"]
        
        # 计算标记的位置
        label_y <- max(subset(longdat_boxplot_Stage_T)$value)-0.2 + offset
        x_p<- groups[i]
        ofset_x_p <- ofset_x+kk
        left_bracket_x <- ofset_x_p[i]
        right_bracket_x <- ofset_x_p[j]
        label_x <- (left_bracket_x + right_bracket_x) / 2
        
        # 根据 p 值决定标记
        if (p_value < 0.001) {
          label <- '***'
        } else if (p_value < 0.01) {
          label <- '**'
        } else if (p_value < 0.05) {
          label <- '*'
        } else {
          label <- 'ns'
        }
        
        # 添加显著性标记
        if(label!="ns"){
        p5 <- p5 + annotate(
          "text",
          x = label_x,
          y = label_y,
          label = label,
          vjust = -0.5,
          size = 1.75
        )
        
        p5 <- p5 + annotate(
          "segment",
          x = left_bracket_x, xend = right_bracket_x,
          y = label_y+0.5, yend = label_y+0.5,
          arrow = arrow(ends = "first", angle = 90, 
                        length = unit(0.005, "npc")),size=0.1
        ) + annotate(
          "segment",
          x = left_bracket_x, xend = right_bracket_x,
          y = label_y+0.5, yend = label_y+0.5,
          arrow = arrow(ends = "last", angle = 90, 
                        length = unit(0.005, "npc"),),size=0.1
        )}
        
        # 更新偏移量
        offset <- offset + 0.3
      }
    }
  }
}
ggsave("output/5-boxplot_Stage_T.pdf", plot = p5, units = "cm",
       width = 16.6, height = 7)
# 8 输出结果保存----
## 8.1 保存工作环境----
save.image("1-correlation.RData.gz",compress = "gzip")

