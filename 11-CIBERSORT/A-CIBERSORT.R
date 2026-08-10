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
required_packages <- c("ggplot2", "ggthemes", "ggpubr","corrplot",'Rmisc',
                       'plyr','Hmisc')
 
# 调用函数进行检查和安装
install_if_missing(required_packages)

## 1.3 加载R包----
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(corrplot)
library(Rmisc)
library(plyr)
library(Hmisc)
library(dplyr)
source('input/A-CIBERSORT_Source.R')

## 1.4 创建文件夹----
dir.create("input/")
dir.create("output/")

## 1.5 定义变量----
treat_name = "STAD"    ## 疾病组命名
con_name = "Control"  ## 对照组命名

# 2 免疫浸润----
## 2.1 读取文件 ----
dat <- data.table::fread("input/Matrix_FPKM.csv",data.table = F)
row.names(dat) <- dat[,1]
dat <- dat[,-1]
group <- data.table::fread("input/PD.csv", data.table = F, header = T) %>%
  dplyr::select("ID","Group")

table_cell <- read.table('input/LM22_input.txt',fill=T,header=T,sep='\t',
                         row.names=1,check.names=F)

## 2.2 CIBERSORT ----
data <- dat[intersect(rownames(dat),rownames(table_cell)),]
table_cell <- table_cell[intersect(rownames(data),rownames(table_cell)),]
# 读取已保存的免疫浸润结果文件
result <- readRDS("output/cibersort_result.Rds")
# 首次做运行以下两行代码
# result <- CIBERSORT(table_cell, data, perm = 1000, QN = T)
# saveRDS(result,file = 'output/cibersort_result.Rds')

cibersort <- as.data.frame(result)[,1:22]
colnames(cibersort)[9] <- "T cells regulatory Tregs"
cibersort <- cibersort[group$ID,] 

# 3 可视化----
## 3.1 成分柱状图 ----
ciber_filter <- cibersort[,colSums(cibersort) > 0]
conNum=nrow(group[group$Group == con_name,])
treatNum=nrow(group[group$Group ==	treat_name,])
colA <- "#a5dff9" # Control组颜色
colB <- "#FDD692" # Disease组颜色

pdf("output/1-Barplot.pdf",16.6,5)
set.seed(500)
randomcoloR::distinctColorPalette(22) 
mycol <- alpha(rainbow(ncol(ciber_filter)), 0.4) #创建彩虹色板（带60%透明度）
par(bty="o", mgp = c(2.5,0.3,0),
    mar = c(4.1,4.1,2.1,10.1),
    tcl=-.25,las = 1,xpd = F)
a1 <- barplot(as.matrix(t(ciber_filter)),space = 0,
              names.arg = rep("",nrow(ciber_filter)),
              yaxt = "n",
              main = bquote(''),
              ylab = "Relative percentage",
              col = mycol,
              border = NA)
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), 
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2], # 这里-20要根据实际出图的图例位置情况调整
       par("usr")[4], 
       legend = colnames(ciber_filter), 
       xpd = T,fill = mycol,cex = 0.8, 
       border = NA, y.intersp = 1,
       x.intersp = 0.2,bty = "n")
par(srt=0,xpd=T)
rect(xleft = 0, ybottom = -0.01, xright = a1[conNum]+(a1[1]-0), 
     ytop = -0.06,col=colA,border = "black")
text(a1[conNum]/2,-0.033,con_name,cex=1)
rect(xleft = a1[conNum]+(a1[1]-0), ybottom = -0.01, 
     xright =a1[length(a1)]+(a1[1]-0) , ytop = -0.06,col=colB,border = "black")
text((a1[length(a1)]+a1[conNum])/2,-0.033,treat_name,cex=1)
dev.off()


## 3.2 免疫细胞分组比较图 ----
#进行显著性检验并标记p值
boxplot_p <- function(data,vs,group){#vs：要进行比较的列名。
  formulas <- sapply(vs,
                     function(x) as.formula(paste(x,"~",group)))
  group <- data[[group]]
  if(length(levels(factor(group)))>2){
    test=lapply(formulas, function(x){kruskal.test(x, data = data)})
  }else{
    test=lapply(formulas, function(x){wilcox.test(x, data = data)})
  }
  test_results <- lapply(test,
                         function(x){ 
                           p.value<-signif(x$p.value, digits=3)
                           pstar <-ifelse(x$p.value<0.05,
                                          ifelse(x$p.value<0.01,"**","*"),"ns")
                           res <- c(p.value,pstar)
                           return(res)
                         })
  
  res <- t(as.data.frame(test_results, check.names = FALSE))
  colnames(res) <- c("p.value","pstar")
  return(res)
}

dat_boxplot <- as.data.frame(cibersort)
identical(group$ID,rownames(dat_boxplot))
dat_boxplot <- cbind(group$Group,dat_boxplot)
colnames(dat_boxplot) <- gsub(" ","_",colnames(dat_boxplot))
colnames(dat_boxplot)[1] <- "Group"

dat_boxplot_p <- boxplot_p(dat_boxplot,colnames(dat_boxplot)[-1],"Group")

dat_boxplot_p <- as.data.frame(dat_boxplot_p)
dat_boxplot_p <- na.omit(dat_boxplot_p)
length(which(dat_boxplot_p$pstar != "ns"))

#输出差异显著的免疫细胞
gene_sig <- paste0(rownames(dat_boxplot_p)[dat_boxplot_p$pstar != "ns"],collapse = "，")
cat(gene_sig, file = "output/报告.txt")

# 绘制显著差异免疫细胞
# ciber_boxplot <- ciber_filter[,dat_boxplot_p$pstar != "ns"]
#保留在分组比较图中全部细胞
ciber_boxplot <- ciber_filter
immu<- rep(colnames(ciber_boxplot),each=nrow(ciber_boxplot)) #组别变量
immu <- factor(immu)
a<-c(group$Group)
boxplot_group <- rep(a,ncol(ciber_boxplot))
boxplot_group <- factor(boxplot_group,levels = c(con_name,treat_name))
value <- c()
for (j in 1:ncol(ciber_boxplot)) { value<-c(value,ciber_boxplot[,j])}
value<-as.numeric(value)
boxplot_dat <- data.frame(immu_cell=immu,group=group,value=value)

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
  theme(legend.key = element_rect(fill = NA,color = NA),
        legend.key.size = unit(c(0.3, 0.3), "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(),
        legend.box.margin = margin(),
        legend.box.spacing = unit(0, "cm"),
        legend.background = element_blank(), legend.spacing = unit(0, "cm"),
        legend.box.background = element_blank())

boxplot <- 
  ggplot(boxplot_dat,aes(x=immu_cell,y=value,fill=boxplot_group)) +
  geom_boxplot(width=0.7*0.47,size=0.3*0.47,outlier.color = NA) +
  # 修改肿瘤对照组名称
  scale_fill_manual(values = c("STAD" = colB, "Control" = colA)) +
  mytheme +
  stat_compare_means(method	= "wilcox.test",
                     symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")),
                     label = "p.signif",size = 2)+
  theme(axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position = "top",
        legend.direction = "horizontal", # 图注方向（水平）
        legend.title = element_blank(),
        legend.text = element_text(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.75)) +
  xlab('') +
  ylab('Infiltration Abundance')
boxplot
ggsave(boxplot,file="output/2-DiffBoxplot.pdf",units = "cm",width = 16.6,height = 7.0)


## 3.3 免疫细胞相关性热图----
# #绘制显著差异免疫细胞
heatmap <- t(ciber_filter[,dat_boxplot_p$pstar != "ns"])
# 保留在分组比较图中全部细胞
# heatmap <- t(ciber_filter)
heatmap <- t(heatmap[,group$ID])
heatmap_col <- colorRampPalette(colors = c("#4DBBD5","white","#E64B35"),
                                space="Lab")

pdf("output/3-HeatMap.pdf",width = 18,height = 16)
corrplot(corr =cor(heatmap,method = 'spearman'),order="AOE",type="upper",
         tl.pos="tp",method="color",
         tl.col = "black",col=heatmap_col(50))
corrplot(corr = cor(heatmap,method = 'spearman'),add=TRUE, type="lower", 
         method="pie",order="AOE", 
         diag=FALSE,tl.pos="n", cl.pos="n",col=heatmap_col(50))
dev.off()

## 3.4 输出相关性结果为表格----
result_con_name <- data.frame(immune_cells1 = character(), 
                              immune_cells2 = character(), cor = numeric(), 
                              p.value = numeric(), stringsAsFactors = FALSE)
# 循环计算并添加结果到数据框
for (i in 1:(ncol(heatmap)-1)) {
  for (j in (i+1):ncol(heatmap)) {
    # 计算相关性和p-value
    dd <- cor.test(heatmap[,i], heatmap[,j],method = 'spearman')
    
    # 将结果添加到数据框
    result_con_name <- rbind(result_con_name, 
                             data.frame(immune_cells1 = colnames(heatmap)[i], 
                                        immune_cells2 = colnames(heatmap)[j], 
                                        cor = dd$estimate, p.value = dd$p.value))
  }
}

write.csv(result_con_name, "output/3-Heatmap.csv", row.names = F)

## 3.5 免疫细胞与关键基因的相关性热图----
gene <- data.table::fread("input/8-LASSO_hubgenes.csv",data.table = F)
hub <- gene$x
corheatmap_hub <- dat[,group$ID]
corheatmap_hub <- corheatmap_hub[hub,]

#合并免疫细胞样本矩阵与目的基因的样本矩阵，前22列为免疫细胞，后面为目的基因
corheatmap_dat <- as.data.frame(ciber_filter[,dat_boxplot_p$pstar != "ns"])
# corheatmap_dat <- as.data.frame(ciber_filter)
identical(row.names(corheatmap_dat),row.names(t(corheatmap_hub)))

corheatmap <- cbind(corheatmap_dat,t(corheatmap_hub))

cor<-data.frame()
for (mm in 10:25) { #修改：目的基因所在列
  cor1<-data.frame(0,0,0,0)
  for (i in 1:9) { #修改：免疫细胞所在列
    c<-rcorr(corheatmap[,i],corheatmap[,mm],type = 'spearman')
    cor1[i,1]<-c$r[2]
    cor1[i,2]<-c$P[2]
    cor1[i,3]<-colnames(corheatmap)[i]
    cor1[i,4]<-colnames(corheatmap)[mm]
  }
  cor<-rbind(cor,cor1)
}
colnames(cor)<-c("correlation","Pvalue","immu_cell","gene")
cor_filter <- cor[cor$Pvalue<0.05,]
cor_filter<-na.omit(cor_filter)
cor_filter<-cor_filter[order(cor_filter$correlation,decreasing = F),]
y=factor(cor_filter$gene)
x=factor(cor_filter$immu_cell)
# a<-(-log10(cor_filter$Pvalue))
cor_filter <- cor_filter %>% 
  dplyr::mutate(Pvalue=case_when(0.03 < Pvalue & Pvalue <= 0.05 ~ "0.03 ~ 0.05",
                                 0.01 < Pvalue & Pvalue <= 0.03 ~ "0.01 ~ 0.03",
                                 Pvalue <= 0.01 ~ "0 ~ 0.01"))

a<-factor(cor_filter$Pvalue,levels = 
            c("0.03 ~ 0.05","0.01 ~ 0.03","0 ~ 0.01"))

corheatmap <- 
  ggplot(cor_filter,aes(x,y)) + 
  geom_point(aes(size=a,color=correlation)) +
  scale_colour_gradient2(high="#E64B35",mid = "white",low="#4DBBD5") + 
  mytheme +
  # theme(plot.title = element_text(hjust = 0.5)) +
  # labs(size="-log10Pvalue",x="",y="") +
  labs(size="P value",x="",y="") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.key = element_blank(),
        legend.margin = margin(l = 0.2,unit = "cm"))
corheatmap
ggsave(corheatmap,file="output/4-CorHeatMap.pdf",unit="cm", width =16,height =12)
write.csv(cor, file = "output/4-Bubblemap.csv")

