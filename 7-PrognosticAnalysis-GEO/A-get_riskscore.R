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
required_packages <- c("dplyr")

#调用函数进行检查和安装
install_if_missing(required_packages)


## 1.3 加载R包----
library(dplyr)


## 1.4 创建文件夹----
dir.create("input/")
dir.create("output/")

## 1.5 定义变量----
treat_name = "STAD"    # 修改：疾病组命名
con_name = "Control"  # 对照组命名

# 2 读取表达矩阵和分组----
ES_GSE_mat <- as.data.frame(data.table::fread("input/GSE84437_Datasets_Matrix.csv"))
rownames(ES_GSE_mat) <- ES_GSE_mat$V1
ES_GSE_mat <- ES_GSE_mat[,-1]
ES_sampleinfo <- as.data.frame(data.table::fread("input/GSE_PD.csv",
                                                 header = T)) %>% 
  dplyr::select("geo_accession","Group")
colnames(ES_sampleinfo) <- c('geo_accession', 'Group')
ES_sampleinfo$Group <- factor(ES_sampleinfo$Group,levels = c(con_name,treat_name))
ES_GSE_mat <- ES_GSE_mat[,ES_sampleinfo$geo_accession]

genes <- read.table("input/8-LASSO_hubgenes.csv", header = T, sep = ',')

# 提取模型基因表达量
# colnames(genes) <- "gene_symbol"
genelist <- ''
for(i in genes$x){
  genelist <- paste0(genelist, '，', i)
}
genelist

ES_boxpgse <- as.data.frame(t(ES_GSE_mat[genes$x,]))

identical(ES_sampleinfo$geo_accession,rownames(ES_boxpgse))


lasso_coef <- readRDS("input/9-LASSO_Coef.rds")
rownames(lasso_coef) <- lasso_coef$genes
lasso_coef$coef <- as.numeric(lasso_coef$coef)

# 根据公式计算风险评分
risk.score <- apply(ES_boxpgse,1,function(x) 
  # 表达加权系数求和计算riskscore
  {crossprod(as.numeric(x),lasso_coef[colnames(ES_boxpgse),]$coef)})
ES_boxpgse$RiskScore = risk.score

ES_boxpgse$group<-c(rep(con_name,table(ES_sampleinfo$Group)[1]),
                    rep(treat_name,table(ES_sampleinfo$Group)[2]))
write.csv(ES_boxpgse,file="output/1-GEO_riskscore.csv")

# 获取高低风险组
risk_group <- ES_boxpgse[ES_boxpgse$group==treat_name,]
high.risk <- risk_group %>%
  filter(RiskScore >= median(risk_group$RiskScore))
low.risk <- risk_group %>%
  filter(RiskScore < median(risk_group$RiskScore))
risk_group1 <- rbind(high.risk,low.risk)
risk_group1$cla <- c(rep('HighRisk',216),rep('LowRisk',215)) # 手动调整

write.table(risk_group1[,ncol(risk_group1),drop = F],
            file="output/2-LASSO_RiskGroup.txt",sep="\t",row.names=T,quote=F)

