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
required_packages <- c("org.Hs.eg.db")

#调用函数进行检查和安装
install_if_missing(required_packages)

## 1.3 加载R包----
library(org.Hs.eg.db)

# 2 处理数据----
## 2.1 输入基因信息----
genes <- data.table::fread("input/8-LASSO_hubgenes.csv",data.table = F)

## 2.2 将基因SYMBOL转化为ENSEMBL----
genesymbol <- clusterProfiler::bitr(genes$x,fromType = "SYMBOL",
                                    toType = "ENSEMBL",OrgDb = org.Hs.eg.db)


# 3 mRNA-Drug----

## 3.1 循环下载所有hub genes相关文件----
for (i in 1:nrow(genesymbol)){
  a <- as.character(genesymbol[i,])
  file = paste("./output/",a[1],".txt",sep = "")
  link = paste("https://ctdbase.org/query.go?chem=&d-1339283-e=1&gene=",
               a[1],"&pathwayqt=equals&taxonqt=equals&chemqt=equals&go=&action",
               "DegreeTypes=increases&actionDegreeTypes=decreases&actionDegree",
               "Types=affects&sort=chemNmSort&type=ixn&actionTypes=ANY&perPage=",
               "50&pathway=&action=Search&taxon=TAXON%3A9606&goqt=",
               "equals&6578706f7274=1&geneqt=equals",sep = "")#下载链接
  download.file(link,file)
}
## 3.2 将所有Drug文件合并----
file_list <- list.files("output/")
file_path <- paste("output/",file_list,sep = "")
Drug_CTD <- data.table::fread(file_path[1],header = T)
Drug_CTD$genesymbol <- genesymbol$SYMBOL[1]
for(i in 2:length(file_path)){
  Drug <- data.table::fread(file_path[i],header = T)
  Drug$genesymbol <- stringr::str_extract(file_list[i],".*(?=\\.)")
  Drug_CTD <- rbind(Drug_CTD,Drug)
}
## 3.3 根据Reference Count的数据进行筛选----
table(Drug_CTD$`Reference Count`)
mRNA_Drug <- Drug_CTD[which(Drug_CTD$`Reference Count`>2),c(10,1)]#可修改阈值
colnames(mRNA_Drug)[1] <- "node1"
colnames(mRNA_Drug)[2] <- "node2"
mRNA_Drug <- mRNA_Drug[!duplicated(mRNA_Drug)]

# 4 输出结果保存----
## 4.1保存csv文件和txt文件
write.table(mRNA_Drug,"output/15-mRNA-Drug_nodes.txt",row.names = F, quote = F,sep = "\t")
Drug_attr <- data.frame("node" = c(unique(mRNA_Drug$node1),unique(mRNA_Drug$node2)),
                        "attribute" = c(rep("mRNA",length(unique(mRNA_Drug$node1))),
                                        rep("Drug",length(unique(mRNA_Drug$node2)))))
write.table(Drug_attr,"output/16-mRNA-Drug_attribute.txt",row.names = F,quote = F,sep = "\t")
table_Drug <- mRNA_Drug
colnames(table_Drug)[1] <- "mRNA"
colnames(table_Drug)[2] <- "Drug"
write.csv(table_Drug,"output/17-TableS5 mRNA-Drug.csv",row.names = F,quote = F)

## 4.2 保存工作环境----
table(Drug_attr$attribute)
save(mRNA_Drug,Drug_attr,file = "output/18-mRNA_Drug.RData")


