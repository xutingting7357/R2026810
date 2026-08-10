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
if(dir.exists("output")!=TRUE){
  dir.create("output")
}

if(dir.exists("intput")!=TRUE){
  dir.create("input")
}

## 1.5 定义变量----
tumor_name <- "STAD" # 肿瘤缩写命名

# 2 读取数据----
# 下载自“仙桃学术”生信工具中的临床意义-预后类-[云]生存曲线(KM图)；“分析数据.xlsx”
# 纯代码可以注释掉，如果是代码+仙桃共同使用的话最好选上用于统一样本
# id <- readxl::read_xlsx("input/sampleid.xlsx")
clinical_old <- readxl::read_xlsx(paste0("input/STAD_rnaseq_clinical_raw.xlsx"))

#去掉没有分组信息的样本
clinical_old <- clinical_old[!is.na(clinical_old$status),]
table(clinical_old$status)

# 主要的预后信息（OS、OS.time）
# 主要的临床变量（age、gender、TMN分期、其他临床stage/分级），也可以用仙桃临床意义 - 临床相关 - [云] 基线资料表板块快速查看可选的临床变量
clinical_new <- data.frame(
  #固定不变
  ID = clinical_old$RNAseq样本编号,
  Group = clinical_old$status,
  Sample = clinical_old$sample,
  OS =  clinical_old$OS,
  OS.time = clinical_old$OS.time.y,
  DSS = clinical_old$DSS,
  DSS.time = clinical_old$DSS.time,
  # PFI = clinical_old$PFI,
  # PFI.time = clinical_old$PFI.time,  ##没有PFI生存时间
  Age = clinical_old$age_at_initial_pathologic_diagnosis.x,
  Gender = clinical_old$gender.x,
  #可选
  Stage_Pathologic = clinical_old$stage_event.pathologic_stage,
  Stage_M = clinical_old$stage_event.tnm_categories.pathologic_categories.pathologic_M,
  Stage_N = clinical_old$stage_event.tnm_categories.pathologic_categories.pathologic_N,
  Stage_T = clinical_old$stage_event.tnm_categories.pathologic_categories.pathologic_T
)

## 2.1 整理肿瘤样本信息----
tumor_sample <- clinical_new[clinical_new$Group == "Tumor",]

## 2.2 整理对照样本信息----
Control <- clinical_new[clinical_new$Group == "Normal",]

# 3 根据仙桃临床信息过滤样本----
# 纯代码可以注释掉，如果是代码+仙桃共同使用的话最好选上用于统一样本
# tumor_sample <- tumor_sample[tumor_sample$sample_id %in% id$sample_id,]

# 对肿瘤和对照分组重命名
tumor_sample$Group <- gsub("Tumor",tumor_name,tumor_sample$Group)

#有时候对比的不一定是正常样本，可能是癌旁，这里统一成Control
Control$Group <- gsub("Normal","Control",Control$Group)

# 去除缺少生存信息的样本
tumor_sample <- tumor_sample %>% 
  filter(!is.na(OS) & !is.na(OS.time) & OS.time > 0)

# 合并肿瘤和对照样本
PD <- rbind(Control,tumor_sample)

# 查看临床变量是否有异常或者可以合并
table(PD$Age)
table(PD$Gender)
table(PD$Stage_M)
table(PD$Stage_N)
table(PD$Stage_T)
table(PD$Stage_Pathologic)

## 3.1 Stage_Pathologic清洗----
#特殊的肿瘤分期信息合并(正则表达式+gsub嵌套)
PD$Stage_Pathologic <- gsub("^(Stage IV)(A|B|C).*", "\\1",
                            gsub("^(Stage III)(A|B|C).*", "\\1",
                                 gsub("^(Stage II)(A|B|C).*",  "\\1",
                                      gsub("^(Stage I)(A|B|C).*", "\\1", PD$Stage_Pathologic))))

# Not Available和Discrepancy样本标记为空
PD$Stage_Pathologic <- gsub("Not Available|Discrepancy","",PD$Stage_Pathologic)

# 或者删掉
# PD <- PD[!grepl("Not Available|Discrepancy", PD$Stage_Pathologic), ] 

# # 合并Stage III,Stage IV样本(样本少的时候选择合并)
# PD$Stage_Pathologic <- gsub("^(Stage III).*|^(Stage IV)", "Stage III&IV", PD$Stage_Pathologic)
table(PD$Stage_Pathologic)

## 3.2 TMN清洗----
#特殊的肿瘤分期信息合并(正则表达式+gsub嵌套)
PD$Stage_T <- gsub("^T4.*", "T4",
                   gsub("^T3.*", "T3",
                        gsub("^T2.*", "T2",
                             gsub("^T1.*", "T1", 
                                  gsub("^T0.*", "T0", PD$Stage_T)))))
PD$Stage_N <- gsub("^N3.*", "N3",
                   gsub("^N2.*", "N2",
                        gsub("^N1.*", "N1", 
                             gsub("^N0.*", "N0", PD$Stage_N))))
PD$Stage_M <- gsub("^M1.*", "M1", 
                   gsub("^M0.*", "M0", PD$Stage_M))

# Not Available和后缀带X的未知分类样本标记为空
PD$Stage_T <- gsub("Discrepancy|Not Available|Tis|*.X$","",PD$Stage_T)
PD$Stage_M <- gsub("Discrepancy|Not Available|*.X$","",PD$Stage_M)
PD$Stage_N <- gsub("Discrepancy|Not Available|*.X$","",PD$Stage_N)

# 或者刪掉
# PD <- PD[!grepl("Not Available|X$", PD$Stage_T) 
#          & !grepl("Not Available|X$", PD$Stage_N) 
#          & !grepl("Not Available|X$", PD$Stage_M), ]

table(PD$Stage_T)
table(PD$Stage_M)
table(PD$Stage_N)

# # 合并T2,T1样本(样本少的时候选择合并)
# PD$Stage_T <- gsub("^T1.*|^T2.*", "T1&T2", PD$Stage_T)
# table(PD$Stage_T)
# 
# # # 合并T3,T4样本(样本少的时候选择合并)
# PD$Stage_T <- gsub("^T3.*|^T4.*", "T3&T4", PD$Stage_T)
# table(PD$Stage_T)

## 3.3 age清洗----
# Not Available的未知分类样本标记为空
PD$Age <- gsub("Not Available","",PD$Age) 

# 或者刪掉
# PD <- PD[!grepl("Not Available", PD$Age), ]

# table(PD$Age)
# # 定义一个年龄分组（老年、中年、青年、少年都可以自己定义，先默认定义老年65岁）
# PD$Age <- case_when(PD$Age == "" ~ "",
#                     PD$Age >= 65 ~ ">=65",
#                     PD$Age <65 ~ "<65")
# table(PD$Age)

# 4 输出结果----
# 将分组信息和临床信息保存
table(PD$Group)
write.csv(PD,"output/PD.csv",row.names = F)
group <- PD[,c("ID","Group")]
colnames(group)=c("ID","status")#重新命名
write.csv(group,"output/sample_group.csv",row.names = F)

