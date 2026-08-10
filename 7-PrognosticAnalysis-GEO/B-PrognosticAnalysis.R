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
required_packages <- c("rms", "dplyr", "stringr","timeROC","pROC",
                       "survival","forestplot","ggDCA","survminer","ggplot2")

#调用函数进行检查和安装
install_if_missing(required_packages)

## 1.3 加载R包----
library(rms)
library(dplyr)
library(stringr)
library(timeROC)
library(survival)
library(forestplot)
library(ggDCA)
library(survminer)
library(ggplot2)
library(tableone)
library(pROC)

# 2 数据处理----

## 2.1 载入数据----
# 加载风险评分数据，使用data.table包的fread函数读取CSV文件，返回数据框格式
riskscores <- data.table::fread("output/1-GEO_riskscore.csv", data.table = F)

# 将风险评分数据的第一列（样本ID）设置为行名，便于后续基于行名的操作
rownames(riskscores) <- riskscores[, 1]

# 删除原第一列（样本ID），因为已设置为行名
riskscores <- riskscores[, -1]

# 加载临床数据，读取包含临床信息的CSV文件
pd <- data.table::fread("input/GSE_PD.csv", data.table = F)

# 将临床数据的第一列重命名为"RNAseq"，通常表示样本的RNA测序ID
colnames(pd)[1] <- "RNAseq"

# 筛选临床数据中"Group"列为"LUAD"（肺腺癌）的样本
pd <- pd[pd$Group == "STAD", ]

# 同步筛选风险评分数据中"group"列为"LUAD"的样本，确保两组数据对应同一疾病亚型
riskscores <- riskscores[riskscores$group == "STAD", ]

# 将临床数据的行名设置为"RNAseq"列（样本ID），与风险评分数据的行名对齐
rownames(pd) <- pd$RNAseq

# 进一步筛选临床数据，仅保留风险评分数据中存在的样本（基于行名匹配）
pd <- pd[pd$RNAseq %in% rownames(riskscores), ]

# 按风险评分数据的行名顺序对临床数据进行排序，确保两个数据集的样本顺序完全一致
pd <- pd[match(rownames(riskscores), pd$RNAseq), ]

## 2.2 对临床信息进行处理----
#对数值型变量分类
#查看样本是否对齐
identical(rownames(riskscores),pd$RNAseq)

# 修改：选取相应临床信息列
data <- cbind(pd[,c(7:8,3:6)],riskscores$RiskScore)
data <- na.omit(data)

# 修改为相应列
colnames(data)[7] <- c("Risk.Score")
data[, "OS.time"] <- as.numeric(data[, "OS.time"])
data[, "OS"] <- as.numeric(data[, "OS"])


#规范临床分期名称
data$Age <- as.numeric(data$Age)  # 转换为数值
data[data$Age <= 60,]$Age <- "<= 60"
data[data$Age > 60,]$Age <- "> 60"

# 3 time_ROC分析----
# 计算时间依赖性AUC（1年、3年、5年）
# 3 time_ROC分析----
data$OS.time <- as.numeric(data$OS.time)
data$OS <- as.numeric(data$OS)
data$Risk.Score <- as.numeric(data$Risk.Score)

time_roc_res <- timeROC(T = data$OS.time,
                        delta = data$OS,
                        marker = data$Risk.Score,
                        cause = 1,
                        weighting = "marginal",
                        times = c(1 * 365, 3 * 365, 5 * 365),
                        ROC = TRUE,
                        iid = TRUE)
# 计算某时间点的 AUC 和 95% CI
get_auc_with_ci <- function(data, time_point) {
  # 提取该时间点的生存状态（1=事件发生，0=删失或未发生）
  status <- ifelse(
    data$OS.time <= time_point & data$OS == 1,  # 在指定时间前发生事件
    1, 0)                                         # 标记为 1（事件）
  roc_obj <- roc( # 使用 pROC 包计算 AUC 和 95% CI
    response = status,
    predictor = data$Risk.Score,
    levels = c(0, 1),
    direction = "<") # 假设风险评分越高，事件风险越高
  # 返回 AUC 和 95% CI
  ci <- ci.auc(roc_obj, method = "delong")
  return(c(
    auc = as.numeric(ci[2]), # AUC 点估计
    lower = as.numeric(ci[1]), # 95% CI 下限
    upper = as.numeric(ci[3]))) # 95% CI 上限
}
# 计算各时间点的 AUC 和 95% CI
auc_ci_1yr <- get_auc_with_ci(data, 1 * 365)
auc_ci_3yr <- get_auc_with_ci(data, 3 * 365)
auc_ci_5yr <- get_auc_with_ci(data, 5 * 365)

pdf("output/1-timeROC.pdf", height = 5, width = 5)
plot(time_roc_res, time = 1 * 365, col = "red", title = FALSE)
plot(time_roc_res, time = 3 * 365, add = TRUE, col = "blue")
plot(time_roc_res, time = 5 * 365, add = TRUE, col = "orange2")
legend("bottomright", #图例文本（显示 AUC 和 95% CI）
       c(paste0("1 year: AUC = ", round(auc_ci_1yr["auc"], 3), 
                " (95%CI:", round(auc_ci_1yr["lower"], 3), "-", round(auc_ci_1yr["upper"], 3), ")"),
         paste0("3 years: AUC = ", round(auc_ci_3yr["auc"], 3), 
                " (95%CI:", round(auc_ci_3yr["lower"], 3), "-", round(auc_ci_3yr["upper"], 3), ")"),
         paste0("5 years: AUC = ", round(auc_ci_5yr["auc"], 3), 
                " (95%CI:", round(auc_ci_5yr["lower"], 3), "-", round(auc_ci_5yr["upper"], 3), ")")),
       col = c("red", "blue", "orange"),
       lty = 1, lwd = 2, bty = "n")
dev.off()


# 4 KM曲线分析----
# 读取风险分组文件（示例路径："output/2-LASSO_RiskGroup.txt"）
riskstatus <- data.table::fread("output/2-LASSO_RiskGroup.txt", data.table = F)

# 合并风险评分与分组信息
riskscores <- riskscores[riskstatus$V1, ]  # 按样本ID对齐
identical(rownames(riskscores), riskstatus$V1)  # 验证行名一致性

# 添加风险等级列
riskscores <- cbind(riskscores, riskstatus$cla)
colnames(riskscores)[19] <- "Risk.Level"  # 重命名列为Risk.Level

# 准备生存分析数据（pd为原始表型数据）
pd <- pd[rownames(riskscores), ]  # 按样本ID对齐
identical(pd$RNAseq,rownames(riskscores))  # 验证一致性

# 构建KM分析数据集（OS时间、状态、风险等级）
kmdata <- cbind(pd[, 7:8], riskscores$Risk.Level)  # 假设6-7列是OS.time和OS
colnames(kmdata) <- c('OS', 'OS.time', 'Risk.Level')

# 数据清洗
kmdata <- na.omit(kmdata)  # 删除缺失值
kmdata$OS.time <- as.numeric(kmdata$OS.time)  # 确保数值类型
kmdata$OS <- as.numeric(kmdata$OS)


## 4.1 画KM曲线----
head(kmdata)
# 按风险等级排序（高风险在前）
kmdata <- kmdata[order(kmdata$Risk.Level, decreasing = T), ]

# 将风险等级转为因子（确保绘图顺序：LowRisk vs HighRisk）
kmdata$Risk.Level <- factor(kmdata$Risk.Level, levels = c("LowRisk", "HighRisk"))

# 拟合生存曲线
fit <- survfit(Surv(OS.time, OS) ~ Risk.Level, data = kmdata)

# Cox比例风险模型（计算HR）
fit_cox <- coxph(Surv(OS.time, OS) ~ Risk.Level, data = kmdata)
summary_cox <- summary(fit_cox)

# 提取HR及置信区间（保留两位小数）
HR_conf <- round(summary_cox$conf.int[1, ], 2)
HR_label <- paste0("HR=", HR_conf[1], "(", HR_conf[3], " - ", HR_conf[4], ")")

# 绘制生存曲线
ggsurvplot(
  fit = fit, 
  data = kmdata,
  fun = "pct",              # 纵轴显示生存百分比
  palette = c("#6AAFE6", "#ED9282"),  # 颜色方案（蓝=低风险，橙红=高风险）
  linetype = 1,             # 实线
  pval = TRUE,              # 显示log-rank检验p值
  censor = TRUE,            # 显示删失标记
  censor.size = 7,          # 删失点大小
  risk.table = F,           # 不显示风险表
  conf.int = TRUE           # 显示置信区间带
) + 
  ggtitle(HR_label)          # 添加HR值标题

# 保存PDF
ggsave("./output/2-KM_curve.pdf", height = 5, width = 5.1)

# 5 Univariate Cox分析----
# 定义需要分析的协变量（与图片中变量一致）
covariates <- c("Risk.Score", "Age", "Gender", "Stage_T","Stage_N")

# 生成单变量Cox回归公式
univ_formulas <- sapply(
  covariates,
  function(x) as.formula(paste('Surv(OS.time, OS) ~ ', paste("`", x, "`", 
                                                             sep = "")))
)
# 自定义函数：从Cox模型对象中提取关键结果（HR、置信区间、p值）
cox_extr <- function(fit){ ## 定义个功能，后面输出Cox分析结果用
  fit_summary <- summary(fit)
  dat_res <- data.frame(
    row.names = rownames(fit_summary$coef),
    p.value = signif(fit_summary$coef[,"Pr(>|z|)"]),
    mean = signif(fit_summary$coef[,"exp(coef)"]),
    lower = signif(fit_summary$conf.int[,"lower .95"]),
    upper = signif(fit_summary$conf.int[,"upper .95"]),
    coef = signif(fit_summary$coef[,"coef"]),
    check.rows = F
  )
  dat_res <- signif(dat_res,digits = 3)
  return(dat_res)
}
res_uni <- lapply(univ_formulas, function(x){coxph(x, data = data)}) %>% #进行单因素cox回归
  lapply(cox_extr) %>% #提取单因素cox回归结果
  c(use.names = F) %>% 
  do.call(what = rbind) %>% 
  dplyr::mutate(Factor = rownames(.), HR = paste(.$mean, " [", .$lower, " - ", .$upper, "]", sep = "")) %>%
  dplyr::select("Factor", "p.value", "HR", "mean", "lower", "upper")


# 准备森林图标签数据
res1 <- res_uni

# 格式化 p.value 列
res1$p.value_formatted <- ifelse(res1$p.value < 0.001, "<0.001", sprintf("%.3f", res1$p.value))

# 重新构建 labeltext 矩阵
labeltext <- as.matrix(res1[, c(1, 7, 3)])  # 使用格式化后的 p.value

# 将格式化后的p.value_formatted重命名为p.value
colnames(labeltext)[2] <- "p.value"

# 添加表头行
labeltext <- rbind(colnames(res1)[c(1, 2, 3)], labeltext)  


## 5.1 根据Univariate Cox分析结果画ForestPlot----
pdf("./output/4-ForestPlot_Uni.pdf", height = 5.3, width = 6.1, onefile = F)
forestplot(
  labeltext,                            # 左侧标签文本
  mean = c(NA, res1$mean),             # HR点估计值（首行留空）
  lower = c(NA, res1$lower),           # 置信区间下限
  upper = c(NA, res1$upper),           # 置信区间上限
  zero = 1,                            # 参考线位置（HR=1）
  
  # 图形样式参数（与多变量森林图一致）
  lwd.ci = 2,                          # 置信区间线宽
  ci.vertices = FALSE,                 # 禁用箭头状端点
  title = "Univariable Cox Regression Analysis",  # 标题
  is.summary = c(T, F, F, F, F, F, F, F, F, F, F, F, F, F),  # 首行为表头
  colgap = unit(7, 'mm'),              # 列间距
  xlab = "HR",                         # X轴标签
  boxsize = 0.1,                       # HR点估计方框大小
  clip = c(0, 7),                      # X轴显示范围（0-7，超出显示箭头）
  col = fpColors(                       # 颜色方案（与多变量图一致）
    box = "#EE9A49",                   # 橙色方框（HR点）
    line = "#BEBEBE",                  # 灰色置信区间线
    zero = "#7F8C8D"                   # 深灰色参考线
  ),
  xticks = c(0, 1, 2, 3, 4, 5, 6, 7),# 设定刻度
)
dev.off()

# 6 multicox分析----

# 去掉单因素P值小于0.1的临床变量
data1 <- data %>% 
  dplyr::select(-c(Gender))

##补充：
data1 <- data1 %>%
  filter(!.$Stage_N %in% c("N1"))

data1 <- data1 %>%
  filter(!.$Stage_T %in% c(" T2"," T3"))
##

##
# 构建多变量Cox比例风险模型
mul_cox <- coxph(Surv(OS.time, OS) ~ Risk.Score + Age + Stage_T + Stage_N, data = data1)

# 获取模型摘要
mul_cox1 <- summary(mul_cox)

# 提取置信区间列名
colnames(mul_cox1$conf.int)  

# 提取HR和置信区间（保留2位小数）
multi1 <- as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))

# 生成可发表的格式化表格
multi2 <- tableone::ShowRegTable(mul_cox,
                                 exp = TRUE,         # 输出HR而非原始系数
                                 digits = 2,         # HR保留2位小数
                                 pDigits = 3,        # p值保留3位小数
                                 printToggle = TRUE, # 打印结果
                                 ciFun = confint)    # 使用Wald置信区间

# 合并结果
result <- cbind(multi1, multi2)

# 添加因素名称列
result <- tibble::rownames_to_column(result, var = "Characteristics")
res1 <- result 
# 重命名列（对应森林图纵轴标签）
colnames(res1) <- c("Factors","coef","lower","upper","HR(95% CI)","P value")

# 类型转换（确保数值格式）
res1$coef <- as.numeric(res1$coef)
res1$lower <- as.numeric(res1$lower)
res1$upper <- as.numeric(res1$upper)

# 构建森林图标签文本（对应图中左侧文字）
labeltext <- rbind(colnames(res1)[c(1,6,5)], res1[,c(1,6,5)])

## 6.1 根据multicox分析结果画ForestPlot----
pdf("./output/6-ForestPlot_Multi.pdf", height = 5.3, width = 6.1, onefile = F)
forestplot(
  labeltext,               # 左侧标签文本（因素名称、P值、HR值）
  mean = c(NA, res1$coef), # HR点估计值（首行留空）
  lower = c(NA, res1$lower), # 置信区间下限
  upper = c(NA, res1$upper), # 置信区间上限
  zero = 1,                # 参考线位置（HR=1）
  
  # 图形样式参数
  lwd.zero = 2,           # 参考线粗细
  lwd.ci = 2,             # 置信区间线粗细
  title = "Multivariable Cox Regression Analysis", # 标题（与图片一致）
  ci.vertices = FALSE,    # 不使用箭头状置信区间端点
  ci.vertices.height = 0.2, # 端点高度
  is.summary = c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F), # 首行为汇总行
  
  # 布局参数
  colgap = unit(7,'mm'),  # 列间距
  xlab = "HR",            # X轴标签
  boxsize = 0.12,         # HR点估计方框大小
  
  # 颜色参数（与图片配色一致）
  col = fpColors(
    box = "#EE9A49",      # 箱线图橙色（对应图片中的HR点）
    line = "#BEBEBE",     # 置信区间灰色线
    zero = "#7F8C8D"      # 参考线深灰色
  ),
  
  # 坐标轴参数
  clip = c(0,7),          # X轴显示范围0-7（超出显示箭头，与图片横轴一致）
  xticks = c(0, 1, 2, 3, 4, 5, 6, 7),# 设定刻度
  lwd.xaxis = 1           # X轴线宽
)

dev.off()

# 7 nomgram----
# 定义数据分布（使用rms包中的datadist函数，存储变量范围、分位数等信息，便于后续模型预测和绘图）
dd <- datadist(data1) 

# 拟合Cox比例风险模型（使用rms包中的cph函数）
coxm <- cph(
  Surv(OS.time, OS) ~ Risk.Score + Age + Stage_T + Stage_N,  # 模型公式：生存时间与状态 ~ 自变量
  x = T,                                   # 保存模型矩阵（用于预测、验证和计算残差）
  y = T,                                   # 保存响应变量（生存时间与状态）
  data = data1,                            # 使用的数据集
  surv = T                                 # 保存基线生存曲线估计（可通过survest()提取）
)

# 设置全局数据分布选项（确保后续函数如nomogram()能识别变量范围）
options(datadist = "dd")

# 查看模型摘要（输出回归系数、HR、置信区间、p值等）
summary(coxm)

# 提取生存函数对象（Survival()函数生成生存概率计算工具）
surv <- Survival(coxm)

# 定义不同时间点的生存概率函数（lp为线性预测值，即模型中的βX）
# 1年生存概率函数
surv1 <- function(x) surv(365 * 1, lp = x)  
# 3年生存概率函数
surv2 <- function(x) surv(365 * 3, lp = x)  
# 5年生存概率函数
surv3 <- function(x) surv(365 * 5, lp = x)  

# 7. 绘制列线图（nomogram）
nom <- nomogram(
  coxm,                                # 输入的Cox模型对象
  fun = list(surv1, surv2, surv3),     # 需要展示的生存概率函数列表
  funlabel = c(
    '1-year probability',              # 第一个生存概率标签
    '3-year probability',              # 第二个生存概率标签
    '5-year probability'               # 第三个生存概率标签
  ),      
  lp = F,                              # 是否在图中显示线性预测值（lp=F表示不显示）
  fun.at = list(
    c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),  # 1-year 刻度
    c(0.3, 0.5, 0.7, 0.9),             # 3-year 刻度
    c(0.1, 0.3, 0.5, 0.7, 0.9)         # 5-year 刻度
  )
  # fun.at = c('0.9', '0.8', '0.7', '0.6', '0.5', '0.3', '0.1') # 设置生存概率刻度点
)

## 7.1 画nomgram图----
pdf(file = "./output/7-Nomogram.pdf",height = 8.2 ,width = 15.6453 )
par(mar=c(2,5,3,2),cex=1.3)# mar 图形空白边界(下，左，上，右)  cex 文本和符号大小
plot(nom,xfrac=0.6)        # 调整变量名称的右对齐位置（0.6表示占总宽度的60%）
dev.off()

# 8 绘制calibration----
## 8.1 Cali_1yr----
# 构建Cox比例风险模型
geocox <- cph(
  Surv(OS.time, OS) ~ Risk.Score + Age + Stage_T + Stage_N,  # 模型公式：生存时间与状态 ~ 风险评分 + 分期
  data = data1,                            # 使用的数据集
  surv = T,                                # 保存基线生存函数（用于后续预测）
  x = T,                                   # 保存模型矩阵（用于验证和诊断）
  y = T,                                   # 保存响应变量（生存时间与状态）
  time.inc = 365 * 1                       # 指定时间间隔为1年（用于生存函数计算）
)

# 模型校准（Bootstrap重抽样）
cal_1 <- calibrate(
  geocox,                                   # 输入的Cox模型对象
  u = 365 * 1,                              # 校准时间点（1年）
  cmethod = 'KM',                           # 使用Kaplan-Meier方法计算实际生存率
  m = 50,                                   # 将预测概率分为50个区间（m越小，区间越多）
  B = 1000                                  # Bootstrap重抽样次数（B越大，结果越稳定）
)

# 绘制校准曲线并保存为PDF 
pdf(
  file = "./output/8-Cali_1yr.pdf",          # 输出文件路径（需确保目录存在）
  height = 8,                               # 图形高度（英寸）
  width = 8                                  # 图形宽度（英寸，正方形布局）
)

# 设置绘图参数
par(
  mar = c(7, 4, 4, 3),                      # 图形边距：下7行，左4行，上4行，右3行
  cex = 1.0                                  # 全局文本缩放因子（1.0为原始大小）
)

# 绘制校准曲线
plot(
  cal_1,                                     # 校准结果对象
  lwd = 2,                                   # 线条宽度（加粗主对角线）
  lty = 1,                                   # 线条类型（1=实线）
  errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),  # 误差线颜色（蓝色，#0076C0）
  xlab = 'Nomogram-Predicted Probability of 1-year OS',    # x轴标签
  ylab = 'Actual 1-year OS (proportion)',                  # y轴标签
  col = c(rgb(192, 98, 83, maxColorValue = 255)),          # 校准点颜色（红色，#C06253）
  xlim = c(0, 1),                            # x轴范围（预测概率0~100%）
  ylim = c(0, 1),                           # y轴范围（实际概率0~100%）
  subtitles = FALSE                         # 禁用所有统计信息
)  


dev.off()

## 8.2 Cali_3yr----
# 构建Cox比例风险模型
geocox_2 <- cph(
  Surv(OS.time, OS) ~ Risk.Score + Age + Stage_T + Stage_N,  # 模型公式：新增年龄和性别变量
  data = data1,                            # 使用的数据集
  surv = T,                               # 保存基线生存函数（用于后续预测）
  x = T,                                  # 保存模型矩阵（用于验证和诊断）
  y = T,                                  # 保存响应变量（生存时间与状态）
  time.inc = 365 * 3                      # 指定时间间隔为3年（用于生存函数计算）
)

# 模型校准（Bootstrap重抽样）
cal_2 <- calibrate(
  geocox_2,                                # 输入的Cox模型对象（含年龄和性别）
  u = 365 * 3,                             # 校准时间点（3年）
  cmethod = 'KM',                          # 使用Kaplan-Meier方法计算实际生存率
  m = 50,                                  # 将预测概率分为50个区间（m越小，区间越精细）
  B = 1000                                 # Bootstrap重抽样次数（增加结果稳定性）
)

# 绘制校准曲线并保存为PDF
pdf(
  file = "./output/9-Cali_3yr.pdf",         # 输出文件路径（需确保目录存在）
  height = 8,                              # 图形高度（英寸）
  width = 8                                # 图形宽度（英寸，正方形布局）
)

# 设置绘图参数
par(
  mar = c(7, 4, 4, 3),                      # 图形边距：下7行，左4行，上4行，右3行
  cex = 1.0                                 # 全局文本缩放因子（1.0为原始大小）
)

# 绘制校准曲线
plot(
  cal_2,                                    # 校准结果对象（3年生存率）
  lwd = 2,                                  # 线条宽度（加粗主对角线）
  lty = 1,                                  # 线条类型（1=实线）
  errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),  # 误差线颜色（蓝色，#0076C0）
  xlab = 'Nomogram-Predicted Probability of 3-year OS',    # x轴标签
  ylab = 'Actual 3-year OS (proportion)',                  # y轴标签
  col = c(rgb(192, 98, 83, maxColorValue = 255)),          # 校准点颜色（红色，#C06253）
  xlim = c(0, 1),                           # x轴范围（预测概率0~100%）
  ylim = c(0, 1),                           # y轴范围（实际概率0~100%）
  subtitles = FALSE                         # 禁用所有统计信息
)  


dev.off()
## 8.3 Cali_5yr----
# 构建Cox比例风险模型
geocox_3 <- cph(
  Surv(OS.time, OS) ~ Risk.Score + Age + Stage_T + Stage_N,  # 模型公式：含风险评分、年龄、性别、分期
  data = data1,                            # 使用的数据集
  surv = T,                               # 保存基线生存函数（用于后续预测）
  x = T,                                  # 保存模型矩阵（用于验证和诊断）
  y = T,                                  # 保存响应变量（生存时间与状态）
  time.inc = 365 * 5                      # 指定时间间隔为5年（用于生存函数计算）
)

# 模型校准（Bootstrap重抽样）
cal_3 <- calibrate(
  geocox_3,                                # 输入的Cox模型对象
  u = 365 * 5,                             # 校准时间点（5年）
  cmethod = 'KM',                          # 使用Kaplan-Meier方法计算实际生存率
  m = 50,                                  # 将预测概率分为50个区间（m越小，区间越精细）
  B = 1000                                 # Bootstrap重抽样次数（增加结果稳定性）
)

# 绘制校准曲线并保存为PDF 
pdf(
  file = "./output/10-Cali_5yr.pdf",         # 输出文件路径（需确保目录存在）
  height = 8,                              # 图形高度（英寸）
  width = 8                                # 图形宽度（英寸，正方形布局）
)

# 设置绘图参数
par(
  mar = c(7, 4, 4, 3),                     # 图形边距：下7行，左4行，上4行，右3行
  cex = 1.0                                 # 全局文本缩放因子（1.0为原始大小）
)

# 绘制校准曲线
plot(
  cal_3,                                    # 校准结果对象（5年生存率）
  lwd = 2,                                  # 线条宽度（加粗主对角线）
  lty = 1,                                  # 线条类型（1=实线）
  errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),  # 误差线颜色（蓝色，#0076C0）
  xlab = 'Nomogram-Predicted Probability of 5-year OS',    # x轴标签
  ylab = 'Actual 5-year OS (proportion)',                  # y轴标签
  col = c(rgb(192, 98, 83, maxColorValue = 255)),          # 校准点颜色（红色，#C06253）
  xlim = c(0, 1),                           # x轴范围（预测概率0~100%）
  ylim = c(0, 1),                             # y轴范围（实际概率0~100%）
  subtitles = FALSE                         # 禁用所有统计信息
)  

dev.off()

# 9 绘制DCA----
mycol <- c('#6679c9','#aacd89','#f4cf72',
           '#df7971','#94c6df','#60a980',
           '#Ef9366','#9d6eba','#e28cd0')
# 拟合Cox模型（用于DCA分析）
f <- coxph(Surv(OS.time, OS) ~ Risk.Score + Age + Stage_T + Stage_N, data = data1)

# DCA分析
model <- dca(
  f,                            # 输入的Cox模型对象
  model.names = "Model"         # 模型名称（显示在图例中）
)

# 绘制DCA曲线（未指定时间）
ggplot(model, lwd = 1, color = mycol)  # lwd: 线条宽度; color: 使用预定义颜色

# 1年DCA分析
model_1 <- dca(
  f,
  model.names = "1-year",       # 模型名称（图例显示为 "1-year"）
  times = 365 * 1               # 时间点：1年（单位与OS.time一致，如天）
)

# 绘制并调整x轴范围（0-0.6，聚焦常用阈值区间）
ggplot(model_1, lwd = 1, color = mycol) + 
  scale_x_continuous(limits = c(0, 0.2))  # 限制x轴显示范围


# 保存为PDF（需确保路径存在）
ggsave(
  "./output/11-DCA_1yr.pdf",    # 输出路径
  height = 8,                   # 图形高度（英寸）
  width = 8                     # 图形宽度（英寸）
)

# 3年DCA分析
model_2 <- dca(
  f,
  model.names = "3-year",       # 模型名称（图例显示为 "3-year"）
  times = 365 * 3               # 时间点：3年
)

# 绘制并调整x轴范围（0-0.9）
ggplot(model_2, lwd = 1, color = mycol) + 
  scale_x_continuous(limits = c(0, 0.9))

# 保存为PDF
ggsave(
  "./output/12-DCA_3yr.pdf",
  height = 8,
  width = 8
)

# 5年DCA分析
model_3 <- dca(
  f,
  model.names = "5-year",       # 模型名称（图例显示为 "5-year"）
  times = 365 * 5               # 时间点：5年
)

# 绘制并调整x轴范围（0-1）
ggplot(model_3, lwd = 1, color = mycol) + 
  scale_x_continuous(limits = c(0, 0.5))

# 保存为PDF
ggsave(
  "./output/13-DCA_5yr.pdf",
  height = 8,
  width = 8
)

# 10 输出结果保存----
## 10.1 保存csv文件----
write.csv(res_uni,"./output/3-U_cox.csv")
write.csv(result,"output/5-M_cox.csv")

## 10.2 保存工作环境----
save.image("1-PrognosticAnalysis.RData.gz",compress = "gzip")

