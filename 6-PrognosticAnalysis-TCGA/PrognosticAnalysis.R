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
riskscores <- data.table::fread("input/7-RiskScore_LASSO.csv", data.table = F)

# 设置行名为第一列
rownames(riskscores) <- riskscores[, 1]

# 删除第一列（已将其设置为行名，避免重复）
riskscores <- riskscores[, -1]

## 2.2 读取并处理临床/表型数据 ----
# 从CSV文件读取临床数据（假设文件路径为`input/PD.csv`）
pd <- data.table::fread("input/PD.csv", data.table = F)

# 重命名第一列为"RNAseq"（假设该列为样本ID，与风险评分数据中的行名对应）
colnames(pd)[1] <- "RNAseq"

# 筛选特定亚型（例如肺腺癌LUAD）
pd <- pd[pd$Group == "STAD", ]

# 筛选存在于风险评分数据中的样本（确保数据一致性）
pd <- pd[pd$RNAseq %in% rownames(riskscores), ]

# 按风险评分数据的行名顺序重新排列临床数据（确保样本顺序完全一致）
pd <- pd[match(rownames(riskscores), pd$RNAseq), ]


## 2.3 对临床信息进行处理----
#对数值型变量分类
#查看样本是否对齐
identical(riskscores$sampleid,pd$RNAseq)

# 修改：选取相应临床信息列
data <- cbind(pd[,c(4:5,8:9,11:13)],riskscores[,1])
data[data == ""] <- NA
data <- na.omit(data)
colnames(data)[8] <- c("Risk.Score")
table(data$Age)
table(data$Gender)
table(data$Stage_M)
table(data$Stage_N)
table(data$Stage_T)

#规范临床分期名称
data[data$Age <= 60,]$Age <- "<= 60"
data[data$Age > 60,]$Age <- "> 60"

dir.create("output")

# 3 time_ROC分析----
## 3.1 计算时间依赖性AUC ----
library(timeROC)  # 加载timeROC包

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
# 读取风险评分文件（假设路径为input/7-RiskScore_LASSO.csv）
riskstatus <- data.table::fread("input/7-RiskScore_LASSO.csv", data.table = F)

# 按样本ID对齐风险评分数据（假设riskstatus$V1为样本ID）
riskscores <- riskscores[riskstatus$V1, ]

# 将风险等级列合并到风险评分数据中
riskscores <- cbind(riskscores, riskstatus$risklevel)

# 设置临床数据（pd）的行名为样本ID（假设第一列为RNAseq ID）
rownames(pd) <- pd[, 1]

# 确保临床数据与风险评分数据的样本顺序一致
pd <- pd[rownames(riskscores), ]

# 验证样本ID一致性（应返回TRUE）
identical(pd$RNAseq,rownames(riskscores))

# 构建KM分析数据集（提取OS时间、OS状态、风险等级）
kmdata <- cbind(pd[, 4:5], riskscores[, 2])  # 假设第4-5列为OS.time和OS
kmdata <- na.omit(kmdata)                   # 删除缺失值
colnames(kmdata) <- c("OS", "OS.time", "Risk.Level")  # 重命名列

## 4.1 画KM曲线----
head(kmdata)
# 按风险等级降序排列（高风险在前）
kmdata <- kmdata[order(kmdata$Risk.Level, decreasing = T), ]

# 将风险等级转为因子（确保绘图顺序：LowRisk在前，HighRisk在后）
kmdata$Risk.Level <- factor(kmdata$Risk.Level, levels = c("LowRisk", "HighRisk"))

# 拟合Kaplan-Meier生存曲线
fit <- survfit(Surv(OS.time, OS) ~ Risk.Level, data = kmdata)

# 输出生存曲线摘要（中位生存时间等）
print(fit)

# 拟合Cox比例风险模型（计算HR）
fit_cox <- coxph(Surv(OS.time, OS) ~ Risk.Level, data = kmdata)
summary_cox <- summary(fit_cox)

# 提取HR和置信区间（保留2位小数）
HR_conf <- round(summary_cox$conf.int[1, ], 2)
HR_label <- paste0("HR=", HR_conf[1], "(", HR_conf[3], " - ", HR_conf[4], ")")

# 验证比例风险假设（p>0.05表示满足假设）
cox.zph(fit_cox)

# 使用ggsurvplot绘制生存曲线（需加载survminer包）
ggsurvplot(
  fit = fit,                 # 生存曲线对象
  data = kmdata,             # 数据集
  fun = "pct",               # 纵轴显示生存百分比
  palette = c("#6AAFE6", "#ED9282"),  # 颜色方案（蓝色=低风险，橙红=高风险）
  linetype = 1,              # 实线
  pval = TRUE,               # 显示log-rank检验p值
  censor = TRUE,             # 显示删失标记
  censor.size = 7,          # 删失点大小
  risk.table = F,            # 不显示风险表
  conf.int = TRUE            # 显示置信区间带
) + 
  ggtitle(HR_label)          # 添加标题（HR值）

# 保存高清PDF
ggsave("./output/2-KM_curve.pdf", height = 5, width = 5.1)

# 5 Univariate Cox分析----
covariates <- c("Risk.Score", "Age","Gender", "Stage_M", "Stage_N", "Stage_T")

# 生成单变量Cox公式
univ_formulas <- sapply(
  covariates,
  function(x) as.formula(paste('Surv(OS.time, OS) ~ ', paste0("`", x, "`")))
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
# 设置输出为PDF文件，指定路径、尺寸（高度8英寸，宽度8.1英寸），且不合并多页
pdf("./output/4-ForestPlot_Uni.pdf", height = 8, width = 8.1, onefile = F)

# 绘制森林图
forestplot(
  labeltext,                      # 左侧标签文本矩阵（包含变量名、p值、HR+置信区间）
  mean = c(NA, res1$mean),        # HR点估计值（首行NA跳过表头）
  lower = c(NA, res1$lower),      # 置信区间下限
  upper = c(NA, res1$upper),      # 置信区间上限
  zero = 1,                       # 参考线位置（HR=1，灰色垂直线）
  lwd.ci = 2,                     # 置信区间线宽（加粗）
  ci.vertices = FALSE,            # 禁用置信区间端点箭头（使用直线端点）
  title = "Univariable Cox Regression Analysis",  # 图表标题
  is.summary = c(T, F, F, F, F, F, F, F, F, F, F, F), 
  colgap = unit(7, 'mm'),         # 列间距（7毫米，优化标签对齐）
  xlab = "HR",                    # X轴标签（Hazard Ratio）
  align = NULL,                   # 不强制对齐列（保持自然间距）
  txt_gp = fpTxtGp(
    ticks = gpar(cex = 0.8),      # 刻度标签字号
    xlab = gpar(cex = 0.8),       # X轴标签字号（0.8倍默认大小）
    title = gpar(cex = 1.2),      # 标题字号（1.2倍默认大小）
    cex = 1.2                     # 全局文本缩放（可能覆盖个别设置）
  ),
  clip = c(0, 10),                 # X轴显示范围（0-10，超出部分显示为箭头）
  boxsize = 0.1,                  # HR点估计方框大小（值越小，方框越紧凑）
  col = fpColors(
    box = "#EE9A49",             # HR点估计方框颜色（橙色）
    line = "#BEBEBE",            # 置信区间线颜色（灰色）
    zero = "#7F8C8D"             # 参考线（HR=1）颜色（深灰色）
  ),
  # xticks = NULL,                  # 不显示默认X轴刻度（自动生成刻度）
  xticks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),# 设定刻度
  lwd.xaxis = 1                   # X轴线宽（1磅）
)

dev.off()

# 6 multicox分析----
data1 <- data
# 保留单因素P值小于0.1的临床变量
# data1 <- data %>% 
#   dplyr::select(-c(Gender))
# 
# #补充：
# data1 <- data1 %>%
#   filter(!.$Stage_N %in% c("N1","N2"))
# 
# data1 <- data1 %>%
#   filter(!.$Stage_T %in% c("T2"))
# 1. 拟合多变量Cox比例风险模型
mul_cox <- coxph(Surv(OS.time, OS) ~ Risk.Score + Age + Gender + Stage_M + Stage_N + Stage_T, data = data1)

# 2. 生成模型摘要
mul_cox1 <- summary(mul_cox)

# 3. 查看置信区间列名
colnames(mul_cox1$conf.int)

# 4. 提取关键结果：HR(exp(coef))、95%置信区间上下限，保留2位小数
multi1 <- as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))

# 5. 生成可发表的格式化表格（tableone包功能）
multi2 <- tableone::ShowRegTable(
  mul_cox,              # Cox回归模型对象
  exp = TRUE,           # 显示风险比（HR）而非原始系数β
  digits = 2,          # HR和置信区间保留2位小数
  pDigits = 3,         # p值保留3位小数
  printToggle = TRUE,  # 打印结果到控制台
  quote = FALSE,       # 输出不添加引号
  ciFun = confint      # 使用confint函数计算置信区间（默认Wald法）
)

# 6. 合并数值型结果（multi1）和格式化表格（multi2）
result <- cbind(multi1, multi2)


result <- tibble::rownames_to_column(result, var = "Characteristics")
res1 <- result
colnames(res1) <- c("Factors","coef","lower","upper","HR(95% CI)","P value")
res1$coef <- as.numeric(res1$coef)
res1$lower <- as.numeric(res1$lower)
res1$upper <- as.numeric(res1$upper)
labeltext <- rbind(colnames(res1)[c(1,6,5)],res1[,c(1,6,5)])

## 6.1 根据multicox分析结果画ForestPlot----
# 设置输出为PDF文件
pdf("./output/6-ForestPlot_Multi.pdf", height = 8, width = 8.1, onefile = F)

# 绘制森林图
forestplot(
  labeltext,                        # 左侧标签文本矩阵（变量名、p值、HR+置信区间）
  mean = c(NA, res1$coef),         # HR点估计值（首行NA跳过表头）
  lower = c(NA, res1$lower),       # 置信区间下限
  upper = c(NA, res1$upper),       # 置信区间上限
  zero = 1,                        # 参考线位置（HR=1）
  lwd.zero = 2,                    # 参考线线宽（2磅，加粗）
  lwd.ci = 2,                      # 置信区间线宽（2磅）
  title = "Multivariable Cox Regression Analysis",  # 图表标题
  ci.vertices = FALSE,             # 禁用置信区间端点箭头（使用直线端点）
  ci.vertices.height = 0.2,        # 端点高度（若启用箭头时控制高度）
  is.summary = c(T, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F),  # 首行（表头）标记为汇总行（加粗）
  colgap = unit(7, 'mm'),           # 列间距（7毫米，优化对齐）
  xlab = "HR",                     # X轴标签（Hazard Ratio）
  align = NULL,                    # 不强制对齐列（自然间距）
  txt_gp = fpTxtGp(                # 文本样式设置
    ticks = gpar(cex = 0.8),       # 刻度标签字号
    xlab = gpar(cex = 0.8),        # X轴标签字号（0.8倍默认大小）
    title = gpar(cex = 1.2),       # 标题字号（1.2倍默认大小）
    cex = 1.2                      # 全局文本缩放（可能覆盖个别设置）
  ),
  boxsize = 0.12,                  # HR点估计方框大小（0.12英寸，控制点大小）
  clip = c(0, 7),                  # X轴显示范围（HR=0-7，超出部分显示为箭头）
  col = fpColors(                  # 颜色方案
    box = "#EE9A49",              # HR点估计方框颜色（橙色）
    line = "#BEBEBE",             # 置信区间线颜色（灰色）
    zero = "#7F8C8D"            # 参考线（HR=1）颜色（深灰色）
  ),
  xticks = c(0, 1, 2, 3, 4, 5, 6, 7),# 设定刻度
  lwd.xaxis = 1                   # X轴线宽（1磅）
)

dev.off()

# 7 nomgram----
# 定义数据分布（用于后续模型预测和绘图）
dd <- datadist(data1)

# 构建Cox比例风险模型（使用rms包的cph函数）
coxm <- cph(Surv(OS.time, OS) ~ Risk.Score + Age + Gender + Stage_M + Stage_N + Stage_T,
            x = T, y = T,  # 保存模型矩阵和响应变量
            data = data1, surv = T)  # surv=T保存基线生存曲线

# 设置全局数据分布选项（确保后续函数识别变量范围）
options(datadist = "dd")

# 查看模型摘要（输出系数、HR、置信区间等）
summary(coxm)

# 提取生存概率计算工具
surv <- Survival(coxm)

# 定义1年、3年、5年生存概率函数（lp为线性预测值）
surv1 <- function(x) surv(365 * 1, lp = x)  # 1年生存概率
surv2 <- function(x) surv(365 * 3, lp = x)  # 3年生存概率
surv3 <- function(x) surv(365 * 5, lp = x)  # 5年生存概率

# 创建列线图对象
nom <- nomogram(
  coxm,  # Cox模型
  fun = list(surv1, surv2, surv3),  # 需要展示的生存概率函数
  funlabel = c('1-year probability', '3-year probability', '5-year probability'),  # 生存概率标签
  lp = F,  # 不显示线性预测值（LP）
  fun.at = c('0.9', '0.8', '0.7', '0.6', '0.5', '0.3', '0.1')  # 生存概率刻度点
)

# 绘制并保存列线图
pdf(file = "./output/7-Nomogram.pdf", height = 8.2, width = 15.6453)
par(mar = c(2, 5, 3, 2), cex = 1.3)  # 设置边距和文本大小
plot(nom, xfrac = 0.6)  # xfrac控制右侧文本对齐位置
dev.off()

# 8 绘制calibration----
## 8.1 Cali_1yr----
# 重新拟合Cox模型（指定时间间隔为1年）
geocox <- cph(Surv(OS.time, OS) ~ Risk.Score + Age + Gender + Stage_M + Stage_N + Stage_T,
              data = data1, surv = T, x = T, y = T, time.inc = 365 * 1)

# 计算校准曲线（Bootstrap重抽样评估模型校准度）
cal_1 <- calibrate(geocox,
                   u = 365 * 1,       # 校准时间点（1年）
                   cmethod = 'KM',    # 使用Kaplan-Meier方法计算实际生存率
                   m = 50,           # 将预测概率分为50个区间
                   B = 1000)         # Bootstrap重抽样次数

# 绘制校准曲线
pdf(file = "./output/8-Cali_1yr.pdf", height = 8, width = 8)
par(mar = c(7, 4, 4, 3), cex = 1.0)  # 增大底部边距防止标签被截断
plot(cal_1,
     lwd = 2, lty = 1,              # 线宽和线型
     errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),  # 误差线颜色（蓝色）
     xlab = 'Nomogram-Predicted Probability of 1-year OS',  # X轴标签
     ylab = 'Actual 1-year OS (proportion)',                # Y轴标签
     col = c(rgb(192, 98, 83, maxColorValue = 255)),        # 校准点颜色（红色）
     xlim = c(0, 1), ylim = c(0, 1),                        # 坐标轴范围
     subtitles = FALSE)  # 不显示副标题
dev.off()


## 8.2 Cali_3yr----
geocox_2 <- cph(Surv(OS.time,OS)~Risk.Score + Age + Gender +Stage_M+Stage_N+Stage_T,
                data=data1,surv=T,x=T,y=T,time.inc = 365*3)
cal_2 <- calibrate(geocox_2,u=365*3,cmethod='KM',m=50,B=1000)
pdf(file = "./output/9-Cali_3yr.pdf",height = 8 ,width = 8 )
par(mar=c(7,4,4,3),cex=1.0)
plot(cal_2,lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 3-year OS',#便签
     ylab='Actual 3-year OS(proportion)',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)),#设置一个颜色
     xlim = c(0,1),ylim = c(0,1), ##x轴和y轴范围
     subtitles = FALSE)  # 禁用所有统计信息
dev.off()

## 8.3 Cali_5yr----
geocox_3 <- cph(Surv(OS.time,OS)~Risk.Score + Age + Gender +Stage_M+Stage_N+Stage_T,
                data=data1,surv=T,x=T,y=T,time.inc = 365*5)

cal_3 <- calibrate(geocox_3,u=365*5,cmethod='KM',m=50,B=1000)
pdf(file = "./output/10-Cali_5yr.pdf",height = 8 ,width = 8 )
par(mar=c(7,4,4,3),cex=1.0)
plot(cal_3,lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 5-year OS',#便签
     ylab='Actual 5-year OS(proportion)',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)),#设置一个颜色
     xlim = c(0,1),ylim = c(0,1), ##x轴和y轴范围
     subtitles = FALSE)  # 禁用所有统计信息
dev.off()

# 9 绘制DCA----
mycol <- c('#6679c9','#aacd89','#f4cf72',
           '#df7971','#94c6df','#60a980',
           '#Ef9366','#9d6eba','#e28cd0')

# 拟合Cox模型（用于DCA分析）
f <- coxph(Surv(OS.time,OS)~Risk.Score + Age + Gender +Stage_M+Stage_N+Stage_T,data=data1)

# DCA分析
model<-dca(f,model.names ="Model")
ggplot(model, lwd = 1,color = mycol)

## 9.1 DCA_1yr----
model_1<-dca(f,model.names ="1-year",times = 365*1)
ggplot(model_1, lwd = 1,color = mycol) + scale_x_continuous (limits = c (0, 0.6))
ggsave("./output/11-DCA_1yr.pdf",height = 8 ,width = 8 )

## 9.2 DCA_3yr----
model_2<-dca(f,model.names ="3-year",times = 365*3)
ggplot(model_2, lwd = 1,color = mycol) + scale_x_continuous (limits = c (0, 0.9))
ggsave("./output/12-DCA_3yr.pdf",height = 8,width = 8)

## 9.3 DCA_5yr----
model_3<-dca(f,model.names ="5-year",times = 365*5)
ggplot(model_3, lwd = 1,color = mycol)+ scale_x_continuous (limits = c (0, 1))
ggsave("./output/13-DCA_5yr.pdf",height = 8,width = 8)


# 10 输出结果保存----
## 10.1 保存csv文件----
write.csv(res_uni,"./output/3-U_cox.csv")
write.csv(result,"output/5-M_cox.csv")

## 10.2 保存工作环境----
save.image("1-PrognosticAnalysis.RData.gz",compress = "gzip")

