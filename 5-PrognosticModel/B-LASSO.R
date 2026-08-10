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
required_packages <- c("glmnet", "survival", "magrittr")

#调用函数进行检查和安装
install_if_missing(required_packages)

## 1.3 加载R包----
for(i in c("glmnet", "survival", "magrittr")){
  library(i, character.only = T)
}

# 2 Lasso分析----

#' lasso分析
#' 
#' @param genes(characters) 选择要筛选的基因
#' @param group(factor or Surv) 选择分组信息，诊断模型是二分类变量，预后模型是生存信息（包括time和event）event需要转变为0和1
#'   \cr 诊断模型默认methods为"binomial"
#'   \cr 预后模型默认methods为"cox"
#' @param data(data.frame) 包含基因和分组信息的数据
#' @param lambda_choose(characters) 'lambda.1se', 'lambda.min'。lambda.min下的模型交叉验证误差最小。lambda.1se给出了最正则化的模型，使得交叉验证误差在最小值的一个标准误差内。
#' @param nfold(numeric) K-fold交叉验证，数据量少可以选5
#' @param iter.time(numeric) 循环次数，默认不进行循环
#' @param sed(numeric) 种子数，默认是500
#' @param save(logic) 是否保存
#' @param path(characters) 路径
#' @param height(numeric) 保存图片的高度
#' @param width(numeric) 保存图片的宽度
#' 
#' @return (data.frame) lasso筛选的基因和系数
#' 
#' @export

## 2.1 定义lasso分析函数----

lasso_analysis <- function(genes, group, data, 
                           lambda_choose = c('lambda.1se','lambda.min'), 
                           nfold=10, iter.times=NULL, sed = 500,
                           save=T, path="./", height= 5, width = 5){
  #cv.glmnet输入变量：x, y
  data <- na.omit(data)
  mat_x <- as.matrix(data[, genes])
  y <- data[, group]
  lasso_methods <- "binomial"
  #group如果为两个：预后lasso
  if(length(group) > 1){
    if(length(group) == 2){
      #判断group变量
      len <- apply(y, 2, function(x) length(unique(x)))
      event_col <- which(len == 2)
      time_col <- which(len != 2)
      if(identical(event_col, integer(0))){
        stop("No event variables were found\n")
      } else{
        if(min(y[, time_col]) > 0){
          cat(paste0("use ", colnames(y)[event_col], " as event variable\n"),
              paste0("use ", colnames(y)[time_col], " as time variable\n"))
          y <- data.matrix(Surv(time = as.double(y[, time_col]),
                                event = as.double(y[, event_col])))
          lasso_methods <- "cox"
        } else {
          cat("Delete time <= 0\n",
              paste0("use ", colnames(y)[event_col], " as event variable\n"),
              paste0("use ", colnames(y)[time_col], " as time variable\n"))
          data <- data[-which(y[, time_col] <= 0),]
          mat_x <- as.matrix(data[, genes])
          y <- data[, group]
          y <- data.matrix(Surv(time = as.double(y[, time_col]),
                                event = as.double(y[, event_col])))
          lasso_methods <- "cox"
        }
      }
    }else{
      stop("The length of group must be 1 or 2\n")
    }
  }
  lambda_choose <- match.arg(lambda_choose)
  fit <- glmnet(mat_x, y, family = lasso_methods)
  #是否循环cv.glmnet
  if(is.null(iter.times)){
    set.seed(sed)
    cvfit <- cv.glmnet(mat_x, y, family = lasso_methods, nfolds = nfold)
    lasso_fit <- glmnet(mat_x, y, family = lasso_methods,
                        lambda = cvfit[[lambda_choose]])
    fea <- rownames(coef(lasso_fit))[coef(lasso_fit)[, 1]!= 0]
    if(is.element("(Intercept)", fea)) {
      fea <- fea[-1] # 去掉截距项
      lasso_coef <- coef(lasso_fit)@x[-1]
    } else {
      fea <- fea
      lasso_coef <- coef(lasso_fit)@x
    }
    dat_res <- data.frame(genes = fea,
                          coef = signif(lasso_coef, digits = 3),
                          seed = sed,
                          stringsAsFactors = F)
  }else{
    #cv.glmnet循环iter.times次，选择cv.glmnet结果的最优解
    lasso_fea_list <- lapply(seq_len(iter.times), function(x){ # 
      set.seed(x)
      cvfit <- cv.glmnet(mat_x, y, family = lasso_methods)
      # 取出最优lambda
      fea <- rownames(coef(cvfit, 
                           s = lambda_choose))[coef(cvfit, s = lambda_choose)[, 1]!= 0]
      if(is.element("(Intercept)", fea)) {
        fea <- fea[-1] # 去掉截距项并排序
      } else {
        fea <- fea
      }
      lasso_res <- data.frame(iteration = x,
                              genes = paste0(fea, collapse = " and "),
                              stringsAsFactors = F)
      return(lasso_res)
    })
    lasso_res <- do.call(rbind, lasso_fea_list)
    sel.iter <- lasso_res[which(lasso_res$genes == names(table(lasso_res$genes))[1]),
                          "iteration"][1]
    #根据最优解重新得到lasso结果并绘图
    set.seed(sel.iter)
    cvfit <- cv.glmnet(mat_x, y, family = lasso_methods, nfolds = nfold)
    lasso_fit <- glmnet(mat_x, y, family = lasso_methods, lambda = cvfit[[lambda_choose]])
    fea <- rownames(coef(lasso_fit))[coef(lasso_fit)[, 1]!= 0]
    if(is.element("(Intercept)", fea)) {
      fea <- fea[-1] # 去掉截距项
      lasso_coef <- coef(lasso_fit)@x[-1]
    } else {
      fea <- fea
      lasso_coef <- coef(lasso_fit)@x
    }
    dat_res <- data.frame(genes = fea,
                          coef = signif(lasso_coef, digits = 3),
                          seed = sel.iter,
                          stringsAsFactors = F)
  }
  #绘图
  par(mfrow = c(1, 2))
  plot(fit)
  plot(cvfit)
  if(save) {
    if(is.null(path)){
      stop("No path")
    }else{
      pdf(paste0(path, "/", "5-LASSO_Track.pdf"), height = height, width = width)
      plot(fit, xvar = "dev", label = T)
      dev.off()
      pdf(paste0(path, "/", "6-LASSO_Coef.pdf"), height = height, width = width)
      plot(cvfit)
      dev.off()
    }
  }
  return(dat_res)
}


## 2.2 载入cox结果数据----
dat_res=readRDS("./output/3-cox_res.rds")
data=readRDS("./output/4-cox_data.rds")

## 2.3 进行lasso分析----
prognosis_lasso <- lasso_analysis(rownames(dat_res), group = c("OS.time", "OS"),
                                  data, lambda_choose = 'lambda.min', 
                                  nfold = 10, iter.times= 10, sed = 500, save = T, 
                                  path = "./output", height = 5, width = 8.1)
dev.off()
# Error in data.frame(genes = fea, coef = signif(lasso_coef, digits = 3), : 
# 参数值意味着不同的行数: 0, 1
# 循环过程中，由于LASSO筛选不到genes导致报错，更改iter.times的参数设置

## 2.4 计算RISK Score----
dat_riskscore <- data 
rownames(dat_riskscore) <- rownames(data)
dat_riskscore <- na.omit(dat_riskscore)
#定义计算RISK Score的函数
risk_score_lasso <- function(prognosis_lasso,dat_riskscore){
  genes <- prognosis_lasso$genes
  lasso_coef <- prognosis_lasso[, "coef"]
  riskscore <- dat_riskscore[, genes] %>% 
    apply(1, function(x) {crossprod(as.numeric(x), lasso_coef)})
  
  risklevel <- ifelse(riskscore<=median(riskscore), "LowRisk", "HighRisk")
  risk <- data.frame(riskscore, risklevel)
  return(risk)
}

# 3 输出结果保存----
## 3.1 保存csv文件----
riskscore_lasso <- risk_score_lasso(prognosis_lasso,dat_riskscore)
riskscore_lasso <- riskscore_lasso[order(riskscore_lasso$risklevel),]
riskscore_lasso$sampleid <- rownames(riskscore_lasso)
write.csv(riskscore_lasso,"./output/7-RiskScore_LASSO.csv")
write.csv(prognosis_lasso$genes,"./output/8-LASSO_hubgenes.csv",row.names = F)

## 3.2 保存工作环境----
saveRDS(prognosis_lasso, file = "./output/9-LASSO_Coef.rds")
save.image("2-Lasso.RData.gz",compress = "gzip")

## 3.3 写入风险评分计算公式----
ps<- ''
ps2 <- ''
for(i in 1:16){
  ge <- prognosis_lasso[i, 1]
  va <- round(prognosis_lasso[i, 2], 3)
  ps <- paste0(ps, '，', ge)
  ps2 <- paste0(ps2, ' + ', ge, ' * (', va, ')')
}
print(ps)
print(ps2)
cat(ps, file = "output/报告.txt", append = TRUE)
cat(ps2, file = "output/报告.txt", append = TRUE)

