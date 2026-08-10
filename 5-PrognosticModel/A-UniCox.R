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
required_packages <- c("survival", "forestplot", "magrittr")

#调用函数进行检查和安装
install_if_missing(required_packages)

## 1.3 加载R包----

for(i in c("survival", "forestplot", "magrittr")){
  library(i, character.only = T)
}

cox_extr <- function(fit){
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

# 2 cox分析----
#' 进行cox分析
#' 
#' @param time(numeric) 生存时间
#' @param event(factor) 生存结局, Alive=0, Dead =1，需要转换为0, 1
#' @param covariates(character) 自变量
#'   \cr 如果有age，age需要转变为numeric
#'   \cr 如果有stage，stage需要将IA, IB, IIA, IIB, ..之类的修改为I, II, III..
#' @param data(data.frame) 输入数据，包含cox所需的变量
#' @param type(character) 单因素cox选择"unicox"，多因素cox选择"mulcox"
#' @param save(logic) TRUE or FALSE, 是否需要保存
#' @param path(character) 保存路径, 默认当前文件夹
#' @param name(character) 文件命名
#' 
#' @return (data.frame) cox分析结果
#' 
#' @export

## 2.1 cox分析函数的定义----
cox_analysis <- function(time, event, covariates, data, type=c("unicox","mulcox"),
                         pvalue=0.05, save=TRUE, path="./", name="unicox.csv"){
  data <- na.omit(data)
  type <- match.arg(type)
  if(!all(data[, event] %in% c(0,1))) stop("event must be either 0 or 1")
  if(type == "unicox") {
    #对covariates中的每个因素都构建一个公式，默认使用OS.time和OS
    formulas <- sapply(covariates,
                       function(x) as.formula(paste0('Surv(', time, ",", event, ')~', 
                                                     paste0("`", x, "`", sep = ""))))
    
    #进行单因素cox回归
    dat_res <- lapply(formulas, function(x){coxph(x, data = data)}) %>% 
      lapply(cox_extr) %>% #提取单因素cox回归结果
      c(use.names = F) %>% 
      do.call(what = rbind) %>% #讲结果合并为数据框 
      dplyr::filter(p.value <= pvalue)
  } else {
    dat_res <- as.formula(paste('Surv(', time, ",", event, ')~',  
                                paste(covariates, collapse = "+"))) %>% #构建公式
      coxph(data=data) %>% #进行多因素cox回归
      cox_extr #提取多因素cox回归结果
  }
  #保存结果
  if(save) {
    if(is.null(path)|is.null(name)){
      stop("No path or name")
    }else{
      write.csv(dat_res, paste0(path, "/", name))
    }
  }
  return(dat_res)
}

#' cox分析可视化 
#' @param data_res(data.frame) cox_analysis函数得到的数据框
#' @param type(character) 单因素cox选择"unicox"，多因素cox选择"mulcox"
#' @param boxcol(character) HR平均值方块的颜色
#' @param height, width(numeric) 图片的高度, 图片的宽度;尺寸按照pdf()标准
#' @param save(logic) TRUE or FALSE, 是否需要保存
#' @param path(character) 保存路径，默认当前文件夹
#' @param name(character) 图片的名字
#' 
#' @return 森林图
#' 
#' @export
cox_plot <- function(dat_res, type=c("unicox", "mulcox"), boxcol="red", 
                     height=4, width=5, save=TRUE, path="./", name="unicox.pdf"){
  type <- match.arg(type)
  dat_res <- dat_res %>% 
    dplyr::mutate(GeneID = rownames(.), 
                  HR = paste(dat_res$mean, " [", dat_res$lower, " - ", 
                             dat_res$upper, "]", sep = "")) %>%
    dplyr::select("GeneID", "p.value", "HR", "mean", "lower", "upper")
  colnames(dat_res)[2] <- "P Value"
  labeltext <- rbind(colnames(dat_res)[1:3], dat_res[, 1:3])
  HR <- rbind(rep(NA, 3), dat_res[, 4:6])
  #画图
  p <- forestplot(labeltext = labeltext,             
                  HR,
                  zero = 1, 
                  lwd.ci = 1,#HR线的宽度
                    title=ifelse(type == "unicox", "Univariable Cox Regression Analysis",
                               "multivariable Cox regression analysis"),
                  colgap = unit(2, 'mm'),
                  xlab = "HR", 
                  
                  #坐标轴和变量的字体大小
                  txt_gp=fpTxtGp(ticks = gpar(cex = 0.8), xlab = gpar(cex = 0.8),
                                 title = gpar(cex = 1.2), cex = 1) ,
                  boxsize = 0.1, 
                  col=fpColors(box = boxcol, line = "black", zero = "black")#修改颜色
  ) 
  if(save) {
    if(is.null(path)|is.null(name)){
      stop("No path or name")
    }else{
      pdf(paste0(path, "/", name), height = height, width = width)
      print(p)
      dev.off()
    }
  }
  print(p)
}  

#sample
## 2.2 读取数据----
dat_tumor <- data.table::fread("input/Matrix_FPKM.csv", data.table = F)
rownames(dat_tumor) <- dat_tumor[,1]
dat_tumor <- dat_tumor[,-1]
pd <- data.table::fread("input/PD.csv", data.table = F)
pd <- pd[pd$Group == "STAD",]
genes <- data.table::fread("input/3-RDEGs.csv", data.table = F)
dat_tumor <- dat_tumor[,pd$ID]
dat_tumor <- t(dat_tumor[genes$x,])

identical(pd$ID,rownames(dat_tumor))
data <- cbind(pd[,4:5],dat_tumor) #不能无脑跑，注意pd和dat_tumor是否对齐！！！
rownames(data) <- rownames(dat_tumor)
covariates <- c(colnames(data)[3:ncol(data)])
dir.create("output")

## 2.3 进行cox分析----
dat_res <- cox_analysis(time = "OS.time", event = "OS",
                        covariates = covariates, data = data,
                        type = "unicox", pvalue = 0.05,
                        save = T, path = "./output", name = "1-UnicoxGenes.csv")

cox_plot(dat_res, type = "unicox", boxcol = "#EE7785", height=10, width = 8.1,
         save = T, path = "./output", name = "2-Unicox_ForestPlot.pdf")
dev.off()
# 3 输出结果保存----
## 3.1 保存工作环境----
saveRDS(dat_res, file = "./output/3-cox_res.rds")
saveRDS(data, file = "./output/4-cox_data.rds")
save.image("1-Unicox.RData.gz",compress = "gzip")
