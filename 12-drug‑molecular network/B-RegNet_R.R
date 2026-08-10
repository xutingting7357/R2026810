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
required_packages <- c("magrittr","ggraph","tidygraph")

#调用函数进行检查和安装
install_if_missing(required_packages)

## 1.3 加载R包----
library(magrittr)
library(ggraph)
library(tidygraph)

# 2 mRNA-Drug网络图----
## 2.1 mRNA-Drug网络图数据处理----
load("output/18-mRNA_Drug.RData")
edges <- mRNA_Drug %>% 
  data.frame(check.names = F) %>% 
  `names<-`(c("from", "to")) %>% 
  dplyr::distinct()
nodes <- Drug_attr %>% 
  data.frame(check.names = F) %>% 
  `names<-`(c("id", "type"))
nodes$type <- factor(nodes$type, levels = c("mRNA", "Drug"))


net.tidy <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
set.seed(123)

## 2.2 绘制mRNA-Drug网络图----
ggraph(net.tidy, layout = "graphopt") + 
  geom_node_point(aes(color = type, fill = type, shape = type, size = type)) + 
  # 点信息
  geom_edge_link(alpha = 1, width = 0.1) +  # 边信息
  geom_node_text(aes(label = id), repel = TRUE, check_overlap = T, size = 5) + 
  # 增加节点的标签，reple避免节点重叠
  scale_size_manual(values = c(6, 3)) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c("#008B45", "#3B4992")) + 
  scale_color_manual(values = c("#008B45", "#3B4992")) + 
  theme_graph() +
  coord_cartesian(clip = "off")
ggsave("output/19-Drug.pdf", width = 8, height = 8, device = cairo_pdf)
