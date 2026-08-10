library(dplyr)

GSEA_output <- read.csv("output/GSEA_MSigDb_result.csv",header = T)
#GSEA_output <- openxlsx::read.xlsx("GSEA_output.xlsx", sheet = 1)
colnames(GSEA_output)
data1 <- GSEA_output %>% 
  dplyr::filter(p.adjust<0.05 & qvalue < 0.25)
data2 <- data1 %>% 
  dplyr::filter(stringr::str_detect(Description, 'MITOPHAGY|JAK|TP53|JNK|PI3K|
                                    NF|WNT|HEDGEHOG|NOTCH|TGF|IL|INTERLEUKIN|INFLAMM|
                                    NECROPTOSIS|AUTOPH|PYROPTOSIS|FERROPTOSIS|IRON|M6A|M5C|
                                    M1A|HYPOXIA|ENDOPLAS|OXIDATIVE|GLYCOLYSIS|FATTY|HISTONE|
                                    APOPTOSIS|AGEING|SENESCENCE|EPITHELIAL|EMT|EXOSOME|METABO|CIRCADIAN'))
write.csv(data2,"output/GSEA明星通路.csv")
write.csv(data1,"output/GSEA_padj_0.05_q_0.25.csv")
