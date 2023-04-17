rm(list = ls())
library(clusterProfiler) #BiocManager::install("clusterProfiler")
library(stringr) #install.packages("stringr")
library(tidyr) #install.packages("tidyr")
library(openxlsx)  #install.packages("openxlsx")
library(ggplot2)  #install.packages("ggplot2")
library(R.utils)   #install.packages("R.utils")
R.utils::setOption("clusterProfiler.download.method",'auto')



getwd()
setwd("D:/Desktop/赵玉琼老师售后分析/")
pic_width = 9.6
pic_height = 6.3

Org = "ssc"    #修改物种 主要为人(human)小鼠(mouse)大鼠(rat)三个物种的，其他的根据实际物种情况进行修???
#OrgDB = "org.Ss.eg.db"   #使用其他物种

dir.create("GO分析")
dir.create("KEGG_PATHWAY分析")

if (Org == "human") {
  Org = "hsa"
  print(Org)
  library(org.Hs.eg.db) #BiocManager::install("org.Hs.eg.db")
  OrgDB = "org.Hs.eg.db"
}else if (Org == "mouse") {
  Org = "mmu"
  print(Org)
  library(org.Mm.eg.db)
  OrgDB = "org.Mm.eg.db" 
}else if (Org == "rat") {
  Org = "rno"
  print(Org)
  library(org.Rn.eg.db)
  OrgDB = "org.Rn.eg.db"
}

name <- "并集_80个蛋白_GO.png"  #设置最后生成的GOKEGG文件名称
workbook <- "并集_80个蛋白.xlsx"  #需要做GOKEGG的数据框，默认读取第一个sheet
data <- read.xlsx(workbook,4)
View(data)
#存在多个元素时需要进行分列，否则无法识别，分列操作运行下方代???
#data <- separate_rows(data, 'Gene.Symbol',sep = "/")

ID <- data$Gene.names


ID <- unique(ID)
ID <- bitr(ID, fromType = "SYMBOL", toType = "ENTREZID",
           OrgDb = OrgDB)
go_MF<-enrichGO(gene = ID$ENTREZID,OrgDb = OrgDB,
                ont = "MF",pvalueCutoff = 0.99, pAdjustMethod = "BH",qvalueCutoff = 0.999,
                maxGSSize = 500, readable = T, keyType = "ENTREZID")
go_CC<-enrichGO(gene = ID$ENTREZID,OrgDb = OrgDB,
                ont = "CC",pvalueCutoff = 0.99, pAdjustMethod = "BH",qvalueCutoff = 0.999,
                 maxGSSize = 500, readable = T, keyType = "ENTREZID")
go_BP<-enrichGO(gene = ID$ENTREZID,OrgDb = OrgDB,
                ont = "BP",pvalueCutoff = 0.99, pAdjustMethod = "BH",qvalueCutoff = 0.999,
                maxGSSize = 500, readable = T, keyType = "ENTREZID")


GOGO <- createWorkbook()
addWorksheet(GOGO, sheetName = "GO_BP")
writeData(GOGO, sheet = "GO_BP", go_BP) 
addWorksheet(GOGO, sheetName = "GO_CC")
writeData(GOGO, sheet = "GO_CC", go_CC) 
addWorksheet(GOGO, sheetName = "GO_MF")
writeData(GOGO, sheet = "GO_MF", go_MF) 
saveWorkbook(GOGO, paste0("GO分析/",name, "_GO.xlsx"))


# 渐变条形???
#平台一部芯片作???
barplot(go_BP, title = "Gene Ontology Biological Process",
        showCategory = 15, cex.names = 0.1, font.size = 9)
ggsave(file = paste0("GO分析/",name, "_GO_BP", ".png"), width = pic_width,
       height = pic_height, dpi = 600, device = "png")
ggsave(file = paste0("GO分析/",name, "_GO_BP", ".pdf"), width = pic_width,
       height = pic_height, dpi = 600, device = "pdf")

barplot(go_CC, title = "Gene Ontology Cellular Component",
        showCategory = 15, cex.names = 0.1, font.size = 9)
ggsave(file = paste0("GO分析/",name, "_GO_CC", ".png"), width = pic_width,
       height = pic_height, dpi = 600, device = "png")
ggsave(file = paste0("GO分析/",name, "_GO_CC", ".pdf"), width = pic_width,
       height = pic_height, dpi = 600, device = "pdf")

barplot(go_MF, title = "Gene Ontology Molecular Function",
        showCategory = 15, cex.names = 0.1, font.size = 9)
ggsave(file = paste0("GO分析/",name, "_GO_MF", ".png"), width = pic_width,
       height = pic_height, dpi = 600, device = "png")
ggsave(file = paste0("GO分析/",name, "_GO_MF", ".pdf"), width = pic_width,
       height = pic_height, dpi = 600, device = "pdf")



#柱状???
if (nrow(go_MF) != 0) {
  if (nrow(go_MF) < 10) {
    go_MF <- go_MF
  }
  else {
    go_MF <- go_MF[1:10, ]
  }
}
if (nrow(go_CC) != 0) {
  if (nrow(go_CC) < 10) {
    go_CC <- go_CC
  }
  else {
    go_CC <- go_CC[1:10, ]
  }
}
if (nrow(go_BP) != 0) {
  if (nrow(go_BP) < 10) {
    go_BP <- go_BP
  }
  else {
    go_BP <- go_BP[1:10, ]
  }
}
go_BP$GO <- rep("BP", nrow(go_BP))
go_CC$GO <- rep("CC", nrow(go_CC))
go_MF$GO <- rep("MF", nrow(go_MF))
GOOOO <- rbind(go_BP, go_CC, go_MF)
GOOOO <- GOOOO[, c(2, 9, 10)]
p <- ggplot(data = GOOOO) + geom_bar(aes(x = Description,
                                         y = Count, fill = GO), stat = "identity") +
  scale_x_discrete(label = function(x) str_wrap(x, 64))+
  scale_fill_manual(values = c("#473C8B", "#EE4000", "#FFA500")) + 
  facet_grid(. ~ GO, scales = "free_x", space = "free") + 
  labs(title = "The Most Enriched GO Terms",x = "GO term", y = "Numbers of Genes") +
  theme(axis.text.x = element_text(angle = 80,hjust = 1), 
        axis.title.x = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
p
ggsave(file = paste0("GO分析/",name, "_GO", ".png"),
       plot = p, width = pic_width, height = pic_height,
       dpi = 600, device = "png")
ggsave(file = paste0("GO分析/",name, "_GO", ".pdf"),
       plot = p, width = pic_width, height = pic_height,
       dpi = 600, device = "pdf")

  
#kegg

  kegg0 <- enrichKEGG(gene = ID$ENTREZID, organism = Org, 
                    keyType = "kegg", pvalueCutoff = 0.99,use_internal_data =T, qvalueCutoff = 0.99,
                    maxGSSize = 5000)
kegg <- setReadable(kegg0, OrgDB, keyType = "ENTREZID")


KEGGdata <- createWorkbook()
addWorksheet(KEGGdata, sheetName = "KEGG_Pathway")
writeData(KEGGdata, sheet = "KEGG_Pathway", kegg)
saveWorkbook(KEGGdata, paste0("KEGG_PATHWAY分析/",name, "_KEGG_Pathway.xlsx"),overwrite = T)

dotplot(kegg, showCategory = 15, title = "KEGG_Pathway", font.size = 9)+
 scale_colour_gradient(low = "red",high = "green")

ggsave(file = paste0("KEGG_PATHWAY分析/",name, "_KEGG_Pathway", ".png"), 
       width = pic_width, height = pic_height, dpi = 600, 
       device = "png")
ggsave(file = paste0("KEGG_PATHWAY分析/",name, "_KEGG_Pathway", ".pdf"), 
       width = pic_width, height = pic_height, dpi = 600, 
       device = "pdf")


