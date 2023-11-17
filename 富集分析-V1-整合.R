library(readxl)
library(clusterProfiler, quietly = T)
library(openxlsx)
library(stringr)
library(ggplot2)
library(stringr)
pic_width = 9.6
pic_height = 6.4
Org = "hsa"#人 hsa 小鼠mmu 大鼠 rno
print(Org)
#library(org.Hs.eg.db)
#library(org.Hs.eg.db)
OrgDB <- 'org.Hs.eg.db'

 library(org.Hs.eg.db)
# library(org.Ss.eg.db)
# library(org.Mm.eg.db)
getwd()
setwd("D:/Desktop/2023.11.13-HWEXM03220231010-中国医学科学院阜外医院-鞠文浩老师-Labelfree定量蛋白质谱项目总结报告-rhj/2. 数据结果/2. 差异统计/")
workbook <- "2.差异定量结果.xlsx"
n <- length(readxl::excel_sheets(workbook))


for (i in 1:1) {
  #setwd("C:/Users/wayen/Desktop")
  print(i)
  name <- readxl::excel_sheets(workbook)[i]
  print(name)
  data <- read_xlsx(workbook, i)
    ID = unlist(strsplit(as.character(data$Genes), "/"),use.names = F)
  #ID <- c("MUT","PDK3","GSTO1","HRG","PDLIM4","VAT1L","APOA4","F9","FMOD")
  ID <- unique(ID)
  ID <- bitr(ID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDB)#SYMBOL,UNIPROT
  go_MF<-enrichGO(gene = ID$ENTREZID,OrgDb = OrgDB,
                  ont = "MF",pvalueCutoff = 0.99, pAdjustMethod = "BH",qvalueCutoff = 0.99,
                  maxGSSize = 5000, readable = T, keyType = "ENTREZID")
  go_CC<-enrichGO(gene = ID$ENTREZID,OrgDb = OrgDB,
                  ont = "CC",pvalueCutoff = 0.99, pAdjustMethod = "BH",qvalueCutoff = 0.99,
                   maxGSSize = 5000, readable = T, keyType = "ENTREZID")
  go_BP<-enrichGO(gene = ID$ENTREZID,OrgDb = OrgDB,
                  ont = "BP",pvalueCutoff = 0.99, pAdjustMethod = "BH",qvalueCutoff = 0.99,
                  maxGSSize = 5000, readable = T, keyType = "ENTREZID")

  #go_simplity_MF <- simplify(go_MF)
  #go_simplity_CC <- simplify(go_CC)
  #go_simplity_BP <- simplify(go_BP)

  folder=paste0("./",name,"_富集分析")
  dir.create(folder)
  setwd(folder)
  GOGO <- createWorkbook()
  addWorksheet(GOGO, sheetName = "GO_BP")
  writeData(GOGO, sheet = "GO_BP", go_BP[go_BP@result$Count>1])
  addWorksheet(GOGO, sheetName = "GO_CC")
  writeData(GOGO, sheet = "GO_CC", go_CC[go_CC@result$Count>1])
  addWorksheet(GOGO, sheetName = "GO_MF")
  writeData(GOGO, sheet = "GO_MF", go_MF[go_MF@result$Count>1])
  saveWorkbook(GOGO, paste0(name, "_GO.xlsx"))


# 渐变条形???
#
#   barplot(go_BP, title = "Gene Ontology Biological Process",
#           showCategory = 15, cex.names = 0.1, font.size = 9)
#   ggsave(file = paste0(name, "_GO_BP", ".png"), width = pic_width,
#          height = pic_height, dpi = 600, device = "png")
#   ggsave(file = paste0(name, "_GO_BP", ".pdf"), width = pic_width,
#          height = pic_height, dpi = 600, device = "pdf")
#
#   barplot(go_CC, title = "Gene Ontology Cellular Component",
#           showCategory = 15, cex.names = 0.1, font.size = 9)
#   ggsave(file = paste0(name, "_GO_CC", ".png"), width = pic_width,
#          height = pic_height, dpi = 600, device = "png")
#   ggsave(file = paste0(name, "_GO_CC", ".pdf"), width = pic_width,
#          height = pic_height, dpi = 600, device = "pdf")
#
#   barplot(go_MF, title = "Gene Ontology Molecular Function",
#           showCategory = 15, cex.names = 0.1, font.size = 9)
#   ggsave(file = paste0(name, "_GO_MF", ".png"), width = pic_width,
#          height = pic_height, dpi = 600, device = "png")
#   ggsave(file = paste0(name, "_GO_MF", ".pdf"), width = pic_width,
#          height = pic_height, dpi = 600, device = "pdf")

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
  theme_go <- theme(axis.text.x = element_text(angle = 80,
                                               hjust = 1), axis.title.x = element_text(hjust = 0.5),
                    plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), panel.background = element_blank(),
                    axis.line = element_line(colour = "black"))
  p <- ggplot(data = GOOOO) + geom_bar(aes(x = Description,
                                           y = Count, fill = GO), stat = "identity") +
    scale_x_discrete(label = function(x) str_wrap(x, 64))+
    scale_fill_manual(values = c("#473C8B", "#EE4000",
                                 "#FFA500")) + facet_grid(. ~ GO, scales = "free_x",
                                                          space = "free") + labs(title = "The Most Enriched GO Terms",
                                                                                 x = "GO term", y = "Numbers of Genes") +
    theme_go
  ggsave(file = paste0(name, "_GO", ".png"),
         plot = p, width = pic_width, height = pic_height,
         dpi = 600, device = "png")
  ggsave(file = paste0(name, "_GO", ".pdf"),
         plot = p, width = pic_width, height = pic_height,
         dpi = 600, device = "pdf")

  
#kegg
  
  kegg0 <- enrichKEGG(gene = ID$ENTREZID, organism = Org,
                      keyType = "kegg", pvalueCutoff = 1,use_internal_data =F,qvalueCutoff = 1,
                      maxGSSize = 500,minGSSize = 2)
  kegg <- setReadable(kegg0, OrgDB, keyType = "ENTREZID")
  # write.xlsx(kegg[kegg@result$Count>1], sheetName = "KEGG_Pathway", file = paste0(name, 
  #                                                            "_KEGG_Pathway", ".xlsx"))
  
  KEGGdata <- createWorkbook()
  addWorksheet(KEGGdata, sheetName = "KEGG_Pathway")
  writeData(KEGGdata, sheet = "KEGG_Pathway", kegg)
  saveWorkbook(KEGGdata, paste0(name, "_KEGG_Pathway.xlsx"))
#kegg@result$Description=gsub("- Mus musculus (house mouse)","",x = kegg@result$Description[1])
  dotplot(kegg, showCategory = 15, title = "KEGG_Pathway", font.size = 9,color ="pvalue")+
   scale_colour_gradient(low = "red", high = "green")
  
  ggsave(file = paste0(name, "_KEGG_Pathway", ".png"), 
         width = pic_width, height = pic_height, dpi = 600, 
         device = "png")
  ggsave(file = paste0(name, "_KEGG_Pathway", ".pdf"), 
         width = pic_width, height = pic_height, dpi = 600, 
         device = "pdf")
  setwd("D:/Desktop/")
}

