rm(list = ls())


library(openxlsx)
library(ggplot2)
library(clusterProfiler)
library(stringr)
library(tidyr)
FC = 1.5
name = " "
dir.create("3. 生信分析")
setwd("D:/Desktop/")
dir.create("1. 火山图分析")
dir.create("2. 聚类分析")
dir.create("3. GO分析")
dir.create("4. KEGG_PATHWAY分析")
setwd("../")
quan_data <- readxl::read_xlsx("1.蛋白定量统计.xlsx")
volcano_data <- quan_data[,14:15]
volcano_data$significant <- rep('no', nrow(volcano_data))
volcano_data$significant[volcano_data[2]< 0.05 & 
                           volcano_data[1] >= FC] <- 'up'
volcano_data$significant[volcano_data[2]< 0.05 &
                           volcano_data[1]<= 1/FC] <- 'down'
volcanoplot <-   ggplot(volcano_data, aes(x = log2(as.numeric(unlist(volcano_data[1]))), y = -log10(as.numeric(unlist(volcano_data[2]))))) +
  geom_point(aes(color = significant), size = 1.5) + scale_color_manual(values = c("blue", "gray", "red"))+
  xlab("log2(Fold Change)") + ylab("-log10(P-value)")+
  labs(title = paste0(name, "_volcano"), hjust = 0.5,size=30) +
  geom_hline(yintercept = 1.3, linetype = 3)+
  geom_vline(xintercept = c(-log2(1/FC), log2(1/FC)), linetype = 3)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 15, colour = 'black'),axis.text.y = element_text(size = 15, colour = 'black'),
        axis.title.x = element_text(vjust = 2,size = 14),axis.title.y = element_text(vjust = 2,size = 14),
        legend.key.size = unit(0.75,'cm'),legend.title=element_blank() ,legend.text = element_text(size = 16))
volcanoplot
ggsave(file = paste0("./3. 生信分析/1. 火山图分析/", name, "_volcano.png"), plot = volcanoplot,
       scale = 1, width = 12, height = 7.5, units = "in")
ggsave(file = paste0("./3. 生信分析/1. 火山图分析/", name, "_volcano.pdf"), plot = volcanoplot,
       scale = 1, width = 12, height = 7.5, units = "in")


library(pheatmap)
Norm <- colnames(quan_data)[grep(pattern = "^NORM_INTENSITY_", colnames(quan_data))]
data_differ <- readxl::read_xlsx("2.差异定量结果.xlsx")
heatmap_data <- data_differ[c('ACCESSION',Norm)]
heatmap_data <- data.frame(heatmap_data)
rownames(heatmap_data) <- heatmap_data[,1]
heatmap_data <- heatmap_data[,-1]
colnames(heatmap_data) <- stringr::str_sub(colnames(heatmap_data),start = 16)

annotation_col <- data.frame(Group= c(rep("Z",3),rep("293T",3)))
rownames(annotation_col) <- colnames(heatmap_data)

par(las = 2)
heatmapplot <- pheatmap::pheatmap(heatmap_data, cluster_rows = TRUE, cluster_cols = F,
                                  scale = "row", annotation_col = annotation_col,
                                  color = colorRampPalette(rev(c("red", "white","blue")))(2000),
                                  show_rownames = F, show_colnames = T,border_color = NA,
                                  #cellwidth =48,#cellheight=4, border_color = NA,
                                  fontsize_row = 2, fontsize_col = 10.5, 
                                  legend = TRUE,legend_breaks = -1.5:1.5,
                                  main = paste0(name, "_Cluster"))
heatmapplot
ggsave(file = paste0("./2. 聚类分析/",name, "_Cluster.png"), plot = heatmapplot,
       scale = 1, width = 2.1*3,height = 3.2*3, units = "in")
ggsave(file = paste0("./2. 聚类分析/",name, "_Cluster.pdf"), plot = heatmapplot,
       scale = 1, width = 2.1*3,height = 3.2*3, units = "in")

Org = "hsa"
OrgDB = "org.Hs.eg.db"
if (Org == "human") {
  Org = "hsa"
  print(Org)
  library(org.Hs.eg.db)
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
}else if (Org == "dog") {
  Org = "cfa"
  print(Org)
  OrgDB = loadDb(file = "~/org.Cf.eg.db")
}else if (Org == "eck12") {
  Org = "eco"
  print(Org)
  library(org.EcK12.eg.db)
  OrgDB = 'org.EcK12.eg.db'
}else if (Org == "bovine") {
  Org = "bta"
  print(Org)
  library(org.Bt.eg.db)
  OrgDB = 'org.Bt.eg.db'
}else if (Org == "Rhesus") {
  Org = "mcc"
  print(Org)
  library(org.Mmu.eg.db)
  OrgDB = 'org.Mmu.eg.db'
}

analysname <- name
ID <- data_differ$`Gene Symbol`
ID <- unique(ID)
ID <- bitr(ID, fromType = "SYMBOL", toType = "ENTREZID",
           OrgDb = OrgDB)
go_MF<-enrichGO(gene = ID$ENTREZID,OrgDb = OrgDB,
                ont = "MF",pvalueCutoff = 0.9999, pAdjustMethod = "BH",qvalueCutoff = 0.999,
                maxGSSize = 5000, readable = T, keyType = "ENTREZID")
go_CC<-enrichGO(gene = ID$ENTREZID,OrgDb = OrgDB,
                ont = "CC",pvalueCutoff = 0.9999, pAdjustMethod = "BH",qvalueCutoff = 0.999,
                maxGSSize = 5000, readable = T, keyType = "ENTREZID")
go_BP<-enrichGO(gene = ID$ENTREZID,OrgDb = OrgDB,
                ont = "BP",pvalueCutoff = 0.9999, pAdjustMethod = "BH",qvalueCutoff = 0.999,
                maxGSSize = 5000, readable = T, keyType = "ENTREZID")
GOGO <- createWorkbook()
addWorksheet(GOGO, sheetName = "GO_BP")
writeData(GOGO, sheet = "GO_BP", go_BP) #[go_BP@result$Count>1]
addWorksheet(GOGO, sheetName = "GO_CC")
writeData(GOGO, sheet = "GO_CC", go_CC) #[go_CC@result$Count>1]
addWorksheet(GOGO, sheetName = "GO_MF")
writeData(GOGO, sheet = "GO_MF", go_MF) #[go_MF@result$Count>1]
saveWorkbook(GOGO, paste0("./3. 生信分析/",analysname, "_GO.xlsx"),overwrite = TRUE)

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
theme_go <- theme(axis.text.x = element_text(angle = 80, hjust = 1), axis.title.x = element_text(hjust = 0.5),
                  plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), panel.background = element_blank(),
                  axis.line = element_line(colour = "black"))
p <- ggplot(data = GOOOO) + geom_bar(aes(x = Description,
                                         y = Count, fill = GO), stat = "identity") +
  scale_x_discrete(label = function(x) str_wrap(x, 64))+
  scale_fill_manual(values = c("#473C8B", "#EE4000","#FFA500")) +
  facet_grid(. ~ GO, scales = "free_x", space = "free") +
  labs(title = "The Most Enriched GO Terms", x = "GO term", y = "Numbers of Genes") +
  theme_go
p
pic_width = 9.6
pic_height = 6.3
ggsave(file = paste0("./3. 生信分析/3. GO分析/",analysname, "_GO", ".png"),
       plot = p, width = pic_width, height = pic_height,
       dpi = 600, device = "png")
ggsave(file = paste0("./3. 生信分析/3. GO分析/",analysname, "_GO", ".pdf"),
       plot = p, width = pic_width, height = pic_height,
       dpi = 600, device = "pdf")


  kegg0 <- enrichKEGG(gene = ID$ENTREZID,organism = Org,use_internal_data = F,pvalueCutoff = 0.99999,qvalueCutoff = 0.9999,maxGSSize = 5000)
kegg <- setReadable(kegg0, OrgDB, keyType = "ENTREZID")
# write.xlsx(kegg[kegg@result$Count>1], sheetName = "KEGG_Pathway", file = paste0(name,
#                                                            "_KEGG_Pathway", ".xlsx"))

KEGGdata <- createWorkbook()
addWorksheet(KEGGdata, sheetName = "KEGG_Pathway")
writeData(KEGGdata, sheet = "KEGG_Pathway", kegg)
saveWorkbook(KEGGdata, paste0("./3. 生信分析/",analysname, "_KEGG_Pathway.xlsx"),overwrite = TRUE)

if(nrow(kegg@result)>15){
  data_kegg = kegg@result[c(1:15),]
}else{
  data_kegg = kegg@result
}
data_kegg <- separate(data = data_kegg, col = GeneRatio, into = c("gene", "inputlist"), sep = "/")
data_kegg <- separate(data = data_kegg, col = BgRatio, into = c("termlist", "Bglist"), sep = "/")
data_kegg$fold1 <- (as.numeric(data_kegg$gene)/as.numeric(data_kegg$inputlist))
data_kegg$fold2 <- (as.numeric(data_kegg$termlist)/as.numeric(data_kegg$Bglist))
data_kegg$fold <- data_kegg$fold1/data_kegg$fold2
keggplot <- ggplot(data_kegg,aes(reorder(Description,fold),fold))+
  geom_point(stat ="identity",aes(size = Count,colour = pvalue))+
  scale_x_discrete(label = function(x) str_wrap(x, 64))+
  scale_colour_gradient(low = "red", high = "green")+
  scale_size_continuous(range = c(4,10.5))+
  xlab("")+ylab("Enrich Factor")+labs(title="KEGG_Pathway",hjust=0.5)+
  coord_flip()+
  theme_bw()+
  theme(
    axis.text.x=element_text(color = "black",size = 10),
    axis.text.y=element_text(color = "black",size = 12),
    axis.title.x = element_text(face="bold",color = "black"))
ggsave(file = paste0("./3. 生信分析/4. KEGG_PATHWAY分析/",analysname, "_KEGG_Pathway.png"),keggplot,
       width = pic_width, height = pic_height, dpi = 600,
       device = "png")
ggsave(file = paste0("./3. 生信分析/4. KEGG_PATHWAY分析/",analysname, "_KEGG_Pathway.pdf"),keggplot,
       width = pic_width, height = pic_height, dpi = 600,
       device = "pdf")

