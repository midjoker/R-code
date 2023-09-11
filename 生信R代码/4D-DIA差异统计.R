# install.packages("devtools")
# library(devtools)  # You need to install this package!
# install_github("ltorgo/DMwR2",ref="master")

library(DMwR2)
library(tibble)
library(openxlsx)
library(ggplot2)
library(clusterProfiler)
library(stringr)
library(dplyr)

######设置调变倍数和物种human mouse rat
FC = 1.2
Org = "rat"

############################################################
#######################缺失值填充###########################
############################################################
data = openxlsx::read.xlsx("Proteins.xlsx",1,rowNames = T)
dat_metaname = colnames(data)[1:4]
dat_meta = data[,(1:4)]
dat = data[,5:ncol(data)]
dat = as.data.frame(apply(dat, 2, as.numeric))
row.names(dat) = row.names(data)

dat = dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,]
dim(dat)
knnOutput <- knnImputation(dat,scale = T,k = 10)  # 使用KNN插值,k为选取的最近邻样本的个数
outputdata <- merge(rownames_to_column(dat_meta,"PG.ProteinAccessions"),
                    rownames_to_column(knnOutput,"PG.ProteinAccessions"),
                                       by = "PG.ProteinAccessions")

write.xlsx(outputdata,"./1. 原始数据/Proteins_imputed.xlsx")


############################################################
#######################差异统计#############################
############################################################

protein_quandata <- knnOutput
logdata <- log2(protein_quandata)
View(logdata)
#输入样本及样本分组信息
sampledata <- openxlsx::read.xlsx("sampegroup.xlsx",3)
View(sampledata)
#输入需要进行比较的组的信息
Compare <- openxlsx::read.xlsx("sampegroup.xlsx",4)
View(Compare)



dir.create("2. 定量统计")
dir.create("3. 生信分析")
setwd("./3. 生信分析")
dir.create("1. 火山图分析")
dir.create("2. 聚类分析")
dir.create("3. GO分析")
dir.create("4. KEGG_PATHWAY分析")
setwd("../")
diffs <- matrix(nrow=nrow(Compare), ncol=4)
ratio_pvalue <- rep(0, nrow(outputdata))
DIFF <- createWorkbook()
quandata<- createWorkbook()


for (i in 1:nrow(Compare)) {
  #筛选出进行组间比较的组所对应的样本
  groupA <- subset(sampledata, select = "SampleID", sampledata[["Group"]] ==  Compare[i, 1])
  groupB <- subset(sampledata, select = "SampleID", sampledata[["Group"]] ==  Compare[i, 2])
  
  type <- factor(c(rep(Compare[i, 1], nrow(groupA)), rep(Compare[i,2], nrow(groupB))), levels = c(Compare[i, 1], Compare[i,2]))
  Pvalue <- rep(0, nrow(protein_quandata))
  Ratio <- rep(0, nrow(protein_quandata))
  meanA <- rep(0, nrow(protein_quandata))
  meanB <- rep(0, nrow(protein_quandata))
  
  
  for (j in 1:nrow(logdata)) {
    if (is.na(mean(as.numeric(logdata[c(groupA$SampleID,groupB$SampleID)][j, ])))) {
      Pvalue[j] = "NA"
    }else{
      WT <- t.test(as.numeric(unlist(logdata[c(groupA$SampleID,groupB$SampleID)][j,])) ~ type, exact = FALSE, 
                   alternative = "two.sided", paired = FALSE, 
                   var.equal = FALSE)
      Pvalue[j] <- WT$p.value
    }
  }
  for (t in 1:nrow(protein_quandata)) {
    meanA[t] <- mean(as.numeric(protein_quandata[t, groupA$SampleID]))
    meanB[t] <- mean(as.numeric(protein_quandata[t, groupB$SampleID]))
    Ratio[t] = meanA[t]/meanB[t]
  }
  #Ratio <- log2(Ratio)
  #Pvalue_adjust <- p.adjust(Pvalue,method = "BH")

  
  quan_data <- cbind(Ratio, Pvalue, protein_quandata[c(groupA$SampleID,groupB$SampleID)])
  colnames(quan_data) <- c(paste0(Compare[i, 1], "_VS_", Compare[i, 2],"_Ratio"), 
                           paste0(Compare[i, 1], "_VS_", Compare[i, 2],"_Pvalue"),
                           colnames(protein_quandata[c(groupA$SampleID,groupB$SampleID)]))
  ratio_pvalue <- cbind(ratio_pvalue,quan_data[c(paste0(Compare[i, 1], "_VS_", Compare[i, 2],"_Ratio"),
                        paste0(Compare[i, 1], "_VS_", Compare[i, 2],"_Pvalue"))])
  #write.xlsx(quan_data,paste0(Compare[i, 1], "_VS_", Compare[i, 2],"_result.xlsx"),row.names = T)
  #dir.create("./1. 火山图分析")
  #setwd("./1. 火山图分析")
  volcano_data <- quan_data[,1:2]
  volcano_data$significant <- rep('no', nrow(volcano_data))
  volcano_data$significant[volcano_data[2]< 0.05 &
                             volcano_data[1]>= FC] <- 'up'
  volcano_data$significant[volcano_data[2]< 0.05 &
                             volcano_data[1]<= 1/FC] <- 'down'
  volcanoplot <-   ggplot(volcano_data, aes(x = log2(as.numeric(unlist(volcano_data[1]))), y = -log10(as.numeric(unlist(volcano_data[2]))))) +
    geom_point(aes(color = significant), size = 1.5) + scale_color_manual(values = c("blue", "gray", "red"))+
    xlab("log2(Fold Change)") + ylab("-log10(P-value)")+
    labs(title = paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_volcano"), hjust = 0.5,size=30) +
    geom_hline(yintercept = 1.3, linetype = 3)+
    geom_vline(xintercept = c(-log2(1/FC), log2(1/FC)), linetype = 3)+
    theme_bw()+
    theme(axis.text.x = element_text(size = 15, colour = 'black'),axis.text.y = element_text(size = 15, colour = 'black'),
          axis.title.x = element_text(vjust = 2,size = 14),axis.title.y = element_text(vjust = 2,size = 14),
          legend.key.size = unit(0.75,'cm'),legend.title=element_blank() ,legend.text = element_text(size = 16))
  volcanoplot
  ggsave(file = paste0("./3. 生信分析/1. 火山图分析/",Compare[i, 1], "_VS_", Compare[i,2], "_volcano.png"), plot = volcanoplot,
         scale = 1, width = 12, height = 7.5, units = "in")
  ggsave(file = paste0("./3. 生信分析/1. 火山图分析/",Compare[i, 1], "_VS_", Compare[i, 2], "_volcano.pdf"), plot = volcanoplot,
         scale = 1, width = 12, height = 7.5, units = "in")
  
  pvalue <- paste0(Compare[i, 1], "_VS_",Compare[i, 2], "_Pvalue")
  ratio <- paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_Ratio")
  data_differ <- quan_data[which(quan_data[pvalue] < 0.05),]
  data_differ <- quan_data[which(((quan_data[ratio] > FC)|(quan_data[ratio] < 1/FC))&quan_data[pvalue] < 0.05),]
  data_differ1 <- data_differ[order(data_differ[,ratio], decreasing = TRUE),]
  
  data_differ <- merge(rownames_to_column(dat_meta,"PG.ProteinAccessions"),
                       rownames_to_column(data_differ1,"PG.ProteinAccessions"),by = "PG.ProteinAccessions",
                       all.x = FALSE)
  
  headerStyle1 <- createStyle(fontName = "Arial", fontSize = 10,fgFill = "#575757",
                              fontColour = "white", halign = "center", 
                              valign = "center", wrapText = TRUE)
  headerStyle2 <- createStyle(fontName = "Arial", fontSize = 10,fgFill = "#4682B4",
                              fontColour = "white", halign = "center", 
                              valign = "center", wrapText = TRUE)
  datastyle <- createStyle(fontName = "Arial", fontSize = 10, wrapText = FALSE)
  name <- paste0(Compare[i, 1], "_VS_", Compare[i, 2])
  
  addWorksheet(DIFF, sheetName = paste0(Compare[i, 1], "_VS_", Compare[i, 2]))
  setColWidths(DIFF, name,cols = c(seq(from = 1, to = ncol(data_differ), by = 1)), widths = c(rep(14, ncol(data_differ))))
  setRowHeights(DIFF, name, rows = 1, heights = 50)
  writeData(DIFF, name, data_differ)
  freezePane(DIFF, sheet = name, firstActiveRow = 2)
  addStyle(DIFF, sheet = name, headerStyle1, rows = 1, cols = 1:5, gridExpand = TRUE)
  addStyle(DIFF, sheet = name, headerStyle2, rows = 1, cols = 6:ncol(data_differ), gridExpand = TRUE)
  addStyle(DIFF, sheet = name, datastyle, rows = 2:nrow(data_differ), cols = 1:ncol(data_differ), gridExpand = TRUE)
  diffs[i,1] <- paste0(Compare[i, 1], "_VS_", Compare[i, 2])
  diffs[i,2] <- nrow(data_differ[data_differ[6]>1,])
  diffs[i,3] <- nrow(data_differ[data_differ[6]<1,])
  diffs[i,4] <- nrow(data_differ)
  
  
  library(pheatmap)
  heatmap_data <- data_differ[c(groupA$SampleID,groupB$SampleID)]
  # rownames(heatmap_data) <- heatmap_data[,1]
  # heatmap_data <- heatmap_data[,-1]
  heatmap_data <- log2(heatmap_data)
  #colnames(heatmap_data) <- stringr::str_sub(colnames(heatmap_data),start = 26)
  
  annotation_colA <- subset(sampledata, sampledata$Group == Compare[i,1])
  annotation_colB <- subset(sampledata, sampledata$Group == Compare[i,2])
  annotation_colX <- rbind(annotation_colA,annotation_colB)
  #annotation_colX$SamplesName <- stringr::str_sub(annotation_colX$SamplesName, start = 26)
  annotation_col <- data.frame(Group= annotation_colX$Group)
  rownames(annotation_col) <- annotation_colX[,1]
  
  par(las = 2)
  heatmapplot <- pheatmap(heatmap_data, cluster_rows = TRUE, cluster_cols = T,
                                    scale = "row", annotation_col = annotation_col,
                                    color = colorRampPalette(rev(c("red", "white","blue")))(2000),
                                    show_rownames = F, show_colnames = T,border_color = NA,
                                    #cellwidth =48,#cellheight=4, border_color = NA,
                                    fontsize_row = 5.2, fontsize_col = 10.5, 
                                    legend = TRUE,legend_breaks = -1.5:1.5,
                                    main = paste0(Compare[i,1], "_VS_", Compare[i, 2], "_Cluster"))
  #setwd("../")
  # dir.create("./2. 聚类分析")
  ggsave(file = paste0("./3. 生信分析/2. 聚类分析/",Compare[i, 1], "_VS_", Compare[i,2], "_Cluster.png"), plot = heatmapplot,
         scale = 1, width = 2.1*3,height = 3.2*3, units = "in")
  ggsave(file = paste0("./3. 生信分析/2. 聚类分析/",Compare[i,1], "_VS_", Compare[i, 2], "_Cluster.pdf"), plot = heatmapplot,
         scale = 1, width = 2.1*3,height = 3.2*3, units = "in")
  
  # dir.create("./3. GO分析")
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
  }
  analysname <- paste0(Compare[i, 1], "_VS_", Compare[i, 2])
  library(tidyr)
  data_differ12 <- separate_rows(data_differ, 'PG.Genes',sep = ";")
  ID <- unique(data_differ12$PG.Genes)
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
  saveWorkbook(GOGO, paste0("./3. 生信分析/3. GO分析/",analysname, "_GO.xlsx"),overwrite = TRUE)
  
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
  
  #dir.create("4. KEGG_PATHWAY分析")
  library(R.utils)
  R.utils::setOption("clusterProfiler.download.method",'auto')
  kegg0 <- enrichKEGG(gene = ID$ENTREZID, organism = Org,
                      keyType = "kegg", pvalueCutoff = 0.9999, use_internal_data = FALSE,qvalueCutoff = 0.99999999,
                      maxGSSize = 5000)
  kegg <- setReadable(kegg0, OrgDB, keyType = "ENTREZID")
  # write.xlsx(kegg[kegg@result$Count>1], sheetName = "KEGG_Pathway", file = paste0(name,
  #                                                            "_KEGG_Pathway", ".xlsx"))
  
  KEGGdata <- createWorkbook()
  addWorksheet(KEGGdata, sheetName = "KEGG_Pathway")
  writeData(KEGGdata, sheet = "KEGG_Pathway", kegg)
  saveWorkbook(KEGGdata, paste0("./3. 生信分析/4. KEGG_PATHWAY分析/",analysname, "_KEGG_Pathway.xlsx"),overwrite = TRUE)
  
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
}
headerstyle3 <- createStyle(fontName = "Arial",fontSize = 11, fgFill = "#808080",
                            fontColour = "white", halign = "center",
                            valign = "center", wrapText = TRUE)
datastyle <- createStyle(fontName = "Arial", fontSize = 10, valign = "center",
                         wrapText = FALSE)
saveWorkbook(DIFF, paste0("./2. 定量统计/", "2. 差异定量结果.xlsx"),overwrite = TRUE)

ratio_pvalue <- ratio_pvalue[,-1]
quandata_C <- merge(rownames_to_column(dat_meta,"PG.ProteinAccessions"),
                  rownames_to_column(ratio_pvalue,"PG.ProteinAccessions"),
                  by = "PG.ProteinAccessions")

quan_dataD <- merge(quandata_C,
                    rownames_to_column(knnOutput,"PG.ProteinAccessions"),
                    by = "PG.ProteinAccessions")



headerStyle1 <- createStyle(fontName = "Arial", fontSize = 10,fgFill = "#575757",
                            fontColour = "white", halign = "center",
                            valign = "center", wrapText = TRUE)
headerStyle2 <- createStyle(fontName = "Arial", fontSize = 10,fgFill = "#4682B4",
                            fontColour = "white", halign = "center",
                            valign = "center", wrapText = TRUE)
datastyle <- createStyle(fontName = "Arial", fontSize = 10, wrapText = FALSE)
addWorksheet(quandata, sheetName = "蛋白定量统计")

setColWidths(quandata, "蛋白定量统计",
             cols = c(seq(from = 1, to = ncol(quan_dataD),
                          by = 1)), widths = c(rep(14, ncol(quan_dataD))))
setRowHeights(quandata, "蛋白定量统计", rows = 1, heights = 50)
writeData(quandata, "蛋白定量统计", quan_dataD)
freezePane(quandata, sheet = "蛋白定量统计", firstActiveRow = 2)
addStyle(quandata, sheet = "蛋白定量统计", headerStyle1, rows = 1, cols = 1:5, gridExpand = TRUE)
addStyle(quandata, sheet = "蛋白定量统计", headerStyle2, rows = 1, cols = 6:ncol(quan_dataD), gridExpand = TRUE)
addStyle(quandata, sheet = "蛋白定量统计", datastyle, rows = 2:nrow(quan_dataD), cols = 1:ncol(quan_dataD), gridExpand = TRUE)
saveWorkbook(quandata, paste0("./2. 定量统计/", "1. 蛋白定量统计.xlsx"),overwrite = TRUE)

colnames(diffs) <- c("组间比较","上调差异","下调差异","汇总")
diffs <- as.data.frame(diffs)
write.xlsx(diffs,"差异结果统计.xlsx")
