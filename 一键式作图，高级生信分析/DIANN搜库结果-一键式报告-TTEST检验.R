library(openxlsx)
library(stringr)
dir.create("2. 数据结果")
setwd("./2. 数据结果")
dir.create("1. 原始数据")
setwd("../")
###########proteins#################
inputfiles <- read.xlsx("parameter.xlsx",2)

Proteins <- read.delim("report.pg_matrix.tsv",header = T)
if (length(colnames(Proteins))-5==length(inputfiles$Sample)) {
  colnames(Proteins)[6:length(colnames(Proteins))] <- inputfiles$Sample
}
orgprotein <- createWorkbook()
headerStyle1 <- createStyle(fontName = "Arial", fontSize = 10,fgFill = "#575757",
                            fontColour = "white", halign = "center", 
                            valign = "center", wrapText = TRUE)
headerStyle2 <- createStyle(fontName = "Arial", fontSize = 10,fgFill = "#4682B4",
                            fontColour = "white", halign = "center", 
                            valign = "center", wrapText = TRUE)
datastyle <- createStyle(fontName = "Arial", fontSize = 10, wrapText = FALSE)

addWorksheet(orgprotein, sheetName = "Proteins")
setColWidths(orgprotein, "Proteins", 
             cols = c(seq(from = 1, to = ncol(Proteins), 
                          by = 1)), widths = c(rep(14, ncol(Proteins))))
setRowHeights(orgprotein, "Proteins", rows = 1, heights = 50)
writeData(orgprotein, "Proteins", Proteins)
freezePane(orgprotein, sheet = "Proteins", firstActiveRow = 2)
addStyle(orgprotein, sheet = "Proteins", headerStyle1, rows = 1, cols = 1:5, gridExpand = TRUE)
addStyle(orgprotein, sheet = "Proteins", headerStyle2, rows = 1, cols = 6:ncol(Proteins), gridExpand = TRUE)
addStyle(orgprotein, sheet = "Proteins", datastyle, rows = 2:(nrow(Proteins)+1), cols = 1:ncol(Proteins), gridExpand = TRUE)
saveWorkbook(orgprotein, "./2. 数据结果/1. 原始数据/Proteins.xlsx", overwrite = TRUE)


###########peptides#################
peptides <- read.delim("report.pr_matrix.tsv",sep = "\t",header = T, stringsAsFactors = FALSE,quote = "")
colnames(peptides)[(length(colnames(peptides))-length(inputfiles$Sample)+1):length(colnames(peptides))] <- inputfiles$Sample
orgpeptides <- createWorkbook()
addWorksheet(orgpeptides, sheetName = "Peptides")
setColWidths(orgpeptides, "Peptides", 
             cols = c(seq(from = 1, to = ncol(peptides), 
                          by = 1)), widths = c(rep(14, ncol(peptides))))
setRowHeights(orgpeptides, "Peptides", rows = 1, heights = 50)
writeData(orgpeptides, "Peptides", peptides)
freezePane(orgpeptides, sheet = "Peptides", firstActiveRow = 2)
addStyle(orgpeptides, sheet = "Peptides", headerStyle1, rows = 1, cols = 1:10, gridExpand = TRUE)
addStyle(orgpeptides, sheet = "Peptides", headerStyle2, rows = 1, cols = 11:ncol(peptides), gridExpand = TRUE)
addStyle(orgpeptides, sheet = "Peptides", datastyle, rows = 2:(nrow(peptides)+1), cols = 1:ncol(peptides), gridExpand = TRUE)
saveWorkbook(orgpeptides, "./2. 数据结果/1. 原始数据/Peptides.xlsx", overwrite = TRUE)


###########差异统计情况#################
proteindata <- Proteins

#######样本鉴定蛋白数量
proteindata1 <- proteindata[,colnames(proteindata) %in% inputfiles$Sample]
peptidedata1 <- peptides[,colnames(peptides) %in% inputfiles$Sample]
Found <- colnames(proteindata)[colnames(proteindata) %in% inputfiles$Sample]
proteinnum <- matrix(nrow=length(Found), ncol=3)
for (i in 1:length(Found)) {
  proteinnum[i,1] <- Found[i]
  proteinnum[i,2] <- sum(proteindata1[,Found[i]] >0,na.rm = T)
  proteinnum[i,3] <- sum(peptidedata1[,Found[i]] >0,na.rm = T)
}
proteinnum <- as.data.frame(proteinnum)
colnames(proteinnum) <- c("样本名称","鉴定蛋白数目","鉴定肽段数量")

###总鉴定蛋白数量
resultdata <- matrix(nrow = 2, ncol = 2)
resultdata[c(1:2), 1] <- c("总蛋白数目","总肽段数目")
resultdata[c(1:2), 2] <- c(nrow(Proteins),nrow(peptides))
colnames(resultdata) <- c('Items', 'Number')
resultdata <- as.data.frame(resultdata)



###数据填充及差异统计
study <- readxl::read_xlsx("parameter.xlsx", sheet = 1)
sample_name <- readxl::read_xlsx("parameter.xlsx", sheet =  3)
search_parameter<-readxl::read_xlsx("parameter.xlsx", sheet =  4)
library(DMwR2)
library(ggplot2)
library(stringr)
library(openxlsx)
library(agricolae)
library(tidyr)
library(AnnotationDbi)
library(clusterProfiler)
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')


samplename <- readxl::read_xlsx("parameter.xlsx", 3)
samplename$No. <- as.character(samplename$No.)
samplegroup <- samplename[3:4]
colnames(samplegroup) <- c("SamplesName", "Groups")
samplegroup <- as.data.frame(samplegroup)
Compare <- readxl::read_xlsx("parameter.xlsx", 4)
Compare <- as.data.frame(Compare)
headers <- readxl::read_xlsx("parameter.xlsx", 7)

reportset <- readxl::read_xlsx("parameter.xlsx", sheet =  6)
FC <- as.numeric(reportset[3,2])
Org = reportset[4,2]
FC1 <- as.numeric(reportset[5,2])
UQPep <- as.numeric(reportset[6,2])

dir.create("3. 生信分析")
setwd("./3. 生信分析")
dir.create("1. 火山图分析")
dir.create("2. 聚类分析")
dir.create("3. GO分析")
dir.create("4. KEGG_PATHWAY分析")
setwd("../")


quan_data <- proteindata
quan_data[(length(colnames(quan_data))-length(inputfiles$Sample)+1):length(colnames(quan_data))] <- 
 lapply(quan_data[(length(colnames(quan_data))-length(inputfiles$Sample)+1):length(colnames(quan_data))],as.numeric)

norm <- colnames(quan_data)[(length(colnames(quan_data))-length(inputfiles$Sample)+1):length(colnames(quan_data))]

nalist <- which(rowSums(is.na(quan_data[,norm]))>0&rowSums(is.na(quan_data[,norm]))<length(quan_data[,norm]))
chisdata <- quan_data[nalist,]

nacol <- length(quan_data[, norm])
nodata <- quan_data[rowSums(is.na(quan_data[, norm])) == nacol,]
protein_quandata <- quan_data[rowSums(is.na(quan_data[, norm])) != nacol,] 
#mindata <- range(protein_quandata[, norm], na.rm = TRUE)
#mindata <- mindata[which.min(mindata)]
#
protein_quandata = protein_quandata[rowSums(is.na(protein_quandata[,norm])) < ncol(protein_quandata[,norm])*0.6,]
protein_quandata[,norm] <- knnImputation(protein_quandata[,norm],scale = T,k = 10)
# protein_quandata[, norm][is.na(protein_quandata[, norm])] <- 1000

imputedprotein <- createWorkbook()
addWorksheet(imputedprotein, sheetName = "Proteins")
setColWidths(imputedprotein, "Proteins", 
             cols = c(seq(from = 1, to = ncol(protein_quandata), 
                          by = 1)), widths = c(rep(14, ncol(protein_quandata))))
setRowHeights(imputedprotein, "Proteins", rows = 1, heights = 50)
writeData(imputedprotein, "Proteins", protein_quandata)
freezePane(imputedprotein, sheet = "Proteins", firstActiveRow = 2)
addStyle(imputedprotein, sheet = "Proteins", headerStyle1, rows = 1, cols = 1:5, gridExpand = TRUE)
addStyle(imputedprotein, sheet = "Proteins", headerStyle2, rows = 1, cols = 6:ncol(protein_quandata), gridExpand = TRUE)
addStyle(imputedprotein, sheet = "Proteins", datastyle, rows = 2:(nrow(protein_quandata)+1), cols = 1:ncol(protein_quandata), gridExpand = TRUE)
saveWorkbook(imputedprotein, "./2. 数据结果/1. 原始数据/Proteins_imputed.xlsx", overwrite = TRUE)



quan_dataA1 <- rep(0,nrow(protein_quandata))
quan_dataA2 <- rep(0,nrow(protein_quandata))
quan_dataB <- rep(0,nrow(protein_quandata))
conven_sample <- c('Protein.Group','Protein.Ids','Protein.Names','Genes','First.Protein.Description')

diffs <- matrix(nrow=nrow(Compare), ncol=4)
quandata <- createWorkbook()
DIFF <- createWorkbook()
chisqdata <- createWorkbook()

for (i in 1:nrow(Compare)) {
  
  
  groupA <- subset(samplegroup, select = "SamplesName", samplegroup$Groups ==  Compare[i, 1])
  groupB <- subset(samplegroup, select = "SamplesName", samplegroup$Groups ==  Compare[i, 2])
  type <- factor(c(rep(Compare[i, 1], nrow(groupA)), rep(Compare[i,2], nrow(groupB))), levels = c(Compare[i, 1], Compare[i,2]))
  
  Ratio <- rep(0, nrow(protein_quandata))
  meanA <- rep(0, nrow(protein_quandata))
  meanB <- rep(0, nrow(protein_quandata))
  
  if(nrow(groupA)>=3 && nrow(groupB)>=3){
    Pvalue <- rep(0, nrow(protein_quandata))
    
    for (j in 1:nrow(protein_quandata)) {
      if (is.na(mean(as.numeric(protein_quandata[c(groupA$SamplesName,groupB$SamplesName)][j, ])))) {
        Pvalue[j] = "NA"
      }else if(mean(as.numeric(protein_quandata[groupA$SamplesName][j,]))==mean(as.numeric(protein_quandata[groupB$SamplesName][j,]))){
        Pvalue[j] = 1
      }else{
        WT <- t.test(as.numeric(unlist(protein_quandata[c(groupA$SamplesName,groupB$SamplesName)][j,])) ~ type, exact = FALSE, 
                     alternative = "two.sided", paired = FALSE, 
                     var.equal = TRUE)
        Pvalue[j] <- WT$p.value
      }
    }
    Pvalue = as.numeric(Pvalue)
    for (t in 1:nrow(protein_quandata)) {
      meanA[t] <- mean(as.numeric(protein_quandata[t, groupA$SamplesName]))
      meanB[t] <- mean(as.numeric(protein_quandata[t, groupB$SamplesName]))
      Ratio[t] = meanA[t]/meanB[t]
    }
    #Ratio <- round(Ratio, 5)
    #Ratio <- signif(Ratio, 5)
    #Pvalue <- round(Pvalue, 5)
    #Pvalue <- signif(Pvalue, 5)
    
    quan_data <- cbind(protein_quandata[conven_sample], Ratio, Pvalue, meanA, meanB,protein_quandata[c(groupA$SamplesName,groupB$SamplesName)])
    colnames(quan_data) <- c(conven_sample,paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_Ratio"), 
                             paste0(Compare[i, 1], "_VS_",Compare[i, 2], "_Pvalue"),
                             paste0('Mean_',Compare[i, 1]),paste0('Mean_',Compare[i, 2]),
                             colnames(protein_quandata[c(groupA$SamplesName,groupB$SamplesName)]))
    
    quan_dataA1 <- cbind(quan_dataA1,quan_data[c(paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_Ratio"), 
                                                 paste0(Compare[i, 1], "_VS_",Compare[i, 2], "_Pvalue"))])
    
    quan_dataA2 <- cbind(quan_dataA2,quan_data[c(paste0('Mean_',Compare[i, 1]),paste0('Mean_',Compare[i, 2]))])
    
    
    
    volcano_data <- quan_data[c(paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_Ratio"), 
                                paste0(Compare[i, 1], "_VS_",Compare[i, 2], "_Pvalue"))]
    volcano_data$significant <- rep('no', nrow(volcano_data))
    volcano_data$significant[volcano_data[2]< 0.05 & 
                               volcano_data[1] >= FC] <- 'up'
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
    
    ggsave(file = paste0("./3. 生信分析/1. 火山图分析/", Compare[i, 1], "_VS_", Compare[i,2], "_volcano.png"), plot = volcanoplot,
           scale = 1, width = 12, height = 7.5, units = "in")
    ggsave(file = paste0("./3. 生信分析/1. 火山图分析/", Compare[i, 1], "_VS_", Compare[i, 2], "_volcano.pdf"), plot = volcanoplot,
           scale = 1, width = 12, height = 7.5, units = "in")
    
    pvalue <- paste0(Compare[i, 1], "_VS_",Compare[i, 2], "_Pvalue")
    ratio <- paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_Ratio")
    data_differ <- quan_data[which(quan_data[pvalue] < 0.05),]
    data_differ <- quan_data[which(((quan_data[ratio] > FC)|(quan_data[ratio] < 1/FC))&quan_data[pvalue] < 0.05),]
    data_differ <- data_differ[order(data_differ[,4], decreasing = TRUE),]
    
    diffs[i,1] <- paste0(Compare[i, 1], "_VS_", Compare[i, 2])
    diffs[i,2] <- nrow(data_differ[data_differ[paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_Ratio")]>1,])
    diffs[i,3] <- nrow(data_differ[data_differ[paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_Ratio")]<1,])
    diffs[i,4] <- nrow(data_differ)
    
    headerStyle1 <- createStyle(fontName = "Arial", fontSize = 10,fgFill = "#575757",
                                fontColour = "white", halign = "center", 
                                valign = "center", wrapText = TRUE)
    headerStyle2 <- createStyle(fontName = "Arial", fontSize = 10,fgFill = "#4682B4",
                                fontColour = "white", halign = "center", 
                                valign = "center", wrapText = TRUE)
    datastyle <- createStyle(fontName = "Arial", fontSize = 10, wrapText = FALSE)
    name <- paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_定量统计")
    
    addWorksheet(DIFF, sheetName = paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_定量统计"))
    setColWidths(DIFF, name,cols = c(seq(from = 1, to = ncol(data_differ), by = 1)), widths = c(rep(14, ncol(data_differ))))
    setRowHeights(DIFF, name, rows = 1, heights = 50)
    writeData(DIFF, name, data_differ)
    freezePane(DIFF, sheet = name, firstActiveRow = 2)
    addStyle(DIFF, sheet = name, headerStyle1, rows = 1, cols = 1:5, gridExpand = TRUE)
    addStyle(DIFF, sheet = name, headerStyle2, rows = 1, cols = 6:ncol(data_differ), gridExpand = TRUE)
    addStyle(DIFF, sheet = name, datastyle, rows = 2:nrow(data_differ), cols = 1:ncol(data_differ), gridExpand = TRUE)
    
    
    library(pheatmap)
    heatmap_data <- data_differ[c('Protein.Group',groupA$SamplesName,groupB$SamplesName)]
    rownames(heatmap_data) <- heatmap_data[,1]
    heatmap_data <- heatmap_data[,-1]
    #colnames(heatmap_data) <- stringr::str_sub(colnames(heatmap_data),start = 16)
    
    annotation_colA <- subset(samplegroup, samplegroup$Groups == Compare[i,1])
    annotation_colB <- subset(samplegroup, samplegroup$Groups == Compare[i,2])
    annotation_colX <- rbind(annotation_colA,annotation_colB)
    #annotation_colX$SamplesName <- stringr::str_sub(annotation_colX$SamplesName, start = 16)
    annotation_col <- data.frame(Group= annotation_colX$Groups)
    rownames(annotation_col) <- annotation_colX[,1]
    
    par(las = 2)
    heatmapplot <- pheatmap::pheatmap(heatmap_data, cluster_rows = TRUE, cluster_cols = T,
                                      scale = "row", annotation_col = annotation_col,
                                      color = colorRampPalette(rev(c("red", "white","blue")))(2000),
                                      show_rownames = F, show_colnames = T,border_color = NA,
                                      #cellwidth =48,#cellheight=4, border_color = NA,
                                      fontsize_row = 6.8, fontsize_col = 10.5, 
                                      legend = TRUE,legend_breaks = -1.5:1.5,
                                      main = paste0(Compare[i,1], "_VS_", Compare[i, 2], "_Cluster"))
    ggsave(file = paste0("./3. 生信分析/2. 聚类分析/",Compare[i, 1], "_VS_", Compare[i,2], "_Cluster.png"), plot = heatmapplot,
           scale = 1, width = 2.1*3,height = 3.2*3, units = "in")
    ggsave(file = paste0("./3. 生信分析/2. 聚类分析/",Compare[i,1], "_VS_", Compare[i, 2], "_Cluster.pdf"), plot = heatmapplot,
           scale = 1, width = 2.1*3,height = 3.2*3, units = "in")
    #units = "px", width = 210*3, height = 320*3, res = 52*3
    
    
    ggsave(file = "volcano.png", plot = volcanoplot,
           scale = 1, width = 12, height = 7.5, units = "in")
    ggsave(file = "Cluster.png", plot = heatmapplot,
           scale = 1, width = 2.1*3,height = 3.2*3, units = "in")
  }
  if(nrow(groupA)<3 || nrow(groupB)<3){
    
    for (t in 1:nrow(protein_quandata)) {
      meanA[t] <- mean(as.numeric(protein_quandata[t, groupA$SamplesName]))
      meanB[t] <- mean(as.numeric(protein_quandata[t, groupB$SamplesName]))
      Ratio[t] = meanA[t]/meanB[t]
    }
    Ratio <- round(Ratio, 5)
    Ratio <- signif(Ratio, 5)
    
    
    quan_data <- cbind(protein_quandata[conven_sample], Ratio,  meanA, meanB,protein_quandata[c(groupA$SamplesName,groupB$SamplesName)])
    colnames(quan_data) <- c(conven_sample,paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_Ratio"), 
                             paste0('Mean_',Compare[i, 1]),paste0('Mean_',Compare[i, 2]),
                             colnames(protein_quandata[Normalized]))
    
    quan_dataA1 <- cbind(quan_dataA1,quan_data[c(paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_Ratio") 
    )])
    
    quan_dataA2 <- cbind(quan_dataA2,quan_data[c(paste0('Mean_',Compare[i, 1]),paste0('Mean_',Compare[i, 2]))])
    
    ratio <- paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_Ratio")
    data_differ <- quan_data[which(quan_data["# Unique Peptides"] >= UQPep),]
    data_differ <- quan_data[which(((quan_data[ratio] > FC1)|(quan_data[ratio] < 1/FC1))&quan_data["# Unique Peptides"] >= UQPep),]
    data_differ <- data_differ[order(data_differ[,4], decreasing = TRUE),]
    
    
    diffs[i,1] <- paste0(Compare[i, 1], "_VS_", Compare[i, 2])
    diffs[i,2] <- nrow(data_differ[data_differ[paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_Ratio")]>1,])
    diffs[i,3] <- nrow(data_differ[data_differ[paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_Ratio")]<1,])
    diffs[i,4] <- nrow(data_differ)
    
    headerStyle1 <- createStyle(fontName = "Arial", fontSize = 10,fgFill = "#575757",
                                fontColour = "white", halign = "center", 
                                valign = "center", wrapText = TRUE)
    headerStyle2 <- createStyle(fontName = "Arial", fontSize = 10,fgFill = "#4682B4",
                                fontColour = "white", halign = "center", 
                                valign = "center", wrapText = TRUE)
    datastyle <- createStyle(fontName = "Arial", fontSize = 10, wrapText = FALSE)
    name <- paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_定量统计")
    
    addWorksheet(DIFF, sheetName = paste0(Compare[i, 1], "_VS_", Compare[i, 2], "_定量统计"))
    setColWidths(DIFF, name,cols = c(seq(from = 1, to = ncol(data_differ), by = 1)), widths = c(rep(14, ncol(data_differ))))
    setRowHeights(DIFF, name, rows = 1, heights = 50)
    writeData(DIFF, name, data_differ)
    freezePane(DIFF, sheet = name, firstActiveRow = 2)
    addStyle(DIFF, sheet = name, headerStyle1, rows = 1, cols = 1:5, gridExpand = TRUE)
    addStyle(DIFF, sheet = name, headerStyle2, rows = 1, cols = 6:ncol(data_differ), gridExpand = TRUE)
    addStyle(DIFF, sheet = name, datastyle, rows = 2:nrow(data_differ), cols = 1:ncol(data_differ), gridExpand = TRUE)
    
    
  }
  
  if(Org != FALSE) {
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
    }else if (Org == "pig") {
      Org = "ssc"
      print(Org)
      library(org.Ss.eg.db)
      OrgDB = 'org.Ss.eg.db'
    }
    analysname <- paste0(Compare[i, 1], "_VS_", Compare[i, 2])
    ID <- data_differ$Genes
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
    
    
    kegg0 <- enrichKEGG(gene = ID[,2], organism = Org,
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
    #keggplot <- dotplot(kegg, showCategory = 15, title = "KEGG_Pathway", font.size = 9)+
    #  scale_colour_gradient(low = "red", high = "green")
    ggsave(file = paste0("./3. 生信分析/4. KEGG_PATHWAY分析/",analysname, "_KEGG_Pathway.png"),keggplot,
           width = pic_width, height = pic_height, dpi = 600,
           device = "png")
    ggsave(file = paste0("./3. 生信分析/4. KEGG_PATHWAY分析/",analysname, "_KEGG_Pathway.pdf"),keggplot,
           width = pic_width, height = pic_height, dpi = 600,
           device = "pdf")
    
  }
  
}
colnames(diffs) <- c("组间比较","上调差异","下调差异","汇总")
diffs <- as.data.frame(diffs)

setwd("./2. 数据结果/")
dir.create("2. 差异统计")

headerstyle3 <- createStyle(fontName = "Arial",fontSize = 11, fgFill = "#808080",
                            fontColour = "white", halign = "center",
                            valign = "center", wrapText = TRUE)
datastyle <- createStyle(fontName = "Arial", fontSize = 10, wrapText = FALSE)
#addWorksheet(DIFF, sheetName = "表头说明")
#setRowHeights(DIFF, "表头说明", row = 1:c(nrow(headers)+1),heights = 28)
#setColWidths(DIFF, "表头说明",
#             cols = c(seq(from = 1, to = ncol(headers),
#                          by = 1)), widths = 28)

#freezePane(DIFF, sheet = "表头说明", firstActiveRow = 2)
#writeData(DIFF, "表头说明", headers)
#addStyle(DIFF,"表头说明", headerstyle3, rows = 1,cols = 1:3, gridExpand = TRUE)
saveWorkbook(DIFF, paste0("./2. 差异统计/", "2. 差异定量结果.xlsx"),overwrite = TRUE)
quan_dataA2 <- quan_dataA2[,unique(names(quan_dataA2), with = FALSE)]
quan_dataC <- cbind(quan_data[conven_sample],quan_dataA1, quan_dataA2, protein_quandata[norm])
quan_dataD <- subset(quan_dataC, select = -c(quan_dataA1, quan_dataA2, quan_dataB))

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

headerstyle3 <- createStyle(fontName = "Arial",fontSize = 11, fgFill = "#808080",
                            fontColour = "white", halign = "center",
                            valign = "center", wrapText = TRUE)
datastyle <- createStyle(fontName = "Arial", fontSize = 10, wrapText = FALSE)
addWorksheet(quandata, sheetName = "表头说明")
setRowHeights(quandata, "表头说明", row = 1:c(nrow(headers)+1),heights = 28)
setColWidths(quandata, "表头说明",
             cols = c(seq(from = 1, to = ncol(headers),
                          by = 1)), widths = 28)

freezePane(quandata, sheet = "表头说明", firstActiveRow = 2)
writeData(quandata, "表头说明", headers)
addStyle(quandata,"表头说明", headerstyle3, rows = 1,cols = 1:3, gridExpand = TRUE)
saveWorkbook(quandata, paste0("./2. 差异统计/", "1. 蛋白定量统计.xlsx"),overwrite = TRUE)

#dir.create("3. 卡方检验")
#setwd("../")
#saveWorkbook(chisqdata, paste0("./2. 数据结果/3. 卡方检验/", "CHISQ_TEST_RESULT.xlsx"),overwrite = TRUE)




if(Org != FALSE) {
  ggsave(file = "~/Software/data/result/GO.png", 
         plot = p, width = pic_width, height = pic_height,
         dpi = 600, device = "png")
  ggsave(file = "~/Software/data/result/KEGG_Pathway.png",plot = keggplot,
         width = pic_width, height = pic_height, dpi = 600,
         device = "png")
}
if(Org == FALSE) {
  file.copy("~/Software/code/GO.png","~/Software/data/result/GO.png")
  file.copy("~/Software/code/KEGG_Pathway.png","~/Software/data/result/KEGG_Pathway.png")
}
