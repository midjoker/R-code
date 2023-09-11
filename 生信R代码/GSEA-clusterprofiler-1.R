library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(openxlsx)
setwd("D:/Desktop/何晶GSEA/蛋白溶液")
data=read.xlsx("D:/Desktop/何晶GSEA/蛋白溶液/1. 蛋白定量统计-蛋白溶液.xlsx",1)
#1.GSEA
cluster=data
colnames(cluster)[1] <- "gene"
FCgenelist <- cluster$`log2(1_VS_2_Ratio)` #numeric vector
names(FCgenelist) <- as.character(cluster$Gene.Symbol)
FCgenelist <- sort(FCgenelist,decreasing=T) #decreasing order
# head(FCgenelist)
library(msigdbr)

#
gene.tx <- bitr(names(FCgenelist),fromType="SYMBOL",toType=c("ENTREZID"),
                OrgDb = org.Mm.eg.db)
colnames(gene.tx)[1] <- "gene"
gene.tx <- merge(gene.tx,cluster,by="gene")
FCgenelist <- gene.tx$`log2(1_VS_2_Ratio)`
names(FCgenelist) <- as.character(gene.tx$ENTREZID) #named vector
FCgenelist <- sort(FCgenelist,decreasing=T) #decreasing orde
# head(FCgenelist)
egseKEGG <- gseKEGG(FCgenelist,organism='mus',keyType="ncbi-geneid",
                    nPerm=1000, minGSSize=5, maxGSSize=500,
                    pvalueCutoff=1, pAdjustMethod = "BH",verbose = TRUE)

kegg<- setReadable(egseKEGG, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
write.csv(kegg,file = "GSEA-KEGG富集结果.csv")
#
library(enrichplot)
upsetplot(egseKEGG)

egseKEGG@result$Description = gsub(pattern = "/",replacement = "_",
                                   x = egseKEGG@result$Description)
p<-ridgeplot(egseKEGG,showCategory = 20)+
  theme(axis.text.y = element_text(size = 10, color="black"))
print(p)
ggsave(filename = "GSEA-picture.png",p,dpi = 600,width = 8,height =10)
# p <- gseaplot2(egseGO,geneSetID=1:15, pvalue_table=F,base_size = 10,rel_heights = c(5,0.5,0.5),
#                 ES_geom = "line")
# print(p)
for (i in c(1:nrow(egseKEGG@result))) {
  p <- gseaplot2(egseKEGG,geneSetID=i, pvalue_table=F,base_size = 11,
                 rel_heights = c(5,1.2,1.2),subplots= 1:3,
                 ES_geom = "line")
  pdf(file = paste0("./GSEA-picture/",egseKEGG@result$Description[i],".pdf"),width = 9,height = 6)
  print(p)
  dev.off()
  
}

