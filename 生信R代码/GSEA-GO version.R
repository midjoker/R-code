library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(openxlsx)
setwd("D:/Desktop/何晶GSEA/细胞上清/A_VS_B")
data=read.xlsx("D:/Desktop/何晶GSEA/细胞上清/1. 蛋白定量统计.xlsx",1)
#1.GSEA
cluster=data
colnames(cluster)[1] <- "gene"
FCgenelist <- cluster$`LOG2(A_VS_B_Ratio)`
names(FCgenelist) <- as.character(cluster$Gene.Symbol)
FCgenelist <- sort(FCgenelist,decreasing=T) #decreasing order
head(FCgenelist)
library(msigdbr)
Human_msigdbr <- msigdbr(species="Mus musculus")
egseGO <- gseGO(FCgenelist, OrgDb=org.Mm.eg.db,
                ont='BP',keyType="SYMBOL",
                nPerm=1000, minGSSize=20, maxGSSize=500,pvalueCutoff = 1,
                pAdjustMethod = "BH", verbose = TRUE, by="fgsea")
# head(egseGO,1);dim(egseGO)
# head(data.frame(egseGO$ID,egseGO$Description))
#egseGO<- setReadable(egseGO, OrgDb = org.Rn.eg.db, keyType="ENTREZID")
write.csv(egseGO,file = "GSEA-GO富集结果.csv")
#

# #
library(enrichplot)
upsetplot(egseGO)
p<-ridgeplot(egseGO,showCategory = 20)+
  theme(axis.text.y = element_text(size = 10, color="black"))
print(p)
ggsave(filename = "GSEA-GO-Top20展示.png",p,dpi = 600,width = 12,height =10)
p <- gseaplot2(egseGO,geneSetID=1:15, pvalue_table=F,base_size = 10,rel_heights = c(5,0.5,0.5),
                ES_geom = "line")
print(p)
pdf(file = "GSEA-GO-linetopn.pdf",width = 9,height = 9)
print(p)
dev.off()

for (i in c(1:nrow(egseGO@result))) {
  p <- gseaplot2(egseGO,geneSetID=i, pvalue_table=F,base_size = 11,
                 rel_heights = c(5,1.2,1.2),subplots= 1:3,
                 ES_geom = "line")
  pdf(file = paste0("./GSEA-GO-picture/",egseGO@result$Description[i],".pdf"),width = 9,height = 6)
  print(p)
  dev.off()
}




# 2. ORA
DEGdata<-cluster[c(1:28),]
DEGdata2 <- subset(DEGdata,DEGdata$logFC<c(-1))
DEgenelist <- as.character(DEGdata2$gene)
library(org.Hs.eg.db)
egoBP <- enrichGO(DEgenelist, OrgDb=org.Hs.eg.db, ont='BP',
                                  pAdjustMethod='BH', pvalueCutoff=0.05, 
                                   minGSSize=5,keyType='SYMBOL')
write.csv(egoBP,file = "~/Documents/尿液外泌体/ORA_mp_vs_hc_down.csv")
#
gene.df <- bitr(DEgenelist,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                               OrgDb = org.Hs.eg.db)
gene.kegg <- bitr_kegg(gene.df$ENTREZID,fromType="ncbi-geneid",
                                            toType="kegg",organism='hsa')
ekegg <- enrichKEGG(gene.df$ENTREZID, organism='hsa',keyType="ncbi-geneid",
                                      pvalueCutoff=0.05,pAdjustMethod='BH',
                                    minGSSize=5,maxGSSize=500,use_internal_data=F)
ekeggx <- setReadable(ekegg,'org.Hs.eg.db','ENTREZID')
write.csv(ekeggx,file = "~/Documents/尿液外泌体/ORA_kegg_vs_hc_down.csv")
#
library(enrichplot)
dotplot(egoBP,showCategory=5)
barplot(ekegg)


