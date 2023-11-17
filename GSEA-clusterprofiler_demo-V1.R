# BiocManager::install("org.Rn.eg.db")
####1、导入需要的包####
library(clusterProfiler)
library(org.Hs.eg.db)#此处为人的库包,小鼠为：org.Mm.eg.db,大鼠为：org.Rn.eg.db
library(ggplot2)
library(openxlsx)
library(org.Mm.eg.db)
#R语言做GSEA的优势

####1.1设置物种等信息####
organism='hsa' #小鼠：mmu，大鼠：rno
OrgDb="org.Hs.eg.db"#人 org.Hs.eg.db #小鼠为:org.Mm.eg.db,大鼠为:org.Rn.eg.db
species="Homo sapiens"# 人：Homo sapiens，小鼠:Mus musculus,大鼠:Rattus norvegicus


#1.GSEA
###2、设置路径，读取数据####
#（注意每个人电脑的路径有所差异）
getwd()
setwd("D:/Desktop/鞠文浩GSEA/")

dir.create("GSEA富集分析")


cluster <- read.xlsx("./1.蛋白定量统计.xlsx",1)
setwd("GSEA富集分析/")
# cluster = cluster[,c(2,11)]
####3、排序####
colnames(cluster)[1] = "gene"
FCgenelist = cluster$实验组_VS_对照组_FC
names(FCgenelist)  = cluster$gene
FCgenelist <- sort(FCgenelist,decreasing=T) #decreasing order
head(FCgenelist)

####4、GO-BP数据库的GSEA分析####
library(msigdbr)
#### 4.1、GO物种为人的库下载####
Human_msigdbr <- msigdbr(species=species)
#### 4.2、GO的GSEA分析####
#设置随机种子确保可重复性
set.seed(123)
egseGO <- gseGO(FCgenelist, OrgDb=OrgDb,
                ont='BP',keyType="SYMBOL",
                nPerm=1000, minGSSize=20, maxGSSize=500,pvalueCutoff = 1,
                pAdjustMethod = "BH", verbose = TRUE, by="fgsea")


#查看结果
head(egseGO,1);dim(egseGO)
head(data.frame(egseGO$ID,egseGO$Description))
enrichplot::goplot(egseGO)


####4.3、导出GO的GSEA结果####
write.csv(egseGO,file = "GOBP_GSEA_.csv",row.names = F)


####4.4导出GO的GSEA图结果####
dir.create("GO_GSEA结果图片")
setwd("./GO_GSEA结果图片")
egseGO@result$Description = gsub(pattern = "/",
                                 replacement = "_",
                                 x = egseGO@result$Description)
#选择导出数量
#1 导出所有显著的
topn_GO = nrow(egseGO@result[egseGO@result$p.adjust<0.05,])
#2导出top N
topn_GO =20

for (i in c(1:topn_GO)) {
  p <- enrichplot::gseaplot2(egseGO ,geneSetID=i, pvalue_table=F,base_size = 11,
                 rel_heights = c(5,1.2,1.2),subplots= 1:3,
                 title = paste0(egseGO@result$Description[i]),
                 ES_geom = "line")
  ggsave(p,filename = paste0(egseGO@result$Description[i],".pdf"),
         width = 7,height = 4.8)
}

####5、KEGG分析####
####5.1基因名称转化#####
gene.tx <- bitr(names(FCgenelist),fromType="SYMBOL",toType=c("ENTREZID"),
                OrgDb = OrgDb)
colnames(gene.tx)[1] <- "gene"
gene.tx <- merge(gene.tx,cluster,by="gene")
FCgenelist <- gene.tx$实验组_VS_对照组_FC
names(FCgenelist) <- as.character(gene.tx$ENTREZID) #named vector
#排序
FCgenelist <- sort(FCgenelist,decreasing=T) #decreasing orde
head(FCgenelist)

####5.2、KEGG数据库的GSEA分析######
set.seed(111)
egseKEGG <- gseKEGG(FCgenelist,organism=organism,keyType="ncbi-geneid",
                    nPerm=1000, minGSSize=5, maxGSSize=500,
                    pvalueCutoff=1, pAdjustMethod = "BH",verbose = TRUE)

kegg2<- setReadable(egseKEGG, OrgDb = OrgDb, keyType="ENTREZID")
head(kegg2)
####5.3、导出KEGG的GSEA结果####
setwd("../")
write.csv(kegg2,file = "KEGG_GSEA_lixin.csv")
# #install.packages("ggridges")
library(ggridges)
#install.packages("ggupset")
library(ggupset)
enrichplot::upsetplot(egseKEGG)
pkegg = enrichplot::ridgeplot(egseKEGG,showCategory = 10,fill = "pvalue")
ggsave("./KEGG_GSEA_ridgeplot.png",plot = pkegg,width = 7,height = 6.5,dpi = 600)
ggsave("./KEGG_GSEA_ridgeplot.pdf",plot = pkegg,width = 7,height = 6.5,dpi = 600)
####5.4、导出GSEA图#####
#设定导出路径，注意如果选择导出的数量较大，建议新建文件夹作为导出的路径
#不建议直接导出到桌面
getwd()
dir.create("KEGG_GSEA结果图片")
setwd("./KEGG_GSEA结果图片")

####5.5、设置导出的GSEA富集到的前多少通路#####
topn_KEGG = nrow(kegg2@result[kegg2@result$pvalue<0.05,])
#替换"/"为"_",通路名称中个别包含"/",保存时无法作为文件名
kegg2@result$Description = gsub(pattern = "/","_",
                                x = kegg2@result$Description)

#导出图片
for (i in c(1:topn_KEGG)) {
  p <- enrichplot::gseaplot2(kegg2,geneSetID=i, pvalue_table=F,base_size = 11,
                 rel_heights = c(5,1.2,1.2),subplots= 1:3,
                 title = paste0(kegg2@result$Description[i]),
                 ES_geom = "line")
  ggsave(p,filename = paste0(kegg2@result$Description[i],".pdf"),
         width = 7,height = 4.8)
}






#####wiki pathway####
set.seed(111)
gse_wiki = clusterProfiler::gseWP(FCgenelist,organism=species,
                                  minGSSize=5, maxGSSize=500,
                                  pvalueCutoff=1, 
                                  pAdjustMethod = "BH",
                                  verbose = TRUE)


gse_wiki2 <- setReadable(gse_wiki, OrgDb = OrgDb, keyType="ENTREZID")
head(gse_wiki2)
####5.3、导出KEGG的GSEA结果####
setwd("../")
write.csv(gse_wiki2,file = "wiki_GSEA_lixin.csv")
# #install.packages("ggridges")

pkegg = enrichplot::ridgeplot(gse_wiki2,showCategory = 15,fill = "pvalue")
pkegg
ggsave("./wiki_GSEA_ridgeplot.png",plot = pkegg,width = 7,height = 6.5,dpi = 600)
ggsave("./wiki_GSEA_ridgeplot.pdf",plot = pkegg,width = 7,height = 6.5,dpi = 600)


pdot = enrichplot::dotplot(gse_wiki2,showCategory = 15,color = "pvalue",
                           title="                     Total\n WiKipathway Enrichment Top 15")
pdot
ggsave("./GSEA富集分析/GSEA_WIKI/wiki_GSEA_dotplot.png",plot = pdot,width = 7,height = 6.5,dpi = 600)
ggsave("./GSEA富集分析/GSEA_WIKI/wiki_GSEA_dotplot.pdf",plot = pdot,width = 7,height = 6.5,dpi = 600)

####5.4、导出GSEA图#####
#设定导出路径，注意如果选择导出的数量较大，建议新建文件夹作为导出的路径
#不建议直接导出到桌面
getwd()
dir.create("wiki_GSEA结果图片")
setwd("./wiki_GSEA结果图片")

####5.5、设置导出的GSEA富集到的前多少通路#####
topn_wiki = nrow(gse_wiki2@result[gse_wiki2@result$pvalue<0.05,])
topn_wiki = 20
#替换"/"为"_",通路名称中个别包含"/",保存时无法作为文件名
gse_wiki2@result$Description = gsub(pattern = "/","_",
                                x = gse_wiki2@result$Description)

#导出图片
for (i in c(1:topn_wiki)) {
  p <- enrichplot::gseaplot2(gse_wiki2,geneSetID=i, pvalue_table=F,base_size = 11,
                             rel_heights = c(5,1.2,1.2),subplots= 1:3,
                             title = paste0(gse_wiki2@result$Description[i]),
                             ES_geom = "line")
  ggsave(p,filename = paste0(gse_wiki2@result$Description[i],".pdf"),
         width = 7,height = 4.8)
}






