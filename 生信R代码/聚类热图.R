#导入R包
#install.packages("pheatmap")
rm(list = ls())
library(ggplot2)
library(pheatmap)
#install.packages("openxlsx")
library(openxlsx)
#读入数据
getwd()
setwd("D:/Desktop/赵玉琼老师售后分析")

data = read.xlsx("并集_80个蛋白.xlsx",4,rowNames = T)
str(data)
View(data)

# rownames(Gene names) = Gene names[,1]　　　#给data1添加列名
# data1 <- data1[,-1]　　　　　　#去除第一列的名称
# data1 <- as.matrix(data1)



#直接画图
pheatmap(data)



#进行列明即样本分组注释
annotation_col = data.frame(group = c(rep("2DM+AS",4),rep("AS",4),rep("Control",4)))
row.names(annotation_col) = colnames(data)

#,rep("敲减_Phos/ Unphos
#",1),rep("FC≥1.5 敲减_vs_对照_Ratio
#",1)))
#row.names(annotation_col) = colnames(data)

p=pheatmap(data,#data1[newOrder[,"Cluster"]==3,]
         #display_numbers = T,
         cluster_rows = T,
         cluster_cols=T,
         scale="row",
         annotation_col = annotation_col,
         #annotation_row = annotation_row,
         annotation_row = NA,
         #annotation_colors = annotation_colors,
         color=colorRampPalette(rev(c("red","white","blue")))(20000),
         show_rownames=F,show_colnames=T,
         cellwidth =20,cellheight=7,
         border_color ="NA",
         fontsize =12,fontsize_row = 7,fontsize_col =10,fontsize_number=12,
         clustering_method='single',
         legend = TRUE,na_col = "grey",drop_levels=T,angle_col = 45,
         annotation_legend = T,#cutree_rows = 3,
         main="")
p
ggsave(filename = "热图.png",p,dpi = 600,width = 10,height =10 )
pdf(file = "热图.pdf",width = 10,height = 10)
print(p)
dev.off()
