#导入R包
# install.packages("pheatmap")
library(pheatmap)
# install.packages("openxlsx")
library(openxlsx)
#读入数据
library(ggplot2)
getwd()
setwd("D:/Desktop/")

data = read.xlsx("2.差异定量结果.xlsx",2,rowNames = T)

str(data)
View(data)




#进行列明即样本分组注释
annotation_col = data.frame(group = c(rep("Exp",5),rep("Con",5)))
#annotation_col = data.frame(group = c(rep("TBM",15),rep("CM",8),rep("NI",5)))
                         #   ,type= c(rep("Cell Lysates",9),rep("Cell Supernatant",9)))
row.names(annotation_col) = colnames(data)
annotation_col$group=factor(annotation_col$group,levels = unique(annotation_col$group))
p=pheatmap(data,#data1[newOrder[,"Cluster"]==3,]
         #display_numbers = T,
         cluster_rows = T,
         cluster_cols=F,
         scale="row",
         annotation_col = annotation_col,
         #annotation_row = annotation_row,
         #annotation_row = NA,
         #annotation_colors = annotation_colors,
         color=colorRampPalette(rev(c("red","white","blue")))(20000),
         show_rownames=T,show_colnames=T,
         cellwidth =40,cellheight=5,
         border_color =NA,
         fontsize =14,fontsize_row = 3,fontsize_col =15,fontsize_number=4,
         clustering_method='ward.D2',
         legend = TRUE,na_col = "grey",drop_levels=T,angle_col = 45,
         annotation_legend = T,#cutree_rows = 3,
         main="")
p
ggsave(filename = "聚类热图.png",p,dpi = 600,width =8,height =10)
pdf(file = "聚类热图.pdf",width =8,height =10)
print(p)
dev.off()

