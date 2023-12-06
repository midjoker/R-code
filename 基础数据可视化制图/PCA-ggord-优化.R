# install.packages("pca3d")
# install.packages('devtools')
# library(devtools)
# install_github('fawda123/ggord')
# library(ggord)

library(stats)
library(pca3d)
library(rgl)
library(openxlsx)
library(ggplot2)
library(ggord)
library(ggrepel)
library(DMwR)
rm(list = ls())
setwd("D:/desktop/丁婷婷/2023.9.21-HWYC04120230710-北京大学第一医院-丁婷婷老师-Label free定量蛋白质谱项目报告-ymm/2. 数据结果/2. 差异统计/")
#setwd("C:/Users/wayen/Desktop/吴艳???-相关性分???/菌门_vs_差异代谢???/主成分分???/PCA数据")
#数据导入
workbook<-"1.蛋白定量统计.xlsx"
# dat<-read.xlsx(workbook,4,check.names=F,encoding = "UTF-8",row.names=1)


dat<-read.xlsx(workbook,3,rowNames = T)#
# dat<-read.csv(workbook,check.names=F,row.names=1)
# #dat <- read.table("PCA-I.txt",header = T,sep = "\t",check.names =F)
# dat = read.csv("PCA.csv",header=T,sep=",",check.names = F,encoding = "UTF-8")
View(dat)
# row.names(dat) = dat[,1]
data <- dat
#View(data)
# data1 = data[rowSums(is.na(data)) < ncol(data)-1,]
# dim(data1)
# knnOutput <- knnImputation(data1,scale = T,k = 3)  # 使用KNN插???,k为选取的最近邻样本的个???
#主成分计
str(data)
pca<-prcomp(t(data),center=T,scale=T,retx = T)
#pca<-prcomp(knnOutput,center=T,scale=T,retx = T)
#pca <- princomp(data,cor = T,scores = T)
#查看
summary(pca)
#提取主成
scores=as.data.frame(pca$x)
scores
#选取2个主成分
newdata2<-scores[,c(1:2)]
View(newdata2)
#write.xlsx(newdata2,"PCA-score.xlsx")
#去读分组数据，同一个excel的第二个sheet
#group <-read.xlsx(workbook,2,check.names=F,encoding = "UTF-8")
#group = dat[,1]

group = data.frame(group = c(rep("HC",6),rep("LN",6),rep("SLE",6)))

#group = c(rep("HT",27),rep("HC",26))
#合并组别及主成分数据
total1<-cbind(newdata2,group)
total1
#row.names(total1) <- colnames(data)
#查看后数据情???
View(total1)
#关联total1，方便条用里面数???
attach(total1)

#绘图

p <- ggord(pca,grp_in= as.factor(total1$group),txt =NULL,ellipse=T,
           arrow=0,alpha_el = 0.4,size = 3,veclsz =1,vec_ext=0,)
p

p1 <- p+geom_text_repel(aes(PC1, PC2, segment.size = .3,label=rownames(total1),color=group$group), 
                        force_pull = 0.6,force = 0.7,max.overlaps = 30,
                         fontface="bold", max.iter = 5e8,
                         box.padding=unit(0.25, "lines"), point.padding=unit(0.25, "lines"),max.time=5)+ theme_classic(base_size = 16)

p1+theme_bw()


p1 <- p+geom_label_repel(aes(PC1, PC2, fill=factor(group$group),  segment.size = .3,
                       label=rownames(total1)), force_pull = 0.6,force = 0.7,max.overlaps = 30,
                       fontface="bold", color="white", max.iter = 5e8,
                   box.padding=unit(0.25, "lines"), point.padding=unit(0.25, "lines"),max.time=5,
                   segment.colour = "grey50")+ theme_classic(base_size = 16)
p1
ggsave(file = "PCA图.png",p1,width =10,height = 8)
pdf(file = "PCA图.pdf",width =10,height = 8)

print(p1) # 打印图片
dev.off()
# p+geom_text(aes(label=row.names(newdata2)),size=4,hjust = 0, nudge_x = 0.05,check_overlap = T)+
#   labs(title = "PCA plot")+
#   theme_bw()+theme(plot.title=element_text(hjust=19.5),
#                    text = element_text(size = 15,color="black"),
#                    axis.text.x = element_text(size = 16,color="black"),
#                    axis.text.y = element_text(size = 16,color="black"),
#                    axis.title.x = element_text(size = 17,color="black"),
#                    axis.title.y = element_text(size = 17,color="black"),
#                    legend.title = element_text(size = 17,color="black"),
#                    legend.text=element_text(size=17),
#                    legend.key.width = unit(0.8, "cm"),
#                    legend.key.height = unit(0.8, "cm"))
dev.off()
#3D PCA
pca3d(pca, group=group$group, show.ellipses=TRUE, show.plane=F)

library(rgl)
colour<-c(rep(2,6),rep(7,6),rep(12,6),rep(17,6),rep(22,6))
pca3d(pca,col=colour,group=group$group,show.ellipses = T,
      show.plane=T,show.labels = T)



