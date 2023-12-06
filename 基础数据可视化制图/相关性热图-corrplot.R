library(corrplot)
library(openxlsx)
setwd("D:/desktop/姚香平老师售后处理-V2/姚香平老师-v2 （文章）/第1步/")
#读取数据
data =read.xlsx("数据提取-剔除5样本.xlsx",sheet = 2,rowNames = T)
View(data)
data1=data[,24:28]
data1=t(data1)
data1=data1[,c(-5)]
#计算相关性系数
cor_data = cor(data1)
#corrplot包函数cor.mtest()函数计算p值：
p.mat <-cor.mtest(data1, conf.level= .95)
#这种方法计算得到的是一个p值的列表：
# class(pdata)
## [1] "list"
# head(pdata,1)

#pheatmap相关性热图
pheatmap(cor_data)

# corrplot相关性热图
p=p.cor <- corrplot(corr = cor_data, # 相关性矩阵
                  method = "square", # 矩形c("circle", "square", "ellipse", "number", "shade", "color", "pie"),
                  type = "full", # 上三角
                  my_color_palette <- colorRampPalette(c("blue", "white", "red"))(100),
                  order = "AOE", # 特征向量角序
                  tl.col = "black", # 标签颜色
                  p.mat = p.mat$p,sig.level = c(0.001,0.01,0.05), # 设置显著性p值为0.05
                  insig = "label_sig", # blank(不显著的设置为空白)
                  pch.cex = 0.9,pch.col = "white" # 设置显著性符号的颜色
)
p
ggsave(file = "姚香平老师-脑脊液样本-5NI-相关性热图.png",p,width =10,height = 10)
pdf(file = "姚香平老师-脑脊液样本-8C7M-相关性热图.pdf",width =10,height = 10)

print(p) 
dev.off()
library(xlsx)
write.xlsx(x = cor_data,sheetName = 'cor_data',file = "cor_data.xlsx")
write.xlsx(x = p.mat$p,sheetName = 'p_data',file = "cor_data.xlsx",append = T)



