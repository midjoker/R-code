##气泡图
# install.packages("egg")
#install.packages("ggplot2")
library(ggplot2)
library(openxlsx)
library(egg)
library(grid)

setwd("C:/Users/wayen/Desktop")
#读取数据
rm(list = ls())
workbook<-"BD_VS_NC_KEGG_Pathway.xlsx"
data<-read.xlsx(workbook,1,check.names=F)
#data <- read.table("1185_KEGG_RichFactor.txt",header = T,sep = "/t",encoding = "UTF-8")
View(data)
str(data)
data = data[1:15,]
#绘制气泡图
p = ggplot(data=data,
       aes(Count,y =reorder(Description,Count)))+ #定义绘图数据，横纵坐标， reorder(Term,FoldEnrichment))
    geom_point(aes(size=Count,          #加气泡，大小根据count,颜色根据pvalue
                 color=-log10(pvalue)))+      
    scale_colour_gradientn(colours = c("red","yellow","green"))+  #定义颜色
    # scale_x_continuous(limits = c(1,15),breaks = c(5,10,15))+ #设置x轴范围
    labs(size="Count",            #加图例以及标题
       color="-log10pvalue",
       x="GO enrichment",
       y="Term",
       title="")+
  # coord_fixed(ratio=1/1)+ #横坐标及纵坐标的比例
    theme_bw()+#主题 
    scale_colour_gradient(low = "green", high = "red")+
  #facet_wrap( ~type , scales = "free", ncol = 1,strip.position = "right")+
 # scale_y_discrete(labels = function(x) str_wrap(x,width = 60))
     theme(axis.text.x = element_text(face="plain",size=11),
           axis.text.y = element_text(face="plain",size=11),
           axis.title = element_text(face="plain",size=12),
          plot.margin=unit(c(1,1,1,1),"cm"))


p #p为ggplot绘制出来的图形的对象

##柱形图
str(data)
p = ggplot(data=data,
           aes(x=Count,y =reorder(Description,pvalue),
               fill=-log10(pvalue)))+ #定义绘图数据，横纵坐标， reorder(Term,FoldEnrichment))
  geom_bar(stat = "identity")+  #加气泡，大小根据count,颜色根据pvalue
  # scale_x_continuous(limits = c(1,15),breaks = c(5,10,15))+ #设置x轴范围
  labs(size="Count",            #加图例以及标题
       color="-log10pvalue",
       x="Count",
       y="Term",
       title="Top 15 of KEGG")+
  coord_flip()+#guides(fill=F)+
  geom_text(data=data,aes(x=Count+0.2,y =reorder(Description,Count)),size=3.3,
            label=format(data$pvalue, digits = 2))+
  #coord_fixed(ratio=1/1)+ #横坐标及纵坐标的比例
  theme_bw()+theme(axis.text.x = element_text(size = 12,angle = 80,hjust = 1))#主题 
  #scale_y_discrete(labels = function(x) str_wrap(x,width = 60))
 # facet_wrap( ~type , scales = "free_x", ncol = 3)


p
# colors =  c("red","blue")
p 
  #geom_area(position='fill')+scale_fill_gradient( low = "blue",high = "red")









library(egg)
library(grid)
#设置图形坐标轴大小
p1 = set_panel_size(p,#对ggplot绘制的对象p进行调整
                    width  = unit(8, "cm"),
                    height = unit(20, "cm"))
#绘制形成设置坐标大小后的图形
grid.draw(p1)


