
rm (list = ls ())
library(ggplot2)
library(openxlsx)
library(ggrepel)
library(dplyr)
setwd("D:/Desktop/")
#windowsFonts(Ariar = windowsFont("Arial"))
workbook<-"1.蛋白定量统计.xlsx"

data<-read.xlsx(workbook,1,rowNames = F)
#workbook1<-"新建 Microsoft Excel 工作表.xlsx"
#data1<-read.xlsx(workbook1, 1, rowNames=F)

#全转录组替换负无穷
#data[is.na(data$log2FC),1] = -Inf
#data$log2FC = as.numeric(data$log2FC)

#注释掉的这段代码实际是可以不用手动输入significant的，可以通过这段代码直接生成。
# data$signifcant = ifelse(data$B_VS_C_Pvalue < 0.05 & abs(log2(data$B_VS_C_Ratio)) >= log2(1.2), 
#                       ifelse(data$B_VS_C_Ratio >= 1.2 ,'Up','Down'),'No')
  data$signifcant = ifelse(data$LN_VS_SLE_Pvalue < 0.05 & abs(log2(data$LN_VS_SLE_FC)) >= log2(2),
                         ifelse(log2(data$LN_VS_SLE_FC)>= log2(2),'Up','Down'),'No')


# data$signifcant = ifelse(data$pvalue < 0.05 & abs(data$LOG2FC) >= 1,
#                          ifelse(data$LOG2FC >= 1 ,'Up','Down'),'No')


table(data$signifcant)
str(data)
 p <- ggplot( data = data, 
              aes(x = log2(LN_VS_SLE_FC), 
                  y = -log10(LN_VS_SLE_Pvalue), 
                 colour=signifcant)) +#                 label = Name
         geom_point(alpha=0.5, size=2) +
         scale_color_manual(values=c("green", "grey","red"))+
         xlim(c(-7, 7)) +ylim(c(0,7))+geom_vline(xintercept=c(-log(2,2),log(2,2)),lty=4,col="grey",lwd=0.8) +
         geom_hline(yintercept = -log10(0.05),lty=5,col="grey",lwd=0.8) +
         labs(x="log2(Fold Change)",
              y="-log10 (Pvalue)",
              title="")  + 
         theme_bw()+theme(panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),
                        legend.title = element_text(color="black", # 修改图例的标题
                                                            #family = "Ariar",
                                                            size=20,
                                                            face="bold"),
                        axis.title.x  = element_text(size = 20,
                                                     #family = "Ariar",
                                                     color = "black",
                                                     face = "bold",
                                                     vjust = 1.9,
                                                     hjust = 0.5,
                                                     angle = 0),
                        axis.title.y  = element_text(size = 20,
                                                     #family = "Ariar",
                                                     color = "black",
                                                     face = "bold",
                                                     vjust = 1.9,
                                                     hjust = 0.5,
                                                     angle = 90),
                                                     
                        legend.text = element_text(color="black", # 设置图例标签文字
                                                           #family = "Ariar",
                                                           size = 20,
                                                           face = "bold"),
                        axis.text.x = element_text(size = 20,  # 修改X轴上字体大小，
                                                           #family = "Ariar", # 类型
                                                           color = "black", # 颜色
                                                           face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                                           vjust = 0.5, # 位置
                                                           hjust = 0.5,
                                                           angle = 0), #角度
                        axis.text.y = element_text(size =20,  # 修改y轴上字体大小，
                                                           #family = "Ariar", # 类型
                                                           color = "black", # 颜色
                                                           face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                                           vjust = 0.5, # 位置
                                                           hjust = 0.5,
                                                           angle = 0), 
                        legend.position="right")
 p
 
 # #保存
 ggsave(filename = "2-A_VS_2-S火山图.png",p,dpi = 300,width = 8,height =6 )
 pdf(file = "2-A_VS_2-S火山图.pdf",width = 8,height = 6)
 print(p)
 dev.off()
 
 
 data$label = ifelse(data$Genes %in% c("ASPN","MMP7","CP","CXCR4","C8B","SAA4","SMC3","C3","C4A","IGLV1-47","C1S"
                                             ,"CALR","IGHV1-18","C2","HLA-DPB1","HLA-DQB1","TYMP","CLIC5","ST3GAL4","LRRC15",
                                             "GATD1","CPNE5","KMO","COLEC11"),data$Genes,"")
 
 
 table(data$label)
 o= p+geom_text_repel(data = data, aes(x = log2(LN_VS_SLE_FC), 
                                       y = -log10(LN_VS_SLE_Pvalue), 
                                       label = label),
                      size = 4,box.padding = unit(1.2, "lines"),
                      max.overlaps = 5000,
                      point.padding = unit(0.5, "lines"), 
                      segment.color = "black", 
                      
                      
                      show.legend = FALSE)
 o
 
 ggsave(filename = "LN_VS_SLE火山图.png",o,dpi = 300,width = 12,height =9 )
 pdf(file = "LN_VS_SLE火山图.pdf",width = 12,height =9)
 print(o)
 dev.off()
 
 
 
 
 
 
 
 
 ###后面是标注内容，此次不需要
 #添加标签代码
 data$label = rep("",nrow(data))
 #data$label = ifelse(data$miRNA %in% c("mmu-miR-132-3p","mmu-miR-466l-3p"),data$miRNA,"")
 

 
top5.1 = data %>%
  filter(病例_VS_对照_FC<0.8333) %>%
        top_n(n=15,wt = -病例_VS_对照_pvalue)  %>%
        select(Genes)
top5.2 = data %>%
  filter(病例_VS_对照_FC>1.2)  %>%
  top_n(n=15,wt = -病例_VS_对照_pvalue)  %>%
  select(Genes)

top5=rbind(top5.1, top5.2)

data$label = ifelse(data$Genes %in% top5[, 1],data$Genes,"")
 #data$label = ifelse(data$`Gene Symbol` != "KNG1",as.character(data$label),"KNG1")
table(data$label)

o= p+geom_text_repel(data = data, aes(x = log2(病例_VS_对照_FC), 
                                      y = -log10(病例_VS_对照_pvalue), 
                                    label = label),
                   size = 4,box.padding = unit(1.2, "lines"),
                   max.overlaps = 5000,
                   point.padding = unit(0.5, "lines"), 
                   segment.color = "black", 
                   show.legend = FALSE)
o
ggsave(filename = "鞠文浩老师标注火山图.png",o,dpi = 300,width =14,height =10 )
 pdf(file = "鞠文浩老师标注火山图.pdf",width = 14,height = 10)
 print(o)
 dev.off()

 
 
 