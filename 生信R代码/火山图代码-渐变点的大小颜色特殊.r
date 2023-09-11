#好看的
rm (list = ls ())
library(ggplot2)
library(openxlsx)
library(ggrepel)
library(dplyr)
setwd("D:/Desktop/魏斌老师售后处理")
#windowsFonts(Ariar = windowsFont("Arial"))
workbook<-"1. 蛋白定量统计.xlsx"
data<-read.xlsx(workbook,1,rowNames = F)
#全转录组替换负无穷
data[is.na(data$`log2(t3_VS_t0_Ratio)`),2] = -Inf
data$`log2(t3_VS_t0_Ratio)` = as.numeric(data$`log2(t3_VS_t0_Ratio)`)
#注释掉的这段代码实际是可以不用手动输入significant的，可以通过这段代码直接生成。
# data$signifcant = ifelse(data$B_VS_C_Pvalue < 0.05 & abs(log2(data$B_VS_C_Ratio)) >= log2(1.2), 
#                       ifelse(data$B_VS_C_Ratio >= 1.2 ,'Up','Down'),'No')
  data$signifcant = ifelse(data$t3_VS_t0_Pvalue < 0.05 & abs(data$`log2(t3_VS_t0_Ratio)`) >= log2(1.2),
                      ifelse(data$`log2(t3_VS_t0_Ratio)` >= log2(1.2) ,'Up(23)','Down(10)'),'No')

# data$signifcant = ifelse(data$pvalue < 0.05 & abs(data$LOG2FC) >= 1,
#                          ifelse(data$LOG2FC >= 1 ,'Up','Down'),'No')


table(data$signifcant)
str(data)
 p <- ggplot(data = data, 
              aes(x = log2(t3_VS_t0_Ratio), 
                 y = -log10(t3_VS_t0_Pvalue), 
                 colour=-log10(t3_VS_t0_Pvalue))) +#                 label = Name
         geom_point(alpha=0.6, aes(size=-log10(data$t3_VS_t0_Pvalue))) +
   scale_size(range = c(2,5))+
         # scale_color_manual(values=c("blue", "grey","red"))+
   scale_color_gradient(low = "blue",high = "red")+
         xlim(c(-8, 11)) +ylim(c(0,5))+geom_vline(xintercept=c(log(0.8333,2),log(1.2,2)),lty=4,col="grey",lwd=0.8) +
         geom_hline(yintercept = -log10(0.05),lty=5,col="grey",lwd=0.8) +
         labs(x="log2(fold change)",
              y="p-value",
              title="",
              color="-log10(pvalue)")  + guides(size = "none")+
         theme_bw()+theme(panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),
                        legend.title = element_text(color="black", # 修改图例的标题
                                                            #family = "Ariar",
                                                            size=16,
                                                            face="bold"),
                        axis.title.x  = element_text(size = 12,
                                                     #family = "Ariar",
                                                     color = "black",
                                                     face = "bold",
                                                     vjust = 1.9,
                                                     hjust = 0.5,
                                                     angle = 0),
                        axis.title.y  = element_text(size = 12,
                                                     #family = "Ariar",
                                                     color = "black",
                                                     face = "bold",
                                                     vjust = 1.9,
                                                     hjust = 0.5,
                                                     angle = 90),
                                                     
                        legend.text = element_text(color="black", # 设置图例标签文字
                                                           #family = "Ariar",
                                                           size = 15,
                                                           face = "bold"),
                        axis.text.x = element_text(size = 12,  # 修改X轴上字体大小，
                                                           #family = "Ariar", # 类型
                                                           color = "black", # 颜色
                                                           face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                                           vjust = 0.5, # 位置
                                                           hjust = 0.5,
                                                           angle = 0), #角度
                        axis.text.y = element_text(size = 13,  # 修改y轴上字体大小，
                                                           #family = "Ariar", # 类型
                                                           color = "black", # 颜色
                                                           face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                                           vjust = 0.5, # 位置
                                                           hjust = 0.5,
                                                           angle = 0), 
                        legend.position="right")
 p
 ggsave(filename = "Model_Control.png",p,dpi = 600,width = 8,height =10 )
 pdf(file = "Model_Control.pdf",width = 8,height = 10)
 print(p)
 dev.off()
 
  
 #添加标签代码
 data$label = rep("",nrow(data))
 #data$label = ifelse(data$miRNA %in% c("mmu-miR-132-3p","mmu-miR-466l-3p"),data$miRNA,"")
 top5.1 = data %>%
   filter(log2(t3_VS_t0_Ratio)< -1) %>%
   top_n(n=10,wt = -t3_VS_t0_Pvalue)  %>%
   select(Gene.Symbol)
 top5.2 = data %>%
   filter(log2(t3_VS_t0_Ratio)> 1)  %>%
   top_n(n=10,wt = -t3_VS_t0_Pvalue)  %>%
   select(Gene.Symbol)
 

 # top5=rbind(top5.1, top5.2) 
 top5=c("GSTP1","HBB",
            "HBA",
            "LOC100624077",
            "COL1A2",
            "COL1A1") 
 
 data$label = ifelse(data$Gene.Symbol %in% top5,data$Gene.Symbol,"")
 #data$label = ifelse(data$`Gene Symbol` != "KNG1",as.character(data$label),"KNG1")
 table(data$label)
 
 o= p+geom_text_repel(data = data, aes(x = log2(t3_VS_t0_Ratio), 
                                       y = -log10(t3_VS_t0_Pvalue), 
                                       label = label),
                      size = 4,box.padding = unit(2, "lines"),
                      max.overlaps = 5000,
                      point.padding = unit(0.5, "lines"), 
                      segment.color = "black", 
                      show.legend = FALSE)
 o
 ggsave(filename = "t3_VS_t0火山图.png",o,dpi = 300,width = 10,height =6 )
 pdf(file = "t3_VS_t0火山图.pdf",width = 10,height = 6)
 print(o)
 dev.off()
 
 
 #data$label = ifelse(data$Pvalue<0.05&abs(log2(data$fc))>=1,data$`Majority protein IDs`,"")
 # c1=(1/2)
 # c2=(-log10(0.01))
 # threshold<-as.factor((data$fc>=2|data$fc<= c1)& (data$Pvalue<0.01))
 # #View(threshold)
 # r1=ggplot(data, aes(x=log2(FC), y=-log10(Pvalue_t.test)))+
 #         geom_point(alpha=1, size=1.5)+ 
 #         xlim(c(-4, 4)) + ylim(c(0,3.5))+
 #         xlab("log2(Fold Change)") + ylab("-log10(p-value)")+
 #         labs(title=" ",hjust=0.5)+geom_hline(yintercept=1.3,linetype=3)+
 #         geom_vline(xintercept=c(-log2(c1),log2(c1)),linetype=3)
 # r1
 # r2=r1+geom_point(aes(color =threshold))
 # 
 # volcano=r2+scale_color_manual(values =  c("blue","black","red"))
 # volcano
 # dev.off()


#火山图+名称

# library(ggplot2)
# library(xlsx)
#  library(ggrepel)
# setwd("c:/users/shenqueying/desktop")
# workbook<-"volcano.xlsx"
# data<-read.xlsx(workbook,1,check.names=F,encoding = 'UTF-8')
# #data$label = ifelse(data$Pvalue<0.05&abs(log2(data$fc))>=1,data$ID,"")
# View(data)
# c1=(1/2)
# c2=(-log10(0.01))
# threshold<-as.factor((data$fc>=2|data$fc<= c1)& (data$Pvalue<0.01))
# r1=ggplot(data, aes(x=log2(fc), y=-log10(Pvalue)))+geom_point(alpha=1, size=1.5)+ xlim(c(-7.6, 7.6)) + ylim(c(0,3.5))+xlab("log2(Fold change)") + ylab("-log10(p_value)")+labs(title="volcano",hjust=0.5)+geom_hline(yintercept=1.3,linetype=3)+geom_vline(xintercept=c(-log2(c1),log2(c1)),linetype=3)+geom_text(label=paste(data$label),colour="black",size=3)
# r2=r1+geom_point(aes(color =significant))
# volcano=r2+scale_color_manual(values =  c("blue","black","red"))
# volcano
# dev.off()

