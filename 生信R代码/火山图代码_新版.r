??ɽͼ
 library(ggplot2)
 library(xlsx)
 setwd("C:/Users/wayenbiotech_103/Downloads")
 workbook<-"蛋白定量统计-火山图分析.xlsx"
 data<-read.xlsx(workbook,5,check.names=F,encoding = 'UTF-8')
 View(data)
 #data$label = ifelse(data$Pvalue<0.05&abs(log2(data$fc))>=1,data$`Majority protein IDs`,"")
 c1=(1/1.5)
 c2=(-log10(0.05))
 threshold<-as.factor((data$fc>=1.5|data$fc<= c1)& (data$Pvalue<0.05))
 #View(threshold)
 r1=ggplot(data, aes(x=log2(fc), y=-log10(Pvalue)))+geom_point(alpha=1, size=1.5)+ xlim(c(-7.6, 7.6)) + ylim(c(0,3.5))+xlab("log2(Fold Change)") + ylab("-log10(p-value)")+labs(title=" ",hjust=0.5)+geom_hline(yintercept=1.3,linetype=3)+geom_vline(xintercept=c(-log2(c1),log2(c1)),linetype=3)
 r2=r1+geom_point(aes(color =significant))
 
 volcano=r2+scale_color_manual(values =  c("blue","black","red"))
 volcano
 dev.off()


#??ɽͼ+????

library(ggplot2)
library(xlsx)
 library(ggrepel)
setwd("C:/Users/wayenbiotech_103/Downloads")
workbook<-"volcano_1.xlsx"
data<-read.xlsx(workbook,1,check.names=F,encoding = 'UTF-8')
#data$label = ifelse(data$Pvalue<0.05&abs(log2(data$fc))>=1,data$ID,"")
View(data)
c1=(1/1.5)
c2=(-log10(0.05))
threshold<-as.factor((data$fc>=1.2|data$fc<= c1)& (data$Pvalue<0.05))
r1=ggplot(data, aes(x=log2(fc), y=-log10(Pvalue)))+geom_point(alpha=1, size=1.5)+ xlim(c(-7.6, 7.6)) + ylim(c(0,3.5))+xlab("log2(Fold change)") + ylab("-log10(p_value)")+labs(title="volcano",hjust=0.5)+geom_hline(yintercept=1.3,linetype=3)+geom_vline(xintercept=c(-log2(c1),log2(c1)),linetype=3)+geom_text(label=paste(data$label),colour="black",size=3)
r2=r1+geom_point(aes(color =significant))
volcano=r2+scale_color_manual(values =  c("blue","black","red"))
volcano
dev.off()

#?ÿ???
library(ggplot2)
library(xlsx)
library(ggrepel)
library(dplyr)
setwd("C:/Users/wayenbiotech_103/Downloads")
workbook<-"volcano_1.xlsx"
data<-read.xlsx(workbook,3,check.names=F,encoding = 'UTF-8')

View(data)
#ע?͵??????δ???ʵ???ǿ??Բ????ֶ?????significant?ģ?????ͨ?????δ???ֱ?????ɡ?
#data$change = ifelse(data$Pvalue < 0.05 & abs(log2(data$fc)) >= 0.584962500721156, 
#                     ifelse(log2(data$fc)> 1 ,'Up','Down'),
#                     'No')
p <- ggplot(data = data, 
            aes(x = log2(data$fc), 
                y = -log10(data$Pvalue), 
                colour=significant,
                label = data$ID)) +
        geom_point(alpha=0.4, size=3.5) +
        scale_color_manual(values=c("blue", "black","red"))+
        xlim(c(-7.6, 7.6)) +geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
        geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
        labs(x="log2(fold change)",
             y="-log10 (p-value)",
             title="")  + 
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5), 
              legend.position="right", 
              legend.title = element_blank())
p
#???ӱ?ǩ????
data$label=ifelse(data$Pvalue < 0.05 & abs(log2(data$fc)) >= 0.584962500721156,as.character(data$ID),"")

p+geom_text_repel(data = data, aes(x = log2(data$fc), 
                                   y = -log10(data$Pvalue), 
                                   label = label),
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)
dev.off()

