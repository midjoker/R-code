library(openxlsx)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tibble)
library(reshape2)
if (dir.exists('4. 生信分析')==FALSE) {
  dir.create('4. 生信分析')
}

setwd("D:/desktop/磷酸化廖老师")

fc=1.2
fc1=1/fc
phospratio <- readxl::read_xlsx("表4 磷酸化比值.xlsx",4)
DistributionPlot <- phospratio[8]
names(DistributionPlot) <- 'y'
DistributionPlot <- arrange(DistributionPlot,y)
DistributionPlot$x <- 1:nrow(DistributionPlot)
DistributionPlot$color <- ifelse(DistributionPlot$y<fc1,'Down',
                                 ifelse(DistributionPlot$y>fc,'Up','No Diffenence'))
displot <- ggplot(data=DistributionPlot,aes(x=x,y=y,fill=color))  + 
  geom_col(width=1) + 
  scale_fill_manual(values=c('Down'="#3399FF", 
                             'No Diffenence'="grey", 
                             'Up'="red"))  + 
  geom_hline(linetype='dashed',yintercept = c(fc1,fc))  + 
  geom_hline(linetype='solid',yintercept = 0)  + 
  labs(x = "Protein numbers",y="Ratio_ASN_vs_BSA",title="Ratio Distribution Plot",fill="Ratio")  + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5,size=20),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.position = "top")

ggsave(file = paste0('4. 生信分析/', "ASN_vs_BSA蛋白磷酸化水平变化.png"), plot = displot,
       width = 8.4, height = 8.4, dpi = 600,
       device = "png")



pathwayplot1 = readxl::read_xlsx("表1 PATHWAY分析结果.xlsx",2)
pheatdata <- log2(pathwayplot1[,-1])
row.names(pheatdata) <- pathwayplot1$Name
pheatplot <- pheatmap(pheatdata,
                      cluster_rows=F,cluster_cols=F,
                      scale="none",
                      #annotation_col = annotation_col, 
                      #annotation_row = NA,
                      annotation_colors = annotation_colors,
                      color=colorRampPalette(rev(c("red","white","blue")))(20000),
                      show_rownames=TRUE,show_colnames=TRUE,
                      #cellwidth =24,cellheight=12,
                      border_color =TRUE,#angle_col = '45',
                      #fontsize_row =10,fontsize_col = 13,
                      na_col = "#DDDDDD",
                      legend_breaks = NA,legend = F,main=" ",width = 200,height=60)
ggsave(file = paste0('4. 生信分析/', "图2 差异磷酸化蛋白热图.png"), plot = pheatplot,
       dpi = 600,width = nrow(pheatdata)*0.12, height = nrow(pheatdata)*0.2, 
       device = "png")


pathwaycount <- data.frame(colSums(!is.na(pathwayplot1[,-1])))
names(pathwaycount) <- 'Number of diff phospho'
proteinscount <- readxl::read_xlsx("PEX100 磷酸化位点通路展示.xlsx")
proteinscount <- data.frame(colSums(!is.na(proteinscount)))
names(proteinscount) <- 'Number of proteins'

proteinscount <- rownames_to_column(proteinscount,var = "name")
pathwaycount <- rownames_to_column(pathwaycount,var = "name")
termplotdata <- merge(pathwaycount,proteinscount,by="name")

termplotdata <- melt(termplotdata)
#termplotdata$name <- factor(termplotdata$name,levels = termplotdata$name)
pathplot <- ggplot(termplotdata, aes(fill=variable, y=value, x=name))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25)+
  scale_fill_manual(values=c("#00AFBB", "#FC4E07", "#E7B800"))+
  xlab("")+ylab("Count")+
  theme_bw()+
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    legend.title=element_text(size=14,face="plain",color="black"),
    legend.background  =element_blank(),
    legend.position = c(0.75,0.94)
  )+
  coord_flip()
ggsave(file = paste0('4. 生信分析/', "图3 信号通路蛋白统计图.png"), plot = pathplot,
       width = 8.4, height = 7.2, dpi = 600,
       device = "png")


