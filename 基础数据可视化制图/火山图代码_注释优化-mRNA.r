
rm (list = ls ())
library(ggplot2)
library(openxlsx)
library(ggrepel)
library(dplyr)
setwd("D:/Desktop/")
#windowsFonts(Ariar = windowsFont("Arial"))
workbook<-"1.���׶���ͳ��.xlsx"

data<-read.xlsx(workbook,1,rowNames = F)
#workbook1<-"�½� Microsoft Excel ������.xlsx"
#data1<-read.xlsx(workbook1, 1, rowNames=F)

#ȫת¼���滻������
#data[is.na(data$log2FC),1] = -Inf
#data$log2FC = as.numeric(data$log2FC)

#ע�͵�����δ���ʵ���ǿ��Բ����ֶ�����significant�ģ�����ͨ����δ���ֱ�����ɡ�
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
                        legend.title = element_text(color="black", # �޸�ͼ���ı���
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
                                                     
                        legend.text = element_text(color="black", # ����ͼ����ǩ����
                                                           #family = "Ariar",
                                                           size = 20,
                                                           face = "bold"),
                        axis.text.x = element_text(size = 20,  # �޸�X���������С��
                                                           #family = "Ariar", # ����
                                                           color = "black", # ��ɫ
                                                           face = "bold", #  faceȡֵ��plain��ͨ��bold�Ӵ֣�italicб�壬bold.italicб��Ӵ�
                                                           vjust = 0.5, # λ��
                                                           hjust = 0.5,
                                                           angle = 0), #�Ƕ�
                        axis.text.y = element_text(size =20,  # �޸�y���������С��
                                                           #family = "Ariar", # ����
                                                           color = "black", # ��ɫ
                                                           face = "bold", #  faceȡֵ��plain��ͨ��bold�Ӵ֣�italicб�壬bold.italicб��Ӵ�
                                                           vjust = 0.5, # λ��
                                                           hjust = 0.5,
                                                           angle = 0), 
                        legend.position="right")
 p
 
 # #����
 ggsave(filename = "2-A_VS_2-S��ɽͼ.png",p,dpi = 300,width = 8,height =6 )
 pdf(file = "2-A_VS_2-S��ɽͼ.pdf",width = 8,height = 6)
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
 
 ggsave(filename = "LN_VS_SLE��ɽͼ.png",o,dpi = 300,width = 12,height =9 )
 pdf(file = "LN_VS_SLE��ɽͼ.pdf",width = 12,height =9)
 print(o)
 dev.off()
 
 
 
 
 
 
 
 
 ###�����Ǳ�ע���ݣ��˴β���Ҫ
 #��ӱ�ǩ����
 data$label = rep("",nrow(data))
 #data$label = ifelse(data$miRNA %in% c("mmu-miR-132-3p","mmu-miR-466l-3p"),data$miRNA,"")
 

 
top5.1 = data %>%
  filter(����_VS_����_FC<0.8333) %>%
        top_n(n=15,wt = -����_VS_����_pvalue)  %>%
        select(Genes)
top5.2 = data %>%
  filter(����_VS_����_FC>1.2)  %>%
  top_n(n=15,wt = -����_VS_����_pvalue)  %>%
  select(Genes)

top5=rbind(top5.1, top5.2)

data$label = ifelse(data$Genes %in% top5[, 1],data$Genes,"")
 #data$label = ifelse(data$`Gene Symbol` != "KNG1",as.character(data$label),"KNG1")
table(data$label)

o= p+geom_text_repel(data = data, aes(x = log2(����_VS_����_FC), 
                                      y = -log10(����_VS_����_pvalue), 
                                    label = label),
                   size = 4,box.padding = unit(1.2, "lines"),
                   max.overlaps = 5000,
                   point.padding = unit(0.5, "lines"), 
                   segment.color = "black", 
                   show.legend = FALSE)
o
ggsave(filename = "���ĺ���ʦ��ע��ɽͼ.png",o,dpi = 300,width =14,height =10 )
 pdf(file = "���ĺ���ʦ��ע��ɽͼ.pdf",width = 14,height = 10)
 print(o)
 dev.off()

 
 
 