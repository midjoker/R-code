library(tidyverse)
library(openxlsx)
library(ggplot2)
library(readxl)
library(gdata)
library(XLConnect)
setwd("D:/Desktop/苗妍老师/")
data1 <- read.xlsx("hsa-miR-877-3p-靶基因预测.xlsx",2)
data2 <- read.xlsx("all.counts.R_45min_Hypothalamus_vs_C_Hypothalamus.edgeR_differential_UP.xlsx",1)

sampleid=data1$diana_microt
sampleid2=data1$mirtarbase
  sampleid3=data1$targetscan
  
list = list(diana_microt=sampleid,mirtarbase =sampleid2,targetscan=sampleid3)


venn <- ggvenn::ggvenn(list,fill_color = c("#fbd9e9","#f5a6bb","#d998b6"),
                       stroke_size = 0.4,text_size = 4,set_name_size = 5,show_elements = F)


venn

dir.create("mRNA_Venn图")
setwd("mRNA_Venn图")
ggsave(filename = "Venn图_数据库交集.png",venn,dpi = 600,width =6,height =6)
ggsave(filename = "Venn图_数据库交集.pdf",venn,dpi = 600,width =6,height =6)


