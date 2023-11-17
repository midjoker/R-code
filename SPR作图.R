library(tidyverse)
library(hrbrthemes)
library(kableExtra)
library(babynames)
library(streamgraph)
library(viridis)
library(DT)
library(plotly)
setwd("D:/Desktop/胡雪晴/")
# Load dataset from github
data <- readxl::read_xlsx("RAW_mGSDMD-Rb1_Kinetic_20231017.xlsx",2)

data1=gather(data,name,fit_Y,-fit_X)
 #color=RColorBrewer::brewer.pal(n = 10,name = "Set3")
data1$name=factor(data1$name,levels = unique(data1$name))
data1 %>%
  ggplot( aes(x=fit_X, y=fit_Y, group=name, color=name)) +
  geom_line(linewidth = 1.5) +
  scale_color_viridis(discrete = TRUE) +
  scale_color_manual(values =c('#0000ff','#ff00ff','#00ced1','#daa520','#8b008b','#32cd32','#FA7F6F'))+
  xlab("Time(s)")+ylab('Response(RU)')+
  labs(
    color=NULL)+
  theme_classic()+
  theme(
    plot.title = element_text(size=14)
  ) +
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    legend.title=element_text(size=14,face="plain",color="black"),
   # plot.caption.position = "plot",
    plot.caption = element_text(hjust = 1,vjust = 0.5,size=11,color="black"),
    plot.tag.position = c(0,1),
    plot.tag = element_text(vjust = 1.5, hjust = -0.2),
    legend.background  =element_blank(),
    legend.position = c(0.84,0.7)
  )


ggsave(file = "RAW_mGSDMD-Rb1_Kinetic_20231017.png", 
       width = 5, height = 4, dpi = 600, 
       device = "png")
ggsave(file = "RAW_mGSDMD-Rb1_Kinetic_20231017.pdf", 
       width = 5, height = 4, dpi = 600,
       device = "pdf")

