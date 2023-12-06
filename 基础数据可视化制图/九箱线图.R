library(ggplot2)
library(openxlsx)
library(ggrepel)



# X = rnorm(800,mean = 0,sd = 20)
# Y = rnorm(800,mean = 0,sd = 30)
setwd("D:/Desktop/ÖÜÑåÄÝ/")
data = read.xlsx("9ÏóÏÞÍ¼_DOX_VS_Control.xlsx",1,rowNames=F)
# X = data$stool_Burn_vs_Control
# Y = data$serum_Burn_vs_Control

#data = data.frame("stool_Burn_vs_Control"=X,"serum_Burn_vs_Control"=Y)

data$color = NA
data$color = ifelse(data[,2] < log(0.5,2)&data[,3] < log(0.5,2),"red",
                    ifelse(data[,2] > log(2,2)&data[,3] > log(2,2),"red",
                           ifelse(abs(data[,2]) >log(2,2) & abs(data[,3]) <log(2,2),"blue",
                                  ifelse(abs(data[,2])<log(2,2)&abs(data[,3])>log(2,2),"green",
                                         ifelse(abs(data[,2])<log(2,2) & abs(data[,3])<log(2,2),"grey","orange")))))
  
data$label = ifelse(data$color=="red",data$SYMBOL,"")

str(data)
p = ggplot(data=data,aes(log2fc_mRNA,log2fc_protein))+
    geom_point(color=data$color,size=1)+
    ylim(c(-15,15))+
    geom_hline(yintercept = c(log(0.5,2),log(2,2)),linetype=6)+
    geom_vline(xintercept = c(log(0.5,2),log(2,2)),linetype=6)+
    theme_bw()+
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15))
  # geom_text_repel(data = data, aes(x = log2fc_mRNA,
  #                                  y =  log2fc_protein),#label
  #                 size = 2.5,box.padding = unit(0.2, "lines"),
  #                 point.padding = unit(0.2, "lines"),
  #                 segment.color = "black",force_pull = 0,
  #                 label = data$label,max.overlaps = 10000,
  #                 show.legend = FALSE)

p
ggsave(p,filename = "9????Í¼.png",width=6,height=6,dpi=600)
  # geom_text_repel(data = data, aes(x = log2fc_Mrna, 
  #                                  y =  log2fc_protein),#label
  #                 size = 4.5,box.padding = unit(0.2, "lines"),
  #                 point.padding = unit(0.2, "lines"), 
  #                 segment.color = "black", 
  #                 show.legend = FALSE)
# plot(X,Y,col=data$color)
# abline(h=10,col="black")
