library(ggplot2)
library(openxlsx)
library(ggrepel)
setwd("D:/desktop/袁国红老师-V5/")
ppp=c("0209d","1207d","0221d","1207_D1","0221H","0209dEcm")
term=c("neuron development","angiogenesis","wound healing","axon guidance","axon extension","collagen fibril organization")
n=1
for (sample in ppp) {

data=read.xlsx("袁国红老师重新统计结果.xlsx",n)

for (yyy in term) {
  
setwd("D:/desktop/袁国红老师-V5/")
data=read.xlsx("袁国红老师重新统计结果.xlsx",n)
data1=read.xlsx("整合-细胞外基质全部蛋白.xlsx",1)
data10=read.xlsx("蛋白统计修改.xlsx",1)

data2=read.xlsx(paste0(sample,"_GO.xlsx"),1)  

list_GO=data2[data2$Description==yyy,"geneID"]

# Split the string by "/"
split_vector <- unlist(strsplit(list_GO, "/"))

# Print the resulting vector
print(split_vector)

intersection_vector_1 <- intersect(split_vector,data1$Gene)
data[,sample]=log2(data[,sample])
data <- data[order(data[,sample],decreasing = T),]
data <- data[data[,sample]!=-Inf,]
data$x=1:nrow(data)
percentile_10 <- quantile(data$x, 0.1)


data_sort=data[data$Gene.Symbol %in% intersection_vector_1, ]

data_sort2=data[data$Gene.Symbol %in% data1$Gene,]
ttt=nrow(data_sort2)
sorted_data <- data_sort[order(data_sort[,sample],decreasing = TRUE), ]

top10=sorted_data$Gene.Symbol[1:10]
top10.data=sorted_data[1:10,]
# write.xlsx(top10.data[,1:2],"top10.xlsx")   #输出排名前十的蛋白表格


p=ggplot(data = data,mapping = aes(x = 1:nrow(data),y=!!sym(sample)))+
  geom_rect(aes(xmin = -Inf, xmax = 0.1*(1:nrow(data)), ymin = -Inf, ymax = Inf),
            fill = "orange", alpha = 0.5)+

  geom_point(color="blue")+
  
  geom_point(data = data_sort2,aes(x = x,y = !!sym(sample)),color="black",size = 0.2, stroke = 1.5)+
  geom_point(data=data_sort[1:10,],aes(x = x,y = !!sym(sample)),color="red",shape = 0 ,size = 2, stroke = 1.5)+
  scale_y_continuous(labels = scales::scientific_format())+
 
  geom_text(aes(x =0.25*nrow(data) , y = 12, label = "90%\nabundance"),color="orange",size=6, vjust = -0.5,inherit.aes = FALSE)+
  geom_text(aes(x =0.75*nrow(data) , y = 26, label = paste0("All proteins\nextracellular matrix(",ttt,"/",nrow(data),")")),color="black",size=4, vjust = -0.5,inherit.aes = FALSE)+
  theme_classic()+
  theme(title = element_text(size = 17))+
  labs(x="Ranked proteins",
       y="Abundance",
       title="Protein abundance distribution")

p

data$label = rep("",nrow(data))



data$label = ifelse(data$Gene.Symbol %in% top10,data$Gene.Symbol,"")

table(data$label)






# Adjust the label positions to center them


# o= p+geom_text_repel(data = data, aes(x = 1:nrow(data),
#                                       y = !!sym(sample),
#                                       label = label),
# 
# 
#                      nudge_x = label_x - data$x,
#                      nudge_y = label_y - data$`!!sym(sample)`,
# size = 4,box.padding = unit(1.2, "lines"),
# max.overlaps = 5000,
# # min.segment.length = 0.1,check_overlap = F,
# 
# point.padding = unit(0.5, "lines"),
# segment.color = "red",
# show.legend = FALSE)
# o

 ggsave(filename = paste0(sample,"_",yyy,"_蛋白丰度图.png"),p,dpi = 300,width = 6,height =6)
 ggsave(filename = paste0(sample,"_",yyy,"_蛋白丰度图.pdf"),p,dpi = 300,width = 6,height =6)


 # ggsave(filename = paste0(sample,"_",yyy,"_蛋白丰度图.png"),o,dpi = 300,width = 10,height =10)
 # ggsave(filename = paste0(sample,"_",yyy,"_蛋白丰度图.pdf"),o,dpi = 300,width = 10,height =10)

}
n=n+1
}
