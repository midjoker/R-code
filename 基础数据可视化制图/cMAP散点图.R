library(readxl)
library(ggplot2)
setwd("D:/Desktop/上海中医药大学-陈颋老师-cMAP报告/过程文件/")
# 从Excel文件读取数据
data <- read_excel("1、camp结果总表.xlsx")
data=data[1:30,]




# 对数据进行排序，选择前20个点


# 创建散点图
scatter_plot <- ggplot(data, aes(x =1:nrow(data) , y =data$`abs(median_tau_score)` ,color=type)) +
  geom_point() +
  labs(title = "Abs_score point plot",
       x = "Number",
       y = "Abs−Score")+
theme_classic()
scatter_plot
# 在散点图上标注排名前20的点
data$label = rep("",nrow(data))

top20 = data %>%
  top_n(n=20,wt = abs(median_tau_score))  %>%
  select(name)

data$label = ifelse(data$name %in% top20$name,data$name,"")

table(data$label)
o= scatter_plot +geom_text_repel(data = data, aes(x =1:nrow(data) , y =abs(median_tau_score), 
                                      label = label),max.overlaps = 5000)
                    
o
ggsave(filename = "绝对值得分Top20散点图.png",o,dpi = 600,width =6,height =5)
ggsave(filename = "绝对值得分Top20散点图.svg",o,dpi = 600,width =6,height =5)
pdf(file = "绝对值得分Top20散点图.pdf",width =6,height =5)
print(o)
dev.off()





