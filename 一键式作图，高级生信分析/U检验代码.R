# 假设数据已经准备好并保存在一个名为data的数据框中
# data应该有70行（每行一个指标）和9列（分组1、2、3各3个重复）
rm(list = ls())

# 加载必要的库
library(stats)

setwd("D:/Desktop/公伟老师/")
data = read.xlsx("公伟老师-总信号校正结果.xlsx",1,rowNames = T)
data=data[,9:16]

str(data)
wilcox_result=c(rep(1,40))
options(digits = 4)
for (i in 1:40) {
  

  wilcox_result[[i]]=wilcox.test(as.numeric(as.character(data[i,1:4])), as.numeric(as.character(data[i,5:8])),exact = F,correct = F)$p.value


}

data2=cbind(wilcox_result,data)
colnames(data2)[1]="MC_vs_NC_pvalue"

write.xlsx(x = data2,file = "公伟老师U检验结2果.xlsx",rowNames = T)
