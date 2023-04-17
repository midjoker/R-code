library(ggplot2)
library(grid)
library(openxlsx)
setwd("D:/Desktop/HWFAA03220210108-中国人民解放军总医院-尹彤老师-PEX100 芯片项目实验报告/3. 实验数据")
bplot = read.xlsx("表4 组间比较结果.xlsx",2,rowNames = T)
View(bplot)
# bplot2=read.xlsx("bplot.xlsx",sheetIndex = 2)
# bplot3=read.xlsx("bplot.xlsx",sheetIndex = 3)
# bplot4=read.xlsx("bplot.xlsx",sheetIndex = 4)

tiff(filename = "2倍差异-Protein Combination.tiff",width = 4508 ,height = 3000,res=700)
bplot$color = ifelse(bplot$敲减_vs_对照_Ratio
>=1.5,"Up",ifelse(bplot$敲减_vs_对照_Ratio <=0.667,"Down","No Difference"))
ggplot(data=bplot,aes(x=reorder(c(584:1),敲减_vs_对照_Ratio),y=敲减_vs_对照_Ratio,fill=color))  + 
  geom_col(width=1) + 
  scale_fill_manual(values=c('Down'="#3399FF", 'No Difference'="grey", 'Up'="red"))  + 
  geom_hline(linetype='dashed',yintercept = c(0.667,1.5))  + 
  geom_hline(linetype='solid',yintercept = 0)  + 
  labs(x = "Posphoprotein numbers",y="Ratio",title="",fill="Ratio")  + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5,size=25),axis.title=element_text(size=12),
        legend.text=element_text(size=15),legend.title=element_text(size=15),legend.position = "top")
dev.off()
