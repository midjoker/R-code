# #安装各种包
 if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 BiocManager::install("cgdsr")

library(cgdsr)
library(DT) 
library(survminer)
library(xlsx)
#设置路径及导入数据
setwd("D:/Desktop/孙洋售后/340生存曲线")
#前期准备
mycgds <- CGDS("http://www.cbioportal.org/")
#测试是否成功
test(mycgds)
#查看所有数据集，挑选某一个数据集，复制其名称如肝癌数据集：'lihc_tcga'
all_tcga_studies =getCancerStudies(mycgds)
DT::datatable(all_tcga_studies)
#将对应数据集导入
mycancerstudy = 'paad_tcga'  
mycaselist = getCaseLists(mycgds,mycancerstudy)[6,1]#第1列第6行代表mRNA数据
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[7,1]#第1列第2行代表mRNA表达量
#指定基因，当有大量基因时候可导入表格，选择基因名那一列
#choose_genes = c("HMGB1","SERPINH1","RAB14","GPC1")
data_genelist = read.xlsx("待分析蛋白(1).xlsx",8,check.names=F,encoding = "UTF-8")
choose_genes = data_genelist$`Gene Symbol`
# choose_genes = "TMC4"
#下载表达量数据
expr = getProfileData(mycgds,choose_genes,mygeneticprofile,mycaselist)
write.csv(expr,"340生存曲线.csv")
sum(is.na(expr))
## 下载临床数据
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]
myclinicaldata = getClinicalData(mycgds,mycaselist)#[11,1]
#选取生存数据，OS代表总生存期，另外有"DFS_STATUS","DFS_MONTHS"代表无进展生存期
choose_columns = c("OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS")
choose_clinicaldata = myclinicaldata[row.names(expr),choose_columns]
#绘制简单生存曲线
library(survival)
data = cbind(expr,choose_clinicaldata)
table(data$OS_STATUS)
table(data$DFS_STATUS)


for (i in 1:length(choose_genes)){
  #设置中位表达量值
   # dat = data[ifelse(is.na(data$HMGB1),F,T),]
   #dat = data[!is.na(data$HMGB1),]
  dat = data
  dat$expression = ifelse(dat[,choose_genes[i]]>quantile(dat[,choose_genes[i]],probs = 0.75),
                          "high","low")
  #DFS的话，dat$DFS_STATUS=="1:DECEASED",DFS为dat$DFS_STATUS=="1:Recurred/Progressed"
  kmfit = survfit(Surv(dat$OS_MONTHS,dat$OS_STATUS=="1:DECEASED")~expression,data = dat)
 
  
  p = ggsurvplot(kmfit,data = dat,pval = T,
             risk.table = T, # 添加风险表
             xlab = "Follow up time(m)", # 指定x轴标签
             legend = c(0.77,0.85), # 指定图例位置
             font.legend = 17,#legend字体大小
             legend.title = "", #设置图例标题，这里设置不显示标题，用空格替代
             legend.labs = c("high", "low"), #指定图例分组标签
             #break.x.by = 200,  # 设置x轴刻度间距
             ggtheme = theme_bw(), #主题设置
             surv.median.line = "hv", #标注中位生存时间
             title=paste0(choose_genes[i],"-","OS survival")) #标题)
  pdf(file = paste0(choose_genes[i], "_survival_OS", ".pdf"))
  print(p,newpage = FALSE)
  dev.off()
  
  
  kmfit = survfit(Surv(dat$DFS_MONTHS,dat$DFS_STATUS=="1:Recurred/Progressed")~expression,data = dat)
  
  
  p = ggsurvplot(kmfit,data = dat,pval = T,
                 risk.table = T, # 添加风险表
                 xlab = "Follow up time(m)", # 指定x轴标签
                 legend = c(0.77,0.85), # 指定图例位置
                 font.legend = 17,#legend字体大小
                 legend.title = "", #设置图例标题，这里设置不显示标题，用空格替代
                 legend.labs = c("high", "low"), #指定图例分组标签
                 #break.x.by = 200,  # 设置x轴刻度间距
                 ggtheme = theme_bw(), #主题设置
                 surv.median.line = "hv", #标注中位生存时间
                 title=paste0(choose_genes[i],"-","DFS survival")) #标题)
  pdf(file = paste0(choose_genes[i], "_survival_DFS", ".pdf"))
  print(p,newpage = FALSE)
  dev.off()
}










