# 需要提前安装好一下的包ggplot2 ggpubr ggsignif RColorBrewer dplyr ggord ggrepel pheatmap openxlsx
# install.packages("openxlsx")
# 
# # Enable the r-universe repo
# options(repos = c(
#   fawda123 = 'https://fawda123.r-universe.dev',
#   CRAN = 'https://cloud.r-project.org'))
# 
# # Install ggord
# install.packages('ggord')

#luminex包为自建包，索取后本地安装
library(luminex)
library(openxlsx)

#清空环境数据
rm(list = ls())
#####################################################################################
#根据数据文件夹路径设置路径
# 赋值+作图统计分析，需要导入分组及对比文件，每个人电脑的路径不同，需要修改为自己路径
setwd("D:/Desktop/王宇晨老师细胞因子分析/")

#导入数据，数据格式参见输入文件，直接在xlsx表格里增加group和组间比较的sheet即可，顺序不能错
data_raw = read.xlsx("表1 原始数据.xlsx",5,rowNames = T)
  
data_nd = read.xlsx("表3 样本浓度.xlsx",1,rowNames = T)
group = read.xlsx("表3 样本浓度.xlsx",2,rowNames = T)
compare_data = read.xlsx("表3 样本浓度.xlsx",3,colNames = F,rowNames = T)
data_dilution =read.xlsx("表3 样本浓度.xlsx",4,rowNames = T)
str(data_nd)
str(data_raw)
str(compare_data)
dim(data_nd)
names = colnames(data_raw)
names_1 = gsub(pattern = "/",replacement = "_",names)
colnames(data_raw) = names_1
colnames(data_nd) = names_1
colnames(data_dilution) = names_1
#执行分析程序，其中各个参数对应不同的分析绘图差异，见#后的注释，excel_name后一般需要修改
Clean_Plot(data_raw = data_raw,
           data_nd = data_nd,
           excel_name = "王宇晨老师_数据及作图",#表格名称可在引号内修改为实际项目特定的名称
           start = 1.3, #最低显著性标识出现在多少倍的最大值位置
           step_increase = 0.13,#显著性标识之间的上下间隔
           y_max = 1.9, #y轴的最大高度，默认为1.8倍的细胞因子最大值 0.9.5新增
           paired = F,#是否进行配对检验
           method  = "t.test", #"wilcox.test"#检验方法
           compare_data = compare_data,
           data_dilution = data_dilution,
           group = group,
           p_width = 6,#小提琴图宽度
           p_height = 7,#小提琴图高度
           plot_pca = T,#是否绘制PCA图
           pca_width = 11,#PCA图宽度，样本较少时设置8左右就行
           pca_height = 10,#PCA图高度
           rm_oor_cytok = T,#是否移除全部为OOR的指标,如果不移除可能导致PCA出问题
           plot_heatmap = T,#是否绘制热图
           hp_width = 10,#热图宽度，根据data_nd的行列设置
           hp_height = 10)#热图高度，根据data_nd的行列设置


#####################################################################################
只赋值
 setwd("D:/Desktop/岳亮/") #设置路径
 data_raw = read.xlsx("表1 原始数据.xlsx",5,rowNames = T)
 data_nd = read.xlsx("表3 样本浓度.xlsx",1,rowNames = T)
 data_dilution =read.xlsx("表3 样本浓度.xlsx",2,rowNames = T)
 Clean(data_raw = data_raw,
       data_nd = data_nd,
       data_dilution = data_dilution,
       excel_name = "赋值后数据" )#表格名称可修改








