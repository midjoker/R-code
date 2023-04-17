#luminex包为自建包，索取后本地安装
library(SFND)
library(openxlsx)

#清空环境数据
rm(list = ls())
#####################################################################################
#根据数据文件夹路径设置路径
# setwd("C:/Users/gpli/Desktop/luminex赋值统计作图分析/")
getwd()
setwd("D:/Desktop/结核杆菌芯片")

#导入数据，数据格式参见输入文件，直接在xlsx表格里增加group和组间比较的sheet即可，顺序不能错
data_nd = read.xlsx("hj.xlsx",1,rowNames = T)
group = read.xlsx("hj.xlsx",2,rowNames = T)
compare_data = read.xlsx("hj.xlsx",3,colNames = F,rowNames = T)

SFND_Plot(data_nd = data_nd,
          excel_name = "min-max统计结果",#表格名称可在引号内修改为实际项目特定的名称
          start = 1.1, #最低显著性标识出现在多少倍的最大值位置
          step_increase = 0.1,#显著性标识之间的上下间隔
          y_max = 1.9, #y轴的最大高度，默认为1.8倍的细胞因子最大值 0.9.5新增
          paired = F,#是否进行配对检验
          method  = "t.test", #"wilcox.test"#检验方法
          compare_data = compare_data,
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

