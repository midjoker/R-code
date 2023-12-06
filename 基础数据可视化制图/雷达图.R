# 安装并加载所需的R包
devtools::install_github("ricardo-bion/ggradar", dependencies=TRUE)
library(ggradar)
pacman::p_load(ggradar)
# 构建示例数据
library(dplyr)
library(scales)
library(tibble)



data_radar <- read.xlsx("表4 组间比较结果.xlsx",2) %>% 
  as_tibble(rownames = "group") %>% 
  mutate_at(vars(-group), rescale) %>% 
  tail(4) %>% 
  select(1:10)
# 查看示例数据
mtcars_radar<- read.xlsx("表4 组间比较结果.xlsx",2)
normalized_mtcars_radar = (mtcars_radar - min_value) / (max_value - min_value)

## # A tibble: 4 x 10
##   group            mpg   cyl  disp    hp  drat    wt   qsec    vs    am
##   <chr>          <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl>
## 1 Ford Pantera L 0.230   1   0.698 0.749 0.673 0.424 0          0     1
## 2 Ferrari Dino   0.396   0.5 0.184 0.435 0.396 0.321 0.119      0     1
## 3 Maserati Bora  0.196   1   0.573 1     0.359 0.526 0.0119     0     1
## 4 Volvo 142E     0.468   0   0.124 0.201 0.622 0.324 0.488      1     1

# 使用ggradar函数绘制雷达图
p=ggradar(mtcars_radar,base.size = 12,
        values.radar = c("0%","25%","50%","75%","100%"),
        legend.title = "group",legend.text.size = 12,
        legend.position = "right")
p
ggsave(filename = "何芹老师雷达图-差异部分.png",p,dpi = 600,width = 15,height =15 )
ggsave(filename = "何芹老师雷达图-差异部分.pdf",p,dpi = 600,width = 30,height =30 )


