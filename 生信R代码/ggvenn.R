library(ggVennDiagram)
library(openxlsx)
library(ggplot2)
library(VennDiagram)

data_A <- read.table("AAA.txt",header = TRUE,sep = "\t")
name_A <- "AAA"
data_B <- read.table("BBB.txt",header = TRUE,sep = "\t")
name_B <- "BBB"


list_A <- data_A[1]
names(list_A) <- name_A
list_B <- data_B[1]
names(list_B) <- name_B
list <- c(list_A,list_B)
venn <- ggvenn::ggvenn(list,c(names(list[1]),names(list[2])),
                       fill_color = c("red","blue"),stroke_size = 0.4,
                       text_size = 4,set_name_size = 5 )
ggsave(paste0("Venn_", name_A, "_VS_", name_B,"_异同分析.png"), venn, dpi = 600)  

# 查看交集详情,并导出结果
inter <- get.venn.partitions(list)
inters <- inter[1:(which(colnames(inter) == '..set..')-1)]
inter$setNames <- apply(inters, 1, function(categories) {
  include <- paste(names(list)[categories], collapse = "|")
})
partresult <- do.call(cbind, lapply(lapply(inter[['..values..']], 
                                           unlist), `length<-`,max(lengths(inter[['..values..']]))))
colnames(partresult) <- inter[['setNames']]
partresult <- as.data.frame(partresult)
partresult
write.csv(partresult,"venn_list.csv",na = "")
