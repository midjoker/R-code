rm(list = ls())
library(openxlsx)
library(ggplot2)
library(stringr)
pic_width = 9.6
pic_height = 8
workbook <- "Y_VS_S差异miRNA靶基因_GO富集分析.XLSX"
name <- str_sub(workbook, end = -6)
BP <- readxl::read_xlsx(workbook, 1)
CC <- readxl::read_xlsx(workbook, 2)
MF <- readxl::read_xlsx(workbook, 3)

BP <- BP[1:10,2:3]
BP$GO <- "BP"
BP$Term <- str_sub(BP$Term,start = 12)
CC <- CC[1:10,2:3]
CC$GO <- "CC"
CC$Term <- str_sub(CC$Term,start = 12)
MF <- MF[1:10,2:3]
MF$GO <- "MF"
MF$Term <- str_sub(MF$Term,start = 12)

GO <- rbind(BP, CC, MF)
#View(GO)
p <- ggplot(GO) + 
  geom_bar(aes(x = Term, y = Count, fill = GO), stat = "identity") + 
  scale_fill_manual(values = c("#473C8B", "#EE4000", "#FFA500")) + 
  facet_grid(. ~ GO, scales = "free_x",  space = "free") + 
  labs(title = "The Most Enriched GO Terms",  x = "GO term", y = "Numbers of Genes") + 
  theme(axis.text.x = element_text(angle = 80, 
                                   hjust = 1), axis.title.x = element_text(hjust = 0.5), 
        plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
p
ggsave(file = paste0(name, "_GO", ".png"),
       plot = p, width = pic_width, height = pic_height,
       dpi = 600, device = "png")
ggsave(file = paste0(name, "_GO", ".pdf"),
       plot = p, width = pic_width, height = pic_height,
       dpi = 600, device = "pdf")


data <- read.xlsx("表4 KEGG PATHWAY分析.xlsx",1)
ggplot(data=data[1:15,],aes(x=`Fold.Enrichment`,y=reorder(Term, `Fold.Enrichment`)))+
  geom_point(aes(size=Count,color=-1*log10(PValue)))+scale_colour_gradient(low="green",high="red")+
  scale_size_continuous(range=c(4,12))+
  labs(color=expression(-log[10](PValue)), x="Fold Enrichment",y="Pathway name",title="KEGG Pathway ")+
  theme(plot.title = element_text(hjust = 0.8))+
  theme_bw()+  #主题
  theme(axis.text.x=element_text(size=15,color = "black"),
        axis.text.y=element_text(size=12,color = "black"),
        plot.margin=unit(c(2,3,3,4),"cm"))
ggsave(file="KEGGDB_vs_DC.png",width = 12, height = 12,dpi = 600,device = "png")
ggsave(file="KEGGDB_vs_DC.pdf",width = 12, height = 12,dpi = 600,device = "pdf") 
