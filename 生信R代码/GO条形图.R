library(ggplot2)
library(openxlsx)
df= read.xlsx("Thymosin β4_VS_Control_KEGG_Pathway.xlsx",3)

df$GO<-as.factor(df$GO)
GO_term_order=factor(as.integer(rownames(df)),labels=df$Description)

#报告中对应配色
ggplot(df) + 
  geom_bar(aes(x = GO_term_order, y = Count, fill = GO), stat = "identity") + 
  scale_fill_manual(values = c("#473C8B", "#EE4000", "#FFA500")) + 
  facet_grid(. ~ GO, scales = "free_x",  space = "free") + 
  labs(title = "The Most Enriched GO Terms",  x = "GO term", y = "Numbers of Genes") + 
  theme(axis.text.x = element_text(angle = 80, 
                                   hjust = 1), axis.title.x = element_text(hjust = 0.5), 
        plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


ggplot(data=df, aes(x=GO_term_order,y=Count, fill=GO)) + 
  geom_bar(stat="identity", width=0.8) + 
  coord_flip() +  
  xlab("GO term") + ylab("Num of Genes") + 
  theme_bw()

COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

#不等距
ggplot(data=df, aes(x=GO_term_order,y=Count, fill=GO)) +
  geom_bar(stat="identity", width=0.8)  + 
  scale_fill_manual(values = COLS) + theme_bw()  +
  xlab("GO term") + ylab("Num of Genes") + labs(title = "The Most Enriched GO Terms")+ 
  facet_wrap(vars(GO), strip.position = "bottom", scales = "free_x")+ 
  theme(axis.line=element_line(color="black"),
        axis.text.x=element_text(angle = 80,vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        panel.spacing = unit(0, "lines"),
        panel.border=element_blank(),
        strip.background=element_blank(),
        strip.switch.pad.wrap=unit(1,"cm"))

#等距1
ggplot(data=df, aes(x=GO_term_order,y=Count, fill=GO)) +
  geom_bar(stat="identity", width=0.8)  + 
  scale_fill_manual(values = COLS) + theme_bw()  +
  xlab("GO term") + ylab("Num of Genes") + labs(title = "The Most Enriched GO Terms")+ 
  facet_grid(cols=vars(GO), scales = "free_x",
             space = "free",switch = "x")+ 
  theme(axis.line=element_line(color="black"),
        axis.text.x=element_text(angle = 80,vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.placement = "outside",panel.spacing = unit(0, "lines"),
        panel.border=element_blank(),strip.background=element_blank() )

#等距2
ggplot(data=df, aes(x=GO_term_order,y=Count, fill=GO)) +
  geom_bar(stat="identity", width=0.8)  + 
  scale_fill_manual(values = COLS) + theme_bw()  +
  xlab("GO term") + ylab("Num of Genes") + labs(title = "The Most Enriched GO Terms")+ 
  facet_grid(cols=vars(GO), scales = "free_x",space = "free",switch = "x")+ 
  theme(axis.line=element_line(color="black"),
        axis.text.x=element_text(angle = 80,vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.spacing = unit(0, "lines"))
ggsave(file="GO.png",width = 9, height = 6,dpi = 600,device = "png")
ggsave(file="GO.pdf",width = 9, height = 6,dpi = 600,device = "pdf") 
