

setwd("D:/Desktop/徐磊SPR/")
# Load dataset from github
data <- readxl::read_xlsx("RAW_LZTS2-P6(61-75)_Kinetic_20231023.xlsx",3)

  ggplot(data, aes(x = fit_X, y = Y))+
  geom_point(size = 3,color=c('#0000ff','#ff00ff','#00ced1','#daa520','#32cd32','#FA7F6F')) +
  # geom_line() +  
    # scale_color_manual(values =)+
  theme_classic()+
  xlab("Concentration(nM)")+ylab('Response(RU)')+
    geom_line()
  ggsave(file = "RAW_LZTS2-P6(61-75)_Kinetic_SPR折线图.tiff", 
         width = 5, height = 4, dpi = 600, )
  ggsave(file = "RAW_LZTS2-P6(61-75)_Kinetic_SPR折线图.pdf", 
         width = 5, height = 4, dpi = 600,
         device = "pdf")
  





