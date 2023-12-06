#' @title luminex caculate
#' @description Assignment, statistical analysis and graphing of luminex data
#' @details Input the original standard concentration, sample measurement results, grouping, comparison between groups, dilution ratio and other information for unified data processing and integration analysis
#' @param data_raw,data_nd,paired,compare_data,data_dilution,group,step_increase,start,y_max,p_width,p_height,excel_name,plot_pca,pca_width,pca_height,rm_oor_cytok,hp_width,hp_height
#' @return Violin plot of cytokines and corresponding statistical analysis of components, etc.
#' @export
#' @import ggplot2 ggpubr ggsignif RColorBrewer dplyr ggord ggrepel pheatmap openxlsx
#' @examples

Clean_Plot = function(data_raw,data_nd,compare_data,data_dilution,group,step_increase=0.15,start=1.1,y_max=1.8,paired=F,method="t.test",p_width=4,p_height=5,plot_pca=F,excel_name="luminex_result",pca_width = 8,pca_height=7,rm_oor_cytok=T,plot_heatmap=T,hp_width = 7,hp_height=7) {
  data_raw[data_raw=="N/A"] = NA
  data_raw = data.frame(data_raw[,-c(1,2)])
  data_raw = apply(data_raw, 2, as.numeric)
  data_raw = data.frame(data_raw[,])
  # str(data_raw)
  x = apply(data_nd,2,function(x) gsub("[*]","",x))
  data_ndraw = data_nd[,-1:-2]
  #step1 OOR
  data_nd[data_nd=="OOR <"] = ""
  data_nd[data_nd=="OOR >"] = ""

  #step 2
  data_nd1 = data.frame(x[,-c(1,2)])
  data_nd1 =apply(data_nd1, 2, as.numeric)
  data_nd1 = data.frame(data_nd1)
  row.names(data_nd1) = row.names(data_nd)
  #str(data_nd1)
  # class(data_nd1)

  #step3
  data_nd2 = data_nd1/data_dilution[,-1:-2]
  data_hb = rbind(data_raw,data_nd2)

  for (i in 1:ncol(data_hb)) {
    #replace_data = replace_na(data_nd1[,i],min(data_nd1[,i]/2,na.rm = T))
    data_ndraw[which(data_ndraw[,i]=="OOR <"),i] = (min(data_hb[,i],na.rm = T)/2)*data_dilution[,-1:-2][which(data_ndraw[,i]=="OOR <"),i]
    data_ndraw[which(data_ndraw[,i]=="OOR >"),i] = max(data_hb[,i],na.rm = T)*2*data_dilution[,-1:-2][which(data_ndraw[,i]=="OOR >"),i]
  }

  # for (i in 1:ncol(data_hb)) {
  #   #replace_data = replace_na(data_nd1[,i],min(data_nd1[,i]/2,na.rm = T))
  #   if (sum(is.na(data_hb[,i])) == 0){data_hb[,i] = data_hb[,i]}
  #   else {data_hb[,i] = replace_na(data_hb[,i],min(data_hb[,i]/2,na.rm = T))}
  # }
  data_result = apply(data_ndraw,2,function(x) gsub("[*]","",x))
  data_result = as.data.frame(apply(data_result, 2, as.numeric))
  #str(data_result)
  if (rm_oor_cytok) {
    data_result = data_result[,!(apply(data_result, 2, var)==0)]
  }
  # data_result = data_result[-c(1:nrow(data_raw)),]

  #
  # data = cbind(group,data_result)
  row.names(data_result) = row.names(data_nd)
  data = merge(group,data_result,by =  "row.names",all = T,sort = F)
  if (length(unique(group[,1]))!=2){
    color = brewer.pal(n = length(unique(group[,1])) , name = "Set3")#S
  }else{
    color = c("#8DD3C7","#FFFFB3")
  }
  compare = c(compare_data) ##
  #compare
  #p_data
  p_data = as.data.frame(matrix(nrow=ncol(data)-2,ncol=ncol(compare_data)*2))

  for (y in seq(1,2*ncol(compare_data),by=2)) {
    colnames(p_data)[y] = paste0(compare[[(y+1)/2]][1],"_vs_",compare[[(y+1)/2]][2],"_Pvalue")
    colnames(p_data)[y+1] = paste0(compare[[(y+1)/2]][1],"_vs_",compare[[(y+1)/2]][2],"_FC")
  }
  row.names(p_data) = colnames(data)[-1:-2]

  folder1=paste0("./",excel_name,"_",Sys.Date())
  dir.create(folder1)
  setwd(folder1)
  folder2=paste0("./","Violin_plot/")
  dir.create(folder2)

  for (i in c(3:ncol(data))){
    data_SOD <- dplyr::select(data,2,colnames(data)[i])
    #data_SOD1 <- data_SOD
    colnames(data_SOD) = c("group","cytokine")
    SOD_mean <- data_SOD %>%
      dplyr::group_by(group) %>%
      dplyr::summarize(
        mean = mean(cytokine),
        sd = sd(cytokine))
    ###################################################################

    colnames(data_SOD) = c("group","cytokine")
    #p = pairwise.t.test(data_SOD$cytokine,data_SOD$group,p.adjust.method = "none")
    # p$p.value

    for (j in c(1:ncol(compare_data))) {
      #
      sd = max(var(data_SOD[data_SOD$group==compare[[j]][1],2]),
               var(data_SOD[data_SOD$group==compare[[j]][2],2]))
      p_data[i-2,j*2-1] = ifelse(sd==0,1,
                                 compare_means(cytokine~group,paired = paired,
                                               data = data_SOD[data_SOD$group==compare[[j]][1]|data_SOD$group==compare[[j]][2],],
                                               method = method)[,"p"][[1]])#
      datax = data_SOD[data_SOD$group==compare[[j]][1]|data_SOD$group==compare[[j]][2],]
      p_data[i-2,j*2] = mean(datax[datax$group==compare[[j]][1],2])/mean(datax[datax$group==compare[[j]][2],2])
    }

    p4 = ggviolin(data_SOD, x = "group", y = "cytokine",fill = "group",add = "boxplot",
                  add.params = list(fill = "group")) +#,add = "boxplot",add.params = list(fill = "white")
      stat_compare_means(label = "p.signif",step.increase = step_increase,comparisons = compare,method = method,vjust = 0.25,paired = paired,
                         #args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.001), symbols = c("****", "***", "**", "*", "ns")),
                         tip.length = 0,size=5,label.y = start*max(data_SOD$cytokine))+
      scale_y_continuous(limits =c(0,y_max*max(data_SOD$cytokine)),expand = c(0,0))+
      theme_bw()+
      labs(title= paste0(colnames(data)[i],"-","boxplot"),x="Group",y="Concentration (pg/ml)")+#
      theme(plot.title = element_text(size = 20,
                                      family = "Times",
                                      colour = "red",
                                      hjust = 0.5),
            panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.ticks.length.x = unit(-0.2,"cm"),
            axis.ticks.length.y = unit(-0.2,"cm"),
            axis.title.y = element_text(size = 15,
                                        family = "Times",
                                        color = "black",
                                        face = "bold",
                                        vjust = 1.9,
                                        hjust = 0.5,
                                        angle = 90),
            axis.title.x  = element_text(size = 15,
                                         family = "Times",
                                         color = "black",
                                         face = "bold",
                                         vjust = 1.9,
                                         hjust = 0.5,
                                         angle = 0),
            legend.title = element_text(color="black", #
                                        family = "Times",
                                        size=15,
                                        face="bold"),
            legend.text = element_text(color="black", #
                                       family = "Times",
                                       size = 15,
                                       face = "bold"),
            axis.text.x = element_text(size = 12,  #
                                       family = "Times", #
                                       color = "black", #
                                       face = "bold", #
                                       vjust = 0.5, #
                                       hjust = 0.5,
                                       angle = 45), #
            axis.text.y = element_text(size = 13,  #
                                       family = "Times", #
                                       color = "black", #
                                       face = "bold", #
                                       vjust = 0.5, #
                                       hjust = 0.5,
                                       angle = 0)) #

    ggsave(filename =paste0(folder2,colnames(data)[i], "_violin_plot",".pdf"),plot = p4,width = p_width,height = p_height)
    ggsave(filename =paste0(folder2,colnames(data)[i], "_violin_plot",".png"),plot = p4,width =p_width,height = p_height)
  }
  row.names(data) = data[,1]
  data_out = cbind(p_data,t(data[,-1:-2]))
  write.xlsx(data_out,file = paste0(excel_name,".xlsx"),sheetName = "compare",rowNames=T)
  if (plot_pca) {
    print("pca")
    row.names(data) = data$Row.names
    pca<-prcomp(data[,-1:-2],center=T,scale=T,retx = T)#
    scores=as.data.frame(pca$x)
    newdata2<-scores[,c(1:2)]
    # total1<-merge(newdata2,group,by = "row.names",)
    p <- ggord(pca,
               grp_in=factor(data[,2],levels = unique(data[,2])),
               txt =NULL ,
               arrow=0,alpha_el = 0.4,size = 3,veclsz =1,vec_ext=0)

    p1 <- p+geom_label_repel(aes(newdata2$PC1,newdata2$PC2, fill=data$group,  segment.size = .5,
                                 label=row.names(newdata2)), fontface="bold", color="white", max.iter = 3e3,
                             box.padding=unit(0.35, "lines"), point.padding=unit(0.25, "lines"),
                             segment.colour = "grey50")+ theme_classic(base_size = 16)
    fold3 = "./PCA/"
    dir.create(fold3)

    ggsave(filename = paste0(fold3,"PCA.pdf"),plot = p1,width =pca_width,height = pca_height)
    ggsave(filename = paste0(fold3,"PCA.png"),plot = p1,width = pca_width,height = pca_height)
  }
  if (plot_heatmap) {
    data_heatmap = as.data.frame(t(data[,-1:-2]))
    colnames(data_heatmap) = data[,1]
    annotation_col = group
    annotation_col$group = factor(annotation_col$group,levels = unique(annotation_col$group))
    p_heatmap=pheatmap(data_heatmap,#data1[newOrder[,"Cluster"]==3,]
               #display_numbers = T,
               cluster_rows = T,
               cluster_cols=F,
               scale="row",
               annotation_col = annotation_col,
               #annotation_row = annotation_row,
               #annotation_row = NA,
               #annotation_colors = annotation_colors,
               color=colorRampPalette(rev(c("red","white","blue")))(200),
               # color=c(colorRampPalette(colors = c("blue"))(length(bk)*1.9/5),
               #         colorRampPalette(colors = c("blue","white"))(length(bk)*0.38/3),
               #         colorRampPalette(colors = c("white","red"))(length(bk)/4),
               #         colorRampPalette(colors = c("red"))(length(bk)/4)),
               show_rownames=T,show_colnames=F,
               cellwidth =300/ncol(data_heatmap),cellheight=350/nrow(data_heatmap),
               border_color ="grey",family = "ARL",
               fontsize =14,fontsize_row = min(350/nrow(data_heatmap),20),fontsize_col = 5,fontsize_number=4,
               clustering_method='ward.D2',
               legend = TRUE,drop_levels=T,angle_col = 45,#na_col = "grey",
               annotation_legend = T,
               main="")
    fold4 = "./Heatmap/"
    dir.create(fold4)
    ggsave(filename = paste0(fold4,"Heatmap.pdf"),plot = p_heatmap,width =hp_width,height = hp_height)
    ggsave(filename = paste0(fold4,"Heatmap.png"),plot = p_heatmap,width = hp_width,height = hp_height)
  }
  setwd("../")
}


