pacman::p_load("pathview")

# box_ggplot
box_ggplot = function(df, x, y, color=NULL, fill=NULL, x_levels = NULL,
                      x_lab = "",y_lab = "",
                      x.axis.rot = 0,x.axis.size=5,y.axis.size=11,
                      legend_hide = T, boxwidth = .5, lwd = 1, dodge_width=.6,
                      manual_colors = NULL, stat_line=F){
  # color pallete info = print(brewer.pal.info)
  if(is.null(x_levels)){
    if(is.null(levels(df[[x]]))){
      x_levels = unique(df[[x]])
    } else {
      x_levels = levels(df[[x]])
    }
  }
  
  if(is.null(color) & is.null(fill)){
    gpl = ggplot(df, aes(x=factor(!! sym(x),levels=x_levels),y= !! sym(y)))
  } else if(is.null(color)) {
    gpl = ggplot(df, aes(x=factor(!! sym(x),levels=x_levels),y= !! sym(y),
                         fill= !! sym(fill)))
    if(is.null(manual_colors)){
      gpl = gpl + scale_fill_manual(values = brewer.pal(length(x_levels), "Paired"))
    } else { gpl = gpl + scale_fill_manual(values = manual_colors) }
  } else {
    gpl = ggplot(df, aes(x=factor(!! sym(x), levels=x_levels),y= !! sym(y),
                         color= !! sym(color)))
    if(is.null(manual_colors)){
      gpl = gpl + scale_color_manual(values = brewer.pal(length(x_levels), "Paired"))
    } else { gpl = gpl + scale_color_manual(values = manual_colors) }}
  
  gpl=gpl+geom_boxplot(width=boxwidth,lwd=lwd,alpha=1,
                       position=position_dodge(width=dodge_width))+
    theme_bw()+theme(axis.text.x=element_text(angle=x.axis.rot,size =x.axis.size),
                     axis.title=element_text(size=10),
                     axis.title.y=element_text(size = y.axis.size))+labs(x=x_lab,y=y_lab)
  
  if(stat_line){
    gpl=gpl+stat_summary(fun.y=median, geom='line',
                         aes(group = !! sym(ifelse(is.null(color), fill, color)),
                             colour = !! sym(ifelse(is.null(color), fill, color))),
                      position=position_dodge(width=dodge_width))
    if(is.null(manual_colors)){
      gpl = gpl + scale_color_manual(values = brewer.pal(length(x_levels), "Paired"))
    } else { gpl = gpl + scale_color_manual(values = manual_colors) }
    
    }
  if(legend_hide){ gpl = gpl + theme(legend.position="none") }
  return(gpl)
}



# TM_PA_ORA_multiple_condition
TM_PA_ORA_multiple_condition = function(df, order,
                                        code_df, cluster_info = NULL,
                                        lgd_wid = 1,x_lab_size=13, y_lab_size=13,
                                        title = NULL, pval_scipen_digit_0=F){
  if(is.null(cluster_info)){
    dcp_order = code_df$Description[length(code_df$Description):1]
    df = df %>% filter(Description %in% dcp_order)
  } else {
    dcp_order_merge = merge(code_df, cluster_info)
    dcp_order_merge = dcp_order_merge[order(dcp_order_merge$cluster,
                                            dcp_order_merge$code),]
    dcp_order = dcp_order_merge[nrow(dcp_order_merge):1,]$Description
    df = df %>% filter(Description %in% dcp_order)
  }
  print(dcp_order)
  gpl = ggplot(data = df, aes(x = factor(cond, levels = order),
                              y = factor(Description ,level = dcp_order),
                              color = `p.adjust`, size = Count)) + 
    geom_point() + 
    theme_bw() + 
    labs(title= title, x ="", y = "") +
    theme(
      axis.text.x = element_text(
        color = "grey20",size=x_lab_size,angle=90,hjust=.5,vjust=.5,face="plain"),
      axis.text.y = element_text(
        color = "grey20",size=y_lab_size,angle=0,hjust=1,vjust=0,face="plain"),  
      axis.title.x = element_text(
        color = "grey20", size= 12,angle=0, hjust = .5, vjust = 0, face ="plain"),
      axis.title.y = element_text(
        color = "grey20", size= 12,angle=90,hjust =.5, vjust =.5, face="plain"),
      plot.title = element_text(
        hjust = 1, # 타이틀 정중앙위치 = 0.5 / 우측정렬 = 1.0
        vjust = 1,   # 필요하면 세로 위치도 조정 가능
        size = 14,   # 타이틀 크기 (원하는 크기로 변경)
        # face = "bold" # 글꼴 스타일 (bold, italic 등)
      )) + 
    theme(#legend.justification = "right",
      #legend.position = "bottom",
      # legend.position = c(-0.1, -0.1),
      legend.key.size = unit(0.3, 'cm'),
      legend.key.width= unit(lgd_wid,'cm'),
      # legend.margin = margin(0, 0, 0, 0),
      # legend.direction = "horizontal",
      # face = 'bold',
      legend.title = element_text(color = 'black', size = 10),
      legend.text= element_text(size = unit(9, 'cm')))
  if(pval_scipen_digit_0){
    gpl = gpl + scale_color_gradient(low = "red", high = "blue",
                                     labels = label_scientific(digits = 0))
  }else{
    gpl = gpl + scale_color_gradient(low = "red", high = "blue")
  }
  
  if(!is.null(cluster_info)){
    for(my_y in nrow(dcp_order_merge) - (which(!duplicated(dcp_order_merge$cluster,
                                                           fromLast =T))-0.5)){
      gpl = gpl + annotate("segment", x = 0, xend = (length(order)+1),
                           y = my_y, yend = my_y, colour = "black")
    }
  }
  
  return(gpl)
}

# TM_PA_ORA_fromTCGA
TM_PA_ORA_fromTCGA = function(df1, df2 = NULL){
  toPlot = data.frame("v1" = df1$Description,
                      "v2" = df1$p.adjust_log,
                      "v3" = as.numeric(matrix(unlist(strsplit(df1$GeneRatio, "/")),
                                               nrow = 2)[1,])/
                        as.numeric(matrix(unlist(strsplit(df1$BgRatio, "/")),
                                          nrow = 2)[1,]))
  
  
  toplot_fun = function(toPlot){
    color = c("#B2B200", "cyan", "green", "yellow"); xlim = NULL
    xAxis <- EDASeq::barplot(toPlot[, 2], horiz = TRUE,
                             col = color[1], main = "GO:Biological Process",
                             xlab = "-log10(FDR)", xlim = xlim)
    labs <- toPlot[, 1]
    text(x = 1, y = xAxis, labs, pos = 4, cex = 1)
    lines(x = toPlot[, 3], y = xAxis, col = "red")
    points(x = toPlot[, 3], y = xAxis, col = "red")
    axis(side = 3, at = pretty(range(0:1)), col = "red")
  }
  
  
  if(is.null(df2)){
    toplot_fun(toPlot)
  } else {
    par(mfrow = c(1,2))
    toPlot2 = data.frame("v1" = df2$Description,
                         "v2" = df2$p.adjust_log,
                         "v3" = as.numeric(matrix(unlist(strsplit(df2$GeneRatio, "/")),
                                                  nrow = 2)[1,])/
                           as.numeric(matrix(unlist(strsplit(df2$BgRatio, "/")),
                                             nrow = 2)[1,]))
    
    toplot_fun(toPlot)
    toplot_fun(toPlot2)
  }
  return(TRUE)
}


# TM_PA_GO_emapplot
TM_PA_GO_emapplot = function( resGO, EMA_n = 10 ){
  
  resGO_PT = enrichplot::pairwise_termsim(resGO)
  gg_p1 = enrichplot::emapplot(resGO_PT, showCategory = EMA_n,
                               min_edge = 0.2,
                               cex_label_category = 0.7,
                               cex_category = 0.7,
                               cex_line = 0.7,
                               repel = FALSE)
  return(gg_p1)
}

# TM_PA_ORA_barplot
TM_PA_ORA_barplot = function(resGO, title = NULL, BAR_n = 10, str_wrap = 50,
                             y = "Description", x, color, scientific_legend = FALSE){
  gg_p1 = bar_ggplot(resGO, title = title, BAR_n = BAR_n, y, x, color, str_wrap=str_wrap)
  return(gg_p1)
}


# TM_PA_ORA_barplot_UpDw
TM_PA_ORA_barplot_UpDw = function(GOresUp, GOresDw,
                                  title_Up = NULL, title_Dw = NULL, BAR_n = 10,
                                  x, color, scientific_legend = FALSE ){
  gg_p1 = bar_ggplot(GOresUp, title = title_Up, BAR_n=BAR_n, x, color, scientific_legend)
  gg_p2 = bar_ggplot(GOresDw, title = title_Dw, BAR_n=BAR_n, x, color, scientific_legend)
  return((gg_p1 + gg_p2))
}



# TM_PA_GSE_bar
TM_PA_GSE_bar = function(GSEres, title ="", x = "Description",
                         # y = "NES",  y is setted with "adj.P.Val".
                         xlab = "Pathway", ylab = "Normalized Enrichment Score",
                         x.size = 10, y.size = 10){
  GSEres$p.adjust_log = -log10(GSEres$adj.P.Val)
  GSEres$Regulation = factor(GSEres$Regulation, levels = c("Up", "Down"))
  
  GSEres$y_ax = ifelse(GSEres$Regulation == "Up",
                       abs(GSEres$p.adjust_log),
                       -abs(GSEres$p.adjust_log))
  GSEres = GSEres %>% #arrange(., p.adjust_log, !! sym(y)) %>% 
    arrange(., y_ax)
  GSEres$x_ax = factor(GSEres[[x]], levels= GSEres[[x]])
  
  #GSEres$Description = paste(GSEres$ID, GSEres$Description, sep = "~")
  g1 = ggplot(GSEres, aes(x_ax, y_ax)) + 
    geom_col(aes(fill=Regulation)) + # facet_grid(Regulation~., scales="free_y") +
    scale_fill_manual(values=c("Up"="#A43029","Down"="#6E95E6")) +
    coord_flip() + theme(axis.text.x = element_text(size = x.size),
                         axis.text.y = element_text(size = y.size)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
    labs(x=xlab, y=ylab, title=title)
  return(g1)
}


# bar_ggplot
bar_ggplot = function(resGO, title = NULL, BAR_n = 10,
                      y = "Description", x, color, scientific_legend=T,str_wrap=50){
  # plot.margin, legend.position, legend.text 는 manual 수정하는 것을 권장.
  resGO[[y]] = factor(resGO[[y]], levels = rev(unique(resGO[[y]])))
  
  g1 = resGO %>% as.data.frame() %>% .[1:BAR_n,] %>% 
    ggplot(aes(x = !! sym(x),
               # y = fct_reorder(!! sym(y), !! sym(x), .desc = FALSE),
               y = !! sym(y),
               fill = !! sym(color) )) + geom_col() +
    scale_fill_gradientn(colours=c("#FD000F", "#4200F5"), 
                         guide=guide_colorbar(reverse=FALSE),
                         label = function(x) format(x,scientific=scientific_legend,
                                                    digits =3)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = str_wrap), 
                     position = "left") +
    theme_bw() + labs(title = title) +
    theme(axis.text.y = element_text(size = 10),
          
          plot.margin = grid::unit(c(0.1, 0.1, 3, 0.1), 'lines'),
          axis.title.y = element_blank(),
          
          # legend.position  = c(-0.1, -0.15),
          #legend.direction = "horizontal",
          
          # legend.key.size = unit(0.3, 'cm'),
          # legend.key.width= unit(1,'cm'),
          
          legend.title = element_text(face = 'bold', vjust = 1,
                                      color = 'black',
                                      size = 8)
          # legend.text= element_text(angle = 15,
          #                           size = unit(7, 'cm'))
    )
  return(g1)
}


# TM_PA_KEGG_pathview
TM_PA_KEGG_pathview = function(geneList, pathway.id, species,
                               study_title, pdf_ver = F){
  
  ori_wd = getwd()
  
  target_dir = paste0(ori_wd, "/sub_data/", study_title, "/pathview/", pathway.id, "/")
  
  if(!dir.exists(target_dir)){
    dir.create(target_dir, recursive = T)
  }
  
  setwd(target_dir)
  mmu03030 <- pathview(gene.data  = geneList,
                       pathway.id = pathway.id,
                       species    = species, kegg.native = pdf_ver,
                       limit      = list(gene=max(abs(geneList)), cpd=1))
  setwd(ori_wd)
  return("Success")
}
