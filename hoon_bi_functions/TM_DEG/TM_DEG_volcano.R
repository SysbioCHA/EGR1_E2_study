pacman::p_load(ggrepel, ggtext)

TM_DEG_volcano = function(DEG_res, x, y, gene_lab = NULL,
                          x_cut = 0.6, y_cut = 0.05, show_genename = F,
                          p_size = 1, l_size=3, title = "", Res4Genes=NULL,
                          x_lab="log2FC", y_lab="-log10(adj.p-value)"){
  args <- lapply(list(x,y), function(x) if (!is.null(x)) sym(x))
  
  de <- DEG_res
  de$diffexpressed <- ""
  de$diffexpressed[de[[x]] >= x_cut & de[[y]] <= y_cut] <- "Up regulated"
  de$diffexpressed[de[[x]] <= -x_cut & de[[y]] <= y_cut] <- "Down regulated"
  
  n_up = sum(de$diffexpressed == "Up regulated")
  n_down = sum(de$diffexpressed == "Down regulated")
  
  de$rank20 = ""
  de$rank20[ sort(de[[y]], index.return =T, decreasing = F)$ix[1:50] ] = "top"
  
  mycolors <- c("#6E95E6", "#A43029", "#BEBEBE")
  names(mycolors) <- c("Down regulated", "Up regulated", "")
  de$delabel <- NA
  
  # de$delabel[de$rank20 != ""] <- de[[gene_lab]][de$rank20 != ""]
  de$logP = -log10(de[[y]])
  
  if(show_genename){
    de$delabel[de$diffexpressed != ""] <- de[[gene_lab]][de$diffexpressed != ""]
    if(!is.null(Res4Genes)){
      de$delabel[ !de$delabel %in% Res4Genes ] = NA
      print(unique(de$delabel))
    }
    
    g1=ggplot(data=de,aes(x=!! sym(x),y=logP)) + theme_bw() + #xlim(-4,4) + ylim(0,5) +
      geom_point(aes(col=diffexpressed),size=p_size) +
      scale_colour_manual(values = mycolors, name = "") +
      geom_text_repel(data = (de %>% filter(!!sym(y) <= y_cut) %>% 
                        filter((!!sym(x) <= -x_cut) | (!!sym(x) >= x_cut))),
                      aes(label=delabel), show.legend = F, size=l_size,
                      box.padding = .5, max.overlaps = Inf)
    
  } else {
    g1 = ggplot(data=de,aes(x=!! sym(x),y=logP)) + theme_bw() + #xlim(-4,4) + ylim(0,5) +
      geom_point(aes(col=diffexpressed),size=p_size) +
      scale_colour_manual(values = mycolors, name = "")
  }
  g1 = g1 + labs(title=paste0(title, "<span style='font-size:14pt'> (",
                              "<span style='color:#6E95E6;font-size:12pt'>",
                              "Down: ",n_down,"</span>  ",
                              "<span style='color:#A43029;font-size:12pt'>",
                              "Up: ",n_up,"</span>)</span>"),
                 x=x_lab, y=y_lab) +
    theme(plot.title = element_markdown(hjust = 0.5, vjust = 0.5),
          axis.text=element_text(size=13), axis.title=element_text(size=14))
  return(g1)
}
