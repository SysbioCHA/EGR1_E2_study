TM_DEG_stats = function(DEG_res, pval_cut = 0.05, logFC_cut = 0.6, id="entrez"){
  
  tmp = DEG_res %>% filter(adj.P.Val <= pval_cut)
  up_de = tmp %>% filter(logFC >= logFC_cut) %>% .[[id]]
  dw_de = tmp %>% filter(logFC <= -logFC_cut) %>% .[[id]]
  return(list(DEG_stats = list(up=length(up_de), down=length(dw_de)), 
              DEGs = list(up=up_de, dw=dw_de)))
}

TM_DEG_stats_barplot = function(stat_res){
  
  ggplot(sapply(stat_res, function(x) sapply(x$DEGs,function(y) length(y))) %>% 
           as.data.frame() %>% rownames_to_column("reg") %>% 
           pivot_longer(col=-reg,names_to="tp",values_to="degs") %>% group_by(tp) %>% 
           mutate('csum'=unlist(by(data=degs,INDICES=tp,FUN=cumsum))) %>% ungroup() %>% 
           mutate('ypos'=ifelse((csum-degs)==0, csum/2, degs/2+(csum-degs))),
         aes(x = factor(tp, levels = rev(unique(tp))), y = degs, fill = reg), ) +
    geom_bar(stat = "identity", width=.7) + labs(x="",y="the number of DEGs") +
    scale_fill_manual(values = c("#6E95E6", "#A43029"), name = "") +
    geom_text(aes(x=tp, y=ypos, label=degs), colour="white", size=3.5) +
    theme(axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),
          axis.title.x=element_text(size=13),legend.position = "none",
          axis.title.y=element_text(size=13)) + coord_flip()
}





