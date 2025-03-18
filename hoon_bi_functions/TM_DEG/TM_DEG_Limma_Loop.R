TM_DEG_Limma_loop = function(pval_cut=0.05, log2fc_cut=0.6,
                             INFO, GE, GENEANNO, loop,
                             Group, case, control,
                             var_cut = 0.5, sample_colname = NULL,
                             symAnno = NULL
){
  
  rt_li = list(cutoff = c(pval_cut, log2fc_cut))
  rt_li$limma = sapply(unique(INFO[[loop]]), function(lp){
    tmp_info = INFO %>% filter(!! sym(loop) == lp)
    print(tmp_info)
    tmp_1 = TM_DEG_Limma(GE, tmp_info, Group = Group, case = case, control = control,
                         var_cut = var_cut, sample_colname= sample_colname)
    tmp_2 = merge(tmp_1, GENEANNO, by = "probeID") %>% 
      group_by(entrez) %>% dplyr::filter(P.Value == min(P.Value)) %>% ungroup()
    return(tmp_2)
  },simplify = F)
  
  rt_li$stats = sapply(unique(INFO[[loop]]), function(lp){
    res = TM_DEG_stats(rt_li$limma[[lp]], pval_cut, log2fc_cut, "entrez")
  }, simplify = F)
  rt_li$stats_bar = TM_DEG_stats_barplot(rt_li$stats)
  
  if(!is.null(symAnno)){
    rt_li$statsSym = sapply(rt_li$stats, function(x){
      DEGs = sapply(x$DEGs, function(y) {
        y = y[y %in% symAnno$entrez]
      }, simplify = F)
      DEG_stats = x$DEG_stats
      return(list(DEGs = DEGs, DEG_stats = DEG_stats))
    }, simplify = F)
    
    rt_li$statsSym_bar = TM_DEG_stats_barplot(rt_li$statsSym)
  }
  
  return(rt_li)
}


