pacman::p_load(clusterProfiler)

#
TM_PA_ORA_reduce = function(genes, org="Hs", keyType = "ENTREZID", ont = "BP",
                            pAdjustMethod="BH", pvalueCutoff=.05, qvalueCutoff=.05,
                            minGSSize = 25, maxGSSize = 500, readable = T,
                            sim_threshold = 0.7){
  res = TM_PA_ORA(genes=genes, org=org, keyType=keyType,ont=ont,
                  pAdjustMethod=pAdjustMethod,pvalueCutoff=pvalueCutoff,
                  qvalueCutoff=qvalueCutoff,minGSSize=minGSSize,maxGSSize=maxGSSize,
                  readable = readable)
  res$GOreduced = TM_PA_utils_reduceSim(resGO=res$GO, org=org, threshold=sim_threshold)
  return(res)
}


#
TM_PA_ORA = function( genes, universe = NULL, org = "Hs", # Hs Mm
                      keyType = 'ENTREZID', # ENSEMBL
                      # One of "BP", "MF", and "CC" subontologies, or "ALL" for all three
                      ont= 'BP',
                      # "BH", "BY", "fdr", "none"
                      pAdjustMethod =  "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                      minGSSize = 10, maxGSSize = 500, readable = T,
                      doGO = T, doGO_sim = T # doKEGG = F
){
  if(org=="Hs"){
    require(org.Hs.eg.db)
    orgdb = org.Hs.eg.db
  } else if(org=="Mm"){
    require(org.Mm.eg.db)
    orgdb = org.Mm.eg.db
  } else {
    print("Wrong 'org'")
    return(F)
  }
  
  perform_sts = TRUE %in% c(doGO)
  if(!perform_sts){
    print("No ORA is allowed.")
    return(FALSE)
  } else {
    res_li = list()
  }

  if(doGO){
    # GO
    pw_GO <- enrichGO(gene = genes,
                      OrgDb = orgdb,
                      ont = ont,
                      keyType = keyType,
                      universe = universe,
                      minGSSize = minGSSize,
                      maxGSSize = maxGSSize,
                      
                      pAdjustMethod = pAdjustMethod,
                      pvalueCutoff = pvalueCutoff,
                      qvalueCutoff = qvalueCutoff,
                      readable = readable)
    
    if(doGO_sim){
      if(nrow(pw_GO) > 1){
        pw_GO = enrichplot::pairwise_termsim(pw_GO)
        pw_GO_prc = clusterProfiler::simplify(pw_GO, cutoff= 0.7,
                             by="p.adjust", select_fun = min, measure = "Wang")
      } else {
        pw_GO_prc = pw_GO
        
      }
    } else {
      pw_GO_prc = pw_GO
    }
    
    pw_GO_prc = pw_GO_prc %>% dplyr::mutate(p.adjust_log = -log10(p.adjust))
    res_li[["GO"]] = pw_GO_prc
    
  }
  
  # if(doKEGG){
  #   # KEGG  
  #   pw_KEGG  = enrichKEGG(gene = genes,
  #                         organism = organism,
  #                         universe = universe,
  #                         
  #                         minGSSize = minGSSize,
  #                         maxGSSize = maxGSSize,
  #                         
  #                         pAdjustMethod = pAdjustMethod,
  #                         pvalueCutoff = pvalueCutoff,
  #                         qvalueCutoff = qvalueCutoff)
  #   pw_KEGG = pw_KEGG %>% dplyr::mutate(p.adjust_log = -log10(p.adjust))
  #   res_li[["KEGG"]] = pw_KEGG
  #   
  # }
  # MKEGG ( KEGG Module is a collection of manually defined function units. In some situation, KEGG Modules have a more straightforward interpretation. )
  # pw_MKEGG <- enrichMKEGG(gene = genes,
  #                         organism = organism,
  #                         
  #                         minGSSize = minGSSize,
  #                         maxGSSize = maxGSSize,
  #                         
  #                         pAdjustMethod = pAdjustMethod,
  #                         pvalueCutoff = pvalueCutoff,
  #                         qvalueCutoff = qvalueCutoff)
  return(res_li)
}
