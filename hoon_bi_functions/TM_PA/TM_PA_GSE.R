pacman::p_load(clusterProfiler)

TM_PA_GSE = function( DEGres, # should be contain logFC, name이 gene, 값은 FC
                      organism = "mmu",
                      keyType = 'ENTREZID', # ENSEMBL
                      
                      # One of "BP", "MF", and "CC" subontologies, or "ALL" for all three
                      ont= 'BP',
                      
                      # "BH", "BY", "fdr", "none"
                      pAdjustMethod =  "BH",
                      pvalueCutoff = 0.05,
                      
                      minGSSize = 10,
                      maxGSSize = 500,
                      readable = T,
                      
                      doGO_sim = T){
  
  geneList = DEGres$logFC
  names(geneList) = DEGres$entrez
  geneList = sort(geneList,decreasing = T)
  
  if(organism == "mmu"){
    require(org.Mm.eg.db)
    org = org.Mm.eg.db
  } else if(organism == "hsa") {
    require(org.Hs.eg.db)
    org = org.Hs.eg.db
  }
  
  # GO
  pw_GO <- gseGO(geneList  = geneList,
                 OrgDb = org,
                 ont = ont,
                 keyType = keyType,
                 
                 minGSSize = minGSSize,
                 maxGSSize = maxGSSize,
                    
                 pAdjustMethod =  pAdjustMethod,
                 pvalueCutoff = pvalueCutoff)
  
  if(doGO_sim){
    if(nrow(pw_GO) > 1){
      pw_GO_prc = clusterProfiler::simplify(pw_GO, cutoff= 0.7,
                           by="p.adjust", select_fun = min, measure = "Wang")
      pw_GO_prc@result$Regulation = ifelse(pw_GO_prc@result$NES > 0, "Up", "Down")
    } else {
      pw_GO_prc = NULL
    }
    
  } else {
    pw_GO_prc = NULL
  }
  pw_GO@result$Regulation = ifelse(pw_GO@result$NES > 0, "Up", "Down")
  
  # KEGG
  pw_KEGG  = gseKEGG(geneList  = geneList,
                     organism = organism,
                        
                     minGSSize = minGSSize,
                     maxGSSize = maxGSSize,
                     
                     pAdjustMethod = pAdjustMethod,
                     pvalueCutoff = pvalueCutoff)
  pw_KEGG@result$Regulation = ifelse(pw_KEGG@result$NES > 0, "Up", "Down")
  # MKEGG
  # pw_MKEGG <- gseMKEGG(geneList = geneList,
  #                      organism = organism,
  #                      
  #                      minGSSize = minGSSize,
  #                      maxGSSize = maxGSSize,
  #                      
  #                      pAdjustMethod = pAdjustMethod,
  #                      pvalueCutoff = pvalueCutoff)
  
  
  return(list("GO"= pw_GO, "GO_sim" =pw_GO_prc, "KEGG"=pw_KEGG))
}


