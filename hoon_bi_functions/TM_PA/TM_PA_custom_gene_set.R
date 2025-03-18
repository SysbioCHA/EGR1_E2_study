pacman::p_load(GSVA, msigdbr)

# TM_PA_custom_geneset_analysis
TM_PA_custom_geneset_analysis = function(mat, method="ssgsea", gs, group_labels,
                                         group_levels = NULL){
# <<<<<<< HEAD
  
  if(method != "gsva"){
    print("The GSVA package, altered there functions. plz, recognize it and fix the codes.")
  }
  
  # old
  # ssgsea_score = gsva(mat, gs, method = method)
  
  # new
  if(method == "gsva"){
    ssgsea_score = gsvaParam(mat, sapply(gs, function(x) as.character(x)), maxDiff=TRUE)
    ssgsea_score = gsva(ssgsea_score)
  }
# =======
#   ssgsea_score = gsva(mat, gs, method = method)
# >>>>>>> parent of 5a66854 (ck)
  
  ssgsea_df = TM_PA_custom_geneset_plot(ssgsea_score, group_labels)
  # ssgsea_df = data.frame(group_labels,
  #                        as.numeric(paste(apply(ssgsea_score, 1, function(x) return(x)))),
  #                        rep(rownames(ssgsea_score), each=length(group_labels))) %>%
  #   setNames(c("group", "Score", "Pathway"))
  
  return(list(mat=ssgsea_score, boxdf = ssgsea_df))
}

# TM_PA_custom_geneset_plot
TM_PA_custom_geneset_plot = function(ssgsea_score, group_labels, sample_label = NULL, 
                                     gs_name= NULL){
  if(!is.null(sample_label)){
    ssgsea_score = subset(ssgsea_score,T,sample_label)
  }
  if(!is.null(gs_name)){
    ssgsea_score = ssgsea_score[gs_name, ]
  }
  
  ssgsea_df = data.frame(group_labels,
                         as.numeric(paste(apply(ssgsea_score, 1, function(x) return(x)))),
                         rep(rownames(ssgsea_score), each=length(group_labels))) %>%
    setNames(c("group", "Score", "Pathway"))
  
  # df2 = subset(df,rownames(df)==gs_name, sample_label)
  # df3 = data.frame(Pathway=gs_name, Score = as.numeric(df2),
  #                  Group = str_split(sample_cond, "_", simplify = T)[,1],
  #                  hours = str_split(sample_cond, "_", simplify = T)[,2])
  return(ssgsea_df)
}

# TM_PA_custom_msigDBr
TM_PA_custom_msigDBr_query = function(x, species = "Mus musculus", category = "H" ){
  # msigdbr_show_species()
  # "Homo sapiens"                   
  # "Mus musculus"
  # H: hallmark gene sets
  # C1: positional gene sets
  # C2: curated gene sets
  # C3: motif gene sets
  # C4: computational gene sets
  # C5: GO gene sets
  # C6: oncogenic signatures
  # C7: immunologic signatures
  res = msigdbr(species = species, category = category) %>% 
    dplyr::select(gs_name, entrez_gene)
  
  # analysis func
  print("Analysis Methods: enricher(ORA), GSEA(GSE)")
  print(">  enricher(gene, TERM2GENE= query_res)")
  print(">  GSEA(geneList, TERM2GENE= query_res)")
  return(res)
}
