pacman::p_load("rrvgo")

# go_clustering_service
go_clustering_service = function(hsGO, all_cover_go){
  
  hsGO_mini = hsGO %>% filter(ID %in% all_cover_go)
  hsGO_mini = enrichplot::pairwise_termsim(hsGO_mini,
                                           showCategory = nrow(hsGO_mini))
  hsGO_tree = enrichplot::treeplot(hsGO_mini, showCategory = nrow(hsGO_mini),
                                   cluster.params = list(n = round(nrow(hsGO_mini)/15)))
  hsGO_cluster = hsGO_tree$data %>% filter(!is.na(label)) %>%
    dplyr::select(group, label) %>% 
    dplyr::rename(cluster = group, Description = label) %>% 
    mutate(cluster = as.numeric(substr(cluster, 9, nchar(.))))
  return(hsGO_cluster)
}


# TM_PA_utils_GO_clustering
TM_PA_utils_GO_clustering = function(PA_ORA_results, hsGO){
  
  all_upgo = sapply(names(PA_ORA_results), function(x){
    ids = PA_ORA_results[[x]]$Up$GO %>% as.data.frame() %>% .$ID
    return(ids)
  }, simplify = F) %>% unlist() %>% unique()
  
  all_dwgo = sapply(names(PA_ORA_results), function(x){
    ids = PA_ORA_results[[x]]$Down$GO %>% as.data.frame() %>% .$ID
    return(ids)
  }, simplify = F) %>% unlist() %>% unique()
  
  hsGO_cluster_up = go_clustering_service(hsGO, all_upgo)
  hsGO_cluster_dw = go_clustering_service(hsGO, all_dwgo)

  return(list(Up=hsGO_cluster_up, Down=hsGO_cluster_dw))
}


# TM_PA_utils_reduceSim
TM_PA_utils_reduceSim = function(resGO, org = "Mm", ont = "BP", threshold = 0.7){
  if(nrow(resGO) < 2){
    print("No enriched terms detected.")
    return(F)
  }
  if(org=="Hs"){
    require(org.Hs.eg.db)
    orgdb = "org.Hs.eg.db"
  } else if(org=="Mm"){
    require(org.Mm.eg.db)
    orgdb = "org.Mm.eg.db"
  } else {
    print("Wrong 'org'")
    return(F)
  }
  
  simMatrix_m <- calculateSimMatrix( resGO$ID, orgdb = orgdb,
                                     ont= ont, method = 'Wang')

  simMatrix_m = simMatrix_m[rownames(simMatrix_m) %in% resGO$ID,
                            colnames(simMatrix_m) %in% resGO$ID]
  scores <- setNames(-log10(resGO$p.adjust), resGO$ID)
  reducedTerms08=reduceSimMatrix(simMatrix_m,scores,threshold=threshold,orgdb=orgdb)
  
  res <- resGO %>% filter(ID %in% unique(reducedTerms08$parent))
  return(list(Res = res, reducedTerms = reducedTerms08))
}


# TM_PA_utils_coder
TM_PA_utils_coder = function(df, order){
  codes = sapply(unique(df$Description), function(x){
    
    
    
    tmp1 = df %>% filter(Description == x) %>% .$cond
    code = (order %in% tmp1)
    code_prc = paste(substr(as.character(code), 1, 1), collapse = "")
    return(code_prc)
  })
  
  coder_res = data.frame(Description = names(codes), code = codes)
  return(coder_res)
}


# TM_PA_utils_mmuhsa_converter
TM_PA_utils_mmuhsa_converter = function(input_v,
                                        input_type = "entrez", input_species = "mouse",
                                        output_type = "entrez", output_species = "human",
                                        study_title = NULL){
  # type: entrez or symbol
  # study_title: dir name
  
  mhg_url = "http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
  
  if(is.null(study_title)){
    mouse_human_genes = fread(mhg_url) %>% 
      setnames(., make.names(colnames(.)))
               
  } else {
    target_dir = paste0(getwd(), "/sub_data/", study_title, "/mmuhsa_map/")
    target_rds = paste0(target_dir, "HOM_MouseHumanSequence.rds")
    
    if(file.exists(target_rds)){
      mouse_human_genes = readRDS(target_rds)
      
    } else {
      dir.create(target_dir)
      mouse_human_genes = fread(mhg_url) %>% 
        setnames(., make.names(colnames(.)))
      saveRDS(mouse_human_genes, target_rds)
      
    }
  }
  
  Split1 <- split(mouse_human_genes, mouse_human_genes$Common.Organism.Name) %>% 
    setNames(., substr(names(.), 1, 5))
  
  select_cols = c("DB.Class.Key","Common.Organism.Name",
                  "NCBI.Taxon.ID","Symbol","EntrezGene.ID")
  
  Split2 = sapply(names(Split1), function(species){
    old_col = sym(paste0(species, "_DB.Class.Key"))
    tmp = Split1[[species]] %>%
      dplyr::select(all_of(select_cols)) %>%
      setNames(paste0(species, "_", names(.))) %>% 
      dplyr::rename( "DB.Class.Key" := !! old_col)
    return(tmp)
  }, simplify = F)
  
  if(input_type == "entrez"){
    input_query = paste0(input_species, "_EntrezGene.ID")
  } else if(input_type == "symbol"){
    input_query = paste0(input_species, "_Symbol")
  }
  if(output_type == "entrez"){
    output_query = paste0(output_species, "_EntrezGene.ID")
  } else if(output_type == "symbol"){
    output_query = paste0(output_species, "_Symbol")
  }
  
  input_to_DBkey = sapply(input_v, function(x){
    ipqry = Split2[[input_species]] %>% .[[input_query]]
    
    class_key=(Split2[[input_species]] %>% 
                 .[ str_to_upper(ipqry) == str_to_upper(x), ])[["DB.Class.Key"]]
    if(identical(class_key, integer(0))){
      # cat("\n"); print(paste0("notice: ", x, " has no class_key, it would be skipped."))
      return(NA)
    } else {
      return(class_key)
    }
  })
  
  input_to_DBkey_status = table(sapply(input_to_DBkey, function(x) length(x)))
  cat("\n"); print("DB key Status. ('1' mean 1:1 matches.) ")
  print(input_to_DBkey_status)
  
  input_to_DBkey_prc = lapply(input_to_DBkey, function(x){
    if(length(x) > 1){
      x = x[1]
    }
    return(x)
  }) %>% unlist()
  
  cat("\n"); print("Multiple DB Key correction!")
  mapped_DF = data.frame(A = names(input_to_DBkey_prc),
                         B = input_to_DBkey_prc) %>% 
    setNames(c(input_query, "DB.Class.Key")) %>% 
    merge(., Split2[[output_species]], by="DB.Class.Key")
  
  cat("\n"); print("Duplicate Check.")
  table(duplicated(mapped_DF[[input_query]]))
  # mapped_DF_rmdup = mapped_DF[-which(duplicated(mapped_DF[[input_query]])),]
  # table(duplicated(mapped_DF_rmdup[[input_query]]))
  
  if(any(grepl("[.]y$",colnames(mapped_DF)))){
    print("Duplicated col detected, it would be corrected ..")
    mapped_DF = subset(mapped_DF,T,-grep("[.]y$",colnames(mapped_DF)))
    tmp_i = grep("[.]x$",colnames(mapped_DF))
    colnames(mapped_DF)[tmp_i] = gsub(".x","",colnames(mapped_DF)[tmp_i])
  }
  return(mapped_DF)
}

#TM_PA_utils_jaccard
TM_PA_utils_jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


#TM_PA_hoonvenn
TM_PA_hoonvenn <- function(input, include, exclude = NULL) {
  if(!is.null(exclude)){
    include = setdiff(include, exclude)
  }
  res = Reduce(intersect, input[include])
  
  if(!is.null(exclude)){
    res = setdiff(res, unique(unlist(input[exclude])))
    
  }
  return (res)
}

