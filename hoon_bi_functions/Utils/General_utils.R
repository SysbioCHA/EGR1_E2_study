rtnP = function(dirs,filename,type="png",recursive=F){
  if(length(dirs)>1){
    dirs = paste0(dirs, collapse = "/")
  }
  
  if(!dir.exists(dirs)){
    if(recursive){
      print(f("Create {dirs}", environment()))
      dir.create(dirs,recursive = T)
    } else {
      print("No such dirs exist.")
      return(FALSE)
    }
  }
  if(length(filename)>1){
    filename = paste0(filename, collapse = "")
  }
  rtn = paste(dirs, paste(filename,type, sep = "."), sep = "/")
  return(rtn)
}

categorize = function(x, map_list){
  tmp1 = AnnotationDbi::unlist2(map_list)
  tmp2 = setNames(names(tmp1), tmp1)
  sapply(x, function(y) {
    if(y %in% names(tmp2)){
      tmp2[y]
    } else {y}
    })
}

color_map = function(x, palette = "Set3"){
  if(length(unique(x)) > 12){
    return(FALSE)
  } else {
    return(setNames(brewer.pal(n=length(unique(x)), name = palette),unique(x)))
  }
}

pairing_GE_INFO = function(GE, INFO){
  if(ncol(GE) == nrow(INFO)){
    if(all(colnames(GE) == rownames(INFO))){
      merge_li = list("GE"=GE, "INFO"=INFO)
    } else {
      return(FALSE)
    }
  } else {
    samps = intersect(colnames(GE), rownames(INFO))
    GE = GE %>% dplyr::select(all_of(samps))
    INFO = INFO %>% dplyr::filter(row.names(.) %in% samps)
    merge_li = list("GE"=GE, "INFO"=INFO)
  }
  return(merge_li)
}

# Ent to Sym
convert_id = function(org = "Hs", from, to, keys, help=F){
  if(help){
    print(keytypes(org.Hs.eg.db))
    return(T)
  }
  
  # org: Hs or Mm
  if(org == "Hs"){
    require(org.Hs.eg.db)
    org = org.Hs.eg.db
  } else if(org == "Mm"){
    require(org.Mm.eg.db)
    org = org.Mm.eg.db
  }
  columns = c(from, to)
  return(AnnotationDbi::select(org, keys=keys, columns=columns, keytype=from))
}




f = function(string, e) return(glue(string, .envir = e))

split_tibble <- function(tibble, column = 'col') {
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}