# TM_sub_data_EGR1_immature_mouse!
pacman::p_load(GEOquery, limma, umap, sva, tidyverse)

# 초기 다운로드 파일
# GSE2195LIST=getGEO(filename ="/home/hoon/R_prjs/EGR1_4th/contents/GSE2195_family.soft")
# saveRDS(GSE2195LIST, "./sub_data/EGR1_study/immature_mouse/GSE2195LIST.rds")
gse2195_prc1 = function(path){
  GSE2195LIST = readRDS(path)
  gse2195li = sapply(GSE2195LIST@gsms, function(x){
    data_title = x@header$title
    data_title = gsub(pattern = " ",replacement = "_", data_title)
    df = x@dataTable@table %>% dplyr::select("ID_REF", "VALUE", "ABS_CALL") %>% 
      setnames(old = colnames(.)[-1],
               new = paste0(data_title, "_", colnames(.)[-1]))
    return(list(df) %>% setNames(data_title))
  }, USE.NAMES = F)
  
  gse2195df = gse2195li %>% purrr::reduce(full_join, by='ID_REF')
  
  colnames(gse2195df) = gsub(pattern = "E2", replacement = "IM", colnames(gse2195df))
  colnames(gse2195df) = gsub(pattern = "hr", replacement = "h", colnames(gse2195df))
  return(gse2195df)
}

# 배치이펙트 및 Outlier 발견 하기 위해 Step13B로 Raw matrix를 보내줄 코드
# (if detected, then P for Present; if not detected, then A for Absent; M is Marginal for borderline detection)
gse2195_PMAprc = function(df){
  tmp1 = df %>% dplyr::select(contains("Value"))
  rownames(tmp1) = df$ID_REF
  
  ProbeCall = df %>% dplyr::select(contains("ABS_CALL")) %>% 
    # mutate(score = apply(., 1, function(x) return((ncol(.) - sum(x =="A"))))) %>% 
    mutate(score = apply(., 1, function(x) return(sum(x =="P")))) %>% 
    mutate(meanv = rowMeans(tmp1)) %>% 
    mutate(ID_REF = df$ID_REF) %>% column_to_rownames("ID_REF") %>% 
    mutate(probeID = rownames(.))
  tmp2 = ProbeCall %>% filter(score == 0)

  tmp3 = tmp1[!rownames(tmp1) %in% rownames(tmp2),] %>% as.matrix()
  tmp3[which(tmp3 <= 0)] <- NaN
  
  SSS = which(sapply(str_split(colnames(tmp3),"_"),function(x) length(x))==3)
  
  GSE2195exp = log2(tmp3) %>% na.omit(.) %>% as.data.frame() %>%
    setnames(old=colnames(.),new=str_split(colnames(.),"_VALUE",simplify=T)[,1]) %>%
    setnames(old=colnames(.)[min(SSS):max(SSS)],
             new=paste("HUN",colnames(.)[min(SSS):max(SSS)],sep="_"))
  
  k = strsplit( colnames(GSE2195exp), split = '_')
  GSE2195info = data.frame(sample = colnames(GSE2195exp),
                           condition = sapply(k, function(x) {paste0(x[2],"_",x[3])}),
                           type = sapply(k, function(x) {x[2]}),
                           hours = sapply(k, function(x) {gsub("hr", "h", x[3])}),
                           batch = sapply(k, function(x) {x[1]})) %>% 
    # mutate(hours = ifelse(type == "AO", "0h", hours)) %>% 
    column_to_rownames('sample')
  
  return(list("GE"=GSE2195exp, "INFO"=GSE2195info, "PROBECALL"=ProbeCall))
}

