pacman::p_load(illuminaMousev2.db)

#@@@@@ TM_EGR1_processing_1 @@@@@#
TM_EGR1_processing_1 = function(load_path){
  
  # Microarray Gene exp : probe id 별 전체 샘플의 signal
  if(is.character(load_path)){
    HT = fread(load_path) %>% dplyr::select(`Probe ID`, contains("Signal")) %>% 
      dplyr::rename(probeID = "Probe ID")
  } else {
    HT = load_path
  }
  
  # HT의 colname을 분석하기 좋게 바꾸고, probe id는 rowname으로 설정
  prc1_GE = data.frame(as.matrix(HT %>% dplyr::select(-probeID)),
                   row.names = HT$probeID) %>% 
    setNames(gsub(pattern = "\\.", replacement = "_", names(.))) %>% 
    setNames(substr(names(.), 1, nchar(names(.))-7)) %>% 
    setNames(sapply(strsplit( names(.), split = '_'), function(x) {
      paste(x[1],paste( x[2],'h',sep = ""),x[3],sep= "_")}))
  
  k = strsplit( colnames(prc1_GE), split = '_')
  prc1_INFO = data.frame(sample = colnames(prc1_GE),
                         condition = substr(colnames(prc1_GE),
                                            1, nchar(colnames(prc1_GE))-2),
                         type = sapply(k, function(x) {x[1]}),
                         hours = sapply(k, function(x) {x[2]})) %>% 
    column_to_rownames('sample')
  
  return(list("GE"=prc1_GE, "INFO"=prc1_INFO))
}


#@@@@@ TM_EGR1_processing_2 @@@@@#
TM_EGR1_processing_2 = function(load_path){
  
  # get Detection data and remove outlier
  HT_rm_outlier = fread(load_path) %>% 
    dplyr::select(`Probe ID`, contains("Signal"), contains("Detection")) %>% 
    dplyr::rename(probeID = "Probe ID") %>% 
    setNames(gsub(pattern = "-",replacement = "_",names(.))) %>% 
    dplyr::select(!contains("KO_3_A")) %>% as.data.frame() # <- KO 3 A.   Sample QC Fail
  Detection = grepl("p_value", colnames(HT_rm_outlier))
  
  # probes display detect p val < 0.05 in least 3 of 17
  isexpr = rowSums(HT_rm_outlier[,Detection] <= 0.05) >= 3
  HT_rm_low = HT_rm_outlier[isexpr,] %>% dplyr::select(!contains("Detection"))
  
  # Probe quality
  x <- illuminaMousev2PROBEQUALITY
  xx <- as.list(x[mappedkeys(x)])
  
  table(unlist(xx))
  xxx = data.frame(t(as.data.frame(xx)))
  probe_quality = data.frame(probeID = rownames(xxx),
                             quality = xxx[[1]],
                             entrez = unlist(
                               mget(x=rownames(xxx),
                                    envir=illuminaMousev2ENTREZREANNOTATED))) %>%
    filter(!quality %in% c("No match", "Bad")) %>% filter(!is.na(entrez))
  
  table(probe_quality$quality); table(duplicated(probe_quality$entrez))
  
  # probe quality에 score를 매김 good- good* -... -perfect*** 순으로, 좋은 등급 == 큰 숫자
  pq_temple = unique(probe_quality$quality)[c(2,6,5,1,4,3)]; print(pq_temple)
  probe_quality$qualscore = sapply(probe_quality$quality,
                                   function(x) return(which(pq_temple == x)))
  probe_quality2 = probe_quality %>% group_by(entrez) %>% filter(qualscore == max(qualscore))  
  
  #entrez ID 당 quality가 가장 높은 probe 선택
  table(duplicated(probe_quality2$entrez))
  
  # Merge Good quality "probe_quality2" / outlier & low expression removed "HT_rm_low"
  HT_final = HT_rm_low %>% inner_join(probe_quality2) %>% group_by(entrez) %>%
    filter(qualscore == max(qualscore)) %>% ungroup() %>% 
    dplyr::select(probeID, entrez, contains("signal"))
  
  prc_resli = TM_EGR1_processing_1(HT_final %>% dplyr::select(-entrez))
  prc_resli$INFO = prc_resli$INFO
  GeneAnno_ProbeToEntrez = HT_final %>% dplyr::select(probeID, entrez)
  
  return(list("GE"=prc_resli$GE,"INFO"=prc_resli$INFO,"GENEANNO"=GeneAnno_ProbeToEntrez))
}


# 필요 앱 : string app,
# 1. File > import > import network from Public Databases
# 2. (설정값) STRING: protein query, Mis musculus, Confidence ( score ) cutoff = 0.7
# 3. Identifiers 넣고 Import
# 
# 4. 우측 사이드의 Nodes 탭에서, Singletons 해제
# 
# 5. Import table from file ( logFC 등 가져오기 )
# - Node table columns 설정
# 
# 6. 좌측 사이드
# - Fill color : Map. > Log2FC > Continuous Mapping
# - Image/Chart~ 는 모두 휴지통
# - Label : Map. > display name > Passthrough mapping
# 
# 7. App > clustermaker cluster network > fast greedy > Attribute: fastgreedycluster & create new clustered network

  