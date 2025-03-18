# main script

# rm(list=ls());gc()
pacman::p_load(dplyr, data.table, RColorBrewer, pathview, dendextend, ggplot2, forcats,
               stringr, factoextra, STRINGdb, rrvgo, patchwork, msigdbr, GSVA,
               cmapR, corrplot, colorspace, future, future.apply,magrittr)
pacman::p_load(sva, colorspace, factoextra, mgu74av2.db) # hcl_palettes(plot = TRUE)

source("./TM_sub_data/TM_sub_data_EGR1.R")
source("./TM_sub_data/TM_sub_data_EGR1_immature_mouse.R")
sapply(list.files("./hoon_bi_functions", recursive = T,full.names = T), function(x){
  source(x); print(x); return(T)})

###########################################################
#      Immature mouse ...
###########################################################
{
  # 0. primary
  gse2195_li=gse2195_PMAprc(
    gse2195_prc1("./sub_data/EGR1_study/immature_mouse/GSE2195LIST.rds"))
  par(mfrow=c(1,2));plotDensities(gse2195_li$GE, legend = F);boxplot(gse2195_li$GE, las=2)
  dev.off()
  # 1. primary: Heatmap
  gse2195_li$hm = TM_EDA_Heatmap(mat = gse2195_li$GE, selection_method = "IQR", width=20,
                                 selection_n = 300, top_annot_df = gse2195_li$INFO %>% 
                                   dplyr::select(-c("condition"))); gse2195_li$hm
  # 2. primary: correction
  gse2195_li$GE_combat=TM_PRC_Combat_batch_correction(GE=gse2195_li$GE,
                                                      batch=gse2195_li$INFO$batch)
  # 3. primary: Heatmap
  gse2195_li$hm2 = TM_EDA_Heatmap(mat = gse2195_li$GE_combat, selection_method = "IQR",
                                  selection_n = 300, top_annot_df = gse2195_li$INFO %>% 
                                    dplyr::select(-c("condition"))); gse2195_li$hm2
  rm(gse2195_li)
}


# 4. second: remove Outlier
gse2195_li_corrected=gse2195_PMAprc((
  gse2195_prc1("./sub_data/EGR1_study/immature_mouse/GSE2195LIST.rds")%>%
    dplyr::select(-(contains("L2_AO_8h")|contains("48h")|contains("72h")))))
gse2195_li_corrected$INFO$hours_num=as.numeric(
  gsub("h","",gse2195_li_corrected$INFO$hours))

gse2195_li_corrected$GENEANNO=data.frame(
  probeID=rownames(gse2195_li_corrected$GE),
  entrez=unlist(mget(x=rownames(gse2195_li_corrected$GE),envir=mgu74av2ENTREZID )),
  symbol=unlist(mget(x=rownames(gse2195_li_corrected$GE),envir=mgu74av2SYMBOL ))) %>% 
  filter(!is.na(entrez))

# 5. second: correction
gse2195_li_corrected$GE_combat = TM_PRC_Combat_batch_correction(
  gse2195_li_corrected$GE, gse2195_li_corrected$INFO$batch)

gse2195_li_corrected$GE_combat_entrez = gse2195_li_corrected$GE_combat %>% 
  mutate(probeID = rownames(.)) %>% merge(., gse2195_li_corrected$GENEANNO) %>% 
  merge(.,gse2195_li_corrected$PROBECALL) %>% group_by(entrez) %>% 
  filter(score==max(score)) %>% filter(meanv==max(meanv)) %>% ungroup() %>% 
  dplyr::select(-contains("ABS_CALL"), -c("score", "meanv", "probeID")) %>%
  column_to_rownames("entrez")

if(!B_source_mode){
  gse2195_li_corrected$hm=draw(
    TM_EDA_Heatmap(mat=gse2195_li_corrected$GE_combat,selection_method="IQR",
                   selection_n=300,top_annot_df=gse2195_li_corrected$INFO %>% 
                     dplyr::select(-c("condition"))),
    annotation_legend_side="left", heatmap_legend_side="left")
}

gse2195_li_corrected$GE_entrez_by_hours=sapply(
  unique(gse2195_li_corrected$INFO$condition), function(x){
    mg_tgs = gse2195_li_corrected$INFO %>% filter(condition == x) %>% rownames(.)
    res=rowMeans(gse2195_li_corrected$GE_combat_entrez %>% dplyr::select(all_of(mg_tgs)))
    return(res)})

gse2195_li_corrected$GE_ent_logFC=sapply(c("2h","4h","24h"),function(t) {
  df = gse2195_li_corrected$GE_entrez_by_hours
  return((df[,paste0("IM_",t)]-df[,paste0("AO_",t)]))}) %>% as.data.frame()

gse2195_li_corrected$Col_fun=setNames(
  genv$h_col_fun(5),unique(gse2195_li_corrected$INFO$hours))

gse2195_li_corrected$cola1_uniq_hrs = HeatmapAnnotation(
  annotation_name_side="right",annotation_name_align=T, show_legend = T,
  annotation_name_gp=gpar(fontsize=8), simple_anno_size=unit(.3,"cm"),
  `E2 treatment` = unique(gse2195_li_corrected$INFO$hours),
  col=list(`E2 treatment`=gse2195_li_corrected$Col_fun),
  annotation_legend_param=list(
    `E2 treatment`=list(at=c(names(gse2195_li_corrected$Col_fun)))))

gse2195_li_corrected$cola1_uniq_hrs2 = HeatmapAnnotation(
  annotation_name_side="right",annotation_name_align=T, show_legend = T,
  annotation_name_gp=gpar(fontsize=8), simple_anno_size=unit(.3,"cm"),
  `E2 treatment` = unique(gse2195_li_corrected$INFO$hours)[c(2,3,5)],
  col=list(`E2 treatment`=gse2195_li_corrected$Col_fun),
  annotation_legend_param=list(
    `E2 treatment`=list(at=c(names(gse2195_li_corrected$Col_fun)[c(2,3,5)]))))

# 7. second: GE entrez ####
dim(gse2195_li_corrected$INFO)
dim(gse2195_li_corrected$GE)
gse2195_li_corrected$GE_entrez = gse2195_li_corrected$GE %>% 
  mutate(probeID = rownames(.)) %>% merge(., gse2195_li_corrected$GENEANNO) %>% 
  merge(.,gse2195_li_corrected$PROBECALL) %>% group_by(entrez) %>% 
  filter(score==max(score)) %>% filter(meanv == max(meanv)) %>% ungroup() %>% 
  dplyr::select(-contains("ABS_CALL"), -c("score", "meanv", "probeID")) %>%
  column_to_rownames("entrez") %>% select(-symbol)
dim(gse2195_li_corrected$GE_entrez)



# 9. merge: immature and EGR1 data ####
# The initial effects of E2 are rapid (4–6 hr) and involve the uptake of fluid resulting from hyperemia and vasodilation of uterine capillaries, 
# which cause ...
dim(prc2_resli$GE_entrez2)

merge_li=list(
  "merge"=merge(prc2_resli$GE_entrez2, gse2195_li_corrected$GE_entrez,
                by='row.names',all=F) %>% column_to_rownames("Row.names"),
  "INFO" =rbind(prc2_resli$INFO %>% mutate(batch = "EGR1"),
                gse2195_li_corrected$INFO) %>% 
    mutate(
      hour_g=ifelse(hours_num %in% c(0,1),"0~1h",
                    ifelse(hours_num %in% c(2,3),"2~3h",
                           ifelse(hours_num %in% c(4,6), "4~6h",
                                  ifelse(hours_num == 8, "8h", "24h"))))) %>% 
    mutate(hour_g=factor(hour_g,levels=c("0~1h", "2~3h", "4~6h", "8h", "24h"))) %>% 
    mutate("E2_hours" = ifelse(type == "AO", "0h", hours)) %>% 
    mutate(
      "E2-elapsed time" = ifelse(
        E2_hours == "0h", "No Estrogen treatment",
        ifelse(E2_hours %in% c("1h","2h","3h"), "Early I",
               ifelse(E2_hours %in% c("4h","6h","8h"), "Early II", "Late")))))

if(!B_source_mode){
  boxplot(merge_li$merge, boxwex=0.7, notch=T, outline=FALSE, las=2)
  par(mfrow=c(1,2))
  plotDensities(prc2_resli$GE_entrez, legend = F)
  plotDensities(gse2195_li_corrected$GE_entrez, legend = F);
  dev.off()
  table(colnames(merge_li$merge) == 
          c(rownames(prc2_resli$INFO), rownames(gse2195_li_corrected$INFO)))
  
  merge_li$hm=draw(TM_EDA_Heatmap(mat = merge_li$merge, selection_method = "IQR",
                                  selection_n = 300, top_annot_df = merge_li$INFO %>% 
                                    dplyr::select(-c("condition","hours_num","hours",
                                                     "E2-elapsed time", "E2_hours"))),
                   annotation_legend_side="left",heatmap_legend_side="left",
                   merge_legend=T)
}

# 10. merge: correction ####
length(colnames(merge_li$merge))
length(merge_li$INFO$batch)

merge_li$GE_combat = TM_PRC_Combat_batch_correction(merge_li$merge, merge_li$INFO$batch)
merge_li$GENEANNO = prc2_resli$GENEANNO_sym %>% 
  filter(entrez %in% rownames(merge_li$GE_combat))

# 11. merge: Heatmap ####
merge_li$Col_fn = setNames(genv$h_col_fun(4), unique(merge_li$INFO$`E2-elapsed time`))
merge_li$cola1=HeatmapAnnotation(
  annotation_name_side = "right", annotation_name_align=T,
  `E2 treatment` = merge_li$INFO$`E2-elapsed time`,
  `Mature status`= ifelse(merge_li$INFO$type %in% c("IM","AO"), "Immature","Mature"),
  `Egr1 genotype` = ifelse(merge_li$INFO$type == "KO", "Knockout", "Wildtype") ,
  col=list(`Egr1 genotype`=c("Wildtype"=genv$wt_col, "Knockout"=genv$ko_col),
           `Mature status`=c("Mature"="#606060", "Immature"=genv$im_col),
           `E2 treatment`=merge_li$Col_fn),
  annotation_legend_param = list(
    `E2 treatment`=list(at=c(names(merge_li$Col_fn)))))

merge_li$rola1=HeatmapAnnotation(
  show_legend = F, Hours=merge_li$INFO$hours_num, col = list(
    Hours = setNames(genv$h_col_fun(length(unique(merge_li$INFO$hours_num))),
                     sort(unique(merge_li$INFO$hours_num)))),
  annotation_name_side="right", annotation_legend_param=list(
    Hours=list(at=c(names(merge_li$INFO$hours_num)))))

if(!B_source_mode){
  merge_li$hm2=draw(TM_EDA_Heatmap(
    mat=merge_li$GE_combat, name="Z score",
    selection_method="IQR",selection_n=300,top_annot_df = merge_li$cola1,do_scale = T,
    column_dend_height=1, bottom_annot_df=merge_li$rola1,
    show_col_names = T,show_row_dend = F, row_names_gp = 10,column_names_rot = 90,
    column_labels=strsplit2(merge_li$INFO$condition,"_")[,2]))
  
  png(rtnP(genv$P2,"SFig2_hm_EGRIM","png"),width =13,height=6,units="in",res=1200)
  draw(merge_li$hm2,annotation_legend_side="left",heatmap_legend_side="left"); dev.off()

  
  par(mfrow=c(1,2))
  boxplot(merge_li$GE_combat,boxwex=0.7,notch=T,outline=FALSE,las=2)
  plotDensities(merge_li$GE_combat,legend=F)

}

# 12-1. merge:  GSVA ####
library(msigdbr); library(GSVA); library(cmapR); library(factoextra); library(corrplot)
if(file.exists("./sub_data/EGR1_study/msigdb/gmts.rds")){
  gsvali_mg = list(gmts = readRDS("./sub_data/EGR1_study/msigdb/gmts.rds"))
} else {
  dir.create("./sub_data/EGR1_study/msigdb/")
  gsvali_mg = list(gmts = sapply(c("H", "C2", "C5", "C6"), function(x) {
    msdb = msigdbr(species = "Mus musculus", category = x) %>% 
      dplyr::select(gs_name, gene_symbol, entrez_gene)
    msdb = sapply(unique(msdb$gs_name), function(x){
      rt = msdb %>% filter(gs_name == x) %>% .[["entrez_gene"]]
      return(rt)})
    return(msdb)}))
  saveRDS(gsvali_mg$gmts, file = "./sub_data/EGR1_study/msigdb/gmts.rds")
}

t_ = gsvali_mg$gmts; gsvali_mg$gmts_filt = list(
  H = sapply(t_$H, function(x) as.character(x))
); rm(t_)

# GSVA
gsvali_mg$INFO = merge_li$INFO %>% filter(type != "AO") %>% select(-batch)
gsvali_mg$GE = merge_li$GE_combat %>% select(all_of(rownames(gsvali_mg$INFO)))

gsvali_mg$gsva_res = future_sapply(gsvali_mg$gmts_filt, function(x){
  gsvapar = gsvaParam(as.matrix(gsvali_mg$GE), x, kcdf = "Gaussian", maxDiff=TRUE,
                      minSize = 10)
  TG = gsva(gsvapar) %>% as.data.frame()
  return(TG)}, simplify = F)


# 전체 샘플의 홀마크 클러스터링
gsvali_mg$Col_fn2=setNames(genv$h_col_fun(8),sort(unique(gsvali_mg$INFO$hours_num)))
gsvali_mg$cola2 = HeatmapAnnotation(
  Group = gsvali_mg$INFO$type,
  `E2 treatment (h)` = gsub("h","",gsvali_mg$INFO$E2_hours),
  # `E2 treatment (term)` = gsvali_mg$INFO$`E2-elapsed time`,
  col=list(Group=genv$Group_col,
           `E2 treatment (h)`=gsvali_mg$Col_fn2),
           # `E2 treatment (term)`=merge_li$Col_fn),
  annotation_legend_param = list(
    `E2 treatment (h)`=list(at=c(names(gsvali_mg$Col_fn2)))))#,
    #`E2 treatment (term)`=list(at=c(names(merge_li$Col_fn)))))

colnames(gsvali_mg$gsva_res$H)

if(!B_source_mode){
  png(rtnP(genv$P2,"SFig3_hm","png"),width=11,height=10, units="in", res=300)

  png(rtnP(genv$P2,"SFig3_hm","png"),width=9,height=10, units="in", res=300)
  draw(TM_EDA_Heatmap(
    mat=gsvali_mg$gsva_res$H, height=18, cluster_row=F,row_title_gp=0,width=12,
    row_names_gp=10, show_col_names=T,show_row_names=T, cluster_col=F, do_scale = F,
    top_annot_df=gsvali_mg$cola2,column_dend_height=1,name = "Z score (GSVA)",
    column_names_rot = 90, col_names_gp = 10,col_split = c(gsvali_mg$INFO$hour_g),
    heatmap_legend_param=list("legend_height"=unit(2,"cm"))),
    annotation_legend_side="left",heatmap_legend_side="left",merge_legend=T)
  dev.off()
  
  tmp=TM_EDA_PCA(gsvali_mg$gsva_res$H,scale=T,center=T,addEllipses=F,
                 label_size=4, pointsize=3.5, legend_label_size=12, show_label="all",
                 colors = genv$Group_col,
                 ellipse.type = "norm", # norm confidence
                 axis_label_size=9,axis_title_label_size=10,
                 group=gsvali_mg$INFO$type); tmp$plot; rm(tmp)
  
  tmp=TM_EDA_PCA(gsvali_mg$gsva_res$H %>% select(!contains("AO")),scale=T,center=T,
                 addEllipses=F,label_size=4, pointsize=3.5, legend_label_size=12, 
                 show_label="all",colors = genv$Group_col,
                 ellipse.type = "norm", # norm confidence
                 axis_label_size=9,axis_title_label_size=10,
                 group=gsvali_mg$INFO$type %>% .[. != "AO"]); tmp$plot; rm(tmp)
}


# 각 시간대별 !!!!!!!!
# apply(gsvali_mg$now$mat, 1, function(x) { return( mad(x) )}) %>% 
#   sort(decreasing = T) %>% names(.)
# tmp = gsvaFilt %>% sort(decreasing = T) %>% names(.)
# now 에 대한 PCA
print(unique(gsvali_mg$INFO$hour_g))
if(T){
  gsvali_mg$now$tp= unique(gsvali_mg$INFO$hour_g)[4] # 바뀌어야 하는 것.
  
  gsvali_mg$now$tp
  gsvali_mg$now$INFO=gsvali_mg$INFO %>% filter(
    `E2-elapsed time`!="No Estrogen treatment") %>%
    filter(hour_g==gsvali_mg$now$tp)
  gsvali_mg$now$mat=gsvali_mg$gsva_res$H %>% select(all_of(rownames(gsvali_mg$now$INFO)))
  gsvali_mg$now$p=TM_EDA_PCA(
    as.matrix(gsvali_mg$now$mat), scale=F, center=T, addEllipses=F, 
    selection_method=NULL, pointsize=3.5, legend_label_size=12,
    label_size=4,show_label="none",colors = genv$Group_col, axis_label_size=9,
    axis_title_label_size=10, group=gsvali_mg$now$INFO$type); gsvali_mg$now$p$plot
  if(F){
    ggsave(rtnP(genv$P2,c("Fig4A_pca",gsvali_mg$now$tp),"png"),
           width=3.5, height=3, bg="white", plot=gsvali_mg$now$p$plot)
  }
  
  p = gsvali_mg$now$p
  for(num in 1:2){
    gsvali_mg$now[[paste('dim_',num,'st_var',sep = "")]] =
      fviz_contrib(p$PCA_res,"var",axes=num)+theme(axis.text.x=element_blank()) +
      ggtitle(paste("Contribution of Genes to PC",num,sep = ""))
    gsvali_mg$now[[paste('dim_',num,'st_ind',sep = "")]] =
      fviz_contrib(p$PCA_res, "ind",axes =num)+
      ggtitle(paste("Contribution of Samples to PC",num,sep = ""))
  }; rm(num, p)
  
  gsvali_mg$now$pca_gene = list()
  for(num in 1:2){
    dim_tmp=gsvali_mg$now[[paste('dim_',num,'st_var',sep="")]]$data %>%
      .[order(.$contrib,decreasing=T),]
    dim_tmp$cum_contrib = 0
    for(k in 1:nrow(dim_tmp)){
      dim_tmp$cum_contrib[k] = sum(dim_tmp$contrib[1:k])}
    cutoff = 100/length(dim_tmp$contrib)
    # PC1에 기여하는 변수들의 기대 기여도는 (C1 * Eig1) / Eig1 임
    # 이때 C1은 각 변수가 PC1에 기여하는 기여도에 해당함.
    # 그러므로 모든 변수들의 PC1에 대한 기여도는 합이 100이 되며,
    # 평균은 sum(C1) / n 이 됨. => 2.040816
    # cutoff = 2.040816
    ak = dim_tmp$name[dim_tmp$contrib >= cutoff] %>% as.character()
    if(24 %in% gsvali_mg$now$INFO$hours_num){
      ak = setdiff(ak, c("HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_MYOGENESIS"))
    } else if(6 %in% gsvali_mg$now$INFO$hours_num){
      ak = setdiff(ak, c("HALLMARK_SPERMATOGENESIS","HALLMARK_HEME_METABOLISM",
                         "HALLMARK_ALLOGRAFT_REJECTION"))}
    gsvali_mg$now$pca_gene[[paste("PC",num,sep="")]]=ak}; rm(dim_tmp,cutoff,ak,k,num)

  # now 에 대한 heatmap 분석
  gsvali_mg$now$cola = HeatmapAnnotation(
    # Group= anno_simple(gsvali_mg$now$INFO$type, col = list(Group=genv$Group_col)$Group,
    #                    pch = c("A", "B", "A", " B", " A"," B"," C"))
    Group = gsvali_mg$now$INFO$type, col=list(Group=genv$Group_col)
    )
              
  # now가 late면 left annotation 적용.
  if(gsvali_mg$now$tp == "24h"){
    gsvali_mg$now$hmCats=list(
      cats = c("proliferation", "immune" ,"proliferation", "pathway","immune", "immune",
               "cellular component", "DNA damage","cellular component", "development", "immune", "proliferation")) %>% 
      append(list(cats_col = setNames(brewer.pal(n=6,name="Set3"), unique(.$cats))))
    gsvali_mg$now$rola=rowAnnotation(
      Category=gsvali_mg$now$hmCats$cats,
      col=list(Category=gsvali_mg$now$hmCats$cats_col),
      annotation_name_side="top", annotation_name_gp=gpar(fontsize=9))
    
    gsvali_mg$now$bota <- HeatmapAnnotation(
      labels = anno_text(c("K2", "K1", "W1", "W2", "I3","I1","I2"),
                         rot = 0, just = "top", gp = gpar(fontsize = 10)),
      which = "column", # 컬럼 라벨에 적용
      height = unit(2, "cm") # 어노테이션 높이 설정
    )
  } else {
    gsvali_mg$now$hmCats = NULL
    gsvali_mg$now$rola = NULL
    gsvali_mg$now$bota = NULL
  }
  
  
  hm=draw(TM_EDA_Heatmap(
    mat=gsvali_mg$now$mat[gsvali_mg$now$pca_gene$PC1,],show_row_names=T,#row_km=2,
    selection_method=NULL,selection_n=NULL,height=7,cluster_row=F,row_title_gp=0,width=5,
    row_names_gp=10,top_annot_df=gsvali_mg$now$cola,column_dend_height=0.5, row_split = c(1,2,1,1,2,2,2,1,2,2,2,1),
    column_names_rot = 0, col_names_gp = 10, 
    # show_col_names = T, column_labels = c("K2", "K1", "W1", "W2", "I3","I1","I2"),
    bottom_annot_df = gsvali_mg$now$bota,
    name="Z score (GSVA)",left_annotation=gsvali_mg$now$rola,do_scale=T, seed =165,
    heatmap_legend_param=list(legend_height=unit(2,"cm"), direction="horizontal")),
    annotation_legend_side="left",heatmap_legend_side="left",merge_legend=T)
  
  if(!B_source_mode){
    png(rtnP(genv$P2,c("Fig4B_",as.character(gsvali_mg$now$tp)),"png"),
        width=10,height=5,units="in",res=1200)
    draw(hm); dev.off()
  }
  
  # 특정 Pathway에 속하는 유전자들의 발현 데이터로 heatmap 분석
  names(gsvali_mg$gmts_filt$H)
  
  # 아래 코드를 T로 한 번 쭉해서 그림 만들고 F로 해서 그림 만들고 해야댐.
  Fig4_Supple_mode = T

  gsvali_mg$now$each_hallmarks_res = sapply(
    c("C:HALLMARK_E2F_TARGETS", "C:HALLMARK_MYC_TARGETS_V1",
      #"C:HALLMARK_MYC_TARGETS_V2",
      "D:HALLMARK_INTERFERON_ALPHA_RESPONSE", "D:HALLMARK_INTERFERON_GAMMA_RESPONSE"
      #"D:HALLMARK_COAGULATION"#,"D:HALLMARK_COMPLEMENT"
      ),
    
    function(gsName){
      categ = strsplit(gsName, ":")[[1]][1]
      gsName = strsplit(gsName, ":")[[1]][2]
      if(categ == "C"){ hmH = 4.5 } else { hmH = 4.5 }
      gss = gsvali_mg$gmts_filt$H[[gsName]]
      
      if((gsName=="HALLMARK_E2F_TARGETS") |
         (gsName=="HALLMARK_INTERFERON_ALPHA_RESPONSE")){
        if(Fig4_Supple_mode){
          top_annot_df = HeatmapAnnotation(
            Group = gsvali_mg$now$INFO$type, col=list(Group=genv$Group_col)
          )
        } else {
          top_annot_df = HeatmapAnnotation(
            Group = gsvali_mg$now$INFO$type, col=list(Group=genv$Group_col), show_legend = F, show_annotation_name = F
          )
        }
        bot_annot = NULL
        # column_labels = c("A", "B", "A", " B", " A"," B"," C")
        # show_col_names = T
      } else {
        top_annot_df = "no"
        bot_annot = gsvali_mg$now$bota
        
        # column_labels = NULL
        # show_col_names = F
      }
      column_labels = c("A", "B", "A", " B", " A"," B"," C")
      show_col_names = T
      
      gsName = gsub("HALLMARK_","",gsName)
      mat = TM_PRC_GE_probe_to_entrez_by_IQR(
        merge_li$GE_combat[intersect( rownames(merge_li$GE_combat), unique(gss)),] %>%
          select(all_of(rownames(gsvali_mg$now$INFO))),
        merge_li$GENEANNO, "entrez", "symbol")
      
      sapply(strsplit(colnames(mat), "_"), function(x) x[3])
      
      
      if(Fig4_Supple_mode){
        cluster_col = T
        bottom_annot_df = NULL
        fig3b_colorder = NULL
      } else {
        bottom_annot_df = bot_annot
        cluster_col = F
        fig3b_colorder = hm@ht_list[["Z score (GSVA)"]]@column_order
      }
      hm = TM_EDA_Heatmap(
        mat=mat, name="Z score",
        row_names_gp = 8, show_row_names=F,
        column_order = fig3b_colorder,
        
        bottom_annot_df = bottom_annot_df, cluster_col = cluster_col,
        
        top_annot_df=top_annot_df, selection_method = NULL,
        do_scale=T, width=5, column_dend_height=0.5, height=hmH,
        heatmap_legend_param=list(legend_height=unit(3,"cm")),
        row_title = gsub("Interferon","IFN",
                         gsub("Alpha","α",
                              gsub("Gamma","γ",
                                   str_to_title(gsub("_", " ",gsName))))),
        row_title_gp=10, row_title_rot = 90)
      
      pca=TM_EDA_PCA(mat, scale=T, center=T,
                     addEllipses=F,label_size=4, pointsize=3.5, legend_label_size=12,
                     show_label="none", colors = genv$Group_col,
                     ellipse.type = "norm", # norm confidence
                     axis_label_size=9, axis_title_label_size=10,
                     group=gsvali_mg$now$INFO$type)
      
      return(list(mat = mat, hm=hm, pca=pca))}, simplify = F)
  
    if(F){
      if(Fig4_Supple_mode){
        png(rtnP(genv$P2,"SFig4_A","png"),width=4,height=5,units="in",res=1200)
        draw(Reduce("%v%", sapply(gsvali_mg$now$each_hallmarks_res, function(x) x$hm)[1:2]))
        dev.off()
        
        png(rtnP(genv$P2,"SFig4_B","png"),width=4,height=5,units="in",res=1200)
        draw(Reduce("%v%", sapply(gsvali_mg$now$each_hallmarks_res, function(x) x$hm)[3:4]))
        dev.off()
      } else {
        png(rtnP(genv$P2,"Fig4C_hm","png"),width=4,height=5,units="in",res=1200)
        draw(Reduce("%v%", sapply(gsvali_mg$now$each_hallmarks_res, function(x) x$hm)[1:2]))
        dev.off()
        
        png(rtnP(genv$P2,"Fig4D_hm","png"),width=4,height=5,units="in",res=1200)
        draw(Reduce("%v%", sapply(gsvali_mg$now$each_hallmarks_res, function(x) x$hm)[3:4]))
        dev.off()
      }
  }
}


# 특정 PW에서 ANova와 사후검정
anova_3grp = function(expr_mat, info){
  tk=as.data.frame(t(expr_mat)) %>% 
    merge(.,info, by="row.names") %>% column_to_rownames("Row.names")
  
  aov_res = sapply(colnames(tk)[1:(ncol(tk)-ncol(info))], function(gn){
    ANOVA_1 <- aov(tk[,gn] ~ as.factor(tk$type))
  }, simplify = F)
  
  aov_pval = sapply(aov_res, function(x) {summary(x)[[1]]$`Pr(>F)`[1]})
  aov_pval.adj = p.adjust(aov_pval, "BH")
  
  signif_aov_res = aov_res[aov_pval.adj < 0.05]
  print(table(aov_pval.adj < 0.05))
  final = sapply(signif_aov_res, function(x) {
    tk_hsd = TukeyHSD(x)
    tmp = tk_hsd$`as.factor(tk$type)`
    if(
      tmp["KO-IM","p adj"] > 0.05 &
      max(tmp["WT-IM","p adj"], tmp["WT-KO","p adj"]) < 0.05 &
      tmp["WT-IM","diff"]*tmp["WT-KO","diff"] > 0
    ){return(tk_hsd)} else {return("Sig but, no interest")}
  }, simplify = F) %>% .[. != "Sig but, no interest"]
  return(final)}

hm@ht_list[["Z score (GSVA)"]]@column_order
{
  # 특정 Geneset의 ANOVA 통과한 유전자들로 HM
  tgtime = 4:1
  mat = rbindlist(sapply(gsvali_mg$now$each_hallmarks_res, function(x){
    x$mat %>% rownames_to_column("Gene")}, simplify=F)[ tgtime ]) %>% 
    .[!duplicated(.),] %>% `rownames<-`(NULL) %>% column_to_rownames("Gene")
  
  dk = anova_3grp(mat, merge_li$INFO)
  print(paste(names(dk), collapse = ", "))
  print(length(names(dk)))
  hm = TM_EDA_Heatmap(
    mat=mat[names(dk),], row_names_gp = 9, show_row_names =T, cluster_col = F, 
    cluster_row = F, column_order = c(3,4,6,7,5,2,1), top_annot_df=gsvali_mg$now$cola,# row_order = rev(names(dk)),
    selection_method = NULL,  name = "Z score",
    row_split = c(rep("IFN α/γ response", 3), rep("E2f/Myc target", 13)),
    do_scale=T, width=5, column_dend_height=0.5, height=9.5, 
    bottom_annot_df = gsvali_mg$now$bota, row_title_gp=9, row_title_rot = 90); hm
  
  if(F){
    png(rtnP(genv$P2, c("Fig4E_","24h"), "png"),
        width=4.5,height=5,units="in",res=1200)
    draw(hm); dev.off(); 
    rm(mat, hm, dk)

  } else {
    draw(hm); rm(mat, hm, dk)
  }
}

if(!B_source_mode){
  # pathway별 GSVA score 변화량
  gsvali_mg$each_HM_line = list(
    HMGS = c(
      #development
      #"HALLMARK_ANGIOGENESIS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
      ##proliferation
      #"HALLMARK_E2F_TARGETS","HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2",
      #"HALLMARK_P53_PATHWAY","HALLMARK_G2M_CHECKPOINT",
      ##dna_damage
      #"HALLMARK_DNA_REPAIR"
      ##immune
      #"HALLMARK_COMPLEMENT",_"HALLMARK_COAGULATION",
      #"HALLMARK_INTERFERON_ALPHA_RESPONSE",_"HALLMARK_INTERFERON_GAMMA_RESPONSE"
      #"HALLMARK_IL6_JAK_STAT3_SIGNALING",_"HALLMARK_INFLAMMATORY_RESPONSE",
      ##Signaling
      "HALLMARK_IL2_STAT5_SIGNALING",'HALLMARK_WNT_BETA_CATENIN_SIGNALING',
      "HALLMARK_HEDGEHOG_SIGNALING",'HALLMARK_MTORC1_SIGNALING',
      "HALLMARK_ESTROGEN_RESPONSE_EARLY","HALLMARK_ESTROGEN_RESPONSE_LATE",
      #"HALLMARK_KRAS_SIGNALING_DN","HALLMARK_KRAS_SIGNALING_UP",
      'HALLMARK_PI3K_AKT_MTOR_SIGNALING',#"HALLMARK_NOTCH_SIGNALING",
      "HALLMARK_TGF_BETA_SIGNALING","HALLMARK_TNFA_SIGNALING_VIA_NFKB"
      ##Cellular_components
      #"HALLMARK_APICAL_JUNCTION",_"HALLMARK_APICAL_SURFACE"
    )
  )
  
  gsvali_mg$each_HM_line$lines = sapply(zipup(
    gsvali_mg$each_HM_line$HMGS,
    c("IL2 STAT5 signaling", "Wnt/β-catenin signalling", "Hedgehog signaling",
      "mTORC1 signaling","E2 response early","E2 response late","PI3K/AKT/mTOR signaling",
      "TGF beta signaling", "TNFA signaling via NFkB")),
    function(x){
      gs=x[1]
      df=TM_PA_custom_geneset_plot(gsvali_mg$gsva_res$H, gsvali_mg$INFO$condition,
                                   rownames(gsvali_mg$INFO), gs) %>% 
        merge(., unique(gsvali_mg$INFO %>% rename(group:=condition)))
      
      ggline(df,x="hour_g",y="Score",color="type",palette=genv$Group_col,add="mean_se") +
        labs(x=NULL,y=x[2]) + stat_compare_means(
          aes(group = type), label = "p.signif", method="anova",size=4, bracket.size=0,
          label.y=max(abs(df$Score))+.05)+lims(y=c(min(df$Score),
                                                   (max(abs(df$Score))+.1)))+
        theme(axis.text.x=element_text(size=11), axis.text.y=element_text(size=11),
              axis.title.x=element_text(size=13),axis.title.y=element_text(size=10),
              legend.position = "none", aspect.ratio=0.6)}, simplify = F)
  
  annotate_figure(ggarrange(plotlist = gsvali_mg$each_HM_line$lines, nrow=3, ncol=3),
                  bottom = textGrob("E2-elapsed time", gp = gpar(cex = 1.3)))
  ggsave(rtnP(genv$P2,"01_Fig2","png"),width=9,height=6,bg="white")
}

