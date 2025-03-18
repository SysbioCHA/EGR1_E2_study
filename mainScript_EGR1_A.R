# main script
rm(list=ls());gc()

# Sys.getenv("GITHUB_PAT")
# Sys.unsetenv("GITHUB_PAT")
# Sys.setenv(GITHUB_PAT = "ghp...")

if (!require("immunedeconv")){
  remotes::install_github('omnideconv/immunedeconv')
}

pacman::p_load(dplyr, data.table, RColorBrewer, pathview, dendextend, ggplot2, forcats,
               stringr, factoextra, STRINGdb, rrvgo, patchwork, msigdbr,VennDiagram,
               GSVA, cmapR, corrplot, colorspace,genefilter, immunedeconv,ggpubr,
               future, future.apply, glue, rstatix, magrittr, colorspace)

library(openxlsx)
source("./TM_sub_data/TM_sub_data_EGR1.R")
source("./TM_sub_data/TM_sub_data_EGR1_immature_mouse.R")
sapply(list.files("./hoon_bi_functions", recursive = T,full.names = T), function(x){
  source(x); print(x); return(T)})

select = dplyr::select

par(family = "sans")

library(ggplot2)

# 기본 테마를 설정
theme_set(
  theme_minimal(base_family = "sans") +
    theme(
      text = element_text(family = "sans"),
      axis.text = element_text(family = "sans"),
      legend.text = element_text(family = "sans"),
      legend.title = element_text(family = "sans")
    )
)

library(grid)
default_gpar <- gpar(fontfamily = "sans")

genv = list(P1 = "./sub_data/EGR1_study", P2 = "./res_plts/", # P2 = "/research_data/egr1_res_plts/0828/",
            P3 = "./egr1_files",
            ko_col = "#F97B22", wt_col = "#7C9070", im_col = "#6C9BCF", ao_col='#adada9',
            Group_col=c("KO"="#F97B22", "WT"="#7C9070", "IM"="#6C9BCF", "AO"="#adada9"),
            UpDwColFn=colorRamp2(c(-2,0,2), c("#2166AC", "#F2F2F2", "#B2182B")),
            h_col_fun=function(x){sequential_hcl(x, palette='Blues 2',rev = T)}) %>% 
  append(list(entrez_annotation=readRDS(rtnP(c(.$P1, "mmuhsa_map"),
                                             "HOM_MouseHumanSequence", "rds"))))

###########################################################
#      Processing & Heatmap & PCA
###########################################################
{
  # The GE_INFO folder must contain the EGR1 dataset named "EGR1_KO_microarray.csv."
  
  # Processing 1
  prc1_resli = TM_EGR1_processing_1(rtnP(c(genv$P1,"GE_INFO"),"EGR1_KO_microarray","csv"))
  sapply(prc1_resli$INFO, function(x) return(length(unique(x))))
  
  # Heatmap
  prc1_resli$hm=TM_EDA_Heatmap(mat=prc1_resli$GE,selection_method="IQR",
                                 selection_n=300, top_annot_df=prc1_resli$INFO %>%
                                   dplyr::select(-c("condition"))); prc1_resli$hm
  rm(prc1_resli)
}

# Processing 2
prc2_resli = TM_EGR1_processing_2(rtnP(c(genv$P1,"GE_INFO"),"EGR1_KO_microarray","csv"))
prc2_resli$INFO$hours_num = as.numeric(gsub("h","",prc2_resli$INFO$hours))

# GENEANNO symbol
prc2_resli$GENEANNO_sym = TM_PA_utils_mmuhsa_converter(
  input_v=unique(prc2_resli$GENEANNO$entrez),input_type="entrez",input_species="mouse",
  output_type = "symbol", output_species = "mouse", study_title = "EGR1_study") %>% 
  select(mouse_EntrezGene.ID,mouse_Symbol)%>% 
  rename("entrez":="mouse_EntrezGene.ID","symbol":="mouse_Symbol")

# Update GE matrix ( row == entrez ), ( replace by mean max )
prc2_resli$GE_entrez = TM_PRC_GE_probe_to_entrez_by_IQR(
  prc2_resli$GE, prc2_resli$GENEANNO, "probeID", "entrez")

# remove EGR1
prc2_resli$GE_entrez2 = prc2_resli$GE_entrez[
  -which(rownames(prc2_resli$GE_entrez) == "13653"),]

# simplify GE by hour
prc2_resli$GE_entrez_by_hours = sapply(unique(prc2_resli$INFO$condition), function(x){
  mg_tgs = rownames(prc2_resli$INFO %>% filter(condition == x))
  res  = rowMeans(prc2_resli$GE_entrez %>% dplyr::select(all_of(mg_tgs)))
  return(res)
})

# simplify INFO by hour
prc2_resli$INFO_by_hours = prc2_resli$INFO[!duplicated(prc2_resli$INFO$condition),]

# IQR filter
prc2_resli$GE_entrez3 = prc2_resli$GE_entrez2 %>%
  mutate(iqrs = apply(., 1, function(y) return(IQR(y))) ) %>%
  filter(iqrs >= quantile(iqrs, prob= c(0.5))) %>% dplyr::select(-iqrs)

# heatmap top annotation
prc2_resli$Col_fun = setNames(genv$h_col_fun(4),unique(prc2_resli$INFO$hours))

prc2_resli$cola1 = HeatmapAnnotation(
  annotation_name_side="right", annotation_name_align=T,
  `E2 treatment` = prc2_resli$INFO$hours, Egr1 = prc2_resli$INFO$type,
  col=list(`E2 treatment`= prc2_resli$Col_fun, Egr1=genv$Group_col),
  annotation_legend_param=list(`E2 treatment`=list(at=unique(prc2_resli$INFO$hours))))

# hm_top_anno$ha1_by_hours
prc2_resli$cola1_by_hours = HeatmapAnnotation(
  annotation_name_side="right",annotation_name_align=T,
  `E2 treatment`=prc2_resli$INFO_by_hours$hours, Egr1=prc2_resli$INFO_by_hours$type,
  col=list(`E2 treatment`= prc2_resli$Col_fun, Egr1=genv$Group_col),
  annotation_legend_param=list(`E2 treatment`=list(at=unique(prc2_resli$INFO$hours))))

# hm_top_anno$ha1_uniq_hrs
prc2_resli$cola1_uniq_hrs = HeatmapAnnotation(
  annotation_name_side="right",annotation_name_align=T,
  annotation_name_gp=gpar(fontsize=8), simple_anno_size=unit(.3,"cm"),
  `E2 treatment` = unique(prc2_resli$INFO_by_hours$hours),
  col=list(`E2 treatment`=prc2_resli$Col_fun),
  annotation_legend_param=list(`E2 treatment`=list(at=unique(prc2_resli$INFO$hours))))

# hm_top_anno$ha1_uniq_hrs_early 0 3 6
prc2_resli$cola1_uniq_hrs_early = HeatmapAnnotation(
  annotation_name_side="right",annotation_name_align=T, show_legend = T,
  annotation_name_gp=gpar(fontsize=8), simple_anno_size=unit(.3,"cm"),
  `E2 treatment`=unique(prc2_resli$INFO_by_hours$hours)[1:3],
  col=list(`E2 treatment`=prc2_resli$Col_fun[1:3]),
  annotation_legend_param=list(
    `E2 treatment`=list(at=unique(prc2_resli$INFO$hours)[1:3])))

# hm_top_anno$ha1_uniq_hrs_early 0 3 6
prc2_resli$cola1_uniq_hrs_early_nozero = HeatmapAnnotation(
  annotation_name_side="right",annotation_name_align=T, show_legend = T,
  annotation_name_gp=gpar(fontsize=8), simple_anno_size=unit(.3,"cm"),
  `E2 treatment`=unique(prc2_resli$INFO_by_hours$hours)[2:3],
  col=list(`E2 treatment`=prc2_resli$Col_fun[2:3]),
  annotation_legend_param=list(
    `E2 treatment`=list(at=unique(prc2_resli$INFO$hours)[2:3])))

# hm_top_anno$ha1_uniq_hrs_nozero 3 6 24 
prc2_resli$cola1_uniq_hrs_nozero = HeatmapAnnotation(
  annotation_name_side="right",annotation_name_align=T, show_legend = T,
  annotation_name_gp=gpar(fontsize=8), simple_anno_size=unit(.3,"cm"),
  `E2 treatment`=unique(prc2_resli$INFO_by_hours$hours)[2:4],
  col=list(`E2 treatment`=prc2_resli$Col_fun[2:4]),
  annotation_legend_param=list(
    `E2 treatment`=list(at=unique(prc2_resli$INFO$hours)[2:4])))

{
  png(rtnP(genv$P2,"Fig1B_heatmap", "png"),width=6,height=4,units="in",res=300)
  draw(TM_EDA_Heatmap(mat=prc2_resli$GE_entrez2,selection_method="IQR",
                      name="Z score", column_dend_height=1,selection_n=1000,
                      show_row_dend=F,top_annot_df=prc2_resli$cola1,width=10,height=5))
  dev.off()
}


###########################################################
#      DEG
###########################################################
## DEG set (main). KO vs WT DEG & Stats
DEG_res=TM_DEG_Limma_loop(.05,.6, prc2_resli$INFO %>% filter(hours_num!=0),
                          prc2_resli$GE,prc2_resli$GENEANNO %>% filter(entrez!="13653"),
                          'hours',"type","KO","WT",var_cut=.4,sample_colname=NULL,
                          symAnno = prc2_resli$GENEANNO_sym)
{
  write.xlsx(sapply(DEG_res$limma,function(x) {
    merge(x,prc2_resli$GENEANNO_sym) %>% 
      dplyr::select(all_of(c("symbol", "probeID", "entrez", "logFC", "P.Value", "adj.P.Val"))) %>% 
      rename('log2 Fold change' = "logFC",
             'p-value' = "P.Value",
             'adj.p-value' = 'adj.P.Val') %>% 
      filter(`adj.p-value` < 0.05)},simplify=F),
    rtnP(genv$P2, "STable1_DEG","xlsx"), rowNames=F)
  print(sapply(DEG_res$stats, function(x){ x[["DEG_stats"]] }))
}

# Entrez genes which have no symbol will be removed.

# DEG volcano
{
  DEG_res$volcanos = sapply(names(DEG_res$limma), function(res){
    TM_DEG_volcano(merge(DEG_res$limma[[res]], prc2_resli$GENEANNO_sym),p_size=1,
                   l_size=3.5,x="logFC",y="adj.P.Val",gene_lab="symbol",show_genename=T,
                   title=res,
                   Res4Genes=fread(rtnP(genv$P3,"Res4Genes","txt"),header=F)$V1)+
      rremove("ylab")+rremove("xlab")}, simplify = F)
  
  DEG_res$volc_figure=annotate_figure(ggarrange(
    plotlist=(DEG_res$volcanos %>% append(list(DEG_res$statsSym_bar))),
    ncol=2,nrow=2, legend="none"),
    bottom = textGrob(expression(log[2]~Fold~Change), gp = gpar(cex = 1.3)),
    left=textGrob(expression(-log[10]~adj.p-value),rot=90,gp=gpar(cex=1.3)))
  
  DEG_res$lgd=Legend(at=c(0,1),labels=c("Up in KO","Down in KO"),type="points",pch=19,
                     background="white",legend_gp=gpar(col=genv$UpDwColFn(c(2,-2))))
  
  ggsave(rtnP(genv$P2,"Fig1C_volcano","png"),plot=DEG_res$volc_figure,width=7,height=7)
  png(rtnP(genv$P2,"Fig1C_lg","png"),width=500,height=100,res=300)
  draw(DEG_res$lgd);dev.off()
}


# Parp1 expression
{ ck = prc2_resli$GE_entrez2["11545",] %>% dplyr::select(contains("24h")) %>% 
    pivot_longer(cols=-NULL,names_to="sample",values_to="Expr") %>% 
    merge(., prc2_resli$INFO, by.x="sample", by.y='row.names') %>% arrange(hours_num)
  plt=box_ggplot(ck, "hours", "Expr",fill="type", stat_line=F, legend_hide=F,lwd=.2,
                 x.axis.size=10,y_lab="Parp1 expression",manual_colors=genv$Group_col)+
    geom_segment(aes(x=0.8,xend=1.2,y=11.5,yend=11.5),linewidth=.1)+
    geom_text(aes(x=1,y=11.53),label="***",size=2.5);plot(plt);rm(ck)
  ggsave(rtnP(genv$P2,"Fig3B_parp1","png"),width=2.3,height=2,dpi=300,plot=plt);rm(plt)}

# E2f1 expression
{ ck = prc2_resli$GE_entrez2["13555",] %>% dplyr::select(contains("24h")) %>% 
    pivot_longer(cols=-NULL,names_to="sample",values_to="Expr") %>% 
    merge(., prc2_resli$INFO, by.x="sample", by.y='row.names') %>% arrange(hours_num)
  plt=box_ggplot(ck, "hours", "Expr",fill="type", stat_line=F, legend_hide=F,lwd=.2,
                  x.axis.size=10,y_lab="E2f1 expression",manual_colors=genv$Group_col)+
      geom_segment(aes(x=0.8,xend=1.2,y=11,yend=11),colour="black",linewidth=.5)+
      geom_text(aes(x=1,y=11.03),label="***",size=2.5);plot(plt);rm(ck)
  ggsave(rtnP(genv$P2,"SFigN_e2f1","png"), width=2.5,height=2,dpi=300,plot=plt);rm(plt)}

# PARylation-related processes.
# https://doi.org/10.1016/j.celrep.2020.108176


###########################################################
#      ORA
###########################################################
# ORA of KO vs WT
Args=list(ora=list(org="Mm",keyType="ENTREZID",ont="BP",pAdjustMethod="BH",readable=T,
                   pvalueCutoff=.05,qvalueCutoff=.05,minGSSize=25,maxGSSize=500,
                   sim_threshold=.7),
          fc=.6,adjp=.05) %>% append(list(OFN=paste0(
            "./sub_data/EGR1_study/PA/ORA2_",.$fc,".rds")))

PA_ORA2_res = list()
if(file.exists(Args$OFN)){
  PA_ORA2_res = list(GO=readRDS(Args$OFN))
} else {
  PA_ORA2_res = list(GO=future_lapply(DEG_res$limma, function(x) {
    genesUp = x %>% dplyr::filter(logFC >= Args$fc & adj.P.Val <= Args$adjp) %>%
      dplyr::select(entrez) %>% as.matrix() %>% as.character()
    genesDw = x %>% dplyr::filter(logFC <= -Args$fc & adj.P.Val <= Args$adjp) %>%
      dplyr::select(entrez) %>% as.matrix() %>% as.character()

    resU = do.call(TM_PA_ORA_reduce, c(list(genesUp), Args$ora))
    resD = do.call(TM_PA_ORA_reduce, c(list(genesDw), Args$ora))
    return(list("Up"=resU, "Down"=resD))}))

  # saveRDS(PA_ORA2_res$GO, file = Args$OFN)
}

PA_ORA2_res$GO_unlist = sapply(c("Up", "Down"), function(z){
  res = sapply(names(PA_ORA2_res$GO), function(x){
    tmp = as.data.frame(PA_ORA2_res$GO[[x]][[z]]$GOreduced$Res) %>%
      mutate(Description = paste0(.$ID,": ", .$Description)) %>% mutate(cond=x)
    if(nrow(tmp) > 0){ tmp = tmp %>% slice_min(n = 20, order_by = p.adjust)
    }
    return(tmp) }, simplify = F) %>% bind_rows(.)
  return(res)}, simplify = F)

PA_ORA2_res$rm_tg = sapply(PA_ORA2_res$GO_unlist, function(regul){
  rm_tg = apply(regul, 1, function(x){
    x_tg = str_split(x["geneID"], "/")[[1]]
    rm_tg = apply(regul, 1, function(y){
      y_tg = str_split(y["geneID"],"/")[[1]]
      if(all(x_tg == y_tg) & (paste0(x["Description"], x["cond"]) !=
                                paste0(y["Description"], y["cond"])) &
         (x["cond"] == y["cond"])){ return(x["Description"])
      } else {return(NA)}}, simplify=F) %>% unlist() %>% unique()
    if(all(is.na(rm_tg))){return(NA)} else {return(rm_tg)}})
  return(rm_tg %>% unlist() %>% as.character() %>% .[!is.na(.)])})

# write.csv(PA_ORA2_res$GO_unlist$Up %>% filter(!Description %in%PA_ORA2_res$rm_tg$Up),
#           file=rtnP(c(genv$P3,"PA"),"wtko_u","csv"),row.names=F)
# write.csv(PA_ORA2_res$GO_unlist$Down %>% filter(!Description %in%PA_ORA2_res$rm_tg$Down),
#           file=rtnP(c(genv$P3,"PA"),"wtko_d","csv"),row.names=F)

PA_ORA2_res$GO_unlist$Up = fread(rtnP(c(genv$P3,"PA"), "wtko_uf0828","csv")) %>% 
  filter(ID != 0) %>% mutate(Description = strsplit(Description, ": ") %>% map_chr(.,2))
PA_ORA2_res$GO_unlist$Down = fread(rtnP(c(genv$P3,"PA"), "wtko_df0828","csv")) %>% 
  filter(ID != 0) %>% mutate(Description = strsplit(Description, ": ") %>% map_chr(.,2))



{ # 3h - 6h
  options(scipen=100)
  pdu = TM_PA_ORA_multiple_condition(
    PA_ORA2_res$GO_unlist$Up %>% filter(cond != "24h"),title="up in KO", y_lab_size=12,
    order=names(PA_ORA2_res$GO),lgd_wid = .5, TM_PA_utils_coder(
      PA_ORA2_res$GO_unlist$Up, names(PA_ORA2_res$GO)) %>% arrange(desc(code)))
  pdw = TM_PA_ORA_multiple_condition(
    PA_ORA2_res$GO_unlist$Down %>% filter(cond != "24h"),title="down in KO",y_lab_size=12,
    order=names(PA_ORA2_res$GO),lgd_wid = .5,TM_PA_utils_coder(
      PA_ORA2_res$GO_unlist$Down, names(PA_ORA2_res$GO)) %>% arrange(desc(code)))
  png(file=rtnP(genv$P2,"Fig2AB","png"),width=3300,height=715,res=300);plot(pdu+pdw);dev.off()
  plot(pdu + pdw); rm(pdu, pdw) }

library(scales)
{ # 24h GO terms
  options(scipen=-3)
  A=rbind(PA_ORA2_res$GO_unlist$Up %>% filter(cond=="24h") %>% mutate(cond="Up"), 
          PA_ORA2_res$GO_unlist$Down %>% filter(cond=="24h") %>% mutate(cond="Down")) %>%
    mutate(Regulation = cond)
  
  p1 = TM_PA_ORA_multiple_condition(A, c("Up","Down"), lgd_wid=.5, title="Egr1 KO 24h",
                                    TM_PA_utils_coder(A, c("Up","Down")), pval_scipen_digit_0=T) +
    theme(axis.text.x = element_text(color = "grey20",size=11,hjust=.5,vjust=.5,angle=0),
          axis.text.y = element_text(color = "grey20",size=11,hjust=1,vjust=0))
  png(file=rtnP(genv$P2,"Fig3A_GO","png"),width=1400,height=850,res=300)
  plot(p1);dev.off()
  plot(p1);rm(A, p1)
  options(scipen=100)
}




# Result 2의 CD를 만들기 위헤, mainScript_EGR1_Network.
# 1. DEG들로 네트워크 그리는거.
# 2. 위에서 만든 wtko_d 의 reproductive와 sex diff 로 네트워크


###########################################################
#      GSVA
###########################################################
if(file.exists("./sub_data/EGR1_study/msigdb/gmts.rds")){
  gsvali = list(gmts=readRDS("./sub_data/EGR1_study/msigdb/gmts.rds"))
} else {
  dir.create("./sub_data/EGR1_study/msigdb/")
  gsvali = list(gmts=sapply(c("H", "C2", "C5", "C6"), function(x) {
    msdb = msigdbr(species = "Mus musculus", category = x) %>% 
      dplyr::select(gs_name, gene_symbol, entrez_gene)
    msdb = sapply(unique(msdb$gs_name), function(x){
      rt = msdb %>% filter(gs_name == x) %>% .[["entrez_gene"]]
      return(rt)}); return(msdb)}))
  saveRDS(gmts, file = "./sub_data/EGR1_study/msigdb/gmts.rds")}

###########################################################
#      Custom Gene set Test ( ORA, GSE )
###########################################################
# Unipep - glycosylation
prc2_resli$GENEANNO_mu_ent_2_hm_sym = TM_PA_utils_mmuhsa_converter(
  input_v=unique(prc2_resli$GENEANNO$entrez),
  input_type = "entrez", input_species = "mouse", output_type = "symbol", 
  output_species = "human", study_title = "EGR1_study")

# write.table(unique(prc2_resli$GENEANNO_mu_ent_2_hm_sym$human_Symbol),
#             "./sub_data/EGR1_study/Glycosyl_Paryl/prc2_ent_to_hm_symbol.txt",
#             col.names = F, row.names = F, quote = F)
# 위 txt를 https://unipep.systemsbiology.net/ (Unipep) 에서
# bulk - search.

# PARylation은 위에 PARylation-related process 분석에 ref 달려 있음
tmp="./sub_data/EGR1_study/Glycosyl_Paryl/"; custom_setTest = list(
  adpribos_mouse_entrez  = TM_PA_utils_mmuhsa_converter(
    input_v=fread(paste0(tmp,"Step6_adp_ry_target.csv"))$`Gene name (primary)`,
    input_type = "symbol", input_species = "human", output_type = "entrez", 
    output_species = "mouse", study_title = "EGR1_study"),
  glycosyl_mouse_entrez  = TM_PA_utils_mmuhsa_converter(
    input_v = unique(fread(paste0(tmp,"unipep_bulk_search.xls"))$Symbol),
    # input_v = fread(paste0(tmp,"Step6_EGR1_glycosylation.txt"))$Gene,
    input_type = "symbol", input_species = "human", output_type = "entrez", 
    output_species = "mouse", study_title = "EGR1_study")); rm(tmp)

custom_setTest$gmt = list(
  "ADPribosylation"=custom_setTest$adpribos_mouse_entrez$mouse_EntrezGene.ID,
  "Glycosylation"  =custom_setTest$glycosyl_mouse_entrez$mouse_EntrezGene.ID)

custom_setTest$gsva_res = TM_PA_custom_geneset_analysis(
  as.matrix(prc2_resli$GE_entrez2), method="gsva", custom_setTest$gmt,
  apply(str_split(colnames(prc2_resli$GE),"_",simplify=T),1,
        function(x)paste(x[1:2],collapse="_")))

custom_setTest$gsva_res$boxdf2 = custom_setTest$gsva_res$boxdf %>%
  mutate(Group = str_split(group, "_", simplify = T)[,1]) %>% 
  mutate(Group = factor(Group, levels = unique(Group))) %>% 
  mutate(hours = str_split(group, "_", simplify = T)[,2]) %>% 
  mutate(hours = factor(as.numeric(substr(hours,1,nchar(hours)-1)))) %>% 
  mutate(Pathway=factor(Pathway, levels=c("ADPribosylation","Glycosylation")))


# WT vs KO 는 GSVA 문서에 따라 limma를 활용한 비교
custom_setTest$gsva_res_de = sapply(c("0h","3h","6h","24h"), function(tm){
  res=TM_CUSTOM_geneset_LIMMA(as.data.frame(custom_setTest$gsva_res$mat),var_cut =NULL,test=T,
                   prc2_resli$INFO,sample_colname=NULL,Group="condition",
                   case=paste0("KO_",tm),control=paste0("WT_",tm)) %>% mutate(hours=tm)
  }, simplify = F)

custom_setTest$ttest = custom_setTest$gsva_res_de %>% bind_rows() %>% 
  rename("p":=P.Value) %>% mutate(hours=as.numeric(str_split(hours,"h",simplify=T)[,1]))

# ADPribosylation -  p
custom_setTest$statest_adp = custom_setTest$gsva_res$boxdf2 %>% 
  filter(Pathway=="ADPribosylation") %>% group_by(hours) %>% t_test(Score~Group) %>%
  add_xy_position(x="hours", stack = T) %>% select(-p) %>% merge(
    .,subset(custom_setTest$ttest,probeID=="ADPribosylation",c(p,hours)),by="hours") %>%
  adjust_pvalue(method="BH") %>% arrange(hours) %>% add_significance("p")

custom_setTest$statest_adp$y.position = c(-0.1, 0.241, 0.284, 0.12)

png(file=rtnP(genv$P2,"Fig3C_GSEA_adp","png"), width=1000,height=800, res = 300)
ggline(custom_setTest$gsva_res$boxdf2 %>% filter(Pathway=="ADPribosylation") %>% 
         mutate(hours=paste0(hours,"h")), x="hours", y="Score", color="Group", 
       palette=genv$Group_col, add = "mean_sd") + 
  stat_pvalue_manual(
    custom_setTest$statest_adp,label="p.signif",remove.bracket = T,
                     label.size=4) +
  labs(x = "E2-elapsed time",y="PARylation score") +
  theme(axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),
        axis.title.x=element_text(size=13), axis.title.y=element_text(size=13)) + 
  guides(color=guide_legend(title="group"))
dev.off()

# Glycosylation
custom_setTest$statest_gly=custom_setTest$gsva_res$boxdf2 %>%
  filter(Pathway=="Glycosylation") %>% group_by(hours) %>% t_test(Score~Group) %>%
  add_xy_position(x="hours", stack = T) %>% select(-p) %>% merge(
    ., subset(custom_setTest$ttest,probeID=="Glycosylation",c(p,hours)),by="hours") %>% 
  adjust_pvalue(method="BH") %>% arrange(hours) %>% add_significance("p")

custom_setTest$statest_gly$y.position = c(0.17, 0.125, 0.05, 0.01)

png(file=rtnP(genv$P2,"Fig3D_GSEA_gly","png"),width=1000,height=800,res=300)
ggline(custom_setTest$gsva_res$boxdf2 %>% filter(Pathway=="Glycosylation") %>% 
         mutate(hours=paste0(hours,"h")), x="hours", y="Score", color="Group", 
       palette=genv$Group_col, add = "mean_sd")+ 
  stat_pvalue_manual(custom_setTest$statest_gly,label="p.signif",remove.bracket = T,
                     label.size=4, vjust=.7)+
  labs(x = "E2-elapsed time",y="Glycosylation score") +
  theme(axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),
        axis.title.x=element_text(size=13),axis.title.y=element_text(size=13)) + 
  guides(color=guide_legend(title="group"))
dev.off()

###########################################################
#      Immune deconvolute
###########################################################
library(xCell) # devtools::install_github('dviraran/xCell')

egr1_immune = list()
egr1_immune$GENEANNO = TM_PA_utils_mmuhsa_converter(
  input_v = rownames(prc2_resli$GE_entrez2), study_title = "EGR1_study",
  input_type="entrez",input_species="mouse",output_type="symbol",
  output_species ="human") %>% dplyr::select(mouse_EntrezGene.ID, human_Symbol) %>% 
  rename(entrez:=mouse_EntrezGene.ID, symbol:=human_Symbol)

egr1_immune$GE = TM_PRC_GE_probe_to_entrez_by_IQR(
  prc2_resli$GE_entrez2, egr1_immune$GENEANNO, "entrez", "symbol")
egr1_immune$GE_unlog = apply(egr1_immune$GE, 2, function(x) return(x))

# 성숙자궁
egr1_immune$res =deconvolute(egr1_immune$GE, "xcell",tumor=F,arrays=T)
egr1_immune$prc1=data.frame(egr1_immune$res[,-1], row.names = egr1_immune$res$cell_type)

egr1_immune$xCellLibRes = xCell::xCellAnalysis(egr1_immune$GE, rnaseq=F) %>% 
  as.data.frame() %>% rownames_to_column("cell_type")

egr1_immune$res %>% filter(cell_type == "Endothelial cell") %>% .$KO_0h_A

egr1_immune$prc2 = egr1_immune$res %>% gather(sample, fraction, -cell_type) %>% 
  inner_join(.,prc2_resli$INFO %>% tibble::rownames_to_column('sample')) %>% arrange(hours) %>%
  mutate(cell_type=factor(cell_type, level=rev(unique(cell_type)))) %>%
  mutate(sample=factor(sample,level=rownames(prc2_resli$INFO)[order(
    prc2_resli$INFO$hours_num)])) %>%
  mutate(condition=factor(condition, level=unique(prc2_resli$INFO$condition[order(
    prc2_resli$INFO$hours_num)]))) %>%
  mutate(hours = factor(hours, levels = c("0h", "3h", "6h", "24h")))

egr1_immune$deco_plt= egr1_immune$res %>% gather(sample, fraction, -cell_type) %>%
  ggplot(aes(x=sample,y=fraction,color=cell_type))+geom_point(size=4) +
  facet_wrap(~cell_type, scales = "free_x", ncol = 3) + coord_flip() + theme_bw() +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),legend.position="none")

# DE pw ( para non-para )
egr1_immune$DE_immunes = sapply(c("0h","3h","6h","24h"), function(hour) {
  tmp_1=TM_DE_pws(egr1_immune$prc1, prc2_resli$INFO, Group = "condition",parametric=T,
                  case=paste0('KO_',hour), control = paste0('WT_', hour)) %>% 
    filter(p.value <= 0.05)
  return(tmp_1)}, simplify = F);egr1_immune$DE_immunes

# immune score
egr1_immune$stat_test_Immune = egr1_immune$prc2 %>% filter(grepl("immune",cell_type)) %>%
  mutate(cell_type=str_to_title(cell_type)) %>% dplyr::rename(Group=type) %>% 
  mutate(hours=factor(hours_num,levels=sort(unique(hours_num)))) %>% group_by(hours)%>% 
  t_test(fraction~Group) %>% 
  add_xy_position(x="hours",stack = T) %>% adjust_pvalue(method = "BH") %>%
  arrange(hours) %>% add_significance("p") %>% mutate(y.position = y.position+0.05) %>% 
  mutate(p = ifelse(p > 0.1, "ns", round(p,3)));

egr1_immune$stat_test_Immune$y.position = c(0.36, 0.39, 0.37, 0.2)+0.005

png(rtnP(genv$P2,"Fig3E_Xcell_immune","png"),width=1000,height=800,res=300)
ggline(egr1_immune$prc2 %>% filter(grepl("immune",cell_type)) %>%
         mutate(cell_type=str_to_title(cell_type)) %>% dplyr::rename(Group=type),
       x="hours",y="fraction",color="Group", palette=genv$Group_col, add = "mean_sd")+
  stat_pvalue_manual(egr1_immune$stat_test_Immune, label="p",remove.bracket = T,
                     label.size=4, vjust = .7) +
  labs(x = "E2-elapsed time",y="Immune score") + 
  theme(axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),
        axis.title.x=element_text(size=13),
        axis.title.y=element_text(size=13)) + 
  guides(color=guide_legend(title="group"))
dev.off()

# stroma score
egr1_immune$stat_test_Stroma = egr1_immune$prc2 %>% filter(grepl("stroma",cell_type)) %>%
  mutate(cell_type=str_to_title(cell_type)) %>% dplyr::rename(Group=type) %>% 
  mutate(hours=factor(hours_num,levels=sort(unique(hours_num)))) %>% group_by(hours)%>% 
  t_test(fraction~Group) %>% 
  add_xy_position(x="hours",stack = T) %>% adjust_pvalue(method = "BH") %>%
  arrange(hours) %>% add_significance("p") %>% mutate(y.position = y.position+0.01)

egr1_immune$stat_test_Stroma$y.position = c(0.085, 0.06, 0.08, 0.035)+0.005

png(rtnP(genv$P2,"Fig3F_Xcell_stroma","png"),width=1000,height=800,res=300)
ggline(egr1_immune$prc2 %>% filter(grepl("stroma",cell_type)) %>%
         mutate(cell_type=str_to_title(cell_type)) %>% dplyr::rename(Group=type),
       x="hours",y="fraction",color="Group", palette=genv$Group_col, add = "mean_sd")+
  stat_pvalue_manual(egr1_immune$stat_test_Stroma,label="p.signif",remove.bracket = T,
                     label.size=4, vjust = .7) +
  labs(x = "E2-elapsed time",y="Stroma score") + 
  theme(axis.text.x=element_text(size=11), axis.text.y=element_text(size=11),
        axis.title.x=element_text(size=13), axis.title.y=element_text(size=13)) + 
  guides(color=guide_legend(title="group"))
dev.off()
