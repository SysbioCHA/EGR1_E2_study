# pacman::p_load(GenomicFeatures, BSgenome.Mmusculus.UCSC.mm10)
source("./hoon_bi_functions/GM/GM_G4_analysis_funcs.R")
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)

"Source Data Setting"
genome <- BSgenome.Mmusculus.UCSC.mm10

"UCSC mm10 Genome을 불러와서, R에서 쓸 수 있는 TxDb로 변환."
iDB2 = GenomicFeatures::makeTxDbFromUCSC(genome = "mm10", tablename="refGene")# need

seqnames <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
              "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
              "chr18","chr19","chr20","chr21","chr22","chrX","chrY")# need

"iDB2를 통해서, 각 Gene에 대한 Genome상의 position 정보."
"iDB2의 각 gene의 chromsomal num, 위치, gene ID => GRanges class"
txdb_idb = GenomicFeatures::genes(iDB2, filter = list(tx_chrom = seqnames))# need

"Extract Sequence of All Entrez Genes" # need
all_entrez_genes = rownames(prc2_resli$GE_entrez)

if(F){
  all_entrez_genes_tss = sapply(all_entrez_genes, function(gene){
    if(gene %in% txdb_idb$gene_id){
      gene_tssup1000 = GenomicFeatures::extractUpstreamSeqs(
        genome, txdb_idb[gene], width = 1000)  # promoter쪽
      
      gene_tssdw500 = GenomicFeatures::extractUpstreamSeqs(
        genome, txdb_idb[gene], width = -500)  # promoter쪽
      
      paste_string = paste0( toString(gene_tssup1000),
                             toString(gene_tssdw500))
      return(paste_string)
    } else {
      return(NA)
    }
  },simplify = F) %>% .[!is.na(.)]
  saveRDS(all_entrez_genes_tss,
          file = "./sub_data/EGR1_study/G4/all_entrez_genes_tss_Up1000_Dw500.rds")
  
} else {
  all_entrez_genes_tss = readRDS(
    "./sub_data/EGR1_study/G4/all_entrez_genes_tss_Up1000_Dw500.rds")
  
}

all_gene_seqs = list2DF(all_entrez_genes_tss) %>% t(.) %>% as.data.frame() %>%
  mutate(EntrezID = rownames(.))
colnames(all_gene_seqs)[1] = "Sequence"

pacman::p_load(data.table, dplyr, hrbrthemes, ggplot2, Biostrings)

###########################################################
#      G4 work!
###########################################################
col_labels = c("seqnames", "start",  "end", "width", "strand", "score" ,"sequence")

if(F){
  all_gene_g4res = apply(all_gene_seqs, 1, function(gene_seq){
    seq = gene_seq["Sequence"]
    geneid = gene_seq["EntrezID"]
    
    # >>>>>> 전체 1500bp
    print(geneid)
    dna <- DNAStringSet(c(seq))
    names(dna)[1] = geneid
    
    # G4hunt cutoff 1.2, 1.1, 1.0;  gene 하나에 대해서는 i=1, 그렇지 않으면
    toto=G4hunt(i=1, k=25, hl=1.5, gen=dna, masked=5)
    tryCatch(expr = {
      error_occured = TRUE
      titi12 = G4huntrefined(toto,gen=dna,i=1)
      error_occured = FALSE
    },
    error = function(e){
      print("Error")
    },
    warning = function(e) print("Warning"),
    finally = {
      if(error_occured){
        titi12 <- data.frame(geneid,0,0,0,0,0,"Error") %>% setNames(col_labels)
      }
    })
    
    # NULL 처리
    if(is.null(titi12)){
      titi12 = data.frame(geneid,0,0,0,0,0,"Emty") %>% setNames(col_labels)
    }
    return(as.data.frame(titi12))
  }) %>% rbindlist(.)
  
  # Gtrack num 추가
  all_gene_g4res$Gtrack_num = apply(all_gene_g4res,1,function(x){
    strand = x["strand"]
    seq = x["sequence"]
    
    if(strand == "+"){
      seq_g4t = gregexpr("(G){3,}",seq)
      if(seq_g4t[[1]][1] == -1){ return(0) }
      
    } else if(strand == "-"){
      seq_g4t = gregexpr("(C){3,}", seq)
      if(seq_g4t[[1]][1] == -1){ return(0) }
      
    } else {
      return(0)
      
    }
    return(length(seq_g4t[[1]]))
  })
  saveRDS(all_gene_g4res, "./sub_data/EGR1_study/G4/all_gene_g4res.rds")
  
} else {
  all_gene_g4res = readRDS("./sub_data/EGR1_study/G4/all_gene_g4res.rds")
  
}

# "empty에 대해 n => 0 으로 바꿔주는 작업 필요
res_df_all = all_gene_g4res %>% dplyr::filter(sequence != "Error") %>% 
  group_by(seqnames) %>% summarise(score = sum(score),
                                   n = ifelse(max(sequence) == "Emty", 0, n()),
                                   Gtrack_num = max(Gtrack_num),
                                   sequence = max(sequence)) %>% ungroup()

res_df_plus = all_gene_g4res %>% dplyr::filter(sequence != "Error") %>%
  dplyr::filter(strand == "+" | strand == 0) %>% group_by(seqnames) %>%
  summarise(score = sum(score), n = ifelse(max(sequence) == "Emty", 0, n()),
            Gtrack_num = max(Gtrack_num),
            sequence = max(sequence)) %>% ungroup()

res_df_minus = all_gene_g4res %>% dplyr::filter(sequence != "Error") %>%
  dplyr::filter(strand == "-"| strand == 0) %>% group_by(seqnames) %>%
  summarise(score = sum(score), n = ifelse(max(sequence) == "Emty", 0, n()),
            Gtrack_num = max(Gtrack_num),
            sequence = max(sequence)) %>% ungroup()

deg_wt3u = DEGs_rawres$WT_3_0 %>% filter(logFC > 1 & adj.P.Val < 0.05) %>% .[["entrez"]]
deg_ko3u = DEGs_rawres$KO_3_0 %>% filter(logFC > 1 & adj.P.Val < 0.05) %>% .[["entrez"]]
deg_wt3u_spe = setdiff(deg_wt3u, deg_ko3u)
deg_ko3u_spe = setdiff(deg_ko3u, deg_wt3u) # x만 갖고 있는 것

df_of_plotting = res_df_all

x = df_of_plotting %>% dplyr::filter(seqnames %in% deg_wt3u_spe)
y = df_of_plotting %>% dplyr::filter(seqnames %in% deg_ko3u_spe)

# Build dataset with different distributions
plot_g4_hist(x, "WT3", y, "KO3", "n",
             length(deg_wt3u_spe), length(deg_ko3u_spe),"PQS number") + 
  plot_g4_dist(x, "WT", y, "KO", "score")
plot_g4_bar(x, "WT", y, "KO", "Gtrack_num", 5) # 한 gene이 G4-track >= 5 한 PQS를 갖는가?



###########################################################
#      create logFC diff file
###########################################################
deg_wt3u
deg_ko3u

the_table = merge(DEGs_rawres$KO_3_0, DEGs_rawres$WT_3_0,
                  by = "entrez", all = T)   # ko가 x임

tt2 = the_table %>% filter( (adj.P.Val.x < 0.05 & abs(logFC.x) > 0.6) |
                              (adj.P.Val.y < 0.05 & abs(logFC.y) > 0.6))

tt2_vector = apply(tt2, 1, function(x){
  pval_x = as.numeric(x["adj.P.Val.x"])
  pval_y = as.numeric(x["adj.P.Val.y"])

  if(pval_x < 0.05 & pval_y < 0.05){
    p = 3
  } else if(pval_x < 0.05) {
    p = 1
  } else if(pval_y < 0.05) {
    p = 2
  } else {
    p = 0
  }

  log_x = as.numeric(x["logFC.x"])
  log_y = as.numeric(x["logFC.y"])

  if(abs(log_x) > 0.6 & abs(log_y) > 0.6){
    l = 3
  } else if(abs(log_x) > 0.6) {
    l = 1
  } else if(abs(log_y) > 0.6) {
    l = 2
  } else {
    l = 0
  }

  return(c(p,l))
}) %>% data.frame(row.names = c("p_value", "logFC")) %>% t()

head(tt2_vector)

tt3 = cbind(tt2, tt2_vector)
tt4 = tt3 %>% mutate(logfc_diff = logFC.x-logFC.y)
# tt4

###########################################################
#      BP_position as freqPQS
###########################################################
group1 = tt4 %>% filter(logfc_diff > 0.6) %>% filter(p_value == 3) %>% .$entrez
group2 = tt4 %>% filter(abs(logfc_diff) < 0.01)  %>% .$entrez
group3 = tt4 %>% filter(logfc_diff < -0.6) %>% filter(p_value == 3) %>%  .$entrez

# 막대그래프용
tt5_g1 = tt4 %>% filter(logfc_diff > 0.6) %>% filter(p_value == 3)
tt5_g2 = tt4 %>% filter(abs(logfc_diff) < 0.01)
tt5_g3 = tt4 %>% filter(logfc_diff < -0.6) %>% filter(p_value == 3)

###############   logFC.X - logFC.y #############
res_df_g1 = all_gene_g4res %>% dplyr::filter(sequence != "Error") %>%
  dplyr::filter(strand == "+" | strand == "-") %>% filter(seqnames %in% group1)
res_df_g2 = all_gene_g4res %>% dplyr::filter(sequence != "Error") %>%
  dplyr::filter(strand == "+" | strand == "-") %>% filter(seqnames %in% group2)
res_df_g3 = all_gene_g4res %>% dplyr::filter(sequence != "Error") %>%
  dplyr::filter(strand == "+" | strand == "-") %>% filter(seqnames %in% group3)

res_df_plus_g1 = all_gene_g4res %>% dplyr::filter(sequence != "Error") %>%
  dplyr::filter(strand == "+") %>% filter(seqnames %in% group1)
res_df_plus_g2 = all_gene_g4res %>% dplyr::filter(sequence != "Error") %>%
  dplyr::filter(strand == "+") %>% filter(seqnames %in% group2)
res_df_plus_g3 = all_gene_g4res %>% dplyr::filter(sequence != "Error") %>%
  dplyr::filter(strand == "+") %>% filter(seqnames %in% group3)

res_df_minus_g1 = all_gene_g4res %>% dplyr::filter(sequence != "Error") %>%
  dplyr::filter(strand == "-") %>% filter(seqnames %in% group1)
res_df_minus_g2 = all_gene_g4res %>% dplyr::filter(sequence != "Error") %>%
  dplyr::filter(strand == "-") %>% filter(seqnames %in% group2)
res_df_minus_g3 = all_gene_g4res %>% dplyr::filter(sequence != "Error") %>%
  dplyr::filter(strand == "-") %>% filter(seqnames %in% group3)


library("rjson")
# Give the input file name to the function.
test1 = pqs_frac_plotting_for_each_genes(res_df_g1,total_len = 2000,
                                         add_stuff = "score", strand_status = "both")
test2 = pqs_frac_plotting_for_each_genes(res_df_g2,total_len = 2000,
                                         add_stuff = "score", strand_status = "both")
test3 = pqs_frac_plotting_for_each_genes(res_df_g3,total_len = 2000,
                                         add_stuff = "score", strand_status = "both")

test1 = pqs_frac_plotting_for_each_genes(res_df_plus_g1,total_len = 2000,
                                         add_stuff = "score")
test2 = pqs_frac_plotting_for_each_genes(res_df_plus_g2,total_len = 2000,
                                         add_stuff = "score")
test3 = pqs_frac_plotting_for_each_genes(res_df_plus_g3,total_len = 2000,
                                         add_stuff = "score")

test4 = pqs_frac_plotting_for_each_genes(res_df_minus_g1,total_len = 2000,
                                         add_stuff = "score")
test5 = pqs_frac_plotting_for_each_genes(res_df_minus_g2,total_len = 2000,
                                         add_stuff = "score")
test6 = pqs_frac_plotting_for_each_genes(res_df_minus_g3,total_len = 2000,
                                         add_stuff = "score")


tk = lapply(test1, function(x){
  return(sum(x$counts))
})
table(tk == 0)
sum(tk == 0)/length(tk)
library("ggpubr")
library(RColorBrewer)
# install.packages("https://cran.r-project.org/src/contrib/Archive/colortools/colortools_0.1.5.tar.gz", repos =NULL, type = "source")
library(colortools)
require(graphics)

# Freq
value = "counts"
yl = 70
leg_y = 70
# Percent
# value = "counts_percent"
# yl = 0.2
# leg_y =0.2

"red, orange, blue"
df_4_plot = test1
my_palette = sequential("blue", alpha = 1, percentage = 100/length(df_4_plot),plot = F)[-1]
for(n in c(1:length(names(df_4_plot)))){
  k = names(df_4_plot)[n]
  if(k %in% c("11504", "16590")){
    print("l")
  }else{
    next
  }
  k_genename = tt5 %>% filter(entrez == k) %>% .$To
  k_df = df_4_plot[[k]]
  if(sum(k_df$counts) == 0){
    next
  }
  # if(k == names(df_4_plot)[1]){
  #   plot(smooth.spline(x=k_df[,"relative_bp"],k_df[,value],df=1000),lwd=2,type='l',xlab=paste('Dist to TSS (bp)',sep=''),ylab="Freq of G4FS",col=my_palette[1],cex.axis=1.4,cex.lab=1.4,font=2,font.lab=2,ylim = c(-3,3))
  # } else {
  #   lines(smooth.spline(x=k_df[,"relative_bp"],y =k_df[,value], df=1000),lwd=2,col=my_palette[n])
  # }
  # if(k == names(df_4_plot)[1]){
  plot(x=k_df[,"relative_bp"],k_df[,value],lwd=2,type='l',xlab=paste('Dist to TSS (bp)',sep=''),ylab="Freq of G4FS",
       # col=my_palette[1],
       col = "blue",main = k_genename,
       cex.axis=1.4,cex.lab=1.4,font=2,font.lab=2,ylim = c(-3,3), xlim = c(-2000, 0))
  print("plot complete")
  # }
  # else {
  #   lines(x=k_df[,"relative_bp"],y =k_df[,value],lwd=2,type = 'l',col=my_palette[n])
  # }
}

{df_4_plot_g1 = data.frame(
  relative_bp = test1[[1]]$relative_bp,
  counts = lapply(test1, function(x) return(x$counts)) %>% list2DF(.) %>% as.data.frame() %>% rowSums()
)
  df_4_plot_g2 = data.frame(
    relative_bp = test2[[1]]$relative_bp,
    counts = lapply(test2, function(x) return(x$counts)) %>% list2DF(.) %>% as.data.frame() %>% rowSums()
  )
  df_4_plot_g3 = data.frame(
    relative_bp = test3[[1]]$relative_bp,
    counts = lapply(test3, function(x) return(x$counts)) %>% list2DF(.) %>% as.data.frame() %>% rowSums()
  )}

# Version 1    smooth.spline
plot(x=df_4_plot_g2[,"relative_bp"],df_4_plot_g2[,value],lwd=1,type='l',
     xlab=paste('Dist to TSS (bp)',sep=''),ylab="Freq of G4FS",col='black',
     ylim = c(-15,15),xlim = c(-1500,0))
abline(h=0, lwd = 0.5)
legend("topleft", # 60 or 0.19
       c("Group1","Group3"),
       col = c("red","black"),lty=1,cex=0.8)
lines(x=df_4_plot_g1[,"relative_bp"],y =df_4_plot_g1[,value],lwd=1,col='red')
lines(x=df_4_plot_g3[,"relative_bp"],y =df_4_plot_g3[,value],lwd=1,col='blue')


lines(smooth.spline(x=df_4_plot_g2[,"relative_bp"],y =df_4_plot_g2[,value],df=100),lwd=2,col='orange')
lines(smooth.spline(x=df_4_plot_g3[,"relative_bp"],y =df_4_plot_g3[,value],df=100),lwd=2,col='blue')
lines(smooth.spline(x=df_4_plot_g3[,"relative_bp"],y =rep(0,length(df_4_plot_g3[,'relative_bp'])),df=10),
      lwd=3,col='black')
title(main="KO specific 3h DEG")
legend(-540, leg_y, # 60 or 0.19
       c("plus strand/Up DEG",
         "minus strand/Up DEG",
         "plus strand/Down DEG",
         "minus strand/Down DEG"),
       col = c("red","blue","orange",'skyblue'), lty = 1, cex = 0.7)





###########################################################
#        TFBS prediction
###########################################################
# https://bioconductor.org/packages/devel/bioc/vignettes/TFBSTools/inst/doc/TFBSTools.html#s4-classes-in-tfbstools
