
# TM_PRC_edger_voom_Expr
TM_PRC_edger_voom_Expr = function(mat, sample, group, QC=T, title = ""){
  xMat = mat %>% dplyr::select(all_of(sample))
  
  deLi = DGEList(xMat)
  deLi$samples$group = group
  
  L = mean(deLi$samples$lib.size) * 1e-6; M = median(deLi$samples$lib.size) * 1e-6
  keep.exprs = filterByExpr(deLi, group = deLi$samples$group)
  deLi_filtered = deLi[keep.exprs, , keep.lib.sizes = F]
  
  # if(QC){
  #   # QC
  #   # non-filtered lcpm
  #   lcpm = cpm(deLi, log = T)
  #   lcpm.cutoff = log2(10/M + 2/L)
  #   par(mfrow= c(1,2))
  #   plot(density(lcpm[,1]), lwd = 2, ylim = c(0, 0.26), las = 2, main="", xlab="")
  #   title(main = title, sub = "A. Expected Count", xlab = "Log-CPM")
  #   abline(v = lcpm.cutoff, lty = 3)
  #   for(i in 2:ncol(deLi)){
  #     den = density(lcpm[, i])
  #     lines(den$x, den$y, lwd = 2)
  #   }
  #   # filtered lcpm
  #   lcpm = cpm(deLi_filtered, log = T)
  #   plot(density(lcpm[,1]), lwd = 2, ylim = c(0, 0.26), las = 2, main="", xlab="")
  #   title(sub = "B. Filtered Expected Count", xlab = "Log-CPM")
  #   abline(v = lcpm.cutoff, lty = 3)
  #   for(i in 2:ncol(deLi)){
  #     den = density(lcpm[, i])
  #     lines(den$x, den$y, lwd = 2)
  #   }
  # }
  
  deLi_filtered = calcNormFactors(deLi_filtered) # TMM vs Upperquantile
  # lcpm = cpm(deLi, log = T) # <- PCA에 쓸 lcpm임.
  
  design = model.matrix(~0 + group)
  colnames(design) = gsub("group", "", colnames(design))
  par(mfrow=c(1,1))
  v = voom(deLi_filtered, design, plot = QC)
  
  return(list(voomExp = as.data.frame(v$E)))}


# TM_PRC_GE_probe_to_entrez_by_meanmax
TM_PRC_GE_probe_to_entrez_by_IQR = function(ori_GE, Geneanno,
                                               before_rowname_lab, after_rowname_lab){
  
  GE_prc = ori_GE %>% as.data.frame()
  # %>% mutate(exp_mean = rowMeans(.))
  GE_prc$IQR = apply(ori_GE, 1, function(x) IQR(x))
  
  GE_prc[[before_rowname_lab]] = rownames(GE_prc)
  GE_prc = merge(GE_prc, Geneanno) %>% dplyr::select(- !! sym(before_rowname_lab)) %>% 
    dplyr::group_by( !! sym(after_rowname_lab) ) %>%
    filter(IQR == max(IQR)) %>% ungroup() %>%
    tibble::column_to_rownames(after_rowname_lab) %>% dplyr::select(-IQR)
  return(GE_prc)
}


# TM_PRC_Combat_batch_correction
TM_PRC_Combat_batch_correction = function(GE, batch){
  if(!(length(batch) == length(colnames(GE)))){
    print("length of batch and GE column number should be same.")
    return(F)
  }
  
  modcombat = model.matrix(~1,
                           data= data.frame('batch' = batch,
                                            row.names = colnames(GE)))
  combat_res = ComBat(dat=GE, batch=batch, mod=modcombat,
                      par.prior=TRUE, prior.plots=FALSE) %>% as.data.frame()
  return(combat_res)
}
