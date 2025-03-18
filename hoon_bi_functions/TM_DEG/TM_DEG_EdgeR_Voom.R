pacman::p_load(dplyr, limma, edgeR, data.table)

TM_DEG_EdgeR_Voom = function(mat, case, control, QC=F, title = ""){
  xMat = mat %>% dplyr::select(all_of(c(case, control)))
  group = factor(c(rep('case', length(case)), rep('control', length(control))),
                 c('case','control'))
  
  deLi = DGEList(xMat)
  deLi$samples$group = group
  
  L = mean(deLi$samples$lib.size) * 1e-6; M = median(deLi$samples$lib.size) * 1e-6
  keep.exprs = filterByExpr(deLi, group = deLi$samples$group)
  deLi_filtered = deLi[keep.exprs, , keep.lib.sizes = F]
  
  if(QC){
    # QC
    # non-filtered lcpm
    lcpm = cpm(deLi, log = T)
    lcpm.cutoff = log2(10/M + 2/L)
    par(mfrow= c(1,2))
    plot(density(lcpm[,1]), lwd = 2, ylim = c(0, 0.26), las = 2, main="", xlab="")
    title(main = title, sub = "A. Expected Count", xlab = "Log-CPM")
    abline(v = lcpm.cutoff, lty = 3)
    for(i in 2:ncol(deLi)){
      den = density(lcpm[, i])
      lines(den$x, den$y, lwd = 2)
    }
    # filtered lcpm
    lcpm = cpm(deLi_filtered, log = T)
    plot(density(lcpm[,1]), lwd = 2, ylim = c(0, 0.26), las = 2, main="", xlab="")
    title(sub = "B. Filtered Expected Count", xlab = "Log-CPM")
    abline(v = lcpm.cutoff, lty = 3)
    for(i in 2:ncol(deLi)){
      den = density(lcpm[, i])
      lines(den$x, den$y, lwd = 2)
    }
  }
  
  deLi_filtered = calcNormFactors(deLi_filtered, method = "upperquartile") # TMM vs Upperquantile
  # lcpm = cpm(deLi, log = T) # <- PCA에 쓸 lcpm임.
  
  design = model.matrix(~0 + group)
  colnames(design) = gsub("group", "", colnames(design))
  
  contr.matrix = makeContrasts(caseVScontrol = case - control,
                               levels = colnames(design))
  
  v = voom(deLi_filtered, design, plot = QC)
  vfit = lmFit(v, design)
  vfit = contrasts.fit(vfit, contrasts = contr.matrix)
  if(QC){
    efit = eBayes(vfit)
    plotSA(efit, main = paste0(title, ": Final model: Mean-variance trend"))
  }
  
  tfit = treat(vfit, lfc = 0.58);
  DEGsTreat = topTreat(tfit, n = Inf)
  return(list(DEGsTreat = DEGsTreat,
              voomExp = as.data.frame(v$E)))}
