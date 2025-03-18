pacman::p_load(dplyr, data.table, limma)

#@@@@@ TM_DEG_Limma @@@@@#
TM_DEG_Limma =  function(Rawcount, clinical,
                         sample_colname = "sample", var_cut=0.5,
                         Group, case, control, padjmethod="BH", verb=F, test=F){
  if(is.null(sample_colname)){
    case_samples = rownames(clinical)[clinical[Group] ==case ] # Group case
    control_samples=rownames(clinical)[clinical[Group]==control] # Group control
  } else {
    case_samples = clinical[sample_colname][ clinical[ Group ] == case ] # Group case
    control_samples=clinical[sample_colname][clinical[Group]==control] # Group control
  }
  
  HT = Rawcount %>% as.data.frame() %>% 
    dplyr::select(all_of(c(case_samples, control_samples)))
  
  if(!is.null(var_cut)){
    print("*******************")
    if(verb){print(paste("Var filter ... Before: ",nrow(HT)))}
  
    HT = ExpressionSet(as.matrix(HT))
    HT = varFilter(HT, var.cutoff = var_cut)
    HT =HT@assayData$exprs
    if(verb){print(paste("Var filter ... After: ",nrow(HT)))}
  }
  
  # HT_group=c(paste("A", rep(case, length(case_samples)), sep = "_"),
  #            paste("B", rep(control, length(control_samples)), sep="_")) %>% as.factor()
  
  HT_group=c(rep("A", length(case_samples)),
             rep("B", length(control_samples))) %>% as.factor()
  
  design = model.matrix(~HT_group+0)
  colnames(design) = levels(HT_group)
  
  fit = lmFit(HT, design)

  constrast.matrix = makeContrasts(contrasts = paste0(levels(HT_group)[1], "-",
                                                       levels(HT_group)[2]),
                                   levels = levels(HT_group))
  
  fit.cont = contrasts.fit(fit, constrast.matrix)
  efit = limma::eBayes(fit.cont)
  result = topTable(efit, n=nrow(efit), adjust.method = padjmethod)
  
  DEGs <- result 
  DEGs[, 'probeID'] <- rownames(DEGs)
  
  
  if(test){
    design = model.matrix(~HT_group)
    colnames(design) = levels(HT_group)
# <<<<<<< HEAD
    # print(design)
    fit <- lmFit(HT, design)
    # print(fit)
    fit <- eBayes(fit, legacy =T)
    # print(fit)
    tt <- topTable(fit, coef=2, number = Inf)
# =======
#     fit <- lmFit(HT, design)
#     fit <- eBayes(fit)
#     tt <- topTable(fit, coef=2, n=Inf)
# >>>>>>> parent of 5a66854 (ck)
    tt[, 'probeID'] <- rownames(tt)
    return(tt)
  } else {
    return(DEGs)
  }
}



TM_CUSTOM_geneset_LIMMA = function(Rawcount, clinical,
                                   sample_colname = "sample", var_cut=0.5,
                                   Group, case, control, padjmethod="BH", verb=F, test=F){
  
  my_ebayes1 = function (fit, proportion = 0.01, stdev.coef.lim = c(0.1, 4), 
                         trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1), 
                         legacy = NULL) 
  {
    if (!is.list(fit)) 
      stop("fit is not a valid MArrayLM object")
    if (is.logical(trend) && trend && is.null(fit$Amean)) 
      stop("Need Amean component in fit to estimate trend")
    eb <- my_ebayes2(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim, 
                     trend = trend, robust = robust, winsor.tail.p = winsor.tail.p, 
                     legacy = legacy)
    fit$df.prior <- eb$df.prior
    fit$s2.prior <- eb$s2.prior
    fit$var.prior <- eb$var.prior
    fit$proportion <- proportion
    fit$s2.post <- eb$s2.post
    fit$t <- eb$t
    fit$df.total <- eb$df.total
    fit$p.value <- eb$p.value
    fit$lods <- eb$lods
    if (!is.null(fit$design) && is.fullrank(fit$design)) {
      F.stat <- classifyTestsF(fit, fstat.only = TRUE)
      fit$F <- as.vector(F.stat)
      df1 <- attr(F.stat, "df1")
      df2 <- attr(F.stat, "df2")
      fit$F.p.value <- pf(fit$F, df1, df2, lower.tail = FALSE)
    }
    fit
  }
  
  
  my_ebayes2 = function (fit, proportion = 0.01, stdev.coef.lim = c(0.1, 4), 
                         trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1), 
                         legacy = NULL) 
  {
    coefficients <- fit$coefficients
    stdev.unscaled <- fit$stdev.unscaled
    sigma <- fit$sigma
    df.residual <- fit$df.residual
    if (is.null(coefficients) || is.null(stdev.unscaled) || 
        is.null(sigma) || is.null(df.residual)) 
      stop("No data, or argument is not a valid lmFit object")
    if (identical(max(df.residual), 0)) 
      stop("No residual degrees of freedom in linear model fits")
    if (!any(is.finite(sigma))) 
      stop("No finite residual standard deviations")
    if (is.logical(trend)) {
      if (trend) {
        covariate <- fit$Amean
        if (is.null(covariate)) 
          stop("Need Amean component in fit to estimate trend")
      }
      else {
        covariate <- NULL
      }
    }
    else {
      if (!is.numeric(trend)) 
        stop("trend should be either a logical scale or a numeric vector")
      if (!identical(length(sigma), length(trend))) 
        stop("If trend is numeric then it should have length equal to the number of genes")
      covariate <- trend
    }
    out <- squeezeVar_1(sigma^2, df.residual, covariate = covariate, 
                        robust = robust, winsor.tail.p = winsor.tail.p, legacy = legacy)
    out$s2.prior <- out$var.prior
    out$s2.post <- out$var.post
    out$var.prior <- out$var.post <- NULL
    out$t <- coefficients/stdev.unscaled/sqrt(out$s2.post)
    df.total <- df.residual + out$df.prior
    df.pooled <- sum(df.residual, na.rm = TRUE)
    df.total <- pmin(df.total, df.pooled)
    out$df.total <- df.total
    out$p.value <- 2 * pt(-abs(out$t), df = df.total)
    var.prior.lim <- stdev.coef.lim^2/median(out$s2.prior)
    out$var.prior <- tmixture.matrix(out$t, stdev.unscaled, 
                                     df.total, proportion, var.prior.lim)
    if (anyNA(out$var.prior)) {
      out$var.prior[is.na(out$var.prior)] <- 1/out$s2.prior
      warning("Estimation of var.prior failed - set to default value")
    }
    r <- rep(1, NROW(out$t)) %o% out$var.prior
    r <- (stdev.unscaled^2 + r)/stdev.unscaled^2
    t2 <- out$t^2
    Infdf <- out$df.prior > 10^6
    if (any(Infdf)) {
      kernel <- t2 * (1 - 1/r)/2
      if (any(!Infdf)) {
        t2.f <- t2[!Infdf]
        r.f <- r[!Infdf]
        df.total.f <- df.total[!Infdf]
        kernel[!Infdf] <- (1 + df.total.f)/2 * log((t2.f + 
                                                      df.total.f)/(t2.f/r.f + df.total.f))
      }
    }
    else kernel <- (1 + df.total)/2 * log((t2 + df.total)/(t2/r + 
                                                             df.total))
    out$lods <- log(proportion/(1 - proportion)) - log(r)/2 + 
      kernel
    out
  }
  
  squeezeVar_1 = function (var, df, covariate = NULL, robust = FALSE, winsor.tail.p = c(0.05, 
                                                                                        0.1), legacy = NULL) 
  {
    print("custom")
    n <- length(var)
    if (identical(n, 0L)) 
      stop("var is empty")
    # if (n < 3L)
    #   print('n<3')
    #   return(list(var.post = var, var.prior = var, df.prior = 0))
    if (length(df) > 1L) 
      var[df == 0] <- 0
    if (is.null(legacy)) {
      dfp <- df[df > 0]
      legacy <- identical(min(dfp), max(dfp))
    }
    if (legacy) {
      if (robust) {
        fit <- fitFDistRobustly(var, df1 = df, covariate = covariate, 
                                winsor.tail.p = winsor.tail.p)
        df.prior <- fit$df2.shrunk
      }
      else {
        fit <- fitFDist(var, df1 = df, covariate = covariate)
        df.prior <- fit$df2
      }
    }
    else {
      fit <- fitFDistUnequalDF1(var, df1 = df, covariate = covariate, 
                                robust = robust)
      df.prior <- fit$df2.shrunk
      if (is.null(df.prior)) 
        df.prior <- fit$df2
    }
    print(df.prior)
    if (anyNA(df.prior)) 
      stop("Could not estimate prior df")
    var.post <- squeezevar_2(var = var, df = df, var.prior = fit$scale, 
                             df.prior = df.prior)
    list(df.prior = df.prior, var.prior = fit$scale, var.post = var.post)
  }
  
  squeezevar_2 = function (var, df, var.prior, df.prior) 
  {
    m <- max(df.prior)
    if (is.finite(m)) 
      return((df * var + df.prior * var.prior)/(df + df.prior))
    n <- length(var)
    if (identical(length(var.prior), n)) {
      var.post <- var.prior
    }
    else {
      var.post <- rep_len(var.prior, length.out = n)
    }
    m <- min(df.prior)
    if (m > 1e+100) 
      return(var.post)
    i <- which(is.finite(df.prior))
    if (length(df) > 1L) 
      df <- df[i]
    df.prior <- df.prior[i]
    var.post[i] <- (df * var[i] + df.prior * var.post[i])/(df + 
                                                             df.prior)
    var.post
  }
  
  
  if(is.null(sample_colname)){
    case_samples = rownames(clinical)[clinical[Group] ==case ] # Group case
    control_samples=rownames(clinical)[clinical[Group]==control] # Group control
  } else {
    case_samples = clinical[sample_colname][ clinical[ Group ] == case ] # Group case
    control_samples=clinical[sample_colname][clinical[Group]==control] # Group control
  }
  
  HT = Rawcount %>% as.data.frame() %>% 
    dplyr::select(all_of(c(case_samples, control_samples)))
  # print(dim(HT))
  if(!is.null(var_cut)){
    print("*******************")
    if(verb){print(paste("Var filter ... Before: ",nrow(HT)))}
    
    HT = ExpressionSet(as.matrix(HT))
    HT = varFilter(HT, var.cutoff = var_cut)
    HT =HT@assayData$exprs
    if(verb){print(paste("Var filter ... After: ",nrow(HT)))}
  }
  
  # HT_group=c(paste("A", rep(case, length(case_samples)), sep = "_"),
  #            paste("B", rep(control, length(control_samples)), sep="_")) %>% as.factor()
  
  HT_group=c(rep("A", length(case_samples)),
             rep("B", length(control_samples))) %>% as.factor()
  
  design = model.matrix(~HT_group+0)
  colnames(design) = levels(HT_group)
  
  fit = lmFit(HT, design)
  
  constrast.matrix = makeContrasts(contrasts = paste0(levels(HT_group)[1], "-",
                                                      levels(HT_group)[2]),
                                   levels = levels(HT_group))
  
  fit.cont = contrasts.fit(fit, constrast.matrix)
  efit = limma::eBayes(fit.cont)
  result = topTable(efit, n=nrow(efit), adjust.method = padjmethod)
  
  DEGs <- result 
  DEGs[, 'probeID'] <- rownames(DEGs)
  
  
  if(test){
    design = model.matrix(~HT_group)
    # print(design)
    colnames(design) = levels(HT_group)
    print(design)
    fit <- lmFit(HT, design)
    print(fit)
    fit <- my_ebayes1(fit, legacy =T)
    print(fit)
    tt <- topTable(fit, coef=2, number = Inf)
    tt[, 'probeID'] <- rownames(tt)
    return(tt)
  } else {
    return(DEGs)
  }
  
}





TM_DE_pws = function(Rawcount, clinical,
                     parametric=F,
                     Group, case, control, padjmethod="BH"){
  clinical = clinical %>% tibble::rownames_to_column("sample")
  dk1 = Rawcount %>% mutate(items = rownames(.)) %>% gather(sample, expression, -items) %>%
    inner_join(., clinical, by = "sample")
  dk4 = sapply(unique(dk1$items), function(x){
    dk2 = dk1 %>% filter(items == x) %>% filter(condition %in% c(case, control)) %>%
      mutate(condition = factor(condition))

    if(parametric){
      dk3 = t.test(expression~condition, dk2)[c("p.value", "stderr")] %>%
        bind_rows(.) %>%  as.data.frame(.)
      tmp = t.test(expression~condition, dk2)[["estimate"]]
      dk3$logFC = tmp[1]-tmp[2]
      rownames(dk3) = x
    } else {
      dk3 = wilcox.test(expression~condition, dk2)[c("p.value")] %>%
        bind_rows(.) %>%  as.data.frame(.)
      rownames(dk3) = x
    }
    
    return(dk3)}, simplify = F) %>% bind_rows(.)
  
  dk4[, 'probeID'] <- rownames(dk4)
  return(dk4)
}



TM_DE_Limma_anova = function(score_matrix, info, target_col){
  info["TG"] = info[target_col]
  
  design_popGen <- model.matrix( ~0 + TG, data = info)
  print(design_popGen)
  VarNames = colnames(design_popGen)
  tmp = c()
  for(i in seq(length(VarNames))){
    targets = seq(i, length(VarNames))
    if(i != length(VarNames)){
      for(tg_i in targets){
        if(i != tg_i){
          tmptmp =c(paste0(VarNames[i], "-", VarNames[tg_i]))
          tmp = c(tmp, tmptmp)
          
        }
      }
    }
  }
  
  contrast_popGen = makeContrasts(contrasts=tmp, levels = design_popGen)
  
  fit1 <- lmFit(subset(score_matrix, T, rownames(tmp_info)), design_popGen)
  fit2 <- contrasts.fit(fit1, contrasts = contrast_popGen)
  fit3 <- eBayes(fit2, trend = T, robust = T)
  
  tmp = topTable(fit3, number = nrow(fit3))
  
  return(tmp)
}
