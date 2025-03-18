

# Utils_TCGA_save_expr_data
Utils_TCGA_save_expr_data = function(x, save_path){
  
  func_res = sapply(x, function(y){
    
    if(file.exists(paste0(save_path, y, "expr.rds"))){
      res = readRDS(paste0(save_path, y, "expr.rds"))
      
    } else {
      proj = paste0("TCGA-", y)
      query.exp=GDCquery(project=proj, legacy = TRUE, data.category = "Gene expression",
                         data.type="Gene expression quantification", file.type="results",
                         experimental.strategy="RNA-Seq", platform = "Illumina HiSeq",
                         sample.type = c("Primary Tumor","Solid Tissue Normal"))
      GDCdownload(query.exp, directory = save_path)
      res <- GDCprepare(query = query.exp,directory = save_path, save = FALSE)
      saveRDS(res, file = paste0(save_path, y, "expr.rds"))
      
    }
    return(res)
    
  })
  
  return(func_res)
}


# Utils_get_genes_groupe
Utils_get_genes_group = function(x, gene, prop=0.2){
  gene_exp = x[gene,]
  qt = quantile(gene_exp, p=c(prop, 1-prop))
  samplesLow = names(gene_exp[gene_exp <= qt[1]])
  samplesHigh = names(gene_exp[gene_exp >= qt[2]])
  return(list(Low=samplesLow, High=samplesHigh))
}



# Utils_get_DEG_EA_res
Utils_get_DEG_res = function(x, gene, prop=0.2, fdr.cut = 0.05, logFC.cut = 1,
                             control_sample = NULL, case_sample = NULL){
  # control_sample, case_sample : sample vector OR NULL
  # control - case
  if(is.null(control_sample) & is.null(case_sample)){
    control_sample = Utils_get_genes_group(x, gene, prop=prop)$Low
    case_sample = Utils_get_genes_group(x, gene, prop=prop)$High
    print("Grouping via Gene expression")
    print(paste("Control: ", length(control_sample)))
    print(paste("Case: ", length(case_sample)))
  } else {
    control_sample = control_sample
    case_sample = case_sample
    print("Grouping via Sample type")
    print(paste("Control: ", length(control_sample)))
    print(paste("Case: ", length(case_sample)))
  }

  dataDEGs <- TCGAanalyze_DEA(mat1 = x[, control_sample], mat2 = x[, case_sample],
                              Cond1type = "Control", Cond2type = "Case",
                              fdr.cut = fdr.cut , logFC.cut = logFC.cut,
                              method = "glmLRT")
  # DEGs table with expression values in normal and tumor samples
  dataDEGsFiltLevel <- TCGAanalyze_LevelTab( FC_FDR_table_mRNA = dataDEGs,
                                             typeCond1 = "Control", typeCond2 = "Case",
                                             TableCond1 = x[, control_sample],
                                             TableCond2 = x[, case_sample])
  return(list("DEG_res"=dataDEGs, "DEG_filt"=dataDEGsFiltLevel,
              "Genelist" = rownames(dataDEGsFiltLevel)))
}



# Utils_get_EA_res
Utils_get_EA_res = function(GeneList){
  ansEA <- TCGAanalyze_EAcomplete(
    TFname = "DEA genes Control Vs Case", RegulonList = GeneList)
  
  return(ansEA)
}




# Utils_get_EA_plot
Utils_get_EA_plot = function(x, filename=NULL, nBar = 10){
  TCGAvisualize_EAbarplot(
    tf = rownames(x$EA_res$ResBP), filename = filename,
    GOBPTab = x$EA_res$ResBP, GOCCTab = x$EA_res$ResCC, 
    GOMFTab = x$EA_res$ResMF, PathTab = x$EA_res$ResPat,
    nRGTab = x$Genelist, nBar = nBar,fig.height = 10, fig.width = 10)
}
