pacman::p_load(dplyr, data.table, ComplexHeatmap, RColorBrewer, circlize, gplots)

#@@@@@ TM_EDA_PCA @@@@@#
TM_EDA_PCA = function(matrix, group, colors, axes = c(1,2),
                      center = F, scale = F, title = "",
                      
                      # varfilter = F,
                      # varfilter_cutoff = 0.7,
                      
                      selection_method = NULL, # IQR VAR
                      selection_n = 0.7,
                      
                      show_label = "none",  # none all
                      label_size = 4, pointsize = 2,
                      legend_label_size = 10,
                      legend_title_label_size = 10,
                      axis_label_size = 7.5,
                      axis_title_label_size = 7.5,
                      addEllipses = T, ellipse.type = "norm" # confidence
){
  if(is.null(selection_method)){
    matrix = matrix
  } else if(selection_method=="IQR"){
    if(selection_n <= 1){
      selection_n = round(nrow(matrix) * selection_n)
    }
    mvgenes = order(apply(matrix, 1, function(x) IQR(x)), decreasing = T)[1:selection_n]
    matrix = matrix[mvgenes, ]
  } else if(selection_method=="VAR"){
    if(selection_n > 1){
      selection_n = 1 - selection_n/nrow(matrix)
    }
    tmp = ExpressionSet(as.matrix(matrix))
    tmp = varFilter(tmp, var.cutoff = selection_n)
    matrix = tmp@assayData$exprs
  }
  
  pca_dt <- prcomp(t(matrix), center = center, scale. = scale)
  res = fviz_pca_ind(pca_dt, axes = axes, col.ind = factor(group), palette = colors,
                     addEllipses = addEllipses, ellipse.type = ellipse.type,
                     invisible="quali", legend.title = "Groups", title = title, 
                     label = show_label, labelsize = label_size, pointsize =pointsize,
                     repel=T, pointshape=16) + 
    labs(x=sprintf("PC%s (%s%s)", axes[1],
                   as.integer(summary(pca_dt)$importance[2,axes[1]] *100),"%"),
         y=sprintf("PC%s (%s%s)", axes[2],
                   as.integer(summary(pca_dt)$importance[2,axes[2]] *100),"%")) +
    theme(title=element_text(size=10), text=element_text(size=legend_label_size),
          legend.title=element_text(size=legend_title_label_size),
          axis.text=element_text(size=axis_label_size),
          axis.title=element_text(size=axis_title_label_size, face  = 'bold'))
  bi_var_plot = fviz_pca_biplot(pca_dt, invisible ="ind")
  return(list("PCA_res" = pca_dt, "plot"=res, 'bi'=bi_var_plot))
  }

#@@@@@ TM_EDA_Heatmap @@@@@#
TM_EDA_Heatmap = function(mat,
                          selection_method = NULL, selection_n = 0.5,
                          top_annot_df = "no", bottom_annot_df = NULL,
                          show_row_names = F, row_names_gp = 12,
                          show_row_dend = T,
                          show_col_names = F, column_names_side = "bottom",
                          column_names_rot= 0, col_names_gp = 12,
                          width = 20, height = 8,
                          cluster_col =T, cluster_row =T,
                          row_km = 1,
                          row_split = NULL, col_split = NULL,
                          row_order = NULL, row_labels = NULL,
                          column_order  = NULL,
                          name="Exp",
                          heatmap_legend_param = NULL, show_legend = T,
                          row_title = character(0), row_title_rot = 0,row_title_gp=13.2,
                          column_labels = NULL,
                          column_title = character(0), col=NULL,
                          left_annotation = NULL, right_annotation = NULL,
                          do_scale = T, seed = 1234, column_dend_height = 2){
  # mat
  # col: SampleID, row: Gene
  if(is.null(row_labels)){ row_labels = rownames(mat) }
  if(is.null(column_labels)){ column_labels = colnames(mat) }
  if(mode(top_annot_df) == "S4"){
    print("Custom top annotation detected.")
  } else if(is.character(top_annot_df)){
    if(top_annot_df == "no"){
      top_annot_df = NULL
    } else {
      return(F)
    }
  } else if(is.data.frame(top_annot_df)){
    if(F %in% c(colnames(mat) == rownames(top_annot_df))){
      print("Error: ")
      print("top_annot_df의 sample 열에는 mat의 sampleIDs가 담겨있어야 함.")
      return(F)
    }
  } else {
    return(F)
  }
  
  # Selection
  if(is.null(selection_method)){
    mat_selected = mat
  } else if(selection_method=="IQR"){
    if(selection_n <= 1){
      selection_n = round(nrow(mat) * selection_n)
    }
    mvgenes = order(apply(mat, 1, function(x) IQR(x)), decreasing = T)[1:selection_n]
    mat_selected = mat[mvgenes, ]
    row_labels = rownames(mat_selected)
  } else if(selection_method=="VAR"){
    if(selection_n > 1){
      selection_n = 1 - selection_n/nrow(mat)
    }
    tmp = ExpressionSet(as.matrix(mat))
    tmp = varFilter(tmp, var.cutoff = selection_n)
    mat_selected = tmp@assayData$exprs
    row_labels = rownames(mat_selected)
  }
  
  # Z scaling
  if(do_scale){
    mat_selected_scaled = t(scale(t(mat_selected)))
  } else {
    mat_selected_scaled = mat_selected
  }
  
  # Top Annotation   ...   set.seed(seed)
  if(is.null(top_annot_df)){
    ha = NULL
  } else if(mode(top_annot_df) == "S4"){
    ha = top_annot_df
  } else {
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector=unlist(mapply(brewer.pal,qual_col_pals$maxcolors,rownames(qual_col_pals)))
    col_li = lapply(top_annot_df, function(x) {
      if(is.numeric(x)){
        return(NA)
      } else {
        return_vec = sample(col_vector, size = length(unique(x)), replace = FALSE)
        names(return_vec) = unique(x)
        return(return_vec)
      } 
    }) %>% .[!is.na(.)]
    ha = HeatmapAnnotation(df = top_annot_df, col = col_li, show_legend = show_legend)
  }
  
  # Heatmap legend params
  if(is.null(heatmap_legend_param)){
    heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                title_position = "topleft")
  } else {
    if(!"legend_height" %in% names(heatmap_legend_param)){
      heatmap_legend_param["legend_height"] = unit(4, "cm")
    }
    
    if(!"title_position" %in% names(heatmap_legend_param)){
      heatmap_legend_param["title_position"] = "topleft"
    }
  }

  # Heatmap
  hm = Heatmap(mat_selected_scaled, top_annotation=ha,
               bottom_annotation=bottom_annot_df,
               
               cluster_rows=cluster_row, cluster_columns=cluster_col, 
               show_row_names=show_row_names, show_column_names=show_col_names,
               show_row_dend=show_row_dend, row_labels=row_labels,
               row_names_gp = gpar(fontsize = row_names_gp),
                
               width = unit(width, "cm"), height = unit(height, "cm"),
               column_names_side = column_names_side,
               column_names_rot = column_names_rot, cluster_row_slices=F,
               name = name, row_split = row_split, column_split = col_split,
               
               row_title = row_title, row_title_rot = row_title_rot,
               row_title_gp = gpar(fontsize = row_title_gp),
               column_title = column_title, cluster_column_slices=F,
               column_labels = column_labels,
               column_names_gp = gpar(fontsize = col_names_gp),
               
               column_dend_height = unit(column_dend_height, "cm"),
               # column_split = c( rep('A', 9), rep('B',8)),
               row_km = row_km, row_order=row_order, column_order = column_order,
               col=col,#colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               left_annotation = left_annotation, right_annotation = right_annotation,
               heatmap_legend_param = heatmap_legend_param)
  return(hm)
}


# draw function
# https://jokergoo.github.io/ComplexHeatmap/reference/draw-HeatmapList-method.html