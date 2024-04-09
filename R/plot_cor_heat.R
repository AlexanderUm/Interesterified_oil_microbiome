#-------------------------------------------------------------------------------
# Function: Plot heat maps
#-------------------------------------------------------------------------------
plot_cor_heat <- function(cor_function_results, 
                          estimate_cutoff = 0, 
                          qval_cutoff = 1, 
                          pval_cutoff = 1,
                          group_col = NULL, 
                          group_col_levels = NULL, 
                          estim_display_cut = 0.5) {
  
  tax.list <- cor_function_results %>%
                filter(abs(estimate) >= estimate_cutoff, 
                       qval <= qval_cutoff, 
                       p.value <= pval_cutoff) %>% 
                group_by(Taxa) %>% 
                slice_max(abs(estimate)) %>% 
                arrange(desc(abs(estimate))) %>% 
                pull(Taxa) 
  
  if (!is.null(group_col)) {
    
    cor_function_results <- cor_function_results %>% 
      mutate(Corre_vector = paste0(Corre_vector, "--", get(group_col)))
    
  } 
  
  heat.df <- cor_function_results %>%
                filter(Taxa %in% tax.list) %>%
                select(c("Taxa", "estimate", "Corre_vector")) %>%
                pivot_wider(id_cols = Taxa, 
                            names_from = Corre_vector, 
                            values_from = estimate) %>%
                column_to_rownames("Taxa") %>%
                .[tax.list, ] %>% 
                replace(is.na(.), 0) %>% 
                as.matrix()
  
  col.lab <- gsub("--.*", "", colnames(heat.df)) %>% 
              setNames(colnames(heat.df))
  
  if (is.null(group_col)) {
    
    
    
    heat.plot <- heat.df %>%
                  Heatmap(., name = "estimate",
                          cluster_columns = FALSE,
                          column_gap = unit(0, "mm"),
                          cluster_rows = FALSE,
                          border = TRUE,
                          rect_gp = gpar(col = "gray", lwd = 2),
                          row_names_max_width = max_text_width(rownames(.), 
                            gp = gpar(fontsize = 10)),
                          # column_names_rot = 45,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            if(abs(.[i, j]) >= estim_display_cut)
                              grid.text(sprintf("%.2f", .[i, j]), 
                                        x, y, gp = gpar(fontsize = 8, col="black"))})
    
  } else {
    
    col.split.factor <- factor(sub(".*--", "", colnames(heat.df)), 
                               levels = group_col_levels)
    
    heat.plot <- heat.df %>% 
      Heatmap(., name = "estimate",
              column_split = col.split.factor,
              cluster_columns = FALSE, 
              cluster_column_slices = FALSE,
              cluster_rows = FALSE,
              column_gap = unit(0, "mm"),
              border = TRUE,
              rect_gp = gpar(col = "gray", lwd = 2),
              row_names_max_width = max_text_width(rownames(.), 
                                                   gp = gpar(fontsize = 10)),
              # column_names_rot = 45,
              column_labels = col.lab, 
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(abs(.[i, j]) >= estim_display_cut)
                  grid.text(sprintf("%.2f", .[i, j]), x, y, gp = gpar(fontsize = 8, 
                                                                      col="black"))}
      )
    
    
  }
  
  out.obj <- list(plot = heat.plot, taxa = tax.list)
  
  return(out.obj)
  
}
