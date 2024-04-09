#-------------------------------------------------------------------------------
# Correlation between a list of data frames (otu tables) and a data frame 
#-------------------------------------------------------------------------------
bulk_corr_matrix <- function(otu_tabs_list, 
                             meta_data, 
                             second_df, 
                             strata_cols_meta = NULL,
                             prev_cut_offs = c(0.25), 
                             corr_method = "spearman", 
                             p_adj_method = "BH") {
  
  #-----------------------------------------------------------------------------
  prev_filt_otu_tab <- function(otu_tab, meta_data, 
                                group_col = FALSE, prev_cutoff) {
    
    if(is.logical(group_col)) {
      
      group.col <- rep("a", nrow(otu_tab))
      
    } else {
      
      group.col <- pull(meta_data, group_col) }
    
    filt.vec <- otu_tab %>%
      mutate(across(where(is.numeric), function(x){ifelse(x == 0, 0, 1)})) %>%
      mutate(group_col = group.col)  %>%
      group_by(group_col) %>%
      summarise(across(where(is.numeric), 
                       function(x){sum(x)/length(x)}), .groups = "drop") %>%
      dplyr::select(-group_col) %>% 
      apply(., 2, function(x){any(x >= prev_cutoff)}, simplify = TRUE)
    
    tab.out <- otu_tab[, filt.vec]
    
    return(tab.out)
    
  }
  #-----------------------------------------------------------------------------
  
  
  # Create data frame for stratification ---------------------------------------
  if(is.null(strata_cols_meta)) {
    
    strata.data <- data.frame(Strata = rep("No Stratification", 
                                           nrow(meta_data)))
    
  } else {  
    
    strata.data <- meta_data %>% 
      dplyr::select(all_of(strata_cols_meta)) %>% 
      unite(Strata, all_of(strata_cols_meta))
  }
  #-----------------------------------------------------------------------------
  
  
  # Create a list of lists that will be used later 
  res.out <- list()
  
  #--- Loop through: otu tables ---
  for (i.otu in names(otu_tabs_list)) {
    
    # --- Loop through: strata ---
    for (i.strata in unique(strata.data$Strata)) {
      
      # Subset data ------------------------------------------------------------
      otu.tab.i.strat <- otu_tabs_list[[i.otu]] %>% 
        .[strata.data$Strata == i.strata, ]
      
      second.df.i.strat <- second_df[strata.data$Strata == i.strata, ]
      
      meta.data.strata <- meta_data[strata.data$Strata == i.strata, ]
      
      
      # Row names check --------------------------------------------------------
      if(!identical(rownames(otu.tab.i.strat), 
                    rownames(meta.data.strata))) {
        stop("Row names are not identical")}
      
      if(!identical(rownames(otu.tab.i.strat), 
                    rownames(second.df.i.strat))) {
        stop("Row names are not identical")}
      #-------------------------------------------------------------------------
      
      # Filter based on prevalence 
      otu.tab.i.strat.f <- prev_filt_otu_tab(otu_tab = otu.tab.i.strat, 
                                             meta_data = meta.data.strata, 
                                             prev_cutoff = min(prev_cut_offs))
      
      # Stop if no taxa left after filtering----------------------------------
      if(ncol(otu.tab.i.strat.f) < 1) {
        paste0("No taxa left after prevalece filtering. \n", 
               "Cut Off :", i.cut.prev, "\n",
               "Otu table: ", i.otu, "\n",
               "Starta :", i.strata) %>% 
          stop(.)}
      #-----------------------------------------------------------------------
      
      #-----------------------------------------------------------------------
      # Correlation loop
      cor.res.comb <- NULL
      
      for (i.second in colnames(second.df.i.strat)) {
        
        cor.res <- NULL 
        
        for (i.tax in colnames(otu.tab.i.strat.f)) {
          
          cor.res <- cor.test(second.df.i.strat[, i.second],
                              otu.tab.i.strat.f[, i.tax],
                              method = corr_method, na.rm=TRUE) %>%
            tidy() %>%
            mutate(Taxa = i.tax) %>%
            bind_rows(cor.res, .)
          
        }
        
        cor.res.comb <- cor.res %>% 
          mutate(Corre_vector = i.second) %>% 
          bind_rows(cor.res.comb, .)
      }
      
      
      # Add strata and otu table names
      cor.res.comb <- cor.res.comb %>% 
        mutate(Strata = i.strata, 
               OTU_table_name = i.otu)
      
      
      # Loop though prevalence cutoff and adjust q value per metabolite separately 
      for (i.prev in prev_cut_offs) {
        
        # Filter based on prevalence 
        sel.taxa <- prev_filt_otu_tab(otu_tab = otu.tab.i.strat.f, 
                                      meta_data = strata.data, 
                                      prev_cutoff = i.prev) %>% 
          colnames(.) 
        
        cor.res.comb.f <-  cor.res.comb %>% 
          group_by(Corre_vector) %>% 
          filter(Taxa %in% sel.taxa) %>% 
          mutate(Prevalence_cutoff = i.prev, 
                 qval = p.adjust(p.value, 
                                 method = p_adj_method)) %>% 
          ungroup()
        
        firs.ls.name <- paste0(i.otu, "--", i.prev)
        
        res.out[[firs.ls.name]][[i.strata]] <- cor.res.comb.f
        
        
      }
      
    }
    
  }
  
  return(res.out)
}

