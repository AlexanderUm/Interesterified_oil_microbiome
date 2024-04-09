#-----------------------------------------------------------------------------
# Function: Extract OTU tables for Maaslin from a phyloseq list
#-----------------------------------------------------------------------------
otu_tab_extract <- function(phy_named_list) {
  
  otu.tabs.out <- list()
  
  for(i.phy in names(phy_named_list)) {
    
    otu.tab <- phy_named_list[[i.phy]] %>% 
      otu_table() %>% 
      as.matrix()
    
    if(taxa_are_rows(phy_named_list[[i.phy]])) {
      
      otu.tabs.out[[i.phy]] <- otu.tab %>% 
        t() %>% 
        as.data.frame()
      
    } else {
      
      otu.tabs.out[[i.phy]] <- otu.tab %>% 
        as.data.frame()
    }
    
  }
  
  return(otu.tabs.out)
  
}