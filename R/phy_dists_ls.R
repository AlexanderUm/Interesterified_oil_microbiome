#-------------------------------------------------------------------------------
# Calculate several dissimilarities distances and return a list with them
#-------------------------------------------------------------------------------\
phy_dist_ls <- function(phylo, 
                        dists =  c("unifrac", "wunifrac", 
                                   "jaccard", "bray")) {
  
  require(phyloseq)
  
  ls.dist <- list()

    for (i.dist in dists) {
      
      if(i.dist == "jaccard") {
        
        dist.inst <- phyloseq::distance(phylo, i.dist, binary = TRUE)
        
      } else {
        
        dist.inst <- phyloseq::distance(phylo, i.dist)
        
      }
      
      ls.dist[[i.dist]] <- dist.inst
    }
    
  return(ls.dist)
  
}
