################################################################################
# Alpha diversity 
#-------------------------------------------------------------------------------
set.seed(549794)


#-------------------------------------------------------------------------------
# Load data
load("output/supp/data.Rdata")

# Custom funciton 
source("R/phy_alpha.R")


#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
ps.inst <- data.ls$PS[[prm.ls$Alpha$Tax_lvl]][[prm.ls$Alpha$Norm]]

alpha.ind <- prm.ls$Alpha$measures

meta <- data.ls$meta

ref.gr <- prm.ls$Data$ref_gr

gr.col <- prm.ls$Data$group_col

res.path <- prm.ls$Data$res_path


# Out folders
dir.create(paste0(res.path, "/diversity/tabs"), recursive = TRUE)
dir.create(paste0(res.path, "/diversity/plots"), recursive = TRUE)

# Results 
res.ls <- list()

################################################################################
# Alpha diversity
################################################################################
alpha.df <- phy_alpha(ps.inst, 
                      measures = alpha.ind ) %>% 
                      bind_cols(., meta) 

alpha.long <- alpha.df %>% 
                  pivot_longer(cols = all_of(alpha.ind), 
                               names_to = "Diversity_index")

# Write results 
write.csv(alpha.df, 
          paste0(res.path, "/diversity/tabs/alpha_per_sample.csv"))


#-------------------------------------------------------------------------------
# Statistical comparison with KW test and Dunn test as post-hoc
#-------------------------------------------------------------------------------
res.kw.all <- NULL

res.dunn.all <- NULL

for (i.ind in alpha.ind) { 
  
  form.inst <- paste0(i.ind, "~", gr.col)
  
  res.kw <-  kruskal.test(as.formula(form.inst), alpha.df) %>% 
             tidy() %>% 
             mutate(Index = i.ind)
  
  res.kw.all <- bind_rows(res.kw.all, res.kw)
  
  if (res.kw$p.value <= 0.05) {
    
    res.dunn.all <- dunnTest(as.formula(form.inst), alpha.df) %>% 
                    .[["res"]] %>% 
                    mutate(Diversity_Index = i, Test = "Dunn Test") %>% 
                    bind_rows(res.dunn.all, .)
    
  }
  
}

# Write results 
write.csv(x = res.kw.all, 
          file = paste0(res.path, "/diversity/tabs/alpha_kw_res.csv"))

res.ls[["alpha"]][["kw_tab"]] <- res.kw.all


if(!is.null(res.dunn.all)) {
  
  write.csv(x = res.dunn.all, 
            file = paste0(res.path, "/diversity/tabs/alpha_Dunn_res.csv"))
  
  res.ls[["alpha"]][["dunn_tab"]] <- res.kw.all
  
}


#-------------------------------------------------------------------------------
# Plot Alpha diversity 
#-------------------------------------------------------------------------------
# Plot alpha diversity 
alpha.plot <- ggplot(alpha.long, 
                         aes(y = value, 
                             x = .data[[gr.col]])) + 
                      geom_jitter(aes(color = .data[[gr.col]]), 
                                      size = 2, 
                                      alpha = 0.75, 
                                      width = 0.15, 
                                      height = 0) +
                      geom_violin(fill = "black", 
                                      alpha = 0.01) +
                      facet_grid("Diversity_index", 
                                     scales = "free") + 
                      theme_bw() + 
                      scale_color_manual(name = "Diet type", 
                                             values = aes.ls$gr_col) + 
                      theme(axis.text.x = element_text(size=8, 
                                                           face = "bold", 
                                                           family = "serif"), 
                            legend.text = element_text(size=8, 
                                                       face = "bold", 
                                                       family = "serif"), 
                            legend.title = element_text(size=10), 
                            axis.title.y = element_blank()) 


res.ls[["alpha"]][["plot"]] <- alpha.plot


#-------------------------------------------------------------------------------
# Write results
saveRDS(res.ls, paste0(res.path, "/supp/res_ls.RDS"))

# Clean environment 
rm(list=ls())

gc()

