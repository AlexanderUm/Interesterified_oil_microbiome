################################################################################
# Beta diversity (overall composition)
#-------------------------------------------------------------------------------
set.seed(5445)


#-------------------------------------------------------------------------------
# Load data
load("output/supp/data.Rdata")
res.ls <- readRDS("output/supp/res_ls.RDS")

# Custom functions
source("R/phy_dists_ls.R")
source("R/phy_pcoa_plot.R")


#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
ps.inst <- data.ls$PS[[prm.ls$Beta$Tax_lvl]][[prm.ls$Beta$Norm]]

dists <- prm.ls$Beta$dists

meta <- data.ls$meta

ref.gr <- prm.ls$Data$ref_gr

gr.col <- prm.ls$Data$group_col

res.path <- prm.ls$Data$res_path

form.prem <- prm.ls$Beta$form

n.perm <- prm.ls$Beta$nperm


################################################################################
# Beta diversity
################################################################################
# Adjust tree root due to warnings
phy_tree(ps.inst) <- ape::multi2di(phy_tree(ps.inst))

# Formula for PERMANOVA
a.form <- paste0("dist.inst ~ ", form.prem)


#-------------------------------------------------------------------------------
# Overall comparison 
#-------------------------------------------------------------------------------
# Calculate distances 
dist.ls <- phy_dist_ls(ps.inst, dists = dists)

# Run Adonis function per distance 
adon.res.ls <- list()

for(i.dist in dists) {
  
  dist.inst <- dist.ls[[i.dist]]
  
  adon.res.ls[[i.dist]] <- adonis2(formula = as.formula(a.form), 
                              data = meta, 
                              by = "terms", 
                              permutations = n.perm, 
                              parallel = 4) %>% 
                          tidy() %>% 
                          mutate(Distance = names(dists[dists == i.dist]),
                                 Formula = a.form) %>% 
                          rbind(., " ")
}

# Write results 
adon.res.ls %>% 
  bind_rows() %>% 
  write.csv(., paste0(res.path, "/diversity/tabs/beta_adonis_comb.csv"))

res.ls[["beta"]][["adonis_comb"]] <- adon.res.ls



#-------------------------------------------------------------------------------
# Paired comparison 
#-------------------------------------------------------------------------------
# List with pairs that will be used
pairs.group <- combn(unique(as.character(meta[[gr.col]])), 2) %>% 
                as.data.frame() %>% 
                as.list()


res.adonis.paired <- list()

for (i.pair in pairs.group) {
  
  ps.sub <- prune_samples(meta[[gr.col]] %in% i.pair, ps.inst) 
  
  phy_tree(ps.sub) <- ape::multi2di(phy_tree(ps.sub))
  
  meta.sub <- filter(meta, meta[[gr.col]] %in% i.pair)
  
  dist.ls.inst <- phy_dist_ls(ps.sub, dists = dists)
  
  for (i.dist in dists) {
    
    dist.inst <- dist.ls.inst[[i.dist]]
    
    obj.name <- paste0(i.dist, ":", paste(i.pair, collapse = "-"))
    
    res.adonis.paired[[obj.name]] <- adonis2(formula = as.formula(a.form), 
                                             data = meta.sub, 
                                             by = "terms", 
                                             permutations = n.perm, 
                                             parallel = 4) %>% 
                                      tidy() %>% 
                                      mutate(Distance = names(dists[dists == i.dist]), 
                                             Pair = obj.name, 
                                             Formula = a.form) %>% 
                                      rbind(., " ")
    
  }
  
}

# Write results 
res.adonis.paired %>% 
        bind_rows() %>% 
        write.csv(., paste0(res.path, 
                            "/diversity/tabs/beta_adonis_paired.csv"))

res.ls[["beta"]][["adonis_pairs"]] <- res.adonis.paired


################################################################################
# Plot ordination
#-------------------------------------------------------------------------------
pcoa.plots.ls <- list() 

for (i.dist in dists) { 
    
    dist.inst <- dist.ls[[i.dist]]
    
    # Make PCoA object
    a.pcoa <- ape::pcoa(dist.inst)
    
    axis.pcoa <- a.pcoa$vectors[, 1:2] %>%
                    as.data.frame() %>%
                    bind_cols(., meta)
    
    # Extract percentage
    var.c <- round((a.pcoa$values$Relative_eig[1:2]*100), 1)
    
    pcoa.p <- ggplot(axis.pcoa,
                     aes(x = Axis.1, 
                         y = Axis.2)) +
              stat_ellipse(geom = "polygon", 
                           aes(fill = .data[[gr.col]]), 
                           alpha = 0.15, level = 0.90, size =2) +
              geom_point(aes(color = .data[[gr.col]]),
                         alpha = 0.75, size = 2) +
              ggtitle(paste0("Distance: ", names(dists[dists == i.dist]))) +
              scale_color_manual(values = aes.ls$gr_col) +
              scale_fill_manual(values = aes.ls$gr_col) +
              coord_fixed() +
              theme_bw() +
              xlab(paste0(colnames(axis.pcoa)[1], " [", var.c[1], "%]")) +
              ylab(paste0(colnames(axis.pcoa)[2], " [", var.c[2], "%]")) +
              labs(fill = "Diet Type",
                   color = "Diet Type") + 
              theme(plot.title = element_text(size=11, face="italic"), 
                    panel.grid = element_blank())
    
    
    pcoa.plots.ls[[i.dist]] <- pcoa.p
    
}


#-------------------------------------------------------------------------------
# Combine plots 
#-------------------------------------------------------------------------------
p.ls.inst.f <- lapply(pcoa.plots.ls, 
                      function(x){x+theme(legend.position="none")})
  
legd.inst <- get_legend(pcoa.plots.ls[[1]])
  
p.grid.1 <- plot_grid(plotlist = p.ls.inst.f, 
                      ncol = 2)
  
p.grid.2 <- plot_grid(p.grid.1, 
                      legd.inst, 
                      ncol = 2, 
                      rel_widths = c(0.8, .2))


res.ls[["beta"]][["plots"]] <- pcoa.plots.ls


#-------------------------------------------------------------------------------
# Combine Alpha Diversity and Beta diversity plots
#-------------------------------------------------------------------------------
alpha.comb <- res.ls$alpha$plot + 
                theme(axis.title.x = element_blank(), 
                      axis.title.y = element_blank(),
                      legend.position="none")

pcoa.comb <- p.grid.2  + 
                theme(legend.text = element_text(size=8, 
                                                 face = "bold", 
                                                 family = "serif"), 
                      legend.title = element_text(size=10))

div.comb.p <- cowplot::plot_grid(alpha.comb, pcoa.comb, 
                                 labels = c("A", "B"), 
                                 rel_widths = c(2.5, 8))

ggsave(paste0(res.path, "/diversity/plots/diversity_comb.png"), 
       div.comb.p, 
       width = 10, 
       height = 5, 
       dpi = 600)


#-------------------------------------------------------------------------------
# Write results
saveRDS(res.ls, paste0(res.path, "/supp/res_ls.RDS"))

# Clean environment 
rm(list=ls())

gc()

