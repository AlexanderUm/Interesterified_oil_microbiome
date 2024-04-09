################################################################################
# Differential abundance 
#-------------------------------------------------------------------------------
set.seed(435975)


#-------------------------------------------------------------------------------
# Load data
load("output/supp/data.Rdata")
res.ls <- readRDS("output/supp/res_ls.RDS")

# Custom functions 
source("R/phy_taxa_filter.R")

#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
tax.lvl <- prm.ls$DA$Tax_lvl

count.norm <- prm.ls$DA$Norm

min.prev <- prm.ls$DA$min.prev

meta <- data.ls$meta

ref.gr <- prm.ls$Data$ref_gr

gr.col <- prm.ls$Data$group_col

res.path <- prm.ls$Data$res_path


# Out folders
dir.create(paste0(res.path, "/DA/tabs"), recursive = TRUE)
dir.create(paste0(res.path, "/DA/plots"), recursive = TRUE)


#-------------------------------------------------------------------------------
# MAASLIN: Extract data
#-------------------------------------------------------------------------------
# Create output directory for maaslin 
out.folder <- paste0(res.path, "/DA/masslin")

maas.out <- list()

res.maas.long <- NULL

for (i.lvl in tax.lvl) {
  
  ps.maas <- data.ls$PS[[i.lvl]][[count.norm]] %>% 
             phy_taxa_filter(., prev = min.prev, group_col = gr.col)
  
 
    # Extract OTU and metadata tables from the phyloseq 
    ps.maas.asv <- ps.maas %>% 
                      otu_table()  %>% 
                      as.matrix()  %>% 
                      as.data.frame() 

    # Prepare maaslin out directory 
    maas.folder <- paste0(out.folder, "/", i.lvl)
      
    dir.create(maas.folder, recursive = TRUE)
      
      maas.res  <- Maaslin2( 
                    input_data =  ps.maas.asv, 
                    input_metadata = meta, 
                    output =  maas.folder, 
                    fixed_effects = gr.col, 
                    reference = paste0(gr.col, ",", ref.gr), 
                    correction = "BH", 
                    cores = 4, 
                    min_abundance = 0, 
                    min_prevalence = 0, 
                    min_variance = 0, 
                    normalization = "CSS", 
                    transform = "LOG", 
                    analysis_method = "LM", 
                    max_significance = 0.1, 
                    plot_heatmap = FALSE)
      
      maas.out[[i.lvl]] <- maas.res
      
      res.maas.long <- rbind(res.maas.long, cbind(maas.res$results, 
                                                  Taxa_Level = i.lvl, 
                                                  Cut_off = min.prev, 
                                                  Reference = ref.gr)) 
      
      # Summarize and write results 
      sum.tab <- maas.res$results %>% 
                    select(feature, value, coef, pval, qval) %>% 
                    pivot_wider(names_from = value, 
                                values_from = c("coef", "pval", "qval"), 
                                names_glue = "{.value}_{value}") 
      
      write.csv(sum.tab, file = paste0(res.path, "/DA/tabs/DA_res_", i.lvl, ".csv"))
      
}


#-------------------------------------------------------------------------------
# Plot Maaslin results 
#-------------------------------------------------------------------------------
res.maas.long.f <- res.maas.long %>% 
                      filter(qval <= 0.2) %>% 
                      mutate(qval_char = ifelse(round(qval, 2) == 0, "q>0.01", 
                                                paste0("q=", round(qval, 2)))) %>% 
                      mutate(qval_posit = ifelse(coef > 0, qval_char, NA), 
                             qval_negat = ifelse(coef < 0, qval_char, NA)) %>% 
                      mutate(Taxa_Level = factor(Taxa_Level, 
                                                 levels = tax.lvl, ), 
                             value = factor(value, levels = c("ICT", "HF", "IHF")), 
                             feature = gsub("__", " ", feature)) 

maas.plot.CT <- ggplot(res.maas.long.f, 
                       aes(x = coef, y = feature, color = value)) + 
                            geom_segment(aes(x = 0, xend = coef, 
                                             y = feature, yend = feature), 
                                         color = "black", alpha = 0.75) + 
                            geom_vline(xintercept = 0, alpha = 0.75) + 
                            geom_point(size = 3.5, alpha = 1) + 
                            geom_text(aes(label = qval_negat), 
                                      hjust = 1.3, 
                                      color = "black", 
                                      fontface = "italic", 
                                      size = 2.5) +
                            geom_text(aes(label = qval_posit), 
                                      hjust = -0.3, 
                                      color = "black", 
                                      fontface = "italic", 
                                      size = 2.5) +
                            facet_grid(Taxa_Level ~ value, 
                                       scales = "free_y", 
                                       space = "free") + 
                            theme_bw() +
                            theme(panel.grid.major.x = element_blank(),
                                  panel.grid.minor.x = element_blank(), 
                                  axis.text.y = element_text(face = "italic")) +
                            ylab("") + 
                            xlab("Coeficient") +
                            scale_color_manual(name = "Diet type", 
                                               values = aes.ls$gr_col) + 
                            scale_x_continuous(limits = c(-10.5, 10.5)) + 
                            theme(legend.text = element_text(size=12, 
                                                             face = "bold", 
                                                             family = "serif"), 
                                  legend.title = element_text(size=15))

ggsave(paste0(res.path, "/DA/plots/DA_res_min_prev", 
              gsub("\\.", "", min.prev), ".png"), 
       maas.plot.CT, dpi = 600, 
       width = 9.5, height = 8)


#-------------------------------------------------------------------------------
# Write results
res.ls[["DA"]][["res"]] <- maas.out

res.ls[["DA"]][["plot"]] <- maas.plot.CT

saveRDS(res.ls, paste0(res.path, "/supp/res_ls.RDS"))


#-------------------------------------------------------------------------------
# Significant taxa as boxplots 
#-------------------------------------------------------------------------------
box.maas <- list()

for (i.lvl in tax.lvl) {
  
  ps.maas <- data.ls$PS[[i.lvl]][["CSS"]]
  
  taxa_names(ps.maas) <- gsub("__", " ", taxa_names(ps.maas)) 
                          
    
  
  maas.res <- res.maas.long.f %>% 
                    filter(Taxa_Level == i.lvl)
    
  maas.res.all <- res.maas.long %>% 
                    filter(Taxa_Level == i.lvl) %>% 
                    mutate(feature = gsub("__", " ", feature)) 
    
    if (nrow(maas.res) > 0) {
      
      # Prepare data for plotting
      
      otu.inst <- ps.maas %>% 
                    prune_taxa(maas.res$feature, .) %>% 
                    otu_table() %>%
                    as.matrix() %>% 
                    as.data.frame()
      
      min.count <- min(otu.inst[otu.inst > 0])
      
      otu.inst <- log(otu.inst + min.count/2)
      
      da.taxa.df <- otu.inst %>% 
                        mutate(Taxa = rownames(.)) %>% 
                        gather(key = "SeqID", value = "Abundance", -Taxa) %>% 
                        left_join(., 
                                  mutate(meta, 
                                         SeqID = rownames(meta), 
                                         Level = i.lvl), 
                                  by = "SeqID") %>% 
                        dplyr::select(all_of(c("SeqID", "Abundance", 
                                               "Taxa", "Level", gr.col)))
      
      # Plot data
      #----------
      # Make data frame for significance levels 
      max.div <- da.taxa.df  %>% 
                    group_by(Taxa) %>% 
                    slice(which.max(Abundance)) 
      
      sig.df <- maas.res.all %>% 
                mutate(Taxa = feature, 
                       p.short = ifelse(round(qval, 2) == 0, "q>0.01", 
                                        paste0("q=", round(qval, 2))), 
                       Start = "CT", 
                       End = value) %>% 
                # filter(p.value <= 0.1) %>% 
                select(c("Taxa", "p.short", "Start", "End")) %>% 
                left_join(., max.div, by = "Taxa") %>% 
                group_by(Taxa) %>% 
                arrange(End, .by_group = TRUE) %>% 
                dplyr::filter(!is.na(DietID)) %>% 
                mutate(y.adj = Abundance*rep(c(1.1, 1.35, 1.55)))
      
      
      dif.abund.tax <- ggplot(da.taxa.df, 
                              aes_string(y = "Abundance", 
                                         x = "DietID")) + 
                      geom_jitter(aes_string(color = "DietID"), 
                                  size = 2, 
                                  alpha = 0.5, 
                                  width = 0.15, 
                                  height = 0) +
                      geom_boxplot(fill = "black",
                                   alpha = 0.1, 
                                   outlier.colour = NA) +
                      geom_signif(data = sig.df,
                                  aes(xmin = Start,
                                      xmax = End,
                                      annotations = p.short,
                                      y_position = y.adj),
                                  textsize = 4, 
                                  vjust = -0.1,
                                  manual = TRUE, margin_top = 1, ) +
                      geom_point(data = sig.df,
                                 aes(x = End, 
                                     y = y.adj*1.15), 
                                 x=NA) +
                      facet_wrap(c("Taxa"), 
                                 scales = "free", 
                                 ncol = 4) + 
                      theme_bw() + 
                      scale_color_manual(name = "Diet type", 
                                         values = aes.ls$gr_col) +
                      theme(legend.position = "none") + 
                      ylab("Abundance (CSS, log)")
                
                    
      
     ggsave(paste0(res.path, "/DA/plots/DA_box_", i.lvl, ".png"), 
            dif.abund.tax, 
            dpi = 600, 
            width = 10, 
            height = ceiling(length(unique(da.taxa.df$Taxa))/4)*2.5)              
    
  }
    
}


#-------------------------------------------------------------------------------
# Clean environment 
rm(list=ls())

gc()
