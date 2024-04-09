################################################################################
# Taxa and SCFA correlations
#-------------------------------------------------------------------------------
set.seed(435975)


#-------------------------------------------------------------------------------
# Load data
load("output/supp/data.Rdata")
res.ls <- readRDS("output/supp/res_ls.RDS")

# Custom functions
source("R/otu_tabs_extract.R")
source("R/bulk_correlation.R")
source("R/plot_cor_heat.R")

#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
ps.inst <- data.ls$PS[[prm.ls$Corr$Tax_lvl]][[prm.ls$Corr$Norm]]

min.prev <- prm.ls$Corr$min_prev

est.plot <- prm.ls$Corr$est_plot

meta <- data.ls$meta

ref.gr <- prm.ls$Data$ref_gr

gr.col <- prm.ls$Data$group_col

res.path <- prm.ls$Data$res_path

# Out folders
dir.create(paste0(res.path, "/Correlations/tabs"), recursive = TRUE)
dir.create(paste0(res.path, "/Correlations/plots"), recursive = TRUE)


#-------------------------------------------------------------------------------
# Data preparation 
#-------------------------------------------------------------------------------
# Read in csv with SCFA
scfa.data <- read.csv("data/SCFA_expr3.csv") %>%
                dplyr::arrange(match(SampleID, meta$SampleID))

scfa.data.cor <- scfa.data[, 4:6] %>% 
                    mutate(rowID = rownames(meta)) %>% 
                    column_to_rownames("rowID") %>% 
                    rename("Acetate" = "Acetic.Acid", 
                           "Propionate" = "Propionic.acid", 
                           "Butyrate" = "Butyric.acid")


# Extract otu tables for correlations 
otu.tabs.ls <- otu_tab_extract(list("Genus" = ps.inst))


cor.all.res <- bulk_corr_matrix(otu_tabs_list = otu.tabs.ls, 
                                second_df = scfa.data.cor, 
                                meta_data = meta, 
                                strata_cols_meta = gr.col, 
                                prev_cut_offs = min.prev)




#-------------------------------------------------------------------------------
# Plot heatmaps
#-------------------------------------------------------------------------------
heat.d <- cor.all.res[[1]] %>% 
                bind_rows() %>% 
                mutate(Taxa = gsub("__", " ", Taxa)) 
  
heat.p <- plot_cor_heat(cor_function_results = heat.d, 
                          estimate_cutoff = est.plot, 
                          qval_cutoff = 1, 
                          pval_cutoff = 1, 
                          group_col = "Strata", 
                          group_col_levels = c("CT", "ICT", "HF", "IHF"), 
                          estim_display_cut = (est.plot-0.05)) 


png(paste0(res.path, "/Correlations/plots/heatmap.png"), 
      width = 22, height = (length(heat.p$taxa)/2 + 5), res = 600, units = "cm")
plot(draw(heat.p$plot))
dev.off()  


#-------------------------------------------------------------------------------
# Heatmap summary 
#-------------------------------------------------------------------------------
heat.d.sum <- heat.d %>% 
                select(Taxa, p.value, estimate, qval, Strata, Corre_vector) %>% 
                pivot_wider(names_from = c("Strata", "Corre_vector"),
                            values_from = c("estimate", "p.value", "qval"), 
                            names_glue = "{Corre_vector}_{Strata}_{.value}") 

write.csv(heat.d.sum, paste0(res.path, "/Correlations/tabs/corr_summary.csv"))


#-------------------------------------------------------------------------------
# Correlations scatter plots 
#-------------------------------------------------------------------------------

scfa.data.cor.p <- scfa.data.cor %>% 
                    mutate(SeqID = rownames(.))

cor.all.long.df <- cor.all.res[[1]] %>% 
                        bind_rows() %>% 
                        mutate(DietID = Strata, 
                               Label = paste0(ifelse(abs(estimate) < 0.01, 
                                                     "est<0.01", 
                                                     paste0("est=", round(estimate, 2)))))

otu.inst.cor <- otu.tabs.ls$Genus


for(i.met in unique(heat.d$Corre_vector)) {
  
  cor.d.inst <- cor.all.long.df %>% 
                  filter(Corre_vector == i.met, 
                         abs(estimate) >= est.plot) %>% 
                  slice_max(abs(estimate), by = Taxa, with_ties = FALSE) %>% 
                  arrange(desc(abs(estimate)))  
  
  otu.long.inst <- otu.inst.cor %>% 
                    select(cor.d.inst$Taxa) %>% 
                    mutate(SeqID = rownames(.)) %>% 
                    pivot_longer(., 
                                 cols = cor.d.inst$Taxa,
                                 names_to = "Taxa") %>% 
                    left_join(., scfa.data.cor.p, by = "SeqID") %>%  
                    left_join(., meta, by = "SeqID") %>% 
                    mutate(Taxa = factor(Taxa, levels = cor.d.inst$Taxa)) %>% 
                    arrange(across(all_of(c("Taxa", gr.col))))

  
  labl.tab.y <- otu.long.inst %>% 
            reframe(labl.y = max(value), 
                    .by = c("Taxa")) 
  
  labl.tab.x <- otu.long.inst %>% 
                  reframe(labl.x = min(.data[[i.met]], na.rm = T), 
                          .by = c(gr.col)) 
    
  cor.labs <- cor.all.long.df %>% 
                    filter(Taxa %in% cor.d.inst$Taxa, 
                           Corre_vector == i.met) %>% 
                    arrange(across(all_of(c("Taxa", gr.col)))) %>% 
                                  left_join(., labl.tab.y, by = "Taxa") %>% 
                                  left_join(., labl.tab.x, by = gr.col) %>% 
                    mutate(Taxa = gsub("__", " ", Taxa))
  
  otu.long.inst <- otu.long.inst %>% 
                        mutate(Taxa = gsub("__", " ", Taxa))
  
  p.cor <- ggplot(otu.long.inst, 
                       aes(x = .data[[i.met]], 
                           y = value, 
                           color = .data[[gr.col]])) + 
                  geom_point(size = 2) + 
                  geom_smooth(color = "gray25", 
                              method = "lm", 
                              se = FALSE, 
                              size = 0.5) + 
                  geom_text(data = cor.labs,
                             aes(y = labl.y,
                                 x = labl.x,
                                 label = Label),
                             color = "black",
                             vjust = 0.75,
                             hjust = -0.01, alpha = 0.5) +
                  scale_color_manual(values = aes.ls$gr_col) +
                  facet_grid(c("Taxa", gr.col), scales = "free") + 
                  theme_classic() + 
                  theme(panel.border = element_rect(colour = "black", 
                                                    fill=NA, 
                                                    size=0.5), 
                        strip.text.y = element_text(size = 6, face = "italic"))
     
  ggsave(paste0(res.path, "/Correlations/plots/cor_scatter_", i.met, ".png"),
         plot = p.cor,
         width = length(unique(otu.long.inst[[gr.col]]))*2.5, 
         height = length(cor.d.inst$Taxa)*1.5, limitsize = FALSE)       
  
}

#-------------------------------------------------------------------------------
# Clean environment 
rm(list=ls())

gc()


