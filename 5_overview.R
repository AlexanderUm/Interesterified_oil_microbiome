################################################################################
# General composition overview  
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
tax.lvl <- prm.ls$overview$Tax_lvl

count.norm <- prm.ls$overview$Norm

meta <- data.ls$meta

ref.gr <- prm.ls$Data$ref_gr

gr.col <- prm.ls$Data$group_col

res.path <- prm.ls$Data$res_path

min.prev <- prm.ls$overview$min_prev_plot

# Out folders
dir.create(paste0(res.path, "/overview/tabs"), recursive = TRUE)
dir.create(paste0(res.path, "/overview/plots"), recursive = TRUE)

#-------------------------------------------------------------------------------
# Tables and heatmaps overviews
#-------------------------------------------------------------------------------
# Prepare data
col.split.v <- meta[[gr.col]] %>% 
                  setNames(rownames(meta))


for(i.lvl in tax.lvl) {
  
  ps.inst <- data.ls$PS[[i.lvl]][[count.norm]]
  
  # Extract OTU table
  otu.inst <- ps.inst %>% 
                phy_taxa_filter(., prev = 0.125) %>% 
                transform_sample_counts(., function(x) x/sum(x)*100) %>% 
                otu_table(.) %>% 
                as.matrix() %>% 
                as.data.frame()
  
  # Adjust row names 
  rownames(otu.inst) <- gsub("__", " ", rownames(otu.inst))
  
  #-----------------------------------------------------------------------------
  # Summary table 
  #-----------------------------------------------------------------------------
  otu.inst <- otu.inst %>% as.matrix()
  
  sum.tab <- data.frame(Mean = rowMeans(otu.inst), 
                        Median = rowMedians(otu.inst), 
                        Min = rowMin(otu.inst), 
                        Max = rowMax(otu.inst), 
                        SD = apply(otu.inst, 1, sd, simplify = TRUE), 
                        Prevalence = apply(otu.inst, 1, 
                                           function(x){length(x[x>0])/length(x)}, 
                                           simplify = TRUE)) %>% 
              arrange(desc(Mean))
  
  write.csv(sum.tab, 
            paste0(res.path, "/overview/tabs/overview_", i.lvl, ".csv"))
  
  
  #-------------------------------------------------------------------------------
  # Taxa table
  #-------------------------------------------------------------------------------
  ps.inst %>% 
    tax_table() %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    'rownames<-'(., gsub("__", " ", rownames(.))) %>% 
    write.csv(., paste0(res.path, "/overview/tabs/taxa_table_", i.lvl, ".csv"))
  
  
  #-----------------------------------------------------------------------------
  # Plot heatmap 
  #-----------------------------------------------------------------------------
  # Sort table by mean
  otu.inst <- otu.inst[rownames(sum.tab), ]
  
  # Plot heatmap
  overv.heat <- Heatmap(otu.inst, 
                        column_split = col.split.v, 
                        show_column_names = FALSE, 
                        cluster_column_slices = FALSE, 
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE,
                        column_gap = unit(0.2, "mm"),
                        name = "Abundance (%)", 
                        rect_gp = gpar(col = "gray", lwd = 0.5),
                        border = TRUE, 
                        row_names_max_width = max_text_width(rownames(otu.inst)))
  
  # Write plot  
  png(paste0(res.path, 
             "/overview/plots/heat_", i.lvl, ".png"), 
      width = 24, height = (nrow(otu.inst)/2 + 4), res = 600, units = "cm")
  plot(overv.heat)
  dev.off()
  
}

#-------------------------------------------------------------------------------
# Clean environment 
rm(list=ls())

gc()

  
