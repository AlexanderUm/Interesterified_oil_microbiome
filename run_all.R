#-------------------------------------------------------------------------------
# Libraries
#-------------------------------------------------------------------------------
lib.to.load <- c("tidyverse", "phyloseq", "broom", 
                 "vegan", "Maaslin2", "metagenomeSeq", 
                 "knitr", "ComplexHeatmap", "qiime2R",
                 "ggsignif", "FSA", "cowplot") 

for (i in lib.to.load) {library(i, character.only = TRUE)}


#-------------------------------------------------------------------------------
# Parameters for all analysis
#-------------------------------------------------------------------------------
prm.ls <- list("Data" = list("q_path" = "data/", 
                             "m_path" = "data/samples_data.csv", 
                             "res_path" = "output",
                             "min_read_tax" = 10, 
                             "tax_lvls" = c("ASV", "Genus", "Family", "Phylum"), 
                             "ref_gr" = "CT", 
                             "group_col" = "DietID", 
                             "seq_ids" = "SeqID"))

prm.ls[["Alpha"]] <- list("measures" = c("Observed", "Shannon", 
                                         "InvSimpson", "PhyloDiverity"),
                          "Tax_lvl" = "ASV", 
                          "Norm" = "Rare")

prm.ls[["Beta"]] <- list("dists" = c("Unweighted UniFrac" = "unifrac", 
                                     "Weighted UniFrac" = "wunifrac", 
                                     "Jaccard" = "jaccard", 
                                     "Bray-Curtis" = "bray"),
                         "form" = "DietID",
                         "Tax_lvl" = "ASV", 
                         "Norm" = "CSS", 
                         "nperm" = 999)

prm.ls[["DA"]] <- list("Tax_lvl" = c("ASV", "Genus", "Family"), 
                       "min.prev" = 0.5, 
                       "Norm" = "Raw")

prm.ls[["Corr"]] <- list("Tax_lvl" = "Genus", 
                         "min_prev" = 0.5, 
                         "est_plot" = 0.75,
                         "Norm" = "CSS")

prm.ls[["overview"]] <- list("Tax_lvl" = c("Phylum", "Family", "Genus", "ASV"), 
                             "Norm" = "Raw", 
                             "min_prev_plot" = 0.125)

#-------------------------------------------------------------------------------
# Write parameters 
#-------------------------------------------------------------------------------
# Create folder 
dir.create(paste0(prm.ls$Data$res_path, "/supp"), recursive = TRUE)

# Write parameters
save(list = c("prm.ls"), 
     file = paste0(prm.ls$Data$res_path, "/supp/prm.Rdata"))


#-------------------------------------------------------------------------------
# Run scripts 
#-------------------------------------------------------------------------------
source("0_data.R")

source("1_diversity.R")

source("2_beta.R")

source("3_DA.R")

source("4_corr.R")

source("5_overview.R")


#-------------------------------------------------------------------------------
# Clean environment 
rm(list=ls())

gc()
