set.seed(3409576)

#-------------------------------------------------------------------------------
# Load parameters 
load("output/supp/prm.Rdata")

# Custom functions
source("R/phy_shorten_tax_names.R")
source("R/phy_css_norm.R")

# List for results
data.ls <- list(PS = list())
res.ls <- list()

#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
q.path <- prm.ls$Data$q_path

meta.path <- prm.ls$Data$m_path

tax.lvls <- prm.ls$Data$tax_lvls

min.r <- prm.ls$Data$min_read_tax

ref.gr <- prm.ls$Data$ref_gr

gr.col <- prm.ls$Data$group_col

seq.id <- prm.ls$Data$seq_ids


#-------------------------------------------------------------------------------
# Aesthetics 
#-------------------------------------------------------------------------------
aes.ls <- list(gr_col = c(HF="#377EB8",  CT="#4DAF4A", 
                          ICT="#984EA3", IHF="#A65628"))


#------------------------------------------------------------------------------
# Data import 
#------------------------------------------------------------------------------
# Import data into phyloseq 
ps <- qza_to_phyloseq(features = paste0(q.path, "asv_table.qza"), 
                      tree = paste0(q.path, "tree/rooted-tree.qza"), 
                      taxonomy = paste0(q.path, "taxonomy_07.qza"))


ps.meta <- read.csv(meta.path, row.names = "X")  %>% 
                    mutate(!!gr.col := factor(get(gr.col), 
                                       levels = c("CT", "ICT", "HF", "IHF"))) %>% 
                    .[sample_names(ps), ]
                    
# Add as metadata
sample_data(ps) <- ps.meta
  
data.ls[["meta"]] <- ps.meta

#-------------------------------------------------------------------------------
# Filter out taxa
#-------------------------------------------------------------------------------
# filter taxa with with less than X reads in total   
ps <- prune_taxa(taxa_sums(ps) >= min.r, ps)

# Remove ASVs: 
ps <- prune_taxa(!tax_table(ps)[, "Genus"] %in% "Mitochondria", ps)

ps <- prune_taxa(!tax_table(ps)[, "Genus"] %in% "Chloroplast", ps)

ps <- prune_taxa(!is.na(tax_table(ps)[, "Phylum"])[, "Phylum"], ps)

ps <- prune_taxa(tax_table(ps)[, "Kingdom"] %in% c("d__Bacteria", 
                                                   "d__Archaea"), ps)


#-------------------------------------------------------------------------------
# Glom to higher taxonomic level & Normalize count 
#-------------------------------------------------------------------------------
for(i.lvl in tax.lvls)  {
  
  if(i.lvl == "ASV") { ps.inst <- ps
  
  } else{ ps.inst <- tax_glom(ps, i.lvl)}
  
  # Adjust taxa names
  taxa_names(ps.inst) <- phy_shorten_tax_names(ps.inst) %>% 
                         make.unique()
  
  data.ls$PS[[i.lvl]] <- list("Raw" = ps.inst, 
                              "Rare" = rarefy_even_depth(ps.inst, rngseed = 347), 
                              "CSS" = phy_css_norm(ps.inst))
  
}


#-------------------------------------------------------------------------------
# Save data 
#-------------------------------------------------------------------------------
save(list = c("data.ls", "prm.ls", "aes.ls"), 
     file = paste0(prm.ls$Data$res_path, "/supp/data.Rdata"))


# Clean environment 
rm(list=ls())

gc()
