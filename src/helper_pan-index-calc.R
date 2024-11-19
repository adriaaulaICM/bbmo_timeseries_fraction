library(tidyverse)
library(phyloseq)
# PAN index calculation script, see https://github.com/GuillemSalazar/EcolUtils
library(EcolUtils)
library(vegan)

readRDS(phylo.rem,file = '../data/robjects/blphy10years.rds')

phylo.raref <- rarefy_even_depth(phylo.rem)

abund.raref <- as(otu_table(phylo.raref), "matrix") %>% as.data.frame() 
tax.raref <- as(tax_table(phylo.raref), 'matrix') %>%
  data.frame() %>% 
  rownames_to_column(var = "asv")


fraction <- decostand(as.numeric(sample_data(phylo.raref)$fraction),
                      "range")

nich.df <- niche.val(abund.raref, fraction, n = 1000)
# write_rds(nich.df, file = '../data/robjects/nicheval.rds')

pan <- nich.df %>% 
  rownames_to_column(var = "asv") %>%
  dplyr::rename( pan.index = observed) 

tax.ts <- tax.raref %>% 
  left_join(pan, by = "asv") %>% 
  mutate(taxsums = taxa_sums(phylo.raref))

pan.df <- tax.ts %>% 
  filter(taxsums > 10) %>% 
  filter(sign != "NON SIGNIFICANT") %>% 
  select(asv, pan.index)
  
tax.new <- data.frame(as(tax_table(phylo.rem), 'matrix')) %>%
  rownames_to_column(var = "asv") %>% 
  left_join(pan.df, by = "asv") %>% 
  column_to_rownames(var = "asv") %>% 
  as.matrix() %>% 
  tax_table()
  
phylo.rem.pan <- merge_phyloseq(phylo.rem, tax.new)

write_rds(phylo.rem.pan, path = "../data/robjects/blphyloseq.rds")
