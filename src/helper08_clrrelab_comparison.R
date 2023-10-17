library(tidyverse)
library(phyloseq)
library(gt)

bl.phy <- readRDS("../data/robjects/blphyloseq.rds") 
source("sourcefiles/params-graphs.R") 
source("sourcefiles/main_functions.R")

sums <- data.frame(asv = taxa_names(bl.phy.0.2),
                   as(tax_table(bl.phy.0.2),'matrix'),
                    sums = taxa_sums(bl.phy.0.2))

sums %>% 
  filter(curated == 'Cya_Synechococcus')  %>% 
  arrange(-sums) %>% 
  View()

bl.phy.0.2 <- subset_samples(bl.phy, fraction == '0.2') %>% 
  prune_taxa(taxa = taxa_sums(.) != 0, .)

physeq.list <- data_transformation(bl.phy.0.2)

bl.phy.relab <- physeq.list$relab
bl.phy.raref <- physeq.list$raref
bl.phy.log <- physeq.list$log
bl.phy.clr <- physeq.list$clr

library(ggforce)

synecos <-  subset_taxa(bl.phy.relab,curated == 'Cya_Synechococcus') 

data.frame(relative_abundance = sample_sums(synecos),
           sample_data(synecos))  %>% 
  select(decimal_date, relative_abundance, Synechococcus) %>% 
  rename( 'FC_synecos' = Synechococcus) %>% 
  gather( key = 'procedence', value = 'scaled.val', -decimal_date) %>% 
  group_by( procedence) %>% 
  mutate( scaled.val = scale(scaled.val)) %>% 
  ggplot( aes(decimal_date, scaled.val)) + 
  geom_line(aes(group = decimal_date), color = 'grey') + 
  geom_point(aes(color = procedence)) + 
  geom_line(aes(color = procedence)) + 
  ylab('Scaled value') + 
  leg.bottom  +
  thedark
  
data.frame(relative_abundance = sample_sums(synecos),
           sample_data(synecos))  %>% 
  mutate( Synechos_relative_fc = Synechococcus / Bacteria_joint) %>% 
  select(decimal_date, relative_abundance, Synechos_relative_fc) %>% 
  gather( key = 'procedence', value = 'scaled.val', -decimal_date) %>% 
  group_by( procedence) %>% 
  mutate( scaled.val = scale(scaled.val)) %>% 
  ggplot( aes(decimal_date, scaled.val)) + 
  geom_line(aes(group = decimal_date), color = 'grey') + 
  geom_point(aes(color = rev(procedence)), show.legend = F) + 
  geom_line(aes(color = rev(procedence)), show.legend = F) + 
  ylab('Scaled value') + 
  leg.bottom  +
  thedark

thing <- psmelt_dplyr(bl.phy.clr)

thing %>% 
  filter(OTU == 'asv1')  %>% 
  select(Date,OTU, Abundance, Synechococcus) %>%
  spread(key = OTU, value = Abundance) %>% 
  gather(key = 'sp', value = 'clr', -Date) %>% 
  mutate(cate = ifelse(str_detect(sp, 'asv'), 'clr transformed data',
                       'flow citometry data')) %>% 
  ggplot( aes(Date, clr)) + 
  geom_point() + 
  geom_line(aes(group = sp)) + 
  facet_wrap(~cate, scales = 'free', ncol = 1)
  