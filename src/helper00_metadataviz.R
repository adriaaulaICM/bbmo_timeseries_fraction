library(tidyverse)
library(scales) 
library(phyloseq)

source('sourcefiles/params-graphs.R')
source('sourcefiles/backbone_params-graphs.R')

metadata <- read_tsv('../data/metadata-raw/metadata_Blanes_compact_may2017.tsv')

labels <- readxl::read_xlsx('../data/metadata-raw/labels_pretty.xlsx',
                            sheet = 2,
                            col_names = c('Variable', 'Measure'))  %>% 
  mutate(all = str_c(Variable, "~", Measure )) %>% 
  mutate(all = ifelse(Variable == 'Salinity', 'Salinity', all))

par.order <-  c("Day_length", "Temperature", "Salinity", "Secchi",
    "NO2", "NO3", "PO4", "Si", "Chla_total", "BP_FC1.55",
    "Bacteria_joint", "HNA", "LNA", "Synechococcus", "Prochlorococcus_FC",
    "Peuk1", "Peuk2", "PNF_Micro", "PNF2_5um_Micro", "PNF_5um_Micro",
    "Cryptomonas", "Micromonas", "HNF_Micro")

# Filtering for specific years for 16S!
metadata <- metadata %>% filter(year > 2003, year < 2010)

# Plotting env and bio separately 
env <- metadata %>% 
  mutate(Temperature = as.numeric(str_replace(Temperature, ",", "."))) %>% 
  select(year, day_of_year, season, month, one_of(par.order[1:8])) %>% 
  gather(key = 'parameter', value = 'val',
         -month, -year, -day_of_year, -season,
         factor_key = T) %>%
  # Outlier value the 2008 , at may
  filter(!(year == 2008 & month == '05')) %>% 
  mutate( parameter = factor(parameter, levels = par.order,
                             labels = labels$all)) 

bio <- metadata %>% 
  mutate(Temperature = as.numeric(str_replace(Temperature, ",", "."))) %>% 
  select(year, month, season, day_of_year, one_of(par.order[9:23])) %>% 
  gather(key = 'parameter', value = 'val',
         -month, -year, -day_of_year, -season,
         factor_key = T) %>%
  # Outlier value the 2008 , at may
  filter(!(year == 2008 & month == '05')) %>% 
  mutate( parameter = factor(parameter, levels = par.order,
                             labels = labels$all)) 


means <- env %>% 
  group_by(parameter, month) %>% 
  summarize(total = mean(val, na.rm = T),
            sd = sd(val, na.rm = T) ) %>% 
  mutate( day_of_year = seq(15, 365, by = 30))


options(scipen=3)

env  %>% 
  ggplot( aes(day_of_year, val)) + 
  geom_point(color = 'gray')  +
  stat_smooth(aes(group = parameter),
              method = "gam",
              color = 'gray48',
              formula = y ~ s(x, k =12, bs = 'cc'),
              se = F, size = 1,
              show.legend = F, alpha = 0.7) + 
  geom_pointrange(data = means,
                  aes( y = total, 
                       ymin = total - sd, ymax = total + sd),
                  color = 'cornflowerblue',
                  size = 0.3) + 
  facet_wrap(~parameter, ncol = 2, 
             labeller = label_parsed,
             scales = 'free_y') + 
  ylab(NULL) +
  xlab('Month') +
  scale_dayyear_shrt + 
  lil.strip + 
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 8))

ggsave('../results/figures/metadata_env_6years.pdf',
       width = 8, 
       height = 9)
  

bio.filt <-   bio  %>% 
  filter(!parameter %in% c('PNF~(5~mu*m)~(cells~ml^{-1})',
                          'Cryptomonads~(cells~ml^{-1})',
                          'italic(Micromonas)-like~(cells~ml^{-1})'))
  
means.bio <- bio.filt %>% 
  group_by(parameter, month) %>% 
  summarize(total = mean(val, na.rm = T),
            sd = sd(val, na.rm = T) ) %>% 
  mutate( day_of_year = seq(15, 365, by = 30))

options(scipen=1)
 
bio.filt  %>% 
  filter(!parameter %in% c('PNF~(5~mu*m)~(cells~ml^{-1})',
                          'Cryptomonads~(cells~ml^{-1})',
                          'italic(Micromonas)-like~(cells~ml^{-1})')) %>% 
  ggplot( aes(day_of_year, val)) + 
  geom_point(color = 'gray')  +
  stat_smooth(aes(group = parameter),
              method = "gam",
              color = 'gray48',
              formula = y ~ s(x, k =12, bs = 'cc'),
              se = F, size = 1,
              show.legend = F, alpha = 0.7) + 
  geom_pointrange(data = means.bio,
                  aes( y = total, 
                       ymin = total - sd, ymax = total + sd), 
                  color = 'cornflowerblue',
                  size = 0.3) + 
  facet_wrap(~parameter, ncol = 3, 
             labeller = label_parsed,
             scales = 'free_y')  + 
  ylab(NULL) +
  xlab('Month') +
  scale_y_continuous(labels = function(x) format(x, scientific = -3)) + 
  scale_dayyear_shrt + 
  lil.strip + 
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 8))
 

ggsave('../results/figures/metadata_bio_6years.pdf',
       width = 8, 
       height = 9)
