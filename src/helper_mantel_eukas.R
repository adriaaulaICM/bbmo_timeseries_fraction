library(tidyverse)
library(phyloseq)
library(vegan)

bl.phy <- readRDS("../data/robjects/blphyloseq.rds") 
source("sourcefiles/params-graphs.R") 
source("sourcefiles/main_functions.R")

geoMean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

bl.phy.relab <- transform_sample_counts(bl.phy, function(x) x / sum(x))
seltax <- filter_taxa(bl.phy.relab, function(x) max(x) >= 0.005, TRUE) %>% 
  taxa_names()

bl.phy <- subset_taxa(bl.phy, taxa_names(bl.phy) %in% seltax)

eukas.otus <- read_tsv('../data/abund-tax-raw/OTU_table99_PN_noMetNuclChar_DEF2.txt') 

colnames(eukas.otus) <- colnames(eukas.otus) %>% 
  str_replace_all(c(
    "BL060704" = "BL060705",
    "BL061010" = "BL061009",
    "BL080311" = "BL080312",
    "BL090521" = "BL090512",
    "BL100412" = "BL100414",
    "BL100413" = "BL100414",
    "BL110704" = "BL110705",
    "BL120518" = "BL120511",
    "BL131204" = "BL131215"))

physeq.list <- data_transformation(bl.phy)

bl.phy.relab <- physeq.list$relab
bl.phy.rar <- physeq.list$raref
bl.phy.log <- physeq.list$log
bl.phy.clr <- physeq.list$clr

# since the list is quite big let's erase it
physeq.list <- NULL


bl.phy.02 <- subset_samples(bl.phy,
                                fraction == '0.2') %>% 
  prune_taxa(taxa = taxa_sums(.) > 0, x = .)

bl.phy.3 <- subset_samples(bl.phy,
                               fraction == '3' ) %>% 
  prune_taxa(taxa = taxa_sums(.) > 0, x = .)

bl.phy.hell.02 <- transform_sample_counts(bl.phy.02, function(x) sqrt( x / sum(x)))
bl.phy.hell.3 <- transform_sample_counts(bl.phy.3, function(x) sqrt( x / sum(x)))

bl.phy.clr.02 <- subset_samples(bl.phy.clr,
                                fraction == '0.2') %>% 
  prune_taxa(taxa = taxa_sums(.) > 0, x = .)

bl.phy.clr.3 <- subset_samples(bl.phy.clr,
                               fraction == '3' ) %>% 
  prune_taxa(taxa = taxa_sums(.) > 0, x = .)

longeuka <- eukas.otus %>% 
  select(-CONSENS_GROUP, -CONSENS_SUPREGROUP) %>% 
  pivot_longer(names_to = 'sample',  values_to = 'count', cols = -OTUId) %>% 
  mutate(type = str_remove(sample, pattern = 'BL.*_'),
         sample = str_remove(sample, pattern = '_.*') )
  
keepab <- longeuka  %>% 
  group_by(sample, type) %>% 
  mutate( relab = count / sum(count)) %>% 
  group_by(OTUId, type) %>% 
  filter(max(relab) > 0.005) %>% 
  select(OTUId, type) %>% 
  distinct(OTUId, type) %>% 
  split(.$type) %>% 
  map(~ .x %>% pull(OTUId))


nano <- longeuka %>% 
  filter(type == 'N' & (OTUId %in% keepab$N),
         sample %in% sample_data(bl.phy.clr.02)$samname ) %>% 
  select(-type) %>% 
  group_by(sample) %>% 
  mutate( count = count + 1, 
          count = log(count / geoMean(count))) %>% 
  pivot_wider(names_from = OTUId, values_from = count, id_cols = sample) 

pico <- longeuka %>% 
  filter(type == 'P' & (OTUId %in% keepab$P),
         sample %in% sample_data(bl.phy.clr.3)$samname ) %>% 
  select(-type) %>% 
  group_by(sample) %>% 
  mutate( count = count + 1, 
          count = log(count / geoMean(count))) %>% 
  pivot_wider(names_from = OTUId, values_from = count, id_cols = sample)

  

# Mantels -----------------------------------------------------------------


comm.nano <- dist(nano %>% column_to_rownames(var = 'sample'),
                  method = 'euclidean')
  
comm.02 <- subset_samples(bl.phy.clr.02,
                          sample_data(bl.phy.clr.02)$samname %in% nano$sample) %>% 
  distance(., method = 'euclidean')


comm.pico <- dist(pico %>% column_to_rownames(var = 'sample'),
                  method = 'euclidean')
  
comm.3 <- subset_samples(bl.phy.clr.3,
                          sample_data(bl.phy.clr.3)$samname %in% pico$sample) %>% 
  distance(., method = 'euclidean')


mantel(comm.02, comm.nano)
mantel(comm.3, comm.pico)

comm.02.forpico <- subset_samples(bl.phy.clr.02,
                          sample_data(bl.phy.clr.02)$samname %in% pico$sample) %>% 
  distance(., method = 'euclidean')

comm.3.fornano <- subset_samples(bl.phy.clr.3,
                          sample_data(bl.phy.clr.3)$samname %in% nano$sample) %>% 
  distance(., method = 'euclidean')

mantel(comm.02.forpico, comm.pico)

mantel(comm.3.fornano, 
       dist(nano %>%
              filter(sample %in% sample_data(bl.phy.clr.3)$samname) %>%
              column_to_rownames(var = 'sample'),
                  method = 'euclidean'))  



# pCCA --------------------------------------------------------------------

nano <- longeuka %>% 
  filter(type == 'N' & (OTUId %in% keepab$N)) %>% 
  select(-type) %>% 
  group_by(sample) %>% 
  mutate( count = sqrt(count / sum(count))) %>% 
  pivot_wider(names_from = OTUId, values_from = count, id_cols = sample) 

pico <- longeuka %>% 
  filter(type == 'P' & (OTUId %in% keepab$P)) %>% 
  select(-type) %>% 
  group_by(sample) %>% 
  mutate( count = sqrt(count / sum(count))) %>% 
  pivot_wider(names_from = OTUId, values_from = count, id_cols = sample)

# here

yichun_pcca <- function(y.FL.tran, x.Euk.tran, envi){
  
  x.dca = cca(x.Euk.tran) # use DCA to decide if CCA or RDA is used
  test=summary(x.dca)
  axes=which(test$cont$importance["Cumulative Proportion",]>=0.7)[1]
  new.x=data.frame(x.dca$CA$u[,1:axes])
  
  y.full = cca(y.FL.tran~.,data=new.x) # cca/rda with full independent variables
  y.red = cca(y.FL.tran~1,data=new.x) # cca/rda with 1 independent variable
  y.rda = step(y.red, scope=list(lower=~1, upper=formula(y.full)))
  new.x=data.frame(x.dca$CA$u[,attr(y.rda$terms,"term.labels")])
  y.full = cca(y.FL.tran~.,data=new.x) # cca/rda with full independent variables
  # y.full.FL = varpart(y.FL.tran,new.x,envi,chisquare=T)
 
  print('I have arrived here') 
  y.full = cca(y.FL.tran,new.x,envi) 
  
  # results <- list( varp = y.full.FL,
  #                  cca = y.full)
  
  return(y.full)
  
  }

pcca_horror <- function( bacs, eukas){
  
  bac.raw <- psmelt_dplyr(bacs) %>% 
    as_tibble() %>% 
    mutate_at(c("Temperature", "Chla_total"), as.double) %>% 
    filter(!is.na(Temperature),
           !is.na(Chla_total),
           samname %in% eukas$sample)
  
  bac.table <- bac.raw %>%
    select(samname, Abundance, OTU) %>% 
    pivot_wider( names_from = 'OTU',
                 values_from = 'Abundance',
                 id_cols = 'samname') %>% 
    select(-samname)
  
  env.table <- bac.raw %>% 
    select(Temperature, Chla_total) %>% 
    distinct() %>% 
    mutate( Temperature = (Temperature - mean(Temperature)) / sd(Temperature), 
            Chla_total = (Chla_total - mean(Chla_total)) / sd(Chla_total)) 
  
  euka.table <- eukas %>% 
    filter(sample %in% bac.raw$samname) %>% 
    ungroup() %>% 
    select(-sample)
  
  ccaobj <- yichun_pcca( bac.table,
                         euka.table,
                         env.table)
  
  return(ccaobj)
  
}

cca02.pico <- pcca_horror(bacs = bl.phy.hell.02, eukas = pico)
anova(cca02.pico)
cca3.pico <- pcca_horror(bacs = bl.phy.hell.3, eukas = pico)
anova(cca3.pico)
cca02.nano <- pcca_horror(bacs = bl.phy.hell.02, eukas = nano)
anova(cca02.nano)
cca3.nano <- pcca_horror(bacs = bl.phy.hell.3, eukas = nano)
anova(cca3.nano)

cca02.pico
cca02.nano
cca3.pico
cca3.nano


  anova(cca02.pico)
  anova(cca02.nano)
  anova(cca3.pico)
  anova(cca3.nano)