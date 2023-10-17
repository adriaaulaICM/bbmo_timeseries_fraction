library(tidyverse)
library(phyloseq)
library(vegan)

bl.phy <- readRDS("../data/robjects/blphyloseq.rds") 
source("sourcefiles/params-graphs.R") 
source("sourcefiles/main_functions.R")

bl.phy.relab <- transform_sample_counts(bl.phy, function(x) x / sum(x))
seltax <- filter_taxa(bl.phy.relab, function(x) max(x) >= 0.005, TRUE) %>% 
  taxa_names()

bl.phy <- subset_taxa(bl.phy, taxa_names(bl.phy) %in% seltax)

eukas.otus <- read_tsv('../data/abund-tax-raw/OTU_table99_PN_noMetNuclChar_DEF2.txt') 

# Mantel data -------------------------------------------------------------
mantel.df <- readxl::read_excel('../data/metadata-raw/CHEMTAX_Blanes_2000_2014_ord.xlsx') %>% 
  select(`Decimal year`, Month, Any, `Day of month`, `PRASINOPHYTES nanog L-1`,
         `DINOPHYTES nanog L-1`, `CRYPTOPHYTES nanog L-1`,
         `PELAGOPHYTES nanog L-1`, `Synechococcus nanog L-1`,
         `DIATOMS nanog L-1`, `Prochlorococcus nanog L-1`,
         `HAPTOPHYTES6+7+8 nanog L-1`) %>% 
  rename( year = Any, month = Month, day = `Day of month`)

bl.data <- sample_data(bl.phy) %>% 
  as(., 'matrix') %>% 
  as_tibble() %>% 
  select(year, month, day, fraction, samname) %>% 
  filter(fraction == '0.2') %>% 
  mutate( year = as.double(year), 
          month = as.double(month), 
          day = as.double(day))


selection.for.mantel <- mantel.df %>% 
  left_join(bl.data, by = c('year', 'month', 'day')) %>% 
  filter(!is.na(samname))


physeq.list <- data_transformation(bl.phy)

bl.phy.relab <- physeq.list$relab
bl.phy.rar <- physeq.list$raref
bl.phy.log <- physeq.list$log
bl.phy.clr <- physeq.list$clr

# since the list is quite big let's erase it
physeq.list <- NULL

bl.phy.02 <- subset_samples(bl.phy,
                                fraction == '0.2' &
                                samname %in% selection.for.mantel$samname) %>% 
  prune_taxa(taxa = taxa_sums(.) > 0, x = .)

bl.phy.3 <- subset_samples(bl.phy,
                               fraction == '3' &
                                samname %in% selection.for.mantel$samname) %>% 
  prune_taxa(taxa = taxa_sums(.) > 0, x = .)


bl.phy.rar.02 <- subset_samples(bl.phy.rar,
                                fraction == '0.2' &
                                samname %in% selection.for.mantel$samname) %>% 
  prune_taxa(taxa = taxa_sums(.) > 0, x = .)

bl.phy.rar.3 <- subset_samples(bl.phy.rar,
                               fraction == '3' &
                                samname %in% selection.for.mantel$samname) %>% 
  prune_taxa(taxa = taxa_sums(.) > 0, x = .)

bl.phy.clr.02 <- subset_samples(bl.phy.clr,
                                fraction == '0.2' &
                                samname %in% selection.for.mantel$samname) %>% 
  prune_taxa(taxa = taxa_sums(.) > 0, x = .)

bl.phy.clr.3 <- subset_samples(bl.phy.clr,
                               fraction == '3' &
                                samname %in% selection.for.mantel$samname) %>% 
  prune_taxa(taxa = taxa_sums(.) > 0, x = .)


# Calculating the distance ------------------------------------------------

# for phyto data
phyto.02.raw <- selection.for.mantel %>% 
  select(samname, `PRASINOPHYTES nanog L-1`: `HAPTOPHYTES6+7+8 nanog L-1`) %>% 
  column_to_rownames(var = 'samname') 

phyto.02 <- phyto.02.raw %>% 
  scale() %>% 
  dist(method = 'euclidean')


# the communities
comm.02 <- distance(bl.phy.clr.02, method = 'euclidean')

mantel(phyto.02, comm.02)

comm.3 <- distance(bl.phy.clr.3, method = 'euclidean')  %>% 
  as.matrix()
  
phyto3 <- phyto.02 %>% 
  as.matrix() %>% 
  .[sample_data(bl.phy.clr.3)$samname, sample_data(bl.phy.clr.3)$samname] %>% 
  as.dist()

mantel(phyto3, comm.3)

# Same with pCCA ----------------------------------------------------------

env02 <- sample_data(bl.phy.clr.02) %>% 
  as(., 'matrix') %>% 
  as_tibble() %>% 
  select(Temperature, Chla_total) %>% 
  mutate_all(., as.double)  %>% 
  mutate_all(., scale)


env3 <- sample_data(bl.phy.clr.3) %>% 
  as(., 'matrix') %>% 
  as_tibble() %>% 
  select(Temperature, Chla_total) %>% 
  mutate_all(., as.double)  %>% 
  mutate_all(., scale)


# New approach ------------------------------------------------------------

y.02.tran <- decostand(as(otu_table(bl.phy.02), 'matrix'),"hellinger")
x.phyto <- decostand(phyto.02.raw, 'hellinger')


cca02 <- cca(y.02.tran, x.phyto, env02)

anova(cca02)


y.3.tran <- decostand(as(otu_table(bl.phy.3), 'matrix'),"hellinger")
x.phy3 <- decostand(phyto.02.raw %>%
                       .[sample_data(bl.phy.clr.3)$samname,],
                     'hellinger')

cca3 <- cca(y.3.tran, x.phy3, env3)
anova(cca3)

# Adding the steps that Yi Chun presented 


yichun_pcca <- function(y.FL.tran, x.Euk.tran, envi){
  
  x.dca = cca(x.Euk.tran) # use DCA to decide if CCA or RDA is used
  test=summary(x.dca)
  axes=which(test$cont$importance["Cumulative Proportion",]>0.7)[1]
  new.x=data.frame(x.dca$CA$u[,1:axes])
  
  y.full = cca(y.FL.tran~.,data=new.x) # cca/rda with full independent variables
  y.red = cca(y.FL.tran~1,data=new.x) # cca/rda with 1 independent variable
  y.rda = step(y.red,scope=list(lower=~1,upper=formula(y.full)))
  new.x=data.frame(x.dca$CA$u[,attr(y.rda$terms,"term.labels")])
  y.full = cca(y.FL.tran~.,data=new.x) # cca/rda with full independent variables
  y.full.FL = varpart(y.FL.tran,new.x,envi,chisquare=T)
  
  y.full = cca(y.FL.tran,new.x,envi) 
  
  results <- list( varp = y.full.FL,
                   cca = y.full)
  
  return(results)
  
  }


pcca02phy <- yichun_pcca(y.FL.tran = y.02.tran,
                         x.Euk.tran = x.phyto,
                         envi = env02)

anova(pcca02phy$cca)

pcca3phy <- yichun_pcca(y.FL.tran = y.3.tran,
                        x.Euk.tran = x.phy3,
                        envi = env3)
anova(pcca3phy$cca)




