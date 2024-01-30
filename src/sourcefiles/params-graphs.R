
#### General #####

library(ggthemes)
library(tidyverse)

#Black and white theme and the text size bigger
### If the work is a presentation, its best to give a base size of...around 13
theme_set(theme_bw(base_size = 13))

# Lets also set the seed for random processes
set.seed(42)


#### Palettes phylogeny (BACTERIA) #####

# Old one 
# superpalette <- c("yellow1", #B_Flavobacteriia
#                   # "#ffce00", #B_Sphingobacteriia
#                   "#ff9a00",  #Other_Bacteroidetes
#                   "limegreen", #Cya_Synechococcus
#                   "springgreen4", #Other_Cyanobacteria
#                   "lightcyan2", #P_Alp_Rhodobacterales
#                   # "#78aaff", #P_Alp_Rhodospirillales
#                   "steelblue1", #P_Alp_Rickettsiales
#                   "#4188ff", #"P_Alp_SAR11_clade
#                   "steelblue4", #P_Alp_Sphingomonadales
#                   "royalblue4", #Alp_Proteobacteria
#                   "#efbbff", #P_Gam_Alteromonadales
#                   "#d896ff", #P_Gam_Cellvibrionales
#                   "#b76bce", #P_Gam_Oceanospirillales
#                   # "#be29ec", #P_Gam_Pseudomonadales
#                   "#660066",  #Gam_Proteobacteria
#                   "#330033", #Other_Proteobacteria
#                   "#f7cfcf", #O_Actinobacteria
#                   "indianred1", #O_Firmicutes
#                   "indianred", #O_Planctomycetes
#                   "#a03e3e", #"O_Verrucomicrobia
#                   "#492a2a") #"Other_Bacteria

# new, from Isabel
superpalette <- 	c("#FDAE61FF", #B_Flavobacteriia
                   "#FEE090FF",  #Other_Bacteroidetes
                   "#65A479FF", #Cya_Synechococcus
                   "#92C5DEFF", #P_Alp_Rhodobacterales
                   "#ABD9E9FF", #P_Alp_Rickettsiales
                   "#4393C3FF", #P_Alp_SAR11_clade
                   "#2166ACFF", # P_Alp_Sphingomonadales
                   "#053061FF", #Alp_Proteobacteria
                   "#EAA9BDFF", #P_Gam_Alteromonadales
                   "#CA699DFF",  #P_Gam_Cellvibrionales
                   "#762A83FF", #P_Gam_Oceanospirilalles
                   "#C2A5CFFF", #Gam_Proteobacteria
                   # "#330033", #Other Proteobacteria
                   "#D6604DFF", #Actinobacteria
                   "#B2182BFF", #O_Firmicutes
                   "#BABABAFF", #O_Planctomycetes
                   "#5C6068FF", #O_Verrucomicrobia
                   "#000000FF") #Other Bacteria

tax.order <-c("Flavobacteriales",
              "Other Bacteroidota",
              "Cyanobiaceae",
              "Alp-Rhodobacterales",
              "Alp-Rickettsiales",
              "Alp-SAR11 clade",
              "Alp-Sphingomonadales",
              "Other Alphaproteobacteria",
              "Gam-Alteromonadaceae",
              "Gam-Cellvibrionales",
              "Gam-Oceanospirillales",
              # "P_Gam_Pseudomonadales",
              "Other Gammaproteobacteria",
              # "Other Proteobacteria",
              "Actinobacteriota",
              "Firmicutes",
              "Planctomycetota",
              "Verrucomicrobiota",
              "Other Bacteria")

tax.order.palette <-c("Flavobacteriales",
              "Other Bacteroidota",
              "Cyanobiaceae",
              "Alp-Rhodobacterales",
              "Alp-Rickettsiales",
              "Alp-SAR11 clade",
              "Alp-Sphingomonadales",
              "Other Alphaproteobacteria",
              "Gam-Alteromonadaceae",
              "Gam-Cellvibrionales",
              "Gam-Oceanospirillales",
              # "P_Gam_Pseudomonadales",
              "Other Gammaproteobacteria",
              # "Other Proteobacteria",
              "Actinobacteriota",
              "Firmicutes",
              "Planctomycetota",
              "Verrucomicrobiota",
              "Other Bacteria")

names(superpalette) <- tax.order.palette


bac.colScale <- scale_color_manual(name = "16S taxonomy",
                                   values = superpalette,
                                   limits = names(superpalette))
bac.fillScale <- scale_fill_manual(name = "16S taxonomy",
                                   values = superpalette,
                                   limits = names(superpalette))

italics.strip.legend <-  theme(legend.text = ggtext::element_markdown(),
                               strip.text.x = ggtext::element_markdown()  ) 

italizyce_synes <- function(df){
  
  df %>% 
    mutate(curated = ifelse(curated == 'Synechococcus',
                            '*Synechococcus*',
                            curated)) %>% 
    mutate(curated = factor(curated, levels = tax.order.palette)) %>% 
    return()
  
}

pretty.legend <- function(gg){
  
  gg + 
    theme( legend.box.margin = margin(0,0,0,0),
           legend.margin = margin(0,0,0,0),
           legend.key.height = unit(0, 'cm'), 
           legend.key.width = unit(0, 'cm')) + 
    guides( fill = guide_legend(ncol = 1,
                                override.aes = list(size = 3)))  
           
}
### Palettes EUKARYA #####

superpalette.euk <- c("olivedrab3", #Alveolata/Dinoflagellata
                      "springgreen4", #Alveolata/MALV-I
                      "yellow",  #Alveolata/MALV-II
                      "#ffce00", #Other MALVs
                      "#ff9a00", #Other Alveolates
                      "#ff4d00", #Amoebozoa
                      "#78aaff" , #Archaeplastida/Chlorodendrophyceae
                      "steelblue1" , #Archaeplastida/Mamiellophyceae
                      "steelblue4" , #Archaeplastida/Prasinophyceae
                      "royalblue4" , #Other Archaeplastida
                      "#b76bce", #Hacrobia/Cryptomonadales
                      "#efbbff", #Hacrobia/Picozoa
                      "#be29ec", #Other Hacrobia
                      "#660066", #Opisthokonta
                      "indianred", #Rhizaria
                      "#a03e3e",  #Stramenopiles/Diatomea
                      "#492a2a", #Stramenopiles/Labyrinthulomycetes
                      "indianred1", #Other Stramenopiles
                      "#666463"  #Other Eukaryotes 
)

tax.order.euk <- c(   
  "Alveolata/Dinoflagellata",
  "Alveolata/MALV-I",
  "Alveolata/MALV-II",
  "Other MALVs",
  "Other Alveolates",
  "Amoebozoa",
  "Archaeplastida/Chlorodendrophyceae",
  "Archaeplastida/Mamiellophyceae",
  "Archaeplastida/Prasinophyceae",
  "Other Archaeplastida",
  "Hacrobia/Cryptomonadales",
  "Hacrobia/Picozoa",
  "Other Hacrobia",
  "Opisthokonta",
  "Rhizaria",
  "Stramenopiles/Diatomea",
  "Stramenopiles/Labyrinthulomycetes",
  "Other Stramenopiles",
  "Other Eukaryotes")
  
names(superpalette.euk) <- tax.order.euk

euk.colScale <- scale_color_manual(name = "18S taxonomy",
                                   values = superpalette.euk,
                                   limits = names(superpalette.euk))
euk.fillScale <- scale_fill_manual(name = "18S taxonomy",
                                   values = superpalette.euk,
                                   limits = names(superpalette.euk))

###### VARIOUS AESTHETICS #####

##### GGPlot ####

# Label related 
leg.bottom <- theme(legend.direction="horizontal",legend.position="bottom")
lab.flip <- theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
labelsx.diag <- theme(axis.text.x = element_text(angle = -45,vjust = 0.5))

# Strip problems
lil.strip <- theme(strip.background = element_blank(),
                   strip.text.x =element_text(margin = margin(.05, 0, .1, 0, "cm")),
                   strip.text = element_text(size = 10),
                   axis.text.y = element_text(size = 10)) 


#### $$ TIME $$ #####

## ~~~ Season~~~

season.order <- c("winter","spring","summer","autumn")

sea.shape <- c(15,16,17,18)
sea.col <- c("#63bfb8","#52b216","#ecee21","#e3a40c")

names(sea.shape) <- season.order
names(sea.col) <- season.order

sea.shapeScale <- scale_shape_manual(name = "Season",
                                     values = sea.shape,
                                     limits = season.order)
sea.colScale <- scale_color_manual(name = "Season",
                                   values = sea.col,
                                   limits = season.order)
sea.fillScale <- scale_fill_manual(name = "Season",
                                   values = sea.col,
                                   limits = season.order)

## ~~~Month~~~~

mon.col <- c(hc_pal()(10)[-5],
             "#d11919","#B37400","#4C4F8B")

month.num <- c("01", "02",  "03", "04",
              "05", "06", "07", "08",
              "09", "10", "11", "12")
               
month.order <- c("jan","feb","mar","apr","may","jun",
                 "jul","aug","sep","oct","nov","dec")

names(mon.col) <- month.num

mon.colScale <- scale_colour_manual(name = "Month",
                                    values=c(hc_pal()(10)[-5],
                                             "#d11919","#B37400","#4C4F8B"),
                                    limits=month.num,
                                    labels = str_to_title(month.order))

mon.fillScale <- scale_fill_manual(name = "Month",
                                    values=c(hc_pal()(10)[-5],
                                             "#d11919","#B37400","#4C4F8B"),
                                    limits=month.num,
                                    labels = str_to_title(month.order))

scale_x_month <- scale_x_discrete(name = 'Month',
                                  labels = str_to_title(month.order))
## ~~~Decade plots~~~ 

darkzone <- data.frame(x = c(2003:2009) + 0.75,
                       x2 = c(2004:2010) + 0.25)

thedark <- geom_rect(data=darkzone,
                     fill="grey",alpha=0.4,
                     aes(xmin=x,
                         xmax=x2),
                     ymin=-Inf,ymax=Inf, inherit.aes = FALSE)

## ~~~~ Themes plots 

theme_month <- list(
)
  
theme_timeseries <- list(
  thedark ,
  scale_x_continuous(breaks = seq(2004, 2010, 1)) ,
  xlab("Year") ,
  ylab("Relative abundance"),
  lil.strip)


#### ~~~~~Location~~~~~~~ ####

# Implemented for gradients analysis, around 20180606
location.order <-  c('Girona','Barcelona',
                     'Tarragona','Palma')

colors.location <- c("green4", "#F8766D", "#C71CFF", "dodgerblue")

names(colors.location) <- location.order

loc.colScale <- scale_color_manual(name = "Location",
                                   values = colors.location,
                                   limits = location.order)

#### Fraction ####

# Implemented for gradients analysis, around 20180606
frac.linetype <- list(scale_linetype_manual(name = 'Fraction',
                                       values = c(1,5),
                                       breaks = c('pico','nano'),
                                       limits = c('pico', 'nano')),
  guides(linetype = F))


