# Scripts for *Seasonal and interannual variability of the free-living and particle-associated bacteria of a coastal microbiome*

Author: Adrià Auladell

Inside `src` there are all the scripts with the following structure: 

```
src/
    data/ # processing of raw files with DADA2 + phyloseq
    analysis/ # all the analysis generating new data or statistics
    figures/  # the visualizations, mainly through ggplot
    utils/ # scripts called from the abovementioned scripts
```

Inside the `data/cleaned` folder there is a phyloseq R object with the sequence data, abundance table and sample metadata. 
To execute most of the analysis, create an R project in the main folder and execute the desired scripts. Most of the `figures` scripts 
require to have executed scripts in the `analysis` folder for the statistics to be saved in the correspondent folders. 

You will need some minimal packages for running the scripts. For most analysis: 

- [tidyverse](https://www.tidyverse.org/): data wrangling, visualization and general purpose package.
- [phyloseq](https://joey711.github.io/phyloseq/): specific microbiome experiment functions. 
- [lomb](https://cran.r-project.org/web/packages/lomb/index.html): seasonality testing.
