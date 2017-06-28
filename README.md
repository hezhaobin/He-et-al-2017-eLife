# Code for He et al elife 2017 paper

## Overview

This repository includes raw data and code for reproducing Figure 3A-D, 4B-D, 5B-C and figure supplements for 3 and 4 in the paper.

The repository is organized into two folders. [Figure 3 and 4](Figure_3_and_4) are combined because they both use data obtained in the _S. cerevisiae_ background. [Figure 5](Figure_5) uses data obtained in the _C. glabrata_ background.

Each folder contains one or more [R Notebook files](http://rmarkdown.rstudio.com/r_notebooks.html). Data used for plotting are stored in the "input" folder therein and the final figures are in "figure_output" folder.

## Machine and software specifications

Running the code requires installation of the latest version of R and Rstudio. Below are the output of the session information used to generate all results in this repository.

```r
> sessionInfo()
R version 3.3.2 (2016-10-31)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X Yosemite 10.10.5

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods  
[9] base     

other attached packages:
 [1] GenomicRanges_1.22.4 GenomeInfoDb_1.6.3   IRanges_2.4.8       
 [4] S4Vectors_0.8.11     edgeR_3.12.0         RColorBrewer_1.1-2  
 [7] data.table_1.10.4    limma_3.26.8         NMF_0.20.6          
[10] Biobase_2.30.0       BiocGenerics_0.16.1  cluster_2.0.5       
[13] rngtools_1.2.4       pkgmaker_0.22        registry_0.3        
[16] plyr_1.8.4           cowplot_0.7.0        ggplot2_2.2.1       

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.10      XVector_0.10.0    knitr_1.16        magrittr_1.5     
 [5] zlibbioc_1.16.0   doParallel_1.0.10 munsell_0.4.3     colorspace_1.2-6 
 [9] xtable_1.8-2      gridBase_0.4-7    foreach_1.4.3     stringr_1.2.0    
[13] tools_3.3.2       grid_3.3.2        gtable_0.2.0      iterators_1.0.8  
[17] lazyeval_0.2.0    digest_0.6.9      tibble_1.3.0      reshape2_1.4.2   
[21] codetools_0.2-15  stringi_1.0-1     scales_0.4.1     
```

