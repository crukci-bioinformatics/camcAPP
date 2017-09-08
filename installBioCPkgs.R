options(repos = c("CRAN" = "http://cran.ma.imperial.ac.uk"))
install.packages(c("devtools","shiny"))
library(devtools)
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("org.Hs.eg.db","Biobase","ggplot2","dplyr","dbplyr","tidyr","gridExtra","heatmap.plus",
           "RColorBrewer","party","survival","knitr","broom","WGCNA","gplots",
           "prostateCancerTaylor","prostateCancerCamcap","prostateCancerStockholm","prostateCancerGrasso","prostateCancerVarambally","GGally","DT","ggthemes","stringr","COSMIC.67"))
