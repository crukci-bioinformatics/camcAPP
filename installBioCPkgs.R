install.packages(c("devtools","shiny"))
library(devtools)
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("ggplot2","dplyr","tidyr","gridExtra","heatmap.plus","RColorBrewer","party","survival","knitr"))
install_github("markdunning/prostateCancerTaylor")