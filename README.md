# Mining human prostate cancer datasets: The 'camcAPP' shiny app

camcAPP is a shiny app that allows various Prostate Cancer datasets to be interrogated through a web-based interface.



## The Datasets

The following datasets are accessible through the interface. Each of which has been made available through Gene Expression Omnibus (GEO). We have re-packaged each dataset into a Bioconductor-compatible package

- [Cambridge](http://www.sciencedirect.com/science/article/pii/S2352396415300712)
    + [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70768)
    + [Bioconductor](http://bioconductor.org/packages/release/data/experiment/html/prostateCancerCamcap.html)
- [Stockholm](http://www.sciencedirect.com/science/article/pii/S2352396415300712)
    + [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70769)
    + [Bioconductor](http://bioconductor.org/packages/release/data/experiment/html/prostateCancerStockholm.html)
- [MSKCC](http://www.nature.com/nature/journal/v487/n7406/full/nature11125.html)
    + [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21034)
    + [Bioconductor](http://bioconductor.org/packages/release/data/experiment/html/prostateCancerTaylor.html)
- [Michigan2005](http://www.sciencedirect.com/science/article/pii/S1535610805003053)
    + [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3325)
    + [Bioconductor](http://bioconductor.org/packages/release/data/experiment/html/prostateCancerVarambally.html)
- [Michigan2012](http://www.nature.com/nature/journal/v487/n7406/full/nature11125.html)
    + [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35988)
    + [Bioconductor](http://bioconductor.org/packages/release/data/experiment/html/prostateCancerGrasso.html)

## Reproducibility

The code to run particular analyses and produce plots can be downloaded from the Shiny interface. Details of the R packages required are given in the file [installBioCPkgs.R](https://raw.githubusercontent.com/crukci-bioinformatics/camcAPP/master/installBioCPkgs.R). 

