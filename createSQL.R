library(dplyr)
library(Biobase)
db <- src_sqlite("camcap.sqlite3", create=TRUE)
  
data(camcap,package = "prostateCancerCamcap")
pd_camcap <- tbl_df(pData(camcap))
fd_camcap <- tbl_df(fData(camcap))
expression <-
  tbl_df(data.frame(ID = as.character(featureNames(camcap)),exprs(camcap))) %>% 
  gather(geo_accession,Expression,-ID)
pd <- pd_camcap
fd <- fd_camcap
expression <- left_join(expression, select(fd, ID, Symbol))

copy_to(db, expression, temporary = FALSE, indexes = list("Symbol"))
copy_to(db, pd, temporary = FALSE, indexes = list("geo_accession"))
copy_to(db, fd, temporary = FALSE, indexes = list("ID"))

copyNumber <- read.delim("data/2014-06-04_UK_OncoSNP_rank3_gene_matrix.txt") %>% 
  gather(Sample, Call, -(EntrezID:position)) %>% 
  tbl_df %>% 
  select(Symbol,Sample,Call)
copy_to(db, copyNumber, temporary = FALSE, indexes = list("Symbol"))


db <- src_sqlite("stockholm.sqlite3", create=TRUE)

data(stockholm,package = "prostateCancerStockholm")
pd_stockholm <- tbl_df(pData(stockholm))
fd_stockholm <- tbl_df(fData(stockholm))
expression <-
  tbl_df(data.frame(ID = as.character(featureNames(stockholm)),exprs(stockholm))) %>% 
  gather(geo_accession,Expression,-ID)
fd <- fd_stockholm
pd <- pd_stockholm

expression <- left_join(expression, select(fd, ID, Symbol))

copy_to(db, expression, temporary = FALSE, indexes = list("Symbol"))
copy_to(db, pd, temporary = FALSE, indexes = list("geo_accession"))
copy_to(db, fd, temporary = FALSE, indexes = list("ID"))

copyNumber <- read.delim("data/2016-01-06_Stockholm_0.9_CNA_threshold_geneXsample_mat.txt") %>% 
  gather(Sample, Call, -(EntrezID:Symbol)) %>% 
  tbl_df %>% 
  select(Symbol,Sample,Call)
copy_to(db, copyNumber, temporary = FALSE, indexes = list("Symbol"))

db <- src_sqlite("taylor.sqlite3",create=TRUE)

data(taylor,package = "prostateCancerTaylor")
pd_taylor <- tbl_df(pData(taylor))
fd_taylor <- tbl_df(fData(taylor))
expression <-
  tbl_df(data.frame(ID = as.character(featureNames(taylor)),log2(exprs(taylor)))) %>% 
  gather(geo_accession,Expression,-ID)
fd <- fd_taylor
pd <- pd_taylor
expression <- left_join(expression, select(fd, ID, Gene))

copy_to(db, expression, temporary = FALSE, indexes = list("Gene"))
copy_to(db, pd, temporary = FALSE, indexes = list("geo_accession"))
copy_to(db, fd, temporary = FALSE, indexes = list("ID"))

copyNumber <- read.delim("data/2014-12-02_cna_matrix_localized_cap_w_clinical.txt") %>% 
  gather(Sample, Call, -(EntrezID:GeneSymbol)) %>% 
  tbl_df %>% 
  mutate(Symbol = GeneSymbol) %>% 
  select(Symbol,Sample,Call)
copy_to(db, copyNumber, temporary = FALSE, indexes = list("Symbol"))


db <- src_sqlite("varambally.sqlite3",create=TRUE)

data(varambally,package = "prostateCancerVarambally")
pd_varambally <- tbl_df(pData(varambally))
fd_varambally <- tbl_df(fData(varambally))
expression <-
  tbl_df(data.frame(ID = as.character(featureNames(varambally)),exprs(varambally))) %>% 
  gather(geo_accession,Expression,-ID)
fd <- fd_varambally
pd <- pd_varambally
expression <- left_join(expression, select(fd, ID, Symbol))

copy_to(db, expression, temporary = FALSE, indexes = list("Symbol"))
copy_to(db, pd, temporary = FALSE, indexes = list("geo_accession"))
copy_to(db, fd, temporary = FALSE, indexes = list("ID"))

db <- src_sqlite("varambally.sqlite3",create=TRUE)

data(varambally,package = "prostateCancerVarambally")
pd_varambally <- tbl_df(pData(varambally))
fd_varambally <- tbl_df(fData(varambally))
expression <-
  tbl_df(data.frame(ID = as.character(featureNames(varambally)),exprs(varambally))) %>% 
  gather(geo_accession,Expression,-ID)
fd <- fd_varambally
pd <- pd_varambally
expression <- left_join(expression, select(fd, ID, Symbol))

copy_to(db, expression, temporary = FALSE, indexes = list("Symbol"))
copy_to(db, pd, temporary = FALSE, indexes = list("geo_accession"))
copy_to(db, fd, temporary = FALSE, indexes = list("ID"))

db <- src_sqlite("grasso.sqlite3",create=TRUE)

data(grasso,package = "prostateCancerGrasso")
pd_grasso <- tbl_df(pData(grasso))
fd_grasso <- tbl_df(fData(grasso))
expression <-
  tbl_df(data.frame(ID = as.character(featureNames(grasso)),exprs(grasso))) %>% 
  gather(geo_accession,Expression,-ID)
fd <- fd_grasso
pd <- pd_grasso
expression <- left_join(expression, select(fd, ID, GENE_SYMBOL))

copy_to(db, expression, temporary = FALSE, indexes = list("GENE_SYMBOL"))
copy_to(db, pd, temporary = FALSE, indexes = list("geo_accession"))
copy_to(db, fd, temporary = FALSE, indexes = list("ID"))


