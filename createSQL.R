library(dplyr)
library(tidyr)
library(Biobase)
library(COSMIC.67)



data(cgc_67, package = "COSMIC.67")

curated.genes <- na.omit(cgc_67[,1])
curated.genes <- c(curated.genes, c('ESR1','AR','AARS','ACAD9','ACVR1','AFG3L2','AKAP11','ANKFY1','ANKRD46','ARF5','ATP6V1B2','BRF2','C18orf25','C8orf58','C8orf76','CCAR2','CDC16','CHMP4C','CHMP7','COL20A1','COPS7A','COX10','COX4I1','CPNE3','CTNS','CUEDC2','CXXC1','CYB5D2','DCTN2','DEF8','DMTN','EBAG9','EIF4G1','ELAC2','ELP3','ERCC3','ERI1','FAM160B2','FIG4','FSCN3','FUK','GALNT12','GHITM','GPR88','GTF2E2','GTF2F2','HDDC2','HDHD2','HMCES','HSBP1','IMPA1','INTS9','KAT6A','KCNK12','KLHDC4','LINC00938','LPPR1','LRRC31','LSM1','MAP3K7','MAPKAP1','MELK','MMS19','MPDU1','MTMR9','MYLK2','NSMCE2','NUDT16P1','PHF11','PIGN','PPP3CC','PQLC1','PTDSS1','PTK2B','R3HCC1','RALBP1','RARS2','RBFA','RCBTB2','RIPK2','RNGTT','RP11-60I3.5','SAT2','SCAF11','SCN4A','SETDB2','SFT2D3','SLC25A11','SORBS3','SPAG7','SPIDR','SUGT1','TCF25','TM2D2','TPD52','TRAPPC2L','TRIM13','WDR59','YEATS2','ZBTB34','ZBTB4','ZNF252P'))
write.table(curated.genes, file="curated.genes.txt",row.names = FALSE,quote=FALSE)

db <- src_sqlite("camcap.sqlite3", create=TRUE)
db.curated <- src_sqlite("camcap.curated.sqlite3", create=TRUE)

data(camcap,package = "prostateCancerCamcap")
pd_camcap <- tbl_df(pData(camcap))
fd_camcap <- tbl_df(fData(camcap))
expression <-
  tbl_df(data.frame(ID = as.character(featureNames(camcap)),exprs(camcap))) %>% 
  gather(geo_accession,Expression,-ID)
pd <- pd_camcap
fd <- fd_camcap
expression.all <- left_join(expression, select(fd, ID, Symbol))

expression <- filter(expression.all, Symbol %in% curated.genes)
  
copy_to(db.curated, expression, temporary = FALSE, indexes = list("Symbol"))
copy_to(db.curated, pd, temporary = FALSE, indexes = list("geo_accession"))
copy_to(db.curated, fd, temporary = FALSE, indexes = list("ID"))

expression <- expression.all

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


