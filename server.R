startLoad <- date()
message("Started loading at:")
print(startLoad)
gc()

library(shiny)
library(ggplot2)

library(tidyr)
library(devtools)
library(gridExtra)
library(heatmap.plus)
library(gplots)
library(RColorBrewer)

library(knitr)
library(Biobase)
library(broom)
library(GGally)

library(DT)



#if(!require(prostateCancerTaylor)) install_github("crukci-bioinformatics/prostateCancerTaylor");library(prostateCancerTaylor)
#if(!require(prostateCancerCamcap)) install_github("crukci-bioinformatics/prostateCancerCamcap");library(prostateCancerCamcap)
#if(!require(prostateCancerStockholm)) install_github("crukci-bioinformatics/prostateCancerStockholm");library(prostateCancerStockholm)

plotDendroAndColors.mod <- function (dendro, colors, groupLabels = NULL, rowText = NULL, 
                                     rowTextAlignment = c("left", "center", "right"), rowTextIgnore = NULL, 
                                     textPositions = NULL, setLayout = TRUE, autoColorHeight = TRUE, 
                                     colorHeight = 0.2, rowWidths = NULL, dendroLabels = NULL, 
                                     addGuide = FALSE, guideAll = FALSE, guideCount = 50, guideHang = 0.2, 
                                     addTextGuide = FALSE, cex.colorLabels = 0.8, cex.dendroLabels = 0.9, 
                                     cex.rowText = 0.8, marAll = c(1, 5, 3, 1), saveMar = TRUE, 
                                     abHeight = NULL, abCol = "red",cutType="k",k=2, ...) 
{
  oldMar = par("mar")
  if (!is.null(dim(colors))) {
    nRows = dim(colors)[2]
  }
  else nRows = 1
  if (!is.null(rowText)) 
    nRows = nRows + if (is.null(textPositions)) 
      nRows
  else length(textPositions)
  if (autoColorHeight) 
    colorHeight = 0.2 + 0.3 * (1 - exp(-(nRows - 1)/6))
  if (setLayout) 
    layout(matrix(c(1:2), 2, 1), heights = c(1 - colorHeight, 
                                             colorHeight))
  par(mar = c(0, marAll[2], marAll[3], marAll[4]))
  plot(dendro, labels = dendroLabels, cex = cex.dendroLabels, 
       ...)
  if(cutType == "k"){
    rect.hclust(dendro,k=k,border="blue")
    abHeight <- NULL
  } else rect.hclust(dendro, h=abHeight,border="blue")
  
  if (addGuide) 
    addGuideLines(dendro, count = if (guideAll) 
      length(dendro$height) + 1
      else guideCount, hang = guideHang)
  if (!is.null(abHeight)) 
    abline(h = abHeight, col = abCol)
  par(mar = c(marAll[1], marAll[2], 0, marAll[4]))
  plotColorUnderTree(dendro, colors, groupLabels, cex.rowLabels = cex.colorLabels, 
                     rowText = rowText, rowTextAlignment = rowTextAlignment, 
                     rowTextIgnore = rowTextIgnore, textPositions = textPositions, 
                     cex.rowText = cex.rowText, rowWidths = rowWidths, addTextGuide = addTextGuide)
  if (saveMar) 
    par(mar = oldMar)
}

library(dplyr)

data(camcap,package = "prostateCancerCamcap")
db_camcap <- src_sqlite("camcap.sqlite3")
#pd_camcap <- tbl_df(pData(camcap))

pd_camcap <- collect(tbl("pd",src=db_camcap))
fd_camcap <- collect(tbl("fd",src=db_camcap))

#fd_camcap <- tbl_df(fData(camcap))
#exp_camcap <-
 # tbl_df(data.frame(ID = as.character(featureNames(camcap)),exprs(camcap)))

exp_camcap <- tbl("expression",src=db_camcap)

#loadStockholm <- function(){

#data(stockholm, package="prostateCancerStockholm")
db_stockholm <- src_sqlite("stockholm.sqlite3")
#pd_stockholm <- tbl_df(pData(stockholm))
#fd_stockholm <- tbl_df(fData(stockholm))
pd_stockholm <- collect(tbl("pd",src=db_stockholm))
fd_stockholm <- collect(tbl("fd",src=db_stockholm))


#exp_stockholm <- tbl_df(data.frame(ID=as.character(featureNames(stockholm)),exprs(stockholm)))
#}
exp_stockholm <- tbl("expression",src=db_stockholm)



#loadTaylor <- function(){

#data(taylor, package="prostateCancerTaylor")
db_taylor <- src_sqlite("taylor.sqlite3")
pd_taylor <- collect(tbl("pd",src=db_taylor))
fd_taylor <- collect(tbl("fd",src=db_taylor))
#exp_taylor <- tbl_df(data.frame(ID=as.character(featureNames(taylor)),log2(exprs(taylor))))
exp_taylor <- tbl("expression",src=db_taylor)

#}

#loadVarambally <- function(){
db_varambally <- src_sqlite("varambally.sqlite3")
#data(varambally,package="prostateCancerVarambally")
pd_varambally <- collect(tbl("pd",src=db_varambally))
fd_varambally <- collect(tbl("fd",src=db_varambally))

#exp_varambally <- tbl_df(data.frame(ID=as.character(featureNames(varambally)),log2(exprs(varambally))))
exp_varambally <- tbl("expression",src=db_varambally)

#}

#loadGrasso <- function(){

#data(grasso,package="prostateCancerGrasso")
db_grasso <- src_sqlite("grasso.sqlite3")
pd_grasso <- collect(tbl("pd",src=db_grasso))
fd_grasso <- collect(tbl("fd",src=db_grasso))

#exp_grasso <- tbl_df(data.frame(ID=as.character(featureNames(grasso)),exprs(grasso)))
exp_grasso <- tbl("expression",src=db_grasso)
#}


## dplyr needs to be the last package loaded, otherwise the 'select' function seems to be over-written

select <- dplyr::select
iclusPal <- brewer.pal(5, "Set1")
gradeCols <- rev(brewer.pal(11, "RdYlGn"))
message("READY FOR INPUT")



#library(org.Hs.eg.db)

#gList <- unique(as.character(as.matrix(select(fd_camcap,Symbol))))
#gList <- union(gList, as.character(as.matrix(select(fd_taylor,Gene))))


#geneTable <- AnnotationDbi:::select(org.Hs.eg.db, keys=gList, keytype = "SYMBOL",columns=c("SYMBOL","GENENAME"))

endLoad <- date()
message("Ended loading at:")
print(endLoad)

shinyServer(function(input, output,session){
  

  
  getCurrentGene <- reactive({
    if(input$inputType == "Single Gene"){
      updateTextInput(session, inputId = "profileBasename", value=paste0(input$currentGene,"-profile"))
      updateTextInput(session, inputId = "survivalBasename", value=paste0(input$currentGene,"-survival"))
      updateTextInput(session, inputId = "copyNumberBasename", value=paste0(input$currentGene,"-copyNumber"))
      updateTextInput(session, inputId = "correlationBasename", value=paste0(input$currentGene,"-versus-",getSecondGene()))
    }
    input$currentGene
  })
  
  output$geneTable <- DT:::renderDataTable(geneTable)
  
  #  getCurrentGene <- reactive({
  #   input$currentGene_cambridge
  #  })
  
  # getCurrentGene <- reactive({
  #    input$currentGene_stockholm
  #  })
  
  #  getCurrentGene <- reactive({
  #    input$currentGene_taylor
  #  })get
  
  
  
  
  getCambridgeVariable <- reactive({
    input$clinvar_boxplot
  })
  
  getStockholmVariable <- reactive({
    input$clinvar_stockholm
  })
  
  getTaylorVariable <- reactive({
    input$clinvar_taylor
  })
  
  getVaramballyVariable <- reactive({
    input$clinvar_varambally
  })
  
  getGrassoVariable <- reactive({
    input$clinvar_grasso
  })
  
  
  
  getCambridgeOverlay <- reactive({
    input$overlay_cambridge
  })
  getStockholmOverlay <- reactive({
    input$overlay_stockholm
  })
  getTaylorOverlay <- reactive({
    input$overlay_taylor
  })
  getVaramballyOverlay <- reactive({
    input$overlay_varambally
  })
  
  getGrassoOverlay <- reactive({
    input$overlay_grasso
  })
  
  
  getCambridgeZ <- reactive({
    input$z_cambridge
  })
  
  getStockholmZ <- reactive({
    input$z_stockholm
  })
  getTaylorZ <- reactive({
    input$z_taylor
  })
  getVaramballyZ <- reactive({
    input$z_varambally
  })
  getGrassoZ <- reactive({
    input$z_grasso
  })
  
  getCambridgeCompositePlot <- reactive({
    input$cambridgeCombPlot
  })
  
  
  getRpDataset <- reactive({
    input$rpDataset
  })
  
  getHeatmapDataset <- reactive({
    input$heatmapDataset
  })
  
  getDistFun <- reactive({
    input$distfun
    
  })
  
  getReordRows <- reactive({
    input$reordRows
  })
  
  getHclustMethod <- reactive({
    input$hclustfun
  })
  
  getScaleMethod <- reactive({
    input$scale
  })
  
  getCorDataset <- reactive({
    dataset <- input$corDataset
    
    covars <- switch(dataset,
                     Cambridge=c("None", "iCluster","Gleason","Sample_Group"),
                     Stockholm=c("None", "iCluster","Gleason"),
                     MSKCC = c("None", "CopyNumberCluster","Gleason"),
                     Michigan2005 = c("None", "Sample_Group"),
                     Michigan2012 = c("None", "Sample_Group")
    )
    
    updateSelectInput(session, inputId="clinvar_cor", choices=covars,selected=covars[1])
    
    dataset
    
    
  })
  
  getSecondGene <- reactive({
    
    gene2 <- input$secondGene
    
  })
  
  getCorType <- reactive({
    input$corType
  })
  
  getProfileBasename <- reactive({
    input$profileBasename
  })
  
  getSurvivalBasename <- reactive({
    input$survivalBasename
    
  })
  
  getHeatmapBasename <- reactive({
    input$heatmapBasename
    
  })
  
  getBoxplotDataset <- reactive({
    dataset <- input$boxplotDataset
    covars <- switch(dataset,
                     Cambridge=c("iCluster","Gleason","Sample_Group"),
                     Stockholm=c("iCluster","Gleason"),
                     MSKCC = c("CopyNumberCluster","Gleason"),
                     Michigan2005 = "Sample_Group",
                     Michigan2012 = "Sample_Group"
    )
    
    updateSelectInput(session, inputId="clinvar_boxplot", choices=covars,selected=covars[1])
    dataset
    
  })
  
  
  
  
  
  getGeneList <- reactive({inFile <- input$file1
  
  if (is.null(inFile))
    return(c("STAT3","ESR1","AR"))
  
  
  #    print(inFile$name)
  genes <- read.delim(inFile$datapath)[,1]
  #   print(dim(genes))
  #   if(input$inputType != "Single Gene") updateTextInput(session, inputId = "profileBasename",value=paste0(basename(input$datapath),"-profile"))
  updateSelectInput(session, inputId = "cambridgeGeneChoice",choices = as.character(genes),selected=as.character(genes)[1])
  updateSelectInput(session, inputId = "survivalGeneChoice",choices = as.character(genes),selected=as.character(genes)[1])
  
  if(input$inputType != "Single Gene"){
    inFile <- input$file1
    updateTextInput(session, inputId = "profileBasename", value=paste0(basename(inFile$name),"-profile"))
    updateTextInput(session, inputId = "survivalBasename", value=paste0(as.character(genes)[1],"-survival"))
    updateTextInput(session, inputId = "heatmapBasename", value=paste0(basename(inFile$name),"-heatmap"))
    updateTextInput(session, inputId = "copyNumberBasename", value=paste0(basename(inFile$name),"-copyNumber"))
  }
  
  genes
  
  })
  
  
  output$profileMapping <- renderPrint({
    
    
    plotType <- input$inputType
    
    
    if(plotType != "Single Gene") {
      
      dataset <- getBoxplotDataset()
      
      genes <- getGeneList()
      
      
      if(dataset == "Cambridge"){
        
        check <- all(genes %in% fd_camcap$Symbol)
        
        if(!check) genediff <- setdiff(genes, fd_camcap$Symbol)
        
      }
      
      else if (dataset == "Stockholm"){
        check <- all(genes %in% fd_stockholm$Symbol)
        if(!check) genediff <- setdiff(genes, fd_stockholm$Symbol)
        
      }
      
      else if(dataset == "MSKCC"){
        check <- all(genes %in% fd_taylor$Gene)
        if(!check) genediff <- setdiff(genes, fd_taylor$Gene)
      }
      
      else if(dataset == "Michigan2005"){
        check <- all(genes %in% fd_varambally$Symbol)
        if(!check) genediff <- setdiff(genes, fd_varambally$Gene)
        
      }
      
      else if (dataset == "Michigan2012"){
        check <- all(genes %in% fd_grasso$GENE_SYMBOL)
        if(!check) genediff <- setdiff(genes, fd_grasso$GENE_SYMBOL)
      }
      
      if(check) msg <- "All genes were sucessfully mapped to the dataset"
      else {
        msg <-paste("Genes that could not be mapped to the dataset:-", genediff)
        validGenes <- as.character(genes[-which(genes %in% genediff)])
        updateSelectInput(session, "cambridgeGeneChoice", choices=validGenes, selected = validGenes[1])
        
      }
      
      msg
    }    
  })
  
  
  
  output$heatmapMapping <- renderPrint({
    
    
    
    dataset <- getHeatmapDataset()
    
    genes <- getGeneList()
    
    
    if(dataset == "Cambridge"){
      
      check <- all(genes %in% fd_camcap$Symbol)
      
      if(!check) genediff <- setdiff(genes, fd_camcap$Symbol)
      
    }
    
    else if (dataset == "Stockholm"){
      check <- all(genes %in% fd_stockholm$Symbol)
      if(!check) genediff <- setdiff(genes, fd_stockholm$Symbol)
      
    }
    
    else if(dataset == "MSKCC"){
      check <- all(genes %in% fd_taylor$Gene)
      if(!check) genediff <- setdiff(genes, fd_taylor$Gene)
    }
    
    else if(dataset == "Michigan2005"){
      check <- all(genes %in% fd_varambally$Symbol)
      if(!check) genediff <- setdiff(genes, fd_varambally$Gene)
      
    }
    
    else if (dataset == "Michigan2012"){
      check <- all(genes %in% fd_grasso$GENE_SYMBOL)
      if(!check) genediff <- setdiff(genes, fd_grasso$GENE_SYMBOL)
    }
    
    if(check) msg <- "All genes were sucessfully mapped to the dataset"
    else {
      msg <-paste("Genes that could not be mapped to the dataset:-", genediff)
      
    }
    
    msg
    
  })
  
  
  
  output$survivalMapping <- renderPrint({
    
    plotType <- input$inputType
    
    
    if(plotType != "Single Gene") {
      
      dataset <- getRpDataset()
      
      genes <- getGeneList()
      
      
      if(dataset == "Cambridge"){
        
        check <- all(genes %in% fd_camcap$Symbol)
        
        if(!check) genediff <- setdiff(genes, fd_camcap$Symbol)
        
      }
      
      else if (dataset == "Stockholm"){
        check <- all(genes %in% fd_stockholm$Symbol)
        if(!check) genediff <- setdiff(genes, fd_stockholm$Symbol)
        
      }
      
      else if(dataset == "MSKCC"){
        check <- all(genes %in% fd_taylor$Gene)
        if(!check) genediff <- setdiff(genes, fd_taylor$Gene)
      }
      
      else if(dataset == "Michigan2005"){
        check <- all(genes %in% fd_varambally$Symbol)
        if(!check) genediff <- setdiff(genes, fd_varambally$Gene)
        
      }
      
      else if (dataset == "Michigan2012"){
        check <- all(genes %in% fd_grasso$GENE_SYMBOL)
        if(!check) genediff <- setdiff(genes, fd_grasso$GENE_SYMBOL)
      }
      
      if(check) msg <- "All genes were sucessfully mapped to the dataset"
      else {
        msg <-paste("Genes that could not be mapped to the dataset:-", genediff)
        validGenes <- as.character(genes[-which(genes %in% genediff)])
        updateSelectInput(session, "survivalGeneChoice", choices=validGenes, selected = validGenes[1])
      }
      
      msg
    }
  })
  
  
  prepareBoxplot <- reactive({
    
    dataset <- getBoxplotDataset()
    
    
    plotType <- input$inputType
    if(plotType == "Single Gene") {
      genes <- getCurrentGene()
    } else  genes <- getGeneList()
    
    if(dataset == "Cambridge"){
      
#      probes <- fd_camcap %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
       data <- collect(exp_camcap)  %>% filter(Symbol %in% genes)
        #gather(geo_accession,Expression,-ID)
      fd <- fd_camcap
      pd <-  mutate(pd_camcap, Gleason=factor(Gleason,levels=c("5=3+2","6=2+4","6=3+3", "7=3+4","7=4+3","8=3+5","8=4+4","9=4+5","9=5+4","10=5+5",NA))) %>% 
        mutate(Time = as.numeric(FollowUpTime), Event = ifelse(BCR=="Y",1,0))
      
    } else if (dataset == "Stockholm"){
      
#      probes <- fd_stockholm %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
      data<- collect(exp_stockholm)  %>% filter(Symbol %in% genes)
#        gather(geo_accession,Expression,-ID)
      fd <- fd_stockholm
      pd <-  mutate(pd_stockholm, Gleason=factor(Gleason,levels=c("5=3+2","6=2+4","6=3+3", "7=3+4","7=4+3","8=3+5","8=4+4","9=4+5","9=5+4","10=5+5",NA))) %>% 
        mutate(Time = as.numeric(FollowUpTime), Event = ifelse(BCR=="Y",1,0))
    }
    
    else if(dataset == "MSKCC"){
      
      #probes <- fd_taylor %>% filter(Gene %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
#      data <- exp_taylor  %>% filter(ID %in% probes) %>% 
 #       gather(geo_accession,Expression,-ID)
      data <- collect(exp_taylor) %>% filter(Gene %in% genes)
      
      
      fd <- fd_taylor %>% mutate(Symbol = Gene)
      pd <- mutate(pd_taylor,Gleason = gsub("4+3", "7=4+3", pd_taylor$Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+3", "8=5+3", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("3+4", "7=3+4", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+3", "7=4+3", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("3+3", "6=3+3", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+5", "9=4+5", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("3+5", "8=3+5", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+4", "8=4+4", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason=factor(Gleason,levels=c("5=3+2","6=2+4","6=3+3", "7=3+4","7=4+3","8=3+5","8=4+4","9=4+5","9=5+4","10=5+5",NA)))
    }
    
    else if(dataset == "Michigan2005"){
      
      #probes <- fd_varambally %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
      #data <- exp_varambally  %>% filter(ID %in% probes) %>% 
       # gather(geo_accession,Expression,-ID)
      data <- collect(exp_varambally) %>% filter(Symbol %in% genes)
      
      fd <- fd_varambally
      pd <- pd_varambally 
    }
    
    else {
      #probes <- fd_grasso %>% filter(GENE_SYMBOL %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
      #data <- exp_grasso  %>% filter(ID %in% probes) %>% 
      #  gather(geo_accession,Expression,-ID)
      data <- collect(exp_grasso) %>% filter(GENE_SYMBOL %in% genes)
      fd <- mutate(fd_grasso, Symbol = GENE_SYMBOL)
      pd <- mutate(pd_grasso, Sample_Group = Group) 
      
    }
    
    
    summary_stats <- data %>% group_by(ID) %>% 
      summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
    
    data <- left_join(data,summary_stats) %>% mutate(Z = (Expression - mean) / sd)
    
    #    mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
    #    mu <- summary_stats$mean[which.max(summary_stats$iqr)]
    #    sd <- summary_stats$sd[which.max(summary_stats$iqr)]
    
    #    data <- filter(data, ID== mostVarProbe) %>%
    #      mutate(Z = (Expression -mu) /sd)
    
    #    data <- left_join(data,pd)
    #    data <- left_join(data, fd)
    
    
    
    mostVarProbes <- left_join(summary_stats,fd) %>% 
      arrange(Symbol,desc(iqr)) %>% 
      distinct(Symbol) %>% 
      select(ID) %>%  as.matrix %>%  as.character
    
    data <- filter(data, ID %in% mostVarProbes)
    data <- left_join(data, select(fd, ID, Symbol))
    data <- left_join(data, pd)
    data    
  
    })
  
  
  filterByGene <-function(dataset, genes){
    
    
    if(dataset == "Cambridge"){
      
      probes <- fd_camcap %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
      data<- collect(exp_camcap)  %>% filter(ID %in% probes) %>% 
        gather(geo_accession,Expression,-ID)
      fd <- fd_camcap
      pd <-  mutate(pd_camcap, Gleason=factor(Gleason,levels=c("5=3+2","6=2+4","6=3+3", "7=3+4","7=4+3","8=3+5","8=4+4","9=4+5","9=5+4","10=5+5",NA))) %>% 
        mutate(Time = as.numeric(FollowUpTime), Event = ifelse(BCR=="Y",1,0))
      
    } else if (dataset == "Stockholm"){
      
      probes <- fd_stockholm %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
      data<- exp_stockholm  %>% filter(ID %in% probes) %>% 
        gather(geo_accession,Expression,-ID)
      fd <- fd_stockholm
      pd <-  mutate(pd_stockholm, Gleason=factor(Gleason,levels=c("5=3+2","6=2+4","6=3+3", "7=3+4","7=4+3","8=3+5","8=4+4","9=4+5","9=5+4","10=5+5",NA))) %>% 
        mutate(Time = as.numeric(FollowUpTime), Event = ifelse(BCR=="Y",1,0))
    }
    
    else if(dataset == "MSKCC"){
      
      probes <- fd_taylor %>% filter(Gene %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
      data <- exp_taylor  %>% filter(ID %in% probes) %>% 
        gather(geo_accession,Expression,-ID)
      fd <- fd_taylor %>% mutate(Symbol = Gene)
      pd <- mutate(pd_taylor,Gleason = gsub("4+3", "7=4+3", pd_taylor$Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+3", "8=5+3", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("3+4", "7=3+4", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+3", "7=4+3", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("3+3", "6=3+3", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+5", "9=4+5", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("3+5", "8=3+5", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+4", "8=4+4", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason=factor(Gleason,levels=c("5=3+2","6=2+4","6=3+3", "7=3+4","7=4+3","8=3+5","8=4+4","9=4+5","9=5+4","10=5+5",NA)))
    }
    
    else if(dataset == "Michigan2005"){
      
      probes <- fd_varambally %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
      data <- exp_varambally  %>% filter(ID %in% probes) %>% 
        gather(geo_accession,Expression,-ID)
      fd <- fd_varambally
      pd <- pd_varambally 
    }
    
    else {
      probes <- fd_grasso %>% filter(GENE_SYMBOL %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
      data <- exp_grasso  %>% filter(ID %in% probes) %>% 
        gather(geo_accession,Expression,-ID)
      fd <- mutate(fd_grasso, Symbol = GENE_SYMBOL)
      pd <- mutate(pd_grasso, Sample_Group = Group) 
      
    }
    
    
    summary_stats <- data %>% group_by(ID) %>% 
      summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
    
    data <- left_join(data,summary_stats) %>% mutate(Z = (Expression - mean) / sd)
    
    #    mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
    #    mu <- summary_stats$mean[which.max(summary_stats$iqr)]
    #    sd <- summary_stats$sd[which.max(summary_stats$iqr)]
    
    #    data <- filter(data, ID== mostVarProbe) %>%
    #      mutate(Z = (Expression -mu) /sd)
    
    #    data <- left_join(data,pd)
    #    data <- left_join(data, fd)
    
    
    
    mostVarProbes <- left_join(summary_stats,fd) %>% 
      arrange(Symbol,desc(iqr)) %>% 
      distinct(Symbol) %>% 
      select(ID) %>%  as.matrix %>%  as.character
    
    data <- filter(data, ID %in% mostVarProbes)
    data <- left_join(data, select(fd, ID, Symbol))
    data <- left_join(data, pd)
    data    
    
    data
    
  }
  
  
  
  output$displayBoxplot <- renderPlot({
    
    plotType <- input$inputType
    
    p1 <- singleBoxplot()
    
    if(plotType != "Single Gene") {
      message("Plotting a gene list....")
      inFile <- input$file1
      
      if(!is.null(inFile)) updateTextInput(session, inputId = "profileBasename", value=paste0(basename(inFile$name),"-profile"))
      else updateTextInput(session, inputId = "profileBasename", value="examplegenes-profile")
      
      if(getCambridgeCompositePlot() == "No"){
        
        message("Not making a composite plot,...")
        geneToPlot <- input$cambridgeGeneChoice
        var <- getCambridgeVariable()
        overlay <- getCambridgeOverlay()
        dataset <- getBoxplotDataset()
        data <- filter(p1$data, Symbol == geneToPlot)
        
        if(dataset == "Cambridge"){
          
          p1 <- switch(var,
                       iCluster = {data %>% 
                           filter(Sample_Group == "Tumour",!is.na(iCluster)) %>% 
                           ggplot(aes(x = iCluster, y = Expression, fill=iCluster)) + geom_boxplot() +  scale_fill_manual(values=iclusPal) 
                       },
                       
                       Gleason = {data  %>% 
                           filter(Sample_Group == "Tumour") %>% 
                           filter(!(is.na(Gleason))) %>% 
                           ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot() + scale_fill_manual(values = c("5=3+2"= gradeCols[1], 
                                                                                                                                  "6=3+3"=gradeCols[2],
                                                                                                                                  "6=2+4"= gradeCols[3], 
                                                                                                                                  "7=3+4" = gradeCols[4],
                                                                                                                                  "7=4+3" = gradeCols[5],
                                                                                                                                  "8=3+5" = gradeCols[6],
                                                                                                                                  "8=4+4"=gradeCols[7],
                                                                                                                                  "9=4+5"=gradeCols[8],
                                                                                                                                  "9=5+4"=gradeCols[9],
                                                                                                                                  "10=5+5"=gradeCols[10]))
                       },
                       Sample_Group ={data %>% mutate(Sample_Group = factor(Sample_Group,levels=c("Benign","Tumour","CRPC"))) %>% 
                           ggplot(aes(x = Sample_Group, y = Expression, fill=Sample_Group)) + geom_boxplot() + scale_fill_manual(values = c("darkseagreen1","darkorange1","firebrick1"))
                         
                       }
                       
          )
        }
        
        else if(dataset == "Stockholm"){
          
          p1 <- switch(var,
                       iCluster = {data %>% 
                           ggplot(aes(x = iCluster, y = Expression, fill=iCluster)) + geom_boxplot() +  scale_fill_manual(values=iclusPal) 
                       },
                       
                       Gleason = {data  %>% 
                           filter(!(is.na(Gleason))) %>% 
                           ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot()  + scale_fill_manual(values = c("5=3+2"= gradeCols[1], 
                                                                                                                                   "6=3+3"=gradeCols[2],
                                                                                                                                   "6=2+4"= gradeCols[3], 
                                                                                                                                   "7=3+4" = gradeCols[4],
                                                                                                                                   "7=4+3" = gradeCols[5],
                                                                                                                                   "8=3+5" = gradeCols[6],
                                                                                                                                   "8=4+4"=gradeCols[7],
                                                                                                                                   "9=4+5"=gradeCols[8],
                                                                                                                                   "9=5+4"=gradeCols[9],
                                                                                                                                   "10=5+5"=gradeCols[10]))
                       }
                       
          )
          
        }
        
        else if(dataset == "MSKCC"){
          
          p1 <- switch(var,
                       CopyNumberCluster = {data %>% 
                           filter(!is.na(Copy.number.Cluster)) %>% 
                           ggplot(aes(x = Copy.number.Cluster, y = Expression, fill=Copy.number.Cluster)) + geom_boxplot()
                       },
                       
                       Gleason = {data  %>% 
                           filter(!(is.na(Gleason))) %>% 
                           ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot() + scale_fill_manual(values = c("5=3+2"= gradeCols[1], 
                                                                                                                                  "6=3+3"=gradeCols[2],
                                                                                                                                  "6=2+4"= gradeCols[3], 
                                                                                                                                  "7=3+4" = gradeCols[4],
                                                                                                                                  "7=4+3" = gradeCols[5],
                                                                                                                                  "8=3+5" = gradeCols[6],
                                                                                                                                  "8=4+4"=gradeCols[7],
                                                                                                                                  "9=4+5"=gradeCols[8],
                                                                                                                                  "9=5+4"=gradeCols[9],
                                                                                                                                  "10=5+5"=gradeCols[10]))
                         
                       }
          )
          
        }
        
        else if(dataset == "Michigan2005"){
          
          
          p1 <- data %>% ggplot(aes(x = Sample_Group, y = Expression, fill=Sample_Group)) + geom_boxplot()
        }
        
        else if(dataset == "Michigan2012"){
          
          p1 <- data %>% ggplot(aes(x = Group, y = Expression, fill=Sample_Group)) + geom_boxplot() 
          
        }    
        
        
      }
      p1 <- p1 + facet_wrap(~Symbol)       
      
    }
    else     updateTextInput(session, "profileBasename", value = paste0(input$currentGene,"-profile"))   
    
    message("Going to print the plot..")
    print(p1)
  }
  
  
  )
  
  
  
  singleBoxplot <- reactive({
    
    dataset <- getBoxplotDataset()
    #dataset <- prepareBoxplot()
    
    plotType <- input$inputType
    if(plotType == "Single Gene") {
      genes <- getCurrentGene()
    } else  genes <- getGeneList()
    
    
#    data <- filterByGene(dataset, genes)
    data <- prepareBoxplot()
    
    
    doZ <- ifelse(getCambridgeZ() == "Yes",TRUE,FALSE)
    
    if(doZ) data <- mutate(data, Expression=Z)
    
    var <- getCambridgeVariable()
    overlay <- getCambridgeOverlay()
    
    
    if(dataset == "Cambridge"){
      
      p1 <- switch(var,
                   iCluster = {data %>% 
                       filter(Sample_Group == "Tumour",!is.na(iCluster)) %>% 
                       ggplot(aes(x = iCluster, y = Expression, fill=iCluster)) + geom_boxplot() +  scale_fill_manual(values=iclusPal) 
                   },
                   
                   Gleason = {data  %>% 
                       filter(Sample_Group == "Tumour") %>% 
                       filter(!(is.na(Gleason))) %>% 
                       ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot() + scale_fill_manual(values = c("5=3+2"= gradeCols[1], 
                                                                                                                              "6=3+3"=gradeCols[2],
                                                                                                                              "6=2+4"= gradeCols[3], 
                                                                                                                              "7=3+4" = gradeCols[4],
                                                                                                                              "7=4+3" = gradeCols[5],
                                                                                                                              "8=3+5" = gradeCols[6],
                                                                                                                              "8=4+4"=gradeCols[7],
                                                                                                                              "9=4+5"=gradeCols[8],
                                                                                                                              "9=5+4"=gradeCols[9],
                                                                                                                              "10=5+5"=gradeCols[10]))
                   },
                   Sample_Group ={data %>% mutate(Sample_Group = factor(Sample_Group,levels=c("Benign","Tumour","CRPC"))) %>% 
                       ggplot(aes(x = Sample_Group, y = Expression, fill=Sample_Group)) + geom_boxplot()
                     
                   }
                   
      )
    }
    
    else if(dataset == "Stockholm"){
      
      p1 <- switch(var,
                   iCluster = {data %>% 
                       ggplot(aes(x = iCluster, y = Expression, fill=iCluster)) + geom_boxplot() +  scale_fill_manual(values=iclusPal) 
                   },
                   
                   Gleason = {data  %>% 
                       filter(!(is.na(Gleason))) %>% 
                       ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot() + scale_fill_manual(values = c("5=3+2"= gradeCols[1], 
                                                                                                                              "6=3+3"=gradeCols[2],
                                                                                                                              "6=2+4"= gradeCols[3], 
                                                                                                                              "7=3+4" = gradeCols[4],
                                                                                                                              "7=4+3" = gradeCols[5],
                                                                                                                              "8=3+5" = gradeCols[6],
                                                                                                                              "8=4+4"=gradeCols[7],
                                                                                                                              "9=4+5"=gradeCols[8],
                                                                                                                              "9=5+4"=gradeCols[9],
                                                                                                                              "10=5+5"=gradeCols[10]))
                   }
                   
      )
      
    }
    
    else if(dataset == "MSKCC"){
      
      p1 <- switch(var,
                   CopyNumberCluster = {data %>% 
                       filter(!is.na(Copy.number.Cluster)) %>% 
                       ggplot(aes(x = Copy.number.Cluster, y = Expression, fill=Copy.number.Cluster)) + geom_boxplot()
                   },
                   
                   Gleason = {data  %>% 
                       filter(!(is.na(Gleason))) %>% 
                       ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot()+ scale_fill_manual(values = c("5=3+2"= gradeCols[1], 
                                                                                                                             "6=3+3"=gradeCols[2],
                                                                                                                             "6=2+4"= gradeCols[3], 
                                                                                                                             "7=3+4" = gradeCols[4],
                                                                                                                             "7=4+3" = gradeCols[5],
                                                                                                                             "8=3+5" = gradeCols[6],
                                                                                                                             "8=4+4"=gradeCols[7],
                                                                                                                             "9=4+5"=gradeCols[8],
                                                                                                                             "9=5+4"=gradeCols[9],
                                                                                                                             "10=5+5"=gradeCols[10]))
                   }
      )
      
    }
    
    else if(dataset == "Michigan2005"){
      
      
      p1 <- data %>% ggplot(aes(x = Sample_Group, y = Expression, fill=Sample_Group)) + geom_boxplot()
    }
    
    else if(dataset == "Michigan2012"){
      
      p1 <- data %>% ggplot(aes(x = Group, y = Expression, fill=Group)) + geom_boxplot()
      
    }
    
    
    if(overlay == "Yes")  p1 <- p1 + geom_jitter(position=position_jitter(width = .05),alpha=0.75) 
    
    p1 +  facet_wrap(~Symbol) 
    
    
  })
  
  
  
  
  
  
  output$cambridgeBoxplotPDF <- downloadHandler(
    filename = function(){
      paste0(input$profileBasename,".",input$profilePlotFormat)
    },
    content = function(file) {
      
      if(input$profilePlotFormat == "pdf") pdf(file, width=12, height=6.3)
      else png(file, width=1200,height=600)
      
      plotType <- input$inputType
      if(plotType == "Single Gene") {
        print(singleBoxplot())
      }
      else{
        p1 <- singleBoxplot()
        
        
        if(getCambridgeCompositePlot() == "Yes"){
          
          print(p1)
          
          
        }
        
        else {
          
          doZ <- ifelse(getCambridgeZ() == "Yes",TRUE,FALSE)
          if(doZ) p1$data <- mutate(p1$data, Expression=Z)
          data <- split(p1$data, factor(p1$data$Symbol))
          
          
          
          
          var <- getCambridgeVariable()
          overlay <- getCambridgeOverlay()
          
          pList <- NULL
          dataset <- getBoxplotDataset()
          
          if(dataset == "Cambridge"){
            
            for(i in 1:length(data)){
              pList[[i]] <- switch(var,
                                   iCluster = {data[[i]] %>% 
                                       filter(Sample_Group == "Tumour",!is.na(iCluster)) %>% 
                                       ggplot(aes(x = iCluster, y = Expression, fill=iCluster)) + geom_boxplot() +  scale_fill_manual(values=iclusPal) 
                                   },
                                   
                                   Gleason = {data[[i]]  %>% 
                                       filter(Sample_Group == "Tumour") %>% 
                                       filter(!(is.na(Gleason))) %>% 
                                       ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot() 
                                   },
                                   Sample_Group ={data[[i]] %>% mutate(Sample_Group = factor(Sample_Group,levels=c("Benign","Tumour","CRPC"))) %>% 
                                       ggplot(aes(x = Sample_Group, y = Expression, fill=Sample_Group)) + geom_boxplot()
                                     
                                   }
                                   
              )
              
              if(overlay == "Yes")  pList[[i]] <- pList[[i]] + geom_jitter(position=position_jitter(width = .05),alpha=0.75) 
              
              pList[[i]] <- pList[[i]] +  facet_wrap(~Symbol) 
              
            }        
          }
          
          else if (dataset == "Stockholm"){
            
            for(i in 1:length(data)){
              pList[[i]]  <- switch(var,
                                    iCluster = {data[[i]] %>% 
                                        ggplot(aes(x = iCluster, y = Expression, fill=iCluster)) + geom_boxplot() +  scale_fill_manual(values=iclusPal) 
                                    },
                                    
                                    Gleason = {data[[i]]  %>% 
                                        filter(!(is.na(Gleason))) %>% 
                                        ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot() 
                                    }
                                    
              )
              if(overlay == "Yes")  pList[[i]] <- pList[[i]] + geom_jitter(position=position_jitter(width = .05),alpha=0.75) 
              
              pList[[i]] <- pList[[i]] +  facet_wrap(~Symbol) 
            }
          }
          
          else if (dataset == "MSKCC"){
            for(i in 1:length(data)){
              pList[[i]] <- switch(var,
                                   CopyNumberCluster = {data[[i]] %>% 
                                       filter(!is.na(Copy.number.Cluster)) %>% 
                                       ggplot(aes(x = Copy.number.Cluster, y = Expression, fill=Copy.number.Cluster)) + geom_boxplot()
                                   },
                                   
                                   Gleason = {data[[i]]  %>% 
                                       filter(!(is.na(Gleason))) %>% 
                                       ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot()
                                   }
              )
              if(overlay == "Yes")  pList[[i]] <- pList[[i]] + geom_jitter(position=position_jitter(width = .05),alpha=0.75) 
              
              pList[[i]] <- pList[[i]] +  facet_wrap(~Symbol) 
            }
          }
          else if (dataset == "Michigan2005"){
            for(i in 1:length(data)){
              pList[[i]] <- data[[i]] %>% ggplot(aes(x = Sample_Group, y = Expression, fill=Sample_Group)) + geom_boxplot()
            }
            if(overlay == "Yes")  pList[[i]] <- pList[[i]] + geom_jitter(position=position_jitter(width = .05),alpha=0.75) 
            
            pList[[i]] <- pList[[i]] +  facet_wrap(~Symbol) 
          }  
          
          else if (dataset == "Michigan2012"){
            for(i in 1:length(data)){
              pList[[i]] <- data[[i]] %>% ggplot(aes(x = Group, y = Expression, fill=Group)) + geom_boxplot()
            }
            if(overlay == "Yes")  pList[[i]] <- pList[[i]] + geom_jitter(position=position_jitter(width = .05),alpha=0.75) 
            
            pList[[i]] <- pList[[i]] +  facet_wrap(~Symbol) 
          }
          lapply(pList,print)
        }
        
        
      }
      dev.off()
    }      
    
  )
  
  
  
  output$anovaResult <- renderPrint({
    dataset <- getBoxplotDataset()  
    plotType <- input$inputType
    if(plotType == "Single Gene") {
      genes <- getCurrentGene()
    } else  genes <- getGeneList()
    
#    data <- filterByGene(dataset,genes)
    
    data <- prepareBoxplot()
    
    var <- getCambridgeVariable()
    
    
    if(dataset == "Cambridge"){
      switch(var,
             iCluster = group_by(data, Symbol)  %>% do(tidy(aov(Expression~iCluster,data=.)))%>% filter(term != "Residuals"),
             
             Gleason = group_by(data, Symbol)  %>% do(tidy(aov(Expression~Gleason,data=.)))%>% filter(term != "Residuals"),
             
             Sample_Group = group_by(data, Symbol)  %>% do(tidy(aov(Expression~Sample_Group,data=.))) %>% filter(term != "Residuals")
             
      )
    }
    
    else if (dataset == "Stockholm"){
      
      
      switch(var,
             iCluster = group_by(data, Symbol)  %>% do(tidy(aov(Expression~iCluster,data=.)))%>% filter(term != "Residuals"),
             
             Gleason = group_by(data, Symbol)  %>% do(tidy(aov(Expression~Gleason,data=.)))%>% filter(term != "Residuals")
      )
      
    }
    
    else if (dataset == "MSKCC"){
      
      switch(var,
             CopyNumberCluster = group_by(data, Symbol)  %>% do(tidy(aov(Expression~Copy.number.Cluster,data=.)))%>% filter(term != "Residuals"),
             Gleason = group_by(data, Symbol)  %>% do(tidy(aov(Expression~Gleason,data=.)))%>% filter(term != "Residuals")
      )
      
    }
    
    else{
      group_by(data, Symbol)  %>% do(tidy(aov(Expression~Sample_Group,data=.)))%>% filter(term != "Residuals")
      
    }
    
    
  }
  )
  
  
  prepareSurvival <- reactive({
    
    library(party)
    library(survival)
    
    dataset <- input$rpDataset
    
    plotType <- input$inputType_survival
    if(plotType == "Single Gene") {
      genes <- getCurrentGene()
    }
    else {
      genes<- getGeneList()
      
    }
    
    if(dataset == "Cambridge"){
      
      data <- collect(exp_camcap)  %>% filter(Symbol %in% genes)
      fd <- fd_camcap
      pd <-  mutate(pd_camcap, Gleason=factor(Gleason,levels=c("5=3+2","6=2+4","6=3+3", "7=3+4","7=4+3","8=3+5","8=4+4","9=4+5","9=5+4","10=5+5",NA))) %>% 
        mutate(Time = as.numeric(FollowUpTime), Event = ifelse(BCR=="Y",1,0))
      
    } else if (dataset == "Stockholm"){
      
      data <- collect(exp_stockholm)  %>% filter(Symbol %in% genes)
      
      fd <- fd_stockholm
      pd <-  mutate(pd_stockholm, Gleason=factor(Gleason,levels=c("5=3+2","6=2+4","6=3+3", "7=3+4","7=4+3","8=3+5","8=4+4","9=4+5","9=5+4","10=5+5",NA))) %>% 
        mutate(Time = as.numeric(FollowUpTime), Event = ifelse(BCR=="Y",1,0))
    }
    
    else if(dataset == "MSKCC"){
      
      data <- collect(exp_taylor)  %>% filter(Gene %in% genes)
      
      fd <- fd_taylor %>% mutate(Symbol = Gene)
      pd <- mutate(pd_taylor,Gleason = gsub("4+3", "7=4+3", pd_taylor$Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+3", "8=5+3", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("3+4", "7=3+4", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+3", "7=4+3", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("3+3", "6=3+3", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+5", "9=4+5", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("3+5", "8=3+5", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+4", "8=4+4", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason=factor(Gleason,levels=c("5=3+2","6=2+4","6=3+3", "7=3+4","7=4+3","8=3+5","8=4+4","9=4+5","9=5+4","10=5+5",NA)))
    }
    
 
    
    summary_stats <- data %>% group_by(ID) %>% 
      summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
    
    data <- left_join(data,summary_stats) %>% mutate(Z = (Expression - mean) / sd)
    
    #    mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
    #    mu <- summary_stats$mean[which.max(summary_stats$iqr)]
    #    sd <- summary_stats$sd[which.max(summary_stats$iqr)]
    
    #    data <- filter(data, ID== mostVarProbe) %>%
    #      mutate(Z = (Expression -mu) /sd)
    
    #    data <- left_join(data,pd)
    #    data <- left_join(data, fd)
    
    
    
    mostVarProbes <- left_join(summary_stats,fd) %>% 
      arrange(Symbol,desc(iqr)) %>% 
      distinct(Symbol) %>% 
      select(ID) %>%  as.matrix %>%  as.character
    
    data <- filter(data, ID %in% mostVarProbes)
    data <- left_join(data, select(fd, ID, Symbol))
    data <- left_join(data, pd)

    combined.data <- data %>% 
      filter(!is.na(Time) & !is.na(Event))
    
    
    combined.data
    
  })
  
  #prepareSurvival <- reactive({
  
  #    dataset <- input$rpDataset
  
  #    plotType <- input$inputType_survival
  #    if(plotType == "Single Gene") {
  #      currentGene <- getCurrentGene()
  #    }
  #    else {
  
  #      currentGene <- getGeneList()
  
  #    }
  
  
  
  #  if(dataset == "MSKCC"){
  
  #    probes <- fd_taylor %>% filter(Gene == currentGene) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
  
  #    taylor<- exp_taylor  %>% filter(ID %in% probes) %>% 
  #      gather(geo_accession,Expression,-ID)
  #    
  #    summary_stats <- taylor %>% group_by(ID) %>% 
  #     summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
  
  #    mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
  #    mu <- summary_stats$mean[which.max(summary_stats$iqr)]
  #    sd <- summary_stats$sd[which.max(summary_stats$iqr)]
  
  #    taylor <- filter(taylor, ID== mostVarProbe)
  
  #    taylor <- full_join(taylor,pd_taylor)
  
  
  #    combined.data <- taylor %>% 
  #      filter(!is.na(Event) & !is.na(Time))
  
  #  } else if (dataset == "Cambridge"){
  
  #    camcap <- getSelectedGeneCambridge()
  
  #    camcap <- full_join(camcap,pd_camcap) %>% 
  #      mutate(Time = FollowUpTime, Event = ifelse(BCR=="Y",1,0))
  
  #    combined.data <- camcap %>% 
  #      filter(!is.na(Time) & !is.na(Event))
  
  
  #  } else{
  #   stockholm <- getSelectedGeneStockholm()
  
  #    stockholm <- full_join(stockholm,pd_stockholm) %>% 
  #      mutate(Time = FollowUpTime, Event = ifelse(BCR=="Y",1,0))
  
  #    combined.data <- stockholm %>% 
  #      filter(!is.na(Time) & !is.na(Event))
  
  
  #   }
  
  #    combined.data
  
  # })
  
  output$rpSummary <- renderTable({
    
    combined.data <- prepareSurvival()
    
    plotType <- input$inputType_survival
    if(plotType == "Single Gene") {
      genes <- getCurrentGene()
    }
    else {
      
      genes <- as.character(getGeneList())
      
    }
    message("Gene list is....")
    print(genes)
    pVals <- NULL
    cOffs <- NULL
    
    
    genes <- genes[genes %in% combined.data$Symbol]
    
    
    for(i in 1:length(genes)){
      
      currentGene <- genes[i]
      message("Current Gene is", currentGene)
      
      if(currentGene %in% combined.data$Symbol) message("Gene was found in data")
      else message("Gene was not found in data")
      
      data <- filter(combined.data, Symbol == currentGene)
      message("Entering rpAnalysis.......")
      results <- rpAnalysis(data)
      message("RP analysis done")
      
      ctree_xfs <- results[[1]]
      newPval <- results[[2]]
      
      pVals[i] <- newPval
      
      if(newPval < 0.05){
        
        ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
        ps  <- signif(ps2[1], digits = 3)
        cOffs[i] <- ps  
      } else cOffs[i] <- NA
      
      
    }
    
    df <- data.frame(Gene = genes,RP_p.value=pVals, RP_Cut.off = cOffs)
    
    df
    
    #if(input$cutoffMethod == "RP"){
    
    # summary <- ifelse(newPval <0.05, ". The partitioning and survival curves will apear below", ". Unfortunately, no significant partitioning of the expression values could be found")
    #  paste("Recursive Partitioning was run on ", currentGene, "from the", getRpDataset(), "dataset",summary)
    #  } else if(input$cutoffMehod == "Median") {
    
    #    print("Selecting expression cut-off based on median expression of the gene")
    #  }
    
    #  else print("Using manual cut-off of; ",as.numeric(input$expCutoff))
    
    
  })
  
  
  rpAnalysis <- function(combined.data){
    surv.xfs <- Surv((as.numeric(as.character(combined.data$Time))/12), combined.data$Event)
    print(surv.xfs)
    combined.data$surv.xfs <- surv.xfs
    ctree_xfs   <- ctree(surv.xfs~Expression,data=combined.data)
    pvalue      <- 1 - ctree_xfs@tree$criterion$maxcriterion
    newPval     <- signif(pvalue, digits = 2)
    list(ctree_xfs, newPval,surv.xfs)
    
  }
  
  output$rpPlot <- renderPlot({
    
    
    combined.data <- prepareSurvival()
    
    plotType <- input$inputType_survival
    if(plotType == "Single Gene") {
      currentGene <- getCurrentGene()
    }
    else {
      
      currentGene <- input$survivalGeneChoice
      
    }
    
    message("Gene for RP plot is ", currentGene)
    data <- filter(combined.data, Symbol == currentGene)
    
    #    results <- rpAnalysis(data)
    
    
    #    ctree_xfs <- results[[1]]
    #    newPval <- results[[2]]
    
    #    ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
    #    ps  <- signif(ps2[1], digits = 3)
    
    
    
    #    combined.data <- prepareSurvival()
    #    surv.xfs <- Surv((as.numeric(as.character(combined.data$Time))/12), combined.data$Event)
    #    combined.data$surv.xfs <- surv.xfs
    
    surv.xfs <- Surv((as.numeric(as.character(data$Time))/12), data$Event)
    data$surv.xfs <- surv.xfs
    if(input$cutoffMethod == "RP"){
      results <- rpAnalysis(data)
      ctree_xfs <- results[[1]]
      newPval <- results[[2]]
      
      if(newPval<0.05) {
        
        
        ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
        ps  <- signif(ps2[1], digits = 3)
        plot(ctree(surv.xfs~Expression, data=data))
      }
      
      else {
        med <- median(data$Expression)
        p <- ggplot(data, aes(x=Expression)) + geom_histogram() + geom_vline(xintercept=med,col="red",lty=2) + ggtitle("Partitioning using median expression")
        print(p)
      }
      
    }
    
    else{
      
      if (input$cutoffMethod == "Median") {
        
        med <- median(combined.data$Expression)
        p <- ggplot(combined.data, aes(x=Expression)) + geom_histogram() + geom_vline(xintercept=med,col="red",lty=2)
        print(p)
      }
      else  ggplot(combined.data, aes(x=Expression)) + geom_histogram() + geom_vline(xintercept=as.numeric(input$expCutOff),col="red",lty=2) + ggtitle("Partitioning using defined cut-off")
      
      
    }
    
    
    
    
  }
  )
  
  
  
  output$survivalPlot <- renderPlot({
    combined.data <- prepareSurvival()
    
    plotType <- input$inputType_survival
    if(plotType == "Single Gene") {
      currentGene <- getCurrentGene()
    }
    else {
      
      currentGene <- input$survivalGeneChoice
      
    }
    
    updateTextInput(session, "survivalBasename",value=paste0(currentGene,"-survival"))
    
    data <- filter(combined.data, Symbol == currentGene)
    updateTextInput(session, "expCutOff",value=round(median(data$Expression),2))
    surv.xfs <- Surv((as.numeric(as.character(data$Time))/12), data$Event)
    
    if(input$cutoffMethod == "RP"){
      
      results <- rpAnalysis(data)
      
      data$surv.xfs <- surv.xfs
      ctree_xfs <- results[[1]]
      newPval <- results[[2]]
      
      if(newPval<0.05) {
        
        
        ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
        ps  <- signif(ps2[1], digits = 3)
        
        
        if(length(ps2)==1) {
          data$geneexp_cp <- data$Expression<=ps2[1]
          nt                       <- table(data$geneexp_cp)
          geneexp.survfit.xfs      <- survfit(surv.xfs~geneexp_cp,data=data)
          plot(geneexp.survfit.xfs, xlab="Time to BCR (years)", ylab="Probability of Freedom from Biochemical Recurrence", main=paste(currentGene,", p=", newPval), col=c(2,4))
          legend("bottomleft", c(paste(currentGene, ">", ps, "n=", nt[[1]]), paste(currentGene, "<=", ps, "n=", nt[[2]])), col=c(2,4), lty=1, lwd=1.5, bty="n")
          newPval2                 <- NA
        }
        
      }
      
      else{
        ps <- round(median(data$Expression),3)
        data$geneexp_cp <- data$Expression<= ps
        nt                       <- table(data$geneexp_cp)
        geneexp.survfit.xfs      <- survfit(surv.xfs~geneexp_cp,data=data)
        
        test <- survdiff(surv.xfs~geneexp_cp,data=data)
        
        newPval <- round(pchisq(test$chisq, df = length(test$n)-1, lower.tail=FALSE),3)
        
        plot(geneexp.survfit.xfs, xlab="Time to BCR (years)", ylab="Probability of Freedom from Biochemical Recurrence", main=paste(currentGene,", p=", newPval), col=c(2,4))
        legend("bottomleft", c(paste(currentGene, ">", ps, "n=", nt[[1]]), paste(currentGene, "<=", ps, "n=", nt[[2]])), col=c(2,4), lty=1, lwd=1.5, bty="n")
        
      }
      
    }
    
    else{
      
      if (input$cutoffMethod == "Median"){
        
        
        ps <- round(median(data$Expression),3)
      } else ps <- as.numeric(input$expCutOff)
      
      data$geneexp_cp <- data$Expression<= ps
      nt                       <- table(data$geneexp_cp)
      geneexp.survfit.xfs      <- survfit(surv.xfs~geneexp_cp,data=data)
      
      test <- survdiff(surv.xfs~geneexp_cp,data=data)
      
      newPval <- round(pchisq(test$chisq, df = length(test$n)-1, lower.tail=FALSE),3)
      
      plot(geneexp.survfit.xfs, xlab="Time to BCR (years)", ylab="Probability of Freedom from Biochemical Recurrence", main=paste(currentGene,", p=", newPval), col=c(2,4))
      legend("bottomleft", c(paste(currentGene, ">", ps, "n=", nt[[1]]), paste(currentGene, "<=", ps, "n=", nt[[2]])), col=c(2,4), lty=1, lwd=1.5, bty="n")
      
      
    }
    
  })
  
  
  
  output$survivalPlotPDF <- downloadHandler(
    filename = function(){
      paste0(input$survivalBasename, "." ,input$survivalPlotFormat)
    },
    content = function(file) {
      
      if(input$survivalPlotFormat == "pdf") pdf(file, width=12, height=6.3)
      else png(file, width=1200,height=600)
      
      plotType <- input$inputType_survival
      if(plotType == "Single Gene") {
        currentGene <- getCurrentGene()
      }
      else {
        
        currentGene <- input$survivalGeneChoice
        
      }
      combined.data <- prepareSurvival()
      
      
      data <- filter(combined.data, Symbol == currentGene)
      updateTextInput(session, "expCutOff",value=round(median(data$Expression),2))
      surv.xfs <- Surv((as.numeric(as.character(data$Time))/12), data$Event)
      
      if(input$cutoffMethod == "RP"){
        
        results <- rpAnalysis(data)
        
        data$surv.xfs <- surv.xfs
        ctree_xfs <- results[[1]]
        newPval <- results[[2]]
        
        if(newPval<0.05) {
          
          
          ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
          ps  <- signif(ps2[1], digits = 3)
          
          
          if(length(ps2)==1) {
            data$geneexp_cp <- data$Expression<=ps2[1]
            nt                       <- table(data$geneexp_cp)
            geneexp.survfit.xfs      <- survfit(surv.xfs~geneexp_cp,data=data)
            plot(geneexp.survfit.xfs, xlab="Time to BCR (years)", ylab="Probability of Freedom from Biochemical Recurrence", main=paste(currentGene,", p=", newPval), col=c(2,4))
            legend("bottomleft", c(paste(currentGene, ">", ps, "n=", nt[[1]]), paste(currentGene, "<=", ps, "n=", nt[[2]])), col=c(2,4), lty=1, lwd=1.5, bty="n")
            newPval2                 <- NA
          }
          
        }
        
        else{
          ps <- round(median(data$Expression),3)
          data$geneexp_cp <- data$Expression<= ps
          nt                       <- table(data$geneexp_cp)
          geneexp.survfit.xfs      <- survfit(surv.xfs~geneexp_cp,data=data)
          
          test <- survdiff(surv.xfs~geneexp_cp,data=data)
          
          newPval <- round(pchisq(test$chisq, df = length(test$n)-1, lower.tail=FALSE),3)
          
          plot(geneexp.survfit.xfs, xlab="Time to BCR (years)", ylab="Probability of Freedom from Biochemical Recurrence", main=paste(currentGene,", p=", newPval), col=c(2,4))
          legend("bottomleft", c(paste(currentGene, ">", ps, "n=", nt[[1]]), paste(currentGene, "<=", ps, "n=", nt[[2]])), col=c(2,4), lty=1, lwd=1.5, bty="n")
          
        }
        
      }
      
      else{
        
        if (input$cutoffMethod == "Median"){
          
          
          ps <- round(median(data$Expression),3)
        } else ps <- as.numeric(input$expCutOff)
        
        data$geneexp_cp <- data$Expression<= ps
        nt                       <- table(data$geneexp_cp)
        geneexp.survfit.xfs      <- survfit(surv.xfs~geneexp_cp,data=data)
        
        test <- survdiff(surv.xfs~geneexp_cp,data=data)
        
        newPval <- round(pchisq(test$chisq, df = length(test$n)-1, lower.tail=FALSE),3)
        
        plot(geneexp.survfit.xfs, xlab="Time to BCR (years)", ylab="Probability of Freedom from Biochemical Recurrence", main=paste(currentGene,", p=", newPval), col=c(2,4))
        legend("bottomleft", c(paste(currentGene, ">", ps, "n=", nt[[1]]), paste(currentGene, "<=", ps, "n=", nt[[2]])), col=c(2,4), lty=1, lwd=1.5, bty="n")
        
        
      }
      
      
      dev.off()
    }
    
  )
  
  prepareHeatmap <- reactive({
    
    
    genes <- getGeneList()
    dataset <- getHeatmapDataset()
    
    if(dataset == "MSKCC"){
      
      taylor <- collect(exp_taylor)  %>% filter(Gene %in% genes)
      
      summary_stats <- taylor %>% group_by(ID) %>% 
        summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
      
      mostVarProbes <- left_join(summary_stats,fd_taylor) %>% 
        arrange(Gene,desc(iqr)) %>% 
        distinct(Gene) %>% 
        select(ID) %>%  as.matrix %>%  as.character
      
      
      
      samples <- filter(pd_taylor, Sample_Group == "prostate cancer") %>% 
        select(geo_accession) %>% as.matrix %>% as.character
      
      
      
      taylor <- filter(taylor, ID %in% mostVarProbes,geo_accession %in% samples)
      geneMatrix <- taylor %>% 
        spread(geo_accession,Expression) %>% select(-Gene) %>% data.frame
      
      geneMatrix <- as.matrix(geneMatrix[,-1])
      
      
      symbols <- filter(fd_taylor, ID %in% mostVarProbes) %>% 
        select(Gene) %>% as.matrix %>% as.character
      
      rownames(geneMatrix) <- symbols
      
      pd <- left_join(taylor,pd_taylor) %>% distinct(geo_accession)
      
      colMatrix <- matrix(nrow = ncol(geneMatrix),ncol = 2)
      grp <- pd$Copy.number.Cluster
      cols <- brewer.pal(9,"Set1")[1:length(levels(factor(as.character(grp))))]
      
      grp <- as.factor(as.character(grp))
      levels(grp) <- cols
      
      colMatrix[,1] <- as.character(grp)
      
      grp <- pd$Gleason
      cols <- brewer.pal(8,"Set2")[1:length(levels(factor(as.character(grp))))]
      
      grp <- as.factor(as.character(grp))
      levels(grp) <- cols
      colMatrix[,2] <- as.character(grp)
      colnames(colMatrix) <- c("Copy Number Cluster", "Gleason")
      
      
    } else if(dataset == "Cambridge"){
      
      camcap <- collect(exp_camcap)  %>% filter(Symbol %in% genes)
      
      summary_stats <- camcap %>% group_by(ID) %>% 
        summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
      
      mostVarProbes <- left_join(summary_stats,fd_camcap) %>% 
        arrange(Symbol,desc(iqr)) %>% 
        distinct(Symbol) %>% 
        select(ID) %>%  as.matrix %>%  as.character
      
      
      
      
      samples <- filter(pd_camcap, Sample_Group == "Tumour") %>% 
        select(geo_accession) %>% as.matrix %>% as.character
      
      camcap <- filter(camcap, ID %in% mostVarProbes,geo_accession %in% samples) %>% select(-Symbol)
      
      geneMatrix <- camcap %>% 
        spread(geo_accession,Expression) %>% data.frame
      
      geneMatrix <- as.matrix(geneMatrix[,-1])
      
      symbols <- filter(fd_camcap, ID %in% mostVarProbes) %>% 
        select(Symbol) %>% as.matrix %>% as.character
      
      rownames(geneMatrix) <- symbols
      
      
      pd <- left_join(camcap,pd_camcap) %>% distinct(geo_accession)
      
      colMatrix <- matrix(nrow = ncol(geneMatrix),ncol = 2)
      grp <- pd$iCluster
      cols <- brewer.pal(5,"Set1")[1:length(levels(factor(as.character(grp))))]
      
      grp <- as.factor(as.character(grp))
      levels(grp) <- cols
      
      colMatrix[,1] <- as.character(grp)
      
      grp <- pd$Gleason
      cols <- brewer.pal(8,"Set2")[1:length(levels(factor(as.character(grp))))]
      
      grp <- as.factor(as.character(grp))
      levels(grp) <- cols
      colMatrix[,2] <- as.character(grp)
      
      colnames(colMatrix) <- c("iCluster", "Gleason")
      
      
      
      
      
    } else{
      
      stockholm <- collect(exp_stockholm)  %>% filter(Symbol %in% genes)
      
      summary_stats <- stockholm %>% group_by(ID) %>% 
        summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
      
      mostVarProbes <- left_join(summary_stats,fd_stockholm) %>% 
        arrange(Symbol,desc(iqr)) %>% 
        distinct(Symbol) %>% 
        select(ID) %>%  as.matrix %>%  as.character
      
      
      
      
      samples <-  select(pd_stockholm,geo_accession) %>% as.matrix %>% as.character
      
      stockholm <- filter(stockholm, ID %in% mostVarProbes,geo_accession %in% samples) 
      
      geneMatrix <- stockholm %>% 
        spread(geo_accession,Expression) %>% select(-Symbol) %>% data.frame
      
      geneMatrix <- as.matrix(geneMatrix[,-1])
      
      symbols <- filter(fd_stockholm, ID %in% mostVarProbes) %>% 
        select(Symbol) %>% as.matrix %>% as.character
      
      rownames(geneMatrix) <- symbols
      
      
      pd <- left_join(stockholm,pd_stockholm) %>% distinct(geo_accession)
      
      colMatrix <- matrix(nrow = ncol(geneMatrix),ncol = 2)
      grp <- pd$iCluster
      cols <- brewer.pal(5,"Set1")[1:length(levels(factor(as.character(grp))))]
      
      grp <- as.factor(as.character(grp))
      levels(grp) <- cols
      
      colMatrix[,1] <- as.character(grp)
      
      grp <- pd$Gleason
      cols <- brewer.pal(8,"Set2")[1:length(levels(factor(as.character(grp))))]
      
      grp <- as.factor(as.character(grp))
      levels(grp) <- cols
      colMatrix[,2] <- as.character(grp)
      colnames(colMatrix) <- c("iCluster", "Gleason")
      
    }
    
    return(list(geneMatrix, colMatrix))
  })
  
  
  
  doClustering <- reactive({
    
    hm <- prepareHeatmap()
    geneMatrix <- hm[[1]]
    
    if(getDistFun() == "Correlation") distfun <- function(x) as.dist(1 - cor(t(x)))
    else distfun <- dist
    
    hclustfun <- function(x) hclust(x,method=getHclustMethod())
    
    
    dMat <- as.dist(distfun(t(geneMatrix)))
    
    clusObj <- hclustfun(dMat)
    
    clusObj
    
  })
  
  output$heatmap<- renderPlot({
    
    hm <- prepareHeatmap()
    geneMatrix <- hm[[1]]
    
    clusObj <- doClustering()
    
    hmcol <- rev(brewer.pal(11 , "RdBu"))
    
    
    if(getDistFun() == "Correlation") distfun <- function(x) as.dist(1 - cor(t(x)))
    else distfun <- dist
    
    hclustfun <- function(x) hclust(x,method=getHclustMethod())
    
    
    scale <- getScaleMethod()
    
    if(getReordRows() == "Yes") heatmap.2(geneMatrix,Colv = as.dendrogram(clusObj),col=hmcol,distfun=distfun,hclustfun = hclustfun,scale=scale,trace="none",cexRow = 0.9)
    else heatmap.2(geneMatrix,Colv = as.dendrogram(clusObj),col=hmcol,distfun=distfun,Rowv=NA,hclustfun = hclustfun,scale=scale,trace="none",cexRow = 0.9)
  }
  )
  
  output$HeatmapPDF <- downloadHandler(
    filename = function(){
      paste0(input$heatmapBasename,".",input$heatmapPlotFormat)
    },
    content = function(file) {
      
      if(input$heatmapPlotFormat == "pdf") pdf(file, width=12, height=6.3)
      else png(file, width=1200,height=600)
      
      hm <- prepareHeatmap()
      geneMatrix <- hm[[1]]
      
      clusObj <- doClustering()
      
      hmcol <- rev(brewer.pal(11 , "RdBu"))
      
      
      if(getDistFun() == "Correlation") distfun <- function(x) as.dist(1 - cor(t(x)))
      else distfun <- dist
      
      hclustfun <- function(x) hclust(x,method=getHclustMethod())
      
      
      scale <- getScaleMethod()
      
      if(getReordRows() == "Yes") heatmap.2(geneMatrix,Colv = as.dendrogram(clusObj),col=hmcol,distfun=distfun,hclustfun = hclustfun,scale=scale,trace="none",cexRow = 0.9)
      else heatmap.2(geneMatrix,Colv = as.dendrogram(clusObj),col=hmcol,distfun=distfun,Rowv=NA,hclustfun = hclustfun,scale=scale,trace="none",cexRow = 0.9)
      
      dev.off()   
      
    }
    
  )
  
  
  
  output$dendrogram <- renderPlot({
    library(WGCNA)
    hm <- prepareHeatmap()
    colMatrix <- hm[[2]]
    clusObj <- doClustering()
    
    
    h <- min(as.numeric(input$hCut),max(clusObj$height))
    k <- input$kGrps
    
    if(input$cutType == "k") newLabs <- cutree(clusObj,k=k)
    else newLabs <-  cutree(clusObj,h=h)
    
    newLabs <- factor(newLabs)
    levels(newLabs) <- rainbow(n=length(unique(newLabs)))
    
    colMatrix <- data.frame(colMatrix, SampleCluster=as.character(newLabs))
    
    plotDendroAndColors.mod(clusObj,colors = as.matrix(colMatrix),abHeight = h,k=k,cutType=input$cutType)    
  }
  
  )
  
  
  
  output$sampleBreakown <- renderPlot({
    
    hm <- prepareHeatmap()
    colMatrix <- hm[[2]]
    clusObj <- doClustering()
    
    dataset <- getHeatmapDataset()
    
    if(input$cutType == "k"){
      kGrps <- cutree(clusObj,k=input$kGrps)
    } else {
      h=min(as.numeric(input$hCut),max(clusObj$height))
      message("cutting at height",h)
      kGrps <- cutree(clusObj,h=h)
    }
    newGrps <- tbl_df(data.frame(geo_accession=names(kGrps), Cluster = factor(kGrps)))
    
    
    if(dataset == "Cambridge"){
      
      new_pheno <- left_join(pd_camcap,newGrps) %>% filter(!is.na(Cluster))
      p0 <- ggplot(new_pheno, aes(x = Cluster,fill=Cluster)) + geom_bar() +  scale_fill_manual(values=as.character(rainbow(n=length(unique(kGrps))))) + coord_flip()
      p1 <- ggplot(new_pheno,aes(x=iCluster,fill=iCluster)) + geom_bar() + facet_wrap(~Cluster,nrow=1) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position="none")  +  scale_fill_manual(values=iclusPal) 
      p2 <- ggplot(new_pheno,aes(x=Gleason,fill=Gleason)) + geom_bar() + facet_wrap(~Cluster,nrow=1) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position="none") 
      grid.arrange(p0,p1,p2)
      
    } else if(dataset == "Stockholm"){
      
      new_pheno <- left_join(pd_stockholm,newGrps) %>% filter(!is.na(Cluster))
      p0 <- ggplot(new_pheno, aes(x = Cluster,fill=Cluster)) + geom_bar() +  scale_fill_manual(values=as.character(rainbow(n=length(unique(kGrps))))) + coord_flip()
      p1 <- ggplot(new_pheno,aes(x=iCluster,fill=iCluster)) + geom_bar() + facet_wrap(~Cluster,nrow=1) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position="none")  +  scale_fill_manual(values=iclusPal) 
      p2 <- ggplot(new_pheno,aes(x=Gleason,fill=Gleason)) + geom_bar() + facet_wrap(~Cluster,nrow=1) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position="none") 
      grid.arrange(p0,p1,p2)
      
    }
    
    else{
      
      new_pheno <- left_join(pd_taylor,newGrps) %>% filter(!is.na(Cluster))
      p0 <- ggplot(new_pheno, aes(x = Cluster,fill=Cluster)) + geom_bar() +  scale_fill_manual(values=as.character(rainbow(n=length(unique(kGrps))))) + coord_flip()
      p1 <- ggplot(new_pheno,aes(x=Copy.number.Cluster,fill=Copy.number.Cluster)) + geom_bar() + facet_wrap(~Cluster,nrow=1) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position="none")  
      p2 <- ggplot(new_pheno,aes(x=Gleason,fill=Gleason)) + geom_bar() + facet_wrap(~Cluster,nrow=1) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position="none") 
      grid.arrange(p0,p1,p2)
      
    }
    
    
  }
  )
  
  getCopyNumberTable <- reactive({
    
    if(input$inputType_cn == "Single Gene") {
      genes <- getCurrentGene()
      
      
    } else  genes <- getGeneList()
    
    
    camcap.cn <- collect(tbl("copyNumber",src=db_camcap))
    stockholm.cn <- collect(tbl("copyNumber",src=db_stockholm))
    taylor.cn <- collect(tbl("copyNumber",src=db_taylor))
    
    camcap.cnts <- filter(camcap.cn, Symbol %in% genes) %>% 
      group_by(Symbol) %>% 
      summarise(NEUTRAL = 100*sum(Call==0)/n(),DEL = -100*sum(Call==-1)/n(), AMP = 100*sum(Call==1)/n())  %>% 
      gather(Event,Percentage,-Symbol) %>% 
      mutate(Cohort = "Cambridge")
    
    stockholm.cnts <- filter(stockholm.cn, Symbol %in% genes) %>% 
      group_by(Symbol) %>% 
      summarise(NEUTRAL = 100*sum(Call==0)/n(),DEL = -100*sum(Call==-1)/n(), AMP = 100*sum(Call==1)/n())  %>% 
      gather(Event,Percentage,-Symbol) %>% 
      mutate(Cohort = "Stockholm")
    
    taylor.cnts <- filter(taylor.cn, Symbol %in% genes) %>% 
      group_by(Symbol) %>% 
      summarise(NEUTRAL = 100*sum(Call==0)/n(),DEL = -100*sum(Call==-1)/n(), AMP = 100*sum(Call==1)/n())  %>% 
      gather(Event,Percentage,-Symbol) %>% 
      mutate(Cohort = "MSKCC")
    
    cn.all <- bind_rows(camcap.cnts,stockholm.cnts,taylor.cnts) %>% 
      mutate(Event = factor(Event, levels=c("DEL","NEUTRAL","AMP")))
    
    cn.all
  }
  
  )
  
  
  output$copyNumberTable <- renderTable({
    
    cn.all <- getCopyNumberTable()
    
    cn.all
    
  })
  
  
  output$copyNumber <- renderPlot({
    
    cn.all <- getCopyNumberTable()
    
    p <- ggplot(cn.all, aes(x = Event, y=Percentage,fill=Event)) + geom_bar(stat="identity") + facet_wrap(Symbol~Cohort) + scale_fill_manual(values=c("dodgerblue4", "grey","firebrick3"))
    print(p)
  })
  
  
  output$copyNumberPDF <- downloadHandler(
    filename = function(){
      paste0(input$copyNumberBasename,".",input$copyNumberPlotFormat)
    },
    content = function(file) {
      
      if(input$correlationPlotFormat == "pdf") pdf(file, width=12, height=6.3)
      else png(file, width=1200,height=600)
      
      if(input$inputType_cn == "Single Gene") {
        genes <- getCurrentGene()
        
        
      } else  genes <- getGeneList()
      
      camcap.cnts <- filter(camcap.cn, Symbol %in% genes) %>% 
        group_by(Symbol) %>% 
        summarise(NEUTRAL = 100*sum(Call==0)/n(),DEL = -100*sum(Call==-1)/n(), AMP = 100*sum(Call==1)/n())  %>% 
        gather(Event,Percentage,-Symbol) %>% 
        
        mutate(Cohort = "Cambridge")
      
      stockholm.cnts <- filter(stockholm.cn, Symbol %in% genes) %>% 
        group_by(Symbol) %>% 
        summarise(NEUTRAL = 100*sum(Call==0)/n(),DEL = -100*sum(Call==-1)/n(), AMP = 100*sum(Call==1)/n())  %>% 
        gather(Event,Percentage,-Symbol) %>% 
        
        mutate(Cohort = "Stockholm")
      
      taylor.cnts <- filter(taylor.cn, Symbol %in% genes) %>% 
        group_by(Symbol) %>% 
        summarise(NEUTRAL = 100*sum(Call==0)/n(),DEL = -100*sum(Call==-1)/n(), AMP = 100*sum(Call==1)/n())  %>% 
        gather(Event,Percentage,-Symbol) %>% 
        
        mutate(Cohort = "MSKCC")
      
      cn.all <- bind_rows(camcap.cnts,stockholm.cnts,taylor.cnts) %>% 
        mutate(Event = factor(Event, levels=c("DEL","NEUTRAL","AMP")))
      
      cn.all
      
      p <- ggplot(cn.all, aes(x = Event, y=Count,fill=Event)) + geom_bar(stat="identity") + facet_wrap(Symbol~Cohort) + scale_fill_manual(values=c("dodgerblue4", "grey","firebrick3"))
      print(p)
      dev.off()
      
    }
    
  )
  
  
  prepareCorrelation <- reactive({
    
    plotType <- input$inputType_correlation
    
    if(plotType == "Single Gene") {
      gene1 <- getCurrentGene()
      gene2 <- input$secondGene
      genes <- c(gene1,gene2)
      
    } else  genes <- getGeneList()
    
    
    dataset <- getCorDataset()
    
    if(dataset == "Cambridge"){
      
      data <- collect(exp_camcap)  %>% filter(Symbol %in% genes)
      fd <- fd_camcap
      pd <-  mutate(pd_camcap, Gleason=factor(Gleason,levels=c("5=3+2","6=2+4","6=3+3", "7=3+4","7=4+3","8=3+5","8=4+4","9=4+5","9=5+4","10=5+5",NA))) %>% 
        mutate(Time = as.numeric(FollowUpTime), Event = ifelse(BCR=="Y",1,0))
      
    } else if (dataset == "Stockholm"){
      
      data <- collect(exp_stockholm)  %>% filter(Symbol %in% genes)
      fd <- fd_stockholm
      pd <-  mutate(pd_stockholm, Gleason=factor(Gleason,levels=c("5=3+2","6=2+4","6=3+3", "7=3+4","7=4+3","8=3+5","8=4+4","9=4+5","9=5+4","10=5+5",NA))) %>% 
        mutate(Time = as.numeric(FollowUpTime), Event = ifelse(BCR=="Y",1,0))
    }
    
    else if(dataset == "MSKCC"){
      
      data <- collect(exp_taylor)  %>% filter(Gene %in% genes)
      fd <- fd_taylor %>% mutate(Symbol = Gene)
      pd <- mutate(pd_taylor,Gleason = gsub("4+3", "7=4+3", pd_taylor$Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+3", "8=5+3", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("3+4", "7=3+4", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+3", "7=4+3", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("3+3", "6=3+3", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+5", "9=4+5", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("3+5", "8=3+5", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason = gsub("4+4", "8=4+4", Gleason,fixed=TRUE)) %>% 
        mutate(Gleason=factor(Gleason,levels=c("5=3+2","6=2+4","6=3+3", "7=3+4","7=4+3","8=3+5","8=4+4","9=4+5","9=5+4","10=5+5",NA)))
    }
    
    else if(dataset == "Michigan2005"){
      
      data <- collect(exp_varambally)  %>% filter(Symbol %in% genes)
      fd <- fd_varambally
      pd <- pd_varambally 
    }
    
    else {
      data <- collect(exp_grasso)  %>% filter(GENE_SYMBOL %in% genes)
      fd <- mutate(fd_grasso, Symbol = GENE_SYMBOL)
      pd <- mutate(pd_grasso, Sample_Group = Group) 
      
    }
    
    
    summary_stats <- data %>% group_by(ID) %>% 
      summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
    
    data <- left_join(data,summary_stats) %>% mutate(Z = (Expression - mean) / sd)
    
    #    mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
    #    mu <- summary_stats$mean[which.max(summary_stats$iqr)]
    #    sd <- summary_stats$sd[which.max(summary_stats$iqr)]
    
    #    data <- filter(data, ID== mostVarProbe) %>%
    #      mutate(Z = (Expression -mu) /sd)
    
    #    data <- left_join(data,pd)
    #    data <- left_join(data, fd)
    
    
    
    mostVarProbes <- left_join(summary_stats,fd) %>% 
      arrange(Symbol,desc(iqr)) %>% 
      distinct(Symbol) %>% 
      select(ID) %>%  as.matrix %>%  as.character
    
    data <- filter(data, ID %in% mostVarProbes)
    data <- left_join(data, select(fd, ID, Symbol))
    data <- left_join(data, pd)
    
  })
  
  
  output$corPlot <- renderPlot({
    
    plotType <- input$inputType_correlation
    
    if(plotType == "Single Gene") {
      gene1 <- getCurrentGene()
      gene2 <- input$secondGene
      genes <- c(gene1,gene2)
      
    } else  genes <- getGeneList()
    
    
    dataset <- getCorDataset()
    
    #cordata <- filterByGene(dataset, genes)
    
    cordata <- prepareCorrelation()
    
    covar <- input$clinvar_cor
    
    if(dataset == "Cambridge"){
      
      if(covar == "iCluster"){
        
        data <- select(cordata, geo_accession,Expression, iCluster, Symbol) %>% 
          filter(iCluster %in% c("clust1","clust2","clust3","clust4","clust5")) %>% 
          #        mutate(Symbol=gsub(genes[1], "Gene1",Symbol)) %>% 
          #       mutate(Symbol=gsub(genes[2], "Gene2",Symbol)) %>% 
          spread(Symbol,Expression)
        #      cor <- round(cor(data$Gene1,data$Gene2,method=getCorType()),3)
        df <- as.data.frame(data)
        p <- ggpairs(df, columns = 3:ncol(df),col="iCluster")
        
        #        ggplot(data, aes(x = Gene1, y=Gene2,col=iCluster)) + geom_point() + xlab(genes[1]) + ylab(genes[2]) +  scale_fill_manual(values=iclusPal) + ggtitle(paste("Correlation = ",cor))
        
      }
      
      else if(covar == "Gleason"){
        
        data <- select(cordata, geo_accession,Expression, Gleason, Symbol) %>% 
          spread(Symbol,Expression)
        df <- as.data.frame(data)
        p <- ggpairs(df, columns = 3:ncol(df),col="Gleason")
        
      }
      
      else if(covar == "Sample_Group"){
        
        data <- select(cordata, geo_accession,Expression, Sample_Group, Symbol) %>% 
          spread(Symbol,Expression)
        df <- as.data.frame(data)
        p <- ggpairs(df, columns = 3:ncol(df),col="Sample_Group")
        
      }
      
      else {
        
        data <- select(cordata, geo_accession,Expression, Gleason, Symbol) %>% 
          spread(Symbol,Expression)
        df <- as.data.frame(data)
        p <- ggpairs(df, columns = 3:ncol(df))
      }
      
      
      
    }
    
    else if (dataset == "Stockholm"){
      
      if(covar == "iCluster"){
        
        data <- select(cordata, geo_accession,Expression, iCluster, Symbol) %>% 
          filter(iCluster %in% c("clust1","clust2","clust3","clust4","clust5")) %>% 
          #        mutate(Symbol=gsub(genes[1], "Gene1",Symbol)) %>% 
          #       mutate(Symbol=gsub(genes[2], "Gene2",Symbol)) %>% 
          spread(Symbol,Expression)
        #      cor <- round(cor(data$Gene1,data$Gene2,method=getCorType()),3)
        df <- as.data.frame(data)
        p <- ggpairs(df, columns = 3:ncol(df),col="iCluster")
        
      }
      
      else if(covar == "Gleason"){
        
        data <- select(cordata, geo_accession,Expression, Gleason, Symbol) %>% 
          spread(Symbol,Expression)
        df <- as.data.frame(data)
        p <- ggpairs(df, columns = 3:ncol(df),col="Gleason")
        
      }
      
      else {
        
        data <- select(cordata, geo_accession,Expression, Gleason, Symbol) %>% 
          spread(Symbol,Expression)
        df <- as.data.frame(data)
        p <- ggpairs(df, columns = 3:ncol(df))
        
      }
      
      
    }
    
    else if (dataset=="MSKCC"){
      
      if(covar == "CopyNumberCluster"){
        
        data <- select(cordata, geo_accession,Expression, Copy.number.Cluster, Symbol) %>% 
          spread(Symbol,Expression)
        df <- as.data.frame(data)
        cor <- round(cor(data$Gene1,data$Gene2,method=getCorType()),3)
        p <- ggpairs(df, columns = 3:ncol(df),col="Copy.number.Cluster")
      }
      else if(covar == "Gleason"){
        
        
        data <- select(cordata, geo_accession,Expression, Gleason, Symbol) %>% 
          spread(Symbol,Expression)
        df <- as.data.frame(data)
        p <- ggpairs(df, columns = 3:ncol(df),col="Gleason")
      }
      
      else {
        
        data <- select(cordata, geo_accession,Expression, Gleason, Symbol) %>% 
          spread(Symbol,Expression)
        df <- as.data.frame(data)
        p <- ggpairs(df, columns = 3:ncol(df))
        
      }
      
      
    }
    
    else if (dataset == "Michigan2005"){
      
      if(covar == "Sample_Group"){
        data <- select(cordata, geo_accession,Expression, Sample_Group, Symbol) %>% 
          spread(Symbol,Expression)
        df <- as.data.frame(data)
        p <- ggpairs(df, columns = 3:ncol(df),col="Sample_Group")
        
      }
      else{
        data <- select(cordata, geo_accession,Expression, Sample_Group, Symbol) %>% 
          spread(Symbol,Expression)
        df <- as.data.frame(data)
        p <- ggpairs(df, columns = 3:ncol(df))
        
      }
      
    }
    
    else if (dataset == "Michigan2012"){
      
      if(covar == "Sample_Group"){
        data <- select(cordata, geo_accession,Expression, Sample_Group, Symbol) %>% 
          spread(Symbol,Expression)
        df <- as.data.frame(data)
        p <- ggpairs(df, columns = 3:ncol(df),col="Sample_Group")
        
      }
      else{
        data <- select(cordata, geo_accession,Expression, Sample_Group, Symbol) %>% 
          spread(Symbol,Expression)
        df <- as.data.frame(data)
        p <- ggpairs(df, columns = 3:ncol(df))
        
      }
      
    }
    
    print(p)
  }
  
  )
  
  
  
  
  
  
  output$CorrelationPDF <- downloadHandler(
    filename = function(){
      paste0(input$correlationBasename,".",input$correlationPlotFormat)
    },
    content = function(file) {
      
      if(input$correlationPlotFormat == "pdf") pdf(file, width=12, height=6.3)
      else png(file, width=1200,height=600)
      
      plotType <- input$inputType_correlation
      if(plotType == "Single Gene") {
        gene1 <- getCurrentGene()
        gene2 <- input$secondGene
        genes <- c(gene1,gene2)
        
      } else  genes <- getGeneList()
      
      
      dataset <- getCorDataset()
      
      #cordata <- filterByGene(dataset, genes)
      cordata <- prepareCorrelation()
      covar <- input$clinvar_cor
      
      if(dataset == "Cambridge"){
        
        if(covar == "iCluster"){
          
          data <- select(cordata, geo_accession,Expression, iCluster, Symbol) %>% 
            filter(iCluster %in% c("clust1","clust2","clust3","clust4","clust5")) %>% 
            #        mutate(Symbol=gsub(genes[1], "Gene1",Symbol)) %>% 
            #       mutate(Symbol=gsub(genes[2], "Gene2",Symbol)) %>% 
            spread(Symbol,Expression)
          #      cor <- round(cor(data$Gene1,data$Gene2,method=getCorType()),3)
          df <- as.data.frame(data)
          p <- ggpairs(df, columns = 3:ncol(df),col="iCluster")
          
          #        ggplot(data, aes(x = Gene1, y=Gene2,col=iCluster)) + geom_point() + xlab(genes[1]) + ylab(genes[2]) +  scale_fill_manual(values=iclusPal) + ggtitle(paste("Correlation = ",cor))
          
        }
        
        else if(covar == "Gleason"){
          
          data <- select(cordata, geo_accession,Expression, Gleason, Symbol) %>% 
            spread(Symbol,Expression)
          df <- as.data.frame(data)
          p <- ggpairs(df, columns = 3:ncol(df),col="Gleason")
          
        }
        
        else if(covar == "Sample_Group"){
          
          data <- select(cordata, geo_accession,Expression, Sample_Group, Symbol) %>% 
            spread(Symbol,Expression)
          df <- as.data.frame(data)
          p <- ggpairs(df, columns = 3:ncol(df),col="Sample_Group")
          
        }
        
        else {
          
          data <- select(cordata, geo_accession,Expression, Gleason, Symbol) %>% 
            spread(Symbol,Expression)
          df <- as.data.frame(data)
          p <- ggpairs(df, columns = 3:ncol(df))
        }
        
        
        
      }
      
      else if (dataset == "Stockholm"){
        
        if(covar == "iCluster"){
          
          data <- select(cordata, geo_accession,Expression, iCluster, Symbol) %>% 
            filter(iCluster %in% c("clust1","clust2","clust3","clust4","clust5")) %>% 
            #        mutate(Symbol=gsub(genes[1], "Gene1",Symbol)) %>% 
            #       mutate(Symbol=gsub(genes[2], "Gene2",Symbol)) %>% 
            spread(Symbol,Expression)
          #      cor <- round(cor(data$Gene1,data$Gene2,method=getCorType()),3)
          df <- as.data.frame(data)
          p <-  ggpairs(df, columns = 3:ncol(df),col="iCluster")
          
        }
        
        else if(covar == "Gleason"){
          
          data <- select(cordata, geo_accession,Expression, Gleason, Symbol) %>% 
            spread(Symbol,Expression)
          df <- as.data.frame(data)
          p <- ggpairs(df, columns = 3:ncol(df),col="Gleason")
          
        }
        
        else {
          
          data <- select(cordata, geo_accession,Expression, Gleason, Symbol) %>% 
            spread(Symbol,Expression)
          df <- as.data.frame(data)
          p <- ggpairs(df, columns = 3:ncol(df))
          
        }
        
        
      }
      
      else if (dataset=="MSKCC"){
        
        if(covar == "CopyNumberCluster"){
          
          data <- select(cordata, geo_accession,Expression, Copy.number.Cluster, Symbol) %>% 
            spread(Symbol,Expression)
          df <- as.data.frame(data)
          cor <- round(cor(data$Gene1,data$Gene2,method=getCorType()),3)
          p <- ggpairs(df, columns = 3:ncol(df),col="Copy.number.Cluster")
        }
        else if(covar == "Gleason"){
          
          
          data <- select(cordata, geo_accession,Expression, Gleason, Symbol) %>% 
            spread(Symbol,Expression)
          df <- as.data.frame(data)
          p <- ggpairs(df, columns = 3:ncol(df),col="Gleason")
        }
        
        else {
          
          data <- select(cordata, geo_accession,Expression, Gleason, Symbol) %>% 
            spread(Symbol,Expression)
          df <- as.data.frame(data)
          p <- ggpairs(df, columns = 3:ncol(df))
          
        }
        
        
      }
      
      else if (dataset == "Michigan2005"){
        
        if(covar == "Sample_Group"){
          data <- select(cordata, geo_accession,Expression, Sample_Group, Symbol) %>% 
            spread(Symbol,Expression)
          df <- as.data.frame(data)
          p <- ggpairs(df, columns = 3:ncol(df),col="Sample_Group")
          
        }
        else{
          data <- select(cordata, geo_accession,Expression, Sample_Group, Symbol) %>% 
            spread(Symbol,Expression)
          df <- as.data.frame(data)
          p <- ggpairs(df, columns = 3:ncol(df))
          
        }
        
      }
      
      else if (dataset == "Michigan2012"){
        
        if(covar == "Sample_Group"){
          data <- select(cordata, geo_accession,Expression, Sample_Group, Symbol) %>% 
            spread(Symbol,Expression)
          df <- as.data.frame(data)
          p <- ggpairs(df, columns = 3:ncol(df),col="Sample_Group")
          
        }
        else{
          data <- select(cordata, geo_accession,Expression, Sample_Group, Symbol) %>% 
            spread(Symbol,Expression)
          df <- as.data.frame(data)
          p <- ggpairs(df, columns = 3:ncol(df))
          
        }
        
      }
      
      print(p)   
      dev.off()   
    }
    
    
  )
  
  
  
  
  
  output$geneProfileScript <- downloadHandler(
    filename = function() {
      paste(input$profileBasename, '.R', sep='')
    },
    content = function(file) {
      cat(file=file,as.name("##Load the required libraries\n"))
      cat(file=file,as.name("library(ggplot2)\n"),append=TRUE)
      cat(file=file,as.name("library(tidyr)\n"),append=TRUE)
      cat(file=file,as.name("library(devtools)\n"),append=TRUE)
      cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
      
      cat(file=file,as.name("library(RColorBrewer)\n"),append=TRUE)
      cat(file=file,as.name("iclusPal <- brewer.pal(5, 'Set1')\n"),append=TRUE)
      
      dataset <- input$boxplotDataset
      
      if(input$inputType == "Single Gene"){
        gene <- input$currentGene
        cat(file=file,as.name(paste0("genes <- '",gene,"'\n")),append=TRUE)
        
      }
      else {
        if(is.null(input$file1)) cat(file=file,as.name("genes <- c('STAT3', 'ESR1','AR')\n"),append=TRUE)
        else{
          inFile <- input$file1
          cat(file=file,as.name(paste0('myfile <- \"' , inFile$name, '\"\n')),append=TRUE)
          cat(file=file,as.name("genes <- read.delim(myfile)[,1]\n"),append=TRUE)
        }
      }
      
      
      if(dataset == "Cambridge"){
        
        cat(file=file,as.name("if(!require(prostateCancerCamcap)) {\n source('http://www.bioconductor.org/biocLite.R')\n biocLite('prostateCancerCamcap')\n}\n"),append=TRUE)
        
        cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(camcap,package = 'prostateCancerCamcap')\n"),append=TRUE)
        cat(file=file,as.name("pd_camcap <- tbl_df(pData(camcap))\n"),append=TRUE)
        cat(file=file,as.name("fd_camcap <- tbl_df(fData(camcap))\n"),append=TRUE)
        cat(file=file,as.name("exp_camcap <- tbl_df(data.frame(ID = as.character(featureNames(camcap)),exprs(camcap)))\n"),append=TRUE)
        cat(file=file,as.name("probes <- fd_camcap %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file,as.name("data<- exp_camcap  %>% filter(ID %in% probes) %>%\n"),append=TRUE)        
        cat(file=file,as.name("gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file,as.name("fd <- fd_camcap\n"),append=TRUE)
        cat(file=file,as.name("pd <-  mutate(pd_camcap, Gleason=factor(Gleason,levels=c('5=3+2','6=2+4','6=3+3', '7=3+4','7=4+3','8=3+5','8=4+4','9=4+5','9=5+4','10=5+5',NA)))\n"),append=TRUE)
        
      }
      
      else if (dataset == "Stockholm"){
        cat(file=file,as.name("if(!require(prostateCancerStockholm)) {\n source('http://www.bioconductor.org/biocLite.R') \n biocLite('prostateCancerStockholm')\n}\n"),append=TRUE)
        
        cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(stockholm,package = 'prostateCancerStockholm')\n"),append=TRUE)
        cat(file=file,as.name("pd_stockholm <- tbl_df(pData(stockholm))\n"),append=TRUE)
        cat(file=file,as.name("fd_stockholm <- tbl_df(fData(stockholm))\n"),append=TRUE)
        cat(file=file,as.name("exp_stockholm <- tbl_df(data.frame(ID = as.character(featureNames(stockholm)),exprs(stockholm)))\n"),append=TRUE)
        
        cat(file=file,as.name("probes <- fd_stockholm %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file, as.name("data<- exp_stockholm  %>% filter(ID %in% probes) %>% \n"),append=TRUE)
        cat(file=file,as.name("gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file,as.name("fd <- fd_stockholm\n"),append=TRUE)
        cat(file=file,as.name("pd <-  mutate(pd_stockholm, Gleason=factor(Gleason,levels=c('5=3+2','6=2+4','6=3+3', '7=3+4','7=4+3','8=3+5','8=4+4','9=4+5','9=5+4','10=5+5',NA))) %>% \n"),append=TRUE)
        
      }
      
      else if(dataset =="MSKCC"){
        
        cat(file=file,as.name("if(!require(prostateCancerTaylor)) {\n source('http://www.bioconductor.org/biocLite.R') \n biocLite('prostateCancerTaylor')\n}\n"),append=TRUE)
        
        cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(taylor,package = 'prostateCancerTaylor')\n"),append=TRUE)
        cat(file=file,as.name("pd_taylor <- tbl_df(pData(taylor))\n"),append=TRUE)
        cat(file=file,as.name("fd_taylor <- tbl_df(fData(taylor))\n"),append=TRUE)
        cat(file=file,as.name("exp_taylor <- tbl_df(data.frame(ID = as.character(featureNames(taylor)),log2(exprs(taylor)))\n"),append=TRUE)
        
        
        cat(file=file,as.name("probes <- fd_taylor %>% filter(Gene %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file, as.name("data <- exp_taylor  %>% filter(ID %in% probes) %>% \n"),append=TRUE)
        cat(file=file, as.name("gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file, as.name("fd <- fd_taylor %>% mutate(Symbol = Gene)\n"),append=TRUE)
        cat(file=file, as.name("pd <- mutate(pd_taylor,Gleason = gsub('4+3', '7=4+3', pd_taylor$Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason = gsub('4+3', '8=5+3', Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason = gsub('3+4', '7=3+4', Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason = gsub('4+3', '7=4+3', Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason = gsub('3+3', '6=3+3', Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason = gsub('4+5', '9=4+5', Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason = gsub('3+5', '8=3+5', Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason = gsub('4+4', '8=4+4', Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason=factor(Gleason,levels=c('5=3+2','6=2+4','6=3+3', '7=3+4','7=4+3','8=3+5','8=4+4','9=4+5','9=5+4','10=5+5',NA)))\n"),append=TRUE)
        
      }
      
      else if(dataset == "Michigan2005"){
        
        cat(file=file,as.name("if(!require(prostateCancerVarambally)) {\nsource('http://www.bioconductor.org/biocLite.R')\n biocLite('prostateCancerVarambally')\n}\n"),append=TRUE)
        
        cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(varambally,package = 'prostateCancerVarambally')\n"),append=TRUE)
        cat(file=file,as.name("pd_varambally <- tbl_df(pData(varambally))\n"),append=TRUE)
        cat(file=file,as.name("fd_varambally <- tbl_df(fData(varambally))\n"),append=TRUE)
        cat(file=file,as.name("exp_varambally <- tbl_df(data.frame(ID = as.character(featureNames(varambally)),exprs(varambally)))\n"),append=TRUE)
        
        
        cat(file=file, as.name("probes <- fd_varambally %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file, as.name("data <- exp_varambally  %>% filter(ID %in% probes) %>% \n"),append=TRUE)
        cat(file=file, as.name("gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file, as.name("fd <- fd_varambally\n"),append=TRUE)
        cat(file=file, as.name("pd <- pd_varambally \n"),append=TRUE)
        
      }
      
      else if(dataset == "Michigan2012"){
        
        cat(file=file,as.name("if(!require(prostateCancerGrasso)) {\n source('http://www.bioconductor.org/biocLite.R')\nbiocLite('prostateCancerGrasso')\n}\n"),append=TRUE)
        
        
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(grasso,package = 'prostateCancerGrasso')\n"),append=TRUE)
        cat(file=file,as.name("pd_grasso <- tbl_df(pData(grasso))\n"),append=TRUE)
        cat(file=file,as.name("fd_grasso <- tbl_df(fData(grasso))\n"),append=TRUE)
        cat(file=file,as.name("exp_grasso <- tbl_df(data.frame(ID = as.character(featureNames(grasso)),exprs(grasso)))\n"),append=TRUE)
        
        
        cat(file=file, as.name("probes <- fd_grasso %>% filter(GENE_SYMBOL %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file, as.name("data <- exp_grasso  %>% filter(ID %in% probes) %>% \n"),append=TRUE)
        cat(file=file, as.name("gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file, as.name("fd <- mutate(fd_grasso, Symbol = GENE_SYMBOL)\n"),append=TRUE)
        cat(file=file, as.name("pd <- mutate(pd_grasso, Sample_Group = Group) \n"),append=TRUE)
        
        
      }
      
      cat(file=file, as.name("summary_stats <- data %>% group_by(ID) %>% \n"),append=TRUE)
      cat(file=file, as.name("summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))\n"),append=TRUE)
      cat(file=file, as.name("data <- left_join(data,summary_stats) %>% mutate(Z = (Expression - mean) / sd)\n"),append=TRUE)
      cat(file=file, as.name("mostVarProbes <- left_join(summary_stats,fd) %>% \n"),append=TRUE)
      cat(file=file, as.name("arrange(Symbol,desc(iqr)) %>% \n"),append=TRUE)
      cat(file=file ,as.name("distinct(Symbol) %>% \n"),append=TRUE)
      cat(file=file, as.name("select(ID) %>%  as.matrix %>%  as.character\n"),append=TRUE)        
      cat(file=file, as.name("data <- filter(data, ID %in% mostVarProbes)\n"),append=TRUE)
      cat(file=file, as.name("data <- left_join(data, select(fd, ID, Symbol))\n"),append=TRUE)
      cat(file=file, as.name("data <- left_join(data, pd)\n"),append=TRUE)
      
      
      doZ <- ifelse(getCambridgeZ() == "Yes",TRUE,FALSE)
      if(doZ){
        cat(file=file, as.name("data <- mutate(data, Expression=Z)\n"),append=TRUE)
        
      }
      var <- getCambridgeVariable()
      overlay <- getCambridgeOverlay()
      
      if(dataset == "Cambridge"){
        
        if(var=="iCluster") cat(file=file,as.name("p <- data %>%  filter(Sample_Group == 'Tumour',!is.na(iCluster)) %>% ggplot(aes(x = iCluster, y = Expression, fill=iCluster)) + geom_boxplot() +  scale_fill_manual(values=iclusPal)\n"),append=TRUE)
        else if(var == "Gleason") cat(file=file, as.name("p <- data  %>% filter(Sample_Group == 'Tumour') %>% filter(!(is.na(Gleason))) %>% ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot()\n"),append=TRUE)
        else if(var == "Sample_Group") cat(file=file, as.name("p <- data %>% mutate(Sample_Group = factor(Sample_Group,levels=c('Benign','Tumour','CRPC'))) %>% ggplot(aes(x = Sample_Group, y = Expression, fill=Sample_Group)) + geom_boxplot()\n"),append=TRUE)
        
      }
      
      else if(dataset == "Stockholm"){
        
        if(var =="iCluster") cat(file=file, as.name("p<- data %>% ggplot(aes(x = iCluster, y = Expression, fill=iCluster)) + geom_boxplot() +  scale_fill_manual(values=iclusPal)\n"),append=TRUE)
        else if(var == "Gleason") cat(file=file, as.name("p <- data  %>% filter(!(is.na(Gleason))) %>%  ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot()\n"),append=TRUE)
        
      }
      
      if(overlay=="Yes"){
        cat(file=file, as.name("p <- p + geom_jitter(position=position_jitter(width = .05),alpha=0.75)\n"),append=TRUE)
      }
      cat(file=file, as.name("p <- p + facet_wrap(~Symbol)\n"),append=TRUE)
      cat(file=file, as.name("print(p)\n"),append=TRUE)
    }
  )
  
  
  output$survivalScript <- downloadHandler(
    filename = function() {
      paste(input$survivalBasename, '.R', sep='')
    },
    content = function(file) {
      
      cat(file=file,as.name("##Load the required libraries\n"))
      cat(file=file,as.name("library(ggplot2)\n"),append=TRUE)
      cat(file=file,as.name("library(tidyr)\n"),append=TRUE)
      cat(file=file,as.name("library(devtools)\n"),append=TRUE)
      cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
      cat(file=file,as.name("library(party)\n"),append=TRUE)
      cat(file=file,as.name("library(survival)\n"),append=TRUE)
      
      dataset <- input$rpDataset
      
      if(input$inputType_survival == "Single Gene"){
        gene <- input$currentGene
      }
      else {
        gene <- input$survivalGeneChoice
      }
      
      cat(file=file,as.name(paste0("genes <- '",gene,"'\n")),append=TRUE)
      
      if(dataset == "Cambridge"){
        
        cat(file=file,as.name("if(!require(prostateCancerCamcap)) {\n source('http://www.bioconductor.org/biocLite.R')\n biocLite('prostateCancerCamcap')\n}\n"),append=TRUE)
        
        cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(camcap,package = 'prostateCancerCamcap')\n"),append=TRUE)
        cat(file=file,as.name("pd_camcap <- tbl_df(pData(camcap))\n"),append=TRUE)
        cat(file=file,as.name("fd_camcap <- tbl_df(fData(camcap))\n"),append=TRUE)
        cat(file=file,as.name("exp_camcap <- tbl_df(data.frame(ID = as.character(featureNames(camcap)),exprs(camcap)))\n"),append=TRUE)
        cat(file=file,as.name("probes <- fd_camcap %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file,as.name("data<- exp_camcap  %>% filter(ID %in% probes) %>%\n"),append=TRUE)        
        cat(file=file,as.name("gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file,as.name("fd <- fd_camcap\n"),append=TRUE)
        cat(file=file, as.name("pd <- mutate(pd_camcap, Time = as.numeric(FollowUpTime), Event = ifelse(BCR=='Y',1,0))\n"),append=TRUE)  
        
      }
      
      else if (dataset == "Stockholm"){
        cat(file=file,as.name("if(!require(prostateCancerStockholm)) {\n source('http://www.bioconductor.org/biocLite.R') \n biocLite('prostateCancerStockholm')\n}\n"),append=TRUE)
        
        cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(stockholm,package = 'prostateCancerStockholm')\n"),append=TRUE)
        cat(file=file,as.name("pd_stockholm <- tbl_df(pData(stockholm))\n"),append=TRUE)
        cat(file=file,as.name("fd_stockholm <- tbl_df(fData(stockholm))\n"),append=TRUE)
        cat(file=file,as.name("exp_stockholm <- tbl_df(data.frame(ID = as.character(featureNames(stockholm)),exprs(stockholm)))\n"),append=TRUE)
        
        cat(file=file,as.name("probes <- fd_stockholm %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file, as.name("data<- exp_stockholm  %>% filter(ID %in% probes) %>% \n"),append=TRUE)
        cat(file=file,as.name("gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file,as.name("fd <- fd_stockholm\n"),append=TRUE)
        cat(file=file, as.name("pd <- mutate(pd_stockholm, Time = as.numeric(FollowUpTime), Event = ifelse(BCR=='Y',1,0))\n"),append=TRUE)  
        
        
      }
      
      else if(dataset =="MSKCC"){
        
        cat(file=file,as.name("if(!require(prostateCancerTaylor)) {\n source('http://www.bioconductor.org/biocLite.R') \n biocLite('prostateCancerTaylor')\n}\n"),append=TRUE)
        
        cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(taylor,package = 'prostateCancerTaylor')\n"),append=TRUE)
        cat(file=file,as.name("pd_taylor <- tbl_df(pData(taylor))\n"),append=TRUE)
        cat(file=file,as.name("fd_taylor <- tbl_df(fData(taylor))\n"),append=TRUE)
        cat(file=file,as.name("exp_taylor <- tbl_df(data.frame(ID = as.character(featureNames(taylor)),log2(exprs(taylor)))\n"),append=TRUE)
        
        
        cat(file=file,as.name("probes <- fd_taylor %>% filter(Gene %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file, as.name("data <- exp_taylor  %>% filter(ID %in% probes) %>% \n"),append=TRUE)
        cat(file=file, as.name("gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file, as.name("fd <- fd_taylor %>% mutate(Symbol = Gene)\n"),append=TRUE)
        cat(file=file, as.name("pd <- pd_taylor\n"),append=TRUE)
      }
      
      cat(file=file, as.name("summary_stats <- data %>% group_by(ID) %>% \n"),append=TRUE)
      cat(file=file, as.name("summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))\n"),append=TRUE)
      cat(file=file, as.name("data <- left_join(data,summary_stats) %>% mutate(Z = (Expression - mean) / sd)\n"),append=TRUE)
      cat(file=file, as.name("mostVarProbes <- left_join(summary_stats,fd) %>% \n"),append=TRUE)
      cat(file=file, as.name("arrange(Symbol,desc(iqr)) %>% \n"),append=TRUE)
      cat(file=file ,as.name("distinct(Symbol) %>% \n"),append=TRUE)
      cat(file=file, as.name("select(ID) %>%  as.matrix %>%  as.character\n"),append=TRUE)        
      cat(file=file, as.name("data <- filter(data, ID %in% mostVarProbes)\n"),append=TRUE)
      cat(file=file, as.name("data <- left_join(data, select(fd, ID, Symbol))\n"),append=TRUE)
      cat(file=file, as.name("data <- left_join(data, pd)\n"),append=TRUE)
      cat(file=file, as.name("data <- data %>% filter(!is.na(Time) & !is.na(Event))\n"),append=TRUE)
      cat(file=file, as.name("surv.xfs <- Surv((as.numeric(as.character(data$Time))/12), data$Event)\n"),append=TRUE)
      cat(file=file, as.name("data$surv.xfs <- surv.xfs\n"),append=TRUE)
      cat(file=file, as.name("ctree_xfs   <- ctree(surv.xfs~Expression,data=data)\n"),append=TRUE)
      cat(file=file, as.name("pvalue <- 1 - ctree_xfs@tree$criterion$maxcriterion\n"),append=TRUE)
      cat(file=file, as.name("newPval <- signif(pvalue, digits = 2)\n"),append=TRUE)
      
      cat(file=file, as.name("if(newPval<0.05) {\n"),append=TRUE)
      cat(file=file, as.name("ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)\n"),append=TRUE)
      cat(file=file, as.name("ps  <- signif(ps2[1], digits = 3)\n"),append=TRUE)
      cat(file=file, as.name("if(length(ps2)==1) {\n"),append=TRUE)
      cat(file=file, as.name("data$geneexp_cp <- data$Expression<=ps2[1]\n"),append=TRUE)
      cat(file=file, as.name("nt <- table(data$geneexp_cp)\n"),append=TRUE)
      cat(file=file, as.name("geneexp.survfit.xfs <- survfit(surv.xfs~geneexp_cp,data=data)\n"),append=TRUE)
      cat(file=file, as.name("plot(geneexp.survfit.xfs, xlab='Time to BCR (years)', ylab='Probability of Freedom from Biochemical Recurrence', main=paste(genes,', p=', newPval), col=c(2,4))\n"),append=TRUE)
      cat(file=file, as.name("legend('bottomleft', c(paste(genes, '>', ps, 'n=', nt[[1]]), paste(genes, '<=', ps, 'n=', nt[[2]])), col=c(2,4), lty=1, lwd=1.5, bty='n')\n"),append=TRUE)
      cat(file=file, as.name("}\n } else {\n"),append=TRUE)
      cat(file=file, as.name("ps <- round(median(data$Expression),3)\n"),append=TRUE)
      cat(file=file, as.name(" data$geneexp_cp <- data$Expression<= ps\n"),append=TRUE)
      cat(file=file, as.name("nt <- table(data$geneexp_cp)\n"),append=TRUE)
      cat(file=file, as.name("geneexp.survfit.xfs <- survfit(surv.xfs~geneexp_cp,data=data)\n"),append=TRUE)
      cat(file=file, as.name("test <- survdiff(surv.xfs~geneexp_cp,data=data)\n"),append=TRUE)
      cat(file=file, as.name("newPval <- round(pchisq(test$chisq, df = length(test$n)-1, lower.tail=FALSE),3)\n"),append=TRUE)
      cat(file=file, as.name("plot(geneexp.survfit.xfs, xlab='Time to BCR (years)', ylab='Probability of Freedom from Biochemical Recurrence', main=paste(genes,', p=', newPval), col=c(2,4))\n"),append=TRUE)
      cat(file=file, as.name("legend('bottomleft', c(paste(genes, '>', ps, 'n=', nt[[1]]), paste(genes, '<=', ps, 'n=', nt[[2]])), col=c(2,4), lty=1, lwd=1.5, bty='n')\n"),append=TRUE)
      cat(file=file, as.name("}\n"),append=TRUE)
      
    }
  )
  
  
  output$correlationScript <- downloadHandler(
    filename = function() {
      paste(input$correlationBasename, '.R', sep='')
    },
    content = function(file) {
      cat(file=file,as.name("##Load the required libraries\n"))
      cat(file=file,as.name("library(ggplot2)\n"),append=TRUE)
      cat(file=file,as.name("library(tidyr)\n"),append=TRUE)
      cat(file=file,as.name("library(devtools)\n"),append=TRUE)
      cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
      cat(file=file,as.name("library(GGally)\n"),append=TRUE)
      cat(file=file,as.name("library(RColorBrewer)\n"),append=TRUE)
      cat(file=file,as.name("iclusPal <- brewer.pal(5, 'Set1')\n"),append=TRUE)
      
      dataset <- input$corDataset
      
      if(input$inputType_correlation == "Single Gene"){
        gene1 <- getCurrentGene()
        gene2 <- input$secondGene
        genes <- c(gene1,gene2)
        
        cat(file=file,as.name(paste0("genes <- c('",gene1,"','",gene2,"')\n")),append=TRUE)
        
      }
      else {
        if(is.null(input$file1)) cat(file=file,as.name("genes <- c('STAT3', 'ESR1','AR')\n"),append=TRUE)
        else{
          inFile <- input$file1
          cat(file=file,as.name(paste0('myfile <- \"' , inFile$name, '\"\n')),append=TRUE)
          cat(file=file,as.name("genes <- read.delim(myfile)[,1]\n"),append=TRUE)
        }
      }
      
      
      if(dataset == "Cambridge"){
        
        cat(file=file,as.name("if(!require(prostateCancerCamcap)) {\n source('http://www.bioconductor.org/biocLite.R')\n biocLite('prostateCancerCamcap')\n}\n"),append=TRUE)
        
        cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(camcap,package = 'prostateCancerCamcap')\n"),append=TRUE)
        cat(file=file,as.name("pd_camcap <- tbl_df(pData(camcap))\n"),append=TRUE)
        cat(file=file,as.name("fd_camcap <- tbl_df(fData(camcap))\n"),append=TRUE)
        cat(file=file,as.name("exp_camcap <- tbl_df(data.frame(ID = as.character(featureNames(camcap)),exprs(camcap)))\n"),append=TRUE)
        cat(file=file,as.name("probes <- fd_camcap %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file,as.name("data<- exp_camcap  %>% filter(ID %in% probes) %>%\n"),append=TRUE)        
        cat(file=file,as.name("gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file,as.name("fd <- fd_camcap\n"),append=TRUE)
        cat(file=file,as.name("pd <-  mutate(pd_camcap, Gleason=factor(Gleason,levels=c('5=3+2','6=2+4','6=3+3', '7=3+4','7=4+3','8=3+5','8=4+4','9=4+5','9=5+4','10=5+5',NA)))\n"),append=TRUE)
        
      }
      
      else if (dataset == "Stockholm"){
        cat(file=file,as.name("if(!require(prostateCancerStockholm)) {\n source('http://www.bioconductor.org/biocLite.R') \n biocLite('prostateCancerStockholm')\n}\n"),append=TRUE)
        
        cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(stockholm,package = 'prostateCancerStockholm')\n"),append=TRUE)
        cat(file=file,as.name("pd_stockholm <- tbl_df(pData(stockholm))\n"),append=TRUE)
        cat(file=file,as.name("fd_stockholm <- tbl_df(fData(stockholm))\n"),append=TRUE)
        cat(file=file,as.name("exp_stockholm <- tbl_df(data.frame(ID = as.character(featureNames(stockholm)),exprs(stockholm)))\n"),append=TRUE)
        
        cat(file=file,as.name("probes <- fd_stockholm %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file, as.name("data<- exp_stockholm  %>% filter(ID %in% probes) %>% \n"),append=TRUE)
        cat(file=file,as.name("gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file,as.name("fd <- fd_stockholm\n"),append=TRUE)
        cat(file=file,as.name("pd <-  mutate(pd_stockholm, Gleason=factor(Gleason,levels=c('5=3+2','6=2+4','6=3+3', '7=3+4','7=4+3','8=3+5','8=4+4','9=4+5','9=5+4','10=5+5',NA))) %>% \n"),append=TRUE)
        
      }
      
      else if(dataset =="MSKCC"){
        
        cat(file=file,as.name("if(!require(prostateCancerTaylor)) {\n source('http://www.bioconductor.org/biocLite.R') \n biocLite('prostateCancerTaylor')\n}\n"),append=TRUE)
        
        cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(taylor,package = 'prostateCancerTaylor')\n"),append=TRUE)
        cat(file=file,as.name("pd_taylor <- tbl_df(pData(taylor))\n"),append=TRUE)
        cat(file=file,as.name("fd_taylor <- tbl_df(fData(taylor))\n"),append=TRUE)
        cat(file=file,as.name("exp_taylor <- tbl_df(data.frame(ID = as.character(featureNames(taylor)),log2(exprs(taylor)))\n"),append=TRUE)
        
        
        cat(file=file,as.name("probes <- fd_taylor %>% filter(Gene %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file, as.name("data <- exp_taylor  %>% filter(ID %in% probes) %>% \n"),append=TRUE)
        cat(file=file, as.name("gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file, as.name("fd <- fd_taylor %>% mutate(Symbol = Gene)\n"),append=TRUE)
        cat(file=file, as.name("pd <- mutate(pd_taylor,Gleason = gsub('4+3', '7=4+3', pd_taylor$Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason = gsub('4+3', '8=5+3', Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason = gsub('3+4', '7=3+4', Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason = gsub('4+3', '7=4+3', Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason = gsub('3+3', '6=3+3', Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason = gsub('4+5', '9=4+5', Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason = gsub('3+5', '8=3+5', Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason = gsub('4+4', '8=4+4', Gleason,fixed=TRUE)) %>% \n"),append=TRUE)
        cat(file=file, as.name("mutate(Gleason=factor(Gleason,levels=c('5=3+2','6=2+4','6=3+3', '7=3+4','7=4+3','8=3+5','8=4+4','9=4+5','9=5+4','10=5+5',NA)))\n"),append=TRUE)
        
      }
      
      else if(dataset == "Michigan2005"){
        
        cat(file=file,as.name("if(!require(prostateCancerVarambally)) {\nsource('http://www.bioconductor.org/biocLite.R')\n biocLite('prostateCancerVarambally')\n}\n"),append=TRUE)
        
        cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(varambally,package = 'prostateCancerVarambally')\n"),append=TRUE)
        cat(file=file,as.name("pd_varambally <- tbl_df(pData(varambally))\n"),append=TRUE)
        cat(file=file,as.name("fd_varambally <- tbl_df(fData(varambally))\n"),append=TRUE)
        cat(file=file,as.name("exp_varambally <- tbl_df(data.frame(ID = as.character(featureNames(varambally)),exprs(varambally)))\n"),append=TRUE)
        
        
        cat(file=file, as.name("probes <- fd_varambally %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file, as.name("data <- exp_varambally  %>% filter(ID %in% probes) %>% \n"),append=TRUE)
        cat(file=file, as.name("gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file, as.name("fd <- fd_varambally\n"),append=TRUE)
        cat(file=file, as.name("pd <- pd_varambally \n"),append=TRUE)
        
      }
      
      else if(dataset == "Michigan2012"){
        
        cat(file=file,as.name("if(!require(prostateCancerGrasso)) {\n source('http://www.bioconductor.org/biocLite.R')\nbiocLite('prostateCancerGrasso')\n}\n"),append=TRUE)
        
        
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(grasso,package = 'prostateCancerGrasso')\n"),append=TRUE)
        cat(file=file,as.name("pd_grasso <- tbl_df(pData(grasso))\n"),append=TRUE)
        cat(file=file,as.name("fd_grasso <- tbl_df(fData(grasso))\n"),append=TRUE)
        cat(file=file,as.name("exp_grasso <- tbl_df(data.frame(ID = as.character(featureNames(grasso)),exprs(grasso)))\n"),append=TRUE)
        
        
        cat(file=file, as.name("probes <- fd_grasso %>% filter(GENE_SYMBOL %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file, as.name("data <- exp_grasso  %>% filter(ID %in% probes) %>% \n"),append=TRUE)
        cat(file=file, as.name("gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file, as.name("fd <- mutate(fd_grasso, Symbol = GENE_SYMBOL)\n"),append=TRUE)
        cat(file=file, as.name("pd <- mutate(pd_grasso, Sample_Group = Group) \n"),append=TRUE)
        
        
      }
      
      cat(file=file, as.name("summary_stats <- data %>% group_by(ID) %>% \n"),append=TRUE)
      cat(file=file, as.name("summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))\n"),append=TRUE)
      cat(file=file, as.name("data <- left_join(data,summary_stats) %>% mutate(Z = (Expression - mean) / sd)\n"),append=TRUE)
      cat(file=file, as.name("mostVarProbes <- left_join(summary_stats,fd) %>% \n"),append=TRUE)
      cat(file=file, as.name("arrange(Symbol,desc(iqr)) %>% \n"),append=TRUE)
      cat(file=file ,as.name("distinct(Symbol) %>% \n"),append=TRUE)
      cat(file=file, as.name("select(ID) %>%  as.matrix %>%  as.character\n"),append=TRUE)        
      cat(file=file, as.name("data <- filter(data, ID %in% mostVarProbes)\n"),append=TRUE)
      cat(file=file, as.name("data <- left_join(data, select(fd, ID, Symbol))\n"),append=TRUE)
      cat(file=file, as.name("data <- left_join(data, pd)\n"),append=TRUE)
      
      covar <- input$clinvar_cor
      
      if(covar == "iCluster"){
        
        cat(file=file, as.name("data <- select(data, geo_accession,Expression, iCluster, Symbol) %>% \n"),append=TRUE)
        cat(file=file, as.name("filter(iCluster %in% c('clust1','clust2','clust3','clust4','clust5')) %>% \n"),append=TRUE)
        cat(file=file, as.name("spread(Symbol,Expression)\n"),append=TRUE)
        cat(file=file, as.name("df <- as.data.frame(data)\n"),append=TRUE)
        cat(file=file, as.name("p <- ggpairs(df, columns = 3:ncol(df),col='iCluster')\n"),append=TRUE)
      }
      
      else if(covar == "Gleason"){
        cat(file=file, as.name("data <- select(data, geo_accession,Expression, Gleason, Symbol) %>% \n"),append=TRUE)
        cat(file=file, as.name("spread(Symbol,Expression)\n"),append=TRUE)
        cat(file=file, as.name("df <- as.data.frame(data)"),append=TRUE)
        cat(file=file, as.name("p <- ggpairs(df, columns = 3:ncol(df),col='Gleason')\n"),append=TRUE)
        
      }
      
      else if(covar == "Sample_Group"){
        cat(file=file, as.name("data <- select(data, geo_accession,Expression, Sample_Group, Symbol) %>% \n"),append=TRUE)
        cat(file=file, as.name("spread(Symbol,Expression)\n"),append=TRUE)
        cat(file=file, as.name("df <- as.data.frame(data)\n"),append=TRUE)
        cat(file=file, as.name("p <- ggpairs(df, columns = 3:ncol(df),col='Sample_Group')\n"),append=TRUE)
        
      }
      else {
        cat(file=file, as.name("data <- select(data, geo_accession,Expression, Gleason, Symbol) %>%\n"),append=TRUE)
        cat(file=file, as.name(" spread(Symbol,Expression)\n"),append=TRUE)
        cat(file=file, as.name("df <- as.data.frame(data)\n"),append=TRUE)
        cat(file=file, as.name("p <- ggpairs(df, columns = 3:ncol(df))\n"),append=TRUE)
        
      }
      cat(file=file,as.name("print(p)\n"),append=TRUE)       
    }
    
    
    
  )
  
  
  
  output$heatmapScript <- downloadHandler(
    filename = function() {
      paste(input$heatmapBasename, '.R', sep='')
    },
    content = function(file) {
      cat(file=file,as.name("##Load the required libraries\n"))
      cat(file=file,as.name("library(ggplot2)\n"),append=TRUE)
      cat(file=file,as.name("library(tidyr)\n"),append=TRUE)
      cat(file=file,as.name("library(devtools)\n"),append=TRUE)
      cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
      cat(file=file,as.name("library(RColorBrewer)\n"),append=TRUE)
      cat(file=file,as.name("iclusPal <- brewer.pal(5, 'Set1')\n"),append=TRUE)
      
      
      if(is.null(input$file1)) cat(file=file,as.name("genes <- c('STAT3', 'ESR1','AR')\n"),append=TRUE)
      else{
        inFile <- input$file1
        cat(file=file,as.name(paste0('myfile <- \"' , inFile$name, '\"\n')),append=TRUE)
        cat(file=file,as.name("genes <- read.delim(myfile)[,1]\n"),append=TRUE)
      }
      
      dataset <- input$heatmapDataset
      
      if(dataset == "MSKCC"){
        cat(file=file,as.name("if(!require(prostateCancerTaylor)) {\n source('http://www.bioconductor.org/biocLite.R') \n biocLite('prostateCancerTaylor')\n}\n"),append=TRUE)
        
        cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(taylor,package = 'prostateCancerTaylor')\n"),append=TRUE)
        cat(file=file,as.name("pd_taylor <- tbl_df(pData(taylor))\n"),append=TRUE)
        cat(file=file,as.name("fd_taylor <- tbl_df(fData(taylor))\n"),append=TRUE)
        cat(file=file,as.name("exp_taylor <- tbl_df(data.frame(ID = as.character(featureNames(taylor)),log2(exprs(taylor))))n"),append=TRUE)
        
        
        cat(file=file, as.name("probes <- fd_taylor %>% filter(Gene %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file, as.name("taylor <- exp_taylor  %>% filter(ID %in% probes) %>% \n"),append=TRUE)
        cat(file=file, as.name("gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file, as.name("summary_stats <- taylor %>% group_by(ID) %>% \n"),append=TRUE)
        cat(file=file, as.name("summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))\n"),append=TRUE)
        cat(file=file, as.name("mostVarProbes <- left_join(summary_stats,fd_taylor) %>% \n"),append=TRUE)
        cat(file=file, as.name("arrange(Gene,desc(iqr)) %>% \n"),append=TRUE)
        cat(file=file, as.name("distinct(Gene) %>% \n"),append=TRUE)
        cat(file=file, as.name("select(ID) %>%  as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file, as.name("samples <- filter(pd_taylor, Sample_Group == 'prostate cancer') %>% \n"),append=TRUE)
        cat(file=file, as.name("select(geo_accession) %>% as.matrix %>% as.character\n"),append=TRUE)
        cat(file=file, as.name("taylor <- filter(taylor, ID %in% mostVarProbes,geo_accession %in% samples)\n"),append=TRUE)
        cat(file=file, as.name("geneMatrix <- taylor %>% \n"),append=TRUE)
        cat(file=file, as.name("spread(geo_accession,Expression) %>% data.frame\n"),append=TRUE)
        cat(file=file, as.name("geneMatrix <- as.matrix(geneMatrix[,-1])\n"),append=TRUE)
        cat(file=file, as.name("symbols <- filter(fd_taylor, ID %in% mostVarProbes) %>% select(Gene) %>% as.matrix %>% as.character\n"),append=TRUE)
        cat(file=file, as.name("rownames(geneMatrix) <- symbols\n"),append=TRUE)
        cat(file=file, as.name("pd <- left_join(taylor,pd_taylor) %>% distinct(geo_accession)\n"),append=TRUE)
        cat(file=file, as.name("colMatrix <- matrix(nrow = ncol(geneMatrix),ncol = 2)\n"),append=TRUE)
        cat(file=file, as.name("grp <- pd$Copy.number.Cluster\n"),append=TRUE)
        cat(file=file, as.name("cols <- brewer.pal(9,'Set1')[1:length(levels(factor(as.character(grp))))]\n"),append=TRUE)
        cat(file=file, as.name("grp <- as.factor(as.character(grp))\n"),append=TRUE)
        cat(file=file, as.name("levels(grp) <- cols\n"),append=TRUE)
        cat(file=file, as.name("colMatrix[,1] <- as.character(grp)\n"),append=TRUE)
        cat(file=file, as.name("grp <- pd$Gleason\n"),append=TRUE)
        cat(file=file, as.name("cols <- brewer.pal(8,'Set2')[1:length(levels(factor(as.character(grp))))]\n"),append=TRUE)
        cat(file=file, as.name("grp <- as.factor(as.character(grp))\n"),append=TRUE)
        cat(file=file, as.name("levels(grp) <- cols\n"),append=TRUE)
        cat(file=file, as.name("colMatrix[,2] <- as.character(grp)\n"),append=TRUE)
        cat(file=file, as.name("colnames(colMatrix) <- c('Copy Number Cluster', 'Gleason')\n"),append=TRUE)
      }  
      else if(dataset == "Cambridge") {
        cat(file=file,as.name("if(!require(prostateCancerCamcap)) {\n source('http://www.bioconductor.org/biocLite.R')\n biocLite('prostateCancerCamcap')\n}\n"),append=TRUE)
        
        cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(camcap,package = 'prostateCancerCamcap')\n"),append=TRUE)
        cat(file=file,as.name("pd_camcap <- tbl_df(pData(camcap))\n"),append=TRUE)
        cat(file=file,as.name("fd_camcap <- tbl_df(fData(camcap))\n"),append=TRUE)
        cat(file=file,as.name("exp_camcap <- tbl_df(data.frame(ID = as.character(featureNames(camcap)),exprs(camcap)))\n"),append=TRUE)
        
        cat(file=file, as.name("probes <- fd_camcap %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file, as.name("camcap <- exp_camcap  %>% filter(ID %in% probes) %>% gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file, as.name("summary_stats <- camcap %>% group_by(ID) %>% summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))\n"),append=TRUE)
        cat(file=file, as.name("mostVarProbes <- left_join(summary_stats,fd_camcap) %>% arrange(Symbol,desc(iqr)) %>% distinct(Symbol) %>%   select(ID) %>%  as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file, as.name("samples <- filter(pd_camcap, Sample_Group == 'Tumour') %>% select(geo_accession) %>% as.matrix %>% as.character\n"),append=TRUE)
        cat(file=file, as.name("camcap <- filter(camcap, ID %in% mostVarProbes,geo_accession %in% samples) \n"),append=TRUE)
        cat(file=file, as.name("geneMatrix <- camcap %>% spread(geo_accession,Expression) %>% data.frame\n"),append=TRUE)
        cat(file=file, as.name("geneMatrix <- as.matrix(geneMatrix[,-1])\n"),append=TRUE)
        cat(file=file, as.name("symbols <- filter(fd_camcap, ID %in% mostVarProbes) %>% select(Symbol) %>% as.matrix %>% as.character\n"),append=TRUE)
        cat(file=file, as.name("rownames(geneMatrix) <- symbols\n"),append=TRUE)
        cat(file=file, as.name("pd <- left_join(camcap,pd_camcap) %>% distinct(geo_accession)\n"),append=TRUE)
        cat(file=file, as.name("colMatrix <- matrix(nrow = ncol(geneMatrix),ncol = 2)\n"),append=TRUE)
        cat(file=file, as.name("grp <- pd$iCluster\n"),append=TRUE)
        cat(file=file, as.name("cols <- brewer.pal(5,'Set1')[1:length(levels(factor(as.character(grp))))]\n"),append=TRUE)
        cat(file=file, as.name(" grp <- as.factor(as.character(grp))\n"),append=TRUE)
        cat(file=file, as.name("levels(grp) <- cols\n"),append=TRUE)
        cat(file=file, as.name("colMatrix[,1] <- as.character(grp)\n"),append=TRUE)
        cat(file=file, as.name("grp <- pd$Gleason\n"),append=TRUE)
        cat(file=file, as.name("cols <- brewer.pal(8,'Set2')[1:length(levels(factor(as.character(grp))))]\n"),append=TRUE)
        cat(file=file, as.name("grp <- as.factor(as.character(grp))\n"),append=TRUE)
        cat(file=file, as.name("levels(grp) <- cols\n"),append=TRUE)
        cat(file=file, as.name("colMatrix[,2] <- as.character(grp)\n"),append=TRUE)
        cat(file=file, as.name("colnames(colMatrix) <- c('iCluster', 'Gleason')\n"),append=TRUE)
      }
      
      else if (dataset=="Stockholm"){
        cat(file=file,as.name("if(!require(prostateCancerStockholm)) {\n source('http://www.bioconductor.org/biocLite.R')\n biocLite('prostateCancerStockholm')\n}\n"),append=TRUE)
        
        cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
        cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
        cat(file=file,as.name("data(stockholm,package = 'prostateCancerStockholm')\n"),append=TRUE)
        cat(file=file,as.name("pd_stockholm <- tbl_df(pData(stockholm))\n"),append=TRUE)
        cat(file=file,as.name("fd_stockholm <- tbl_df(fData(stockholm))\n"),append=TRUE)
        cat(file=file,as.name("exp_stockholm <- tbl_df(data.frame(ID = as.character(featureNames(stockholm)),exprs(stockholm)))\n"),append=TRUE)
        
        cat(file=file, as.name("probes <- fd_stockholm %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n"),append=TRUE)
        cat(file=file, as.name("stockholm <- exp_stockholm  %>% filter(ID %in% probes) %>% gather(geo_accession,Expression,-ID)\n"),append=TRUE)
        cat(file=file, as.name("summary_stats <- stockholm %>% group_by(ID) %>% summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))\n"),append=TRUE)
        cat(file=file, as.name("mostVarProbes <- left_join(summary_stats,fd_camcap) %>% arrange(Symbol,desc(iqr)) %>% distinct(Symbol) %>%   select(ID) %>%  as.matrix %>%  as.character\n"),append=TRUE)
        #      cat(file=file, as.name("stockholm <- filter(stockholm, ID %in% mostVarProbes,geo_accession %in% samples) \n"),append=TRUE)
        cat(file=file, as.name("geneMatrix <- stockholm %>% spread(geo_accession,Expression) %>% data.frame\n"),append=TRUE)
        cat(file=file, as.name("geneMatrix <- as.matrix(geneMatrix[,-1])\n"),append=TRUE)
        cat(file=file, as.name("symbols <- filter(fd_stockholm, ID %in% mostVarProbes) %>% select(Symbol) %>% as.matrix %>% as.character\n"),append=TRUE)
        cat(file=file, as.name("rownames(geneMatrix) <- symbols\n"),append=TRUE)
        cat(file=file, as.name("pd <- left_join(camcap,pd_camcap) %>% distinct(geo_accession)\n"),append=TRUE)
        cat(file=file, as.name("colMatrix <- matrix(nrow = ncol(geneMatrix),ncol = 2)\n"),append=TRUE)
        cat(file=file, as.name("grp <- pd$iCluster\n"),append=TRUE)
        cat(file=file, as.name("cols <- brewer.pal(5,'Set1')[1:length(levels(factor(as.character(grp))))]\n"),append=TRUE)
        cat(file=file, as.name(" grp <- as.factor(as.character(grp))\n"),append=TRUE)
        cat(file=file, as.name("levels(grp) <- cols\n"),append=TRUE)
        cat(file=file, as.name("colMatrix[,1] <- as.character(grp)\n"),append=TRUE)
        cat(file=file, as.name("grp <- pd$Gleason\n"),append=TRUE)
        cat(file=file, as.name("cols <- brewer.pal(8,'Set2')[1:length(levels(factor(as.character(grp))))]\n"),append=TRUE)
        cat(file=file, as.name("grp <- as.factor(as.character(grp))\n"),append=TRUE)
        cat(file=file, as.name("levels(grp) <- cols\n"),append=TRUE)
        cat(file=file, as.name("colMatrix[,2] <- as.character(grp)\n"),append=TRUE)
        cat(file=file, as.name("colnames(colMatrix) <- c('iCluster', 'Gleason')\n"),append=TRUE)
        
        
      }
      
      if(getDistFun() == "Correlation") cat(file=file, as.name("distfun <- function(x) as.dist(1 - cor(t(x)))\n"),append=TRUE)
      else cat(file=file, as.name("distfun <- dist\n"),append=TRUE)
      
      hclustfun <- getHclustMethod()
      
      cat(file=file, as.name(paste0("hclustfun <- function(x) hclust(x,method='",hclustfun,"')\n")),append=TRUE)
      
      cat(file=file,as.name("dMat <- as.dist(distfun(t(geneMatrix)))\n"),append=TRUE)
      
      cat(file=file, as.name("clusObj <- hclustfun(dMat)\n"),append=TRUE)
      
      cat(file=file, as.name("hmcol <- rev(brewer.pal(11 , 'RdBu'))\n"),append=TRUE)  
      
      if(getDistFun() == "Correlation") cat(file=file, as.name("distfun <- function(x) as.dist(1 - cor(t(x)))\n"),append=TRUE)
      else cat(file=file, as.name("distfun <- dist\n"),append=TRUE)
      
      scale <- getScaleMethod()
      
      cat(file=file, as.name(paste0("scale <- '", scale,"'\n")),append=TRUE)
      
      if(getReordRows() == "Yes"){
        
        cat(file=file, as.name("heatmap.2(geneMatrix,Colv = as.dendrogram(clusObj),col=hmcol,distfun=distfun,hclustfun = hclustfun,scale=scale,trace='none',cexRow = 0.9)\n"),append=TRUE)
        
      }
      else cat(file=file, as.name("heatmap.2(geneMatrix,Colv = as.dendrogram(clusObj),col=hmcol,distfun=distfun,Rowv=NA,hclustfun = hclustfun,scale=scale,trace='none',cexRow = 0.9)\n"),append=TRUE)
      
    }
    
  )
  
}
)