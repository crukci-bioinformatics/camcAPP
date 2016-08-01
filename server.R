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

db_camcap <- src_sqlite("camcap.sqlite3")

pd_camcap <- collect(tbl("pd",src=db_camcap))
fd_camcap <- collect(tbl("fd",src=db_camcap))


exp_camcap <- tbl("expression",src=db_camcap)

db_stockholm <- src_sqlite("stockholm.sqlite3")
pd_stockholm <- collect(tbl("pd",src=db_stockholm))
fd_stockholm <- collect(tbl("fd",src=db_stockholm))


exp_stockholm <- tbl("expression",src=db_stockholm)



db_taylor <- src_sqlite("taylor.sqlite3")
pd_taylor <- collect(tbl("pd",src=db_taylor))
fd_taylor <- collect(tbl("fd",src=db_taylor))
exp_taylor <- tbl("expression",src=db_taylor)


db_varambally <- src_sqlite("varambally.sqlite3")
pd_varambally <- collect(tbl("pd",src=db_varambally))
fd_varambally <- collect(tbl("fd",src=db_varambally))

exp_varambally <- tbl("expression",src=db_varambally)


db_grasso <- src_sqlite("grasso.sqlite3")
pd_grasso <- collect(tbl("pd",src=db_grasso))
fd_grasso <- collect(tbl("fd",src=db_grasso))

exp_grasso <- tbl("expression",src=db_grasso)


## dplyr needs to be the last package loaded, otherwise the 'select' function seems to be over-written

select <- dplyr::select
iclusPal <- brewer.pal(5, "Set1")
gradeCols <- rev(brewer.pal(11, "RdYlGn"))
message("READY FOR INPUT")




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
  
  getSecondGene <- reactive({
    
    gene2 <- input$secondGene
    
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
  
  
  getAllGenes <- reactive({
   ###Get the names of all the genes that could currently be used for analysis
    ##Gene name selected from drop-down on first page
    ##Genes from gene list currently uploaded
    ##Gene from correlation page
    
    gene1 <- getCurrentGene()
    genes <- getGeneList()
    gene2 <- getSecondGene()
    
    genelist <- unique(c(gene1,gene2,genes))
    
    genelist
  }
  )
    
  
  
  getCambridgeVariable <- reactive({
    input$clinvar_boxplot
  })
  
  
  getCambridgeZ <- reactive({
    input$z_cambridge
  })
  
  getCambridgeOverlay <- reactive({
    input$overlay_cambridge
  })
  
  getCambridgeCompositePlot <- reactive({
    input$cambridgeCombPlot
  })
  
  
  
  getDataset <- reactive({
    dataset <- input$theDataset
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
  
  
  ##########################################

    
  
  prepareExpressionMatrix <- reactive({
    
    message("Preparing expression matrix.....")
    
    
    dataset <- getDataset()
    genes <- getAllGenes()
    
    
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
  
  
  
  
  
  
  ##########################################
  
  prepareBoxplot <- reactive({
    
   data <- prepareExpressionMatrix() 
    
   plotType <- input$inputType
   
   if(plotType == "Single Gene") {
     genes <- getCurrentGene()
   } else {
     
     if (input$cambridgeCombPlot == "Yes"){
       genes <- getGeneList()
     }
     else genes <- input$cambridgeGeneChoice
   }
     
     

   
   
   #    data <- filterByGene(dataset, genes)
   data <- data %>% 
     filter(Symbol %in% genes)
   
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
   
   p1 <- p1 +  facet_wrap(~Symbol) 
   
   
   theme <- input$boxplotTheme
   
   p1 <- switch(theme,
                ggplot2 = p1,
                bw = p1 + theme_bw(),
                classic = p1 + theme_classic(),
                minimal = p1 + theme_minimal(),
                light = p1 + theme_light()
   )
   
                
   p1
   
   
    
  }
  )
  
  
  
  
  output$displayBoxplot <- renderPlot({
    
    message("Preparing the boxplot")
    plotType <- input$inputType
    
    data <- prepareExpressionMatrix() 
    message("Have data been re-calculated?")


    
    if(plotType == "Single Gene") {
      genes <- getCurrentGene()
    } else {
      
      if (input$cambridgeCombPlot == "Yes"){
        genes <- getGeneList()
      }
      else genes <- input$cambridgeGeneChoice
    }
    
    
    
    
    
    #    data <- filterByGene(dataset, genes)
    data <- data %>% 
      filter(Symbol %in% genes)
    
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
    
    p1 <- p1 +  facet_wrap(~Symbol) 
    
    
    theme <- input$boxplotTheme
    
    p1 <- switch(theme,
                 ggplot2 = p1,
                 bw = p1 + theme_bw(),
                 classic = p1 + theme_classic(),
                 minimal = p1 + theme_minimal(),
                 light = p1 + theme_light()
    )
    
    
    p1
    
    
    
    p1  
    
  }
  )
  
  
  output$anovaResult <- renderPrint({
    dataset <- getDataset()  
    plotType <- input$inputType
    
    
    if(plotType == "Single Gene") {
      genes <- getCurrentGene()
    } else  genes <- getGeneList()
    
    #    data <- filterByGene(dataset,genes)
    
    data <- prepareExpressionMatrix()
    
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
  
  output$cambridgeBoxplotPDF <- downloadHandler(
    filename = function(){
      paste0(input$profileBasename,".",input$profilePlotFormat)
    },
    content = function(file) {
      
      if(input$profilePlotFormat == "pdf") pdf(file, width=12, height=6.3)
      else png(file, width=1200,height=600)
      


        p1 <- prepareBoxplot()
        
        print(p1)
        dev.off()
      }


  )
  
  
  ##########################################
  
  

  
  output$rpSummary <- renderTable({
    
    combined.data <- prepareExpressionMatrix()
    
    combined.data <- filter(combined.data, !is.na(Event), !is.na(Time))
    
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
    
    
    combined.data <- prepareExpressionMatrix()
    
    combined.data <- filter(combined.data, !is.na(Event), !is.na(Time))
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
    combined.data <- prepareExpressionMatrix()
    combined.data <- filter(combined.data, !is.na(Event), !is.na(Time))
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
      combined.data <- prepareExpressionMatrix()
      combined.data <- filter(combined.data, !is.na(Event), !is.na(Time))
      
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
  
  
  ##########################################
  
  
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
  
  
  prepareHeatmap <- reactive({
    
    
    genes <- getGeneList()
    data <- prepareExpressionMatrix()
    
    if(dataset == "MSKCC"){
    
      taylor <- data
    
      samples <- filter(pd_taylor, Sample_Group == "prostate cancer") %>% 
        select(geo_accession) %>% as.matrix %>% as.character
      
      
      
      taylor <- filter(taylor,geo_accession %in% samples)
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
      
    }
    
    else if(dataset == "Cambridge"){
      
      camcap <- data
      
      summary_stats <- camcap %>% group_by(ID) %>% 
        summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
      

      
      
      
      
      samples <- filter(pd_camcap, Sample_Group == "Tumour") %>% 
        select(geo_accession) %>% as.matrix %>% as.character
      
      camcap <- filter(camcap, geo_accession %in% samples) %>% select(-Symbol)
      
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
      
      stockholm <- data
    
      samples <-  select(pd_stockholm,geo_accession) %>% as.matrix %>% as.character
      
      stockholm <- filter(stockholm,geo_accession %in% samples) 
      
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
  
  
  ##########################################
  
  
  getCopyNumberTable <- reactive({
    
      genes <- getCurrentGene()
      genes <- c(genes,getGeneList())
    
    message("Retrieving copy-number data....")
    
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
    
    if(input$inputType_cn == "Single Gene") {
      genes <- getCurrentGene()
      
      
    } else  genes <- getGeneList()
    
    cn.all <- filter(cn.all, Symbol %in% genes)
    
  })
  
  
  output$copyNumber <- renderPlot({
    
    cn.all <- getCopyNumberTable()
    if(input$inputType_cn == "Single Gene") {
      genes <- getCurrentGene()
      
      
    } else  genes <- getGeneList()
    
    cn.all <- filter(cn.all, Symbol %in% genes)
    
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
  
  
  
 
}
)