library(shiny)
library(ggplot2)

library(tidyr)
library(devtools)
library(gridExtra)
library(heatmap.plus)
library(RColorBrewer)
library(party)
library(survival)
library(knitr)
library(Biobase)
library(broom)


#if(!require(prostateCancerTaylor)) install_github("crukci-bioinformatics/prostateCancerTaylor");library(prostateCancerTaylor)
#if(!require(prostateCancerCamcap)) install_github("crukci-bioinformatics/prostateCancerCamcap");library(prostateCancerCamcap)
#if(!require(prostateCancerStockholm)) install_github("crukci-bioinformatics/prostateCancerStockholm");library(prostateCancerStockholm)
library(dplyr)

data(camcap,package = "prostateCancerCamcap")
pd_camcap <- tbl_df(pData(camcap))
fd_camcap <- tbl_df(fData(camcap))
exp_camcap <- tbl_df(data.frame(ID=as.character(featureNames(camcap)),exprs(camcap)))

data(stockholm, package="prostateCancerStockholm")
pd_stockholm <- tbl_df(pData(stockholm))
fd_stockholm <- tbl_df(fData(stockholm))
exp_stockholm <- tbl_df(data.frame(ID=as.character(featureNames(stockholm)),exprs(stockholm)))

data(taylor, package="prostateCancerTaylor")
pd_taylor <- tbl_df(pData(taylor))
fd_taylor <- tbl_df(fData(taylor))
exp_taylor <- tbl_df(data.frame(ID=as.character(featureNames(taylor)),log2(exprs(taylor))))

data(varambally,package="prostateCancerVarambally")
pd_varambally <- tbl_df(pData(varambally))
fd_varambally <- tbl_df(fData(varambally))
exp_varambally <- tbl_df(data.frame(ID=as.character(featureNames(varambally)),log2(exprs(varambally))))


data(grasso,package="prostateCancerGrasso")
pd_grasso <- tbl_df(pData(grasso))
fd_grasso <- tbl_df(fData(grasso))
exp_grasso <- tbl_df(data.frame(ID=as.character(featureNames(grasso)),exprs(grasso)))



## dplyr needs to be the last package loaded, otherwise the 'select' function seems to be over-written

select <- dplyr::select
iclusPal <- brewer.pal(5, "Set1")

message("READY FOR INPUT")

shinyServer(function(input, output){
  
  getCurrentGene <- reactive({
    input$currentGene
  })
  
#  getCurrentGene <- reactive({
 #   input$currentGene_cambridge
#  })

 # getCurrentGene <- reactive({
#    input$currentGene_stockholm
#  })
  
#  getCurrentGene <- reactive({
#    input$currentGene_taylor
#  })
  
  
    
  getCambridgeVariable <- reactive({
    input$clinvar_cambridge
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
    input$corDataset
  })
  
  getSecondGene <- reactive({
    input$secondGene
    
  })
  
  getCorType <- reactive({
    input$corType
  })
  
  getGeneList <- reactive({inFile <- input$file1
  
    if (is.null(inFile))
      return(c("STAT3","ESR1","AR"))
    
    
    print(inFile$name)
    genes <- read.delim(inFile$datapath)[,1]
    print(dim(genes))
    genes
    #read.csv("GraphPad Course Data/diseaseX.csv")
  })
  
  
  getSelectedGeneCambridge <-reactive({
    
    currentGene <- getCurrentGene()
    probes <- fd_camcap %>% filter(Symbol == currentGene) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
    
    camcap<- exp_camcap  %>% filter(ID %in% probes) %>% 
      gather(geo_accession,Expression,-ID)
    
    summary_stats <- camcap%>% group_by(ID) %>% 
      summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
    
    mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
    mu <- summary_stats$mean[which.max(summary_stats$iqr)]
    sd <- summary_stats$sd[which.max(summary_stats$iqr)]
    
    camcap <- filter(camcap, ID== mostVarProbe) %>%
      mutate(Z = (Expression -mu) /sd)
    
    camcap <- left_join(camcap,pd_camcap)
    camcap <- left_join(camcap, fd_camcap)

    
    camcap
    
  })
  
  getSelectedGeneListCambridge <- reactive({
    genes <- getGeneList()
    
    probes <- fd_camcap %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
    
    #      taylor<- exp_taylor  %>% filter(ID %in% probes) %>% 
    #       gather(geo_accession,Expression,-ID)
    camcap <- exp_camcap  %>% filter(ID %in% probes) %>% 
      gather(geo_accession,Expression,-ID)
    
    summary_stats <- camcap %>% group_by(ID) %>% 
      summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
    
    camcap <- left_join(camcap,summary_stats) %>% mutate(Z = (Expression - mean) / sd)
    
    mostVarProbes <- left_join(summary_stats,fd_camcap) %>% 
      arrange(Symbol,desc(iqr)) %>% 
      distinct(Symbol) %>% 
      select(ID) %>%  as.matrix %>%  as.character
  
    camcap <- filter(camcap, ID %in% mostVarProbes)
    camcap <- left_join(camcap, select(fd_camcap, ID, Symbol))
    camcap <- left_join(camcap, pd_camcap)
    camcap    
  })
  
  output$boxplotCambridge <- reactivePlot(function(){
    
    plotType <- input$inputType_cambridge
    if(plotType == "Single Gene") {
      camcap <- getSelectedGeneCambridge()
    } else  camcap <- getSelectedGeneListCambridge()
        
    doZ <- ifelse(getCambridgeZ() == "Yes",TRUE,FALSE)
    
    if(doZ) camcap <- mutate(camcap, Expression=Z)
    
    var <- getCambridgeVariable()
    overlay <- getCambridgeOverlay()
  
    p1 <- switch(var,
                 iCluster = {camcap %>% 
                   filter(Sample_Group == "Tumour",!is.na(iCluster)) %>% 
                   ggplot(aes(x = iCluster, y = Expression, fill=iCluster)) + geom_boxplot() +  scale_fill_manual(values=iclusPal) 
                   },
                 
                 Gleason = {camcap  %>% 
                   filter(Sample_Group == "Tumour") %>% 
                   filter(!(is.na(Gleason))) %>% 
                   ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot() 
                 }
                 )
    
    if(overlay == "Yes")  p1 <- p1 + geom_jitter(position=position_jitter(width = .05),alpha=0.75) 
    
      p1 +  facet_wrap(~Symbol) 
  }
  )


getSelectedGeneStockholm <- reactive({
  
  currentGene <- getCurrentGene()
  
  message(paste("Plotting gene", currentGene))
  probes <- fd_stockholm %>% filter(Symbol == currentGene) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
  
  stockholm <- exp_stockholm  %>% filter(ID %in% probes) %>% 
    gather(geo_accession,Expression,-ID)
  
  summary_stats <- stockholm %>% group_by(ID) %>% 
    summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
  
  mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
  mu <- summary_stats$mean[which.max(summary_stats$iqr)]
  sd <- summary_stats$sd[which.max(summary_stats$iqr)]
  
  stockholm <- filter(stockholm, ID== mostVarProbe) %>%
    mutate(Z = (Expression -mu) /sd)
  
  stockholm <- full_join(stockholm,pd_stockholm)
  

  
  
})
  
output$boxplotStockholm <- reactivePlot(function(){

  currentGene <- getCurrentGene()
  
  stockholm <- getSelectedGeneStockholm()
  
  doZ <- ifelse(getStockholmZ() == "Yes",TRUE,FALSE)
  
  if(doZ) stockholm <- mutate(stockholm, Expression=Z)
  
  var <- getStockholmVariable()
  overlay <- getStockholmOverlay()
  
  p1 <- switch(var,
               iCluster = {stockholm %>% 
                   filter(!is.na(iCluster)) %>% 
                   ggplot(aes(x = iCluster, y = Expression, fill=iCluster)) + geom_boxplot() +  scale_fill_manual(values=iclusPal) +  ggtitle(currentGene)
               },
               
               Gleason = {stockholm  %>% 
                   filter(!(is.na(Gleason))) %>% 
                   ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot() +  ggtitle(currentGene)
               }
  )
  
  if(overlay == "Yes")  p1 <- p1 + geom_jitter(position=position_jitter(width = .05),alpha=0.75)    
  
  p1
  
}

  
)


getSelectedGeneTaylor <- reactive({
  currentGene <- getCurrentGene()
  message(paste("Plotting gene", currentGene))
  probes <- fd_taylor %>% filter(Gene == currentGene) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
  
  taylor <- exp_taylor  %>% filter(ID %in% probes) %>% 
    gather(geo_accession,Expression,-ID)
  
  summary_stats <- taylor %>% group_by(ID) %>% 
    summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
  
  mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
  mu <- summary_stats$mean[which.max(summary_stats$iqr)]
  sd <- summary_stats$sd[which.max(summary_stats$iqr)]
  
  taylor <- filter(taylor, ID== mostVarProbe) %>%
    mutate(Z = (Expression -mu) /sd)
  
  taylor <- full_join(taylor,pd_taylor)
  
  
  
})


output$boxplotTaylor <- reactivePlot(function(){
  currentGene <- getCurrentGene()
  
  taylor <- getSelectedGeneTaylor()
  
  doZ <- ifelse(getTaylorZ() == "Yes",TRUE,FALSE)
  
  if(doZ) taylor <- mutate(taylor, Expression=Z)
  
  var <- getTaylorVariable()
  overlay <- getTaylorOverlay()
  
  p1 <- switch(var,
               CopyNumberCluster = {taylor %>% 
                   filter(!is.na(Copy.number.Cluster)) %>% 
                   ggplot(aes(x = Copy.number.Cluster, y = Expression, fill=Copy.number.Cluster)) + geom_boxplot()  +  ggtitle(currentGene)
               },
               
               Gleason = {taylor  %>% 
                   filter(!(is.na(Gleason))) %>% 
                   ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot() +  ggtitle(currentGene)
               }
  )
  
  if(overlay == "Yes")  p1 <- p1 + geom_jitter(position=position_jitter(width = .05),alpha=0.75)    
  
  p1
  
}


)


output$boxplotVarambally <- reactivePlot(function(){
  currentGene <- getCurrentGene()
  
  message(paste("Plotting gene", currentGene))
  probes <- fd_varambally %>% filter(Symbol == currentGene) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
  
  varambally <- exp_varambally  %>% filter(ID %in% probes) %>% 
    gather(geo_accession,Expression,-ID)
  
  summary_stats <- varambally %>% group_by(ID) %>% 
    summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
  
  mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
  mu <- summary_stats$mean[which.max(summary_stats$iqr)]
  sd <- summary_stats$sd[which.max(summary_stats$iqr)]
  
  varambally <- filter(varambally, ID== mostVarProbe) %>%
    mutate(Z = (Expression -mu) /sd)
  
  varambally <- full_join(varambally,pd_varambally)
  
  doZ <- ifelse(getVaramballyZ() == "Yes",TRUE,FALSE)
  
  if(doZ) varambally <- mutate(varambally, Expression=Z)
  
  var <- getVaramballyVariable()
  overlay <- getVaramballyOverlay()
  
  p1 <- varambally %>% ggplot(aes(x = Sample_Group, y = Expression, fill=Sample_Group)) + geom_boxplot()  +  ggtitle(currentGene)
               
  
  
  if(overlay == "Yes")  p1 <- p1 + geom_jitter(position=position_jitter(width = .05),alpha=0.75)    
  
  p1
  
}


)


output$boxplotGrasso <- reactivePlot(function(){
  currentGene <- getCurrentGene()
  
  message(paste("Plotting gene", currentGene))
  probes <- fd_grasso %>% filter(GENE_SYMBOL == currentGene) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
  
  grasso <- exp_grasso  %>% filter(ID %in% probes) %>% 
    gather(geo_accession,Expression,-ID)
  
  summary_stats <- grasso %>% group_by(ID) %>% 
    summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
  
  mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
  mu <- summary_stats$mean[which.max(summary_stats$iqr)]
  sd <- summary_stats$sd[which.max(summary_stats$iqr)]
  
  grasso <- filter(grasso, ID== mostVarProbe) %>%
    mutate(Z = (Expression -mu) /sd)
  
  grasso <- full_join(grasso,pd_grasso)
  
  doZ <- ifelse(getGrassoZ() == "Yes",TRUE,FALSE)
  
  if(doZ) grasso <- mutate(grasso, Expression=Z)
  
  var <- getGrassoVariable()
  overlay <- getGrassoOverlay()
  
  p1 <- grasso %>% ggplot(aes(x = Group, y = Expression, fill=Group)) + geom_boxplot()  +  ggtitle(currentGene)
  
  
  
  if(overlay == "Yes")  p1 <- p1 + geom_jitter(position=position_jitter(width = .05),alpha=0.75)    
  
  p1
  
}


)




output$anovaCambridge <- renderPrint({
  
  plotType <- input$inputType_cambridge
  if(plotType == "Single Gene") {
    camcap <- getSelectedGeneCambridge()
  } else  camcap <- getSelectedGeneListCambridge()

  var <- getCambridgeVariable()
  switch(var,
         iCluster = group_by(camcap, Symbol)  %>% do(tidy(aov(Expression~iCluster,data=.))),
         
         Gleason = group_by(camcap, Symbol)  %>% do(tidy(aov(Expression~Gleason,data=.)))
         
  )
  
  
}
)


output$anovaStockholm <- renderPrint({
  
  currentGene <- getCurrentGene()
  
  stockholm <- getSelectedGeneStockholm()
  stockholm <- full_join(stockholm,pd_stockholm)
  
  var <- getStockholmVariable()
  switch(var,
         iCluster = summary(aov(lm(Expression~iCluster,stockholm))),
         
         Gleason = summary(aov(lm(Expression~Gleason,stockholm)))
         
  )
  
  
}
)


output$anovaTaylor <- renderPrint({
  
  currentGene <- getCurrentGene()
  
  taylor <- getSelectedGeneTaylor()
  
  taylor <- full_join(taylor,pd_taylor)
  
  var <- getTaylorVariable()
  switch(var,
         CopyNumberCluster = summary(aov(lm(Expression~Copy.number.Cluster,taylor))),
         Gleason = summary(aov(lm(Expression~Gleason,taylor)))
  )
  
  
}
)

prepareSurvival <- reactive({

  dataset <- getRpDataset()
  currentGene <- getCurrentGene()
  
  if(dataset == "MSKCC"){
    
    probes <- fd_taylor %>% filter(Gene == currentGene) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
    
    taylor<- exp_taylor  %>% filter(ID %in% probes) %>% 
      gather(geo_accession,Expression,-ID)
    
    summary_stats <- taylor %>% group_by(ID) %>% 
      summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
    
    mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
    mu <- summary_stats$mean[which.max(summary_stats$iqr)]
    sd <- summary_stats$sd[which.max(summary_stats$iqr)]
    
    taylor <- filter(taylor, ID== mostVarProbe)
    
    taylor <- full_join(taylor,pd_taylor)
    
    
    combined.data <- taylor %>% 
      filter(!is.na(Event) & !is.na(Time))
    
  } else if (dataset == "Cambridge"){
    
    camcap <- getSelectedGeneCambridge()
    
    camcap <- full_join(camcap,pd_camcap) %>% 
      mutate(Time = FollowUpTime, Event = ifelse(BCR=="Y",1,0))
    
    combined.data <- camcap %>% 
      filter(!is.na(Time) & !is.na(Event))
    
    
  } else{
    stockholm <- getSelectedGeneStockholm()
    
    stockholm <- full_join(stockholm,pd_stockholm) %>% 
      mutate(Time = FollowUpTime, Event = ifelse(BCR=="Y",1,0))
    
    combined.data <- stockholm %>% 
      filter(!is.na(Time) & !is.na(Event))
    
    
  }
  
  combined.data
  
})

output$rpSummary <- renderPrint({
  
  combined.data <- prepareSurvival()
  currentGene <- getCurrentGene()
  results <- rpAnalysis(combined.data)
  ctree_xfs <- results[[1]]
  newPval <- results[[2]]
  summary <- ifelse(newPval <0.05, ". The partitioning and survival curves will apear below", ". Unfortunately, no significant partitioning of the expression values could be found")
  paste("Recursive Partitioning was run on ", currentGene, "from the", getRpDatase(), "dataset",summary)
})


rpAnalysis <- function(combined.data){
  surv.xfs <- Surv((as.numeric(as.character(combined.data$Time))/12), combined.data$Event)
  combined.data$surv.xfs <- surv.xfs
  ctree_xfs   <- ctree(surv.xfs~Expression,data=combined.data)
  pvalue      <- 1 - ctree_xfs@tree$criterion$maxcriterion
  newPval     <- signif(pvalue, digits = 2)
  list(ctree_xfs, newPval,surv.xfs)
  
}

output$rpPlot <- reactivePlot(function(){

  combined.data <- prepareSurvival()
  surv.xfs <- Surv((as.numeric(as.character(combined.data$Time))/12), combined.data$Event)
  combined.data$surv.xfs <- surv.xfs
  
  results <- rpAnalysis(combined.data)
  ctree_xfs <- results[[1]]
  newPval <- results[[2]]

      if(newPval<0.05) {
      
      
      ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
      ps  <- signif(ps2[1], digits = 3)

      par(mfrow=c(2,1))
      plot(ctree(surv.xfs~Expression, data=combined.data))
    }
      
  
}
)



output$survivalPlot <- reactivePlot(function(){
  combined.data <- prepareSurvival()
  results <- rpAnalysis(combined.data)
  surv.xfs <- Surv((as.numeric(as.character(combined.data$Time))/12), combined.data$Event)
  combined.data$surv.xfs <- surv.xfs
  ctree_xfs <- results[[1]]
  newPval <- results[[2]]
  currentGene <- getCurrentGene()
  
    
    if(newPval<0.05) {
      
      
      ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
      ps  <- signif(ps2[1], digits = 3)
      
      
      if(length(ps2)==1) {
        combined.data$geneexp_cp <- combined.data$Expression<=ps2[1]
        nt                       <- table(combined.data$geneexp_cp)
        geneexp.survfit.xfs      <- survfit(surv.xfs~geneexp_cp,data=combined.data)
                  plot(geneexp.survfit.xfs, xlab="Time to BCR (years)", ylab="Probability of Freedom from Biochemical Recurrence", main=paste(currentGene,", p=", newPval), col=c(2,4))
                 legend("bottomleft", c(paste(currentGene, ">", ps, "n=", nt[[1]]), paste(currentGene, "<=", ps, "n=", nt[[2]])), col=c(2,4), lty=1, lwd=1.5, bty="n")
        newPval2                 <- NA
      }
      
    }

})




  
  output$heatmap<- reactivePlot(function(){
    
    genes <- getGeneList()
    dataset <- getHeatmapDataset()
    
    if(dataset == "MSKCC"){
      
      probes <- fd_taylor %>% filter(Gene %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
      
#      taylor<- exp_taylor  %>% filter(ID %in% probes) %>% 
 #       gather(geo_accession,Expression,-ID)
      taylor <- exp_taylor  %>% filter(ID %in% probes) %>% 
        gather(geo_accession,Expression,-ID)
      
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
        spread(geo_accession,Expression) %>% data.frame
  
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
      
      
    } else if(dataset == "Cambridge"){
      
      probes <- fd_camcap %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
      
      #      taylor<- exp_taylor  %>% filter(ID %in% probes) %>% 
      #       gather(geo_accession,Expression,-ID)
      camcap <- exp_camcap  %>% filter(ID %in% probes) %>% 
        gather(geo_accession,Expression,-ID)
      
      summary_stats <- camcap %>% group_by(ID) %>% 
        summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
      
      mostVarProbes <- left_join(summary_stats,fd_camcap) %>% 
        arrange(Symbol,desc(iqr)) %>% 
        distinct(Symbol) %>% 
        select(ID) %>%  as.matrix %>%  as.character
      
      
      
      
      samples <- filter(pd_camcap, Sample_Group == "Tumour") %>% 
                select(geo_accession) %>% as.matrix %>% as.character
      
      camcap <- filter(camcap, ID %in% mostVarProbes,geo_accession %in% samples) 
      
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
      
      
      


        
    } else{
      
      probes <- fd_stockholm %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
      
      #      taylor<- exp_taylor  %>% filter(ID %in% probes) %>% 
      #       gather(geo_accession,Expression,-ID)
      stockholm <- exp_stockholm  %>% filter(ID %in% probes) %>% 
        gather(geo_accession,Expression,-ID)
      
      summary_stats <- stockholm %>% group_by(ID) %>% 
        summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
      
      mostVarProbes <- left_join(summary_stats,fd_stockholm) %>% 
        arrange(Symbol,desc(iqr)) %>% 
        distinct(Symbol) %>% 
        select(ID) %>%  as.matrix %>%  as.character
      
      
      
      
      samples <-  select(pd_stockholm,geo_accession) %>% as.matrix %>% as.character
      
      stockholm <- filter(stockholm, ID %in% mostVarProbes,geo_accession %in% samples) 
      
      geneMatrix <- stockholm %>% 
        spread(geo_accession,Expression) %>% data.frame
      
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
      
      
    }
    
    
    hmcol <- brewer.pal(11 , "RdBu")
    
    if(getDistFun() == "Correlation") distfun <- function(x) as.dist(1 - cor(t(x)))
    else distfun <- dist
    
    hclustfun <- function(x) hclust(x,method=getHclustMethod())
    scale <- getScaleMethod()
    
    if(getReordRows() == "Yes") heatmap.plus(geneMatrix,ColSideColors = colMatrix,col=hmcol,distfun=distfun,hclustfun = hclustfun,scale=scale)
    else heatmap.plus(geneMatrix,ColSideColors = colMatrix,col=hmcol,distfun=distfun,Rowv=NA,hclustfun = hclustfun,scale=scale)
  }
  )
  


  output$corPlot <- reactivePlot(function(){
    
    dataset <- getCorDataset()
    currentGene <- getCurrentGene()
    secondGene <- getSecondGene()
    genes <- c(currentGene, secondGene)
    
    if(dataset == "MSKCC"){

      probes <- fd_taylor %>% filter(Gene %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
      
      #      taylor<- exp_taylor  %>% filter(ID %in% probes) %>% 
      #       gather(geo_accession,Expression,-ID)
      taylor <- exp_taylor  %>% filter(ID %in% probes) %>% 
        gather(geo_accession,Expression,-ID)
      
      summary_stats <- taylor %>% group_by(ID) %>% 
        summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
      
      mostVarProbes <- left_join(summary_stats,fd_taylor) %>% 
        arrange(Gene,desc(iqr)) %>% 
        distinct(Gene) %>% 
        select(ID) %>%  as.matrix %>%  as.character
      
      
      
      
      samples <- filter(pd_taylor, Sample_Group == "prostate cancer") %>% 
        select(geo_accession) %>% as.matrix %>% as.character
      
      taylor <- filter(taylor, ID %in% mostVarProbes,geo_accession %in% samples) 
      cor_data <- left_join(taylor,select(fd_taylor,ID,Gene)) %>% mutate(Gene = ifelse(Gene == currentGene,"Gene1","Gene2")) %>% 
        select(geo_accession,Expression,Gene) %>% spread(Gene,Expression)
      cor_data <- left_join(cor_data, select(pd_taylor, geo_accession,Copy.number.Cluster)) %>% mutate(Group = Copy.number.Cluster)
      
    } else if(dataset == "Cambridge"){
    

    
    probes <- fd_camcap %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
    
    #      taylor<- exp_taylor  %>% filter(ID %in% probes) %>% 
    #       gather(geo_accession,Expression,-ID)
    camcap <- exp_camcap  %>% filter(ID %in% probes) %>% 
      gather(geo_accession,Expression,-ID)
    
    summary_stats <- camcap %>% group_by(ID) %>% 
      summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
    
    mostVarProbes <- left_join(summary_stats,fd_camcap) %>% 
      arrange(Symbol,desc(iqr)) %>% 
      distinct(Symbol) %>% 
      select(ID) %>%  as.matrix %>%  as.character
    
    
    
    
    samples <- filter(pd_camcap, Sample_Group == "Tumour") %>% 
      select(geo_accession) %>% as.matrix %>% as.character
    
    camcap <- filter(camcap, ID %in% mostVarProbes,geo_accession %in% samples) 
    cor_data <- left_join(camcap,select(fd_camcap,ID,Symbol)) %>% mutate(Gene = ifelse(Symbol == currentGene,"Gene1","Gene2")) %>% 
      select(geo_accession,Expression,Gene) %>% spread(Gene,Expression)
    
    cor_data <- left_join(cor_data, select(pd_camcap, geo_accession,iCluster)) %>% mutate(Group = iCluster)
                     
    } else {
      
      probes <- fd_stockholm %>% filter(Symbol %in% genes) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
      
      #      taylor<- exp_taylor  %>% filter(ID %in% probes) %>% 
      #       gather(geo_accession,Expression,-ID)
      stockholm <- exp_stockholm  %>% filter(ID %in% probes) %>% 
        gather(geo_accession,Expression,-ID)
      
      summary_stats <- stockholm %>% group_by(ID) %>% 
        summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
      
      mostVarProbes <- left_join(summary_stats,fd_stockholm) %>% 
        arrange(Symbol,desc(iqr)) %>% 
        distinct(Symbol) %>% 
        select(ID) %>%  as.matrix %>%  as.character
      

      stockholm <- filter(stockholm, ID %in% mostVarProbes) 
      cor_data <- left_join(stockholm,select(fd_stockholm,ID,Symbol)) %>% mutate(Gene = ifelse(Symbol == currentGene,"Gene1","Gene2")) %>% 
        select(geo_accession,Expression,Gene) %>% spread(Gene,Expression)
      cor_data <- left_join(cor_data, select(pd_stockholm, geo_accession,iCluster)) %>% mutate(Group = iCluster)
      
    }
        
    cor <- round(with(cor_data, cor(Gene1,Gene2,method=getCorType())),3)
    ggplot(cor_data, aes(x = Gene1, y=Gene2,col=Group)) + geom_point() + xlab(currentGene) + ylab(secondGene) + ggtitle(paste0("Correlation = ", cor))
                 
  }
    )
  
output$cambridgeProfileScript <- downloadHandler(
  filename = function() {
    paste(input$outfile, '.R', sep='')
  },
  content = function(file) {
    inFile <- input$file1
    cat(file=file,as.name("##Load the required libraries\n"))
    cat(file=file,as.name("library(ggplot2)\n"),append=TRUE)
    cat(file=file,as.name("library(tidyr)\n"),append=TRUE)
    cat(file=file,as.name("library(devtools)\n"),append=TRUE)
    
    cat(file=file,as.name("###Install camcap dataset if not present\n"),append=TRUE)
    cat(file=file,as.name("if(!require(prostateCancerCamcap)) install_github('crukci-bioinformatics/prostateCancerCamcap')\n"),append=TRUE)
    
    cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
    cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
    cat(file=file,as.name("data(camcap,package = 'prostateCancerCamcap')\n"),append=TRUE)
    cat(file=file,as.name("pd_camcap <- tbl_df(pData(camcap))\n"),append=TRUE)
    cat(file=file,as.name("fd_camcap <- tbl_df(fData(camcap))\n"),append=TRUE)
    cat(file=file,as.name("exp_camcap <- tbl_df(data.frame(ID=as.character(featureNames(camcap)),exprs(camcap)))\n"),append=TRUE)

    
    cat(file=file,as.name("library(RColorBrewer)\n"),append=TRUE)
    cat(file=file,as.name("iclusPal <- brewer.pal(5, 'Set1')\n"),append=TRUE)
    
    
    currentGene <- getCurrentGene()
    doZ <- ifelse(getCambridgeZ() == "Yes",TRUE,FALSE)
    var <- getCambridgeVariable()
    overlay <- getCambridgeOverlay()
    
    cat(file=file,as.name(paste0("currentGene <-'", currentGene,"'\n")),append=TRUE)
    
    cat(file=file, as.name(paste0("probes <- fd_camcap %>% filter(Symbol == \'",currentGene,"\') %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n")),append=TRUE)
    
    cat(file=file,as.name("camcap<- exp_camcap  %>% filter(ID %in% probes) %>% gather(geo_accession,Expression,-ID)\n"),append=TRUE)
    cat(file=file,as.name("summary_stats <- camcap%>% group_by(ID) %>% summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))\n"),append=TRUE)
    cat(file=file,as.name("mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])\n"),append=TRUE)
    cat(file=file,as.name("mu <- summary_stats$mean[which.max(summary_stats$iqr)]\n"),append=TRUE)
    cat(file=file,as.name("sd <- summary_stats$sd[which.max(summary_stats$iqr)]\n"),append=TRUE)
    
    cat(file=file,as.name("camcap <- filter(camcap, ID== mostVarProbe) %>% mutate(Z = (Expression -mu) /sd)\n"),append=TRUE)
    cat(file=file,as.name("camcap <- full_join(camcap,pd_camcap)\n"),append=TRUE)
    
    if(doZ) cat(file=file,as.name("camcap <- mutate(camcap, Expression=Z)\n"),append=TRUE)
    
    if(var == "iCluster"){
      
      cat(file=file,as.name("p1 <- camcap %>% filter(Sample_Group == 'Tumour',!is.na(iCluster)) %>% ggplot(aes(x = iCluster, y = Expression, fill=iCluster)) + geom_boxplot() + ggtitle(currentGene)+ scale_fill_manual(values=iclusPal)\n"),append=TRUE)
      
    } else if(var == "Gleason"){
      
      cat(file=file,as.name("p1 <- camcap %>% filter(Sample_Group == 'Gleason',!is.na(iCluster)) %>% ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot() +  ggtitle(currentGene)\n"))
    }
    
    
    if(overlay=="Yes") cat(file=file,as.name("p1 <- p1 + geom_jitter(position=position_jitter(width = .05),alpha=0.75) \n"),append=TRUE)
    
    cat(file=file,as.name("p1"),append=TRUE)
   }
)


output$stockholmProfileScript <- downloadHandler(
  filename = function() {
    paste(input$outfile, '.R', sep='')
  },
  content = function(file) {
    inFile <- input$file1
    cat(file=file,as.name("##Load the required libraries\n"))
    cat(file=file,as.name("library(ggplot2)\n"),append=TRUE)
    cat(file=file,as.name("library(tidyr)\n"),append=TRUE)
    cat(file=file,as.name("library(devtools)\n"),append=TRUE)
    
    cat(file=file,as.name("###Install stockholm dataset if not present\n"),append=TRUE)
    cat(file=file,as.name("if(!require(prostateCancerStockholm)) install_github('crukci-bioinformatics/prostateCancerStockholm')\n"),append=TRUE)
    
    cat(file=file,as.name("library(dplyr)\n"),append=TRUE)
    cat(file=file,as.name("###Convert into data convenient for dplyr\n"),append=TRUE)
    cat(file=file,as.name("data(stockholm,package = 'prostateCancerStockholm')\n"),append=TRUE)
    cat(file=file,as.name("pd_stockholm <- tbl_df(pData(stockholm))\n"),append=TRUE)
    cat(file=file,as.name("fd_stockholm <- tbl_df(fData(stockholm))\n"),append=TRUE)
    cat(file=file,as.name("exp_stockholm <- tbl_df(data.frame(ID=as.character(featureNames(stockholm)),exprs(stockholm)))\n"),append=TRUE)
    
    
    cat(file=file,as.name("library(RColorBrewer)\n"),append=TRUE)
    cat(file=file,as.name("iclusPal <- brewer.pal(5, 'Set1')\n"),append=TRUE)
    
    
    currentGene <- getCurrentGene()
    doZ <- ifelse(getStockholmZ() == "Yes",TRUE,FALSE)
    var <- getStockholmVariable()
    overlay <- getStockholmOverlay()
    
    cat(file=file,as.name(paste0("currentGene <-'", currentGene,"'\n")),append=TRUE)
    
    cat(file=file, as.name(paste0("probes <- fd_stockholm %>% filter(Symbol == \'",currentGene,"\') %>% select(ID) %>% unique %>% as.matrix %>%  as.character\n")),append=TRUE)
    
    cat(file=file,as.name("stockholm <- exp_stockholm  %>% filter(ID %in% probes) %>% gather(geo_accession,Expression,-ID)\n"),append=TRUE)
    cat(file=file,as.name("summary_stats <- stockholm %>% group_by(ID) %>% summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))\n"),append=TRUE)
    cat(file=file,as.name("mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])\n"),append=TRUE)
    cat(file=file,as.name("mu <- summary_stats$mean[which.max(summary_stats$iqr)]\n"),append=TRUE)
    cat(file=file,as.name("sd <- summary_stats$sd[which.max(summary_stats$iqr)]\n"),append=TRUE)
    
    cat(file=file,as.name("stockholm <- filter(stockholm, ID== mostVarProbe) %>% mutate(Z = (Expression -mu) /sd)\n"),append=TRUE)
    cat(file=file,as.name("stockholm <- full_join(stockholm,pd_stockholm)\n"),append=TRUE)
    
    if(doZ) cat(file=file,as.name("stockholm <- mutate(stockholm, Expression=Z)\n"),append=TRUE)
    
    if(var == "iCluster"){
      
      cat(file=file,as.name("p1 <- stockholm %>% filter(!is.na(iCluster)) %>% ggplot(aes(x = iCluster, y = Expression, fill=iCluster)) + geom_boxplot() + ggtitle(currentGene)+ scale_fill_manual(values=iclusPal)\n"),append=TRUE)
      
    } else if(var == "Gleason"){
      
      cat(file=file,as.name("p1 <- stockholm %>% ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot() +  ggtitle(currentGene)\n"))
    }
    
    
    if(overlay=="Yes") cat(file=file,as.name("p1 <- p1 + geom_jitter(position=position_jitter(width = .05),alpha=0.75) \n"),append=TRUE)
    
    cat(file=file,as.name("p1"),append=TRUE)
  }
)


  
}
)