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



#if(!require(prostateCancerTaylor)) install_github("crukci-bioinformatics/prostateCancerTaylor");library(prostateCancerTaylor)
#if(!require(prostateCancerCamcap)) install_github("crukci-bioinformatics/prostateCancerCamcap");library(prostateCancerCamcap)
#if(!require(prostateCancerStockholm)) install_github("crukci-bioinformatics/prostateCancerStockholm");library(prostateCancerStockholm)


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

## dplyr needs to be the last package loaded, otherwise the 'select' function seems to be over-written
library(dplyr)
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
  
  
  getCambridgeOverlay <- reactive({
    input$overlay_cambridge
  })
  getStockholmOverlay <- reactive({
    input$overlay_stockholm
  })
  getTaylorOverlay <- reactive({
    input$overlay_taylor
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

  
  getRpDataset <- reactive({
    input$rpDataset
  })
  
  getHeatmapDataset <- reactive({
    input$heatmapDataset
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
  
  output$boxplotCambridge <- reactivePlot(function(){
    
    currentGene <- getCurrentGene()
    message(paste("Plotting gene", currentGene))
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
    
    camcap <- full_join(camcap,pd_camcap)
    
    doZ <- ifelse(getCambridgeZ() == "Yes",TRUE,FALSE)
    
    if(doZ) camcap <- mutate(camcap, Expression=Z)
    
    var <- getCambridgeVariable()
    overlay <- getCambridgeOverlay()
    

    
    p1 <- switch(var,
                 iCluster = {camcap %>% 
                   filter(Sample_Group == "Tumour",!is.na(iCluster)) %>% 
                   ggplot(aes(x = iCluster, y = Expression, fill=iCluster)) + geom_boxplot() +  scale_fill_manual(values=iclusPal) + ggtitle(currentGene)
                   },
                 
                 Gleason = {camcap  %>% 
                   filter(Sample_Group == "Tumour") %>% 
                   filter(!(is.na(Gleason))) %>% 
                   ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot() +  ggtitle(currentGene)
                 }
                 )
    
    if(overlay == "Yes")  p1 <- p1 + geom_jitter(position=position_jitter(width = .05),alpha=0.75)    

      p1
  }
  )


output$boxplotStockholm <- reactivePlot(function(){
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


output$boxplotTaylor <- reactivePlot(function(){
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




output$anovaCambridge <- renderPrint({
  
  currentGene <- getCurrentGene()
  
  probes <- fd_camcap %>% filter(Symbol == currentGene) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
  
  cambridge <- exp_camcap  %>% filter(ID %in% probes) %>% 
    gather(geo_accession,Expression,-ID)
  
  summary_stats <- cambridge %>% group_by(ID) %>% 
    summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
  
  mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
  mu <- summary_stats$mean[which.max(summary_stats$iqr)]
  sd <- summary_stats$sd[which.max(summary_stats$iqr)]
  
  cambridge <- filter(cambridge, ID== mostVarProbe)
  cambridge <- full_join(cambridge,pd_camcap)
  
  var <- getCambridgeVariable()
  switch(var,
         iCluster = summary(aov(lm(Expression~iCluster,cambridge))),
         
         Gleason = summary(aov(lm(Expression~Gleason,cambridge)))
         
  )
  
  
}
)


output$anovaStockholm <- renderPrint({
  
  currentGene <- getCurrentGene()
  
  probes <- fd_stockholm %>% filter(Symbol == currentGene) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
  
  stockholm<- exp_stockholm  %>% filter(ID %in% probes) %>% 
    gather(geo_accession,Expression,-ID)
  
  summary_stats <- stockholm%>% group_by(ID) %>% 
    summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
  
  mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
  mu <- summary_stats$mean[which.max(summary_stats$iqr)]
  sd <- summary_stats$sd[which.max(summary_stats$iqr)]
  
  stockholm <- filter(stockholm, ID== mostVarProbe)
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
  
  probes <- fd_taylor %>% filter(Gene == currentGene) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
  
  taylor<- exp_taylor  %>% filter(ID %in% probes) %>% 
    gather(geo_accession,Expression,-ID)
  
  summary_stats <- taylor%>% group_by(ID) %>% 
    summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
  
  mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
  mu <- summary_stats$mean[which.max(summary_stats$iqr)]
  sd <- summary_stats$sd[which.max(summary_stats$iqr)]
  
  taylor <- filter(taylor, ID== mostVarProbe)
  taylor <- full_join(taylor,pd_taylor)
  
  var <- getTaylorVariable()
  switch(var,
         CopyNumberCluster = summary(aov(lm(Expression~Copy.number.Cluster,taylor))),
         Gleason = summary(aov(lm(Expression~Gleason,taylor)))
  )
  
  
}
)

output$rpPlot <- reactivePlot(function(){
  message("I am doing survival now....")
  currentGene <- getCurrentGene()
  dataset <- getRpDataset()
  
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
    
      probes <- fd_camcap %>% filter(Symbol == currentGene) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
      
      camcap<- exp_camcap  %>% filter(ID %in% probes) %>% 
        gather(geo_accession,Expression,-ID)
      
      summary_stats <- camcap %>% group_by(ID) %>% 
        summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
      
      mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
      mu <- summary_stats$mean[which.max(summary_stats$iqr)]
      sd <- summary_stats$sd[which.max(summary_stats$iqr)]
      
      camcap <- filter(camcap, ID== mostVarProbe)
      
      camcap <- full_join(camcap,pd_camcap) %>% 
        mutate(Time = FollowUpTime, Event = ifelse(BCR=="Y",1,0))
      
      combined.data <- camcap %>% 
        filter(!is.na(Time) & !is.na(Event))
    
    
  } else{
    
      probes <- fd_stockholm %>% filter(Symbol == currentGene) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
      
      stockholm <- exp_stockholm  %>% filter(ID %in% probes) %>% 
        gather(geo_accession,Expression,-ID)
      
      summary_stats <- stockholm %>% group_by(ID) %>% 
        summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
      
      mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
      mu <- summary_stats$mean[which.max(summary_stats$iqr)]
      sd <- summary_stats$sd[which.max(summary_stats$iqr)]
      
      stockholm <- filter(stockholm, ID== mostVarProbe)
      
      stockholm <- full_join(stockholm,pd_stockholm) %>% 
        mutate(Time = FollowUpTime, Event = ifelse(BCR=="Y",1,0))
      
      combined.data <- camcap %>% 
        filter(!is.na(Time) & !is.na(Event))
    
  }
  
  pvalList       <- NA
  pvalList2      <- NA
  splitList      <- NA
  splitList2     <- NA
  highIsGoodList <- NA
  accList <- NA

  i <- 1
  
    if(nrow(combined.data) > 0){
      
      surv.xfs <- Surv((as.numeric(as.character(combined.data$Time))/12), combined.data$Event)
      combined.data$surv.xfs <- surv.xfs
      ctree_xfs   <- ctree(surv.xfs~Expression,data=combined.data)
      pvalue      <- 1 - ctree_xfs@tree$criterion$maxcriterion
      newPval     <- signif(pvalue, digits = 2)
      pvalList[i] <- newPval
      ps3         <- NA
      highIsGood  <- NA
      accList[i] <- as.character(currentGene)
      
      if(newPval<0.05) {
        
        
        ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
        ps  <- signif(ps2[1], digits = 3)

        par(mfrow=c(2,1))
        plot(ctree(surv.xfs~Expression, data=combined.data))
        
        if(length(ps2)==1) {
          combined.data$geneexp_cp <- combined.data$Expression<=ps2[1]
          nt                       <- table(combined.data$geneexp_cp)
          geneexp.survfit.xfs      <- survfit(surv.xfs~combined.data$geneexp_cp)
#          plot(geneexp.survfit.xfs, xlab="Time to BCR (years)", ylab="Probability of Freedom from Biochemical Recurrence", main=paste(currentGene,", p=", newPval), col=c(2,4))
 #         legend("bottomleft", c(paste(currentGene, ">", ps, "n=", nt[[1]]), paste(currentGene, "<=", ps, "n=", nt[[2]])), col=c(2,4), lty=1, lwd=1.5, bty="n")
          newPval2                 <- NA
        }
        
      }
      
    }
})



output$survivalPlot <- reactivePlot(function(){
  message("I am doing survival now....")
  currentGene <- getCurrentGene()
  dataset <- getRpDataset()
  
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
    
    probes <- fd_camcap %>% filter(Symbol == currentGene) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
    
    camcap<- exp_camcap  %>% filter(ID %in% probes) %>% 
      gather(geo_accession,Expression,-ID)
    
    summary_stats <- camcap %>% group_by(ID) %>% 
      summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
    
    mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
    mu <- summary_stats$mean[which.max(summary_stats$iqr)]
    sd <- summary_stats$sd[which.max(summary_stats$iqr)]
    
    camcap <- filter(camcap, ID== mostVarProbe)
    
    camcap <- full_join(camcap,pd_camcap) %>% 
      mutate(Time = FollowUpTime, Event = ifelse(BCR=="Y",1,0))
    
    combined.data <- camcap %>% 
      filter(!is.na(Time) & !is.na(Event))
    
    
  } else{
    
    probes <- fd_stockholm %>% filter(Symbol == currentGene) %>% select(ID) %>% unique %>% as.matrix %>%  as.character
    
    stockholm <- exp_stockholm  %>% filter(ID %in% probes) %>% 
      gather(geo_accession,Expression,-ID)
    
    summary_stats <- stockholm %>% group_by(ID) %>% 
      summarise(mean=mean(Expression,na.rm=TRUE),sd=sd(Expression,na.rm=TRUE),iqr=IQR(Expression,na.rm=TRUE))
    
    mostVarProbe <- as.character(summary_stats$ID[which.max(summary_stats$iqr)])
    mu <- summary_stats$mean[which.max(summary_stats$iqr)]
    sd <- summary_stats$sd[which.max(summary_stats$iqr)]
    
    stockholm <- filter(stockholm, ID== mostVarProbe)
    
    stockholm <- full_join(stockholm,pd_stockholm) %>% 
      mutate(Time = FollowUpTime, Event = ifelse(BCR=="Y",1,0))
    
    combined.data <- camcap %>% 
      filter(!is.na(Time) & !is.na(Event))
    
  }
  
  pvalList       <- NA
  pvalList2      <- NA
  splitList      <- NA
  splitList2     <- NA
  highIsGoodList <- NA
  accList <- NA
  
  i <- 1
  
  if(nrow(combined.data) > 0){
    
    surv.xfs <- Surv((as.numeric(as.character(combined.data$Time))/12), combined.data$Event)
    combined.data$surv.xfs <- surv.xfs
    ctree_xfs   <- ctree(surv.xfs~Expression,data=combined.data)
    pvalue      <- 1 - ctree_xfs@tree$criterion$maxcriterion
    newPval     <- signif(pvalue, digits = 2)
    pvalList[i] <- newPval
    ps3         <- NA
    highIsGood  <- NA
    accList[i] <- as.character(currentGene)
    
    if(newPval<0.05) {
      
      
      ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
      ps  <- signif(ps2[1], digits = 3)
      
      
      if(length(ps2)==1) {
        combined.data$geneexp_cp <- combined.data$Expression<=ps2[1]
        nt                       <- table(combined.data$geneexp_cp)
        geneexp.survfit.xfs      <- survfit(surv.xfs~combined.data$geneexp_cp)
                  plot(geneexp.survfit.xfs, xlab="Time to BCR (years)", ylab="Probability of Freedom from Biochemical Recurrence", main=paste(currentGene,", p=", newPval), col=c(2,4))
                 legend("bottomleft", c(paste(currentGene, ">", ps, "n=", nt[[1]]), paste(currentGene, "<=", ps, "n=", nt[[2]])), col=c(2,4), lty=1, lwd=1.5, bty="n")
        newPval2                 <- NA
      }
      
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
      
      
      
    }
    hmcol <- brewer.pal(11 , "RdBu")
    heatmap.plus(geneMatrix,ColSideColors = colMatrix,col=hmcol)
    
  }
  )
  


output$downloadScript <- downloadHandler(
  filename = function() {
    paste(input$outfile, '.R', sep='')
  },
  content = function(file) {
    inFile <- input$file1
    
    cat(file=file,as.name(paste0('myfile <- \"' , inFile$name, '\"\n')))
    cat(file=file,as.name(paste0('sep <- \'', input$sep,'\'','\n')),append=TRUE)
    cat(file=file,as.name(paste0('quote <- \'', input$quote,'\'','\n')),append=TRUE)
    cat(file=file,as.name(paste('header <- ', input$header,'\n')),append=TRUE)
    cat(file=file,as.name(paste('skip <- ', input$skip,'\n')),append=TRUE)
    cat(file=file,as.name("data <- read.csv(myfile, header=header, sep=sep, quote=quote,skip=skip)\n"),append=TRUE)
    
    cat(file=file,as.name("head(data)\n"),append=TRUE)
  
    cat(file=file,as.name(paste("datacol <- ", input$dataCol,'\n')),append=TRUE)
    cat(file=file,as.name("X <- data[,datacol]\n"),append=TRUE)
    cat(file=file,as.name("summary(X)\n"),append=TRUE)
    cat(file=file,as.name("boxplot(X,horizontal=TRUE)\n"),append=TRUE)
    
    cat(file=file,as.name("colnames(data)[datacol] <- 'X'\n"),append=TRUE)
    cat(file=file, as.name("library(ggplot2)\n"),append=TRUE)
    cat(file=file, as.name("ggplot(data, aes(x=X)) + geom_histogram(aes(y=..density..),binwidth=.5,colour='black', fill='white')+ stat_function(fun=dnorm,color='red',arg=list(mean=mean(data$X), sd=sd(data$X)))\n"),append=TRUE)
    
    cat(file=file,as.name(paste0('alternative <- \'', input$alternative,'\'','\n')),append=TRUE)
    cat(file=file,as.name(paste("mu <- ", input$mu,'\n')),append=TRUE)
    cat(file=file,as.name("t.test(X,mu=mu,alternative=alternative)\n"),append=TRUE)
    cat(file=file,as.name("sessionInfo()\n"),append=TRUE)
    #formatR::tidy_source(source=file,output = file)
  }
)


output$downloadMarkdown <- downloadHandler(
  filename = function() {
    paste(input$outfile, '.Rmd', sep='')
  },
  content = function(file) {
    inFile <- input$file1
    script <- gsub(".Rmd", ".R",file)
    cat(file=script,as.name(paste0('myfile <- \"' , inFile$name, '\"\n')))
    cat(file=script,as.name(paste0('sep <- \'', input$sep,'\'','\n')),append=TRUE)
    cat(file=script,as.name(paste0('quote <- \'', input$quote,'\'','\n')),append=TRUE)
    cat(file=script,as.name(paste('header <- ', input$header,'\n')),append=TRUE)
    cat(file=script,as.name(paste('skip <- ', input$skip,'\n')),append=TRUE)
    cat(file=script,as.name("data <- read.csv(myfile, header=header, sep=sep, quote=quote,skip=skip)\n"),append=TRUE)
    cat(file=script,as.name("head(data)\n"),append=TRUE)
    
    cat(file=script,as.name(paste("datacol <- ", input$dataCol,'\n')),append=TRUE)
    cat(file=script,as.name("X <- data[,datacol]\n"),append=TRUE)
    cat(file=script,as.name("summary(X)\n"),append=TRUE)
    cat(file=script,as.name("boxplot(X,horizontal=TRUE)\n"),append=TRUE)
    cat(file=script,as.name("colnames(data)[datacol] <- 'X'\n"),append=TRUE)
    cat(file=script, as.name("library(ggplot2)\n"),append=TRUE)
    cat(file=script, as.name("ggplot(data, aes(x=X)) + geom_histogram(aes(y=..density..),binwidth=.5,colour='black', fill='white')+ stat_function(fun=dnorm,color='red',arg=list(mean=mean(data$X), sd=sd(data$X)))\n"),append=TRUE)
    
    cat(file=script,as.name(paste0('alternative <- \'', input$alternative,'\'','\n')),append=TRUE)
    cat(file=script,as.name(paste("mu <- ", input$mu,'\n')),append=TRUE)
    cat(file=script,as.name("t.test(X,mu=mu,alternative=alternative)\n"),append=TRUE)
    cat(file=script,as.name("sessionInfo()\n"),append=TRUE)
    knitr:::spin(hair=script,knit = FALSE)
    rmd <- readLines(file)
    
    cat(file = file, paste(input$title, "\n=======================\n"))
    cat(file=file, as.name(paste("###", input$name, "\n")),append=TRUE)    
    cat(file=file, as.name(paste("### Report Generated at: ", as.character(Sys.time()), "\n")),append=TRUE)    
    
    for(i in 1:length(rmd)){
      cat(file=file, as.name(paste(rmd[i], "\n")),append=TRUE)
      
    }
    
    #    formatR::tidy_urce(file,output = file)
  }
)


  
}
)