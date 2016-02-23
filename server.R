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

library(dplyr)

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


iclusPal <- brewer.pal(5, "Set1")

message("READY FOR INPUT")

shinyServer(function(input, output){
  
  getCambridgeGene <- reactive({
    input$currentGene_cambridge
  })

  getStockholmGene <- reactive({
    input$currentGene_stockholm
  })
  
  getTaylorGene <- reactive({
    input$currentGene_taylor
  })
  
  
    
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

  
  
  
  
  
  output$boxplotCambridge <- reactivePlot(function(){
    
    currentGene <- getCambridgeGene()
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
  currentGene <- getStockholmGene()
  
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
  currentGene <- getTaylorGene()
  
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
  
  currentGene <- getCambridgeGene()
  
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
  
  currentGene <- getStockholmGene()
  
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
  
  currentGene <- getTaylorGene()
  
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
  
  var <- getStockholmVariable()
  switch(var,
         CopyNumberCluster = summary(aov(lm(Expression~Copy.Number.Cluster,taylor))),
         
         Gleason = summary(aov(lm(Expression~Gleason,taylor)))
         
  )
  
  
}
)

output$rpCambridge <- renderPrint({
    
    genes <- data()
    print(genes)
    
    pvalList       <- NA
    pvalList2      <- NA
    splitList      <- NA
    splitList2     <- NA
    highIsGoodList <- NA
    accList <- NA
    gg <- alist()

    for(i in 1:length(genes)){
      combined.data <- taylor %>% filter(Gene %in% genes[i]) %>% 
        filter(!is.na(Event) & !is.na(Time))
      
      
      if(nrow(combined.data) > 0){
        surv.xfs <- Surv((combined.data$Time/12), as.numeric(combined.data$Event))
        combined.data$surv.xfs <- surv.xfs

        ctree_xfs   <- ctree(surv.xfs~Expression,data=combined.data)
        pvalue      <- 1 - ctree_xfs@tree$criterion$maxcriterion
        newPval     <- signif(pvalue, digits = 2)
        pvalList[i] <- newPval
        print(newPval)
        ps3         <- NA
        highIsGood  <- NA
        accList[i] <- as.character(combined.data$GB_ACC[1])
        
        if(newPval<0.05) {

          ps2 <- party:::cutpoints_list(ctree_xfs@tree, variableID=1)
          ps  <- signif(ps2[1], digits = 3)
          
          if(length(ps2)==1) {
            combined.data$geneexp_cp <- combined.data$Expression<=ps2[1]
            nt                       <- table(combined.data$geneexp_cp)
            geneexp.survfit.xfs      <- survfit(surv.xfs~combined.data$geneexp_cp)
            gg[[i]] <- ggsurv(geneexp.survfit.xfs) + ggtitle(genes[i])
            newPval2                 <- NA
            
          }
          
          if(length(ps2)==2) {
            if(ps2[1]>ps2[2]) {
              ps3                       <- round(ps2[2], digits=3)
              combined.data$geneexp_cp1 <- combined.data$Expression<=ps2[1]
              combined.data$geneexp_cp2 <- combined.data$Expression<=ps2[2]
              combined.data$geneexp_cp3 <- combined.data$geneexp_cp1+combined.data$geneexp_cp2
              nt                        <- table(combined.data$geneexp_cp3)
              geneexp.survfit.xfs       <- survfit(surv.xfs~combined.data$geneexp_cp3)
              pvalue2                   <- 1 - ctree_xfs@tree$left[[3]][[2]]
              newPval2                  <- signif(pvalue2, digits=2)
              
            } else {
              ps3                       <- round(ps2[2], digits=3)
              combined.data$geneexp_cp1 <- combined.data$geneexp<=ps2[1]
              combined.data$geneexp_cp2 <- combined.data$geneexp<=ps2[2]
              combined.data$geneexp_cp3 <- combined.data$geneexp_cp1+combined.data$geneexp_cp2
              nt                        <- table(combined.data$geneexp_cp3)
              geneexp.survfit.xfs       <- survfit(surv.xfs~combined.data$geneexp_cp3)
              pvalue2                   <- 1 - ctree_xfs@tree$right[[3]][[2]]
              newPval2                  <- signif(pvalue2, digits=2)
              
              
            }
          }
          
          combined.data$geneexp_cp <- combined.data$geneexp <= ps2[1]
          coxphRegressionModel     <- coxph(surv.xfs~combined.data$geneexp_cp)
          hazardRatio              <- exp(coxphRegressionModel$coefficients)
          highIsGood               <- hazardRatio > 1
          
        } else {
          ps       <- NA
          newPval2 <- NA
        }
        
        
        
        
        
        pvalList2[i]      <- newPval2
        splitList[i]      <- ps
        splitList2[i]     <- ps3
        highIsGoodList[i] <- highIsGood
        
      }   else {
        accList[i] <- NA
        pvalList[i] <- NA
        pvalList2[i]      <- NA
        splitList[i]      <- NA
        splitList2[i]     <- NA
        highIsGoodList[i] <- NA
      }
      
      
    }
    
    Accession  <- accList
    PValue1    <- pvalList
    PValue2    <- pvalList2
    CutOff1    <- splitList
    CutOff2    <- splitList2
    HighIsGood <- highIsGoodList
    
    xfsList <- data.frame(GeneName= as.character(genes), Accession = accList,CutOff1, CutOff2, PValue1, PValue2, HighIsGood)
    kable(xfsList)
    
    
    
  })
  



output$rp_plotCambridge <- reactivePlot(function(){
  
  currentGene <- getCambridgeGene()
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
  

  combined.data <- camcap %>% 
    filter(!is.na(BCR) & !is.na(FollowUpTime)) %>% mutate(BCR=ifelse(BCR=="Y",1,0))
    
  
  
  pvalList       <- NA
  pvalList2      <- NA
  splitList      <- NA
  splitList2     <- NA
  highIsGoodList <- NA
  accList <- NA
  gg <- alist()
  
  i <- 1
  
    if(nrow(combined.data) > 0){
      
      surv.xfs <- Surv((as.numeric(as.character(combined.data$FollowUpTime))/12), combined.data$BCR)
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
          plot(geneexp.survfit.xfs, xlab=xlabel, ylab=ylabel, main=paste(geneid,", p=", newPval), col=c(2,4))
          newPval2                 <- NA
          pCount <- pCount +1 
        }
        
      }
})



  
  output$heatmap<- reactivePlot(function(){
    
    genes <- data()
    
    taylor <- taylor %>% filter(Gene %in% genes)
    
    
    varProbes <- taylor %>% group_by(Gene) %>% 
      summarise(Probe = Probe[which.max(IQR)])
    
    taylor <- inner_join(taylor, varProbes,by="Probe") %>% rename(Gene = Gene.x) %>% select(-c(Gene.y))
    
    geneMatrix <- taylor %>% filter(Gene %in% genes) %>% 
      filter(Sample_Group %in% c("prostate cancer","normal adjacent benign prostate")) %>% 
    select(Expression, Sample,Gene) %>% 
    spread(Sample, Expression) %>% data.frame
    
    rownames(geneMatrix) <- geneMatrix$Gene
    geneMatrix <- as.matrix(geneMatrix[,-1])
    
    colMatrix <- matrix(nrow = ncol(geneMatrix),ncol = 2)
    grp <- taylor$Sample_Group[match(colnames(geneMatrix),taylor$Sample)]
    cols <- brewer.pal(9,"Set1")[1:length(levels(factor(as.character(grp))))]
    
    grp <- as.factor(as.character(grp))
    levels(grp) <- cols
    
    colMatrix[,1] <- as.character(grp)

    grp <- taylor$Gleason[match(colnames(geneMatrix),taylor$Sample)]
    cols <- brewer.pal(8,"Set2")[1:length(levels(factor(as.character(grp))))]
    
    grp <- as.factor(as.character(grp))
    levels(grp) <- cols
    colMatrix[,2] <- as.character(grp)
        
    heatmap.plus(geneMatrix,ColSideColors = colMatrix)
    
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