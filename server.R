library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(devtools)
library(gridExtra)
library(heatmap.plus)
library(RColorBrewer)
library(party)
library(survival)
library(knitr)
if(!require(prostateCancerTaylor)) install_github("markdunning/prostateCancerTaylor");library(prostateCancerTaylor)

taylor <- dplyrConvert()

##http://www.r-statistics.com/2013/07/creating-good-looking-survival-curves-the-ggsurv-function/
ggsurv <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                   cens.col = 'red', lty.est = 1, lty.ci = 2,
                   cens.shape = 3, back.white = F, xlab = 'Time',
                   ylab = 'Survival', main = ''){
  
  library(ggplot2)
  strata <- ifelse(is.null(s$strata) ==T, 1, length(s$strata))
  stopifnot(length(surv.col) == 1 | length(surv.col) == strata)
  stopifnot(length(lty.est) == 1 | length(lty.est) == strata)
  
  ggsurv.s <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = ''){
    
    dat <- data.frame(time = c(0, s$time),
                      surv = c(1, s$surv),
                      up = c(1, s$upper),
                      low = c(1, s$lower),
                      cens = c(0, s$n.censor))
    dat.cens <- subset(dat, cens != 0)
    
    col <- ifelse(surv.col == 'gg.def', 'black', surv.col)
    
    pl <- ggplot(dat, aes(x = time, y = surv)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(col = col, lty = lty.est)
    
    pl <- if(CI == T | CI == 'def') {
      pl + geom_step(aes(y = up), color = col, lty = lty.ci) +
        geom_step(aes(y = low), color = col, lty = lty.ci)
    } else (pl)
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  
  ggsurv.m <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = '') {
    n <- s$strata
    
    groups <- factor(unlist(strsplit(names
                                     (s$strata), '='))[seq(2, 2*strata, by = 2)])
    gr.name <-  unlist(strsplit(names(s$strata), '='))[1]
    gr.df <- vector('list', strata)
    ind <- vector('list', strata)
    n.ind <- c(0,n); n.ind <- cumsum(n.ind)
    for(i in 1:strata) ind[[i]] <- (n.ind[i]+1):n.ind[i+1]
    
    for(i in 1:strata){
      gr.df[[i]] <- data.frame(
        time = c(0, s$time[ ind[[i]] ]),
        surv = c(1, s$surv[ ind[[i]] ]),
        up = c(1, s$upper[ ind[[i]] ]),
        low = c(1, s$lower[ ind[[i]] ]),
        cens = c(0, s$n.censor[ ind[[i]] ]),
        group = rep(groups[i], n[i] + 1))
    }
    
    dat <- do.call(rbind, gr.df)
    dat.cens <- subset(dat, cens != 0)
    
    pl <- ggplot(dat, aes(x = time, y = surv, group = group)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(aes(col = group, lty = group))
    
    col <- if(length(surv.col == 1)){
      scale_colour_manual(name = gr.name, values = rep(surv.col, strata))
    } else{
      scale_colour_manual(name = gr.name, values = surv.col)
    }
    
    pl <- if(surv.col[1] != 'gg.def'){
      pl + col
    } else {pl + scale_colour_discrete(name = gr.name)}
    
    line <- if(length(lty.est) == 1){
      scale_linetype_manual(name = gr.name, values = rep(lty.est, strata))
    } else {scale_linetype_manual(name = gr.name, values = lty.est)}
    
    pl <- pl + line
    
    pl <- if(CI == T) {
      if(length(surv.col) > 1 && length(lty.est) > 1){
        stop('Either surv.col or lty.est should be of length 1 in order
             to plot 95% CI with multiple strata')
      }else if((length(surv.col) > 1 | surv.col == 'gg.def')[1]){
        pl + geom_step(aes(y = up, color = group), lty = lty.ci) +
          geom_step(aes(y = low, color = group), lty = lty.ci)
      } else{pl +  geom_step(aes(y = up, lty = group), col = surv.col) +
          geom_step(aes(y = low,lty = group), col = surv.col)}
    } else {pl}
    
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  pl <- if(strata == 1) {ggsurv.s(s, CI , plot.cens, surv.col ,
                                  cens.col, lty.est, lty.ci,
                                  cens.shape, back.white, xlab,
                                  ylab, main)
  } else {ggsurv.m(s, CI, plot.cens, surv.col ,
                   cens.col, lty.est, lty.ci,
                   cens.shape, back.white, xlab,
                   ylab, main)}
  pl
}

shinyServer(function(input, output){
  
  data <- reactive({inFile <- input$file1
  
                    if (is.null(inFile))
                    return(c("STAT3","ESR1","AR"))
  
  
                    print(inFile$name)
                    genes <- read.delim(inFile$datapath)[,1]
                    print(dim(genes))
                    genes
                    #read.csv("GraphPad Course Data/diseaseX.csv")
  })
  

  
  output$geneList <- renderPrint({
    genes <- data()
    genes
#    dput(df, file="data.rda")
  }
  )
  
  
  

  

output$rp <- renderPrint({
    
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
  



output$rp_plot <- reactivePlot(function(){
  
  genes <- data()
  
  pvalList       <- NA
  pvalList2      <- NA
  splitList      <- NA
  splitList2     <- NA
  highIsGoodList <- NA
  accList <- NA
  gg <- alist()
  
  # variable for how many plots have been created
  pCount <- 1
  
  
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
          gg[[pCount]] <- ggsurv(geneexp.survfit.xfs) + ggtitle(genes[i])
          newPval2                 <- NA
          pCount <- pCount +1 
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
        
      }
      
      if(newPval>=0.05) {
        ps       <- NA
        newPval2 <- NA
      }
      
      pvalList2[i]      <- newPval2
      splitList[i]      <- ps
      splitList2[i]     <- ps3
      highIsGoodList[i] <- highIsGood
      
    } else {
      accList[i] <- NA
      pvalList[i] <- NA
      pvalList2[i]      <- NA
      splitList[i]      <- NA
      splitList2[i]     <- NA
      highIsGoodList[i] <- NA
    }
    
  }
  
  pp <- do.call(grid.arrange,c(gg,ncol=2))
  pp   
})


  output$boxplot<- reactivePlot(function(){
    
    genes <- data()
    
    p1 <- taylor %>% filter(Gene %in% genes) %>% 
      filter(Sample_Group %in% c("prostate cancer","normal adjacent benign prostate")) %>% 
      ggplot(aes(x = Sample_Group, y = Expression, fill=Sample_Group)) + geom_boxplot() + facet_wrap(~Gene)
    
    p2 <- taylor %>% filter(Gene %in% genes) %>% 
      filter(Sample_Group %in% c("prostate cancer","normal adjacent benign prostate")) %>% 
      filter(!(is.na(Gleason))) %>% 
      ggplot(aes(x = Gleason, y = Expression, fill=Gleason)) + geom_boxplot() + facet_wrap(~Gene)
    
    grid.arrange(p1,p2)

##    if(input$showMu) p <- p + geom_vline(xintercept = mu,lty=2,col="red")
 #   print(p)
    
  }
  )
  
  output$heatmap<- reactivePlot(function(){
    
    genes <- data()
    
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
  

  output$mapping <- renderPrint({
    
    genes <- data()
    
    taylor %>% filter(Gene %in% genes) %>% 
    group_by(Gene) %>% 
    summarise(Number.Of.Probes = length(unique(Probe)))
    
  })
  
  output$summary <- renderPrint({
    genes <- data()
    
    taylor %>% filter(Gene %in% genes) %>% 
      filter(Sample_Group %in% c("prostate cancer","normal adjacent benign prostate")) %>% 
      group_by(Gene,Sample_Group) %>% 
      summarise(Mean = mean(Expression))

#    myDF <- taylor %>% filter(Gene %in% genes) %>% 
 #     filter(Sample_Group %in% c("prostate cancer","normal adjacent benign prostate"))
    
#    tmpRes <- lapply(split(myDF, unique(as.character(myDF$Gene))), function(x) t.test(x$Expression~x$Sample_Group))
#    testRes<- do.call("rbind",lapply(tmpRes, function(x) data.frame(Benign = x$estimate[1], Cancer = x$estimate[1], CI = paste('(',round(x$conf.int[[1]], 
 #                                                                                                              digits = 4),', ',
#                                                                                                     round(x$conf.int[[2]], 
#                                                                                                           digits = 4), ')',
 #                                                                                                    sep=""),pvalue = round(x$p.value, digits = 4),statistic = x$stat)))
      
  #  testRes
      
    
  })
  




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


output$thecode <- renderPrint({
  
  inFile <- input$file1
  
  print(as.name(paste0('myfile <- \"' , inFile$name,'\"')))
  
  print(as.name(paste0('sep <- \'', input$sep,'\'')))
  print(as.name(paste0('quote <- \'', input$quote,'\'')))
  print(as.name(paste('header <- ', input$header)))
  print(as.name(paste('skip <- ', input$skip)))
  print(as.name("data <- read.csv(myfile, header=header, sep=sep, quote=quote,skip=skip)"))
  
  #dump <- dput(data)
  #print(as.name(paste("data <-", capture.output(dput(data)))))
  print(as.name("head(data)"))
  
  print(as.name(paste("datacol <- ", input$dataCol)))
  print(as.name("X <- data[,datacol]"))
  print(as.name("summary(X)"))
  print(as.name("boxplot(X,horizontal=TRUE)"))
  
  print(as.name("colnames(data)[datacol] <- 'X'"))
  print(as.name("library(ggplot2)"))
  print(as.name("ggplot(data, aes(x=X)) + geom_histogram(aes(y=..density..),binwidth=.5,colour='black', fill='white')+ stat_function(fun=dnorm,color='red',arg=list(mean=mean(data$X), sd=sd(data$X)))"))
  
  
  print(as.name(paste0('alternative <- \'', input$alternative,'\'')))
  print(as.name(paste("mu <- ", input$mu)))
  print(as.name("t.test(X,mu=mu,alternative=alternative)"))
  print(as.name("sessionInfo()"))
}
)
  
  
}
)