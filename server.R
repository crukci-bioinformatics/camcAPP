library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(devtools)
library(gridExtra)
library(heatmap.plus)
library(RColorBrewer)


if(!require(prostateCancerTaylor)) install_github("markdunning/prostateCancerTaylor");library(prostateCancerTaylor)

taylor <- dplyrConvert()

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
  
#  output$plot <- renderPlot({
#    plot(data(), xlab="X", ylab="Y", ylim=c(-300,800))
#    if(input$line) {
#      abline(lm(Y ~ X, data=data()), col="dark blue")
#    }
#    if(input$means) {
#      abline(v = mean(data()[,1]), lty="dotted")
#      abline(h = mean(data()[,2]), lty="dotted")
#    } 
#    if(input$ant) {
#      model = lm(Y ~ X, data=data())
#      txt = paste("The equation of the line is:\nY = ",
#                  round(coefficients(model)[1],0)," + ",
#                  round(coefficients(model)[2],3),"X + error")
      
#      boxed.labels(50,600,labels=txt,bg="white", cex=1.25)
#    }    
    
#  })
 #
  
  output$geneList <- renderPrint({
    genes <- data()
    genes
#    dput(df, file="data.rda")
  }
  )
  


  
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