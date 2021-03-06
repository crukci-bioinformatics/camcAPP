

library(shiny)


shinyUI(navbarPage("Explore Prostate Cancer Datasets", id = "nav",
                  
    

    
    
        tabPanel("Data Input",
                 sidebarLayout(
                   sidebarPanel(
                     tags$head(includeScript("google-analytics.js")),           
                   
                     #),
  #                   radioButtons("inputType", "Use Single or Gene List as input?", choices=c("Single Gene","Gene List"),selected="Single Gene"),
                      img(src="camcAPP_logo_Final_(no_background).png",width="90%"),
                      fileInput('file1', 'Gene List',
                               accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),helpText("Your gene list must tab-delimited, with gene names in the first column"),
                      a("Click Here to downlad an example gene list",href="https://raw.githubusercontent.com/crukci-bioinformatics/camcAPP/master/example-genes.txt"),
                      helpText("If no gene list is uploaded, the genes AR, ESR1, HES6 MELK and STAT3 will be used"),
   #                   textInput("currentGene", "Gene of Interest", value="A1BG"),
                      helpText("If you want to analyse a single-gene, see the Quick Analysis tab"),
                      selectInput("theDataset","Choose a Dataset",choices=c("Cambridge","Stockholm","MSKCC", "Michigan2005","Michigan2012"),selected = "Cambridge"),
                      h2("User Guide"),
                      tags$a(href="https://github.com/crukci-bioinformatics/camcAPP/raw/master/supplementary/supplementary.pdf",class="btn",icon("download"),"Download User Guide"),
                      h2("Citation"),
                      strong("If you use any of the images generated in a publication or presentation, please cite:"),
                      em("Dunning et al (2017). Mining the human prostate cancer datasets: The camcApp Shiny app. EBiomedicine. http://dx.doi.org/10.1016/j.ebiom.2017.02.022"),br(),
                      strong("Please also cite the relevant publication for the dataset;"),br(),
                      h3("Cambridge and Stockholm"),
                      em("Integration of copy number and transcriptomics provides risk stratification in prostate cancer: A discovery and validation cohort study. Ross-Adams et al. (2015) doi:10.1016/j.ebiom.2015.07.017"),
                      h3("MSKCC"),
                      em("Integrative genomic profiling of human prostate cancer. Taylor et al. (2010) doi:10.1016/j.ccr.2010.05.026"),
                      h3("Michigan2005"),
                      em("Integrative genomic and proteomic analysis of prostate cancer reveals signatures of metastatic progression. Varambally et al. (2005) doi:10.1016/j.ccr.2005.10.001"),
                      h3("Michigan2012"),
                      em("The Mutational Landscape of Lethal Castrate Resistant Prostate Cancer. Grasso et al. (2012) doi:10.1038/nature11125")
                   )
                  ,
                   mainPanel(
                     h3("About this app............."),
                     img(src="cruk-cambridge-institute.jpg",width="90%"), br(),a("cruk.cam.ac.uk",href="www.cruk.cam.ac.uk"),
                     helpText("This app was developed by"),a("Cancer Research Uk Cambridge Institute",href="http://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core"),
                     helpText("Source code available on github;"),code("https://github.com/crukci-bioinformatics/camcAPP"),
                     h3("Features"),
                     h2("Gene Profile"),
                     helpText("Produces boxplots to visualise the distribution of the selected genes. For each dataset, you can choose which clinical variable to group the samples on"),
                     helpText("When choosing Cambridge or Stockholm, you will have the option to display the expression in the five different subtypes identified by Ross-Adams et al (2015). These subtypes were shown to have significantly different outcomes"),
                     img(src="km-from-paper.png",width="100%"),
                     helpText("If multiple microarray probes are found for the gene, the probe with the highest inter-quartile range (IQR) will be picked"),
                     helpText("An ANOVA analysis will also be performed to assess whether there are different expression levels in the groups you have chosen"),
                     helpText("The boxplot can be exported as a pdf or png image. An R script can be downloaded, allowing you to repeat the analysis or tweak as you wish"),
                     h2("Survival Analysis"),
                     helpText("You can perform Recursive Partitioning on a selected gene in a dataset with survival information (Cambridge, Stockholm and MSKCC). This analysis will determine if there are sub-groups of samples with significantly different expression level"),
                     helpText("If samples in the dataset can be allocated into different groups based on the expression of the gene, a Kaplan-Meir plot will be displayed. Otherwise, the median expression level of the gene will be used to assign samples to high and low expression groups"),
                     h2("Gene Correlation"),
                     helpText("You can plot one gene against another in a specified dataset. Points on the plot are coloured according to sample group. The correlation is computed and displayed. "),
                     h2("Heatmap"),
                     helpText("The uploaded gene list can be used to generate a heatmap from the chosen dataset. Control is given over the distance metric and clustering method. Samples can be partitioned into different groups based on the clustering, and the composition of each group can be interrogated"),
                     h2("Copy Number"),
                     helpText("For datasets with Copy number information (Cambridge, Stockholm and MSKCC), the frequency of alterations in different clinical covariates is displayed. A heatmap can also be generated"),
                     helpText("We are very grateful to Emilie Lalonde from University of Toronto for supplying the data for these plots"),
                     h2("Images"),
                     helpText("Spinning Wait Icons by Andrew Davidson http://andrewdavidson.com/articles/spinning-wait-icons/")

                  )
                 )
        ),
      



        tabPanel("Gene Profile",
                 sidebarLayout(
                   sidebarPanel(
                     conditionalPanel(
                       condition="($('html').hasClass('shiny-busy'))",
                       p("Shiny is computing something.."),
                       img(src="wait30trans.gif")
                     ),
  #                   selectInput("boxplotDataset","Choose a Dataset",choices=c("Cambridge","Stockholm","MSKCC", "Michigan2005","Michigan2012"),selected = "Cambridge"),
                     selectInput("clinvar_boxplot", "Choose a Clinical Covariate",choices=c("iCluster","Gleason","Sample_Group"),selected="iCluster"),
                     helpText("The covariates you can plot will be different for the various datasets"),
                     radioButtons("z_cambridge","Z-Score transform?",choices=c("Yes","No"),selected="Yes"),
                    helpText("The z-score transformation is recommended to put the expression values for each gene onto comparable scales"),
                     radioButtons("overlay_cambridge","Overlay individual points?",choices=c("Yes","No"),selected="Yes"),
                     h2("Gene List plotting options"),
                     helpText("You can choose whether to plot all genes in the gene list on the same plot"),
                     radioButtons("cambridgeCombPlot","Composite plot?",choices=c("Yes","No"),selected="Yes"),
                     helpText("If No is selected above, a particular gene from the list can be displayed"),
                     selectInput("cambridgeGeneChoice","Gene to plot", choices=c("STAT3","ESR1","AR","HES6","MELK"),selected="STAT3"),
                     h2("Output options"),
                     textInput("profileBasename", label = "What to call the output files",value="boxplot"),
                     selectInput("boxplotTheme", "Pick a plot style",choices=c("ggplot2","bw","classic","minimal","light","Wall Street Journal","Economist", "Excel","solarized","stata","calc","dark","fivethirtyeight","tufte"),selected="ggplot2"),
                    helpText("For more information on the different plot styles see the documentation for the", a("ggthemes", href="https://github.com/jrnold/ggthemes"),"package"),
                     radioButtons("profilePlotFormat", "File format for plots", choices=c("pdf","png"), selected="pdf"),
                     helpText("PDF can be imported into Illustrator (or similar) for editing. PNG plots are suitable for presentation"),
                     helpText("PDF dimensions are measured in inches, and PNG dimensions are measured in pixels"),
                     textInput("profileWidth", label="Width of plot ",value = 1200),
                     textInput("profileHeight", label="Height of plot ",value=600),
                     helpText("Plots and R scripts will have the extension pdf (/png) and R respectively"),
                     downloadButton("cambridgeBoxplotPDF","Export current profile(s)...."),
                     downloadButton("geneProfileScript","Download R script...."),
                     helpText("If you are using a gene list as input for the boxplots and have de-selected the composite plot option each gene will be plotted on a separate page")

                   ),
                   mainPanel(
   #                  helpText("Integration of copy number and transcriptomics provides risk stratification in prostate cancer: A discovery and validation cohort study"),
    #                 a("Ross-Adams et al. (2015) doi:10.1016/j.ebiom.2015.07.017",href="http://www.ebiomedicine.com/article/S2352-3964%2815%2930071-2/abstract"),
     #                helpText("These data were downloaded from GEO: GSE70768 and imported using the GEOquery Bioconductor package"),
      #               a("GSE70768",href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70768"),
       #              helpText("Also available as the "), a("prostateCancerCamcap",href="http://bioconductor.org/packages/devel/data/experiment/html/prostateCancerCamcap.html"), helpText(" Bioconductor package"),
                     h1("Gene Profile"),
                     verbatimTextOutput("profileMapping"),
#                     plotOutput("displayBoxplot",width = 1200,height=600),
                     imageOutput("displayBoxplot"),
                     h1("ANOVA analysis"),
                     helpText("Here we show the results of an ANOVA (analysis of variance) analysis to assess whether there are changes in expression level between the defined groups"),
                     DT::dataTableOutput("anovaResult")

                   )
                   
                 )
                 
        ),
    

    tabPanel("Survival",
             sidebarLayout(
                sidebarPanel(

                 conditionalPanel(
                    condition="($('html').hasClass('shiny-busy'))",
                    p("Shiny is computing something.."),
                    img(src="wait30trans.gif")
                  ),
#                  radioButtons("inputType_survival", "Use Single or Gene List as input?", choices=c("Single Gene","Gene List"),selected="Single Gene"),
#                  radioButtons("cutoffMethod","Type of partitioning?",choices=c("RP","Median","Manual"),selected="RP"),

                  helpText("You can select which gene to display the results for"),
                  selectInput("survivalGeneChoice","Gene to plot", choices=c("STAT3","ESR1","AR","HES6","MELK"),selected="STAT3"),
                  h2("Output options"),
                  textInput("survivalBasename", label = "What to call the output files",value="survival"),
                  radioButtons("survivalPlotFormat", "File format for plots", choices=c("pdf","png"), selected="pdf"),
                  helpText("PDF can be imported into Illustrator (or similar) for editing. PNG plots are suitable for presentation"),
                  helpText("PDF dimensions are measured in inches, and PNG dimensions are measured in pixels"),
                  textInput("survivalWidth", label="Width of plot",value=12),
                  textInput("survivalHeight", label="Height of plot",value=6.3),
                  downloadButton("survivalPlotPDF","Export K-M plot...."),
                  downloadButton("survivalScript","Download R script...."),
                  h2("Citation for recursive partitioning (RP)"),
                  em("[1] Torsten Hothorn, Kurt Hornik and Achim Zeileis (2006). Unbiased  Recursive Partitioning: A Conditional Inference Framework. Journal of  Computational and Graphical Statistics, 15(3), 651--674.")
                ),
                mainPanel(

                  verbatimTextOutput("survivalWarn"),
                  helpText("A recursive partitioning (RP) analysis [1] is performed to determine if the samples can be split into groups based on the expression data from your chosen gene(s). An RP p-value < 0.05 indicates a significant split. The p-value from RP and cut-off corresponding to a split are shown in the table below"),
                  helpText("If no cut-off can be found with RP, the samples will be divided according to median expression level in the plots below"),
                  DT::dataTableOutput("rpSummary"),
                  helpText("A histogram of expression level will be shown with a line to indicate the median expression level or RP cut-off"),
#                  plotOutput("rpPlot"),
                  imageOutput("rpPlot"),
                  helpText("The grouping of samples found by RP, or using median expression level, is used to construct a Kaplan-Meier plot"),
 #                 plotOutput("survivalPlot"),helpText("The Kaplan-Meier plot is a useful way of summarising survival data. There is one curve for each group. Each curve starts at 100% probability of survival. The probability of freedom from biochemical recurrence is shown on the y axis and the time (in years) is shown on the x axis. The curve drops each time there is an 'event'. A cross is shown on each curve where a 'censoring'' event takes place. This is where someone drops out of the study for a reason not related to the study, e.g. the study ends before an event has occurred. These subjects are no longer included in any calculations. The lower the survival curve the worse prognosis the patients in that group have.")
                  imageOutput("survivalPlot"),helpText("The Kaplan-Meier plot is a useful way of summarising survival data. There is one curve for each group. Each curve starts at 100% probability of survival. The probability of freedom from biochemical recurrence is shown on the y axis and the time (in years) is shown on the x axis. The curve drops each time there is an 'event'. A cross is shown on each curve where a 'censoring'' event takes place. This is where someone drops out of the study for a reason not related to the study, e.g. the study ends before an event has occurred. These subjects are no longer included in any calculations. The lower the survival curve the worse prognosis the patients in that group have.")
                )
                
             )
    ),
    
    tabPanel("Gene Correlation",
             sidebarLayout(
               sidebarPanel(
                conditionalPanel(
                   condition="($('html').hasClass('shiny-busy'))",
                   p("Shiny is computing something.."),
                   img(src="wait30trans.gif")
                 ),
                 radioButtons("inputType_correlation", "Plot correlations with a single gene, or all Pairwise Correlations?", choices=c("Single Gene","All Pairwise"),selected="Single Gene"),
                  radioButtons("corType", "Type of correlation to compute", selected="Pearson",choices=c("Pearson","Spearman")),
                  selectInput("correlationGeneChoice","Gene to plot", choices=c("STAT3","ESR1","AR","HES6","MELK"),selected="STAT3"),
#                 helpText("If you have selected Single Gene mode (above), now select a second gene to correlate with"),
#                 selectInput("secondGene","Type a Gene Symbol",choices="A2M",selected = "A2M"),
                 selectInput("clinvar_cor", "Choose a Clinical Covariate to colour by...",choices=c("None", "iCluster","Gleason","Sample_Group"),selected="None"),
                 h2("Output options"),
                 textInput("correlationBasename", label = "What to call the output files",value="correlation"),
                selectInput("corTheme", "Pick a plot style",choices=c("ggplot2","bw","classic","minimal","light","Wall Street Journal","Economist", "Excel","solarized","stata","calc","dark","fivethirtyeight","tufte"),selected="ggplot2"),
helpText("For more information on the different plot styles see the documentation for the", a("ggthemes", href="https://github.com/jrnold/ggthemes"),"package"),
                 radioButtons("correlationPlotFormat", "File format for plots", choices=c("pdf","png"), selected="pdf"),
                  helpText("PDF can be imported into Illustrator (or similar) for editing. PNG plots are suitable for presentation"),
                  helpText("PDF dimensions are measured in inches, and PNG dimensions are measured in pixels"),
                  textInput("correlationWidth", label="Width of plot",value=12),
                  textInput("correlationHeight", label="Height of plot",value=6.3),
                 downloadButton("CorrelationPDF","Export Correlation plot ...."),
                 downloadButton("correlationScript","Download R script....")
               ),
               mainPanel(
                plotOutput("corPlot",width = "100%")
                # imageOutput("corPlot")
               )
             )
    ),
    
    tabPanel("Heatmap",
            sidebarLayout(
             sidebarPanel(
              conditionalPanel(
                 condition="($('html').hasClass('shiny-busy'))",
                 p("Shiny is computing something.."),
                 img(src="wait30trans.gif")
               ),
                 radioButtons("distfun","Method to calulate distances",choices=c("Euclidean","Correlation"),selected="Euclidean"),
                 radioButtons("hclustfun", "Method of hierachical clustering",choices=c("ward","single","complete","average","mcquitty","median","centroid"),selected="complete"),
                 radioButtons("reordRows", "Re-order Rows?", choices = c("Yes","No"),selected = "Yes"),
                 radioButtons("scale","Scaling?",choices = c("row", "column", "none"),selected="row"),
                 radioButtons("cutType","Cut the dendrogram into k groups, or h(eight)",choices=c("k","h"),selected = "k"),
                 sliderInput("kGrps","Select k groups from the dendrogram",min=2,max = 10,value=2),
                 sliderInput("hCut","Select a height to cut the dendrogram",value = 3.8,min=1.8, max=3.8,step=0.1),
                 h2("Output options"),
                 textInput("heatmapBasename", label = "What to call the output files",value="example"),
                 radioButtons("heatmapPlotFormat", "File format for plots", choices=c("pdf","png"), selected="pdf"),
                 helpText("PDF can be imported into Illustrator (or similar) for editing. PNG plots are suitable for presentation"),
                 helpText("PDF dimensions are measured in inches, and PNG dimensions are measured in pixels"),
                 textInput("heatmapWidth", label="Width of plot",value=12),
                 textInput("heatmapHeight", label="Height of plot",value=12),
                 downloadButton("HeatmapPDF","Export heatmap..."),
                 downloadButton("heatmapScript","Download R script....")
               ),
               mainPanel(
                 helpText("Construcing a heatmap from the gene list you uploaded in the Analysis Parameters tab. If you haven't uploaded a gene list, an example gene list of three genes will be used"),
                 verbatimTextOutput("heatmapMapping"),
                 #plotOutput("heatmap",width = 1200,height=600),
                 imageOutput("heatmap"),
                 helpText("Sample Clustering"),
                 #plotOutput("dendrogram",width = 1200,height=600),
                 imageOutput("dendrogram"),
                 helpText("Select the number of clusters, k, from the slider"),
                 #plotOutput("sampleBreakown")
                 imageOutput("sampleBreakdown")
               )
               
             )
    ),
    tabPanel("Copy Number",
      sidebarLayout(
       sidebarPanel(
        conditionalPanel(
           condition="($('html').hasClass('shiny-busy'))",
           p("Shiny is computing something.."),
           img(src="wait30trans.gif")
         ),
#         radioButtons("inputType_cn", "Use Single or Gene List as input?", choices=c("Single Gene","Gene List"),selected="Single Gene"),                 h2("Output options"),
         textInput("copyNumberBasename", label = "What to call the output files",value="example"),
         selectInput("cnPlotType", "Type of plot to show", choices=c("Frequency", "Frequency by Dataset", "Heatmap"),selected="Frequency"),
         selectInput("clinvar_cn","Choose a Clinical Covariate to group by...",choices=c("iCluster","Gleason","Sample_Group"),selected="None"),
         helpText("Frequency: Overall Frequency of alterations in Cambridge, Stockholm and MSKCC"),
         helpText("Frequency in Dataset", "Frequency of alterations in a given covariate of interest in the chosen dataset"),
         helpText("Heatmap: Heatmap using the dataset that is currently selected"),
         selectInput("cnTheme", "Pick a plot style",choices=c("ggplot2","bw","classic","minimal","light","Wall Street Journal","Economist", "Excel","solarized","stata","calc","dark","fivethirtyeight","tufte"),selected="ggplot2"),
helpText("For more information on the different plot styles see the documentation for the", a("ggthemes", href="https://github.com/jrnold/ggthemes"),"package"),         
radioButtons("copyNumberPlotFormat", "File format for plots", choices=c("pdf","png"), selected="pdf"),
         helpText("PDF can be imported into Illustrator (or similar) for editing. PNG plots are suitable for presentation"),
         helpText("PDF dimensions are measured in inches, and PNG dimensions are measured in pixels"),
         textInput("copyNumberWidth", label="Width of plot",value=12),
         textInput("copyNumberHeight", label="Height of plot",value=6.3),
         downloadButton("copyNumberPDF","Export plot....")
        ),
        mainPanel(
          verbatimTextOutput("cnWarn"),
          helpText("The Proportion of amplifications and deletions will be shown for your chosen gene(s)."),
          plotOutput("copyNumber",width="100%"),
          #imageOutput("copyNumber"),
          DT::dataTableOutput("copyNumberTable")
        )
        
      )

    ),

tabPanel("Quick Analysis",
         sidebarLayout(
           sidebarPanel(
            conditionalPanel(
               condition="($('html').hasClass('shiny-busy'))",
               p("Shiny is computing something.."),
               img(src="wait30trans.gif")
             ),
             textInput("currentGene", "Gene of Interest", value="A1BG"),
             selectInput("quickDataset","Choose a Dataset",choices=c("Cambridge","Stockholm","MSKCC", "Michigan2005","Michigan2012"),selected = "Cambridge"),
             radioButtons("displayType", "Plot to display", choices=c("Boxplot","Survival","Copy Number"),selected="Boxplot"),
            actionButton("check", "Check Gene Name"),
            actionButton("go","Go!"),
             h2("Boxplot options"),
             selectInput("quick_clinvar_boxplot", "Choose a Clinical Covariate",choices=c("iCluster","Gleason","Sample_Group"),selected="iCluster"),
             helpText("The covariates you can plot will be different for the various datasets"),
             radioButtons("quick_z_cambridge","Z-Score transform?",choices=c("Yes","No"),selected="Yes"),
             radioButtons("quick_overlay_cambridge","Overlay individual points?",choices=c("Yes","No"),selected="Yes"),
             h2("Output options"),
             selectInput("quickTheme", "Pick a plot style",choices=c("ggplot2","bw","classic","minimal","light","Wall Street Journal","Economist", "Excel","solarized","stata","calc","dark","fivethirtyeight","tufte"),selected="ggplot2"),
            helpText("For more information on the different plot styles see the documentation for the", a("ggthemes", href="https://github.com/jrnold/ggthemes"),"package"),
             textInput("geneBasename", label = "What to call the output files",value="example"),
             radioButtons("genePlotFormat", "File format for plots", choices=c("pdf","png"), selected="pdf"),
             helpText("PDF can be imported into Illustrator (or similar) for editing. PNG plots are suitable for presentation"),
             helpText("PDF dimensions are measured in inches, and PNG dimensions are measured in pixels"),
             textInput("geneWidth", label="Width of plot",value=12),
             textInput("geneHeight", label="Height of plot",value=6.3),
            downloadButton("genePDF","Export Current plot...."),
             downloadButton("geneScript","Download R script....")
           ),
           mainPanel(
             verbatimTextOutput("geneCheck"),
             plotOutput("quickPlot",width="100%"),
             DT::dataTableOutput("quickTable")
           )
           
         )
         
)

                   
)

)
                   
