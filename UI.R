

library(shiny)


shinyUI(navbarPage("Explore Prostate Cancer Datasets", id = "nav",
                  
    #    tabPanel("MSKCC",
     #            sidebarLayout(
      #             sidebarPanel(
       #              selectInput("currentGene","Type a Gene Symbol",choices=keys(revmap(org.Hs.egSYMBOL)),selected = "")
        #           ),
         #          mainPanel(
                  #   plotOutput("boxplotMSKCC"),verbatimTextOutput("anovaMSKCC")
          #         )
                   
           #      )
                 
            #     ),
        tabPanel("Data Input",
                 sidebarLayout(
                   sidebarPanel(
                     
                   
                     #),
  #                   radioButtons("inputType", "Use Single or Gene List as input?", choices=c("Single Gene","Gene List"),selected="Single Gene"),
                      fileInput('file1', 'Gene List',
                               accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),helpText("Your gene list must tab-delimited, with gene names in the first column"),
                      helpText("If no gene list is uploaded, the genes ESR1, AR and STAT3 will be used"),
   #                   textInput("currentGene", "Gene of Interest", value="A1BG"),
                      helpText("If you want to analyse a single-gene, see the Quick Analysis tab"),
                      selectInput("theDataset","Choose a Dataset",choices=c("Cambridge","Stockholm","MSKCC", "Michigan2005","Michigan2012"),selected = "Cambridge")
                   )
                  ,
                   mainPanel(
                     h3("About this app............."),
                     img(src="cruk-cambridge-institute.jpg",width=350,height=77), br(),a("cruk.cam.ac.uk",href="www.cruk.cam.ac.uk"),
                     helpText("This app was developed by"),a("Cancer Research Uk Cambridge Institute Bioinformatics Core",href="http://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core"),
                     helpText("Source code available on github;"),code("https://github.com/crukci-bioinformatics/camcAPP"),
                     h3("Features"),
                     h2("Profile a specified gene across different datasets"),
                     helpText("Choose a gene from the 'Type a Gene Symbol' box and the first five tabs (Cambridge Profile ... Michigan-2012) will show the expression of this gene in different datasets. For each dataset, you can choose which clinical variable to group the samples on"),
                     helpText("When choosing Cambridge or Stockholm, you will have the option to display the expression in the five different subtypes identified by Ross-Adams et al (2015). These subtypes were shown to have significantly different outcomes"),
                     img(src="km-from-paper.png"),
                     helpText("If multiple microarray probes are found for the gene, the probe with the highest inter-quartile range (IQR) will be picked"),
                     helpText("An ANOVA analysis will also be performed to assess whether there are different expression levels in the groups you have chosen"),
                     helpText("You can also upload a gene list and the boxplots will be displayed for each gene individually [CURRENTLY CAMBRIDGE ONLY]"),
                     helpText("An R script can be downloaded, allowing you to repeat the analysis or tweak as you wish"),
                     h2("Survival Analysis"),
                     helpText("You can perform Recursive Partitioning on a selected gene in a dataset with survival information (Cambridge, Stockholm and MSKCC)"),
                     helpText("If samples in the dataset can be allocated into different groups based on the expression of the gene, a Kaplan-Meir plot will be displayed"),
                     helpText("Otherwise, no plots will be displayed"),
                     h2("Gene Correlation"),
                     helpText("You can plot one gene against another in a specified dataset. Points on the plot are coloured according to sample group. The correlation is computed and displayed. "),
                     h2("Heatmap"),
                     helpText("The uploaded gene list can be used to generate a heatmap from the chosen dataset. Control is given over the distance metric and clustering method")
                     
                     
                  )
                 )
        ),
      



        tabPanel("Gene Profile",
                 sidebarLayout(
                   sidebarPanel(
                     conditionalPanel(
                       condition="($('html').hasClass('shiny-busy'))",
                       p("Shiny is computing something.."),
                       img(src="http://i.imgur.com/tbIq2nD.gif")
                     ),
  #                   selectInput("boxplotDataset","Choose a Dataset",choices=c("Cambridge","Stockholm","MSKCC", "Michigan2005","Michigan2012"),selected = "Cambridge"),
                     selectInput("clinvar_boxplot", "Choose a Clinical Covariate",choices=c("iCluster","Gleason","Sample_Group"),selected="iCluster"),
                     helpText("The covariates you can plot will be different for the various datasets"),
                     radioButtons("z_cambridge","Z-Score transform?",choices=c("Yes","No"),selected="Yes"),
                     radioButtons("overlay_cambridge","Overlay individual points?",choices=c("Yes","No"),selected="Yes"),
                     h2("Gene List plotting options"),
                     helpText("You can choose whether to plot all genes in the gene list on the same plot"),
                     radioButtons("cambridgeCombPlot","Composite plot?",choices=c("Yes","No"),selected="Yes"),
                     helpText("If No is selected above, a particular gene from the list can be displayed"),
                     selectInput("cambridgeGeneChoice","Gene to plot", choices=c("STAT3","ESR1","AR"),selected="STAT3"),
                     h2("Output options"),
                     textInput("profileBasename", label = "What to call the output files",value="boxplot"),
                     selectInput("boxplotTheme", "Pick a plot style",choices=c("ggplot2","bw","classic","minimal","light"),selected="ggplot2"),
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
                     plotOutput("displayBoxplot",width = 1200,height=600),
                     h1("ANOVA analysis"),
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
                    img(src="http://i.imgur.com/tbIq2nD.gif")
                  ),
#                  radioButtons("inputType_survival", "Use Single or Gene List as input?", choices=c("Single Gene","Gene List"),selected="Single Gene"),
                  radioButtons("cutoffMethod","Type of partitioning?",choices=c("RP","Median","Manual"),selected="RP"),
                  textInput("expCutOff", "Cut-off for partitioning",value = 6),
                  helpText("If you working on a gene list, you can select which gene to display the results for"),
                  selectInput("survivalGeneChoice","Gene to plot", choices=c("STAT3","ESR1","AR"),selected="STAT3"),
                  h2("Output options"),
                  textInput("survivalBasename", label = "What to call the output files",value="survival"),
                  radioButtons("survivalPlotFormat", "File format for plots", choices=c("pdf","png"), selected="pdf"),
                  helpText("PDF can be imported into Illustrator (or similar) for editing. PNG plots are suitable for presentation"),
                  helpText("PDF dimensions are measured in inches, and PNG dimensions are measured in pixels"),
                  textInput("survivalWidth", label="Width of plot",value=12),
                  textInput("survivalHeight", label="Height of plot",value=6.3),
                  downloadButton("survivalPlotPDF","Export K-M plot as pdf...."),
                  downloadButton("survivalScript","Download R script....")
                ),
                mainPanel(
                  verbatimTextOutput("survivalMapping"),
                  helpText("A recursive partitioning (RP) analysis is performed to determine if the samples can be split into groups based on the expression data from your chosen gene(s)."),
                  DT::dataTableOutput("rpSummary"),
                  plotOutput("rpPlot"),
                  plotOutput("survivalPlot"),helpText("The Kaplan-Meier plot is a useful way of summarising survival data. There is one curve for each group. Each curve starts at 100% probability of survival. The probability of freedom from biochemical recurrence is shown on the y axis and the time (in years) is shown on the x axis. The curve drops each time there is an 'event'. A cross is shown on each curve where a 'censoring'' event takes place. This is where someone drops out of the study for a reason not related to the study, e.g. the study ends before an event has occurred. These subjects are no longer included in any calculations. The lower the survival curve the worse prognosis the patients in that group have.")
                )
                
             )
    ),
    
    tabPanel("Gene Correlation",
             sidebarLayout(
               sidebarPanel(
                conditionalPanel(
                   condition="($('html').hasClass('shiny-busy'))",
                   p("Shiny is computing something.."),
                   img(src="http://i.imgur.com/tbIq2nD.gif")
                 ),
                 radioButtons("inputType_correlation", "Plot correlations with a single gene, or all Pairwise Correlations?", choices=c("Single Gene","All Pairwise"),selected="Single Gene"),
                  selectInput("correlationGeneChoice","Gene to plot", choices=c("STAT3","ESR1","AR"),selected="STAT3"),
#                 helpText("If you have selected Single Gene mode (above), now select a second gene to correlate with"),
#                 selectInput("secondGene","Type a Gene Symbol",choices="A2M",selected = "A2M"),
                 selectInput("clinvar_cor", "Choose a Clinical Covariate to colour by...",choices=c("None", "iCluster","Gleason","Sample_Group"),selected="None"),
                 h2("Output options"),
                 textInput("correlationBasename", label = "What to call the output files",value="correlation"),
                 radioButtons("correlationPlotFormat", "File format for plots", choices=c("pdf","png"), selected="pdf"),
                  helpText("PDF can be imported into Illustrator (or similar) for editing. PNG plots are suitable for presentation"),
                  helpText("PDF dimensions are measured in inches, and PNG dimensions are measured in pixels"),
                  textInput("correlationWidth", label="Width of plot",value=12),
                  textInput("correlationHeight", label="Height of plot",value=6.3),
                 downloadButton("CorrelationPDF","Export Correlation plot as pdf...."),
                 downloadButton("correlationScript","Download R script....")
               ),
               mainPanel(
                 plotOutput("corPlot",width = 1200,height=600)
               )
             )
    ),
    
    tabPanel("Heatmap",
            sidebarLayout(
             sidebarPanel(
              conditionalPanel(
                 condition="($('html').hasClass('shiny-busy'))",
                 p("Shiny is computing something.."),
                 img(src="http://i.imgur.com/tbIq2nD.gif")
               ),
                 radioButtons("distfun","Method to calulate distances",choices=c("Euclidean","Correlation"),selected="Euclidean"),
                 radioButtons("hclustfun", "Method of hierachical clustering",choices=c("ward","single","complete","average","mcquitty","median","centroid"),selected="complete"),
                 radioButtons("reordRows", "Re-order Rows?", choices = c("Yes","No"),selected = "Yes"),
                 radioButtons("scale","Scaling?",choices = c("row", "column", "none"),selected="row"),
                 radioButtons("cutType","Cut the dendrogram into k groups, or h(eight)",choices=c("k","h"),selected = "k"),
                 sliderInput("kGrps","Select k groups from the dendrogram",min=2,max = 7,value=2),
                 textInput("hCut","Select a height to cut the dendrogram",value = 10),
                 h2("Output options"),
                 textInput("heatmapBasename", label = "What to call the output files",value="example"),
                 radioButtons("heatmapPlotFormat", "File format for plots", choices=c("pdf","png"), selected="pdf"),
                 helpText("PDF can be imported into Illustrator (or similar) for editing. PNG plots are suitable for presentation"),
                 helpText("PDF dimensions are measured in inches, and PNG dimensions are measured in pixels"),
                 textInput("heatmapWidth", label="Width of plot",value=12),
                 textInput("heatmapHeight", label="Height of plot",value=12),
                 downloadButton("HeatmapPDF","Export heatmap as pdf...."),
                 downloadButton("heatmapScript","Download R script....")
               ),
               mainPanel(
                 helpText("Construcing a heatmap from the gene list you uploaded in the Analysis Parameters tab. If you haven't uploaded a gene list, an example gene list of three genes will be used"),
                 verbatimTextOutput("heatmapMapping"),
                 plotOutput("heatmap",width = 1200,height=600),
                 helpText("Sample Clustering"),
                 plotOutput("dendrogram",width = 1200,height=600),
                 helpText("Select the number of clusters, k, from the slider"),
                 plotOutput("sampleBreakown")
               )
               
             )
    ),
    tabPanel("Copy Number",
      sidebarLayout(
       sidebarPanel(
        conditionalPanel(
           condition="($('html').hasClass('shiny-busy'))",
           p("Shiny is computing something.."),
           img(src="http://i.imgur.com/tbIq2nD.gif")
         ),
#         radioButtons("inputType_cn", "Use Single or Gene List as input?", choices=c("Single Gene","Gene List"),selected="Single Gene"),                 h2("Output options"),
         textInput("copyNumberBasename", label = "What to call the output files",value="example"),
         selectInput("cnPlotType", "Type of plot of show", choices=c("Frequency", "Heatmap"),selected="Frequency"),
         selectInput("cnTheme", "Pick a plot style",choices=c("ggplot2","bw","classic","minimal","light"),selected="ggplot2"),
         radioButtons("copyNumberPlotFormat", "File format for plots", choices=c("pdf","png"), selected="pdf"),
         helpText("PDF can be imported into Illustrator (or similar) for editing. PNG plots are suitable for presentation"),
         helpText("PDF dimensions are measured in inches, and PNG dimensions are measured in pixels"),
         textInput("copyNumberWidth", label="Width of plot",value=12),
         textInput("copyNumberHeight", label="Height of plot",value=6.3),
         downloadButton("copyNumberPDF","Export plot...."),
         downloadButton("copyNumberScript","Download R script....")
        ),
        mainPanel(
          helpText("The Proportion of amplifications and deletions will be shown for your chosen gene(s)."),
          plotOutput("copyNumber",width=1200,height=600),
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
               img(src="http://i.imgur.com/tbIq2nD.gif")
             ),
             textInput("currentGene", "Gene of Interest", value="A1BG"),
             selectInput("quickDataset","Choose a Dataset",choices=c("Cambridge","Stockholm","MSKCC", "Michigan2005","Michigan2012"),selected = "Cambridge"),
             radioButtons("displayType", "Plot to display", choices=c("Boxplot","Survival","Copy Number"),selected="Boxplot"),
            actionButton("go","Go!"),
             h2("Boxplot options"),
             selectInput("quick_clinvar_boxplot", "Choose a Clinical Covariate",choices=c("iCluster","Gleason","Sample_Group"),selected="iCluster"),
             helpText("The covariates you can plot will be different for the various datasets"),
             radioButtons("quick_z_cambridge","Z-Score transform?",choices=c("Yes","No"),selected="Yes"),
             radioButtons("quick_overlay_cambridge","Overlay individual points?",choices=c("Yes","No"),selected="Yes"),
             h2("Output options"),
             selectInput("quickTheme", "Pick a plot style",choices=c("ggplot2","bw","classic","minimal","light"),selected="ggplot2"),
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

             plotOutput("quickPlot",width=1200,height=600),
             DT::dataTableOutput("quickTable")
           )
           
         )
         
)

                   
)

)
                   
