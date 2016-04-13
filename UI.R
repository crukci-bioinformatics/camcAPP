

library(shiny)
library(org.Hs.eg.db)

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
        tabPanel("Analysis Parameters",
                 sidebarLayout(
                   sidebarPanel(
                     fluidRow(
                       h1("Gene Selector"),
                     selectInput("currentGene","Type a Gene Symbol",choices=keys(revmap(org.Hs.egSYMBOL)),selected = "A1BG")

                     #fluidRow(
                       
#                      column(8,DT::dataTableOutput("geneTable"))
                      #textInput("currentGene", "Type a Gene Symbol",value="A1BG")
                     ),
                   
                     #),
                      fileInput('file1', 'Gene List',
                               accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),helpText("Your gene list must tab-delimited, with gene names in the first column"),
                      helpText("If no gene list is uploaded, the genes ESR1, AR and STAT3 will be used"),
radioButtons("inputType", "Use Single or Gene List as input?", choices=c("Single Gene","Gene List"),selected="Single Gene")
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
                     selectInput("boxplotDataset","Choose a Dataset",choices=c("Cambridge","Stockholm","MSKCC", "Michigan2005","Michigan2012"),selected = "Cambridge"),
                     selectInput("clinvar_boxplot", "Choose a Clinical Covariate",choices=c("iCluster","Gleason","Sample_Group"),selected="iCluster"),
                     helpText("The covariates you can plot will be different for the various datasets"),
                     radioButtons("z_cambridge","Z-Score transform?",choices=c("Yes","No"),selected="Yes"),
                     radioButtons("overlay_cambridge","Overlay individual points?",choices=c("Yes","No"),selected="Yes"),
                     h2("Gene List plotting options"),
                     helpText("You can choose whether to plot all genes in the gene list on the same plot"),
                     radioButtons("cambridgeCombPlot","Composite plot?",choices=c("Yes","No"),selected="Yes"),
                     helpText("Choosing a composite plot will display all genes in the gene list side-by-side. If No is selected, a particular gene from the list can be selected"),
                     selectInput("cambridgeGeneChoice","Gene to plot", choices=c("STAT3","ESR1","AR"),selected="STAT3"),
                     h2("Output options"),
                     textInput("profileBasename", label = "What to call the output files",value="sausage"),
                     helpText("Plots and R scripts will have the extension pdf and R respectively"),
                     downloadButton("cambridgeBoxplotPDF","Export current plot as pdf...."),
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
                     plotOutput("displayBoxplot",width = 1200,height=600),
                     h1("ANOVA analysis"),
                     verbatimTextOutput("anovaResult")

                   )
                   
                 )
                 
        ),
    

    tabPanel("Survival",
             sidebarLayout(
                sidebarPanel(
                  radioButtons("inputType_survival", "Use Single or Gene List as input?", choices=c("Single Gene","Gene List"),selected="Single Gene"),
                  selectInput("rpDataset","Choose a Dataset",choices=c("Cambridge","Stockholm","MSKCC"),selected = "MSKCC"),
                  radioButtons("cutoffMethod","Use Recursive Partitioning to choose a cut-off?",choices=c("RP","Median","Manual"),selected="RP"),
                  textInput("expCutOff", "Cut-off for partitioning",value = 6),
                  selectInput("survivalGeneChoice","Gene to plot", choices=c("STAT3","ESR1","AR"),selected="STAT3"),
                  textInput("survivalBasename", label = "What to call the output files",value="sausage"),
                  h2("Output options"),
                  downloadButton("survivalPlotPDF","Export current plot as pdf...."),
                  downloadButton("survivalProfileScript","Download R script....")
                ),
                mainPanel(
                  helpText("A recursive partitioning (RP) analysis is performed to determine if the samples can be split into groups based on the expression data from your chosen gene(s)."),
                  tableOutput("rpSummary"),
                  plotOutput("rpPlot"),
                  plotOutput("survivalPlot"),helpText("The Kaplan-Meier plot is a useful way of summarising survival data. There is one curve for each group. Each curve starts at 100% probability of survival. The probability of freedom from biochemical recurrence is shown on the y axis and the time (in years) is shown on the x axis. The curve drops each time there is an 'event'. A cross is shown on each curve where a 'censoring'' event takes place. This is where someone drops out of the study for a reason not related to the study, e.g. the study ends before an event has occurred. These subjects are no longer included in any calculations. The lower the survival curve the worse prognosis the patients in that group have.")
                )
                
             )
    ),
    
    tabPanel("Gene Correlation",
             sidebarLayout(
               sidebarPanel(
                 selectInput("secondGene","Type a Gene Symbol",choices=keys(revmap(org.Hs.egSYMBOL)),selected = "A2M"),
                 selectInput("corDataset","Choose a Dataset",choices=c("Cambridge","Stockholm","MSKCC"),selected = "Cambridge"),
                 radioButtons("corType","Type of correlation to calculate",choices=c("spearman","pearson"),selected="pearson")
               ),
               mainPanel(
                 plotOutput("corPlot")
               )
             )
    ),
    
    tabPanel("Heatmap",
            sidebarLayout(
             sidebarPanel(
                 selectInput("heatmapDataset","Choose a Dataset",choices=c("Cambridge","Stockholm","MSKCC"),selected = "MSKCC"),
                 radioButtons("distfun","Method to calulate distances",choices=c("Euclidean","Correlation"),selected="Euclidean"),
                 radioButtons("hclustfun", "Method of hierachical clustering",choices=c("ward","single","complete","average","mcquitty","median","centroid"),selected="complete"),
                 radioButtons("reordRows", "Re-order Rows?", choices = c("Yes","No"),selected = "Yes"),
                 radioButtons("scale","Scaling?",choices = c("row", "column", "none"),selected="row"),
                 radioButtons("cutType","Cut the dendrogram into k groups, or h(eight)",choices=c("k","h"),selected = "k"),
                 sliderInput("kGrps","Select k groups from the dendrogram",min=2,max = 7,value=2),
                 textInput("hCut","Select a height to cut the dendrogram",value = 10),
                 h2("Output options"),
                 textInput("heatmapBasename", label = "What to call the output files",value="sausage"),
                 downloadButton("HeatmapPDF","Export heatmap as pdf...."),
                 downloadButton("heatmapScript","Download R script....")
               ),
               mainPanel(
                 helpText("Construcing a heatmap from the gene list you uploaded in the Analysis Parameters tab. If you haven't uploaded a gene list, an example gene list of three genes will be used"),
                 plotOutput("heatmap",width = 1200,height=600),
                 helpText("Sample Clustering"),
                 plotOutput("dendrogram",width = 1200,height=600),
                 helpText("Select the number of clusters, k, from the slider"),
                 plotOutput("sampleBreakown")
               )
               
             )
    )    
                   
)

)
                   
