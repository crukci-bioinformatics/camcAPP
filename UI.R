####shiny::runGitHub("OneSidedT","markdunning")

library(shiny)

shinyUI(pageWithSidebar(
  
  headerPanel("Explore the Taylor Prostate Cancer dataset...."),
  
  sidebarPanel(
    h2("Data Import Parameters"),
    fileInput('file1', 'Gene List',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),helpText("Your gene list must tab-delimited, with gene names in the first column")
    ),
  
  mainPanel(
    tabsetPanel(
#      tabPanel("Plot",plotOutput("plot")),
      
      tabPanel("The data", h2("The gene list provided was.."),verbatimTextOutput("geneList"), h2("Details of the probe mapping"),verbatimTextOutput("mapping")),
      
#      tabPanel("Boxplot",plotOutput("boxplot")),
#      tabPanel("Histogram",plotOutput("histogram")),

#      tabPanel("Summary Statistics",
#               h4("Screen output in R"),
#               verbatimTextOutput("summary")),
tabPanel("Data Distribution", helpText("The boxplot and histogram of the data are shown below"),
         plotOutput("boxplot"),
         verbatimTextOutput("summary")
         ),

tabPanel("Heatmap", helpText("A heatmap showing the clustering of the samples according to your selected genes"),
         plotOutput("heatmap")
),

tabPanel("Recursive Partitioning", helpText("Recursive Partitioning....."),

         verbatimTextOutput("rp"),
         plotOutput("rp_plot")
)
,

      tabPanel("R code",
               helpText("You will be able to re-run this analysis in R by downloading the R code below"),
               strong("The input file that you are analysing must be in your R working directory in order for the script to run"),
               h4("Code Preview"),
               verbatimTextOutput("thecode"),
               downloadLink('downloadScript', 'Download R Script'),
               br(),
               downloadLink('downloadMarkdown', 'Download R Markdown'),
               br(),
  #             downloadLink('downloadPDF', 'Download HTML Report')
  helpText("We recommend RStudio to run the R code and compile reports"),
  img(src="https://www.rstudio.com/wp-content/uploads/2014/03/blue-125.png"), br(),a("RStudio",href="https://www.rstudio.com/"),br(),
  helpText("In order to compile the report in RStudio, you will need to install the ggplot2, rmarkdown and knitr packages"),br(),
  code("install.packages(c('knitr','ggplot2','rmarkdown'))")
              )
    )
  )
  
  ))