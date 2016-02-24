####shiny::runGitHub("OneSidedT","markdunning")

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
                     selectInput("currentGene","Type a Gene Symbol",choices=keys(revmap(org.Hs.egSYMBOL)),selected = "A1BG"),
                     fileInput('file1', 'Gene List',
                               accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),helpText("Your gene list must tab-delimited, with gene names in the first column")
                   ),
                   mainPanel(
                     h3("About this app.............")
                  )
                 )
        ),

        tabPanel("Cambridge Profile",
                 sidebarLayout(
                   sidebarPanel(

                     selectInput("clinvar_cambridge", "Choose a Clinical Covariate",choices=c("iCluster","Gleason"),selected="iCluster"),
                     radioButtons("z_cambridge","Z-Score transform?",choices=c("Yes","No"),selected="Yes"),
                     radioButtons("overlay_cambridge","Overlay individual points?",choices=c("Yes","No"),selected="Yes")
                   ),
                   mainPanel(
                     plotOutput("boxplotCambridge"),verbatimTextOutput("anovaCambridge")
                   )
                   
                 )
                 
        ),
    
      tabPanel("Stockholm Profile",
               sidebarLayout(
                 sidebarPanel(
                   selectInput("clinvar_stockholm", "Choose a Clinical Covariate",choices=c("iCluster","Gleason"),selected="iCluster"),
                   radioButtons("z_stockholm","Z-Score transform?",choices=c("Yes","No"),selected="Yes"),
                   radioButtons("overlay_stockholm","Overlay individual points?",choices=c("Yes","No"),selected="Yes")
                 ),
                 mainPanel(
                   plotOutput("boxplotStockholm"),verbatimTextOutput("anovaStockholm")
                 )
                 
               )
               
      ),
    
      tabPanel("MSKCC Profile",
               sidebarLayout(
                 sidebarPanel(
                   selectInput("clinvar_taylor", "Choose a Clinical Covariate",choices=c("CopyNumberCluster","Gleason"),selected="CopyNumberCluster"),
                   radioButtons("z_taylor","Z-Score transform?",choices=c("Yes","No"),selected="Yes"),
                   radioButtons("overlay_taylor","Overlay individual points?",choices=c("Yes","No"),selected="Yes")
                 ),
                 mainPanel(
                   plotOutput("boxplotTaylor"),verbatimTextOutput("anovaTaylor")
                 )
                 
               )
               
      ),
    
    tabPanel("Survival",
             sidebarLayout(
                sidebarPanel(
                  selectInput("rpDataset","Choose a Dataset",choices=c("Cambridge","Stockholm","MSKCC"),selected = "MSKCC")
                ),
                mainPanel(
                  plotOutput("rpPlot"),helpText("The plot above is a recursive partitioning plot for the first probe of the BIRC5 gene"),
                  plotOutput("survivalPlot"),helpText("The Kaplan-Meier plot is a useful way of summarising survival data. There is one curve for each group. Each curve starts at 100% probability of survival. The probability of freedom from biochemical recurrence is shown on the y axis and the time (in years) is shown on the x axis. The curve drops each time there is an 'event'. A cross is shown on each curve where a 'censoring'' event takes place. This is where someone drops out of the study for a reason not related to the study, e.g. the study ends before an event has occurred. These subjects are no longer included in any calculations. The lower the survival curve the worse prognosis the patients in that group have.")
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
                 radioButtons("scale","Scaling?",choices = c("row", "column", "none"),selected="row")
               ),
               mainPanel(
                 plotOutput("heatmap")
               )
               
             )
    )    
                   
)

)
                   
