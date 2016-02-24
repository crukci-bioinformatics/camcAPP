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
                  plotOutput("rPlot")
                )
                
             )
    ),
    
    tabPanel("Heatamp",
            sidebarLayout(
             sidebarPanel(
                 selectInput("heatmapDataset","Choose a Dataset",choices=c("Cambridge","Stockholm","MSKCC"),selected = "MSKCC")
               ),
               mainPanel(
                 plotOutput("heatmap")
               )
               
             )
    )    
                   
)

)
                   
