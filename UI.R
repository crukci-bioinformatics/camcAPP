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
        tabPanel("Cambridge",
                 sidebarLayout(
                   sidebarPanel(
                     selectInput("currentGene_cambridge","Type a Gene Symbol",choices=keys(revmap(org.Hs.egSYMBOL)),selected = "A1BG"),
                     selectInput("clinvar_cambridge", "Choose a Clinical Covariate",choices=c("iCluster","Gleason"),selected="iCluster"),
                     radioButtons("z_cambridge","Z-Score transform?",choices=c("Yes","No"),selected="Yes"),
                     radioButtons("overlay_cambridge","Overlay individual points?",choices=c("Yes","No"),selected="Yes")
                   ),
                   mainPanel(
                     plotOutput("boxplotCambridge"),verbatimTextOutput("anovaCambridge")
                   )
                   
                 )
                 
        ),
    
    tabPanel("Stockholm",
             sidebarLayout(
               sidebarPanel(
                 selectInput("currentGene_stockholm","Type a Gene Symbol",choices=keys(revmap(org.Hs.egSYMBOL)),selected = "A1BG"),
                 selectInput("clinvar_stockholm", "Choose a Clinical Covariate",choices=c("iCluster","Gleason"),selected="iCluster"),
                 radioButtons("z_stockholm","Z-Score transform?",choices=c("Yes","No"),selected="Yes"),
                 radioButtons("overlay_stockholm","Overlay individual points?",choices=c("Yes","No"),selected="Yes")
               ),
               mainPanel(
                 plotOutput("boxplotStockholm"),verbatimTextOutput("anovaStockholm")
               )
               
             )
             
    ),
    
    tabPanel("MSKCC",
             sidebarLayout(
               sidebarPanel(
                 selectInput("currentGene_taylor","Type a Gene Symbol",choices=keys(revmap(org.Hs.egSYMBOL)),selected = "A1BG"),
                 selectInput("clinvar_taylor", "Choose a Clinical Covariate",choices=c("CopyNumberCluster","Gleason"),selected="CopyNumberCluster"),
                 radioButtons("z_taylor","Z-Score transform?",choices=c("Yes","No"),selected="Yes"),
                 radioButtons("overlay_taylor","Overlay individual points?",choices=c("Yes","No"),selected="Yes")
               ),
               mainPanel(
                 plotOutput("boxplotTaylor"),verbatimTextOutput("anovaTaylor")
               )
               
             )
             
    )
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
) 
                   
                   
                   
                   
)