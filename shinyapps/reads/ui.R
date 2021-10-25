# download public RNA-Seq data from ARCHS4
# needs to run GSEinfo.R script to generate GSE info file
#library(shiny)
#library(DT) # for renderDataTable
#library(shinyBS) # for popup figures
#library(shinyjs)

# Define UI ----
shiny::shinyUI(
  shiny::fluidPage(
    shinyjs::useShinyjs(),
    shiny::fluidRow(
      shiny::column(12,
             offset = 0,
             shiny::titlePanel("iDEP-READS: Uniformly Processed Public RNA Sequencing Data"),
             h4(
               "Download counts data for 23,419 human and mouse datasets from",
               a("ARCHS4 v10", href = "http://amp.pharm.mssm.edu/archs4/help.shiny::HTML"),
               "and 29,662 datasets from",
               a("DEE2 ", href = "http://dee2.io/"),
               "for 9 model organisms. Click",
               actionLink("Statistics", "here"),
               "to see the number of datasets and samples by species.",
               shiny::br(),
               "To begin, select a species and a source below. Then, select a row on the table."
             ),
      ),
    ),
    
    #             "HiSeq 2000 and HiSeq 2500 raw data from NCBI's SRA database are
    #            aligned with Kallisto against human GRCh38 or mouse GRCm38 cDNA reference.")
    
    
    shiny::HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />'), # a solid line
    
   # shinybusy::add_busy_spinner(spin = "breeding-rhombus", position="top-left", timeout=200, color="#8A2BE2", margin= c(50,50)),

   
   
    shiny::sidebarLayout(
      
      # Sidebar with a slider input
      shiny::sidebarPanel( width=2,
                    
                    shiny::radioButtons(
                      inputId = "selectedSpecies", label = "Select a Species", choices = "speciesChoice",
                      selected = NULL
                    ),
                    shiny::fluidRow(
                      shiny::column(12, offset = 0, h5("Upload 'counts' file to ", a("iDEP", href = "http://bioinformatics.sdstate.edu/idep/"), "for analysis")),
                      
                      shiny::column(12, offset = 0, 
                             downloadButton("downloadSearchedData", "Download Gene-Level Counts")),
                      shiny::column(1, downloadButton("downloadSearchedDataTranscript", "Download Transcript-Level Counts")),
                      
                    ),
                    #bottom of sidebar
                    shiny::fluidRow(
                      shiny::conditionalPanel("input.selectedSpecies != 'ARCHS4_Human' & input.selectedSpecies != 'ARCHS4_Mouse'",

                      shiny::column(12,shiny::br(), textOutput("selected_selectedSpecies")),
                      
                      shiny::column(12, downloadButton("downloadSearchedDataGeneInfo", "Gene info")),
                      shiny::column(12, downloadButton("downloadSearchedDataTxInfo", "Transcript info")),
                      shiny::column(12, downloadButton("downloadSearchedDataSummaryMeta", "Summary MetaData")),
                      shiny::column(11, downloadButton("downloadSearchedDataQCmat", "QC Matrix"))
                      
                      )
                    )
      ),
      
      # Show a plot of the generated distribution
      
      shiny::mainPanel(width=10,
                DT::dataTableOutput("SearchData"),
                
                shiny::column(6, shiny::tableOutput("samples1"))

    ),
    

   ),
    # GSEinfo table
    shiny::HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />'), # a solid line
    
    
    shiny::fluidRow(
      shiny::column(3, h4(textOutput("selectedDataset"))),
      
    ),
    shiny::conditionalPanel(
      "input.selectedSpecies != 'ARCHS4_Human' & input.selectedSpecies != 'ARCHS4_Mouse' & input.GSEID != 'null'",
      #shiny::fluidRow(shiny::column(3, offset = 9, textOutput("selected_selectedSpecies"), )),
      shiny::fluidRow(
        shiny::conditionalPanel(
          "input.selectedSpecies != 'ARCHS4_Human' & input.selectedSpecies != 'ARCHS4_Mouse'& input.GSEID != 'null'",
          shiny::column(8, ""),
       #   shiny::column(1, offset = 1, downloadButton("downloadSearchedDataGeneInfo", "Gene info")),
         # shiny::column(1, downloadButton("downloadSearchedDataQCmat", "QC Matrix")),
          
        )
      ),
    ),
    shiny::fluidRow(
      shiny::column(
        5,
        shiny::tableOutput("samples")
        # ,h4("Loading data and R packages ... ...")
        # ,shiny::HTMLOutput('DoneLoading')
      ),
      shiny::conditionalPanel(
        "input.selectedSpecies != 'ARCHS4_Human' & input.selectedSpecies != 'ARCHS4_Mouse' & input.GSEID != 'null'",
        #shiny::column(2, downloadButton("downloadSearchedDataSummaryMeta", "Summary MetaData"))
      )
    ),
    
    # table for dataset counts by species and source
    shinyBS::bsModal("modalExample1021", "Data Set Statistics and Details", "Statistics",
            size = "large",
            h5("ARCHS4 v6 data downloaded May 27, 2021. DEE2 metadata was downloaded August 4, 2021. DEE2 expression data is downloaded via API."),
            shiny::tableOutput("stats")
    ),
    tags$head(includeScript("ga.js")) # tracking usage
  )
)
