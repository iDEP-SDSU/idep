# download public RNA-Seq data from ARCHS4
# needs to run GSEinfo.R script to generate GSE info file
library(shiny)
library(DT) # for renderDataTable
library(shinyBS) # for popup figures
library(shinyjs)

# Define UI ----
shinyUI(
  fluidPage(
    fluidRow(
      column(12,
             offset = 0,
             titlePanel("iDEP-READS: Uniformly Processed Public RNA Sequencing Data"),
             h4(
               "Download counts data for 23,419 human and mouse datasets from",
               a("ARCHS4 v10", href = "http://amp.pharm.mssm.edu/archs4/help.html"),
               "and 29,662 datasets from",
               a("DEE2 ", href = "http://dee2.io/"),
               "for 9 model organisms. Click",
               actionLink("Statistics", "here"),
               "to see the number of datasets and samples by species.",
               br(),
               "To begin, select a species and a source below."
             ),
      ),
    ),
    
    #             "HiSeq 2000 and HiSeq 2500 raw data from NCBI's SRA database are
    #            aligned with Kallisto against human GRCh38 or mouse GRCm38 cDNA reference.")
    
    
    HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />'), # a solid line
    
    sidebarLayout(
      
      # Sidebar with a slider input
      sidebarPanel( width=2,
                    
                    radioButtons(
                      inputId = "selectedSpecies", label = "Select a Species", choices = "speciesChoice",
                      selected = NULL
                    ),
                    fluidRow(
                      column(12, offset = 0, h5("Upload 'counts' file to ", a("iDEP", href = "http://bioinformatics.sdstate.edu/idep/"), "for analysis")),
                      
                      column(12, offset = 0, 
                             downloadButton("downloadSearchedData", "Download Gene-Level Counts")),
                      column(2, downloadButton("downloadSearchedDataTranscript", "Download Transcript-Level Counts")),
                      
                    ),
      ),
      
      # Show a plot of the generated distribution
      
      mainPanel(width=10,
                DT::dataTableOutput("SearchData")
      )
    ),
    
    # GSEinfo table
    HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />'), # a solid line
    
    
    fluidRow(
      column(3, h4(textOutput("selectedDataset"))),
      
    ),
    conditionalPanel(
      "input.selectedSpecies != 'ARCHS4_Human' & input.selectedSpecies != 'ARCHS4_Mouse'", #& !is.null(input.selectedDataset)",
      fluidRow(column(3, offset = 9, textOutput("selected_selectedSpecies"), )),
      fluidRow(
        conditionalPanel(
          "input.selectedSpecies != 'ARCHS4_Human' & input.selectedSpecies != 'ARCHS4_Mouse'",
          column(8, ""),
          column(1, offset = 1, downloadButton("downloadSearchedDataGeneInfo", "Gene info")),
          column(1, downloadButton("downloadSearchedDataQCmat", "QC Matrix")),
          
        )
      ),
    ),
    fluidRow(
      column(
        5,
        tableOutput("samples")
        # ,h4("Loading data and R packages ... ...")
        # ,htmlOutput('DoneLoading')
      ),
      conditionalPanel(
        "input.selectedSpecies != 'ARCHS4_Human' & input.selectedSpecies != 'ARCHS4_Mouse'",
        column(1, offset = 4, downloadButton("downloadSearchedDataTxInfo", "Transcript info")),
        column(2, downloadButton("downloadSearchedDataSummaryMeta", "Summary MetaData"))
      )
    ),
    
    # table for dataset counts by species and source
    bsModal("modalExample1021", "Data Set Statistics and Details", "Statistics",
            size = "large",
            h5("ARCHS4 v6 data downloaded May 27, 2021. DEE2 metadata was downloaded August 4, 2021. DEE2 expression data is downloaded via API."),
            tableOutput("stats")
    ),
    tags$head(includeScript("ga.js")) # tracking usage
  )
)
