# download public RNA-Seq data from ARCHS4
# needs to run GSEinfo.R script to generate GSE info file 
library(shiny)
library(DT) # for renderDataTable
library(shinyBS) # for popup figures

# Define UI ----
shinyUI( fluidPage(
  titlePanel("iDEP-reads: Uniformlly processed public RNA-Seq data")
         ,h5("Read counts data for 5,470 human and mouse datasets from",
             a("ARCHS4 v6", href="http://amp.pharm.mssm.edu/archs4/help.html"), 
             "and 12,670 datasets from",
             a("DEE2 ", href="http://dee2.io/"), "for 9 model organisms.",actionButton("Statistics", "Dataset statistics"))         
#             "HiSeq 2000 and HiSeq 2500 raw data from NCBI's SRA database are 
 #            aligned with Kallisto against human GRCh38 or mouse GRCm38 cDNA reference.")
         ,fluidRow( column(9, h5("Select your species, search,  and click on a dataset to see more information before download.") )
                    #,column(3, selectInput("selectedSpecies", "", 
                    # choices = list("Human_ARCHS4"= "Human_ARCHS4", "Mouse_ARCHS4" = "Mouse_ARCHS4"), selected = "Human_ARCHS4"))
                    ,column(3, selectInput("selectedSpecies", label = NULL,"Human_ARCHS4",width='100%') )
         )

         ,DT::dataTableOutput('SearchData')
         ,HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />') # a solid line    
         ,br()
         ,fluidRow( column(2, h4( textOutput('selectedDataset')))
                    ,column(3, downloadButton('downloadSearchedData', 'Gene-level counts'))
                   ,column(3,downloadButton('downloadSearchedDataTranscript', 'Transcript-level counts') )                   
                    ,column(4, h5("Edit and uplodad file to ",a("iDEP", href="http://bioinformatics.sdstate.edu/idep/"), "for analysis") )                            
         )

      ,conditionalPanel("input.selectedSpecies != 'ARCHS4_Human' & input.selectedSpecies != 'ARCHS4_Mouse'",
                        fluidRow(
                          column(8, "" )
                          ,column(2, downloadButton('downloadSearchedDataGeneInfo', 'Gene info') )
                          ,column(2, downloadButton('downloadSearchedDataTxInfo', 'Transcript info'))
                        )
        )
                         
         ,br()
         ,tableOutput('samples' ) 
        # ,h4("Loading data and R packages ... ...")
        # ,htmlOutput('DoneLoading')
        ,bsModal("modalExample1021", "Data set statistics and details", "Statistics", size = "large",
          h5( "ARCHS4 v6 data downloaded June 4, 2019. DEE2 metadata was downloaded June 4, 2019. DEE2 expression data is downloaded via API.") 
          ,tableOutput('stats') 
        )
        ,tags$head(includeScript("ga.js")) # tracking usage  
)

)