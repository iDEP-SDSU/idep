# download public RNA-Seq data from ARCHS4
# needs to run GSEinfo.R script to generate GSE info file 
library(shiny)
library(DT) # for renderDataTable
library(shinyBS) # for popup figures

# Define UI ----
shinyUI( fluidPage(
  titlePanel("Public RNA-Seq and ChIP-Seq data")
         ,h5("Download gene-level read counts data for 7000+ human and mouse datasets from",
             a("ARCHS4 ", href="http://amp.pharm.mssm.edu/archs4/help.html"), 
             "on Nov. 5, 2018.", 
             "HiSeq 2000 and HiSeq 2500 raw data from NCBI's SRA database are 
             aligned with Kallisto against human GRCh38 or mouse GRCm38 cDNA reference.")
#  ,h5("Human datasets:",textOutput('humanNsamplesOutput'))
#  ,h5("Mouse datasets:",textOutput('mouseNsamplesOutput') )

         ,fluidRow( column(9, h5("Search and click on a dataset to see more information before download.") )
                    ,column(3, selectInput("selected.species.archs4", "", 
                     choices = list("Human"= "human", "Mouse" = "mouse"), selected = "human")     )
         )
         
         
         ,DT::dataTableOutput('SearchData')
         ,br()
         ,fluidRow( column(4, textOutput('selectedDataset'))
                    ,column(3, downloadButton('downloadSearchedData', 'Download'))
                    ,column(5, a("Return to iDEP", href="http://bioinformatics.sdstate.edu/idep/"))
                            
         )        

         ,tableOutput('samples' ) 
         
         

         
)
)