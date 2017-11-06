# download public RNA-Seq data from ARCHS4
# needs to run GSEinfo.R script to generate GSE info file 
library(shiny)
library(DT) # for renderDataTable
library(shinyBS) # for popup figures

# Define UI ----
shinyUI( fluidPage(
  titlePanel("Download processed public RNA-Seq and ChIP-Seq data")
  ,h5("Gene-level read counts data was downloaded from",a("ARCHS4 ", href="http://amp.pharm.mssm.edu/archs4/help.html"), "on Sept 1, 2017.")
  ,h5("HiSeq 2000 and HiSeq 2500 raw data from NCBI's SRA database are aligned with Kallisto against human GRCh38 or mouse GRCm38 cDNA reference.")
  ,h5("Human datasets:",a("HiSeq 2000,", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL11154" )," ", a("HiSeq 2500", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16791" ),textOutput('humanNsamplesOutput'))
  ,h5("Mouse datasets:",a("HiSeq 2000,", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13112" )," ", a("HiSeq 2500", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL17021" ),textOutput('mouseNsamplesOutput') )

,br()
,actionButton("SearchGSEB", "Search available datasets")
,br(),br()
  ,textInput("SearchGSE", h5("Enter GEO series number: "), value = "GSE53280") #GSE62734"
  ,tableOutput('samples' ) 
  ,downloadButton('downloadSearchedData', 'Download read counts matrix')                   
,bsModal("modalExample10", "Converted data (Most variable genes on top)", "SearchGSEB", size = "large", DT::dataTableOutput('SearchGSE'))
          ,br(),br()
		  ,h5("The downloaded read counts data can be analyzed with ", a("iDEP.", href="http://ge-lab.org:3838/idep/"), "The column names often need to be changed to define sample groups.")
		  ,br(),br(),h5( a("Contact us, ",href="mailto:xijin.ge@sdstate.edu?Subject=Downloader" )
         , "or visit our",a(" homepage.", href="http://ge-lab.org/") )
)
)