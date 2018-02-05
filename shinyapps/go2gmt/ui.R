# This app converts gene-GO id mappings to a GMT file. 
library(shiny)
   

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Convert Gene Ontology files to GMT file"),
  
  sidebarLayout(
    sidebarPanel(
      h5("This app converts gene id to GO id mappings into a GMT file containing one GO id per line followed by all gene ids")  
      , fileInput('files', 'Upload gene id to Gene Ontology mapping file'
                  
      )
      
      ),
    
    mainPanel(
      textOutput("stat")
	  ,br(),br()
	  ,downloadButton('download.File',"Download combined file")
      
    )
  )
  )


