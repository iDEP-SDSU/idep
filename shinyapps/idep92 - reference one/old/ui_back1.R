
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(
  
#  fluidPage( 
  navbarPage("ITA",
             tabPanel("Load Data",
    
    titlePanel("Uploading Files"),
    sidebarLayout(
      sidebarPanel(
        fileInput('file1', 'Choose file to upload',
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    'text/tab-separated-values',
                    'text/plain',
                    '.csv',
                    '.tsv'
                  )
        ),
		actionButton("goButton", "Use example data"),
        tags$hr(),
		selectInput("selectOrg", label = NULL,"Best matching species",width='100%'), 
		tableOutput('species' )

      ),
      mainPanel(
        tableOutput('contents')
      )
    )
             ) #plots
			 

    ,tabPanel("Distribution",
             sidebarLayout(
               sidebarPanel(
			     numericInput("lowFilter", label = h5("Only keep genes with expression at least"), value = 1),
                 radioButtons("transform", "Log Transformation (recommended)",
                              c("Yes"=TRUE, "No"=FALSE)),
				numericInput("logStart", label = h5("Add stabilizing constant"), value = 1)
                 
               ),
               mainPanel(
                 plotOutput("EDA")
               )
             )       
    )
	        ,tabPanel("Processed Data",
             sidebarLayout(
               sidebarPanel(
             ),
               mainPanel(
                 tableOutput("debug")
               )
             )       
    )
    ,tabPanel("Heatmap",
              sidebarLayout(
                sidebarPanel(
   					sliderInput("nGenes", label = h4("Most variable genes to include"), min = 10, max = 2000, value = 200,step=20) 
					,downloadButton('downloadData', 'Download heatmap data')
					),
                mainPanel(
                  plotOutput("heatmap")
				 
                )
              )       
    )
     ,tabPanel("PCA",
              sidebarLayout(
                sidebarPanel(
 					sliderInput("nGenes1", label = h4("Most variable genes to include"), min = 40, max = 2000, value = 200,step=20) 
				),
                mainPanel(
                  plotOutput("PCA")
                )
              )       
    )
     ,tabPanel("Pathways",
              sidebarLayout(
                sidebarPanel(
				actionButton("goButton2", "Re-Analyze"),
				tags$hr(),
				htmlOutput("selectGO1"),
				
				sliderInput("minSetSize", label = h5("Min gene set size"), min = 5, max = 30, value = 15,step=1), 
			    sliderInput("maxSetSize", label = h5("Max gene set size"), min = 1000, max = 2000, value = 2000,step=100) ,
				numericInput("pathwayPvalCutoff", label = h5("Pathway signifiance cutoff (FDR)"), value = 0.1),
				numericInput("nPathwayShow", label = h5("Number of top pathway to show"), value = 30)
				
				),
                mainPanel(
                  plotOutput("PGSEAplot")
                )
              )       
    )	
  )# Navibar
 # ) # fluid page
  

)
