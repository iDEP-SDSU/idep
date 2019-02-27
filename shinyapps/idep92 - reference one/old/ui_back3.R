
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library("shinyAce") # for showing text files, code
library(shinyBS) # for popup figures


shinyUI(
  
#  fluidPage( 
  navbarPage("IEA",
             tabPanel("Load Data",
    
    titlePanel("Uploading Files"),
    sidebarLayout(
      sidebarPanel(
		 radioButtons("dataFileFormat", label = "Choose data type", choices = list("RNA-seq raw counts data" = 1, "Normalized expression values(FPKM, microarray etc)" = 2),selected = 2)

	  ,fileInput('file1', 'Choose file to upload',
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    'text/tab-separated-values',
                    'text/plain',
                    '.csv',
                    '.tsv'
                  )
        )

		,br(),br(),
		actionButton("goButton", "Use example data"),
        tags$hr(),
		selectInput("selectOrg", label = NULL,"Best matching species",width='100%'), 
		tableOutput('species' )
		

      ),
      mainPanel(
        tableOutput('contents')
		#dataTableOutput('contents')
		,htmlOutput('fileFormat')
      )
    )
             ) #plots
			 

    ,tabPanel("Pre-Process",
             sidebarLayout(
               sidebarPanel(
			     conditionalPanel("input.dataFileFormat == 2"
				 ,numericInput("lowFilter", label = h5("Only keep genes with expression at least"), value = 1)              
				 ,radioButtons("transform", "Log Transformation",c("No"=FALSE,"Yes"=TRUE) )
				,numericInput("logStart", label = h5("Constant c for started log: log(x+c)"), value = 1)
				,br(),h4( textOutput("text.transform") )
				)
				,conditionalPanel("input.dataFileFormat == 1",
				numericInput("minCounts", label = h5("Removing lowly expressed genes. Only keep those with this level ( in counts per million) in at least 2 samples"), value = 1) 
 				,radioButtons("CountsTransform", "Transform counts data for clustering & PCA.",  c("VST: variance stabilizing transform"=2, "rlog: regularized log (N<10) "= 3,"Started log: log2(x+c)" = 1),selected = 1 )
				,numericInput("countsLogStart", label = h5("Constant c for started log, if chosen above"), value = 4)
 				,radioButtons("CountsDEGMethod", "Detect differntial expression using:",  c("limma-voom"=2,"limma-trend"=1,"DESeq2"= 3) )
				)
                 
				 ,downloadButton('downloadProcessedData', 'Download processed data') 

               ),
               mainPanel(
                 plotOutput("EDA")
               )
             )       
    )
	#        ,tabPanel("Processed Data",
     #        sidebarLayout(
     #          sidebarPanel(
     #        ),
     #          mainPanel(
     #            tableOutput("debug")
     #          )
      #       )       
    #)
    ,tabPanel("Heatmap",
              sidebarLayout(
                sidebarPanel(
   					sliderInput("nGenes", label = h4("Most variable genes to include"), min = 10, max = 4000, value = 1000,step=100) 
					,downloadButton('downloadData', 'Download heatmap data')
					),
                mainPanel(
                  plotOutput("heatmap")
				 
                )
              )       
    )
	,tabPanel("k-Means",
              sidebarLayout(
                sidebarPanel(
				#numericInput("nClusters", label = h4("Number of Clusters (often <15) "), value = 6)
				sliderInput("nClusters", label = h4("Number of Clusters"), min = 3, max = 20, value = 6,step=1) 
   				,sliderInput("nGenesKNN", label = h4("Most variable genes to include "), min = 10, max = 6000, value = 1000,step=20) 
				
				,actionButton("showMotifKmeans", "Promoter analysis of each cluster")
				,br(),br(),downloadButton('downloadDataKNN', 'Download K-means data')
				,h5("Database for enrichment analysis")
				,htmlOutput("selectGO3"),tags$style(type='text/css', "#selectGO3 { width:100%;   margin-top:-9px}")

				),
                mainPanel(
                  plotOutput("KNN")
				  ,br(),br(),br(),br()
				  ,h4("Enriched pathways for each cluster")
				  ,tableOutput("KNN_GO")
				  ,bsModal("modalExample2", "Enriched TF binding motifs in promoters of Kmeans clusters", "showMotifKmeans", size = "large"
				   ,radioButtons("radioPromoterKmeans", label = NULL, choices = list("Upstream 300bp as promoter" = 300, "Upstream 600bp as promoter" = 600),selected = 300)
				   ,tableOutput("KmeansPromoter"))
                )
              )       
    )
     ,tabPanel("PCA",
              sidebarLayout(
                sidebarPanel(
				radioButtons("PCA_MDS", "Methods", c("Principal Component Analysis"=1, "Pathway Analysis of PCA rotation" =2, 
				"Multidimensional Scaling"=3))
				,downloadButton('downloadPCAData', 'Download Coordinates')
 					#sliderInput("nGenes1", label = h4("Most variable genes to include"), min = 40, max = 2000, value = 200,step=20) 
				),
                mainPanel(
                  plotOutput("PCA")
                )
              )       
    )
	     ,tabPanel("DEGs",
              sidebarLayout(
                sidebarPanel(
				h5("Identifying Differential Expressed Genes (DEGs) using  limma.")
				,numericInput("limmaPval", label = h5("FDR cutoff"), value = 0.05)
				,tags$style(type='text/css', "#limmaPval { width:100%;   margin-top:-12px}")
				,numericInput("limmaFC", label = h5("Min. fold change"), value = 2)
				,tags$style(type='text/css', "#limmaFC { width:100%;   margin-top:-12px}")
				,h4( textOutput("text.limma") )
				,actionButton("showVenn", "Venn Diagram of all lists")
				 ,br(),br(),downloadButton('download.DEG.data', 'Download all gene lists')
				 ,hr()
				 ,htmlOutput("listComparisons")
				 ,h5("Database for enrichment analysis")
				 ,htmlOutput("selectGO2")
				,tags$style(type='text/css', "#selectGO2 { width:100%;   margin-top:-9px}")
				 ,downloadButton('download.selectedHeatmap.data', "Download heatmap data" )
				 ,br(),br(),actionButton("showMotif", "TF binding motifs in promoter")
				 #,radioButtons("radio.promoter", label = NULL, choices = list("Upstream 300bp as promoter" = 300, "Upstream 600bp as promoter" = 600),selected = 300)

				,h4("Top Genes"),tableOutput('geneList')
				, width = 4),
				
                mainPanel(
                  plotOutput("selectedHeatmap")
				   ,h4("Enriched GO terms in up (A) or down (B) -regulated genes")
				   ,tableOutput("geneListGO")
				   
				   #,h4("Enriched motif in promoters")
				   #,tableOutput("DEG.Promoter")
				   ,bsModal("modalExample1", "Enriched TF binding motifs in promoters of DEGs", "showMotif", size = "large"
				   ,radioButtons("radio.promoter", label = NULL, choices = list("Upstream 300bp as promoter" = 300, "Upstream 600bp as promoter" = 600),selected = 300)
				   ,tableOutput("DEG.Promoter"))
				   ,bsModal("modalExample", "Venn Diagram", "showVenn", size = "large",plotOutput("vennPlot"))

				   
				  
                )
              )       
    )
   ,tabPanel("Pathways",
              sidebarLayout(
                sidebarPanel(
				actionButton("goButton2", "Run/Re-Run")
				,tags$hr()
				,htmlOutput("listComparisonsPathway")
				,htmlOutput("selectGO1")
				,sliderInput("minSetSize", label = h5("Min gene set size"), min = 5, max = 30, value = 15,step=1) 
			    ,sliderInput("maxSetSize", label = h5("Max gene set size"), min = 1000, max = 2000, value = 2000,step=100) 
				,numericInput("pathwayPvalCutoff", label = h5("Pathway signifiance cutoff (FDR)"), value = 0.1)
				,numericInput("nPathwayShow", label = h5("Number of top pathway to show"), value = 30)
				,downloadButton('download.PGSEAplot.data', 'Download pathway data')
				),
                mainPanel(
                  plotOutput("PGSEAplot")
                )
              )       
    )

	
#   ,tabPanel("Code",
#      fluidRow(
#       column(12,
#          "R codes for reference",
#            sourceCode1 <- aceEditor("code"
#  , value = paste(readLines("ui.R"), collapse="\n")
#  , mode = "r"
#  , theme = "ambience"
#  , height = "400px"
#  , readOnly = TRUE
#  ) )))
  
  )# Navibar
 # ) # fluid page
  

)
