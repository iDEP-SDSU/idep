
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
  navbarPage("iDEP 0.14",
             tabPanel("Load Data",
	h5("Integrated Differential Expression and Pathway analysis (iDEP) of genomics data.  See ",
	a(" documentation", href="https://idepsite.wordpress.com/"), "and",
a(" Demo.", href="https://ge600.wordpress.com/2017/02/28/using-idep-to-analyze-rna-seq-counts-data-v2/"),
   
 "Based on annotation of 69 metazoa and 42 plant genomes in Ensembl BioMart as of 11/15/2016."
            ," Additional  data from"
			,a("KEGG, ", href="www.genome.jp/kegg/")
			,a("Reactome, ", href="http://www.reactome.org/")
			,a("MSigDB (human),", href="https://doi.org/10.1093/bioinformatics/btr260") 
         ,a("GSKB (mouse)", href="http://biorxiv.org/content/early/2016/10/24/082511") 
         ,"and",a("  araPath (arabidopsis).", href="https://doi.org/10.1093/bioinformatics/bts421") 
         ," For feedbacks or data contributions (genes and GO mapping of any species), please"
         ,a("contact us, ",href="mailto:xijin.ge@sdstate.edu?Subject=ShinyGO" )
         , "or visit our",a(" homepage.", href="http://ge-lab.org/")
		 , "Send us suggestions or any error message to help us improve iDEP."
         ),  
    titlePanel("Uploading Files"),
          
    sidebarLayout(
      sidebarPanel(
		radioButtons("dataFileFormat", label = "Choose data type", choices = list("RNA-seq raw counts data (recommended)" = 1, "Normalized expression values (RNA-seq FPKM, microarray, etc.)" = 2),selected = 1)

		,fileInput('file1', 'Choose file to upload (<50Mb)',
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    'text/tab-separated-values',
                    'text/plain',
                    '.csv',
                    '.tsv'
					
					
                  )
        ),
		actionButton("goButton", "Use demo data"),
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
				numericInput("minCounts", label = h5("Removing genes with less than this many reads in all samples: "), value = 10) 
 				,radioButtons("CountsTransform", "Transform counts data for clustering & PCA.",  c("VST: variance stabilizing transform"=2, "rlog: regularized log (only for N<10) "= 3,"Started log: log2(x+c)" = 1),selected = 1 )
				,conditionalPanel("input.CountsTransform == 1",numericInput("countsLogStart", label = h5("Constant c for started log"), value = 4))
 				#,radioButtons("CountsDEGMethod", "Detect differntial expression using:",  c("limma-voom"=2,"limma-trend"=1,"DESeq2"= 3) )
				)
                 
				 ,downloadButton('downloadProcessedData', 'Download processed data') 

               ),
               mainPanel(
			   conditionalPanel("input.dataFileFormat == 1", plotOutput("totalCounts") )
                  , plotOutput("EDA")
               )
             )       
    )

    ,tabPanel("Heatmap",
              sidebarLayout(
                sidebarPanel(
   					sliderInput("nGenes", label = h4("Most variable genes to include in hierarchical clustering"), min = 10, max = 4000, value = 1000,step=50) 
					,downloadButton('downloadData', 'Download heatmap data')
					,br(),br()
					,actionButton("showCorrelation", "Show correlation matrix")

					),
                mainPanel(				  
  				  plotOutput("heatmap")
				  ,bsModal("modalExample8", "Correlation matrix using all genes", "showCorrelation", size = "large",plotOutput("correlationMatrix"))
				 
                )
              )       
    )
	,tabPanel("k-Means",
              sidebarLayout(
                sidebarPanel(
				#numericInput("nClusters", label = h4("Number of Clusters (often <15) "), value = 6)
				sliderInput("nClusters", label = h4("Number of Clusters"), min = 2, max = 20, value = 4,step=1) 
   				,sliderInput("nGenesKNN", label = h4("Most variable genes to include "), min = 10, max = 6000, value = 2000,step=100) 
				
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
				h5("Identifying Differential Expressed Genes (DEGs)")
				,conditionalPanel("input.dataFileFormat == 1",
			    #  radioButtons("CountsDEGMethod", "Method to detect DEGs:",  c("limma-voom"=2,"limma-trend"=1,"DESeq2"= 3) )
				selectInput("CountsDEGMethod", "Method to detect DEGs from read counts:", choices = list("DESeq2"= 3,"limma-voom"=2,"limma-trend"=1), selected = 2)				
				)				
				,fluidRow(
				 column(5,numericInput("limmaPval", label = h5("FDR cutoff"), value = 0.1,min=1e-20,max=1,step=.05)  )
				 ,column(7, numericInput("limmaFC", label = h5("Min fold change"), value = 2,min=1,max=100,step=0.5) )
				) # fluidRow
				
				
				,tags$style(type='text/css', "#limmaPval { width:100%;   margin-top:-12px}")
				,tags$style(type='text/css', "#limmaFC { width:100%;   margin-top:-12px}")
				,fluidRow(
				   column(6,actionButton("showVenn", "Venn Diagram") )
					, column(6, downloadButton('download.DEG.data', 'All lists') )
				) # fluidRow
				 #,hr()
				 ,HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />') # a solid line
				 ,htmlOutput("listComparisons")
				 ,h5("Genesets for enrichment analysis")
				 ,htmlOutput("selectGO2")
				,tags$style(type='text/css', "#selectGO2 { width:100%;   margin-top:-9px}")
				 ,downloadButton('download.selectedHeatmap.data', "Download heatmap data" )
				 ,br(),br(),fluidRow(
				 column(6,actionButton("showVolcano", "Volcano Plot"))
				 ,column(6, actionButton("showScatter", "Scatter Plot") ) )
				 ,br(),actionButton("showMotif", "TF binding motifs in promoter")
				 #,radioButtons("radio.promoter", label = NULL, choices = list("Upstream 300bp as promoter" = 300, "Upstream 600bp as promoter" = 600),selected = 300)
				,br(),br(),h4( textOutput("text.limma") )
				
				, width = 4),
				
                mainPanel(
                  plotOutput("selectedHeatmap")
				   ,h4("Enriched gene sets in up (A) or down (B) -regulated genes")
				   ,tableOutput("geneListGO")
				   ,h4("Top Genes"),tableOutput('geneList')
				   #,h4("Enriched motif in promoters")
				   #,tableOutput("DEG.Promoter")
				   ,bsModal("modalExample1", "Enriched TF binding motifs in promoters of DEGs", "showMotif", size = "large"
				   ,radioButtons("radio.promoter", label = NULL, choices = list("Upstream 300bp as promoter" = 300, "Upstream 600bp as promoter" = 600),selected = 300)
				   ,tableOutput("DEG.Promoter"))
				   ,bsModal("modalExample", "Venn Diagram", "showVenn", size = "large",plotOutput("vennPlot"))
				   ,bsModal("modalExample4", "Volcano plot", "showVolcano", size = "large",plotOutput("volcanoPlot"))
				   ,bsModal("modalExample5", "Scatter plot", "showScatter", size = "large",plotOutput("scatterPlot"))
				  
                )
              )       
    )
   ,tabPanel("Pathways",
              sidebarLayout(
                sidebarPanel(
				htmlOutput("listComparisonsPathway")
				,tags$style(type='text/css', "#listComparisonsPathway { width:100%;   margin-top:-12px}")
			,selectInput("pathwayMethod", label = "Select method and genesets:", choices = list("GAGE" = 1,  "GSEA (preranked fgsea)" =3,"PGSEA" = 2, "PGSEA w/ all samples" =4, "ReactomePA" = 5), selected = 1) #
				,tags$style(type='text/css', "#pathwayMethod { width:100%;   margin-top:-12px}")
				,htmlOutput("selectGO1")
				,tags$style(type='text/css', "#selectGO1 { width:100%;   margin-top:-12px}")
                ,h5("Genesets size filter:")
				,fluidRow( 
				column(6,numericInput("minSetSize", label = "Min", min = 5, max = 30, value = 15,step=1) 
			           ),
				column(6,numericInput("maxSetSize", label = "Max", min = 1000, max = 2000, value = 2000,step=100) 
						)
				) # fluidRow
				,numericInput("pathwayPvalCutoff", label = h5("Pathway signifiance cutoff (FDR)"), value = 0.1,min=1e-20,max=1,step=.05)
				,tags$style(type='text/css', "#pathwayPvalCutoff { width:100%;   margin-top:-12px}")
				,numericInput("nPathwayShow", label = h5("Number of top pathways to show"), value = 30, min=5,max=100,step=5)
				,tags$style(type='text/css', "#nPathwayShow { width:100%;   margin-top:-12px}")
				,checkboxInput("absoluteFold", label = "Use absolute values of fold changes for GSEA and GAGE", value = FALSE)
				#,actionButton("examinePathway", "Examine individual pathways")
				,conditionalPanel("input.pathwayMethod == 2",downloadButton('download.PGSEAplot.data', 'Download PGSEA pathway data'))
				),
                mainPanel(	  
					
				  conditionalPanel("input.pathwayMethod == 2",
                  plotOutput("PGSEAplot") )
				  ,conditionalPanel("input.pathwayMethod == 1",
                  tableOutput("gagePathway") )
				  ,conditionalPanel("input.pathwayMethod == 3",
                  tableOutput("fgseaPathway") )
				   ,conditionalPanel("input.pathwayMethod == 5",
                  tableOutput("ReactomePAPathway") )
				  ,conditionalPanel("input.pathwayMethod == 4",
				  plotOutput("PGSEAplotAllSamples") )
				   ,conditionalPanel("input.pathwayMethod == 1 | input.pathwayMethod == 3" 
				    ,htmlOutput("listSigPathways")
					,downloadButton('downloadSelectedPathwayData', 'Download data on selected pathway')  )
					,conditionalPanel(" (input.pathwayMethod == 1 | input.pathwayMethod == 3 ) & input.selectGO == 'KEGG'",imageOutput("KeggImage"))
					,conditionalPanel("(input.pathwayMethod == 1 | input.pathwayMethod == 3 ) & input.selectGO != 'KEGG'",plotOutput("selectedPathwayHeatmap"))
				  
			#	  ,bsModal("modalExample3", "Examine Pathways", "examinePathway", size = "large"
			#	    ,htmlOutput("listSigPathways")
			#		,downloadButton('downloadSelectedPathwayData', 'Download pathway data')
			#		,conditionalPanel("input.selectGO == 'KEGG'"
					
			#		,imageOutput("KeggImage"))
			#		,plotOutput("selectedPathwayHeatmap")
					#,tableOutput("selectedPathwayData1")
			#		)
				   
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
   ,tags$head(includeScript("google_analytics.js")) # tracking usage  
  )# Navibar
 # ) # fluid page
   

)
