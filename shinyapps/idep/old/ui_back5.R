
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
  navbarPage("iDEP",
             tabPanel("Load Data",
	h5("Integrated Differential Expression and Pathway analysis (iDEP) of genomics data.  ",
	a(" See demo and ask questions here. ", href="https://ge600.wordpress.com/2017/02/28/using-idep-to-analyze-rna-seq-counts-data-v2/")),
         h5("Based on annotation of 69 metazoa and 42 plant genomes in Ensembl BioMart as of 11/15/2016."
            ," Additional  data from",a("MSigDB (human),", href="https://doi.org/10.1093/bioinformatics/btr260") 
         ,a("GSKB (mouse)", href="http://biorxiv.org/content/early/2016/10/24/082511") 
         ,"and",a("  araPath (arabidopsis).", href="https://doi.org/10.1093/bioinformatics/bts421") 
         , "Built with R and", a("Shiny!",href="http://shiny.rstudio.com/")
         ," For feedbacks or data contributions (genes and GO mapping of any species), please"
         ,a("contact us, ",href="mailto:xijin.ge@sdstate.edu?Subject=ShinyGO" )
         , "or visit our",a(" homepage.", href="http://ge-lab.org/")
         ),    
    titlePanel("Uploading Files"),
    sidebarLayout(
      sidebarPanel(
		 radioButtons("dataFileFormat", label = "Choose data type", choices = list("RNA-seq raw counts data (recommended)" = 1, "Normalized expression values (RNA-seq FPKM, microarray, etc.)" = 2),selected = 1)

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

		
		,actionButton("goButton", "Use example data"),
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
				numericInput("minCounts", label = h5("Removing genes with small counts in all samples:"), value = 10) 
 				,radioButtons("CountsTransform", "Transform counts data for clustering & PCA.",  c("VST: variance stabilizing transform"=2, "rlog: regularized log (N<10) "= 3,"Started log: log2(x+c)" = 1),selected = 1 )
				,numericInput("countsLogStart", label = h5("Constant c for started log, if chosen above"), value = 4)
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
				selectInput("CountsDEGMethod", "Method to detect DEGs:", choices = list("DESeq2"= 3,"limma-voom"=2,"limma-trend"=1), selected = 1)
				
				)				
				
				
				,numericInput("limmaPval", label = h5("FDR cutoff"), value = 0.1)
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

				   
				  
                )
              )       
    )
   ,tabPanel("Pathways",
              sidebarLayout(
                sidebarPanel(
				actionButton("goButton2", "Run/Re-Run")
				,tags$hr()
				,htmlOutput("listComparisonsPathway")
				,selectInput("pathwayMethod", label = h5("Method and geneset database for pathway analysis"), choices = list("GAGE" = 1,  "GSEA (preranked fgsea)" =3,"PGSEA" = 2, "PGSEA w/ all samples" =4, "ReactomePA" = 5), selected = 1) #
				,htmlOutput("selectGO1")
				,sliderInput("minSetSize", label = h5("Min gene set size"), min = 5, max = 30, value = 15,step=1) 
			    ,sliderInput("maxSetSize", label = h5("Max gene set size"), min = 1000, max = 2000, value = 2000,step=100) 
				,numericInput("pathwayPvalCutoff", label = h5("Pathway signifiance cutoff (FDR)"), value = 0.2)
				,numericInput("nPathwayShow", label = h5("Number of top pathway to show"), value = 30)
				,checkboxInput("absoluteFold", label = h5("Use absolute values of fold changes in GSEA and GAGE analysis." ), value = FALSE)
				,actionButton("examinePathway", "Examine pathways")
				
				,downloadButton('download.PGSEAplot.data', 'Download PGSEA pathway data')
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
				 # ,htmlOutput("listSigPathways")
				 #	,downloadButton('downloadSelectedPathwayData', 'Download data for selected pathway')
				 #	,plotOutput("selectedPathwayHeatmap")
					#,tableOutput("selectedPathwayData1") 
				  ,bsModal("modalExample3", "Examine Pathways", "examinePathway", size = "large"
				    ,htmlOutput("listSigPathways")
					,downloadButton('downloadSelectedPathwayData', 'Download pathway data')
					,plotOutput("selectedPathwayHeatmap")
					#,tableOutput("selectedPathwayData1")
					)
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
