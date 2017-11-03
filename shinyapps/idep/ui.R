# iDEP user interface
library(shiny)
library("shinyAce") # for showing text files, code
library(shinyBS) # for popup figures
library(plotly)
iDEPversion = "iDEP.44"
# 0.38 Gene ID conversion, remove redudancy;  rlog option set to blind=TRUE
# 0.39 reorganized code. Updated to Bioconductor 3.5; solved problems with PREDA 9/8/17
# 0.40 moved libraries from the beginning to different places to save loading time
   # change colors for heatmap 
   # fix error with voom; also changed method to TMM, see https://www.bioconductor.org/help/workflows/RNAseq123/#transformations-from-the-raw-scale
   
shinyUI(
  
  navbarPage(iDEPversion,
			id='navBar',
            tabPanel("Load Data",
	     	titlePanel("Upload Files"),
    sidebarLayout(
      sidebarPanel(
	    actionButton("goButton", "Click here to load demo data"),
		    tags$head(tags$style("#goButton{color: red;
                                 font-size: 16px;
                                 font-style: italic;
                                 }"
                         )
              ),	
        h5(" and click all the tabs for results!",  style = "color:red"),	
		br(),
		a(h5("Reset all",align = "right"), href="http://ge-lab.org/idep/"), 
		radioButtons("dataFileFormat", label = "1. Choose data type", choices = list("Read counts data (recommended)" = 1, "Normalized expression values (RNA-seq FPKM, microarray, etc.)" = 2),selected = 1)

		,fileInput('file1', '2. Upload expression data (CSV or text)',
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    'text/tab-separated-values',
                    'text/plain',
                    '.csv',
                    '.tsv'					
                  )
        )		
		,fileInput('file2', h5('Optional: Upload sample info file (CSV or text)'),
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    'text/tab-separated-values',
                    'text/plain',
                    '.csv',
                    '.tsv'					
                  )
        )
		
		# tags$hr(),
		#actionButton("goButton2", "Change Species"),
		,strong("3. Verify guessed species. Change if neccessary."),
		selectInput("selectOrg", label = NULL,"Best matching species",width='100%'), 
				
		 conditionalPanel("input.selectOrg == 'NEW'",
			fileInput('gmtFile', 'Upload a geneset .GMT file for enrichment analysis (optional)',
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    'text/tab-separated-values',
                    'text/plain',
                    '.csv',
                    '.tsv',
					'.gmt'					
                  )
        )
		 ),
		#strong("4. Click the tabs for results!"),
tableOutput('species' ),
		a(h5("?",align = "right"), href="https://idepsite.wordpress.com/data-format/",target="_blank")
                                                                                       # new window
      ),
      mainPanel(
      tableOutput('sampleInfoTable')
     ,tableOutput('contents')
		
		,h5("Integrated Differential Expression and Pathway analysis (iDEP) of transcriptomic data.  See ",
			a(" documentation", href="https://idepsite.wordpress.com/"), "and",
			a(" manuscript.", href="http://biorxiv.org/content/biorxiv/early/2017/06/09/148411.full.pdf"),
   			"Based on annotation of",a( " 69 metazoa and 42 plant genomes ",href="https://idepsite.wordpress.com/species/") ,"in Ensembl BioMart as of 6/4/2017."
            ," Additional  data from"
			,a("KEGG, ", href="www.genome.jp/kegg/")
			,a("Reactome, ", href="http://www.reactome.org/")
			,a("MSigDB (human),", href="https://doi.org/10.1093/bioinformatics/btr260") 
			,a("GSKB (mouse)", href="http://biorxiv.org/content/early/2016/10/24/082511") 
			,"and",a("  araPath (arabidopsis).", href="https://doi.org/10.1093/bioinformatics/bts421") 
			," For feedbacks or data contributions (genes and GO mapping of any species), please"
			,a("contact us, ",href="mailto:xijin.ge@sdstate.edu?Subject=iDEP" )
			, "or visit our",a(" homepage.", href="http://ge-lab.org/")
			, "Send us suggestions or any error message to help improve iDEP."
         )
		,h3("Loading R packages ... ...")
		,htmlOutput('fileFormat')
      )
    )
             ) #plots
			 
###############################################################################################################################
 
    ,tabPanel("Pre-Process",
             sidebarLayout(
               sidebarPanel(
			     conditionalPanel("input.dataFileFormat == 2"
					,numericInput("lowFilter", label = h5("Only keep genes with expression at least:"), value = 1)
					,radioButtons("transform", "Log Transformation",c("No"=FALSE,"Yes"=TRUE) )
					,numericInput("logStart", label = h5("Constant c for started log: log(x+c)"), value = 1)
					,br(),h4( textOutput("text.transform") )
				)
				,conditionalPanel("input.dataFileFormat == 1",
					 
					strong("Keep genes with minimal counts per million (CPM) in at least n libraries:")
					,fluidRow(
						column(6, numericInput("minCounts", label = h5("Min. CPM"), value = 0.5) )
						,column(6, numericInput("NminSamples", label = h5("n libraries"), value = 1) )
					) # fluidRow
					,tags$style(type='text/css', "#minCounts { width:100%;   margin-top:-12px}")
					,tags$style(type='text/css', "#NminSamples { width:100%;   margin-top:-12px}")
							
					,radioButtons("CountsTransform", "Transform counts data for clustering & PCA.",  c("VST: variance stabilizing transform"=2, 
					"rlog: regularized log (slow) "= 3,"edgeR's logCPM with prior count" = 1),selected = 1 )
					,conditionalPanel("input.CountsTransform == 1",
						fluidRow(
							column(5, h5("Prior count:")  )
							,column(7, numericInput("countsLogStart", label = NULL, value = 1) )
						)
					        
					
					)
 				#,radioButtons("CountsDEGMethod", "Detect differntial expression using:",  c("limma-voom"=2,"limma-trend"=1,"DESeq2"= 3) )
				
				)
				,actionButton("genePlot1", "Barplot for one or more genes")
				,br(),br()
				,actionButton("examineDataB", "Search processed data")
                 ,br(),br()
				 ,downloadButton('downloadProcessedData', 'Download processed data') 
				 ,br(),br()
				 ,textOutput('nGenesFilter')
				,tags$head(tags$style("#nGenesFilter{color: blue;
											 font-size: 16px;
											 font-style: italic;
											 }"
									 )
					)
			,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/pre-process/",target="_blank")
               ),
               mainPanel(
			   conditionalPanel("input.dataFileFormat == 1", plotOutput("totalCounts") )
                  , plotOutput("EDA")				  
				  ,bsModal("modalExample10", "Converted data (Most variable genes on top)", "examineDataB", size = "large", DT::dataTableOutput('examineData'))
				  ,bsModal("modalExample1021", "Search for genes", "genePlot1", size = "large", 
					textInput("geneSearch", "Enter full or partial gene ID:", "HOXA"),
					checkboxInput("genePlotBox", label = "Show individual samples", value = FALSE),
					plotOutput("genePlot"),
					conditionalPanel("input.genePlotBox == 0", checkboxInput("useSD", label = "Use standard deviation instead of standard error", value = FALSE))
					)
	
               )
             )       
    )

	
	
###############################################################################################################################
 
    ,tabPanel("Heatmap",
              sidebarLayout(
                sidebarPanel(
   					sliderInput("nGenes", label = h4("Most variable genes to include:"), min = 10, max = 6000, value = 1000,step=50) 
					,actionButton("showStaticHeatmap", "Interactive heatmap")
					,br()
					,actionButton("showCorrelation", "Correlation matrix")			
					,HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />') # a solid line
					,strong("Customize hierarchical clustering (Default values work well):")
					,fluidRow(
						column(3, h5("Color")  )
						,column(9, selectInput("heatColors1", label = NULL,"green-black-red",width='100%') )
					)
					,fluidRow(
						column(4, h5("Distance")  )
						,column(8, selectInput("distFunctions", label = NULL,"Correlation",width='100%') )
					)					
					,fluidRow(
						column(4, h5("Linkage")  )
						,column(8, selectInput("hclustFunctions", label = NULL,"average",width='100%') )
					)
					,fluidRow(
						column(8, h5("Cut-off Z score")  )
						,column(4, numericInput("heatmapCutoff", label = NULL, value = 4,min=2,step=1) )
					)					
					,checkboxInput("geneNormalize", "Normalize by gene", value = FALSE)
					,checkboxInput("sampleNormalize", "Normalize by sample", value = FALSE)
					,checkboxInput("noSampleClustering", "Do not re-order or cluster samples", value = FALSE)
					,htmlOutput('listFactorsHeatmap')

					,downloadButton('downloadData', 'Download heatmap data')
					,br(),a(h5("?",align = "right"), href="https://idepsite.wordpress.com/heatmap/",target="_blank")
						),
					mainPanel(				  
						plotOutput("heatmap1")

						# ,verbatimTextOutput("event")
						,bsModal("modalExample8", "Correlation matrix using all genes", "showCorrelation", size = "large",plotOutput("correlationMatrix"))
						,bsModal("modalExample28", "Heatmap with hierarchical clustering tree", "showStaticHeatmap", size = "large",
						sliderInput("nGenesPlotly", label = h4("Most variable genes to include:"), min = 10, max = 6000, value = 50,step=50),
						h4("Mouse over to see gene names. To zoom, click and drag up or downward and release."),
						plotlyOutput("heatmap",width = "100%", height = "800px"))
					 
					)
					)       
					)

###############################################################################################################################
 
	,tabPanel("k-Means",
              sidebarLayout(
                sidebarPanel(
				#numericInput("nClusters", label = h4("Number of Clusters (often <15) "), value = 6)
   				sliderInput("nGenesKNN", label = h4("Most variable genes to include "), min = 10, max = 6000, value = 2000,step=100) 
				,sliderInput("nClusters", label = h4("Number of Clusters"), min = 2, max = 20, value = 4,step=1) 
				,actionButton("NClusters", "How many clusters?")
				,br(),br(),actionButton("showMotifKmeans", "Promoter analysis of each cluster")
				,br(),br(),downloadButton('downloadDataKmeans', 'Download K-means data')
				,h5("Pathway database")
				,htmlOutput("selectGO3"),tags$style(type='text/css', "#selectGO3 { width:100%;   margin-top:-9px}")
				
						 
			,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/k-means/",target="_blank")
				),
                mainPanel(
                  plotOutput("KmeansHeatmap")
				  ,br(),br(),br(),br()
				  ,h4("Enriched pathways for each cluster")
				  ,tableOutput("KmeansGO")
				  ,bsModal("modalExample2", "Enriched TF binding motifs in promoters of Kmeans clusters", "showMotifKmeans", size = "large"
				   ,radioButtons("radioPromoterKmeans", label = NULL, choices = list("Upstream 300bp as promoter" = 300, "Upstream 600bp as promoter" = 600),selected = 300)
				   ,tableOutput("KmeansPromoter"))
				  ,bsModal("modalExample9", "Determining the number of clusters (k)", "NClusters", size = "large",
				  h5("Following the elbow method, one should choose k so that adding another cluster does not substantially reduce the within groups sum of squares."
			,a("Wikipedia", href="https://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set"))
				  ,plotOutput("KmeansNclusters"))

                )
              )       
    )
	
###############################################################################################################################
 
     ,tabPanel("PCA",
              sidebarLayout(
                sidebarPanel(
				radioButtons("PCA_MDS", "Methods", c("Principal Component Analysis"=1, "Pathway Analysis of PCA rotation" =2, 
				"Multidimensional Scaling"=3))
				,downloadButton('downloadPCAData', 'Download Coordinates')
 					#sliderInput("nGenes1", label = h4("Most variable genes to include"), min = 40, max = 2000, value = 200,step=20) 
			,br(),br()
			,conditionalPanel("input.PCA_MDS != 2" # only show if PCA or MDS (not pathway)
				,htmlOutput('listFactors')
				,htmlOutput('listFactors2')
			)
			,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/pca/",target="_blank")				
				),
                mainPanel(
                  plotOutput("PCA")
                )
              )       
    )

###############################################################################################################################
 
	     ,tabPanel("DEGs",
              sidebarLayout(
                sidebarPanel(
				h5("Identifying Differential Expressed Genes (DEGs)"),
				conditionalPanel("input.dataFileFormat == 1",
				selectInput("CountsDEGMethod", "Method:", choices = list("DESeq2"= 3,"limma-voom"=2,"limma-trend"=1), selected = 3)				
				,tags$style(type='text/css', "#CountsDEGMethod { width:100%;   margin-top:-12px}")
				)	
				,conditionalPanel("input.dataFileFormat == 2", h5("Using the limma package")				)				
				,fluidRow(
				 column(5,numericInput("limmaPval", label = h5("FDR cutoff"), value = 0.1,min=1e-5,max=1,step=.05)  )
				 ,column(7, numericInput("limmaFC", label = h5("Min fold change"), value = 2,min=1,max=100,step=0.5) )
				) # fluidRow
				,tags$style(type='text/css', "#limmaPval { width:100%;   margin-top:-12px}")
				,tags$style(type='text/css', "#limmaFC { width:100%;   margin-top:-12px}")
				,actionButton("modelAndComparisons", "Select factors and comparisons")
				,tags$head(tags$style("#modelAndComparisons{color: blue;}"))				
				,br(),br() 
				,fluidRow(
				   column(6,actionButton("showVenn", "Venn Diagram") )
					, column(6, downloadButton('download.DEG.data', 'All lists') )
				) # fluidRow
				 #,hr()
				 ,HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />') # a solid line
				 ,htmlOutput("listComparisons")
				 ,h5("Pathway database")
				 ,htmlOutput("selectGO2")
				,tags$style(type='text/css', "#selectGO2 { width:100%;   margin-top:-9px}")
				 ,br()
				,fluidRow(
					column(5,actionButton("showVolcano", "Volcano Plot"))
					,column(5, actionButton("showScatter", "Scatter Plot") ) 
				)
			  
				  ,br(),actionButton("showMotif", "TF binding motifs in promoters")
				  ,tags$style(type='text/css', "#showMotif { width:100%;   margin-top:-12px}")
				 #,radioButtons("radio.promoter", label = NULL, choices = list("Upstream 300bp as promoter" = 300, "Upstream 600bp as promoter" = 600),selected = 300)
				 ,br(),br(),downloadButton('download.selectedHeatmap.data', "Download gene list & data" )
				  ,tags$style(type='text/css', "#download.selectedHeatmap.data { width:100%;   margin-top:-12px}")				 
				,h5("Also try",  a("ShinyGO", href="http://ge-lab.org:3838/go/") )	
				,br(),h4( textOutput("text.limma") )										
			,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/degs/",target="_blank")
				, width = 4),
				
                mainPanel(
                  plotOutput("selectedHeatmap")
				   ,h4("Enriched pathways in differentially expressed genes:")
				   ,tableOutput("geneListGO")
				   ,h4("Top Genes"),tableOutput('geneList')
				   #,h4("Enriched motif in promoters")
				   #,tableOutput("DEG.Promoter")
				   ,bsModal("modalExample1", "Enriched TF binding motifs in promoters of DEGs", "showMotif", size = "large"
				   ,radioButtons("radio.promoter", label = NULL, choices = list("Upstream 300bp as promoter" = 300, "Upstream 600bp as promoter" = 600),selected = 300)
				   ,tableOutput("DEG.Promoter"))
				   ,bsModal("modalExample", "Venn Diagram", "showVenn", size = "large",
						htmlOutput('listComparisonsVenn')
						,plotOutput("vennPlot"))
				   ,bsModal("modalExample21", "Model and comparisons", "modelAndComparisons", size = "large",	
					fluidRow(
						 column(6, htmlOutput('listFactorsDE'))
						,column(6, htmlOutput('listBlockFactorsDE') ) 
					)				   
						
						
						,br(),br()
						,htmlOutput('listModelComparisons')
						,br(),br()
						,actionButton("submitModelButton", "Submit & re-calculate",style="float:center")
						,tags$head(tags$style("#submitModelButton{color: blue;font-size: 16px;}"))
						,br(),br(),h5("Close this window to see results.")
						)
				   ,bsModal("modalExample4", "Volcano plot", "showVolcano", size = "large",
						checkboxInput("volcanoPlotBox", label = "Show interactive version w/ gene symbols", value = FALSE)
						,conditionalPanel("input.volcanoPlotBox == 0",	plotOutput("volcanoPlot") )
						,conditionalPanel("input.volcanoPlotBox == 1",plotlyOutput("volcanoPlotly",width = "550px", height = "550px") )
					)
				   ,bsModal("modalExample5", "Scatter plot", "showScatter", size = "large",
						checkboxInput("scatterPlotBox", label = "Show interactive version w/ gene symbols", value = FALSE)
						,conditionalPanel("input.scatterPlotBox == 0",	plotOutput("scatterPlot") )
						,conditionalPanel("input.scatterPlotBox == 1",plotlyOutput("scatterPlotly",width = "550px", height = "550px") )
				   )
				   
				   
				#  ,bsModal("modalExample25", "Interactive Scatter plot", "showScatterPlotly", size = "large",						 	plotlyOutput("scatterPlotly",width = "550px", height = "550px"))
				  # ,bsModal("modalExample24", "Interactive Volcano plot", "showVolcanoPlotly", size = "large",plotlyOutput("volcanoPlotly",width = "550px", height = "550px"))
                )  
              )       
    )

	
###############################################################################################################################
 
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
							
			,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/pathways/",target="_blank")		
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
					,downloadButton('downloadSelectedPathwayData', 'Download expression data for selected pathway')  )
					,conditionalPanel(" (input.pathwayMethod == 1 | input.pathwayMethod == 3 ) & input.selectGO == 'KEGG'",imageOutput("KeggImage"))
					,conditionalPanel("(input.pathwayMethod == 1 | input.pathwayMethod == 3 ) & input.selectGO != 'KEGG'",plotOutput("selectedPathwayHeatmap")) 
			   
                )
              )       
    )


###############################################################################################################################
   ,tabPanel("Genome", 
              sidebarLayout(
                sidebarPanel(
				h5("This interactive map shows DEGs on the genome. 
				 Red and blue dots represent up- or down-regulated genes, respectively.
				 Mouse over to see gene symbols. Click and drag to zoom in.") 
				,htmlOutput("listComparisonsGenome")
				,tags$style(type='text/css', "#listComparisonsPathway { width:100%;   margin-top:-12px}")	
				,fluidRow(
				 column(6,numericInput("limmaPvalViz", label = h5("Filters: FDR "), value = 0.1,min=1e-5,max=1,step=.05)  )
				 ,column(6, numericInput("limmaFCViz", label = h5("Fold change"), value = 2,min=1,max=100,step=0.5) )
				) # fluidRow	
				,tags$style(type='text/css', "#limmaPvalViz { width:100%;   margin-top:-12px}")
				,tags$style(type='text/css', "#limmaFCViz { width:100%;   margin-top:-12px}")					
				,br(), h5("To identify genomic regions significatly enriched with up- or down-regulated genes, we can use "
				,a("PREDA.", href="https://academic.oup.com/bioinformatics/article/27/17/2446/224332/PREDA-an-R-package-to-identify-regional-variations"),
				"Very slow (5 mins), but may be useful in studying cancer or other diseases that might involve chromosomal gain or loss."	)

				,actionButton("runPREDA", "Run PREDA (5 mins)")
				
				
				#,h4("Chr. Regions")
				#,tableOutput("chrRegionsList")
				),
                mainPanel(	
					plotlyOutput("genomePlotly",height = "700px")
					,bsModal("modalExample111", "Differentially expressed genomic loci", "runPREDA", size="large"
						,fluidRow( 
						column(3,numericInput("RegionsPvalCutoff", label = h5("Min. FDR"), value = 0.01,min=1e-20,max=1,step=.05) ),
						column(3,numericInput("StatisticCutoff", label = h5("Min. Statistic"), min = .2, max = 1.5, value = .5,step=.1) ),
						#column(6,numericInput("nPermutations", label = h5("Permutations"), min = 500, max = 3000, value = 1000,step=100) 
						#)
						
						#) # fluidRow
					#,fluidRow(
						column(3,actionButton("showRegions", "Significant Loci")),
						column(3,actionButton("showGenesRegions", "Genes"))
					) #fluidRow
					     ,tags$style(type='text/css', "#showRegions { width:100%; margin-top: 40px;}")
						,tags$style(type='text/css', "#showGenesRegions { width:100%; margin-top: 40px;}")
					,plotOutput("genomePlot")
					)
					,bsModal("modalExample15", "Diff. expressed regions", "showRegions", size = "large",downloadButton('downloadRegions', 'Download'),dataTableOutput("chrRegions"))
					,bsModal("modalExample16", "Diff. expressed regions", "showGenesRegions", size = "large",downloadButton('downloadGenesInRegions', 'Download'),dataTableOutput("genesInChrRegions"))

                )
              )       
    )
	
 #  ,tabPanel("Code",
 #     fluidRow(
 #      column(12,
 #         "R codes for reference",
 #           sourceCode1 <- aceEditor("code"
  #, value = paste(readLines("server.R"), collapse="\n")
  #, mode = "r"
  #, theme = "ambience"
  #, height = "400px"
  #, readOnly = TRUE
 # ) )))
 
###############################################################################################################################
# 
,tabPanel("Bicluster",
              sidebarLayout(
                sidebarPanel(
					h5("Biclustering can discover genes correlated on subset of samples. Only useful when  sample size is large(>10). Uses methods implemented in the biclust R package. ")
					,numericInput("nGenesBiclust", label = h5("Most variable genes to include "), min = 10, max = 2000, value = 1000) 
					,selectInput("biclustMethod", "Method:", choices = list( #"QUBIC"= "BCQU()"
																			#,"runibic"= "BCUnibic()"
																			"BCCC" = "BCCC()"
																			,"BCXmotifs"="BCXmotifs()"
																			, "BCPlaid" = "BCPlaid()"
																			,"BCSpectral" = "BCSpectral()"
																			,"BCBimax" = "BCBimax()"
																			,"BCQuest" = "BCQuest()" 
																			), selected = "BCCC()")				
	
					,htmlOutput('listBiclusters')
					,h5("Enrichment database")
					,htmlOutput("selectGO4"),tags$style(type='text/css', "#selectGO4 { width:100%;   margin-top:-9px}")							
					,downloadButton('download.biclust.data', 'Download all biclusters')
					,br(),br(),textOutput('biclusterInfo')
					,tags$head(tags$style("#biclusterInfo{color: blue;
											 font-size: 16px;
											 font-style: italic;
											 }"
									 )
					)
					,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/biclustering/",target="_blank")
				),
				mainPanel(	
					plotOutput('biclustHeatmap')
					,h3("Enriched gene sets in selected bicluster")
					,tableOutput('geneListBclustGO')
					,br(),br()
					,h3("Genes in this cluster")
					,tableOutput('geneListBicluster')
				
				))
				

) 

###############################################################################################################################
# 
,tabPanel("Network",
              sidebarLayout(
                sidebarPanel(
					h5("Identify co-expression networks and sub-modules using",a( "WGCNA.", href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559",target="_blank" )  ,"Only useful when  sample size is large(>15). ")
					,numericInput("nGenesNetwork", label = h5("Most variable genes to include "), min = 10, max = 1000, value = 1000) 
				,fluidRow(
				 column(6, numericInput("mySoftPower", label = h5("Soft Threshold"), min = 1, max = 20, value = 6))
				 ,column(6, numericInput("minModuleSize", label = h5("Min. Module Size"), min = 10, max = 100, value = 20)  )
				) # fluidRow	
				,tags$style(type='text/css', "#mySoftPower { width:100%;   margin-top:-12px}")
				,tags$style(type='text/css', "#minModuleSize { width:100%;   margin-top:-12px}")				
					,actionButton("chooseSoftThreshold", "Choose soft threshold")	
					,actionButton("showModuleHeatmap", "Heatmap")					
					,downloadButton('download.WGCNA.Module.data',"Download all Modules")
				,HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />') # a solid line
				,htmlOutput('listWGCNA.Modules')
				,fluidRow(
				 column(6, numericInput("edgeThreshold", label = h5("Edge Threshold"), min = 0, max = 1, value =.5))
				 ,column(6, numericInput("topGenesNetwork", label = h5("Top genes"), min = 10, max = 1000, value = 50)  )
				) # fluidRow	
				,tags$style(type='text/css', "#mySoftPower { width:100%;   margin-top:-12px}")
				,tags$style(type='text/css', "#minModuleSize { width:100%;   margin-top:-12px}")	
				,h5("Enrichment database")
				,htmlOutput('selectGO5')
				,downloadButton('downloadSelectedModule',"Download selected network")	
					,h5("The network file can be imported to",a(" VisANT", href="http://visant.bu.edu/",target="_blank"),
					" or ", a("Cytoscape.",href="http://www.cytoscape.org/",target="_blank" )  )				
					),
				mainPanel(	
					
					plotOutput('modulePlot')
					,br(),br()
					,h3("Enriched gene sets in selected module")					
					,tableOutput('networkModuleGO')
					,bsModal("modalExample112", "Choose soft threshold", "chooseSoftThreshold", size="large",plotOutput('softPower') )
					,bsModal("modalExample116", "Heatmap of identified modules", "showModuleHeatmap", size="large",plotOutput('networkHeatmap') )
	

				
				))
				

) 
###############################################################################################################################
 
,tabPanel("R",
      fluidRow(
       column(12,
     htmlOutput('RsessionInfo')
 ) ))

  #,tags$head(includeScript("ga.js")) # tracking usage  
  )# Navibar

)
