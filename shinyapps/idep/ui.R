# iDEP user interface
library(shiny,verbose=FALSE)
library("shinyAce",verbose=FALSE) # for showing text files, code
library(shinyBS,verbose=FALSE) # for popup figures
library(plotly,verbose=FALSE)
iDEPversion = "iDEP.68"
# 0.38 Gene ID conversion, remove redudancy;  rlog option set to blind=TRUE
# 0.39 reorganized code. Updated to Bioconductor 3.5; solved problems with PREDA 9/8/17
# 0.40 moved libraries from the beginning to different places to save loading time
   # change colors for heatmap 
   # fix error with voom; also changed method to TMM, see https://www.bioconductor.org/help/workflows/RNAseq123/#transformations-from-the-raw-scale
 #   
shinyUI(
  
  navbarPage(
iDEPversion,
			id='navBar',
            tabPanel("Load Data",
	     	titlePanel("Upload Files "),
    sidebarLayout(
      sidebarPanel(
	    actionButton("goButton", "Click here to load demo data"),
		    tags$head(tags$style("#goButton{color: red;
                                 font-size: 16px;
                                 font-style: italic;
                                 }"
                         )
              )	
		,h5(" and just click the tabs for some magic!",  style = "color:red")
		,p(HTML("<div align=\"right\"> <A HREF=\"javascript:history.go(0)\">Reset</A></div>" ))
		,radioButtons("dataFileFormat", label = "1. Choose data type", choices = list("Read counts data (recommended)" = 1, 
													"Normalized expression values (RNA-seq FPKM, microarray, etc.)" = 2,
													"Fold-changes with FDRs from CuffDiff or any other program" =3),selected = 1)
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
		,fileInput('file2', h5('Optional: Upload an experiment design file(CSV or text)'),
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    'text/tab-separated-values',
                    'text/plain',
                    '.csv',
                    '.tsv'					
                  )
        )
		#,div(downloadButton('downloadSampleInfoData', 'Example'), style="float:right")
		#,downloadButton('downloadSampleInfoData', 'Example')
		#,tags$style(type='text/css', "#downloadSampleInfoData { width:100%;   margin-top:-22px}")
		#,br(),br()
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
		#,h3("Dear users, Thank you for using iDEP. Please help make it better by sending us criticisms, suggestions or feature/functionality requests. We also answer your questions on how to use iDEP to analyze your data.",a("Just email us.",href="mailto:Xijin.Ge@SDSTATE.EDU?Subject=iDEP suggestions"))
		,h5("Integrated Differential Expression and Pathway analysis (iDEP) of transcriptomic data.  See ",
			a(" documentation", href="https://idepsite.wordpress.com/",target="_blank"), "and",
			a(" manuscript.", href="http://biorxiv.org/content/biorxiv/early/2017/06/09/148411.full.pdf",target="_blank"),
   			"Based on annotation of",a( " 163 animal and 45 plant genomes ",href="https://idepsite.wordpress.com/species/",target="_blank") ,"in Ensembl BioMart as of 12/15/2017."
 			,a("STRING-db ", href="https://string-db.org/",target="_blank")
			,"offer API access to protein interaction networks and annotations for 115 archaeal, 1678 bacterial, and 238 eukaryotic species."
			," Additional  data from"
			,a("KEGG, ", href="www.genome.jp/kegg/",target="_blank")
			,a("Reactome, ", href="http://www.reactome.org/",target="_blank")
			,a("MSigDB (human),", href="https://doi.org/10.1093/bioinformatics/btr260",target="_blank") 
			,a("GSKB (mouse)", href="http://biorxiv.org/content/early/2016/10/24/082511",target="_blank") 
			,"and",a("  araPath (arabidopsis).", href="https://doi.org/10.1093/bioinformatics/bts421",target="_blank")

			
			," For feedbacks or data contributions (genes and GO mapping of any species), please"
			,a("contact us, ",href="mailto:xijin.ge@sdstate.edu?Subject=iDEP" )
			, "or visit our",a(" homepage.", href="http://ge-lab.org/",target="_blank")
			, "Send us suggestions or any error message to help improve iDEP."
			,a("Email",href="mailto:Xijin.Ge@SDSTATE.EDU?Subject=iDEP suggestions")
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
			     conditionalPanel("input.dataFileFormat == 2",
				 strong("Only keep genes above this level in at least n samples:" )
				 ,fluidRow(
					column(6, numericInput("lowFilter", label = h5(" Min. level"), value = -1000))
					,column(6, numericInput("NminSamples2", label = h5("n samples"), value = 1) )					
				)
					,tags$style(type='text/css', "#lowFilter { width:100%;   margin-top:-12px}")
					,tags$style(type='text/css', "#NminSamples2 { width:100%;   margin-top:-12px}")			
					,radioButtons("transform", "Log Transformation",c("No"=FALSE,"Yes"=TRUE) )
					,numericInput("logStart", label = h5("Constant c for started log: log(x+c)"), value = 1)
					,tags$style(type='text/css', "#logStart { width:100%;   margin-top:-12px}")
					,textOutput("textTransform") 
						,tags$head(tags$style("#textTransform{color: blue;
								 font-size: 16px;
								 font-style: italic;
								 }"
						 ) )
					
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
					"rlog: regularized log (slow) "= 3,"Log transformatoin: log2(x+c)" = 1),selected = 1 )
					,conditionalPanel("input.CountsTransform == 1",
						fluidRow(
							column(5, h5("Constant c:")  )
							,column(7, numericInput("countsLogStart", label = NULL, value = 4) )
						)
					        
					
					)
 				#,radioButtons("CountsDEGMethod", "Detect differntial expression using:",  c("limma-voom"=2,"limma-trend"=1,"DESeq2"= 3) )
				
				)
				,selectInput("missingValue", label = "Missing values inputation:",choices = list("Gene median"="geneMedian",
																							"Treat as zero"="treatAsZero", 
																							"Median within sample groups"="geneMedianInGroup"),
																							selected = "geneMedian")
				,actionButton("genePlot1", "Barplot for one or more genes")
				,br(),br()
				,actionButton("examineDataB", "Search processed data")
                 ,br(),br()
				 ,checkboxInput("noIDConversion", "Do not convert gene IDs to Ensembl.", value = FALSE)
				 ,downloadButton('downloadProcessedData', 'Processed data') 
				 ,conditionalPanel("input.dataFileFormat == 1",downloadButton('downloadConvertedCounts', 'Converted counts data') )

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
					,actionButton("showSampleTree", "Sample Tree")						
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
					,checkboxInput("geneCentering", "Center genes (substract mean)", value = TRUE)					
					,checkboxInput("geneNormalize", "Normalize genes (divide by SD)", value = FALSE)
					,checkboxInput("sampleCentering", "Center samples (substract mean)", value = FALSE)	
					,checkboxInput("sampleNormalize", "Normalize samples(divide by SD)", value = FALSE)
					,checkboxInput("noSampleClustering", "Do not re-order or cluster samples", value = FALSE)
					,htmlOutput('listFactorsHeatmap')

					,downloadButton('downloadData', 'Heatmap data')
					,downloadButton('downloadHeatmap1', 'High-resolution figure')
						
					,br(),a(h5("?",align = "right"), href="https://idepsite.wordpress.com/heatmap/",target="_blank")
						),
					mainPanel(				  

						plotOutput("heatmap1")
						
						
						# ,verbatimTextOutput("event")
						,bsModal("modalExample8", "Correlation matrix using top 75% genes", "showCorrelation", size = "large"
							,downloadButton("downloadCorrelationMatrix")
							,checkboxInput("labelPCC", "Label w/ Pearson's correlation coefficients", value = TRUE)
							,plotOutput("correlationMatrix"))
						,bsModal("modalExample228", "Hierarchical clustering tree.", "showSampleTree", size = "large"
							,h4("Using genes whose maximum expression level is at top 75%. Data is transformed and clustered as specified in the main page. ")
							,plotOutput("sampleTree"))

						,bsModal("modalExample28", "Heatmap with hierarchical clustering tree", "showStaticHeatmap", size = "large",
						sliderInput("nGenesPlotly", label = h4("Most variable genes to include:"), min = 10, max = 6000, value = 50,step=50),
						h5("Mouse over to see gene names. To zoom, click and drag up or downward and release."),
						plotlyOutput("heatmapPlotly",width = "100%", height = "800px"))
					 
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
				,actionButton("KmeansReRun", "Re-Run")
				,actionButton("NClusters", "How many clusters?")
				,actionButton("showGeneSD", "Gene SD distribution")
				,actionButton("geneTSNE", "t-SNE map")
				,selectInput("kmeansNormalization", h5("Normalize by gene:"), choices = list("Mean center"="geneMean","Standardization"= "geneStandardization","L1 Norm"= "L1Norm"), selected = "geneMean")	
				,tags$style(type='text/css', "#kmeansNormalization { width:100%;   margin-top:-9px}")
				,actionButton("showMotifKmeans", "Enriched TF binding motifs")
				,br(),downloadButton('downloadDataKmeans', 'K-means data')
				,downloadButton('downloadKmeansHeatmap', 'High-resolution figure')
				,HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />') # a solid line
				,h5("Pathway database")
				,htmlOutput("selectGO3")
				,tags$style(type='text/css', "#selectGO3 { width:100%;   margin-top:-9px}")
				,checkboxInput("removeRedudantSets", "Remove redudant genesets", value = TRUE)
				,actionButton("ModalEnrichmentPlotKmeans", "Visualize enrichment")
				,downloadButton('downloadKmeansGO',"Enrichment details")						 
			,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/k-means/",target="_blank")
				),
                mainPanel(
                  plotOutput("KmeansHeatmap",inline=TRUE)
				  ,br()
				  ,h4("Enriched pathways for each cluster")
				  ,tableOutput("KmeansGO")
				  ,bsModal("modalExample2", "Enriched TF binding motifs in promoters of Kmeans clusters", "showMotifKmeans", size = "large"
				   ,radioButtons("radioPromoterKmeans", label = NULL, choices = list("Upstream 300bp as promoter" = 300, "Upstream 600bp as promoter" = 600),selected = 300)
				   ,tableOutput("KmeansPromoter"))
				  ,bsModal("modalExample9", "Determining the number of clusters (k)", "NClusters", size = "large",
				  h5("Following the elbow method, one should choose k so that adding another cluster does not substantially reduce the within groups sum of squares."
					,a("Wikipedia", href="https://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set",target="_blank"))
				  ,plotOutput("KmeansNclusters"))
				  ,bsModal("modalExample229", "t-SNE plot of genes", "geneTSNE", size = "large",
				  h5("We use the dimension reduction algorith "
					,a("t-SNE", href="https://lvdmaaten.github.io/tsne/",target="_blank"), "to map the top genes. Examine the distribution can help choose the nubmer of clusters in k-Means. ")
				  ,checkboxInput("colorGenes", "Color genes by the results of k-Means", value = TRUE)
				  ,actionButton("seedTSNE", "Re-calculate using different random numbers")
				  ,plotOutput("tSNEgenePlot"))
				  ,bsModal("modalExample233", "Distribution of variations among genes", "showGeneSD", size = "large", plotOutput('distributionSD'))
				  ,bsModal("ModalEnrichmentPlotKmeans1", "Visualize enrichment", "ModalEnrichmentPlotKmeans", size="large"
						,h5("Gene sets closer on the tree share more genes. Sizes of dot correspond to adjuested Pvalues")
						,downloadButton('enrichmentPlotKmeans4Download',"Figure"),plotOutput('enrichmentPlotKmeans')
				   
				   )
                )
              )       
    )
	
###############################################################################################################################
 
     ,tabPanel("PCA",
              sidebarLayout(
                sidebarPanel(
				radioButtons("PCA_MDS", "Methods", c("Principal Component Analysis"=1, 
				"Multidimensional Scaling"=3,"t-SNE"=4, "Pathway Analysis of PCA rotation" =2 ))



			,conditionalPanel("input.PCA_MDS == 2" # only show if PCA or MDS (not pathway)
				,htmlOutput('selectGO6')
			)
			,conditionalPanel("input.PCA_MDS != 2" # only show if PCA or MDS (not pathway)
				,htmlOutput('listFactors')
				,htmlOutput('listFactors2')
			)
			
			
			
			
				,br()
			,downloadButton('downloadPCAData', 'Coordinates')
 					#sliderInput("nGenes1", label = h4("Most variable genes to include"), min = 40, max = 2000, value = 200,step=20) 

			,downloadButton('downloadPCA', 'High-resolution figure')
			,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/pca/",target="_blank")				
				),
                mainPanel(
				
                  plotOutput("PCA", inline=TRUE)
				  
				  ,conditionalPanel("input.PCA_MDS == 4", # only show if PCA or MDS (not pathway)
					actionButton("tsneSeed2", "Re-calculate t-SNE"),br(),br())
				
				  ,br(),br()
				 ,conditionalPanel("input.PCA_MDS == 1 |input.PCA_MDS == 2 " # only show if PCA or MDS (not pathway)
					,htmlOutput('PCA2factor')
					)
                )
              )       
    )

###############################################################################################################################
 
	     ,tabPanel("DEG1",
              sidebarLayout(
                sidebarPanel(
				h5("Identifying Differential Expressed Genes (DEGs). See next tab for details."),
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
				,actionButton("modelAndComparisons", "Select factors & comparisons")
				,tags$head(tags$style("#modelAndComparisons{color: blue;font-size: 15px;}"))				
				,br(),br() 
				,fluidRow(
				   column(6,actionButton("showVenn", "Venn Diagram") )
				#, column(6,actionButton("showDEGstats", "Barplot"))
						)
					,br(),fluidRow(					
					column(8, downloadButton('downloadGeneListsGMT', 'Gene lists') ) )
					,br(),fluidRow(
					 column(8, downloadButton('download.DEG.data', 'Gene lists & Data') )
					) # fluidRow
									,br(),h4( textOutput("textLimma") )
				,tags$head(tags$style("#textLimma{color: blue;font-size: 15px;}"))	
			,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/degs/",target="_blank")
				, width = 4),
				
                mainPanel(
					plotOutput('sigGeneStats')
					,br(),br(),h4("Numbers of differentially expressed genes for all comparisons. \"B-A\" means B vs. A. Interaction terms starts with \"I:\" ")
					,tableOutput('sigGeneStatsTable')


				   #,h4("Enriched motif in promoters")
				   #,tableOutput("DEG.Promoter")

				   ,bsModal("modalExample", "Venn Diagram", "showVenn", size = "large",
						checkboxInput("UpDownRegulated", label = "Split gene lists by up- or down-regulation", value = FALSE)
						,htmlOutput('listComparisonsVenn')
						,plotOutput("vennPlot"))
				   ,bsModal("modalExample21", "Build model and/or select comparisons", "modelAndComparisons", size = "large",

					fluidRow(
						 column(6, htmlOutput('listFactorsDE'))
						,column(6, htmlOutput('listBlockFactorsDE') ) 
					)
			
					,fluidRow(
						 column(6, htmlOutput('selectReferenceLevels1'))
						,column(6, htmlOutput('selectReferenceLevels2') ) 
					)
					,fluidRow(
						 column(6, htmlOutput('selectReferenceLevels3'))
						,column(6, htmlOutput('selectReferenceLevels4') ) 
					)
					,fluidRow(
						 column(6, htmlOutput('selectReferenceLevels5'))
						,column(6, htmlOutput('selectReferenceLevels6') ) 
					)						
					,htmlOutput('listInteractionTerms')
					,textOutput('experimentDesign')
					,tags$head(tags$style("#experimentDesign{color: red;font-size: 16px;}"))		
						,htmlOutput('listModelComparisons')
						,actionButton("submitModelButton", "Submit & re-calculate",style="float:center")
						,tags$head(tags$style("#submitModelButton{color: blue;font-size: 20px;}"))
						,h5("Close this window to see results.")
			,a(h5("More info on DESeq2 experiment design",align = "right"), href="http://rpubs.com/ge600/deseq2",target="_blank")
						)
				   #,bsModal("modalExample56", "Summary of differentially expressed genes", "showDEGstats", size = "large",

					#)				   
				   
				#  ,bsModal("modalExample25", "Interactive Scatter plot", "showScatterPlotly", size = "large",						 	plotlyOutput("scatterPlotly",width = "550px", height = "550px"))
				  # ,bsModal("modalExample24", "Interactive Volcano plot", "showVolcanoPlotly", size = "large",plotlyOutput("volcanoPlotly",width = "550px", height = "550px"))
                )  
              )       
    )

	
###############################################################################################################################
	     ,tabPanel("DEG2",
              sidebarLayout(
                sidebarPanel(
				 h5("Examine the results of DEGs for each comparison")
				 ,htmlOutput("listComparisons")
				 ,br()
				,actionButton("showVolcano", "Volcano Plot")
					,actionButton("showMAplot", "MA Plot")  

				
			     ,actionButton("showScatter", "Scatter Plot") 
				 ,actionButton("showMotif", "TF binding motifs in promoters")
				 # ,tags$style(type='text/css', "#showMotif { width:100%;   margin-top:-12px}")
				 ,downloadButton('download.selectedHeatmap.data', "Gene list & data" )
				  ,tags$style(type='text/css', "#download.selectedHeatmap.data { width:100%;   margin-top:-12px}")				 
				 ,downloadButton('downloadSelectedHeatmap',"High-resolution figure")
				  #,radioButtons("radio.promoter", label = NULL, choices = list("Upstream 300bp as promoter" = 300, "Upstream 600bp as promoter" = 600),selected = 300)
				,HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />') # a solid line
				 ,h5("Enrichment analysis")

				 
				 ,htmlOutput("selectGO2")
				,tags$style(type='text/css', "#selectGO2 { width:100%;   margin-top:-9px}")
				 ,actionButton("ModalEnrichmentPlot", "Enrichment tree")
				 ,actionButton("ModalEnrichmentNetwork", "Enrichment network")			 
				 ,downloadButton('downloadGOTerms', "Enrichment details" )
				 ,actionButton("STRINGdb_GO", "Enrichment using STRING API")
				,tags$head(tags$style("#STRINGdb_GO{color: blue}"))	
				 ,h5("Also try",  a("ShinyGO", href="http://ge-lab.org:3838/go/",target="_blank") )	
			,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/degs/",target="_blank")				
				
				, width=4),
				
				mainPanel(
					#h4("Expression pattern of DEGs for selected comparison:")
                   plotOutput("selectedHeatmap", inline=TRUE)
				   ,br()
				   ,h4("Enriched pathways in DEGs for the selected comparison:")				    
				   ,tableOutput("geneListGO")
				   ,h4("Top Genes for selected comparison:"),tableOutput('geneList')
				   ,bsModal("modalExample4", "Volcano plot", "showVolcano", size = "large",
						checkboxInput("volcanoPlotBox", label = "Show interactive version w/ gene symbols", value = FALSE)
						,conditionalPanel("input.volcanoPlotBox == 0",	downloadButton('downloadVolcanoPlot',"Figure" ),plotOutput("volcanoPlot")
						 )
						,conditionalPanel("input.volcanoPlotBox == 1",plotlyOutput("volcanoPlotly",width = "550px", height = "550px") )
					)
				   ,bsModal("modalExample5", "Scatter plot", "showScatter", size = "large",
						checkboxInput("scatterPlotBox", label = "Show interactive version w/ gene symbols", value = FALSE)
						,conditionalPanel("input.scatterPlotBox == 0",	downloadButton('downloadScatterPlot',"Figure" ),	plotOutput("scatterPlot") )
						,conditionalPanel("input.scatterPlotBox == 1",plotlyOutput("scatterPlotly",width = "550px", height = "550px") )
				   )
				   ,bsModal("ModalEnrichmentPlot1", "Visualize enrichment", "ModalEnrichmentPlot", size="large"
						,h5("Gene sets closer on the tree share more genes. Sizes of dot correspond to adjuested Pvalues")
						,downloadButton('enrichmentPlotDEG24Download',"Figure")
						,plotOutput('enrichmentPlotDEG2')
				   
				   )
				   ,bsModal("ModalEnrichmentPlot2", "Visualize enrichment", "ModalEnrichmentNetwork", size="large"
						,h5("Connected gene sets share more genes. Color of node correspond to adjuested Pvalues.")
						,checkboxInput("enrichmentNetworkInteractive", label = "Interactive version", value = FALSE)
						,conditionalPanel("input.enrichmentNetworkInteractive==0" 
								,h5("To change network layout, zoom in or out from your browser. Download repeatedly for different layout.")
								,downloadButton("enrichmentNetworkPlot4Download", "High-resolution figure")
								,plotOutput('enrichmentNetworkPlot'))
						,conditionalPanel("input.enrichmentNetworkInteractive==1" ,plotlyOutput('enrichmentNetworkPlotly',width = "900px", height = "800px"))				   
				   )
				   ,bsModal("modalExample5555", "M-A plot", "showMAplot", size = "large",
						checkboxInput("MAPlotBox", label = "Show interactive version w/ gene symbols", value = FALSE)
						,conditionalPanel("input.MAPlotBox == 0",	downloadButton('downloadMAPlot',"Figure" ),	plotOutput("MAplot") )
						,conditionalPanel("input.MAPlotBox == 1",plotlyOutput("MAplotly",width = "550px", height = "550px") )
				   )
				   ,bsModal("modalExampleSTRING", "Enrichment and network visualization using STRING API", "STRINGdb_GO", size = "large"
						,h5("iDEP tries to match your species with the 115 archaeal, 1678 bacterial, and 238 eukaryotic species in the",
							a(" STRING server", href="https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html",target="_blank"),
								" and send the DEGs. If it is running, please wait until it finishes.")
						,htmlOutput("STRINGDB_species_stat") 
						,tags$head(tags$style("#STRINGDB_species_stat{color: blue;font-size: 15px;}"))						
						, selectizeInput('speciesName', label=NULL,choices = " ",
							   multiple = TRUE, options = list(maxItems = 1,							 
								 placeholder = 'Species name (e.g. Homo sapiens)',
								 onInitialize = I('function() { this.setValue(""); }')
							   )
							 )
						,textOutput('STRINGDB_mapping_stat')
						,tags$head(tags$style("#STRINGDB_mapping_stat{color: blue;font-size: 15px;}"))
						,br()
				,actionButton("ModalPPI","PPI network of DEGs"),br(),br()
						,selectInput("STRINGdbGO", label="Functional Enrichment", choices = list("GO Biological Process"= "Process"
																			  ,"GO Cellular Component"= "Component"
																			  ,"GO Molecular Function"= "Function"
																			  ,"KEGG" = "KEGG"
																			  ,"Pfam" = "Pfam"
																			  ,"InterPro" = "InterPro")
																, selected = "Process")							
						#,actionButton("submit2STRINGdb", "Submit")
						,downloadButton("STRING_enrichmentDownload")
						,tableOutput("stringDB_GO_enrichment")
   				   )
				   ,bsModal("ModalExamplePPI", "Analyze top DEGs on protein interaction networks (PPIs)", "ModalPPI", size = "large"
						,h5("By sending top DEGs (ranked by fold-change) to the STRING website, 
						    iDEP is retrieving a sub-network, calculating PPI enrichment, 
						  and generating custom URLs to the STRING website containing your genes. This can take 5 minutes. Patience will pay off! ")
						,sliderInput("nGenesPPI", label = h5("Top up- or down-regulated genes to include:"), min = 0, max = 400, value = 100,step=10) 
						,htmlOutput("stringDB_network_link")
						,tags$head(tags$style("#stringDB_network_link{color: blue; font-size: 15px;}"))
						,br(),br(),h5("Interactions among proteins encoded by top up-regulated proteins",align = "center")
						,plotOutput("stringDB_network1")		 	   
				   
				   )

				   	,bsModal("modalExample1", "Enriched TF binding motifs in promoters of DEGs", "showMotif", size = "large"
				   ,radioButtons("radio.promoter", label = NULL, choices = list("Upstream 300bp as promoter" = 300, "Upstream 600bp as promoter" = 600),selected = 300)
				   ,tableOutput("DEG.Promoter"))
				   
                )  
              )       
    )								
###############################################################################################################################
 
   ,tabPanel("Pathway",
              sidebarLayout(
                sidebarPanel(
				htmlOutput("listComparisonsPathway")
				,tags$style(type='text/css', "#listComparisonsPathway { width:100%;   margin-top:-12px}")
			,selectInput("pathwayMethod", label = "Select method:", choices = list("GAGE" = 1, 
																								"GSEA (preranked fgsea)" =3,
																								"PGSEA" = 2, 
																								"PGSEA w/ all samples" =4, 
																								"ReactomePA" = 5
																								), selected = 3) #
				,tags$style(type='text/css', "#pathwayMethod { width:100%;   margin-top:-12px}")
				#,h5("Select genesets (use KEGG to show pathway diagrams w/ fold-change):")
				,htmlOutput("selectGO1")
				,tags$style(type='text/css', "#selectGO1 { width:100%;   margin-top:-12px}")
                #,h5("Filter genesets by size:")
				,fluidRow( 
				column(6,numericInput("minSetSize", label = h5("Geneset size: Min."), 
						min = 5, max = 30, value = 15,step=1) ),
				column(6,numericInput("maxSetSize", label = h5("Max."), 
					min = 1000, max = 2000, value = 2000,step=100) )
				) # fluidRow
				,tags$style(type='text/css', "#minSetSize { width:100%;   margin-top:-12px}")
				,tags$style(type='text/css', "#maxSetSize { width:100%;   margin-top:-12px}")
				,numericInput("pathwayPvalCutoff", label = h5("Pathway signifiance cutoff (FDR)"), value = 0.2,min=1e-20,max=1,step=.05)
				,tags$style(type='text/css', "#pathwayPvalCutoff { width:100%;   margin-top:-12px}")
				,numericInput("nPathwayShow", label = h5("Number of top pathways to show"), value = 30, min=5,max=100,step=5)
				,tags$style(type='text/css', "#nPathwayShow { width:100%;   margin-top:-12px}")
				,checkboxInput("absoluteFold", label = "Use absolute values of fold changes for GSEA and GAGE", value = FALSE)
				,numericInput("GenePvalCutoff", label = h5("Remove genes with big FDR before pathway analysis:"), value = 1,min=1e-20,max=1,step=.05)
				,tags$style(type='text/css', "#GenePvalCutoff { width:100%;   margin-top:-12px}")

				,conditionalPanel("input.pathwayMethod == 1 | input.pathwayMethod == 2 | input.pathwayMethod == 3| input.pathwayMethod == 4" 
					,actionButton("ModalEnrichmentPlotPathway", "Pathway tree") 
					,actionButton("ModalEnrichmentNetworkPathway", "Pathway network")
#					,actionButton("ModalExaminePathways", "Gene expression by pathway")
					,downloadButton('downloadPathwayListData', "Pathway list w/ genes")
					
				)


				,conditionalPanel("input.pathwayMethod == 4", downloadButton("PGSEAplotAllSamples.Download", "High-resolution figure") )
				,conditionalPanel("input.pathwayMethod == 2", downloadButton("PGSEAplot.Download", "High-resolution figure") )

				#,actionButton("examinePathway", "Examine individual pathways")
				#,conditionalPanel("input.pathwayMethod == 2",downloadButton('download.PGSEAplot.data', 'Download PGSEA pathway data'))
			,h5("* Warning! The many combinations can lead to false positives in pathway analyses.")				
			,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/pathways/",target="_blank")		
				),
                mainPanel(	  
					
				  conditionalPanel("input.pathwayMethod == 2",
				  	 h5("Red and blue indicates activated and suppressed pathways, respectively. ")				  
					 ,plotOutput("PGSEAplot", inline=TRUE)

				  )
				  ,conditionalPanel("input.pathwayMethod == 1",
                  tableOutput("gagePathway") )
				  ,conditionalPanel("input.pathwayMethod == 3",
                  tableOutput("fgseaPathway") )
				   ,conditionalPanel("input.pathwayMethod == 5",
                  tableOutput("ReactomePAPathway") )
				  ,conditionalPanel("input.pathwayMethod == 4",
				  h5("Red and blue indicates activated and suppressed pathways, respectively.  ") 
				  ,plotOutput("PGSEAplotAllSamples", inline=TRUE)
				  )
				  
				  #,bsModal("ModalExaminePathway", "Gene expression on pathways", "ModalExaminePathways", size="large"				  
				     ,conditionalPanel("input.pathwayMethod == 1 | input.pathwayMethod == 2
											| input.pathwayMethod == 3|input.pathwayMethod == 4" 
				      ,htmlOutput("listSigPathways")
					  ,downloadButton('downloadSelectedPathwayData', 'Expression data for genes in selected pathway')  )
					  ,conditionalPanel(" (input.pathwayMethod == 1 | input.pathwayMethod == 2| 
											input.pathwayMethod == 3 |input.pathwayMethod == 4) & input.selectGO == 'KEGG'"
							,h5("Red and green represent up- and down-regulated genes, respectively.")
							,imageOutput("KeggImage", width = "100%", height = "100%"))
					  ,conditionalPanel(" (input.pathwayMethod == 1 | input.pathwayMethod == 2 
						|input.pathwayMethod == 3 |input.pathwayMethod == 4) & input.selectGO != 'KEGG'",
						plotOutput("selectedPathwayHeatmap")) 
					#)

			,bsModal("ModalEnrichmentPlotPathway1", "Significant pathways", "ModalEnrichmentPlotPathway", size="large"
						,h5("Gene sets closer on the tree share more genes. Sizes of dot correspond to adjuested Pvalues")
						,downloadButton('enrichmentPlotPathway4Download',"Figure")
						,plotOutput('enrichmentPlotPathway')
				   
				   )
				   ,bsModal("ModalEnrichmentPlotPahtway2", "Significant pathways", "ModalEnrichmentNetworkPathway", size="large"
						,h5("Connected gene sets share more genes. Color of node correspond to adjuested Pvalues.")
						,checkboxInput("enrichmentNetworkInteractivePathway", label = "Interactive version", value = FALSE)
						,conditionalPanel("input.enrichmentNetworkInteractivePathway==0"
								,h5("To change network layout, zoom in or out from your browser. Download repeatedly for different layout.")
								,downloadButton("enrichmentNetworkPlotPathway4Download","High-resolution figure")
								,plotOutput('enrichmentNetworkPlotPathway'))
						,conditionalPanel("input.enrichmentNetworkInteractivePathway==1" ,plotlyOutput('enrichmentNetworkPlotlyPathway',width = "900px", height = "800px"))				   
				   )				   
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
				,a("PREDA.", href="https://academic.oup.com/bioinformatics/article/27/17/2446/224332/PREDA-an-R-package-to-identify-regional-variations",target="_blank"),
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
					,selectInput("biclustMethod", "Method:", choices = list( 
																			"BCCC" = "BCCC()"
																			,"QUBIC"= "BCQU()"
																			,"runibic"= "BCUnibic()"																		
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
					,numericInput("nGenesNetwork", label = h5("Most variable genes to include (<3001) "), min = 10, max = 3000, value = 1000) 
				,fluidRow(
				 column(6, numericInput("mySoftPower", label = h5("Soft Threshold"), min = 1, max = 20, value = 5))
				 ,column(6, numericInput("minModuleSize", label = h5("Min. Module Size"), min = 10, max = 100, value = 20)  )
				) # fluidRow	
				,tags$style(type='text/css', "#mySoftPower { width:100%;   margin-top:-12px}")
				,tags$style(type='text/css', "#minModuleSize { width:100%;   margin-top:-12px}")				
					,actionButton("chooseSoftThreshold", "Choose soft threshold")	
					,actionButton("showModuleHeatmap", "Heatmap")					

					,textOutput('moduleStatistics')
					,tags$head(tags$style("#moduleStatistics{color: blue;
											 font-size: 14px;
											 font-style: italic;
											 }"
							  ))
													
				,HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />') # a solid line
				,htmlOutput('listWGCNA.Modules')
				,fluidRow(
				 column(6, numericInput("edgeThreshold", label = h5("Edge Threshold"), min = 0, max = 1, value =.4, step=.1))
				 ,column(6, numericInput("topGenesNetwork", label = h5("Top genes"), min = 10, max = 2000, value = 10, step = 10)  )
				) # fluidRow	
				,tags$style(type='text/css', "#edgeThreshold { width:100%;   margin-top:-12px}")
				,tags$style(type='text/css', "#topGenesNetwork { width:100%;   margin-top:-12px}")
				,actionButton("networkLayout", "Change network layout",style="float:center")				
				,h5("Enrichment database:")
				,htmlOutput('selectGO5')
				,downloadButton('download.WGCNA.Module.data',"Download all modules")
				,downloadButton('downloadSelectedModule',"Download network for selected module")	
					,h5("The network file can be imported to",a(" VisANT", href="http://visant.bu.edu/",target="_blank"),
					" or ", a("Cytoscape.",href="http://www.cytoscape.org/",target="_blank" )  )	
					,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/network/",target="_blank")					
					),
				mainPanel(						
					plotOutput('modulePlot')
					,br(),br()
					,plotOutput('moduleNetwork')
					,h3("Enriched pathways among all genes in selected module")					
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
	 ,h4("Changes")
	 ,h5("iDEP v0.62 2/5/2018:  Add  tree and networks to visualize 
		overlaps between enriched gene sets, in K-means, DEG2, and pathway tags.
		Downloads of pathway analysis results and high-resolution figures.")
	 ,h5("2/6/2018: Fixed errors caused by gene symbol matching for unknown species. More user control of hierarchical clustering tree")
	 ,h5("2/9/2018: V 0.65 Added API access to STRINGdb website on the DEG2 tab. Supports thousands of bacterial species")
	 ,h5("2/10/2018: V 0.66 Improved API access to STRINGdb, by adding automatic species matching.")
	 ,h5("2/11/2018: V 0.67 Tested with larger dataset of 259 samples. Changed figure configurations.")
	 ,h5("2/14/2018: V 0.68 Fixed Pathview loading code. Connected Pathview to PGSEA.")
	 ,h5("In loving memory of my parents.")
 ) ))

  ,tags$head(includeScript("ga.js")) # tracking usage  
  )# Navibar

)
