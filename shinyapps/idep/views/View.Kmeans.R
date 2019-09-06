library('R6')

View.Kmeans <- R6Class("View.Kmeans")

###############################################################################
########################		Main Structure		###########################
###############################################################################

# Side Bar

View.Kmeans$set("public", "sidebarPanel", 
	function(){
		sidebarPanel(
			self$ClusterConfiguration(),
			actionButton("btn_Kmeans_rerun", "Re-Run")
			self$PopupTriggers(),
			self$DownloadPlotSection(),
			# a solid line as divider
			HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />'),
			self$PathwayEnrichment()
		)
	}
)

# Main panel

# Besid main plot, Kmeans page contains 5 different pop-out pages:
# 1. Enriched TF binding motifs in promoters of Kmeans clusters, aka PopPromoters
# 2. Determining the number of clusters (k), aka PopNumberOfClusters
# 3. TSNE plot of genes, aka PopTSNEPlot
# 4. Distribution of variations among genes, aka PopGeneDistribution
# 5. Visualize enrichment, aka PopVisualizeEnrichment

View.Kmeans$set("public", "mainPanel",
	function(){
		mainPanel(
			plotOutput("Kmeans_Heatmap", inline=TRUE),
			br(),
			h4("Enriched pathways for each cluster"),
			tableOutput("tbl_Kmeans_GO"),
			self$PopPromoters(), 
			self$PopNumberOfClusters(),
			self$PopTSNEPlot(),
			self$PopGeneDistribution(),
			self$PopVisualizeEnrichment()
		)
	}
)


###############################################################################
###################			Component Functions				###################
###############################################################################

# Configuration which will affect heatmap change
View.Kmeans$set("public", "ClusterConfiguration",
	function(){
		wellPanel(
			sliderInput(
				"num_Kmeans_GenesKNN", 
				label = h4("Most variable genes to include "),
				min   = 0, 
				max   = 12000, 
				value = 2000,
				step  = 100
			),
			sliderInput(
				"num_Kmeans_Culsters",
				label = h4("Number of Clusters"), 
				min   = 2, 
				max   = 20, 
				value = 4,
				step  = 1 
			),
			selectInput(
				"select_Kmeans_Normalization", 
				h5("Normalize by gene:"), 
				choices = list(
					"Mean center" = "geneMean",
					"Standardization" = "geneStandardization",
					"L1 Norm"	= "L1Norm"
				),
				selected = "Standardization"
			),
			tags$style(type='text/css', "#select_Kmeans_Normalization { width:100%;   margin-top:-9px}")
		)
	}
)

# Buttons trigger pop out windows
View.Kmeans$set("public", "PopupTriggers",
	function(){
		wellPanel(
			actionButton("btn_Kmeans_showNumberOfClusters", "How many clusters?"),
			actionButton("btn_Kmeans_ShowGeneSD", "Gene SD distribution"),
			actionButton("btn_Kmeans_ShowTSNE", "t-SNE map")
		)
	}
)

# Download section
View.Kmeans$set("public", "DownloadPlotSection",
	function(){
		downloadButton('download_Kmeans_KmeansData', 'K-means data'),
        downloadButton('download_Kmeans_Heatmap', 'High-resolution figure')
	}
)

# Enrichment part
View.Kmeans$set("public", "PathwayEnrichment",
	function(){
		h5("Pathway database"),
        htmlOutput("select_Kmeans_PathwayDatabase"),
        tags$style(type='text/css', "#select_Kmeans_PathwayDatabase { width:100%;   margin-top:-9px}"),
		checkboxInput("is_Kmeans_RemoveRedudantSets", "Remove redudant genesets", value = TRUE),
		actionButton("btn_Kmeans_ShowEnrichmentPlot", "Visualize enrichment"),
		downloadButton("download_Kmeans_EnrichementPlot", "Enrichment details" ),
		a( 	
			h5("?",align = "right"), 
			href="https://idepsite.wordpress.com/k-means/",
			target="_blank"
		)
	}
)

# Pop out windows definitions
View.Kmeans$set("public", "PopPromoters",
	function(){
		bsModal(
			"modal_Kmeans_ShowPromoters",
			"Enriched TF binding motifs in promoters of Kmeans clusters",
			"btn_Kmeans_showMotifKmeans",
			size = "large",
			radioButtons("select_PromoterKmeans", 
                        label    = NULL, 
                        choices  = list("Upstream 300bp as promoter" = 300, 
                                        "Upstream 600bp as promoter" = 600),
                        selected = 300),
          	tableOutput("tbl_KmeansPromoter")
		)
	}
)


View.Kmeans$set("public", "PopNumberOfClusters",
	function(){
		bsModal(
			"modal_Kmeans_ShowNumberOfClusters",
			"Determining the number of clusters (k)", 
			"btn_Kmeans_showNumberOfClusters",
			size = "large",
			h5("Following the elbow method, one should choose k so that adding another cluster 
         		does not substantially reduce the within groups sum of squares." ,
          		a(	"Wikipedia", 
				  	href="https://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set",
          			target="_blank"
				)
			),
			plotOutput("Kmeans_Nclusters")
		)
	}
)


View.Kmeans$set("public", "PopTSNEPlot",
	function(){
		bsModal(
			"modal_Kmeans_TSNEPlot",
			"t-SNE plot of genes",
			"btn_Kmeans_ShowTSNE",
			size = "large",
			h5("We use the dimension reduction algorith ",
				a("t-SNE", href="https://lvdmaaten.github.io/tsne/",target="_blank"), 
				"to map the top genes. Examine the distribution can help choose the nubmer of clusters in k-Means. "
			),
			checkboxInput("is_Kmeans_ColorGenes", "Color genes by the results of k-Means", value = TRUE),
			actionButton("btn_Kmeans_Recalculate", "Re-calculate using different random numbers"),
			plotOutput("Kmeans_tSNEgenePlot")
		)
	}
)

View.Kmeans$set("public", "PopGeneDistribution",
	function(){
		bsModal(
			"modal_Kmeans_GeneDistribution",
			"Distribution of variations among genes",
			"btn_Kmeans_ShowGeneSD",
			size = "large",
			downloadButton('download_Kmeans_GeneDistribution',"Download Figure"),
			plotOutput("Kmeans_GeneDistribution")
		)
	}
)

View.Kmeans$set("public", "PopVisualizeEnrichment",
	function(){
		bsModal(
			"modal_Kmeans_VisualizeEnrichment",
			"Visualize enrichment",
			"btn_Kmeans_ShowEnrichmentPlot",
			size = "large",
			h5("Gene sets closer on the tree share more genes. Sizes of dot correspond to adjuested Pvalues"),
			downloadButton('download_Kmeans_EnrichmentPlot',"Download Figure"),
			plotOutput('Kmeans_EnrichmentPlot')
		)
	}
)


