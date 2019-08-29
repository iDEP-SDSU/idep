library('R6')

View.Kmeans <- R6Class("View.Kmeans")

###############################################################################
########################		Main Structure		###########################
###############################################################################

# Side Bar

View.Heatmap$set("public", "sidebarPanel", 
	function(){
		sidebarPanel(
			# not finish yet
		)
	}
)

View.Heatmap$set("public", "mainPanel",
	function(){
		mainPanel(
			plotOutput("KmeansHeatmap", inline=TRUE),
			br(),
			self$EnrichedPathways()
			br(),
			self$PopMotifInPromoters(),
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










