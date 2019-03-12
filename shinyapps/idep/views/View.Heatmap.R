library('R6')

View.Heatmap <- R6Class("View.Heatmap")


###############################################################################
########################		Main Structure		###########################
###############################################################################

# Side Bar
# Side bar contains 3 major panel:
#	1. allow user select how many sample included in the heat map
#	2. contains buttons for pop up panels
#	3. heatmap setting and downloading function
View.Heatmap$set( "public", "sidebarPanel", 
	function(){
		sidebarPanel(
			self$SampleSelection(),
			self$ActionButtons(),
			HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />'), # a solid line
			self$HeatmapSettingAndDownload()			
		)
	}
)

# Main panel
View.Heatmap$set("public", "mainPanel",
	function(){
		mainPanel(
			plotOutput("Heatmap_MainPlot"),
			self$PopShowCorrelation(),
			self$PopShowSampleTree(),
			self$PopShowStaticHeatmap(),
			self$PopShowGeneSDHeatmap()
		)
	}
)




###############################################################################
###################			Component Functions				###################
###############################################################################



