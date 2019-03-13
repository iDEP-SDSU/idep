library('R6')

View.Heatmap <- R6Class("View.Heatmap")


###############################################################################
########################		Main Structure		###########################
###############################################################################

# Side Bar
# Side bar contains 3 major panel:
#	1. allow user select how many samples are included in the heat map
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
			#plotOutput("Heatmap_MainPlot"),
			#self$PopShowCorrelation(),
			#self$PopShowSampleTree(),
			#self$PopShowStaticHeatmap(),
			#self$PopShowGeneSDHeatmap()
		)
	}
)


###############################################################################
###################			Component Functions				###################
###############################################################################

# a slider input
View.Heatmap$set("public", "SampleSelection",
	function(){
		wellPanel(
			sliderInput("nGenes", 
                label = h4("Most variable genes to include:"), 
                min   = 0, 
                max   = 12000, 
                value = 1000,
                step  = 100
			) 
		)
	}
)

# pop up panel triggers
View.Heatmap$set("public", "ActionButtons",
	function(){
		wellPanel(
			strong("More plots:"),
			br(),
			actionButton("btn_ShowGeneSDHeatmap", "Gene SD distribution"),
			actionButton("btn_ShowStaticHeatmap", "Interactive heatmap"),
			br(),
			actionButton("btn_ShowCorrelation", "Correlation matrix"),
			actionButton("btn_ShowSampleTree", "Sample Tree")
		)
	}
)


View.Heatmap$set("public", "HeatmapSettingAndDownload",
	function(){
		wellPanel(
			strong("Customize hierarchical clustering (Default values work well):"),
			fluidRow(
				column(3, h5("Color")  ),
				column(9, selectInput("heatColors1", label = NULL,"green-black-red",width='100%') )
			),
			fluidRow(
			 	column(4, h5("Distance")  ),
			 	column(8, selectInput("distFunctions", label = NULL,"Correlation",width='100%') )
			),
			fluidRow(
			  	column(4, h5("Linkage")  ),
			  	column(8, selectInput("hclustFunctions", label = NULL,"average",width='100%') )
			),
			fluidRow(
			  	column(8, h5("Cut-off Z score")  ),
			  	column(4, numericInput("heatmapCutoff", label = NULL, value = 4,min=2,step=1) )
			),
			checkboxInput("geneCentering", "Center genes (substract mean)", value = TRUE),
			checkboxInput("geneNormalize", "Normalize genes (divide by SD)", value = FALSE),
			checkboxInput("sampleCentering", "Center samples (substract mean)", value = FALSE),
			checkboxInput("sampleNormalize", "Normalize samples(divide by SD)", value = FALSE),
			checkboxInput("noSampleClustering", "Do not re-order or cluster samples", value = FALSE),
			htmlOutput('listFactorsHeatmap'),
			downloadButton('downloadData', 'Heatmap data'),
			downloadButton('downloadHeatmap1', 'High-resolution figure'),
			br(),
			a(h5("?",align = "right"), href="https://idepsite.wordpress.com/heatmap/",target="_blank")
		)
	}
)



