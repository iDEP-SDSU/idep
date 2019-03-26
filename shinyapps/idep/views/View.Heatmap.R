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
View.Heatmap$set("public", "sidebarPanel", 
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
			self$PopShowGeneSDHeatmap(),
			self$PopShowStaticHeatmap(),
			self$PopShowCorrelation(),
			self$PopShowSampleTree()	
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
			sliderInput("num_Heatmap_IncludeGeneCount", 
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
			actionButton("btn_Heatmap_ShowGeneSDHeatmap", "Gene SD distribution"),
			actionButton("btn_Heatmap_ShowStaticHeatmap", "Interactive heatmap"),
			br(),
			actionButton("btn_Heatmap_ShowCorrelationMatrix", "Correlation matrix"),
			actionButton("btn_Heatmap_ShowSampleTree", "Sample Tree")
		)
	}
)

# main heatmap settings
View.Heatmap$set("public", "HeatmapSettingAndDownload",
	function(){
		wellPanel(
			strong("Customize hierarchical clustering (Default values work well):"),
			fluidRow(
				column(3, h5("Color")  ),
				column(9, selectInput("select_Heatmap_MainPlot_HeatColor", label = NULL,"green-black-red",width='100%') )
			),
			fluidRow(
			 	column(4, h5("Distance")  ),
			 	column(8, selectInput("select_Heatmap_MainPlot_DistanceFun", label = NULL,"Correlation",width='100%') )
			),
			fluidRow(
			  	column(4, h5("Linkage")  ),
			  	column(8, selectInput("select_Heatmap_MainPlot_HClustFun", label = NULL,"average",width='100%') )
			),
			fluidRow(
			  	column(8, h5("Cut-off Z score")  ),
			  	column(4, numericInput("num_Heatmap_HeatmapCutoff", label = NULL, value = 4,min=2,step=1) )
			),
			checkboxInput("is_Heatmap_GeneCentering", "Center genes (substract mean)", value = TRUE),
			checkboxInput("is_Heatmap_GeneNormalize", "Normalize genes (divide by SD)", value = FALSE),
			checkboxInput("is_Heatmap_SampleCentering", "Center samples (substract mean)", value = FALSE),
			checkboxInput("is_Heatmap_SampleNormalize", "Normalize samples(divide by SD)", value = FALSE),
			checkboxInput("is_Heatmap_NoSampleClustering", "Do not re-order or cluster samples", value = FALSE),
			conditionalPanel(condition = "output.showSelect_Heatmap_FactorsHeatmap",
				selectInput("select_Heatmap_FactorsHeatmap", label="Sample color bar:",choices= c("Sample_Name"), selected = "Sample_Name")
			),
			downloadButton('download_Heatmap_DownloadHeatmapData', 'Heatmap data'),
			downloadButton('download_Heatmap_DownloadEpsFormatPlot', 'High-resolution figure'),
			br(),
			a(h5("?",align = "right"), href="https://idepsite.wordpress.com/heatmap/",target="_blank")
		)
	}
)

# PopShowGeneSDHeatmap
View.Heatmap$set("public", "PopShowGeneSDHeatmap",
	function(){
		bsModal(
			"modal_Heatmap_ShowGeneSDHeatmap",
			"Distribution of variations among genes",
			"btn_Heatmap_ShowGeneSDHeatmap",
			size="large",
			downloadButton("download_Heatmap_PopShowGeneSDHeatmap", "Figure"),
			plotOutput("Heatmap_PopShowGeneSDHeatmap")
		)
	}
)

# PopShowStaticHeatmap
View.Heatmap$set("public", "PopShowStaticHeatmap",
	function(){
		bsModal(
			"modal_Heatmap_ShowStaticHeatmap",
			"Heatmap with hierarchical clustering tree",
			"btn_Heatmap_ShowStaticHeatmap",
			size = "large",
			sliderInput("num_Heatmap_PlotlyIncludeGeneCount", 
               label = h4("Most variable genes to include:"),
               min   = 0, 
               max   = 12000, 
               value = 50,
               step  = 100
			),
           	h5("Mouse over to see gene names. To zoom, click and drag up or downward and release."),
           	plotlyOutput("Heatmap_HeatmapPlotly", width = "100%", height = "800px")
		)
	}
)

# PopShowCorrelation
View.Heatmap$set("public", "PopShowCorrelation",
	function(){
		bsModal(
			"modal_Heatmap_ShowCorrelation",
			"Correlation matrix using top 75% genes",
			"btn_Heatmap_ShowCorrelationMatrix",
			size = "large",
			downloadButton("download_Heatmap_CorrelationMatrixData","Data"),
          	downloadButton('download_Heatmap_CorrelationMatrixPlot', 'Figure'),
			checkboxInput("isLabelWithPCC", "Label w/ Pearson's correlation coefficients", value = TRUE),
			plotOutput("Heatmap_CorrelationMatrix")
		)
	}
)


# PopShowSampleTree
View.Heatmap$set("public", "PopShowSampleTree",
	function(){
		bsModal(
			"modal_Heatmap_ShowSampleTree",
			"Hierarchical clustering tree.",
			"btn_Heatmap_ShowSampleTree",
			size = "large",
			h4("Using genes with maximum expression level at the top 75%. Data is transformed 
           		and clustered as specified in the main page. "),
          	plotOutput("Heatmap_SampleTree"),
          	downloadButton('download_Heatmap_SampleTree', 'Figure')
		)
	}
)







