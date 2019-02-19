library('R6')
library(shiny)
library(shinyBS)

View.PreProcess <- R6Class("View.PreProcess")


###############################################################################
########################		Main Layout			###########################
###############################################################################

# Main Layout of this page includes 2 conditional panels.

# If data is not ready, the message panel will show up.
# If data is ready, the work panel will show up.

# Major UI design is in WorkPanel.

View.PreProcess$set("public", "mainPanel",
	function(){
		fluidPage(
			self$PreRequestNotFitMessagePanel(),
			self$WorkPanel()
		)
	}
)


###############################################################################
###################			Component Functions			#######################
###############################################################################

View.PreProcess$set("public", "PreRequestNotFitMessagePanel",
	function(){
		conditionalPanel(condition = "output.DataSource==null || output.DataSourceType==null",
			fluidPage(
				h3("Pre process requires load data and select data source type first")
			)
		)
	}
)


###################				Work Panel					###################
# The main UI of Preprocessing tab. 
# This pannel contains two part: 
#	1. A collaspe panel, which allows user to adjust values
#	2. A result panel, which includes multiple plot

View.PreProcess$set("public", "WorkPanel",
	function(){
		conditionalPanel(condition = "output.DataSource!=null && output.DataSourceType!=null",
			fluidPage(
				bsCollapse(id = "clsp_PreProcessingSetting", open="None",
					bsCollapsePanel(title="Click Here for Pre-Process Settings",
						style="info",
						self$PreprocessSettingsPanel()
					)
				),
				self$PreprocessPlotPanel()
			)
		)
	}
)

###################				Setting Panel				###################
# Pre process setting panel.
# This panel have several conditional panels.
# These panels will contains controls and settings only work for given data type

View.PreProcess$set("public", "PreprocessSettingsPanel",
	function(){
		wellPanel(
			self$ConPanel_ReadCountOnlySettings(),
			self$ConPanel_FPKMDdataOnlySettings(),
			
			actionButton("btn_PreProcess_ChangeSettings", TxtLibrary$btn_label_ChangePreprocessSetting)
		)
	}
)

###################				Output Panel				###################
View.PreProcess$set("public", "PreprocessPlotPanel",
	function(){
		fluidRow(			
			fluidRow(
				column(width = 5,
      				plotlyOutput("PreProcess_ReadCount")
    			),
				column(width = 5, offset = 1,
					plotlyOutput("PreProcess_DistTransform")
				)
			),
			fluidRow(
				column(width = 5,
      				plotlyOutput("PreProcess_DensityTransform")
    			),
				column(width = 5, offset = 1,
					plotlyOutput("PreProcess_ScatterPlot")
			 	)
			)
		)
	}
)

###################			Read Count Sepcial Settings		###################
View.PreProcess$set("public", "ConPanel_ReadCountOnlySettings",
	function(){
		conditionalPanel(condition = "output.DataSourceType == 1",
			strong("Keep genes with minimal counts per million (CPM) in at least n libraries:"),
          	fluidRow(
            	column(6, numericInput("numMinCounts", label = h5("Min. CPM"), value = 0.5) ),
            	column(6, numericInput("numNMinSamples", label = h5("n libraries"), value = 1) ),
	
				radioButtons("selectCountsTransform", 
					"Transform counts data for clustering & PCA.",  
					c(	"VST: variance stabilizing transform" = 2, 
						"rlog: regularized log (slow) "       = 3,
						"EdgeR: log2(CPM+c)"                  = 1
					),                          
					selected = 1 
				),
			
				conditionalPanel("input.selectCountsTransform == 1",
					fluidRow(
						column(5, h5("Pseudo count c:")  ),
						column(7, numericInput("numCountsLogStart", label = NULL, value = 4) )
					)                
				)
			)
		)
	}

)


###################			FPKM Data Sepcial Settings		###################
View.PreProcess$set("public", "ConPanel_FPKMDdataOnlySettings",
	function(){
		conditionalPanel(condition = "output.DataSourceType == 2",
			strong("Only keep genes above this level in at least n samples:" ),
			fluidRow(
            	column(6, numericInput("lowFilter", label = h5(" Min. level"), value = -1000)),
            	column(6, numericInput("NminSamples2", label = h5("n samples"), value = 1) )
				          
          	)
		)
	}
)



