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
			self$ShareSettingsPanel(),
			
			actionButton("btn_PreProcess_ChangeSettings", TxtLibrary$btn_label_ChangePreprocessSetting)
		)
	}
)

###################				Output Panel				###################
View.PreProcess$set("public", "PreprocessPlotPanel",
	function(){
		fluidRow(			
			fluidRow(
				conditionalPanel(condition = "output.DataSourceType == 1",
					column(width = 5,
      					plotlyOutput("PreProcess_ReadCount")
    				)
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
            	column(6, numericInput("numNMinSamplesInCountCase", label = h5("n libraries"), value = 1) ),
	
				radioButtons("selectCountsTransform", 
					"Transform counts data for clustering & PCA.",  
					c(	
						"EdgeR: log2(CPM+c)"                  = 1,
						"VST: variance stabilizing transform" = 2, 
						"rlog: regularized log (slow) "       = 3
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
            	column(6, numericInput("numMinFPKM", label = h5(" Min. level"), value = -1000)),
            	column(6, numericInput("numNMinSampleInFPKMCase", label = h5("n samples"), value = 1) ),
				radioButtons("isApplyLogTransFPKM", "Log Transformation",c("No"=FALSE,"Yes"=TRUE) ),
          		numericInput("numFPKMLogStart", label = h5("Constant c for started log: log(x+c)"), value = 1),
          		textOutput("textTransform") 
			)
		)
	}
)

###################			Shared  Preprocess Settings		###################
View.PreProcess$set("public", "ShareSettingsPanel",
	function(){
		fluidPage(
			selectInput("selectMissingValueImputationMethod", 
        		label   = "Missing values imputation:",
				choices = list(	"Gene median" = "geneMedian",
								"Treat as zero" = "treatAsZero", 
								"Median within sample groups" = "geneMedianInGroup"),
        	    selected = "geneMedian"),
        	actionButton("btn_GenePlot1", "Plot one or more genes"),
        	br(),br(),
        	actionButton("btn_ExamineDataB", "Search processed data"),
        	br(),br(),
        	checkboxInput("isNoIDConversion", "Do not convert gene IDs to Ensembl.", value = FALSE),
        	downloadButton('downloadProcessedData', 'Processed data'),
        	conditionalPanel("output.DataSourceType == 1", 
        	   downloadButton('downloadConvertedCounts', 'Converted counts data') ),
        	downloadButton('downloadEDAplot', 'High-resolution figure'),  ## this need change name later
        	br(),br(),
        	textOutput('nGenesFilter'),
        	tags$head(tags$style("#nGenesFilter{color: blue;
        	               font-size: 16px;
        	               font-style: italic;
        	               }")),                 
        	textOutput("readCountsBias"),
        	tags$head(tags$style("#readCountsBias{color: red;
        	         font-size: 16px;
        	         font-style: italic;
        	         }" ) ),
        	a(h5("?",align = "right"), href="https://idepsite.wordpress.com/pre-process/",target="_blank")
		)
	}
)






















