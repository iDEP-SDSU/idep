library('R6')
library(shiny)
library(shinyBS)

View.PreProcess <- R6Class("View.PreProcess")


###############################################################################
###################				Main Panel					###################
###############################################################################

# Main Panel Layout of this page includes 2 conditional panels.
# If data is not ready, the message panel will show up.
# If data is ready, the ResultPanel panel will show up.

View.PreProcess$set("public", "mainPanel",
	function(){
		mainPanel(
			self$PreRequestNotFitMessagePanel(),
			self$ResultPanel()
		)
	}
)

###################				sidebar Panel				###################
# Pre process setting panel.
# This panel have several conditional panels.
# These panels will contains controls and settings only work for given data type

View.PreProcess$set("public", "sidebarPanel",
	function(){
		sidebarPanel(
			self$GuessSpeciesPanel(),
			self$ConPanel_FPKMDdataOnlySettings(),
			self$ConPanel_ReadCountOnlySettings(),
			self$ShareSettingsPanel()
		)
	}
)




###############################################################################
###################			Component Functions				###################
###############################################################################


###################				Warning Panel				###################

View.PreProcess$set("public", "PreRequestNotFitMessagePanel",
	function(){
		conditionalPanel(condition = "!output.DataSource ",
			fluidPage(
				h3("Pre process requires load data and select data source type first")
			)
		)
	}
)


###################				Work Panel					###################
# The main UI of Preprocessing result. 
View.PreProcess$set("public", "ResultPanel",
	function(){
		conditionalPanel(condition = "output.DataSource!=null",
			conditionalPanel(condition = "input.dataFileFormat == 1",
				plotlyOutput("PreProcess_ReadCount")
			),
			hr(),
			plotlyOutput("PreProcess_DistTransform"),
			hr(),
			plotlyOutput("PreProcess_DensityTransform"),
			hr(),
			plotlyOutput("PreProcess_ScatterPlot")
		)
	}
)



###################			Read Count Related Settings		###################
View.PreProcess$set("public", "ConPanel_ReadCountOnlySettings",
	function(){
		conditionalPanel(condition = "input.dataFileFormat == 1",
			strong("Keep genes with minimal counts per million (CPM) in at least n libraries:"),
          	fluidRow(
            	column(6, numericInput("numMinCounts", label = h5("Min. CPM"), value = 0.5) ),
            	column(6, numericInput("numNMinSamplesInCountCase", label = h5("n libraries"), value = 1) )
			),
			tags$style(type='text/css', "#numMinCounts { width:100%;   margin-top:-12px}"),
			tags$style(type='text/css', "#numNMinSamplesInCountCase { width:100%;   margin-top:-12px}"),

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
	}

)


###################			FPKM Data Related Settings		###################
View.PreProcess$set("public", "ConPanel_FPKMDdataOnlySettings",
	function(){
		conditionalPanel(condition = "input.dataFileFormat == 2",
			strong("Only keep genes above this level in at least n samples:" ),
			fluidRow(
            	column(6, numericInput("numMinFPKM", label = h5(" Min. level"), value = -1000)),
            	column(6, numericInput("numNMinSampleInFPKMCase", label = h5("n samples"), value = 1) )
			),
			tags$style(type='text/css', "#numMinFPKM { width:100%;   margin-top:-12px}"),
			tags$style(type='text/css', "#numNMinSampleInFPKMCase { width:100%;   margin-top:-12px}"),
			radioButtons("isApplyLogTransFPKM", "Log Transformation",c("No"=FALSE,"Yes"=TRUE) ),
          	numericInput("numFPKMLogStart", label = h5("Constant c for started log: log(x+c)"), value = 1),
			tags$style(type='text/css', "#numFPKMLogStart { width:100%;   margin-top:-12px}"),
          	textOutput("textTransform"),
			tags$head( 
				tags$style("#textTransform{color: blue;
                    font-size: 16px;
                     font-style: italic;}"
            	)
			)
		)
	}
)

###################			Shared  Preprocess Settings		###################
View.PreProcess$set("public", "ShareSettingsPanel",
	function(){
		wellPanel(
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
        	conditionalPanel("input.dataFileFormat == 1", 
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

###################			Guess Species Panel				###################
View.PreProcess$set("public", "GuessSpeciesPanel",
	function(){
		wellPanel(
			strong("Verify guessed species. Change if neccessary."),
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
			tableOutput('PreProcess_tblSpecies')
		)
	}
)




















