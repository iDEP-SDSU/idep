library('R6')
library(shiny)
library(shinyBS)

View.LoadData <- R6Class("View.LoadData")


###############################################################################
########################		Main Layout			###########################
###############################################################################

# Main Layout of this page includes 3 conditional panels and 2 popout panel.

# These 3 panels are controlled by output.DataSource and output.DataSourceType.
# By the conditions, these 3 conditional panels will exclude each other.
# 	1. ConPanel_DataSource shows when data source is not picked. This panel helps  
# 		user select the data source they want use.
# 	2. ConPanel_DataSourceType shows when data source is picked, but data type is
#		not specified. This panel guild user to pick the data source type.
#	3. ConPanel_ViewData shows when both data source and data type are picked. 
#		This panel give user a preview about the uploaded data.

# Pop_DownloadPublicData carries 'reads' app in legacy code. Which allows user
# 	to search and import data from public source. 

# Pop_UploadClientData is part of load data UI in legacy code. We seperate it 
#	into a popout window to make the main UI more clear.


View.LoadData$set("public", "mainPanel",
	function(){
		fluidPage(
			self$ConPanel_DataSource(),
			self$ConPanel_DataSourceType(),
			self$ConPanel_ViewData(),

			self$Pop_DownloadPublicData(),
			self$Pop_UploadClientData()
		)
	}
)









###############################################################################
###################			Component Functions			#######################
###############################################################################


View.LoadData$set("public", "ConPanel_DataSource",
	function(){
		conditionalPanel(condition = "output.DataSource==null",
			fluidPage(
				h5(TxtLibrary$LoadData_Help_Message),
				actionButton("btn_LoadData_DemoData", TxtLibrary$btn_label_LoadDemoData),
				actionButton("btn_LoadData_SearchDataFromPublic", TxtLibrary$btn_label_SearchPublicData),
				actionButton("btn_LoadData_UploadFileFromClient", TxtLibrary$btn_label_UploadClientData)
			)
		)
	}
)

View.LoadData$set("public", "ConPanel_DataSourceType",
	function(){
		conditionalPanel(condition = "output.DataSource!=null && output.DataSourceType==null",
			fluidPage(
				radioButtons(   
					"dataFileFormat", 
					label = TxtLibrary$`Choose data type`,
					choiceNames = list(
						TxtLibrary$`Read counts data (recommended)`,
						TxtLibrary$`Normalized expression values (RNA-seq FPKM, microarray, etc.)`,
						TxtLibrary$`Fold-changes and corrected P values from CuffDiff or any other program`
					), 
					choiceValues = list(
						1, 
						2,
						3
					),
					selected = 1
				),
				conditionalPanel(
					"input.dataFileFormat == 3",
					checkboxInput("isNoFDR", TxtLibrary$`Fold-changes only, no corrected P values`, value = FALSE)
				),
				actionButton("btn_LoadData_ConfirmDataSourceType", TxtLibrary$btn_label_Confirm)
			)
		)
	}

)

View.LoadData$set("public", "ConPanel_ViewData",
	function(){
		conditionalPanel(condition = "output.DataSource!=null && output.DataSourceType!=null",
			fluidPage(
				h3("Success!"),
				h5(TxtLibrary$ViewData_Help_Message),
				actionButton("btn_LoadData_MoveToPreprocess", TxtLibrary$btn_label_MoveToPreProcess),
				actionButton("btn_Reset_Data", TxtLibrary$btn_label_ResetData, onclick="javascript:history.go(0)"),
				br(),
				conditionalPanel(condition="output.tbl_TestDesign!=null",
					h3("Test Design"),
					tableOutput('tbl_TestDesign'),
					br()
				),
				h3("Top Rows"),
				tableOutput('tbl_RawDataTop20')
			)
		)
	}
)

View.LoadData$set("public", "Pop_DownloadPublicData",
	function(){
		bsModal(
			"modalSearchPublicDataTab", 
			"Public RNA-Seq and ChIP-Seq data", 
			"btn_LoadData_SearchDataFromPublic", 
			size="large",
			fluidPage(			
				fluidRow( 
					column(	9, 
							h5("Search and click on a dataset to see more information before download.") 
					),
        	       	column(	3, 
							selectInput(
								"selected.species.archs4", 
								"", 
        	        			choices = list("Human"= "human", "Mouse" = "mouse"), 
								selected = "human"
							)
					)
				),
				DT::dataTableOutput('SearchData'),
				HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />'),
				br(),
				fluidRow( 	
					column(4, 	textOutput('selectedDataset')),
	               	column(3, 	actionButton('btn_LoadData_UseSelectedPublicData', 'Use This Sample') )	
	     		),
				br(),
				tableOutput('ViewData_SelectFromPublic_Samples'),
				h4("Loading data and R packages ... ..."),
				htmlOutput('DoneLoading')
			)
		)
	}
)


View.LoadData$set("public", "Pop_UploadClientData",
	function(){
        bsModal(
			"modalUploadUserDataTab", 
			"Upload user data", 
			"btn_LoadData_UploadFileFromClient", 
			size="large",
            fluidPage(
				fileInput(
					'fileUploadedData', 
					TxtLibrary$`Upload expression data (CSV or text)`,
					accept = c(
						'text/csv',
						'text/comma-separated-values',
						'text/tab-separated-values',
						'text/plain',
						'.csv',
						'.tsv'		  
					) 
				),
                fileInput(
					'fileUploadedTestDesign', 
					h5(TxtLibrary$`Optional: Upload an experiment design file(CSV or text)`),
                    accept = c(
                        'text/csv',
                        'text/comma-separated-values',
                        'text/tab-separated-values',
                        'text/plain',
                        '.csv',
                        '.tsv'		  
                    )
                ),
				actionButton("btn_LoadData_UseUploadedData", TxtLibrary$btn_label_UseUploadedData)
            )
        )
	}
)


