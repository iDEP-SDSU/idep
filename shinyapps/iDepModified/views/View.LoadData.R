library('R6')
library(shiny)

View.LoadData <- R6Class("View.LoadData")

View.LoadData$set("public", "mainPanel",
	function(){
		mainPanel(
			self$ConPanel_DataSource(),
			self$ConPanel_ViewData(),

			self$Pop_DownloadPublicData(),
			self$Pop_UploadClientData()
		)
	}
)

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

View.LoadData$set("public", "ConPanel_ViewData",
	function(){
		conditionalPanel(condition = "output.DataSource!=null",
			fluidPage(
				h5(TxtLibrary$ViewData_Help_Message),
				actionButton("btn_Reset_Data", TxtLibrary$btn_label_ResetData)
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
                a(TxtLibrary$`New! Analyze public RNA-seq data`, href="http://bioinformatics.sdstate.edu/reads/"),
                fileInput('uploadedDataFile', h5(TxtLibrary$`Optional: Upload an experiment design file(CSV or text)`),
                    accept = c(
                        'text/csv',
                        'text/comma-separated-values',
                        'text/tab-separated-values',
                        'text/plain',
                        '.csv',
                        '.tsv'		  
                    )
                )
            )
        )
	}
)


