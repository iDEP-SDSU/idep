library('R6')

View.LoadData <- R6Class("View.LoadData")

View.LoadData$set("puiblic", "mainPanel",
	function(){
		mainPanel(
			self$ConPanel_DataSource(),
			self$ConPanel_ViewData(),

			self$Pop_DownloadPublicData(),
			self$Pop_UploadClientData()
		)
	}
)

View.LoadData$set("puiblic", "ConPanel_DataSource",
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

View.LoadData$set("puiblic", "ConPanel_ViewData",
	function(){
		conditionalPanel(condition = "output.DataSource!=null",
			
		)
	}
)



