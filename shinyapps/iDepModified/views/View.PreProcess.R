library('R6')
library(shiny)
library(shinyBS)

View.PreProcess <- R6Class("View.PreProcess")
View.PreProcess$set("public", "ShowAdvanceSetting", FALSE)


View.PreProcess$set("public", "mainPanel",
	function(){
		fluidPage(
			self$ErrorMessagePanel(),
			self$WorkPanel()
		)
	}
)


View.PreProcess$set("public", "ErrorMessagePanel",
	function(){
		conditionalPanel(condition = "output.DataSource==null",
			fluidPage(
				h3("Pre process requires load data first")
			)
		)
	}
)

View.PreProcess$set("public", "WorkPanel",
	function(){
		conditionalPanel(condition = "output.DataSource!=null",
			fluidPage(
				#self$InstructionPanel(),
				bsCollapse(id = "clsp_PreProcessingSetting", multiple=TRUE, open="None",
					bsCollapsePanel(title="Click for Pre-Process Settings",
						self$PreprocessSettingsPanel()
					)
				),
				self$PreprocessPlotPanel()
			)
		)
	}
)

View.PreProcess$set("public", "InstructionPanel",
	function(){
		wellPanel(
			h5("some description text"),
			actionButton("btn_PreProcess_ChangeSettings", TxtLibrary$btn_label_ChangePreprocessSetting)
		)
	}
)


View.PreProcess$set("public", "PreprocessSettingsPanel",
	function(){
		wellPanel(
			h5("some setting here")
		)
	}
)


View.PreProcess$set("public", "PreprocessPlotPanel",
	function(){
		wellPanel(			
			fluidRow(
				column(width = 5,
      				plotlyOutput("p1")
    			),
				column(width = 5, offset = 1,
					plotlyOutput("p2")
				)
			),
			fluidRow(
				column(width = 5,
      				plotlyOutput("p3")
    			),
				column(width = 5, offset = 1,
					plotlyOutput("p4")
			 	)
			)
		)
	}
)


