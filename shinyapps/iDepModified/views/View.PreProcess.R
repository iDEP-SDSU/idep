library('R6')
library(shiny)
library(shinyBS)

View.PreProcess <- R6Class("View.PreProcess")
View.PreProcess$set("public", "ShowAdvanceSetting", FALSE)


View.PreProcess$set("public", "mainPanel",
	function(){
		fluidPage(
			self$InstructionPanel(),
			bsCollapse(id = "collaspePrepareData", multiple=TRUE, open="Prepare Data",
				bsCollapsePanel(title="Prepare Data",style="warining",
					self$PreprocessSettingsPanel()
				)
			),
			self$PreprocessPlotPanel()
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
		
	}
)


View.PreProcess$set("public", "PreprocessPlotPanel",
	function(){
		wellPanel(
			fluidRow(
				column(width = 5,
      				plotlyOutput("p1", width ="400", height ="300")
    			),
				column(width = 5, offset = 1,
					plotOutput("p2")
				)
			),
			plotOutput("p3", width ="400", height ="300"),
			plotOutput("p4", width ="400", height ="200")
		)
	}
)


