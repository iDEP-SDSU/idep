# New version UI


source('server.config')
source('views/UIManager.R')
library(plotly)


LoadDataView <- View.LoadData$new()
PreProcessView <- View.PreProcess$new()




shinyUI(
  	fluidPage(
  		titlePanel("iDep"),
		navlistPanel(
			id = "iDepNav",
			"Prepare Data",
			widths = c(2, 10),
			tabPanel("Load Data",
				width = "80%",
				mainPanel(LoadDataView$mainPanel(), width=10)
			),
			tabPanel("Pre Process",
				width = "80%",
				mainPanel(PreProcessView$mainPanel(), width=10)
			),
			"Analysis",
			tabPanel("some analysis method")
		)
	)
)
