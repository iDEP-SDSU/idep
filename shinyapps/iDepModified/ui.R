# New version UI


source('server.config')
source('views/UIManager.R')
library(plotly)


LoadDataView <- View.LoadData$new()
PreProcessView <- View.PreProcess$new()




shinyUI(
#  	fluidPage(
#	
#  		titlePanel("iDep"),
#  	
#		bsCollapse(id = "collaspeMainPanel", multiple=TRUE, open="Prepare Data",
#			bsCollapsePanel(title="Prepare Data",style="warining",
#				bsCollapse(id = "collaspePrePareData", multiple=TRUE, open="Load Data",
#					bsCollapsePanel(	title="Load Data",style="warining", 
#										LoadDataView$mainPanel() ),
#					bsCollapsePanel(	title="Pre Process",style="warining", 
#										PreProcessView$mainPanel() )
#				)
#			)
#		)
#	)


	navlistPanel(
		id = "iDepNav",
		"Prepare Data",
		widths = c(2, 10),
		tabPanel("Load Data",
			mainPanel(LoadDataView$mainPanel())
		),
		tabPanel("Pre Process",
			mainPanel(PreProcessView$mainPanel())
		),
		"Analysis",
		tabPanel("some analysis method")
		
	)
)
