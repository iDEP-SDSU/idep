# New version UI


source('server.config')
source('views/UIManager.R')


LoadDataView <- View.LoadData$new()




shinyUI(
  fluidPage(
	
  	titlePanel("iDep"),
  	
  	navlistPanel(
  		"Prepare Data",
  		widths = c(2, 10),
  		tabPanel("Load Data",
  				 mainPanel(
  				 )),
  		tabPanel("Pre Process",
  				 mainPanel()),
  		"Analysis",
  		tabPanel("some analysis method"
  		)
  		
  	)
  )
)