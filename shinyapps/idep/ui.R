# New version UI


source('server.config')
source('views/UIManager.R')
library(plotly)


LoadDataView <- View.LoadData$new()
PreProcessView <- View.PreProcess$new()
HeatmapView <- View.Heatmap$new()
ReportView <- View.Report$new()


shinyUI(
	navbarPage(
  		CONFIG_SERVER_VERSION,
		id='navBar',

		#==========================================
		#				Load Data Tab
		#==========================================
		tabPanel("Load Data",
			sidebarLayout(
				LoadDataView$sidebarPanel(),
				LoadDataView$mainPanel()
			) 
		),

		#==========================================
		#				Pre-process
		#==========================================
		tabPanel("Pre-Process", 
			sidebarLayout(
				PreProcessView$sidebarPanel(),
				PreProcessView$mainPanel()
			) 
		),

		#==========================================
		#				Heatmap
		#==========================================
		tabPanel("Heatmap", 
			sidebarLayout(
				HeatmapView$sidebarPanel(),
				HeatmapView$mainPanel()
			) 
		),

		
		#==========================================
		#				Report
		#==========================================
		tabPanel("Report", 
			sidebarLayout(
				ReportView$sidebarPanel(),
				ReportView$mainPanel()
			) 
		)
	)
)



