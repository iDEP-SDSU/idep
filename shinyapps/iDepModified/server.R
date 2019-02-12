library(shiny)
library(shinyBS)
library(shinyjs)


source('server.config')
source('controllers/ControllerManager.R')

ReactVars <- reactiveValues()

LoadDataCtrl <- Ctl.LoadData$new()

shinyServer(
	function(input, output, session) {

		############################################################################
		#   					1.0  Read data Main UI
		############################################################################





		############################################################################
		#						1.1	Download public data UI
		############################################################################
		output$ViewData_SelectFromPublic_Samples <- LoadDataCtrl$RenderSampleTable(input)
		output$SearchData <- LoadDataCtrl$SearchGSEIDs(input)
		output$humanNsamplesOutput <- LoadDataCtrl$RenderHumanNsampleOutput(input)
		output$mouseNsamplesOutput <- LoadDataCtrl$RenderMouseNsampleOutput(input)
		output$selectedDataset <- LoadDataCtrl$RenderSelectedDataset(input)
		output$DoneLoading <- LoadDataCtrl$RenderInitDoneUI()






		observeEvent(input$btn_LoadData_UseSelectedPublicData,{
			LoadDataCtrl$EventHandler_UseSelectedPublicData(input, output, session, ReactVars)
		})

		observeEvent(input$btn_LoadData_DemoData,{
			LoadDataCtrl$EventHandler_UseDemoData(input, output, session, ReactVars)
		})

		observeEvent(input$btn_label_UseUploadedData,{
			LoadDataCtrl$EventHandler_UseUploadedFiles(input, output, session, ReactVars)
		})

		observeEvent(input$btn_LoadData_MoveToPreprocess,{
			updateNavbarPage(session, "iDepNav", selected = "Pre Process")
		})

		output$tbl_TestDesign <- renderTable(ReactVars$RawTestDesign)
		output$tbl_RawDataTop20 <- renderTable(ReactVars$RawData[1: min(20, nrow(ReactVars$RawData)) ,])

	}
)
