library(shiny)
library(shinyBS)
library(shinyjs)
library(plotly)


source('server.config')
source('controllers/ControllerManager.R')

ReactVars <- reactiveValues()

LoadDataCtrl <- Ctl.LoadData$new()
PreProcessCtrl <- Ctl.PreProcess$new()

shinyServer(
	function(input, output, session) {

		############################################################################
		#   					1.1  Read data
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

		observeEvent(input$btn_LoadData_UseUploadedData,{
			LoadDataCtrl$EventHandler_UseUploadedFiles(input, output, session, ReactVars)
		})

		observeEvent(input$btn_LoadData_ConfirmDataSourceType,{
			LoadDataCtrl$EventHandler_ConfirmDataSrouceType(input, output, session, ReactVars)
		})

		observeEvent(input$btn_LoadData_MoveToPreprocess,{
			updateNavbarPage(session, "iDepNav", selected = "Pre Process")
		})

		output$tbl_TestDesign <- renderTable(ReactVars$RawTestDesign)
		output$tbl_RawDataTop20 <- renderTable(ReactVars$RawData[1: min(20, nrow(ReactVars$RawData)) ,])



		############################################################################
		#   					1.2  Pre Process
		############################################################################



		output$PreProcess_ReadCount <- renderPlotly({
			PreProcessCtrl$getTotalReadCountsData(input, output, session, ReactVars)
		})
		output$PreProcess_DistTransform <- renderPlotly({
			PreProcessCtrl$testPlotly()
		})
		output$PreProcess_DensityTransform <- renderPlotly({
			PreProcessCtrl$testPlotly()
		})
		output$PreProcess_ScatterPlot <- renderPlotly({
			PreProcessCtrl$testPlotly()
		})
	}
)
