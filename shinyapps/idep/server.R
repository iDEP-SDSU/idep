library(shiny)
library(shinyBS)
library(shinyjs)
library(plotly)


source('server.config')
source('controllers/ControllerManager.R')

ReactVars <- reactiveValues()
RegularVars <- list()

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

		observeEvent(input$fileUploadedData$datapath,{
			LoadDataCtrl$EventHandler_UseUploadedFiles(input, output, session, ReactVars)
		})

		observeEvent(input$btn_Reset_Data,{
			ReactVars$RawTestDesign <- NULL
			ReactVars$RawData <- NULL
			RegularVars <- list()
		})

		output$txt_PackageLoadedMessage <- renderUI(LoadDataCtrl$GetPackageLoadedMessage())
		outputOptions(output, "txt_PackageLoadedMessage", suspendWhenHidden = FALSE)

		output$tbl_TestDesign <- renderTable(ReactVars$RawTestDesign)
		output$tbl_RawDataTop20 <- renderTable(ReactVars$RawData[1: min(20, nrow(ReactVars$RawData)) ,])



		############################################################################
		#   					1.2  Pre Process
		############################################################################

		observe({  updateSelectInput(session, "selectOrg", choices = QuerySpeciesListFromConvertDBOrgInfo() ) })

		ReactVars$PreProcessResult <- reactive({
			PreProcessCtrl$PreProcessResult(input, session, ReactVars)
		})

		ReactVars$ConvertedIDResult <- reactive({
			PreProcessCtrl$ConvertedIDResult(ReactVars$PreProcessResult()$dat, input$selectOrg )
		})

		output$PreProcess_ReadCount <- renderPlotly({
			PreProcessCtrl$GetTotalReadCountsPlot(ReactVars$PreProcessResult()$rawCount)
		})

		output$PreProcess_DistTransform <- renderPlotly({
			PreProcessCtrl$GetTransformedDataBoxPlot(ReactVars$PreProcessResult()$dat)
		})
		

		output$PreProcess_DensityTransform <- renderPlotly({
			PreProcessCtrl$GetTransformedDataDensityPlot(ReactVars$PreProcessResult()$dat)
		})

		output$PreProcess_ScatterPlot <- renderPlotly({
			PreProcessCtrl$GetTransformedDataScatterPlot(ReactVars$PreProcessResult()$dat)
		})

		output$tblSpecies <- renderTable({
			PreProcessCtrl$GetGuessSpeciesResult(ReactVars$ConvertedIDResult())
		}, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T) )

		

		############################################################################
		#						0.0		Test R Markdown Report
		############################################################################



	}
)
