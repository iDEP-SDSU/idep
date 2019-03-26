library(shiny)
library(shinyBS)
library(shinyjs)
library(plotly)


source('server.config')
source('controllers/ControllerManager.R')


ReactVars <- reactiveValues()
RegularVars <- list()

shinyServer(
	function(input, output, session) {

		############################################################################
		#   					1.1  		Read data
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
		#   					1.2  		Pre Process
		############################################################################

		observe({  updateSelectInput(session, "selectOrg", choices = PreProcessCtrl$InitChoiceSelectOrgUI() ) })

		ReactVars$PreProcessResult <- reactive({
			PreProcessCtrl$PreProcessResult(input, session, ReactVars)
		})

		PreprocessSampleInfoResult <- reactive({
			PreProcessCtrl$RawSampleInfoPreprocess(ReactVars$RawTestDesign, ReactVars$PreProcessResult()$dat)
		})

		ConvertedIDResult <- reactive({
			geneNames <- rownames(ReactVars$PreProcessResult()$dat)
			PreProcessCtrl$ConvertedIDResult(geneNames, input$selectOrg )
		})

		ConvertedTransformedData <- reactive({
			PreProcessCtrl$GetConvertedTransformedData(input, ConvertedIDResult(), ReactVars$PreProcessResult())
		})

		ConvertedRawReadcountData <- reactive({
			PreProcessCtrl$GetConvertedRawReadcountData(input, ConvertedIDResult(), ReactVars$PreProcessResult())
		})

		GetConvertedPvals <- reactive({
			PreProcessCtrl$GetConvertedPvals(input, ConvertedIDResult(), ReactVars$PreProcessResult())
		})

		AllGeneInfo <- reactive({
			PreProcessCtrl$GetAllGeneInfomation(ConvertedIDResult(), input)
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

		output$PreProcess_tblSpecies <- renderTable({
			PreProcessCtrl$GetGuessSpeciesResult(ConvertedIDResult())
		}, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)


		############################################################################
		#   					1.3  		Heatmap
		############################################################################

		# Heatmap: UI init
		observe({
			HeatmapCtrl$RefreshUI_Select_Heatmap_FactorsHeatmap(session, output, ReactVars, PreprocessSampleInfoResult())
		})

		observe({
			updateSelectInput(
				session,
				"select_Heatmap_MainPlot_DistanceFun",
				choices = HeatmapCtrl$InitSelectDistFunctionChoices()
			)
		})

		observe({
			updateSelectInput(
				session,
				"select_Heatmap_MainPlot_HeatColor",
				choices = HeatmapCtrl$InitSelectHeatColorChoices()
			)
		})

		observe({
			updateSelectInput(
				session,
				"select_Heatmap_MainPlot_HClustFun",
				choices = HeatmapCtrl$InitSelectClusterFunctionChoices()
			)
		})

		# Heatmap: Main Panel
		output$Heatmap_MainPlot <- renderPlot(
			{
				HeatmapCtrl$GetMainHeatmap(
					input,
					ReactVars,
					ReactVars$PreProcessResult(),
					PreprocessSampleInfoResult()
				)
			},
			height = 900
		)

		# Heatmap: Side Bar
		output$download_Heatmap_DownloadHeatmapData <- downloadHandler(
			filename = "heatmap.csv",
			content = function(file) {
				HeatmapCtrl$SaveHeatmapDataInFile(
					file, input,
					ReactVars$PreProcessResult(), PreprocessSampleInfoResult() 
				)
			}
		)

		output$download_Heatmap_DownloadEpsFormatPlot <- downloadHandler(
			filename = "heatmap.eps",
			content = function(file) {
				HeatmapCtrl$SaveMainPlotEpsInTempFile(
					file, input, ReactVars,
					ReactVars$PreProcessResult(), PreprocessSampleInfoResult()
				)
			}
		)

		# Heatmap PopShowStaticHeatmap
		output$Heatmap_HeatmapPlotly <- renderPlotly({
			HeatmapCtrl$GetMainHeatmapPlotly(
				input,
				ReactVars,
				ReactVars$PreProcessResult(),
				PreprocessSampleInfoResult(),
				AllGeneInfo()
			)
		})

		# Heatmap PopShowGeneSDHeatmap
		output$Heatmap_PopShowGeneSDHeatmap <- renderPlot({
			HeatmapCtrl$GetGeneSDHeatmap(input, ConvertedTransformedData())		
		}, height = 600, width = 800, res=120 )

		output$download_Heatmap_PopShowGeneSDHeatmap <- downloadHandler(
			filename = "gene_SD_distribution.eps",
			content = function(file){
				HeatmapCtrl$SaveGeneSDPlotEpsInTempFile(
					file,
					input,
					ConvertedTransformedData()
				)
			}
		)


		# Heatmap PopShowCorrelation
		output$Heatmap_CorrelationMatrix <- renderPlot({
			HeatmapCtrl$GetCorrelationMatrixPlot(input, ReactVars$PreProcessResult())
		})

		output$download_Heatmap_CorrelationMatrixPlot <- downloadHandler(
			filename = "correlation_matrix.eps",
			content = function(file){
				HeatmapCtrl$SaveCorrelationMatrixPlotEpsInTempFile(
					file,
					input,
					ReactVars$PreProcessResult()
				)
			}
		)

		output$download_Heatmap_CorrelationMatrixData <- downloadHandler(
			filename = "correlationMatrix.csv",
			content = function(file){
				HeatmapCtrl$SaveCorrelationMatrixPlotDataInFile(
					file,					
					ReactVars$PreProcessResult()
				)
			}
		)
		
		# Heatmap 
		output$Heatmap_SampleTree <- renderPlot({
			HeatmapCtrl$GetSampleTreePlot(input, ReactVars$PreProcessResult())
		})

		############################################################################
		#						0.0		Test R Markdown Report
		############################################################################



	}
)
