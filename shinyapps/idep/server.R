library(shiny)
library(shinyBS)
library(shinyjs)
library(plotly)
library(DT,verbose=FALSE) 		# for renderDataTable


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

		# PreProcess: UI init
		observe({
			updateSelectInput(
				session,
				"selectOrg",
				choices = PreProcessCtrl$InitChoiceSelectOrgUI()
			)
		})

		# PreProcess: Reactive variables
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

		ConvertedPvals <- reactive({
			PreProcessCtrl$GetConvertedPvals(input, ConvertedIDResult(), ReactVars$PreProcessResult())
		})

		AllGeneInfo <- reactive({
			PreProcessCtrl$GetAllGeneInfomation(ConvertedIDResult(), input)
		})

		# Preprocess: Plots, Tables and Downloads
		# Main
		output$PreProcess_ReadCount <- renderPlotly({
			ReactVars$plotly_PreProcess_ReadCount <- PreProcessCtrl$GetTotalReadCountsPlot(ReactVars$PreProcessResult()$rawCount)
			ReactVars$plotly_PreProcess_ReadCount
		})

		output$PreProcess_DistTransform <- renderPlotly({
			ReactVars$plotly_PreProcess_DistTransform <- PreProcessCtrl$GetTransformedDataBoxPlot(ReactVars$PreProcessResult()$dat)
			ReactVars$plotly_PreProcess_DistTransform
		})


		output$PreProcess_DensityTransform <- renderPlotly({
			ReactVars$plotly_PreProcess_DensityTransform <- PreProcessCtrl$GetTransformedDataDensityPlot(ReactVars$PreProcessResult()$dat)
			ReactVars$plotly_PreProcess_DensityTransform
		})

		output$PreProcess_ScatterPlot <- renderPlotly({
			ReactVars$plotly_PreProcess_ScatterPlot <- PreProcessCtrl$GetTransformedDataScatterPlot(ReactVars$PreProcessResult()$dat)
			ReactVars$plotly_PreProcess_ScatterPlot
		})

		# Side bar
		output$PreProcess_tblSpecies <- renderTable({
			PreProcessCtrl$GetGuessSpeciesResult(ConvertedIDResult())
		}, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)

		output$download_PreProcess_ProcessedData <- downloadHandler(
			filename = "Processed_Data.csv",
			content = function(file){
				PreProcessCtrl$SaveConvetedTransformedDataInTempFile(
					input, file, AllGeneInfo(),
					ConvertedTransformedData(), ConvertedIDResult()
				)
			}
		)


		output$download_PreProcess_ConvertedCounts <- downloadHandler(
			filename = "Converted_Counts_Data.csv",
			content = function(file){
				PreProcessCtrl$SaveConvetedReadCountDataInTempFile(
					input, file, AllGeneInfo(),
					ConvertedRawReadcountData(), ConvertedIDResult()
				)
			}
		)

		output$download_PreProcess_EDAplot <- downloadHandler(
			filename = "EDA.zip",
      		content = function(file) {
				PreProcessCtrl$SaveAllPlotsInTempFile(
					file,
					ReactVars$plotly_PreProcess_ReadCount,
					ReactVars$plotly_PreProcess_DistTransform,
					ReactVars$plotly_PreProcess_DensityTransform,
					ReactVars$plotly_PreProcess_ScatterPlot
				)
			}
		)
		# Pop: Plot one or more gene

		output$PreProcess_SingleGenePlot <- renderPlot({
			PreProcessCtrl$GetSingleGenePlot(input, ConvertedTransformedData(), AllGeneInfo())
		})

		output$download_PreProcess_PlotSingleGenes <- downloadHandler(
			filename = "genePlot.eps",
			content = function(file) {
				PreProcessCtrl$SaveSingleGenesPlotEpsInTempFile(
					input, ConvertedTransformedData(),
					AllGeneInfo(), file
				)
			}
		)

		# Pop: Search Processed data

		output$PreProcess_tbl_DT_ConvertedTransformedData <- DT::renderDataTable({
			PreProcessCtrl$GetDataTableOfConvetedTransformedData(input, AllGeneInfo(), ConvertedTransformedData())
		})

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
					ReactVars$PreProcessResult(),
					PreprocessSampleInfoResult(),

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
		
		# Heatmap Pop Sample Tree

		output$Heatmap_SampleTree <- renderPlot({
			HeatmapCtrl$GetSampleTreePlot(input, ReactVars$PreProcessResult())
		})

		output$download_Heatmap_SampleTree <- downloadHandler(
			filename = "sample_tree.eps",
			content = function(file){
				HeatmapCtrl$SaveSampleTreePlotEpsInTempFile(
					file,
					input,
					ReactVars$PreProcessResult()
				)
			}
		)

		############################################################################
		#						0.0		Test R Markdown Report
		############################################################################

		# Since the file name is changing (file format is also changing) for different 
		# user selection. We have to use a dynamic filename.
		# For some reason, we cannot put this dynamic filename function into controller. 
		# Putting it into controller will make the function lost activation. 
		# Therefore, we have to put all naming function in the 'downloadHandler' directly. 
		# 
		# Another possible solution is to pack whole downloadHandler into controller. 
		# This haven't being tested.

		output$download_Report_MainReport <- downloadHandler(
			filename = function(){
					fileType = gsub( "_document", "", input$selectReportOutputFormat )

					if(fileType == 'html'){
						return("iDepReport.html")
					}

					if(fileType == 'pdf'){
						return("iDepReport.pdf")
					}

					return("iDepReport.docx")
				},
			content = function(file){
				ReportCtrl$SaveReportInTempFile(
					file,
					input,
					ReactVars,
					ReactVars$PreProcessResult(),
					PreprocessSampleInfoResult(),
					ConvertedIDResult(),
					AllGeneInfo(),
					ConvertedTransformedData(),
					NULL, # ConvertedRawReadcountData()
					NULL  # ConvertedPvals()
				)
			}
		)

	}
)
