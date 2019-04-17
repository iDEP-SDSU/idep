library('R6')

Ctl.Heatmap <- R6Class("Ctl.Heatmap")

###############################################################################
###################			Load/Reload UI Functions		###################
###############################################################################

Ctl.Heatmap$set("public", "InitSelectDistFunctionChoices", 
	function(){
		distFuns <- LogicManager$UtilFuns$DistanceFuns
		distFunNames <- names(distFuns)
		return(setNames(distFunNames, distFunNames))
	}
)

Ctl.Heatmap$set("public", "InitSelectHeatColorChoices",
	function(){
		heatColors <- LogicManager$Display$HeatColors
		colorNames <- names(heatColors)
		return(setNames(colorNames, colorNames))
	}
)

Ctl.Heatmap$set("public", "InitSelectClusterFunctionChoices",
	function(){
		clusterFuns <- LogicManager$UtilFuns$HierarchicalClusteringFuns
		clusterFunsNames <- names(clusterFuns)
		return(setNames(clusterFunsNames,clusterFunsNames))
	}
)

Ctl.Heatmap$set("public", "RefreshUI_Select_Heatmap_FactorsHeatmap",
	function(session, output, internalVarList, sampleInfo){
		if(is.null(sampleInfo)){
			output$showSelect_Heatmap_FactorsHeatmap = renderText('FALSE')
			internalVarList$showSelect_Heatmap_FactorsHeatmap = FALSE
			outputOptions(output, "showSelect_Heatmap_FactorsHeatmap", suspendWhenHidden = FALSE)
		}else{
			output$showSelect_Heatmap_FactorsHeatmap = renderText('TRUE')
			internalVarList$showSelect_Heatmap_FactorsHeatmap = TRUE
			outputOptions(output, "showSelect_Heatmap_FactorsHeatmap", suspendWhenHidden = FALSE)

			updateSelectInput(
				session, 
				"select_Heatmap_FactorsHeatmap", 
				choices = c(colnames(sampleInfo), "Sample_Name"),
				selected = "Sample_Name"
			)
		}
	}
)


###############################################################################
###################			Result Ouput Functions			###################
###############################################################################

###################			GetMainHeatmap Function			###################
#	Generate the heatmap we can see in heat map tab
#	There're two steps:
#		1. Cut data based on user setting
#		2. Generate heatmap
Ctl.Heatmap$set("public", "GetMainHeatmap",
	function(input, internalVarList, preprocessedResult, preprocessedSampleInfo){
		if(is.null(preprocessedResult)){
			return(NULL)
		}

		if(is.null(preprocessedResult$dat)){
			return(NULL)
		}

		dat <- preprocessedResult$dat
		geneCount <- input$num_Heatmap_IncludeGeneCount
		numHeatmapCutoff <- input$num_Heatmap_HeatmapCutoff
		isGeneCentering <- input$is_Heatmap_GeneCentering
		isGeneNormalize <- input$is_Heatmap_GeneNormalize
		isSampleCentering <- input$is_Heatmap_SampleCentering
		isSampleNormalize <- input$is_Heatmap_SampleNormalize
		
		sampInfo <- preprocessedSampleInfo
		isHaveSelectFactorHeatmap <- internalVarList$showSelect_Heatmap_FactorsHeatmap
		isSampleClustering <- !input$is_Heatmap_NoSampleClustering
		selectFactorsHeatmap <- input$select_Heatmap_FactorsHeatmap
		selectedDistFunction <- input$select_Heatmap_MainPlot_DistanceFun
		selectedhclustFunction <- input$select_Heatmap_MainPlot_HClustFun		
		selectedHeatColor <- input$select_Heatmap_MainPlot_HeatColor

		withProgress(
			message=sample(LogicManager$DB$Quotes,1), 
			detail ="Runing hierarchical clustering ", 
			{
				cuttedData <- LogicManager$Heatmap$CutData(dat, geneCount, numHeatmapCutoff,
					isGeneCentering, isGeneNormalize, isSampleCentering, isSampleNormalize)

				incProgress(1/2, "Generate Plot")

				LogicManager$Heatmap$GenerateHeatmap(cuttedData, sampInfo, geneCount, 
					isHaveSelectFactorHeatmap, isSampleClustering, selectFactorsHeatmap, 
					selectedDistFunction, selectedhclustFunction, selectedHeatColor)

				incProgress(1, "Done")
			}
		)
	}
)

# This function allow user to download the data set used for generate main heat map
# It simply saved get cutted data into given file
Ctl.Heatmap$set("public", "SaveHeatmapDataInFile",
	function( file, input, preprocessedResult){
		# pre check
		if(is.null(preprocessedResult)){
			return(NULL)
		}

		if(is.null(preprocessedResult$dat)){
			return(NULL)
		}

		# get required data
		dat <- preprocessedResult$dat
		geneCount <- input$num_Heatmap_PlotlyIncludeGeneCount
		numHeatmapCutoff <- input$num_Heatmap_HeatmapCutoff
		isGeneCentering <- input$is_Heatmap_GeneCentering
		isGeneNormalize <- input$is_Heatmap_GeneNormalize
		isSampleCentering <- input$is_Heatmap_SampleCentering
		isSampleNormalize <- input$is_Heatmap_SampleNormalize

		# start process
		withProgress(
			message=sample(LogicManager$DB$Quotes,1), 
			{
				incProgress(1/5, "Prepare Data...")
				cuttedData <- LogicManager$Heatmap$CutData(dat, geneCount, numHeatmapCutoff,
					isGeneCentering, isGeneNormalize, isSampleCentering, isSampleNormalize)

				incProgress(4/5, "Prepare File...")
				write.csv(cuttedData, file)
				incProgress(1, "Done")
			}
		)
	}
)

# This function allow user to download a high resolution main heat map use eps format
Ctl.Heatmap$set("public", "SaveMainPlotEpsInTempFile",
	function(file, input, internalVarList, preprocessedResult, preprocessedSampleInfo){
		cairo_ps(file, width = 10, height = 15)
		self$GetMainHeatmap(input, internalVarList, preprocessedResult, preprocessedSampleInfo)
		dev.off()
	}
)

###################		GetMainHeatmapPlotly Function			###################
#	Generate an plotly heatmap that can be interactived with
#	There're two steps:
#		1. Cut data based on user setting
#		2. Cluster Gene (and sample)
#		2. Generate heatmap
Ctl.Heatmap$set("public", "GetMainHeatmapPlotly",
	function(input, internalVarList, preprocessedResult, preprocessedSampleInfo, allGeneInfo){
		if(is.null(preprocessedResult)){
			return(NULL)
		}

		if(is.null(preprocessedResult$dat)){
			return(NULL)
		}

		dat <- preprocessedResult$dat
		geneCount <- input$num_Heatmap_PlotlyIncludeGeneCount
		numHeatmapCutoff <- input$num_Heatmap_HeatmapCutoff
		isGeneCentering <- input$is_Heatmap_GeneCentering
		isGeneNormalize <- input$is_Heatmap_GeneNormalize
		isSampleCentering <- input$is_Heatmap_SampleCentering
		isSampleNormalize <- input$is_Heatmap_SampleNormalize
		
		sampInfo <- preprocessedSampleInfo
		isHaveSelectFactorHeatmap <- internalVarList$showSelect_Heatmap_FactorsHeatmap
		isSampleClustering <- !input$is_Heatmap_NoSampleClustering
		selectFactorsHeatmap <- input$select_Heatmap_FactorsHeatmap
		selectedDistFunction <- input$select_Heatmap_MainPlot_DistanceFun
		selectedhclustFunction <- input$select_Heatmap_MainPlot_HClustFun		
		selectedHeatColor <- input$select_Heatmap_MainPlot_HeatColor

		withProgress(
			message=sample(LogicManager$DB$Quotes,1), 
			detail ="Runing hierarchical clustering ", 
			{
				incProgress(1/5, "Prepare Data...")
				cuttedData <- LogicManager$Heatmap$CutData(dat, geneCount, numHeatmapCutoff,
					isGeneCentering, isGeneNormalize, isSampleCentering, isSampleNormalize)

				incProgress(1/2, "Clustering...")	
				clusteredOrder <- LogicManager$Heatmap$ClusterGeneAndSample(cuttedData, isSampleClustering)

				incProgress(3/5, "Generate Plot...")
				p <- LogicManager$Heatmap$GenerateHeatmapPlotly(cuttedData, clusteredOrder, selectedHeatColor, allGeneInfo)

				incProgress(1, "Done")
				return(p)
			}
		)
	}
)



###################		GetGeneSDHeatmap Functions				###################
#	Generate the plot shows in popup tab: GeneSDHeatmap
Ctl.Heatmap$set("public", "GetGeneSDHeatmap",
	function(input, ConvertedTransformedData){
		withProgress(message="Calculating SD distribution", {
			
			geneCount = input$num_Heatmap_PlotlyIncludeGeneCount
			incProgress(1/5, "Prepare Data...")
			CutResult <- LogicManager$Heatmap$CutData_SD(geneCount, ConvertedTransformedData)
			
			
			incProgress(3/5, "Generate Plot...")
			p <- LogicManager$Heatmap$GenerateSDHeatmapPlot(CutResult$SDs, CutResult$Cutoff)
			
			incProgress(1)
		})
		return(p)
	}
)
#	Download high resolution plot for GeneSDHeatmap
Ctl.Heatmap$set("public", "SaveGeneSDPlotEpsInTempFile",
	function(file, input, ConvertedTransformedData){
		cairo_ps(file, width = 6, height = 4)
		self$GetGeneSDHeatmap(input, ConvertedTransformedData)
		dev.off()
	}
)


###################		GetCorrelationMatrixPlot Functions			###################
#	Generate the matrix plot
Ctl.Heatmap$set("public", "GetCorrelationMatrixPlot",
	function(input, Reactive_PreProcessResult){
		if(is.null(Reactive_PreProcessResult)){
			return(NULL)
		}

		transformedData <- Reactive_PreProcessResult$dat

		if(is.null(transformedData)){
			return(NULL)
		}

		dat <- LogicManager$Heatmap$CutDataForCorrelationMatrixPlot(transformedData)

		p <- LogicManager$Heatmap$GenerateCorrelationPlot(dat, input$isLabelWithPCC)

		return(p)
	}
)


#	Download high resolution plot for CorrelationMatrixPlot
Ctl.Heatmap$set("public", "SaveCorrelationMatrixPlotEpsInTempFile",
	function(file, input, Reactive_PreProcessResult){
		cairo_ps(file, width = 6, height = 4)
		self$GetCorrelationMatrixPlot(input, Reactive_PreProcessResult)
		dev.off()
	}
)

# 	Download CorrelationMatrixPlot data 
Ctl.Heatmap$set("public", "SaveCorrelationMatrixPlotDataInFile",
	function(file, Reactive_PreProcessResult){
				if(is.null(Reactive_PreProcessResult)){
			return(NULL)
		}

		transformedData <- Reactive_PreProcessResult$dat

		if(is.null(transformedData)){
			return(NULL)
		}

		dat <- LogicManager$Heatmap$CutDataForCorrelationMatrixPlot(transformedData)

		write.csv(dat, file)
	}
)

###################		SampleTree Tab related Functions			###################	
#	Generate Sample Tree plot
Ctl.Heatmap$set("public", "GetSampleTreePlot",
	function(input, Reactive_PreProcessResult){

		# prepare data
		transformedData <- Reactive_PreProcessResult$dat
		isGeneCentering <- input$is_Heatmap_GeneCentering
		isGeneNormalize <- input$is_Heatmap_GeneNormalize
		isSampleCentering <- input$is_Heatmap_SampleCentering
		isSampleNormalize <- input$is_Heatmap_SampleNormalize
		
		dat <- LogicManager$Heatmap$CutDataForSampleTreePlot(transformedData, isGeneCentering, isGeneNormalize, isSampleCentering, isSampleNormalize)
		

		# generate plot
		selectedDistFunction <- input$select_Heatmap_MainPlot_DistanceFun
		selectedhclustFunction <- input$select_Heatmap_MainPlot_HClustFun		
		
		p <- LogicManager$Heatmap$GenerateSampLeTreePlot(dat, selectedhclustFunction, selectedDistFunction)


		# return plot

		return(p)
	}
)

# download eps plot of sample tree
Ctl.Heatmap$set("public", "SaveSampleTreePlotEpsInTempFile",
	function(file, input, Reactive_PreProcessResult){
		cairo_ps(file, width = 8, height = 6)
		self$GetSampleTreePlot(input, Reactive_PreProcessResult)
		dev.off()
	}
)

