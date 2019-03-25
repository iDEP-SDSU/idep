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


Ctl.Heatmap$set("public", "SaveHeatmapDataInFile",
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
		withProgress(
			message=sample(LogicManager$DB$Quotes,1), 
			{
				incProgress(1/5, "Prepare Data...")
				cuttedData <- LogicManager$Heatmap$CutData(dat, geneCount, numHeatmapCutoff,
					isGeneCentering, isGeneNormalize, isSampleCentering, isSampleNormalize)
				incProgress(4/5, "Prepare File...")
				wrtie.csv(cuttedData, file)
				incProgress(1, "Done")
			}
		)
	}
)

Ctl.Heatmap$set("public", "SaveEpsPlotInTempFile",
	function(file, input, internalVarList, preprocessedResult, preprocessedSampleInfo){
		cairo_ps(file, width = 10, height = 15)
		self$GetMainPlot(input, internalVarList, preprocessedResult, preprocessedSampleInfo)
		dev.off()
	}
)

Ctl.Heatmap$set("public", "GetGeneSDHeatmap",
	function(input, ConvertedTransformedData){
		withProgress(message="Calculating SD distribution", {
			
			geneCount = input$num_Heatmap_PlotlyIncludeGeneCount

			incProgress(1/5, "Prepare Data...")
			SDs <- CutData_SD(geneCount, ConvertedTransformedData)
			Cutoff=sort(SDs,decreasing=TRUE)[geneCount]

			incProgress(3/5, "Generate Plot...")
			p <- GetGeneSDHeatmap(input, ConvertedTransformedData)
			incProgress(1)
		})
		return(p)
	}
)

