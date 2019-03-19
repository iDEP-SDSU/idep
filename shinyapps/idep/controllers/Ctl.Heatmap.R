library('R6')

Ctl.Heatmap <- R6Class("Ctl.Heatmap")


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

Ctl.Heatmap$set("public", "GetMainPlot",
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

Ctl.Heatmap$set("public", "SaveHeatmapDataInTempFile",
	function(){
		#
	}
)

Ctl.Heatmap$set("public", "SaveEpsPlotInTempFile",
	function(file, input, internalVarList, preprocessedResult, preprocessedSampleInfo){
		cairo_ps(file, width = 10, height = 15)
		self$GetMainPlot(input, internalVarList, preprocessedResult, preprocessedSampleInfo)
		dev.off()
	}
)