library('R6')

Ctl.Kmeans <- R6Class("Ctl.Kmeans")

###############################################################################
###################			Load/Reload UI Functions		###################
###############################################################################

Ctl.Kmeans$set("public", "InitPathwayDatabaseSelection",
	function(){
		
	}
)


###############################################################################
###################			Reactive Variables      		###################
###############################################################################

# Calculate Kmeans() reactive variable
Ctl.Kmeans$set("public", "GetKmeansReactiveVar",
    function(input, Reactive_ConvertedTransformedData ){
        withProgress(message="Converting data ... ", {            
            # if no converted transformed data, return null
            if(is.null(Reactive_ConvertedTransformedData)){
                return(NULL)
            }

            GeneCount <- input$num_Kmeans_GenesKNN
            NormalizationMethod <- input$select_Kmeans_Normalization
            RerunSeed <- input$btn_Kmeans_rerun
            NumberOfCluster <- input$num_Kmeans_Culsters

            incProgress(0.3, detail = paste("Calc Kmeans ... "))
            
            result <- LogicManager$Kmeans$CalcKmeansCluster(
                Reactive_ConvertedTransformedData, 
                GeneCount, 
                NormalizationMethod, 
                RerunSeed,
                NumberOfCluster
            )
            
            incProgress(1, detail = paste("Done"))
            
            return(result)
        })
    }
)

# Calculate KmeansDataWithGeneInfo() reactive variable. 
Ctl.Kmeans$set("public", "GetKmeansWithGeneInfo",
	function(input, Reactive_Kmeans, Reactive_AllGeneInfo){
		
		x <- Reactive_Kmeans$x
		bar <- Reactive_Kmeans$bar
		allGeneInfo <- Reactive_AllGeneInfo
		selectedOrg <- input$selectOrg

		return(LogicManager$Kmeans$MergeGenInfoWithClusterResult( x, bar, allGeneInfo, selectedOrg))
	}
)


###############################################################################
###################			Result Ouput Functions			###################
###############################################################################

Ctl.Kmeans$set("public", "GetMainHeatmapPlot",
	function(input, Reactive_Kmeans){
		if( is.null(Reactive_Kmeans) ){
			return(NULL)
		}

		withProgress(message="Creating heatmap", {
			dat <- Reactive_Kmeans$x
			bar <- Reactive_Kmeans$bar
			colorOption <- input$select_Heatmap_MainPlot_HeatColor

			incProgress(1/2, "Generate Plot")

			LogicManager$Kmeans$GetKmeansClusterHeatmapWithGeneBar(
				dat, bar, colorOption
			)
			incProgress(1, "Done")
		})
	}
)


## Why this function using different logic 
Ctl.Kmeans$set("public", "GetNclusterPlot", 
	function(input, Reactive_ConvertedTransformedData ){
        withProgress(message="Converting data ... ", {
            GeneCount <- input$num_Kmeans_GenesKNN
            NormalizationMethod <- input$select_Kmeans_Normalization
            
            incProgress(0.3, detail = paste("Calc Kmeans ... "))

            wss_data <- LogicManager$Kmeans$CalcKmeansNCluster(
                Reactive_ConvertedTransformedData,
                GeneCount,
                NormalizationMethod
            )

            incProgress(0.6, detail = paste("Generate plot ... "))

            LogicManager$Kmeans$PlotWithinGroupSquareSum(wss_data)

            incProgress(1, "Done")
        }
	}
)

# GetTSNEGenePlot 
# The actual business logic will calculate tsne and then plot the result
Ctl.Kmeans$set("public", "GetTSNEGenePlot", 
	function(input, Reactive_Kmeans){
		withProgress(message="Runing t-SNE algorithm", {
			isolate({
				Cluster <- Reactive_Kmeans$bar
				train <- as.data.frame( cbind(Cluster, Reactive_Kmeans$x) )

				train = unique(train)
				Cluster = train$Cluster	
				seed <- input$btn_Kmeans_Recalculate
				colorGenes <- input$colorGenes

				incProgress(1/3, "Calculate t-SNE and generate plot")
				
				LogicManager$Kmeans$CalculateTSNEAndGeneratePlot( train, Cluster, seed, colorGenes )
				
				incProgress(1, "Done")
			})
		})
	}
)


#	Generate the plot shows in popup tab: GeneSDHeatmap
Ctl.Kmeans$set("public", "GetGeneSDHeatmap",
	function(input, ConvertedTransformedData){
		withProgress(message="Calculating SD distribution", {
			
			geneCount = input$num_Kmeans_GenesKNN
			incProgress(1/5, "Prepare Data...")
			CutResult <- LogicManager$Heatmap$CutData_SD(geneCount, ConvertedTransformedData) ## Yes. Kmeans one use exactly same logic and code as Heatmap one
			
			
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


# download eps plot of main heat map
Ctl.Kmeans$set("public", "SaveMainHeatmapPlotEpsInTempFile",
	function(file, input, Reactive_Kmeans){
		cairo_ps(file, width = 8, height = 6)
		self$GetMainHeatmapPlot(input, Reactive_Kmeans)
		dev.off()
	}
)



