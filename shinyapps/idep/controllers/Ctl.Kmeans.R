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

###############################################################################
###################			Result Ouput Functions			###################
###############################################################################

Ctl.Kmeans$set("public", "RenderMainHeatmapPlot",
	function(Reactive_Kmeans, input){
		if( is.null(Reactive_Kmeans) ){
			return(NULL)
		}

		withProgress(message="Creating heatmap", {
			dat <- Reactive_Kmeans$SortedData
			bar <- Reactive_Kmeans$SortedIndex
			colorOption <- input$select_Heatmap_MainPlot_HeatColor

			incProgress(1/2, "Generate Plot")

			LogicManager$Kmeans$GetKmeansClusterHeatmapWithGeneBar(
				dat, bar, colorOption
			)
			incProgress(1, "Done")
		})
	}
)

# download eps plot of main heat map
Ctl.Heatmap$set("public", "SaveMainHeatmapPlotEpsInTempFile",
	function(file, input, Reactive_Kmeans){
		cairo_ps(file, width = 8, height = 6)
		self$RenderMainHeatmapPlot(Reactive_Kmeans, input)
		dev.off()
	}
)



