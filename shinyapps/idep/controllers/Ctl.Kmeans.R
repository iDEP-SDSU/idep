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
	function(){

	}
)





