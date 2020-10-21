library('R6')

Ctl.Kmeans <- R6Class("Ctl.Kmeans")

###############################################################################
###################			Load/Reload UI Functions		###################
###############################################################################

Ctl.Kmeans$set("public", "InitPathwayDatabaseSelection",
    function(input, ConvertedIDResult, ConvertedTransformedData){
        renderUI({
            tem = input$selectOrg
            if (is.null(input$fileUploadedData)&& input$btn_LoadData_DemoData == 0 ){ 
                selectInput(
                    "select_Kmeans_PathwayDatabase", 
                    label = NULL,
                    choices = list("All available gene sets" = "All", "GO Biological Process" = "GOBP","GO Molecular Function" = "GOMF","GO Cellular Component" = "GOCC",
                                        "KEGG metabolic pathways" = "KEGG"), selected = "GOBP")  
            }else { 
                selectInput(
                    "select_Kmeans_PathwayDatabase", 
                    label=NULL,
                    choices=LogicManager$DB$gmtCategory(ConvertedIDResult, ConvertedTransformedData, input$selectOrg,input$gmtFile),
                    selected = "GOBP" 
                )   
            } 
        })
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
        })
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
Ctl.Kmeans$set("public", "SaveGeneSDPlotEpsInTempFile",
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

# the table displaying Kmeans GO data
Ctl.Kmeans$set("public", "GetKmeansGoTableData",
	function(input, Reactive_KmeansGOData){
        is_Kmeans_RemoveRedudantSets = input$is_Kmeans_RemoveRedudantSets

		return(
            LogicManager$Kmeans$CalculateKmeansGoTableData(Reactive_KmeansGOData, is_Kmeans_RemoveRedudantSets)
        )
	}
)

# Enrichment plot
Ctl.Kmeans$set("public", "GetEnrichmentPlot",
    function(input, Reactive_Kmeans, Reactive_KmeansGOData){

        selectedGO <- input$selectGO3
        selectedOrg <- input$selectOrg
        gmtFile <- input$gmtFile

       	##################################  
        # these are needed to make it responsive to changes in parameters
        tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
        if( !is.null(input$dataFileFormat) ) 
            if(input$dataFileFormat== 1)  
                {  tem = input$minCounts ; tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
        if( !is.null(input$dataFileFormat) )
            if(input$dataFileFormat== 2) 
                { tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
        tem = input$KmeansReRun; 
        tem = input$nGenesKNN;
        tem = input$kmeansNormalization
        tem = input$nClusters
        tem = input$removeRedudantSets
        ####################################
        

        return( LogicManager$Kmeans$GenerateEnrichmentPlot(selectedGO, selectedOrg, gmtFile, Reactive_Kmeans, Reactive_KmeansGOData ))
    }
)

# Kmeans Promoter Table
Ctl.Kmeans$set("public", "GetPromoterTable",
    function(input, Reactive_Kmeans){

        nClusters <- input$nClusters
        selectOrg <- input$selectOrg
        selectGO2 <- input$selectGO2
        promoterBP <- input$select_PromoterKmeans

        ##################################  
        # these are needed to make it responsive to changes in parameters
        tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
        if( !is.null(input$dataFileFormat) ) 
            if(input$dataFileFormat== 1)  
                {  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
        if( !is.null(input$dataFileFormat) )
            if(input$dataFileFormat== 2) 
                { tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
        tem = input$KmeansReRun
        ####################################

        withProgress(message="Promoter analysis", {
            result <-  LogicManager$Kmeans$PromoterAnalysis(
                nClusters,
                selectOrg,
                selectGO2,
                promoterBP,
                Reactive_Kmeans
            )
            incProgress(1, detail = paste("Done")) 
	    }) #progress

    }

)
