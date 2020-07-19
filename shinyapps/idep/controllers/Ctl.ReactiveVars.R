library('R6')
library(shiny)
library(shinyBS)
library(plotly)
library(processx)


Ctl.ReactiveVars <- R6Class("Ctl.ReactiveVars")



## ReactVars$PreProcessResult
Ctl.ReactiveVars$set("public", "PreProcessResult",
	function(input, session, storeVariableList){
		withProgress(message="Reading and pre-processing", {

			# Pre checking
			if(is.null(storeVariableList$RawData)){
				#if no raw data, then return null
				return(NULL)
			}



			if(input$dataFileFormat == 1){
				incProgress(1/3, "Pre-processing counts data")
				minCount = input$numMinCounts
				minSample = input$numNMinSamplesInCountCase
				logStart = input$numCountsLogStart
				incProgress(1/2,"transforming raw counts")
			}else{
				# actrually should be input.dataFileFormat == 2
				# However, when input.dataFileFormat == 3, minCound and minSample won't use.
				# So we simply it into 'if else' statement.
				incProgress(1/3,"Pre-processing data")
				minCount = input$numMinFPKM
				minSample = input$numNMinSampleInFPKMCase
				logStart = input$numFPKMLogStart
			}

			#if preprocess is not trigger, do it.
			PreProcessResult <- 
				LogicManager$PreProcessing$RawDataPreprocess(
					storeVariableList$RawData, 
					input$selectMissingValueImputationMethod,
					input$dataFileFormat,
					minCount,
					minSample,
					as.numeric(input$selectCountsTransform),
					logStart,
					input$isApplyLogTransFPKM,
					input$isNoFDR			## This parm is got in load data tab.
				)
			incProgress(1, "Done.")
			return(PreProcessResult)
		})
	}
)






# the gene set used in PCA page
Ctl.ReactiveVars$set("public", "GeneSetPCA", 
	function(input, ConvertedIDResult, ConvertedTransformedData){
		selectOrg <- input$selectOrg
		gmtFile <- input$gmtFile
		GO <- input$select_PCA_GO 
        minGeneSetSize <- input$num_PCA_minGeneSetSize ##this is in pathway ???
        maxGeneSetSize <- input$num_PCA_maxGeneSetSize ##this is in pathway ???
		return(
			DB.Manager$QueryGeneSetsFromPathway(
				ConvertedIDResult,
                ConvertedTransformedData,
                GO,
                selectedOrg,
                c(minGeneSetSize, maxGeneSetSize)
			)
		)

	}
)

# the PCAdata set downloaded in PCA page
Ctl.ReactiveVars$set("public", "PCAdata", 
	function(input, PreProcessResult){
        if (is.null(input$file1)&& input$goButton == 0)   return(NULL)

        PreprocessResultData <- PreProcessResult$data
        TSNESeed <- input$btn_PCA_RecalcTSNE

        return(
            LogicManager$PCA$GeneratePCAData(PreprocessResultData, TSNESeed)
        )
        
	}
)



# exactly same logic as "GetGeneSetPCA"
# but the 'select GO is input$selectGO' in original code
Ctl.ReactiveVars$set("public", "GeneSets", 
	function(input, ConvertedIDResult, ConvertedTransformedData){
		selectOrg <- input$selectOrg
		gmtFile <- input$gmtFile
		GO <- input$selectGO
        minGeneSetSize <- input$minGeneSetSize 
        maxGeneSetSize <- input$maxGeneSetSize 
		return(
			LogicManager$Pathway$GetGeneSetByGOOption(
				ConvertedIDResult,
                ConvertedTransformedData,
                GO,
                selected
			)
		)
	}
)


# Calculate Kmeans() reactive variable
Ctl.ReactiveVars$set("public", "KmeansReactiveVar",
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
Ctl.ReactiveVars$set("public", "KmeansWithGeneInfo",
	function(input, Reactive_Kmeans, Reactive_AllGeneInfo){
		
		x <- Reactive_Kmeans$x
		bar <- Reactive_Kmeans$bar
		allGeneInfo <- Reactive_AllGeneInfo
		selectedOrg <- input$selectOrg

		return(LogicManager$Kmeans$MergeGenInfoWithClusterResult( x, bar, allGeneInfo, selectedOrg))
	}
)



# Caculate Kmeans GO data
Ctl.ReactiveVars$set("public", "KmeansGoData",
    function(input, Reactive_Kmeans, Reactive_ConvertedIDResult, Reactive_AllGeneInfo, Reactive_GeneSets){
        minFDR = 0.01
        selectedOrg <- input$selectOrg
        gmtFile <- input$gmtFile
        nCluster <- input$num_Kmeans_Culsters
        GO <- input$select_Kmeans_PathwayDatabase
        is_Kmeans_RemoveRedudantSets <- input$is_Kmeans_RemoveRedudantSets
        return(
            LogicManager$Kmeans$GetKmeansGoData(
                minFDR, selectedOrg, gmtFile, nCluster, GO, is_Kmeans_RemoveRedudantSets,
                Reactive_Kmeans, Reactive_ConvertedIDResult, Reactive_AllGeneInfo,
                Reactive_GeneSets
            )
        ) 
    }
)

