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










Ctl.ReactiveVars$set("public", "GetGeneSet", 
	function(input, ConvertedIDResult, ConvertedTransformedData){

		selectOrg <- input$selectOrg
		gmtFile <- input$gmtFile
		GO <- input$ ## the select GO 1 in pathway analysis

		return(
			LogicManager$Pathway$GetGeneSetByGOOption(
				ConvertedIDResult
			)
		)
	}
)

Ctl.ReactiveVars$set("public", "GetGeneSetPCA", 
	function(input, ConvertedIDResult, ConvertedTransformedData){
		### this one need use GO 6
		
		selectOrg <- input$selectOrg
		gmtFile <- input$gmtFile
		GO <- input$selectGO6

	}
)


# Caculate Kmeans GO data
Ctl.ReactiveVars$set("public", "GetKmeansGoData",
    function(input, Reactive_Kmeans, Reactive_ConvertedIDResult, Reactive_AllGeneInfo){
        stop('this function need confirm with Dr. Ge')
        withProgress(message=sample(quotes,1), detail ="GO Enrichment", {

            pp=0
            minFDR = 0.01
            selectedOrg <- input$selectOrg
            gmtFile <- input$gmtFile
            nCluster <- input$num_Kmeans_Culsters
			GO <- input$select_Kmeans_PathwayDatabase

			overlap <- LogicManager$DB$findOverlap(query, )			

            return(LogicManager$Kmeans$) 
        })
    }
)

