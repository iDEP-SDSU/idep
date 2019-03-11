library('R6')
library(shiny)
library(shinyBS)
library(plotly)
source('businessLogic/LogicManager.R')

Ctl.PreProcess <- R6Class("Ctl.PreProcess")
LogicManager <- Logic.Manager$new()

Ctl.PreProcess$set("public", "PreProcessResult",
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
		})
		
		return(PreProcessResult)
	}
)

# Preprocessing tab, first plot
Ctl.PreProcess$set("public", "GetTotalReadCountsPlot",
	function(rawCount){
		# check data 
		if(is.null(rawCount)){
			return(NULL)
		}

		# calculate group
		groups = as.factor( LogicManager$PreProcessing$DetectGroups(colnames(rawCount)) )

		# fetch the plot
		return(LogicManager$Display$GetReadcountBarPlot(rawCount, groups))
	}
)


# Preprocessing tab, second plot
Ctl.PreProcess$set("public", "GetTransformedDataBoxPlot",
	function(transformedData){
		#
		# this function should be called after LogicManager$PreProcessing$RawDataPreprocess() function has been called. 
		#

		if(is.null(transformedData)){
			# if we cannot get result, then exit with null
			return(NULL)
		} 

		# Get transfered data
		transformedData <- as.data.frame(transformedData)

		# calculate group
		groups = as.factor( LogicManager$PreProcessing$DetectGroups(colnames(transformedData)) )

		# fetch the plot
		return(LogicManager$Display$GetTransformedDataBoxPlot(transformedData, groups))
	}
)

# Preprocessing tab, third plot
Ctl.PreProcess$set("public", "GetTransformedDataDensityPlot",
	function(transformedData){
		#
		# this function should be called after LogicManager$PreProcessing$RawDataPreprocess() function has been called. 
		#

		if(is.null(transformedData)){
			# if we cannot get result, then exit with null
			return(NULL)
		} 

		# Get transfered data
		transformedData <- as.data.frame(transformedData)

		# calculate group
		groups = as.factor( LogicManager$PreProcessing$DetectGroups(colnames(transformedData)) )

		# fetch the plot
		return(LogicManager$Display$GetTransformedDataDensityPlot(transformedData, groups))
	}
)

# Preprocessing tab, fourth plot
Ctl.PreProcess$set("public", "GetTransformedDataScatterPlot",
	function(transformedData){
		#
		# this function should be called after LogicManager$PreProcessing$RawDataPreprocess() function has been called. 
		#

		if(is.null(transformedData)){
			# if we cannot get result, then exit with null
			return(NULL)
		} 

		# Get transfered data
		transformedData <- as.data.frame(transformedData)

		# fetch the plot
		return(LogicManager$Display$GetTransformedDataScatterPlot(transformedData))
	}
)

Ctl.PreProcess$set("public", "GetGuessSpeciesResult",
	function(convertedIDResult){
		if(is.null(convertedIDResult)){
			return(as.data.frame("ID not recognized."))
		} 
		else { 
			return(convertedIDResult$speciesMatched)
		}
	}
)


Ctl.PreProcess$set("public", "ConvertedIDResult",
	function(transformedData, selectedOrg){
		withProgress(
			message = something we didnt finished(), 
			detail="Converting gene IDs", 
			{
				Converted code here
				incProgress(1, detail = paste("Done"))
			}
		)
	}

)