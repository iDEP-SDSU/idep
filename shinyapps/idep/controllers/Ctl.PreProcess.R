library('R6')
library(shiny)
library(shinyBS)
library(plotly)
source('businessLogic/LogicManager.R')

Ctl.PreProcess <- R6Class("Ctl.PreProcess")
LogicManager <- Logic.Manager$new()

Ctl.PreProcess$set("public","testPlotly",
	function(){
		x <- c(1:100)
		random_y <- rnorm(100, mean = 0)
		data <- data.frame(x, random_y)
		
		p <- plot_ly(data, x = ~x, y = ~random_y, type = 'scatter', mode = 'lines')
		return(p)
	}
)



Ctl.PreProcess$set("public","testPlotlyData",
	function(){
		x <- c(1:100)
		random_y <- rnorm(100, mean = 0)
		data <- data.frame(x, random_y)
		
		
		return(data)
	}
)


Ctl.PreProcess$set("public", "testData",
	function(input, session, storeVariableList){
	
		# Pre checking
		if(is.null(storeVariableList$RawData)){
			#if no raw data, then return null
			return(NULL)
		}



			if(storeVariableList$DataSourceType == 1){
				minCount = input$numMinCounts
				minSample = input$numNMinSamplesInCountCase
				logStart = input$numCountsLogStart
			}else{
				# actrually should be output.DataSourceType == 2
				# However, when output.DataSourceType == 3, minCound and minSample won't use.
				# So we simply it into 'if else' statement.

				minCount = input$numMinFPKM
				minSample = input$numNMinSampleInFPKMCase
				logStart = input$numFPKMLogStart
			}
	
			#if preprocess is not trigger, do it.
			PreProcessResult <- 
				LogicManager$PreProcessing$RawDataPreprocess(
					storeVariableList$RawData, 
					input$selectMissingValueImputationMethod,
					storeVariableList$DataSourceType,
					minCount,
					minSample,
					as.numeric(input$selectCountsTransform),
					logStart,
					input$isApplyLogTransFPKM,
					input$isNoFDR			## This parm is got in load data tab.
				)
			
			return(PreProcessResult)
	}
)

Ctl.PreProcess$set("public", "testPlot",
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
# Preprocessing tab, first plot
Ctl.PreProcess$set("public", "getTotalReadCountsData",
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
Ctl.PreProcess$set("public", "getTransformedDataBoxPlot",
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
Ctl.PreProcess$set("public", "getTransformedDataDensityPlot",
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
Ctl.PreProcess$set("public", "getTransformedDataScatterPlot",
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
