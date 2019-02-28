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

Ctl.PreProcess$set("public", "testPlot",
	function(){
		p <- plot(rnorm(100))
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


# Preprocessing tab, first plot
Ctl.PreProcess$set("public", "getTotalReadCountsData",
	function(input, output, session, storeVariableList){
		
		# Pre checking
		if(is.null(storeVariableList$RawData)){
			#if no raw data, then return null
			return(NULL)
		}

		if(storeVariableList$DataSourceType != 1){
			#if the data is not read count, return null
			return(NULL)
		}

		if( !storeVariableList$PreProcessDone ){
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
			storeVariableList$PreProcessResult <- 
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
			storeVariableList$PreProcessDone <- TRUE
			
			if(is.null(storeVariableList$PreProcessResult)){
				# if we cannot get result, then exit with null
				return(NULL)
			} 
		}

		# Get Read Count
		PreProcessedReadCount <- storeVariableList$PreProcessResult$rawCount

		# calculate group
		groups = as.factor( LogicManager$PreProcessing$DetectGroups(colnames(PreProcessedReadCount)) )

		# fetch the plot
		return(LogicManager$Display$GetReadcountBarPlot(PreProcessedReadCount, groups))
	}
)


# Preprocessing tab, second plot
Ctl.PreProcess$set("public", "getTransDataBoxPlot",
	function(input, output, session, storeVariableList){
		#
		# this function should be called after LogicManager$PreProcessing$RawDataPreprocess() function has been called. 
		#

		if(is.null(storeVariableList$PreProcessResult)){
			# if we cannot get result, then exit with null
			return(NULL)
		} 

		# Get transfered data
		TransferedData <- as.data.frame(storeVariableList$PreProcessResult$dat)

		# calculate group
		groups = as.factor( LogicManager$PreProcessing$DetectGroups(colnames(TransferedData)) )

		# fetch the plot
		return(LogicManager$Display$GetTransedDataBoxPlot(TransferedData, groups))
	}
)

# Preprocessing tab, third plot
Ctl.PreProcess$set("public", "getTransDataDensityPlot",
	function(input, output, session, storeVariableList){
		#
		# this function should be called after LogicManager$PreProcessing$RawDataPreprocess() function has been called. 
		#

		if(is.null(storeVariableList$PreProcessResult)){
			# if we cannot get result, then exit with null
			return(NULL)
		} 

		# Get transfered data
		TransferedData <- as.data.frame(storeVariableList$PreProcessResult$dat)

		# calculate group
		groups = as.factor( LogicManager$PreProcessing$DetectGroups(colnames(TransferedData)) )

		# fetch the plot
		return(LogicManager$Display$GetTransedDataDensityPlot(TransferedData, groups))
	}
)

# Preprocessing tab, fourth plot
Ctl.PreProcess$set("public", "getTransDataScatterPlot",
	function(input, output, session, storeVariableList){
		#
		# this function should be called after LogicManager$PreProcessing$RawDataPreprocess() function has been called. 
		#

		if(is.null(storeVariableList$PreProcessResult)){
			# if we cannot get result, then exit with null
			return(NULL)
		} 

		# Get transfered data
		TransferedData <- as.data.frame(storeVariableList$PreProcessResult$dat)

		# fetch the plot
		return(LogicManager$Display$GetTransedDataScatterPlot(TransferedData))
	}
)
