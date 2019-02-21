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

Ctl.PreProcess$set("public", "getTotalReadCountsData",
	function(input, output, session, storeVariableList ){
		
		if(is.null(storeVariableList$RawData)){
			#if no raw data, then return null
			return(NULL)
		}

		if(storeVariableList$DataSourceType != 1){
			#if the data is not read count, return null
			return(NULL)
		}

		if(is.null(storeVariableList$PreProcessResult)){
			#if preprocess is not trigger, do it.
			storeVariableList$PreProcessResult <- LogicManager$PreProcessing$RawDataPreprocess()  ### parm not right yet. 
		}

		
	}
)