library('R6')
library(shiny)
library(shinyBS)
library(plotly)

Ctl.PreProcess <- R6Class("Ctl.PreProcess")

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