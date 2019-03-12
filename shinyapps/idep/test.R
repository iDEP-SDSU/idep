library(shiny)
library(plotly)

ui <- fluidPage(
	plotlyOutput("plot"),
	verbatimTextOutput("event")
)

server <- function(input, output) {

	# renderPlotly() also understands ggplot2 objects!
	output$plot <- renderPlotly({
		plot_ly(mtcars, x = ~mpg, y = ~wt)
	})

	output$event <- renderPrint({
		d <- event_data("plotly_hover")
		if (is.null(d)) "Hover on a point!" else d
	})
}

shinyApp(ui, server)
