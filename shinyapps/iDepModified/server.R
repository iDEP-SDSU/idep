server <- function(input, output) {
	observeEvent(input$btn_LoadData_SearchDataFromPublic,{
		output$showOther <- renderText(1)
		outputOptions(output, "showOther", suspendWhenHidden = FALSE)
	})
	
  output$distPlot <- renderPlot({
    hist(rnorm(input$obs))
  })
}


shinyApp(ui, server)