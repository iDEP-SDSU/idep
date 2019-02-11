





server <- function(input, output) {
	
    output$ViewData_SelectFromPublic_Samples <- LoadDataCtrl$RenderSampleTable(input) 
    output$SearchData <- LoadDataCtrl$SearchGSEIDs(input)
    output$humanNsamplesOutput <- LoadDataCtrl$RenderHumanNsampleOutput(input)
    output$mouseNsamplesOutput <- LoadDataCtrl$RenderMouseNsampleOutput(input)
    output$selectedDataset <- LoadDataCtrl$RenderSelectedDataset(input)
    output$DoneLoading <- LoadDataCtrl$RenderInitDoneUI()

}


shinyApp(ui, server)