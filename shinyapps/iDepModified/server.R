library(shiny)


source('server.config')
source('controllers/ControllerManager.R')

ReactVars <- reactiveValues()

LoadDataCtrl <- Ctl.LoadData$new()

shinyServer(server)

server <- function(input, output, session) {
    
    ################################################################
    #   Read data
    ################################################################

    output$ViewData_SelectFromPublic_Samples <- LoadDataCtrl$RenderSampleTable(input) 
    output$SearchData <- LoadDataCtrl$SearchGSEIDs(input)
    output$humanNsamplesOutput <- LoadDataCtrl$RenderHumanNsampleOutput(input)
    output$mouseNsamplesOutput <- LoadDataCtrl$RenderMouseNsampleOutput(input)
    output$selectedDataset <- LoadDataCtrl$RenderSelectedDataset(input)
    output$DoneLoading <- LoadDataCtrl$RenderInitDoneUI()

    observeEvent(input$btn_LoadData_UseSelectedPublicData,{
        LoadDataCtrl$EventHandler_UseSelectedPublicData(input, output, session)
    })

    observeEvent(input$LoadData_uploadedDataFile,{
        LoadDataCtrl$EventHandler_UploaedDataFileChanged(input, output, session)
    })
}
