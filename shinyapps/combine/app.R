# Extract and combine files
library(shiny)
# Define UI for application that draws a histogram
ui <- fluidPage(
   titlePanel("Extract read-count or FPKM from multiple files"),
   sidebarLayout(
      sidebarPanel(
        h5("This app extracts a certain column from multiple files with the same structure. It can be used to extract read-count or FPKM/RPKM data from multiple files
           and combine into one file using gene ID. Only support CSV files. Gene IDs can be of different order. Files can have different number of rows/genes. ")
        , numericInput("colID",
                     "Column index for read-count or FPKM to be extracted:",
                     value = 2, min = 1, max=100)
        , numericInput("geneIDcol",
                       "Column index for gene IDs",
                       value = 1, min = 1, max=100)
      ,checkboxInput("hasHeader", "Has column header", value = TRUE)

       , fileInput('files', 'Upload multiple CSV files',
                   accept = c(
                     'text/csv',
                     '.csv'
                   )
                  ,multiple = TRUE
        )

        ),

      mainPanel(
        downloadButton('download.File',"Download combined file")
        ,tableOutput("data1")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    combinedData <- reactive({
     if (is.null(input$files))   return(NULL)

     ci = input$colID-1
     n = length(input$files[,1]) # number of files

     withProgress(message="combining files by gene IDs ", {

     #read the first file
     data = read.csv(input$files[1, 'datapath']
                     ,header=input$hasHeader
                     ,row.names=input$geneIDcol)
     colnames(data)[ci] = input$files$name[1]
     data <- data[,ci,drop=F]
      incProgress(1/n, message= input$files$name[1] )

     #adding files by merging with gene ID
     for(i in 2:n) {
       new1 = read.csv(input$files[i, 'datapath']
                       ,header=input$hasHeader
                       , row.names=input$geneIDcol)
       colnames(new1)[ci] = input$files$name[i]
       new1 <- new1[,ci,drop=F]
       data <- merge(data,new1, by="row.names", all=TRUE)
       row.names(data) <- data[,1]
       data <- data[,-1]

       incProgress(i/n, message= input$files$name[i] )
     }

     # missing values as 0
     data [is.na(data)] <- 0
     colnames(data) = gsub(".*/|\\..*","",colnames(data))

     })

     return(data)
   })

   output$data1 <- renderTable({
     if (is.null(input$files))   return(NULL)
     isolate({
        combinedData()[1:30,]
     })
   },include.rownames=TRUE,
      striped=TRUE,
      bordered = TRUE,
      width = "auto",hover=T)

   output$download.File <- downloadHandler(
     filename = function() {"combined_data.csv"},
     content = function(file) {
       write.csv(combinedData(), file)
     }
   )
}

# Run the application
shinyApp(ui = ui, server = server)
