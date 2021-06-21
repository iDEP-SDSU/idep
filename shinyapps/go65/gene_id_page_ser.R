####################################################
# Author: Eric Tulowetzke, eric.tulowetzke@jacks.sdstate.edu
# Lab: Ge Lab
# R version 4.0.5
# Project: iDEP v93
# File: gene_id_page_ser.R
# Purpose of file:server logic for second tab i.e. Gene ID Examples
#  Allow users view example of database
# Start data: 06-06-2022 (mm-dd-yyyy)
# Data last modified: 06-06-2021, 12:46 PM CST (mm-dd-yyyy,TIME) 
# to help with github merge 
#######################################################

firstTime <- TRUE #Initialize the page at start up in sever logic
MAX_WIDTH_COL <- 150 #Determine the starting with of table columns
EXAMPLE_GENES_WCOL <- 800

#################################################################
# FUNCTION : getExample 
# DESCRIPTION : Gives user example of our IDs, 
# INPUT ARGS : Organism picked by user, user's database option,
#  and path to database
# OUTPUT ARGS : data frame of genes
# IN/OUT ARGS :
# RETURN : returnDf
# Implementation notes : shiny is only used when this is used in shiny app
#################################################################
getExample <- function(userSpecie = NULL, userIDtype = NULL,
                       path2Database = NULL, shiny = FALSE) {
  returnDf = NULL
  convert <- dbConnect(RSQLite::SQLite(), path2Database)
  query4IDtype <- paste('SELECT * FROM idIndex WHERE idType =',
                        shQuote(userIDtype))
  userIDtypeNum <- dbGetQuery(convert, query4IDtype)
  query4IDmap <- paste('SELECT * FROM mapping WHERE species =',
                       as.numeric(userSpecie),
                       'AND idType =', as.numeric(userIDtypeNum$id))
  userIDdf <- dbGetQuery(convert, query4IDmap)
  RSQLite::dbDisconnect(convert)
  returnDf <- userIDdf[-c(3,4)]
  colnames(returnDf) <- c('User_ID', 'Ensembl_ID')
  if (shiny) {incProgress(1)}
  return(returnDf)
} # end of function


#################################################################
# FUNCTION : getExampleDfID
# DESCRIPTION : Gives user example of our IDs, 
# INPUT ARGS : Organism picked by user
# OUTPUT ARGS : data frame of genes
# IN/OUT ARGS :
# RETURN : returnDf
# Implementation notes : shiny is only used when this is used in shiny app
#################################################################
getExampleDfID <- function(userSpecie = NULL, path2Database = NULL,
                           shiny = FALSE) {
  allSample <- feather::read_feather(path2Database)
  returnDF <- allSample[allSample$index == userSpecie,]
  returnDF <- returnDF[, -c(1)]
  colnames(returnDF) <- c('Database', 'Example_genes')
  if (shiny) {incProgress(1)}
  return(returnDF)
}# end of getExampleDfID


############################################
#Purpose: server logic for second tab i.e. Gene ID Examples
############################################
geneIDPage <- function(input, output, session, orgInfo, path) {
  if (firstTime == TRUE) {
    #load packages
    libs <- c('RSQLite','feather') 
    lapply(libs, library, character.only = TRUE)
    #set up input and paths at start up
    SPECIE_LIST <- unique(c("Human", sort(orgInfo$name2)))
    updateSelectizeInput(session = session, inputId = "userSpecie",
                         choices = SPECIE_LIST, server = TRUE)
    PATH <- paste0(path, 'convertIDs.db')
    PATH2 <- paste0(path, '/feather/example_of_id.feather')
    default <- getExampleDfID(userSpecie = SPECIE_LIST[1],
                              path2Database = PATH2)
    
    output$tableDefault <- renderReactable({
      reactable::reactable(data = default,
                           columns = list(Database = colDef(maxWidth = MAX_WIDTH_COL),
                                     Example_genes = colDef(maxWidth = EXAMPLE_GENES_WCOL)),
                           searchable = TRUE, bordered = TRUE, defaultPageSize = 4,
                           highlight = TRUE, resizable = TRUE, minRows = 5, showPageSizeOptions = TRUE,
                           pageSizeOptions = c(4, 10, 25, 50, 100))
    })#end of tableDefault
    shinyjs::hide(id = "downloadIDPage")
    firstTime <- FALSE
  }#end of if 
  
  observeEvent(input$userSpecie, { # update userIDtype when userSpecie changes
    ID_TYPE_LIST <- getExampleDfID(userSpecie = input$userSpecie,
                                   path2Database = PATH2)
    ID_TYPE_FILTER <- c("None", sort(ID_TYPE_LIST$Database))
    
    updateSelectizeInput(session = session, inputId = "userIDtype",
                         choices = ID_TYPE_FILTER, server = TRUE)
  }) # end of observeEvent
  
  observeEvent(input$submitIDPage, {
    #Decide on what to pass to function
    if (input$userIDtype == "None") {#user just gives species
      shinyjs::hide(id = "downloadIDPage")
      getExampleSer <- shiny::reactive({
        withProgress(message = 'Work be done...', value = 0, {
          result <- getExampleDfID(userSpecie = input$userSpecie,
                                   path2Database = PATH2, shiny = TRUE)
        })#end of withProgress
        return(result)
      })#end of reactive
      col <- list(Database = colDef(maxWidth = MAX_WIDTH_COL), 
                  Example_genes = colDef(maxWidth = EXAMPLE_GENES_WCOL))
    } else if (input$userIDtype != "None") { #if user doesn't give genelist
      ## and give id type
      getExampleSer <- shiny::reactive({
        withProgress(message = 'Work be done...', value = 0, {
          userSpecieNum <- orgInfo$id[orgInfo$name2 == input$userSpecie]
          result <- getExample(userSpecie = userSpecieNum,
                               path2Database = PATH,
                               userIDtype = input$userIDtype, shiny = TRUE)
        })#end of withProgress
        return(result)
      })#end of reactive
      col <- list(User_ID = colDef(maxWidth = MAX_WIDTH_COL),
                  Ensembl_ID = colDef(maxWidth = MAX_WIDTH_COL))
      
      output$downloadIDPage <- downloadHandler(
        filename = function() {
          paste0(input$userIDtype, "_mapped_to_Ensembl_ID.csv")
        },
        content = function(file) {
          downloadFile <- getExampleSer()
          colnames(downloadFile) <- c(paste(input$userIDtype),
                                      'Ensembl_ID')
          write.csv(downloadFile, file)
        }
      )#end of downloadIDPage
      shinyjs::show(id = "downloadIDPage")
    }#end of if/else
    
    res <- getExampleSer()
    shinyjs::hide(id = "tableDefault")
    output$tableResult <- renderReactable({
      reactable::reactable(data = res,
                           columns = col,
                           searchable = TRUE, defaultPageSize = 4,
                           highlight = TRUE, resizable = TRUE, minRows = 5, showPageSizeOptions = TRUE,
                           pageSizeOptions = c(4, 10, 25, 50, 100))
    })#end of tableResult
    shinyjs::show(id = "tableResult")
  }) #end of observeEvent
  
  observeEvent(input$resetIDPage, {
    shinyjs::hide(id = "tableResult")
    shinyjs::hide(id = "downloadIDPage")
  }) #end of observeEvent
}# end of GeneIDPage 