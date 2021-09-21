# download public RNA-Seq data from ARCHS4
# need to run Reads_DB_ARCHS4_DEE2.R to create GEO.db

# run in showcase mode
# runApp('~/shinyapps/reads', display.mode = "showcase")
options(shiny.maxRequestSize = 100 * 1024^2)

library(shiny)
library(DT) # for renderDataTable
library("rhdf5")
library(RSQLite)
library(getDEE2)
library(dplyr)
  
  #dataPath <- "../../data/readCounts"
   dataPath <- "C:/Users/bdere/OneDrive/Documents/idep-master/idep-master/data/readCounts"

  destination_fileH <- paste(dataPath, "/human_matrix_v10.h5", sep="")
  destination_fileM <- paste(dataPath, "/mouse_matrix_v10.h5", sep="")
  destination_fileH_transcript <- paste(dataPath, "/human_transcript_v10.h5", sep="")
  destination_fileM_transcript <- paste(dataPath, "/mouse_transcript_v10.h5", sep="")
  destination_file_transcript <- ""
  GEOdbFile <- paste(dataPath, "/GEO.db", sep="") # ~/idep-master/idep-master/shinyapps/reads
  
  
  sqlite <- dbDriver("SQLite")
  convert <- dbConnect(sqlite, GEOdbFile, flags = SQLITE_RO) # read only mode
  
  #dbListTables(convert)
  # dbListFields(convert, 'GSEinfo')
  DEE2Species <- c(
    "athaliana",
    "celegans",
    "dmelanogaster",
    "drerio",
    "ecoli",
    "hsapiens",
    "mmusculus",
    "rnorvegicus",
    "scerevisiae"
  )
  
  
  names(DEE2Species) <- c(
    "Arabidopsis",
    "Worm",
    "Fly",
    "Zebrafish",
    "E coli",
    "Human",
    "Mouse",
    "Rat",
    "Yeast"
  )
  
  # get a list of species
  orgInfo <- dbGetQuery(convert, paste("select distinct species from GSEinfo "))
  orgInfo <- sort(orgInfo[, 1]) # convert data frame to vector, and sort
  speciesChoice <- setNames(as.list(orgInfo), orgInfo)



# Define server logic ----
server <- function(input, output, session) {
  
  GSEID <- "null"
  
  #updattes species botton selection
  reactive(input$selectedSpecies,{
     updateRadioButtons(session, "selectedSpecies", choices = speciesChoice, selected = "ARCHS4_Mouse")
    })
  #disables downlaod button when species changes
  observeEvent(input$selectedSpecies,{
    shinyjs::disable("downloadSearchedData")
    shinyjs::disable("downloadSearchedDataTranscript") 
    shinyjs::disable("downloadSearchedDataTxInfo")
    shinyjs::disable("downloadSearchedDataGeneInfo") 
    shinyjs::disable("downloadSearchedDataQCmat")
    shinyjs::disable("downloadSearchedDataSummaryMeta")
  })
 
  #enables download buttons after row is selected
  observeEvent(input$SearchData_rows_selected,{
    shinyjs::enable("downloadSearchedData")
    shinyjs::enable("downloadSearchedDataTranscript") 
    shinyjs::enable("downloadSearchedDataTxInfo")
    shinyjs::enable("downloadSearchedDataGeneInfo") 
    shinyjs::enable("downloadSearchedDataQCmat")
    shinyjs::enable("downloadSearchedDataSummaryMeta")  
  }) 
  
  dataset.info <- reactive({

    dataset.info <- dbGetQuery(convert, "select * from GSEinfo")
    dataset.info$GSEID <- as.character(dataset.info$GSEID)
    return(dataset.info)
  })
  
  # retrieve sample info and counts data
  Search <- reactive({

    if (is.null(input$SearchData_rows_selected)) {
      return(NULL)
    }

    # row selected by clicking
    iy <- which(dataset.info()$Species == input$selectedSpecies)
    ix <- iy[input$SearchData_rows_selected]
    
    GSEID <- dataset.info()$GSEID[ix]

    GSEID <- gsub(" ", "", GSEID)
    
    querySTMT <- paste("select * from sampleInfo where gse = '",
                       GSEID, "' AND species = '", input$selectedSpecies, "'",
                       sep = ""
    )
    
    results <- dbGetQuery(convert, querySTMT)
    
    selectedSpecies <- names(sort(table(results[, 5]), decreasing = T))[1]

    if (dim(results)[1] == 0) {
      return(NULL)
    } else if (grepl("ARCHS4", selectedSpecies)) {
      cat("archs4 selected")

      withProgress(message = "Parsing ARCHS4 file ...", {
        samp <- results[, 1]
        if (selectedSpecies == "ARCHS4_Human") {
          destination_file <- destination_fileH
          destination_file_transcript <- destination_fileH_transcript
        } else if (selectedSpecies == "ARCHS4_Mouse") {
          destination_file <- destination_fileM
          
          destination_file_transcript <- destination_fileM_transcript
        }
        
        # Identify columns to be extracted
        samples <- h5read(destination_file, "meta/samples/geo_accession")
        sample_locations <- which(samples %in% samp)
        
        # extract gene expression from compressed data
        genes <- h5read(destination_file, "meta/genes")
        expression <- h5read(destination_file, "data/expression", index = list(sample_locations, 1:length(genes$genes)))
        expression <- t(expression)
        sample_title <- h5read(destination_file, "meta/samples/title")[sample_locations]
        H5close()
        
        incProgress(.2)
        
        # add row and column names
        rownames(expression) <- paste(" ", genes$genes)
        colnames(expression) <- paste(samples[sample_locations], sample_title, sep = " ")
        
        expression <- expression[, order(colnames(expression))] # order by colname
        
        incProgress(.3)
        
        # extract transcript level expression
        if (file.exists(destination_file_transcript)) {
          cat("archs4 transcript exists ")
          
          samples <- h5read(destination_file_transcript, "meta/samples/geo_accession")
          sample_locations <- which(samples %in% samp)
          
          transcripts <- h5read(destination_file_transcript, "meta/transcripts")
          incProgress(.1)
          transcriptCounts <- h5read(destination_file_transcript,
                                     "data/expression",
                                     index = list(sample_locations, 1:length(transcripts$transcripts))
          )
          
          rownames(transcriptCounts) <- paste(samples[sample_locations], sample_title, sep = " ")
          colnames(transcriptCounts) <- paste(" ", transcripts$transcripts)
          transcriptCounts <- t(transcriptCounts)
          transcriptCounts <- transcriptCounts[, order(colnames(transcriptCounts))]
        } else { # transcript missing
          
          cat("no archs4 transcript")
          transcriptCounts <- NULL
        }
        
        # sample information table
        a <- results %>% select(gse, gsm, tissue, sample_title)
        results <- a
        remove(a)
        
        
        results <- results[order(results[, 2]), ] # sort results
        
        colnames(results) <- c("GEO ID", "Sample ID", "Tissue", "Sample Title") # rename columns
        results <- results[, -1] # Remove GSE number
        
        incProgress(1)
        
        if (dim(results)[1] > 100) results <- results[1:100, ] # only show 100 results
      })
      return(list(
        info = results,
        counts = expression,
        transcriptCounts = transcriptCounts
      ))
    } else { # DEE2 data
      if (grepl("DEE2_", selectedSpecies)) {
        cat("DEE2 selected")
        withProgress(message = "Downloading expression data from DEE2 server... This can take 5 minutes. ", {
          selectedSpecies <- gsub("DEE2_", "", selectedSpecies) # remove DEE2
          selectedSpecies <- DEE2Species[selectedSpecies] # species code
          SRRlist <- results$SRR_accession
          
          
          
          # download data using DEE2 API
          tmp <- tempfile()
          len <- length(SRRlist)
          if (len <= 500) {
            data1 <- getDEE2(selectedSpecies, SRRvec=SRRlist, outfile = tmp) #outfile = "myfile.zip")
            geneCounts <- as.data.frame(data1@assays@data@listData$counts)
            tmp <- paste(tmp, ".zip", sep="")
            geneInfo <- loadGeneInfo(tmp)
            TranscriptInfo <- loadTxInfo(tmp)
            # These stopped working for some reason
            QCmat <- loadQcMx(tmp)
            SummaryMeta <- loadSummaryMeta(tmp)
            
          } else {
            if (len > 500) {
            # for large samples, we can only read 500 at a time
            iter <- floor(len / 500) + 1
            start <- 1
            end <- 500
            for (i in 1:iter) {
              data1_prime <- getDEE2(selectedSpecies, SRRvec = SRRlist[start:end])#, outfile = "myfile.zip")
              geneInfo <- loadGeneInfo("myfile.zip")
              
              data_chunk <- as.data.frame(data1_prime@assays@data@listData$counts)
              
              if (i == 1) {
                # initiale data frame with data_chunk dimensions
                df <- data.frame(matrix(nrow = dim(data_chunk[1]), ncol = 0))
                geneInfo <- loadGeneInfo(tmp)
                TranscriptInfo <- loadTxInfo(tmp)
                QCmat <- loadQcMx(tmp)
                SummaryMeta <- loadSummaryMeta(tmp)
              }
              df <- cbind(df, data_chunk)
              
              # increment start value
              start <- start + 500
              
              if (i == iter) { # change range on final iteration
                end <- end + (len %% 500)
              } else if (i < iter) { # all iteration before final iteration
                end <- end + 500
              }
            }
            
            # rename df
            geneCounts <- df
            }
          }
          
          ##### Transcript Counts for large number of samples
          if (len <= 500) {
            tc <- getDEE2(selectedSpecies, SRRlist, counts = "TxCounts")
            transcriptCounts <- tc@assays@data@listData$counts
          } else if (len > 500) {
            iter <- floor(len / 500) + 1
            start <- 1
            end <- 500
            for (i in 1:iter) {
              data_chunk <- as.data.frame(getDEE2(selectedSpecies, SRRvec = SRRlist[start:end], counts = "TxCounts")@assays@data@listData$counts)
              
              if (i == 1) { # initiale data frame with data_chunk dimensions
                df <- data.frame(matrix(nrow = dim(data_chunk[1]), ncol = 0))
              }
              df <- cbind(df, data_chunk)
              
              # increment start value
              start <- start + 500
              if (i == iter) { # change range on final iteration
                end <- end + (len %% 500)
              } else if (i < iter) { # all iteration before final iteration
                end <- end + 500
              }
            }
            
            # rename df
            transcriptCounts <- df
          }
          
          incProgress(.5)
          ix <- match(colnames(geneCounts), results$SRR_accession)
          colnames(geneCounts) <- paste(colnames(geneCounts), results$sample_title[ix])
          
          
          ix <- match(colnames(transcriptCounts), results$SRR_accession)
          colnames(transcriptCounts) <- paste(colnames(transcriptCounts), results$sample_title[ix])
          
          results <- results[, c(-4, -5)]
          incProgress(1)
        })
        
        return(list(
          info = results,
          counts = geneCounts,
          transcriptCounts = transcriptCounts,
          geneInfo = geneInfo,
          transcriptInfo = TranscriptInfo,
          QCmatrix = QCmat,
          SummaryMeta = SummaryMeta
        ))
      }
    }
  })
  
  
  
  # selected species
  output$selected_selectedSpecies <- renderText({
    paste("Additional information: ", input$selectedSpecies)
  })
  
  output$samples <- renderTable(
    {
      if (is.null(input$SearchData_rows_selected)) {
        return(NULL)
      }
      if (is.null(Search())) {
        return(as.matrix("No dataset found!"))
      }
      Search()$info
    },
    bordered = TRUE
  )
  
  output$downloadSearchedData <- downloadHandler(
    filename = function() {
      paste(selectedGSEID(), "_gene_level.csv", sep = "")
    },
    content = function(file) {
      write.csv(Search()$counts, file)
    }
  )
  
  output$downloadSearchedDataTranscript <- downloadHandler(
    filename = function() {
      paste(selectedGSEID(), "_transcript_level.csv", sep = "")
    },
    content = function(file) {
      write.csv(Search()$transcriptCounts, file)
    }
  )
  
  output$downloadSearchedDataTxInfo <- downloadHandler(
    filename = function() {
      paste(selectedGSEID(), "_transcript_info.csv", sep = "")
    },
    content = function(file) {
      write.csv(Search()$transcriptInfo, file)
    }
  )
  output$downloadSearchedDataGeneInfo <- downloadHandler(
    filename = function() {
      paste(selectedGSEID(), "_gene_info.csv", sep = "")
    },
    content = function(file) {
      write.csv(Search()$geneInfo, file)
    }
  )
  output$downloadSearchedDataQCmat <- downloadHandler(
    filename = function() {
      paste(selectedGSEID(), "_QC_matrix.csv", sep = "")
    },
    content = function(file) {
      write.csv(Search()$QCmat, file)
    }
  )
  output$downloadSearchedDataSummaryMeta <- downloadHandler(
    filename = function() {
      paste(selectedGSEID(), "_summary_meta.csv", sep = "")
    },
    content = function(file) {
      write.csv(Search()$SummaryMeta, file)
    }
  )
  # search GSE IDs
  output$SearchData <- DT::renderDataTable(
    {
      if (is.null(dataset.info())) {
        return(NULL)
      }
      if (is.null(input$selectedSpecies)) {
        return(NULL)
      }
      dataset.info()[which(dataset.info()$Species == input$selectedSpecies), ]
    },
    selection = "single",
    options = list(pageLength = 5) # only 5 rows shown
  )
  
  selectedGSEID <- reactive({
    if (is.null(input$SearchData_rows_selected)) {
      return(NULL)
    }
    # indices for a certain species
    iy <- which(dataset.info()$Species == input$selectedSpecies)
    ix <- iy[input$SearchData_rows_selected]
    
    return(dataset.info()$GSEID[ix])
  })
  
  output$selectedDataset <- renderText({
    if (is.null(input$SearchData_rows_selected)) {
      return(NULL)
    }
    return(
      # paste( selectedGSEID() )
      paste("Preview:", selectedGSEID())
    )
  })
  output$stats <- renderTable(
    {
      withProgress(message = "Parsing ARCHS4 file ...", {
        # datasets
        GSEs <- dbGetQuery(convert, "select  species, count(GSEID) from GSEinfo GROUP BY species")
        incProgress(0.3)
        GSMs <- dbGetQuery(convert, "select  species, count(gsm) from sampleInfo GROUP BY species")
        incProgress(0.3)
        stats <- merge(GSEs, GSMs, by.x = "Species", by.y = "species")
        
        stats$Source <- stats$Species
        stats$Source <- gsub("_.*", "", stats$Source)
        stats$Species <- gsub(".*_", "", stats$Species)
        stats <- stats[, c(4, 1:3)]
        colnames(stats)[3:4] <- c("#Datasets", "#Samples")
        stats$Source[which(duplicated(stats$Source))] <- ""
      })
      return(stats)
    },
    bordered = TRUE
  )
  output$DoneLoading <- renderUI({
    i <- "<h4>Done. Ready to search.</h4>"
    
    
    HTML(paste(i, collapse = "<br/>"))
  })
  
  # dbDisconnect(convert)
}
