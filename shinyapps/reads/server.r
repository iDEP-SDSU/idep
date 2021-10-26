# download public RNA-Seq data from ARCHS4
# need to run Reads_DB_ARCHS4_DEE2.R to create GEO.db

# run in showcase mode
# runApp('~/shinyapps/reads', display.mode = "showcase")
options(shiny.maxRequestSize = 100 * 1024^2)

# library(shiny)
# library(DT) # for renderDataTable
# library("rhdf5")
 library(RSQLite)
# library(getDEE2)
library(dplyr)
# library("shinybusy")

dataPath <- "../../data/readCounts"
#dataPath <- "C:/Users/bdere/OneDrive/Documents/idep-master/idep-master/data/readCounts" #path on device

destination_fileH <- paste(dataPath, "/human_matrix_v10.h5", sep = "")
destination_fileM <- paste(dataPath, "/mouse_matrix_v10.h5", sep = "")
destination_fileH_transcript <- paste(dataPath, "/human_transcript_v10.h5", sep = "")
destination_fileM_transcript <- paste(dataPath, "/mouse_transcript_v10.h5", sep = "")
destination_file_transcript <- ""
GEOdbFile <- paste(dataPath, "/GEO.db", sep = "") # ~/idep-master/idep-master/shinyapps/reads


sqlite <- RSQLite::dbDriver("SQLite")
convert <- RSQLite::dbConnect(sqlite, GEOdbFile, flags = SQLITE_RO) # read only mode

# dbListTables(convert)
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
orgInfo <- RSQLite::dbGetQuery(convert, paste("select distinct species from GSEinfo "))
orgInfo <- sort(orgInfo[, 1]) # convert data frame to vector, and sort
speciesChoice <- setNames(as.list(orgInfo), orgInfo)






# Define server logic ----
server <- function(input, output, session) {
  GSEID <- "null"

  # updates species botton selection
  shiny::reactive(input$selectedSpecies, {
    updateRadioButtons(session, "selectedSpecies", choices = speciesChoice, selected = "ARCHS4_Human")
  })
  # disables download button when species changes
  observeEvent(input$selectedSpecies, {
    shinyjs::disable("downloadSearchedData")
    shinyjs::disable("downloadSearchedDataTranscript")
    shinyjs::disable("downloadSearchedDataTxInfo")
    shinyjs::disable("downloadSearchedDataGeneInfo")
    shinyjs::disable("downloadSearchedDataQCmat")
    shinyjs::disable("downloadSearchedDataSummaryMeta")
  })

  # enables download buttons after row is selected
  observeEvent(input$SearchData_rows_selected, {
    shinyjs::enable("downloadSearchedData")
    shinyjs::enable("downloadSearchedDataTranscript")
    shinyjs::enable("downloadSearchedDataTxInfo")
    shinyjs::enable("downloadSearchedDataGeneInfo")
    shinyjs::enable("downloadSearchedDataQCmat")
    shinyjs::enable("downloadSearchedDataSummaryMeta")
  })

  dataset.info <- shiny::reactive({
    #shinybusy::show_modal_spinner(
      #spin = "orbit",
      ##text = "Loading",
      #color = "#000000"
    #)  
    withProgress(message = "Loading...",{
    dataset.info <- RSQLite::dbGetQuery(convert, "select GSEID, Species, Samples, Title, Summary from GSEinfo")
    dataset.info$GSEID <- as.character(dataset.info$GSEID)
    dataset.info$Link <- paste("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", dataset.info$GSEID, sep = "")
    dataset.info$"NCBI GEO Link" <- paste0("<a href='", dataset.info$Link, "'>", dataset.info$GSEID, "</a>")
    #shinybusy::remove_modal_spinner()
    })
    return(dataset.info)
  })

  # retrieve sample info and counts data
  Search <- shiny::reactive({
    #shinybusy::remove_modal_spinner()
    
    
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

    results <- RSQLite::dbGetQuery(convert, querySTMT)

    selectedSpecies <- names(sort(table(results[, 5]), decreasing = T))[1]

    if (dim(results)[1] == 0) {
      return(NULL)
    } else if (grepl("ARCHS4", selectedSpecies)) {
      cat("archs4 selected")

      shiny::withProgress(message = "Parsing ARCHS4 file ...", {
        #shinybusy::show_modal_spinner(
          #spin = "orbit",
          #text = "Loading",
          #color = "#000000"
        #)
        
        
        samp <- results[, 1]
        if (selectedSpecies == "ARCHS4_Human") {
          destination_file <- destination_fileH
          destination_file_transcript <- destination_fileH_transcript
        } else if (selectedSpecies == "ARCHS4_Mouse") {
          destination_file <- destination_fileM

          destination_file_transcript <- destination_fileM_transcript
        }

        # Identify columns to be extracted
        samples <- rhdf5::h5read(destination_file, "meta/samples/geo_accession")
        sample_locations <- which(samples %in% samp)

        # extract gene expression from compressed data
        genes <- rhdf5::h5read(destination_file, "meta/genes")
        expression <- rhdf5::h5read(destination_file, "data/expression", index = list(sample_locations, 1:length(genes$genes)))
        expression <- t(expression)
        sample_title <- rhdf5::h5read(destination_file, "meta/samples/title")[sample_locations]
        rhdf5::H5close()

        shiny::incProgress(.2)

        # add row and column names
        rownames(expression) <- paste(" ", genes$genes)
        colnames(expression) <- paste(samples[sample_locations], sample_title, sep = " ")

        expression <- expression[, order(colnames(expression))] # order by colname

        shiny::incProgress(.3)

        # extract transcript level expression
        if (file.exists(destination_file_transcript)) {
          cat("archs4 transcript exists ")

          samples <- rhdf5::h5read(destination_file_transcript, "meta/samples/geo_accession")
          sample_locations <- which(samples %in% samp)

          transcripts <- rhdf5::h5read(destination_file_transcript, "meta/transcripts")
          shiny::incProgress(.1)
          transcriptCounts <- rhdf5::h5read(destination_file_transcript,
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
        a <- results %>% dplyr::select(gse, gsm, tissue, sample_title)
        results <- a
        remove(a)


        results <- results[order(results[, 2]), ] # sort results

        colnames(results) <- c("GEO ID", "Sample ID", "Tissue", "Sample Title") # rename columns
        results <- results[, -1] # Remove GSE number

        shiny::incProgress(1)

        if (dim(results)[1] > 100) results <- results[1:100, ] # only show 100 results
        
        #shinybusy::remove_modal_spinner()
        
      })
      return(list(
        info = results,
        counts = expression,
        transcriptCounts = transcriptCounts
      ))
    } else { # DEE2 data
      if (grepl("DEE2_", selectedSpecies)) {
        #shinybusy::show_modal_spinner(
          #spin = "orbit",
          #text = "Loading",
          #color = "#000000"
        #)
        #shinybusy::remove_modal_spinner()
        
        
        cat("DEE2 selected")
        shiny::withProgress(message = "Downloading expression data from DEE2 server... This can take 5 minutes. ", {
          
          #shinybusy::show_modal_spinner(
            #spin = "orbit",
            #text = "Loading",
            #color = "#000000"
          #)
          
          selectedSpecies <- gsub("DEE2_", "", selectedSpecies) # remove DEE2
          selectedSpecies <- DEE2Species[selectedSpecies] # species code
          SRRlist <- results$SRR_accession

          # download data using DEE2 API
          tmp <- tempfile()
          len <- length(SRRlist)
          if (len <= 500) {
            cat("small sample #")
            data1 <- getDEE2::getDEE2(selectedSpecies, SRRvec = SRRlist, outfile = tmp) # outfile = "myfile.zip")
            geneCounts <- as.data.frame(data1@assays@data@listData$counts)
            tmp <- paste(tmp, ".zip", sep = "")
            geneInfo <- getDEE2::loadGeneInfo(tmp)
            TranscriptInfo <- getDEE2::loadTxInfo(tmp)
            QCmat <- getDEE2::loadQcMx(tmp)
            SummaryMeta <- getDEE2::loadSummaryMeta(tmp)
          } else {
            if (len > 500) {
              cat("large sample #")
              # for large samples, we can only read 500 at a time
              iter <- floor(len / 500) + 1
              start <- 1
              end <- 500
              for (i in 1:iter) {
                data1_prime <- getDEE2::getDEE2(selectedSpecies, SRRvec = SRRlist[start:end], outfile = tmp)
                tmp <- paste(tmp, ".zip", sep = "")
                # geneInfo <- loadGeneInfo(tmp)

                data_chunk <- as.data.frame(data1_prime@assays@data@listData$counts)

                if (i == 1) {
                  # initiale data frame with data_chunk dimensions
                  df <- data.frame(matrix(nrow = dim(data_chunk[1]), ncol = 0))
                  geneInfo <- getDEE2::loadGeneInfo(tmp)
                  TranscriptInfo <- getDEE2::loadTxInfo(tmp)
                  QCmat <- getDEE2::loadQcMx(tmp)
                  SummaryMeta <- getDEE2::loadSummaryMeta(tmp)
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
            tc <- getDEE2::getDEE2(selectedSpecies, SRRlist, counts = "TxCounts")
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

          shiny::incProgress(.5)
          ix <- match(colnames(geneCounts), results$SRR_accession)
          colnames(geneCounts) <- paste(colnames(geneCounts), results$sample_title[ix])


          ix <- match(colnames(transcriptCounts), results$SRR_accession)
          colnames(transcriptCounts) <- paste(colnames(transcriptCounts), results$sample_title[ix])

          results <- results[, c(-4, -5)]
          shiny::incProgress(1)
          #shinybusy::remove_modal_spinner()
          
        })
        
        #shinybusy::remove_modal_spinner()

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
  output$selected_selectedSpecies <- shiny::renderText({
    paste("Additional information: ", input$selectedSpecies)
  })

  output$samples <- shiny::renderTable(
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

  output$samples1 <- shiny::renderTable(
    {
      if (is.null(input$SearchData_rows_selected)) {
        return(NULL)
      }
      if (is.null(Search())) {
        return(as.matrix("No dataset found!"))
      }
      if (nrow(Search()$info) < 8) {
        return(t(Search()$info))
      }
      sampletab <- as.data.frame(t(Search()$info))
      return(sampletab[, 1:8])
      # DT::datatable(st)

      # , row.names = c("<b>GSM</b>","Tissue","Sample Title")), escape = FALSE)

      ## sampletable<-t(Search()$info)


      # DT::datatable(sampletable)
    },
    bordered = TRUE,
    rownames = TRUE,
    colnames = FALSE
  )


  output$downloadSearchedData <- shiny::downloadHandler(
    filename = function() {
      paste(selectedGSEID(), "_gene_level.csv", sep = "")
    },
    content = function(file) {
      utils::write.csv(Search()$counts, file)
    }
  )

  output$downloadSearchedDataTranscript <- shiny::downloadHandler(
    filename = function() {
      paste(selectedGSEID(), "_transcript_level.csv", sep = "")
    },
    content = function(file) {
      utils::write.csv(Search()$transcriptCounts, file)
    }
  )

  output$downloadSearchedDataTxInfo <- shiny::downloadHandler(
    filename = function() {
      paste(selectedGSEID(), "_transcript_info.csv", sep = "")
    },
    content = function(file) {
      utils::write.csv(Search()$transcriptInfo, file)
    }
  )
  output$downloadSearchedDataGeneInfo <- shiny::downloadHandler(
    filename = function() {
      paste(selectedGSEID(), "_gene_info.csv", sep = "")
    },
    content = function(file) {
      utils::write.csv(Search()$geneInfo, file)
    }
  )
  output$downloadSearchedDataQCmat <- shiny::downloadHandler(
    filename = function() {
      paste(selectedGSEID(), "_QC_matrix.csv", sep = "")
    },
    content = function(file) {
      utils::write.csv(Search()$QCmat, file)
    }
  )
  output$downloadSearchedDataSummaryMeta <- shiny::downloadHandler(
    filename = function() {
      paste(selectedGSEID(), "_summary_meta.csv", sep = "")
    },
    content = function(file) {
      utils::write.csv(Search()$SummaryMeta, file)
    }
  )
  # search GSE IDs
  output$SearchData <- DT::renderDataTable(
    rownames = FALSE,
    {
      if (is.null(dataset.info())) {
        return(NULL)
      }
      if (is.null(input$selectedSpecies)) {
        return(NULL)
      }


      dataset.info()[which(dataset.info()$Species == input$selectedSpecies), c(1, 3, 4, 7)]
    },
    selection = "single",
    escape = FALSE,
    options = list(pageLength = 5), # only 5 rows shown
  )

  selectedGSEID <- shiny::reactive({
    #shinybusy::show_modal_spinner(
      #spin = "orbit",
      #text = "Loading",
      #color = "#000000"
    #)
    if (is.null(input$SearchData_rows_selected)) {
      return(NULL)
    }
    # indices for a certain species
    iy <- which(dataset.info()$Species == input$selectedSpecies)
    ix <- iy[input$SearchData_rows_selected]

    #shinybusy::remove_modal_spinner()

    return(dataset.info()$GSEID[ix])
  })

  # shows: Preview: GSE######
  output$selectedDataset <- shiny::renderText({
    if (is.null(input$SearchData_rows_selected)) {
      return(NULL)
    }
    return(
      # paste( selectedGSEID() )
      paste("Preview:", selectedGSEID())
    )
  })


  output$stats <- shiny::renderTable(
    {
      #shinybusy::show_modal_spinner(
        #spin = "orbit",
        #text = "Loading",
        #color = "#000000"
      #)
      shiny::withProgress(message = "Parsing ARCHS4 file ...", {
        # datasets
        GSEs <- RSQLite::dbGetQuery(convert, "select  species, count(GSEID) from GSEinfo GROUP BY species")
        shiny::incProgress(0.3)
        GSMs <- RSQLite::dbGetQuery(convert, "select  species, count(gsm) from sampleInfo GROUP BY species")
        shiny::incProgress(0.3)
        stats <- merge(GSEs, GSMs, by.x = "Species", by.y = "species")

        stats$Source <- stats$Species
        stats$Source <- gsub("_.*", "", stats$Source)
        stats$Species <- gsub(".*_", "", stats$Species)
        stats <- stats[, c(4, 1:3)]
        colnames(stats)[3:4] <- c("#Datasets", "#Samples")
        stats$Source[which(duplicated(stats$Source))] <- ""
      })
      #shinybusy::remove_modal_spinner()

      return(stats)
    },
    bordered = TRUE
  )
  output$DoneLoading <- shiny::renderUI({
    #shinybusy::show_modal_spinner(
      #spin = "orbit",
      #text = "Loading",
      #color = "#000000"
    #)
    i <- "<h4>Done. Ready to search.</h4>"


    shiny::HTML(paste(i, collapse = "<br/>"))
    
    #shinybusy::remove_modal_spinner()
  })

  # dbDisconnect(convert)
}
