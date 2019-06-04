# download public RNA-Seq data from ARCHS4
# needs to run GSEinfo.R script to generate GSE info file
library(shiny)
library(DT) # for renderDataTable
library("rhdf5")
library(RSQLite)
#library(getDEE2)

dataPath = "../../data/readCounts/"
destination_fileH = "../../data/readCounts/human_matrix.h5"
destination_fileM = "../../data/readCounts/mouse_matrix.h5"
destination_fileH_transcript = "../../data/readCounts/human_transcript.h5"
destination_fileM_transcript = "../../data/readCounts/mouse_transcript.h5"
GEOdbFile = "../../data/readCounts/GEO.db"
source("getDEE2_v2.r")

DEE2Species = c("athaliana"
				,"celegans"
				,"dmelanogaster"
				,"drerio"
				,"ecoli"
				,"hsapiens"
				,"mmusculus"
				,"rnorvegicus"
				,"scerevisiae")
names( DEE2Species ) = c("Arabidopsis"
				,"Worm"
				,"Fly"
				,"Zebrafish"
				,"E coli"
				,"Human"
				,"Mouse"
				,"Rat"
				,"Yeast")

sqlite  <- dbDriver("SQLite")
convert <- dbConnect( sqlite, GEOdbFile, flags=SQLITE_RO)  #read only mode

# get a list of species
orgInfo <- dbGetQuery(convert, paste("select distinct species from GSEinfo " ))
orgInfo <- sort( orgInfo[,1] ) # convert data frame to vector, and sort
speciesChoice <- setNames(as.list( orgInfo ), orgInfo )

# Define server logic ----
server <- function(input, output, session) {
  # populate species
  observe({  updateSelectInput(session, "selectedSpecies", choices = speciesChoice )      })

  dataset.info <- reactive({
    dataset.info <- dbGetQuery(convert, "select * from GSEinfo")
    dataset.info$GSEID = as.character(dataset.info$GSEID)
    return(dataset.info)
  })
  
  # retrieve sample info and counts data
  Search <- reactive({
    if (is.null(input$SearchData_rows_selected))   return(NULL)


     # row selected by clicking
    iy = which( dataset.info()$Species == input$selectedSpecies )
    ix = iy[input$SearchData_rows_selected]

    GSEID =  dataset.info()$GSEID[ix]
  
    GSEID = gsub(" ","",GSEID)

 	querySTMT <- paste( "select * from sampleInfo where sample_series_id = '",
                         GSEID, "' AND species = '", input$selectedSpecies, "'", sep="")     

	results <- dbGetQuery(convert, querySTMT)

    selectedSpecies <- names(sort(table(results[,5]),decreasing=T))[1]

	if( dim(results)[1] == 0  ) 
       return(NULL)
    else if ( grepl("ARCHS4", selectedSpecies ) ) {    # ARCHS4
    #sample ids
    withProgress(message = "Parsing ARCHS4 file ...", {
    samp = results[,1]   # c("GSM1532588", "GSM1532592" )
    if( selectedSpecies == "ARCHS4_Human" ) {
      destination_file = destination_fileH
      destination_file_transcript = destination_fileH_transcript
    } else  if( selectedSpecies == "ARCHS4_Mouse" ) { 
       destination_file = destination_fileM
       destination_file_transcript = destination_fileM_transcript
    }
    # Identify columns to be extracted
    samples = h5read(destination_file, "meta/Sample_geo_accession")
    sample_locations = which(samples %in% samp)

    # extract gene expression from compressed data
    genes = h5read(destination_file, "meta/genes")
    expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
    #tissue = h5read(destination_file, "meta/Sample_source_name_ch1")
    sample_title = h5read(destination_file, "meta/Sample_title")[sample_locations]
    H5close()
    incProgress(.2)
    rownames(expression) <-paste(" ",genes)
    colnames(expression) <- paste( samples[sample_locations], sample_title, sep=" ")
    expression <- expression[,order(colnames(expression))]
    incProgress(.3)
    
    # extract transcript level expression    
    if( file.exists(destination_file_transcript) ) { 
    
    samples = h5read(destination_file_transcript, "meta/Sample_geo_accession")
    sample_locations = which(samples %in% samp)
    transcripts = h5read(destination_file_transcript, "meta/transcripts")
    incProgress(.1)  
    transcriptCounts = h5read(destination_file_transcript, "data/expression", index=list(1:length(transcripts), sample_locations))
    rownames(transcriptCounts) <-paste(" ",transcripts)
    colnames(transcriptCounts) <- paste( samples[sample_locations], sample_title, sep=" ")
    transcriptCounts <- transcriptCounts[,order(colnames(transcriptCounts))]
    } else {
      transcriptCounts = NULL
    }
    #sample information table
    results = results[ ,c(4,1:3)]

    results = results[order(results[,2]),]
    colnames(results) <- c( "GEO ID","Sample ID","Tissue","Sample Title")
    results = results[, -1] # Remove GSE number
    incProgress(1)
    if(dim(results)[1]>100) results = results[1:100,]
    })
    return( list(info = results, 
                 counts = expression, 
                 transcriptCounts = transcriptCounts ) )

    } else  {   # DEE2 data

     withProgress(message = "Downloading expression data from DEE2 server ... This can take 5 minutes. ", {   

     selectedSpecies <- gsub("DEE2_", "", selectedSpecies )
     selectedSpecies <- DEE2Species[ selectedSpecies ] # species code
     SRRlist <- results$SRR_accession
     # download data using DEE2 API
     data1 <- getDEE2(selectedSpecies, SRRlist )
     incProgress(.5)
     geneCounts <- data1$GeneCounts
     ix = match( colnames(geneCounts), results$SRR_accession )
     colnames(geneCounts) = paste( colnames(geneCounts),  results$sample_title[ix] )
     
     transcriptCounts <- data1$TxCounts
     ix = match( colnames(transcriptCounts), results$SRR_accession )
     colnames(transcriptCounts) = paste( colnames(transcriptCounts),  results$sample_title[ix] )
     
     results <- results[ , c(-4,-5)]
     incProgress(1)

    })
     return( list(info = results, 
                  counts = geneCounts, 
                  transcriptCounts = transcriptCounts, 
                  geneInfo = data1$GeneInfo, 
                  transcriptInfo = data1$TxInfo  ) )
    }

 })
  
output$samples <- renderTable({
  if (is.null(input$SearchData_rows_selected))   return(NULL)
    if (is.null(Search() )  )   return(as.matrix("No dataset found!"))
    Search()$info
  },bordered = TRUE)

output$downloadSearchedData <- downloadHandler(
  
    filename = function() { paste(selectedGSEID(),"_gene_level.csv",sep="")},
    content = function(file) {
      write.csv( Search()$counts, file )	    }
  )

output$downloadSearchedDataTranscript <- downloadHandler(
  
    filename = function() { paste(selectedGSEID(),"_transcript_level.csv",sep="")},
    content = function(file) {
      write.csv( Search()$transcriptCounts, file )	    }
  )

output$downloadSearchedDataTxInfo <- downloadHandler(
  
    filename = function() { paste(selectedGSEID(),"_transcript_info.csv",sep="")},
    content = function(file) {
      write.csv( Search()$transcriptInfo, file )	    }
  )
output$downloadSearchedDataGeneInfo <- downloadHandler(
  
    filename = function() { paste(selectedGSEID(),"_gene_info.csv",sep="")},
    content = function(file) {
      write.csv( Search()$geneInfo, file )	    }
  )
# search GSE IDs
output$SearchData <- DT::renderDataTable({
    if( is.null( dataset.info())) return(NULL)
  if( is.null( input$selectedSpecies)) return(NULL) 
	     dataset.info()[which( dataset.info()$Species == input$selectedSpecies)    ,]
	
  }, selection = 'single'
	   ,options = list(  pageLength = 5 ) # only 5 rows shown
	)

selectedGSEID <- reactive({
  if (is.null(input$SearchData_rows_selected))   return(NULL)
  # indices for a certain species
  iy = which( dataset.info()$Species == input$selectedSpecies )
  ix = iy[input$SearchData_rows_selected]
  return(   dataset.info()$GSEID[ix]  )

})

output$selectedDataset <- renderText({
  if (is.null(input$SearchData_rows_selected))   return(NULL)
  return( 
    paste( selectedGSEID() ) 
    #paste("Selected:",selectedGSEID() ) 
    )
 
})
output$stats <- renderTable({
    withProgress(message = "Parsing ARCHS4 file ...", {
	# datasets
	GSEs <- dbGetQuery(convert, "select  species, count(GSEID) from GSEinfo GROUP BY species" )
    incProgress(0.3)
	GSMs <- dbGetQuery(convert, "select  species, count(GSMs) from sampleInfo GROUP BY species" )
    incProgress(0.3)
	stats <- merge(GSEs, GSMs, by.x = "Species", by.y = "species")

	stats$Source = stats$Species
	stats$Source <- gsub("_.*", "", stats$Source)
	stats$Species <- gsub(".*_", "", stats$Species)
	stats <- stats[, c(4,1:3)]
	colnames(stats)[3:4] <- c("#Datasets", "#Samples")
    stats$Source[ which( duplicated(stats$Source)  )] <- ""
    })
    return(stats)
  },bordered = TRUE)
output$DoneLoading <- renderUI({
  i = "<h4>Done. Ready to search.</h4>"

  
  HTML(paste(i, collapse='<br/>') )
})


}
