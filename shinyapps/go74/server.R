####################################################
# Author: Steven Ge Xijin.Ge@sdstate.edu
# co-author: Eric Tulowetzke, eric.tulowetzke@jacks.sdstate.edu
# Lab: Ge Lab
# R version 4.0.5
# Project: ShinyGO v74
# File: server.R
# Purpose of file:main server logic of app
# Start data: NA (mm-dd-yyyy)
# Data last modified: 09-2-2021
#######################################################
server <- function(input, output, session){
  options(warn=-1)
  
  observe({  updateSelectizeInput(session, "selectOrg", choices = speciesChoice, selected = speciesChoice[1] )      })
  
  # for gene ID example
  observe({  updateSelectizeInput(session, "userSpecieIDexample", choices = speciesChoice, selected = speciesChoice[1] )      })  
  # load demo data when clicked
  observe({ 
    if( input$useDemo1 ) {
      updateTextInput(session, 'input_text', value = ExampleGeneList1 )
    }
   # if( input$useDemo2 ) {
   #   updateTextInput(session, 'input_text', value = ExampleGeneList2 )
   # }    
  })
  
  # update species for STRING-db related API access
  
  # tried to solve the double reflashing problems	 
  #https://stackoverflow.com/questions/30991900/avoid-double-refresh-of-plot-in-shiny
  observe({  	updateSelectizeInput(session, "speciesName", choices = sort(STRING10_species$official_name) ) 	})
  #click_saved <- reactiveValues(GO = NULL)
  #observeEvent(eventExpr = input$selectGO, handlerExpr = { click_saved$GO <- input$selectGO })
  
  # this defines an reactive object that can be accessed from other rendering functions
  converted <- reactive({
    if (input$goButton == 0)    return()
    
    convertID(input$input_text,input$selectOrg );
    
  } )
  
  
  geneInfoLookup <- reactive({
    if (input$goButton == 0)    return()
    geneInfo(converted(),input$selectOrg )   # uses converted gene ids thru converted() call
    
  } )
  
  detailedGeneInfoLookup <- reactive({
    if (input$goButton == 0)    return()
    geneInfoDetails(converted(),input$selectOrg )   # uses converted gene ids thru converted() call
    
  } )
  # this defines an reactive object that can be accessed from other rendering functions
  converted_background <- reactive({
    if (input$goButton == 0 | is.null(input$input_text_b))    return() 
    if (nchar(input$input_text_b) > 10)    return() 
    
    convertID(input$input_text_b,input$selectOrg );
    
  } )
  geneInfoLookup_background <- reactive({
    if (input$goButton == 0 | nchar(input$input_text_b) > 10)    return()
    if(is.null( converted_background() )) return()
    geneInfo(converted_background(),input$selectOrg )   # uses converted gene ids thru converted() call
    
  } )
  
  significantOverlaps <- reactive({
    if (input$goButton == 0 | is.null( input$selectGO) ) return()
    tem = input$maxTerms
    tem = input$minFDR
    tem = input$selectOrg
    tem = input$selectGO
    isolate({ 
      withProgress(message= sample(quotes,1),detail="enrichment analysis", {
        #gene info is passed to enable lookup of gene symbols
        tem = geneInfoLookup(); tem <- tem[which( tem$Set == "List"),]
        temb = geneInfoLookup_background(); 
        if(class(temb) == "data.frame")
          temb <- temb[which( temb$Set == "List"),]  	  
        FindOverlap( converted(), tem, input$selectGO, input$selectOrg, input$minFDR, input$maxTerms, 
                     converted_background(), temb )
      })
    })
  })
    
  output$species <-renderTable({
    if (input$goButton == 0)    return()
    tem = input$selectGO; tem=input$selectOrg; tem=input$minFDR
    isolate( {  #tem <- convertID(input$input_text,input$selectOrg );
      withProgress(message="Converting gene IDs", {
        tem <- converted()
        incProgress(1, detail = paste("Done"))	  })
      
      if( is.null(tem)) {as.data.frame("ID not recognized.")} else {
        tem$speciesMatched }
      
    }) # avoid showing things initially
  }, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)
 
  output$showGeneIDs4Species <-renderTable({
    if (input$userSpecieIDexample == 0)    return()
      withProgress(message="Retrieving gene IDs (2 minutes)", {
          geneIDs <- showGeneIDs(species = input$userSpecieIDexample, nGenes = 10)
        incProgress(1, detail = paste("Done"))	  })
      geneIDs
  }, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)
  
 output$orgInfoTable <- DT::renderDataTable({

     df <- orgInfo[, c("ensembl_dataset", "name", "totalGenes")]
     colnames(df) <- c("Ensembl/STRING-db ID", "Name (Assembly)", "Total Genes")
     row.names(df) <- NULL
     df
  })
  
  promoterData <-reactive({
    if (input$goButton == 0)    return()
    tem = input$radio
    tem = input$selectOrg
    isolate( {
      myMessage ="Promoter analysis";
      withProgress(message= sample(quotes,1),detail=myMessage, {
        tem <- promoter( converted(),input$selectOrg,input$radio )
        incProgress(1, detail = paste("Done"))	  })
      
      if( is.null(tem)) { return( as.data.frame("ID not recognized.") )} else {
        return(tem) }
      
    }) # avoid showing things initially
  })
  
  output$promoter <-renderTable({
    if (input$goButton == 0)    return()
    tem = input$radio; tem = input$selectOrg
    isolate( {
      promoterData()
    }) # avoid showing things initially
  }, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)
  
  output$downloadPromoter <- downloadHandler(
    filename = function() {"promoterMotif.csv"},
    content = function(file) {
      write.csv(promoterData(), file, row.names=FALSE)
    }
  )
  
  conversionTableData <- reactive({
    if (input$goButton == 0)    return()   # still have problems when geneInfo is not found!!!!!
    tem = input$selectGO; tem=input$selectOrg; tem=input$minFDR
    isolate( {

    withProgress(message= sample(quotes,1),detail="Looking up gene Info", {
      tem <- converted();
      incProgress(0.1)
      tem2 <- geneInfoLookup()
      incProgress(0.3)
      tem3 <- detailedGeneInfoLookup()   
      incProgress(0.6)   
      if( is.null(tem)) {as.data.frame("ID not recognized.")} else {
        if(dim(tem2)[1] == 1) { return( tem$conversionTable) }
        else { # if gene info is not available
          merged <- merge(tem$conversionTable,tem2,by='ensembl_gene_id')
          if(!is.null(tem3))
             merged <- merge(merged, tem3, by='ensembl_gene_id', all.x = TRUE)
          merged <- subset(merged,select=c(User_input,symbol,ensembl_gene_id,entrezgene_id, 
                                           gene_biotype,Species,chromosome_name,start_position, Description  ))
          
          tem3 <- as.data.frame(tem$originalIDs); colnames(tem3) = "User_input"
          merged <- merge(merged, tem3, all=T)
          merged$ensembl_gene_id[which(is.na(merged$ensembl_gene_id))] <- "Not mapped"
          chrName <- suppressWarnings( as.numeric( as.character(merged$chromosome_name) ))
          merged <- merged[order( merged$gene_biotype,
                                  chrName, 
                                  merged$start_position                                                
          ), ]
          merged$start_position <- merged$start_position/1e6
          colnames(merged) <- c("Pasted","Symbol", "Ensembl Gene ID",  "Entrez",
                                "Gene Type", "Species", "Chr", "Position (Mbp)", "Description" )
          i = 1:dim(merged)[1]
          merged = cbind(i,merged)


        }
      }
      incProgress(0.9)
      return(merged)
      })
    }) # avoid showing things initially
  })
  
  output$conversionTable <-renderTable({
    if (input$goButton == 0)    return()   # still have problems when geneInfo is not found!!!!!
    
    isolate( {
      df <- conversionTableData()

      # first see if it is Ensembl gene ID-----------------------
      ix <- grepl("ENS", df$'Ensembl Gene ID')  
      if(sum(ix) > 0) { # at least one has http?
        tem <- paste0("<a href='http://www.ensembl.org/id/", 
                      df$'Ensembl Gene ID',
                      "' target='_blank'>",
                      df$'Ensembl Gene ID', 
                      "</a>" )
        # only change the ones with URL
        df$'Ensembl Gene ID'[ix] <- tem[ix]
      }
      # first see if it is Ensembl gene ID-----------------------
      ix <- !is.na(as.numeric( df$Entrez) ) 
      if(sum(ix) > 0) { # at least one has http?
        tem <- paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/", 
                      df$Entrez,
                      "' target='_blank'>",
                      df$Entrez, 
                      "</a>" )
        # only change the ones with URL
        df$Entrez[ix] <- tem[ix]
      
     }
     return(df)
      
    }) # avoid showing things initially
  }, digits = 4,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T,
    sanitize.text.function = function(x) x)
  
  
  output$downloadGeneInfo <- downloadHandler(
    filename = function() {"geneInfo.csv"},
    content = function(file) {
      write.csv(conversionTableData(), file, row.names=FALSE)
    }
  )
  
  output$EnrichmentTable <-renderTable({
    if (input$goButton == 0  )    return(NULL)
    tem <- input$input_text_b; # just to make it re-calculate if user changes background

    myMessage = "Analyzing genes."

    if(is.null(significantOverlaps() )  ) return(NULL)
    # this solves an error when there is no significant enrichment

    if(ncol(significantOverlaps()$x ) ==1 ) return(significantOverlaps()$x)	
    
    withProgress(message= sample(quotes,1),detail=myMessage, {
      pathways <- significantOverlaps()$x;

      if(input$SortPathways == "Sort by FDR")
           pathways <- pathways[order(pathways[, 1]), ] 
      if(input$SortPathways == "Sort by Fold Enrichment")
           pathways <- pathways[order(pathways[, 4], decreasing = TRUE), ] 
      if(input$SortPathways == "Sort by Genes")
           pathways <- pathways[order(pathways[, 2], decreasing = TRUE), ]  
      if(input$SortPathways == "Sort by Category Name")
           pathways <- pathways[order( pathways[, 5]), ]  

      pathways$Pathway <- hyperText(pathways$Pathway, pathways$URL )
      
      pathways <- pathways[, -7]
      #rownames(pathways) <- NULL
      #pathways[, 1] <- as.numeric( format(pathways[, 1], scientific = TRUE, digits = 3 ) )
      pathways[, 4] <- as.character( round(pathways[, 4], 1))
      pathways[, 2] <- as.character(pathways[, 2]) # convert total genes into character 10/21/19
      pathways[, 3] <- as.character(pathways[, 3]) # convert total genes into character 10/21/19
      colnames(pathways)[5] <- "Pathways (click for details)"
      incProgress(1, detail = paste("Done"))	  })

    if(dim(pathways)[2] >1 ) pathways[,2] <- as.character(pathways[,2])

    if(dim(pathways)[2] ==1 ) return(pathways) else return(pathways[,1:5])  # If no significant enrichment found x only has 1 column.
  }, digits = -1, 
     spacing="s", 
     striped=TRUE,
     bordered = TRUE, 
     width = "auto",
     hover=TRUE, 
     sanitize.text.function = function(x) x )
  
  significantOverlaps2 <- reactive({
    if (input$goButton == 0  )    return()
    tem <- input$input_text_b; # just to make it re-calculate if user changes background
    
    tem <- significantOverlaps();
    if(dim(tem$x)[2] ==1 ) return(NULL)
    tem <- tem$x;
    colnames(tem) = c("adj.Pval","nGenesList","nGenesCategor","Fold", "Pathways","URL","Genes")
    tem$Pathways <- gsub(".*'_blank'>|</a>", "", tem$Pathways) # remove URL
    tem$Direction ="Diff"	
    tem
  })
  
  # duplicate of the above, with the word wrapping. This is for use in the network plot
  significantOverlaps3 <- reactive({
    if (input$goButton == 0  )    return()
    tem <- input$input_text_b; # just to make it re-calculate if user changes background
    
    tem <- significantOverlaps();
    if(dim(tem$x)[2] ==1 ) return(NULL)
    tem <- tem$x;
    colnames(tem) = c("adj.Pval","nGenesList","nGenesCategor","Fold","Pathways","URL","Genes")
    tem$Pathways <- gsub(".*'_blank'>|</a>", "", tem$Pathways) # remove URL
    if(input$wrapTextNetwork)
      tem$Pathways <- wrap_strings( tem$Pathways ) # wrap long pathway names using default width of 30 10/21/19
    
    tem$Direction ="Diff"	
    tem
  })
  
  # duplicate of the above for word wrapping in static networkplot
  significantOverlaps4 <- reactive({
    if (input$goButton == 0  )    return()
    tem <- input$input_text_b; # just to make it re-calculate if user changes background
    tem <- significantOverlaps();
    if(dim(tem$x)[2] ==1 ) return(NULL)
    tem <- tem$x;
    colnames(tem) = c("adj.Pval","nGenesList","nGenesCategor", "Fold", "Pathways","URL","Genes")
    tem$Pathways <- gsub(".*'_blank'>|</a>", "", tem$Pathways) # remove URL

    if(input$wrapTextNetworkStatic)
      tem$Pathways <- wrap_strings( tem$Pathways ) # wrap long pathway names using default width of 30 10/21/19
    
    tem$Direction ="Diff"	
    tem
  })
  
  output$GOTermsTree <- renderPlot({
    if(input$goButton == 0) return(NULL)
    
    if(is.null(significantOverlaps2() ) ) return(NULL)
    enrichmentPlot(significantOverlaps2(), 56  )
    
  }, height=770, width=1000)
  
  output$GOTermsTree4Download <- downloadHandler(
    filename = "GO_terms_Tree.tiff",
    content = function(file) {
      tiff(file, width = 10, height = 6, units = 'in', res = 300, compression = 'lzw');
      enrichmentPlot(significantOverlaps2(), 45  )
      dev.off()
    })
  
  output$enrichmentNetworkPlot <- renderPlot({
    if(is.null(significantOverlaps4())) return(NULL)
    
    enrichmentNetwork(significantOverlaps4(),layoutButton = input$layoutButtonStatic, edge.cutoff = input$edgeCutoff )
    
  }, height=900)
  
  output$enrichmentNetworkPlotDownload <- downloadHandler(
    filename = "enrichmentPlotNetworkPathway.tiff",
    content = function(file) {
      tiff(file, width = 12, height = 12, units = 'in', res = 300, compression = 'lzw')
      enrichmentNetwork(significantOverlaps4(),layoutButton = input$layoutButton, edge.cutoff = input$edgeCutoff )
      dev.off()
    })
  
  # note the same code is used twice as above. They need to be updated together!!!	  
  output$enrichmentNetworkPlotInteractive <- renderVisNetwork({
    if(is.null(significantOverlaps3())) return(NULL)

    g <- enrichmentNetwork(significantOverlaps3(),layoutButton = input$layoutButton, edge.cutoff = input$edgeCutoff )
    data1 <- toVisNetworkData(g)
    
    # Color codes: https://www.rapidtables.com/web/color/RGB_Color.html
    data1$nodes$shape <- "dot"
    # remove the color change of nodes
    #data1$nodes <- subset(data1$nodes, select = -color)
    
    data1$nodes$size <- 5 + data1$nodes$size^2 
    visNetwork(nodes = data1$nodes, edges = data1$edges, height = "700px", width = "700px")%>% 
      visIgraphLayout(layout = "layout_with_fr") %>%
      visNodes( 
        color = list(
          #background = "#32CD32",
          border = "#000000",
          highlight = "#FF8000"
        ),
        font = list(
          color = "#000000",
          size = 20
        ),
        borderWidth = 1,
        shadow = list(enabled = TRUE, size = 10)
      )  %>%
      visEdges(
        shadow = FALSE,
        color = list(color = "#A9A9A9", highlight = "#FFD700")
      ) %>% visExport(type = "jpeg", 
                      name = "export-network", 
                      float = "left", 
                      label = "Export as an image (only what's visible on the screen!)", 
                      background = "white", 
                      style= "") 
  })	
  
  output$enrichmentNetworkPlotInteractiveDownload <- downloadHandler(
    filename = "enrichmentPlotNetwork.html",
    content = function(file) {
      #jpeg(file, width = 12, height = 12, units = 'in', res = 300, compression = 'lzw')
      
      g <- enrichmentNetwork(significantOverlaps3(),layoutButton = input$layoutButton, edge.cutoff = input$edgeCutoff )
      data1 <- toVisNetworkData(g)
      
      # Color codes: https://www.rapidtables.com/web/color/RGB_Color.html
      data1$nodes$shape <- "dot"
      # remove the color change of nodes
      #data1$nodes <- subset(data1$nodes, select = -color)
      
      data1$nodes$size <- 5 + data1$nodes$size^2 
      g2 <-
        visNetwork(nodes = data1$nodes, edges = data1$edges, height = "700px", width = "700px")%>% 
        visIgraphLayout(layout = "layout_with_fr") %>%
        visNodes( 
          color = list(
            #background = "#32CD32",
            border = "#000000",
            highlight = "#FF8000"
          ),
          font = list(
            color = "#000000",
            size = 20
          ),
          borderWidth = 1,
          shadow = list(enabled = TRUE, size = 10)
        )  %>%
        visEdges(
          shadow = FALSE,
          color = list(color = "#A9A9A9", highlight = "#FFD700")
        )   %>% 
        visSave(file = file, background = "white")
      
    })
  
  output$downloadNodes <- downloadHandler(
    filename = function() {"network_nodes.csv"},
    content = function(file) {
      g <- enrichmentNetwork(significantOverlaps3(),layoutButton = input$layoutButton, edge.cutoff = input$edgeCutoff )
      data1 <- toVisNetworkData(g)
      data1$nodes$shape <- "dot"
      data1$nodes$size <- 5 + data1$nodes$size^2 
      
      write.csv(data1$nodes, file, row.names=FALSE)
    }
  )
  output$downloadEdges <- downloadHandler(
    filename = function() {"network_edges.csv"},
    content = function(file) {
      g <- enrichmentNetwork(significantOverlaps3(),layoutButton = input$layoutButton, edge.cutoff = input$edgeCutoff )
      data1 <- toVisNetworkData(g)
      data1$nodes$shape <- "dot"
      data1$nodes$size <- 5 + data1$nodes$size^2 
      
      write.csv(data1$edges, file, row.names=FALSE)
    }
  )   
  
  output$downloadEnrichment <- downloadHandler(
    filename = function() {"enrichment.csv"},
    content = function(file) {
      write.csv(significantOverlaps()$x, file, row.names=FALSE)
    }
  )
  
  
  #----------------------------------------------------
  # STRING-db functionality
  # find Taxonomy ID from species official name 
  findTaxonomyID <- reactive({
    if (input$goButton == 0  )    return(NULL)
    
    if(!is.null(input$speciesName) ) { # if species name is entered
      ix = match(input$speciesName, STRING10_species$official_name)
    } else if( input$selectGO != "ID not recognized!" )
    { # if no species is entered, try to resolve species using existing info 	
      codedNames = sapply(STRING10_species$compact_name,shortSpeciesNames )
      ix = match( gsub("_.*","", converted()$species[1,1] ), codedNames)
      if(input$selectOrg != speciesChoice[[1]]) {  # if species is entered
        selectedSpecies = findSpeciesById(input$selectOrg)[1,1]
        ix = match( gsub("_.*","", selectedSpecies ), codedNames)				
      }
      
    } else return(NULL) 
    
    if(length(ix) == 0 | is.na(ix) ) return(NULL) 
    return(STRING10_species$species_id[ix])
  })
  
  
  STRINGdb_geneList <- reactive({
    
    if (input$goButton == 0  )    return(NULL)						   
    library(STRINGdb,verbose=FALSE)
    tem=input$selectOrg; 
    
    ####################################
    
    if( is.null(conversionTableData()) ) return(NULL) # this has to be outside of isolate() !!!
    #if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
    NoSig = as.data.frame("No significant enrichment found.")
    taxonomyID = findTaxonomyID()
    
    if(is.null( taxonomyID ) ) return(NULL)
    
    isolate({
      withProgress(message=sample(quotes,1), detail ="Mapping gene ids (5 minutes)", {
        
        #Intialization
        string_db <- STRINGdb$new( version=STRING_DB_VERSION, species=taxonomyID,
                                   score_threshold=0, input_directory="" )
        
        # using expression data
        genes <- conversionTableData()
        colnames(genes)[3]=c("gene")
        genes$lfc = 1
        mapped <- string_db$map(genes,"gene", removeUnmappedRows = TRUE )
        
        incProgress(1/4,detail = paste("up regulated")  )
        up= subset(mapped, lfc>0, select="STRING_id", drop=TRUE )
        
        incProgress(1/2, detail ="Down regulated")
        down= subset(mapped, lfc<0, select="STRING_id", drop=TRUE )		
        
        mappingRatio = nrow(mapped)/ nrow(genes)
        if(nrow(mapped) == 0) return(NULL) else
          return( list(up=up, down=down, ratio=mappingRatio, geneTable=mapped ) )
        incProgress(1)
      })#progress
    }) #isolate						   
    
  })
  
  output$STRINGDB_species_stat <- renderUI({
    
    tem =""
    if(is.null(input$speciesName) && !is.null(findTaxonomyID() ) ) {
      ix = match(findTaxonomyID(), STRING10_species$species_id )
      if(length(ix) !=0 && !is.na(ix) ) 
        tem = paste(tem, "If ",STRING10_species$official_name[ix], "is NOT the correct species, change below:")		
    }  else
      tem = paste(tem, " Select species below:")		
    
    
    return( HTML(tem) )
    
  }) 
  
  output$STRINGDB_mapping_stat <- renderText({
    
    if (input$goButton == 0  )    return(NULL)	
    
    if( is.null(STRINGdb_geneList() ) ) return("No genes mapped by STRINGdb. Please enter or double-check species name above.")
    if(! is.null(STRINGdb_geneList() ) ) { 
      tem=paste0( 100*round(STRINGdb_geneList()$ratio,3), "% genes mapped by STRING web server.")
      if(STRINGdb_geneList()$ratio <0.3 ) tem = paste(tem, "Warning!!! Very few gene mapped. Double check if the correct species is selected.")
      return( tem  )
    }
  }) 
  
  stringDB_GO_enrichmentData <- reactive({
    if (input$goButton == 0  ) return(-1)
    taxonomyID = findTaxonomyID()
    if(is.null(taxonomyID)) return(NULL)
    library(STRINGdb,verbose=FALSE)
      withProgress(message=sample(quotes,1), detail ="Enrichment analysis", {
        tem = input$STRINGdbGO
        #Intialization
        string_db <- STRINGdb$new( version=STRING_DB_VERSION, species=taxonomyID,
                                   score_threshold=0, input_directory="" )
        
        # using expression data
        genes <- conversionTableData()
        minGenesEnrichment <- 1
        if(is.null(genes) ) {
          return(-2) 
        } else if(dim(genes)[1] <= minGenesEnrichment ) {
          return(-2) # if has only few genes
        } else {
          # GO
          ids = STRINGdb_geneList()$up
          if( length(ids) <= minGenesEnrichment || is.null(ids)) {
            return(-2)
          }	
          incProgress(1/3)
          result <- string_db$get_enrichment(ids, category = input$STRINGdbGO, methodMT = "fdr", iea = TRUE)
          if(nrow(result) == 0 || is.null(result)) {
            return (-2)
          } else {
            if(min(result$fdr) > input$minFDR) {
              return (-2)
            }  else {
              result <- result[which(result$fdr < input$minFDR),]
              incProgress(1, detail = paste("Done")) 
              return(result)
            }#end of check minFDR
          }# check results 
        }# end of check genes if 
      })#progress
  })#end of stringDB_GO_enrichmentData
  
  output$stringDB_GO_enrichment <- renderTable({
    result <- stringDB_GO_enrichmentData()
    if (result == -1) {
      return(NULL)
    } else if(result == -2) {
      return(as.data.frame("No significant enrichment found."))
    } else {
      result <- dplyr::select(result,
                              c('fdr','number_of_genes','term',
                                'description'))
      colnames(result) <- c('FDR','nGenes','GO terms or pathways',
                            'Description')
      if(nrow(result) > 30) {
        result <- result[1:30,] 
      }
      return(result)
    } # end of if else 
  },
  digits = 4,
  spacing="s",
  include.rownames=F,
  striped=TRUE,
  bordered = TRUE,
  width = "auto",
  hover=T) #renderTable
  
  output$STRING_enrichmentDownload <- downloadHandler(
    filename = function() {
      paste0("STRING_enrichment",input$STRINGdbGO,".csv")
    },
    content = function(file) {
      write.csv(stringDB_GO_enrichmentData(), file)
    }
  ) #downloadHandler
  
  # observeEvent(input$STRING, {
  #   taxonomyID = findTaxonomyID()
  #   if(is.null(taxonomyID)) {
  #     return(NULL)		
  #   } else {
  #     stringDB_GO_enrichmentData(input = input, 
  #                                output = output,
  #                                taxonomyID = taxonomyID) 
  #   } #end of if else 
  # })# end of observeEvent

  
  output$stringDB_network1 <- renderPlot({
    library(STRINGdb)
    if (input$goButton == 0  )    return(NULL)		
    
    
    tem = input$STRINGdbGO
    tem = input$nGenesPPI
    taxonomyID = findTaxonomyID( )
    if(is.null( taxonomyID ) ) return(NULL)		
    ####################################
    
    if(is.null(STRINGdb_geneList() ) ) return(NULL)
    
    isolate({
      withProgress(message=sample(quotes,1), detail ="Enrichment analysis", {
        #Intialization
        string_db <- STRINGdb$new( version=STRING_DB_VERSION, species=taxonomyID,
                                   score_threshold=0, input_directory="" )
        # only up regulated is ploted		
        ngenes1 = input$nGenesPPI
        if(ngenes1 <2) ngenes1 = 2
        for( i in c(1:1) ) {
          incProgress(1/2,detail = paste("Plotting network")  )
          
          
          ids = STRINGdb_geneList()[[i]]
          if(length(ids)> ngenes1 )  # n of genes cannot be more than 400
            ids <- ids[1:ngenes1]
          incProgress(1/3  )
          string_db$plot_network( ids,add_link=FALSE)
          
          
        }
        
      })#progress
    }) #isolate						   
  }, width = 1000, height=600)
  
  output$stringDB_network_link <- renderUI({
    library(STRINGdb,verbose=FALSE)
    
    tem = input$STRINGdbGO
    tem = input$nGenesPPI
    taxonomyID = findTaxonomyID( )
    if(is.null( taxonomyID ) ) return(NULL)		
    
    ####################################
    if(is.null(STRINGdb_geneList() ) ) return(NULL)
    
    isolate({
      withProgress(message=sample(quotes,1), detail ="PPI Enrichment and link", {
        #Intialization
        string_db <- STRINGdb$new( version=STRING_DB_VERSION, species=taxonomyID,
                                   score_threshold=0, input_directory="" )
        # upregulated
        ids = STRINGdb_geneList()[[1]]
        
        ngenes1 = input$nGenesPPI
        if(ngenes1 <2) ngenes1 = 2
        
        if(length(ids)> ngenes1 )  # n of genes cannot be more than 400
          ids <- ids[1:ngenes1]
        incProgress(1/4  )
        link1 = string_db$get_link( ids)
        
        
        tem = paste( "<a href=\"", link1, "\" target=\"_blank\"> Click here for an interactive and annotated network </a>"  )
        #	Pval1 = string_db$get_ppi_enrichment( ids)
        #    tem2 = paste("<h5> PPI enrichment P value: ")  
        #	tem2 = paste0(tem2, sprintf("%-3.2e",Pval1[1]))
        #	tem2 = paste(tem2, ".</h5>  <h5> Small P value indicates more PPIs among your proteins than background. </h5>" )
        #	tem = paste(tem2,tem )
        return(HTML(tem))	
        
        incProgress(1  )
        
      })#progress
    }) #isolate			
    
  }) 
  
  
  output$selectGO1 <- renderUI({   # gene set for pathway analysis
    if(input$goButton == 0 ) return(NULL)
    
    choices=gmtCategory(converted(), input$selectOrg)
    if(length(choices) > 12) {   # more than 12 categories in human and mouse, we default to GOBP
      selected = "GOBP" } else { # otherwise all gene sets
        selected = "All"
      }
    
    selectInput("selectGO", label=NULL,
                choices = choices,
                selected = selected )    	
    
    
  })	
  
  output$tableDetail <-renderTable({
    if (input$goButton == 0)    return()
    
    tem <- significantOverlaps(); tem$x
    
  }, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)
  
  output$grouping <-renderTable({
    if (input$goButton == 0)    return()
    myMessage = "Just a minute. Matching your genes with level 2 and level 3 Gene Ontology biological process terms.
	       This can take up to 1 minute as we have to glue together a large number of gene names. "
    withProgress(message=sample(quotes,1), detail =myMessage , {
      
      tem <- significantOverlaps()
      
      incProgress(1, detail = paste("Done"))	  })
    tem$groupings
  }, digits = 1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)
  
  output$downloadGrouping <- downloadHandler(
    filename = function() {"GO_Gropus.csv"},
    content = function(file) {
      write.csv(significantOverlaps()$groupings, file, row.names=FALSE)
    }
  )
  
  
  output$genomePlot <- renderPlot({
    if (input$goButton == 0  )    return()
    tem=input$selectOrg; 
    isolate( {
      x = geneInfoLookup()
      converted1 = converted()
      
      #chromosomes
      if((sum(!is.na( x$chromosome_name) ) >= minGenes && length(unique(x$chromosome_name) ) > 2 ) && length(which(x$Set == "List") ) > minGenes )
      {
        freq = table( x$chromosome_name,x$Set );
        #freq <- as.matrix(freq[which(nchar(row.names(freq))<3   ),])# remove unmapped chromosomes    Removed. Chr. VII in yeast
        freq <- as.matrix(freq[which( freq[,1]/colSums(freq)[1] >.01),])
        if(dim(freq)[2] >1 && dim(freq)[1]>1 && dim(freq)[1]<100) { # some organisms do not have fully seuqence genome: chr. names: scaffold_99816
          freq <-freq[order( as.numeric(row.names(freq) )), ]
          #freq <- freq[which(freq[,2]>0), ] # remove chromosomes with no genes
          
          
          tem <- subset(x, select = c(chromosome_name,start_position) )
          chrLengthTable = aggregate(start_position~chromosome_name, data=tem,max )
          
          allUserGenes <- x[which(x$Set == "List"),]
          allUserGenes <- merge(allUserGenes, converted1$conversionTable, by = 'ensembl_gene_id'  )
          allUserGenes$preferedIDs = allUserGenes$User_input;
          if(length(unique(allUserGenes$symbol) )/dim(allUserGenes)[1] >.7  ) allUserGenes$preferedIDs = allUserGenes$symbol;
          par(mfrow=c(dim(freq)[1],1))
          for( i in 1:dim(freq)[1] ) {
            #if(freq[i,2] >0)
            {
              par(mar=c(0,0,0,0))
              plot(.1,.1,axes=F,col="white",xlab="",ylab="",xlim=c(0,1), ylim=c(0,1))
              chr = rownames(freq)[i]
              ix = match(chr, chrLengthTable$chromosome_name)
              chrLength = chrLengthTable[ix,2]
              a1 <- allUserGenes[which(allUserGenes$chromosome_name == chr),]
              # if most of the genes have gene symbol, show gene symbol
              
              a1$start_position = a1$start_position/chrLength
              y1 = .50 # vertical position, from 0 - 1, relative to bottom left.
              text(0,y1,"I" );text(1,y1,"I" ); # start and end
              text(0,y1+.2, paste("Chr:",chr,sep=""),cex=2);
              if(chrLength > 1e6) {
                text(1,y1+.2, paste(round(chrLength/1e6,1),"Mb",sep=""),cex=2)
              } else {
                text(1,y1+.2, paste(round(chrLength/1e3,0),"Kb",sep=""),cex=2)
              }
              
              segments(0,y1+.01,1,y1+.01,col="blue")
              sapply(1:dim(a1)[1], function (i) text(a1$start_position[i], y1+.03, "|"))
              if(dim(a1)[1] <100 && freq[i,2] >0)  # if more genes, do not show symbol
                sapply(1:dim(a1)[1], function (i) text(a1$start_position[i], y1, a1$preferedIDs[i], offset=0,srt=90, pos=2,cex=1.5))
              
            }}
          
        }
      }
      
      
    } )# isolate
  }, height = 3000, width = 1000)
  
  # barplots using R base graphics
  output$genePlot <- renderPlot({
    if (input$goButton == 0  )    return()
    tem=input$selectOrg; 
    isolate( {
      withProgress(message="Ploting gene characteristics", {
        x = geneInfoLookup()
        x2 = x[which(x$gene_biotype == "protein_coding"),]  # only coding for some analyses
        
        
        # background genes	--------------------------------------------------------------   
        xB = geneInfoLookup_background()
        convertedB = converted_background()	   
        if(!is.null(xB) && 
           !is.null(convertedB) && 
           length( convertedB$IDs) < maxGenesBackground + 1) { # if more than 30k genes, ignore background genes.
          
          x <- x[ x$Set == "List", ] # remove background from selected genes
          xB <- xB[ xB$Set == "Genome", ] # remove Genome genes from background
          xB$Set <- "Background"
          x <- rbind(x, xB)
          x2 <- x[which(x$gene_biotype == "protein_coding"),]  # only coding for some analyses
        }
        # end background genes
        
        
        if(dim(x)[1]>=minGenes) # only making plots if more than 20 genes
        { # only plot when there 10 genes or more   # some columns have too many missing values
          par(mfrow=c(4,1))
          par(mar=c(8,6,8,2))
          #chromosomes
          if( sum(!is.na( x$chromosome_name) ) >= minGenes && length(unique(x$chromosome_name) ) > 2 && length(which(x$Set == "List") ) > minGenes )
          {
            freq = table( x$chromosome_name,x$Set );
            freq <- as.matrix(freq[which(nchar(row.names(freq))<10   ),])# remove unmapped chromosomes
            if(dim(freq)[2] >1 && dim(freq)[1]>1 ) { # some organisms do not have fully seuqence genome: chr. names: scaffold_99816
              Pval = chisq.test(freq)$p.value
              sig = paste("Distribution of query genes on chromosomes \nChi-squared test P=",formatC(Pval, digits=2, format="G") )
              
              if( Pval <PvalGeneInfo2)  sig = paste(sig," ***" ) else 
                if( Pval <PvalGeneInfo1)  sig = paste(sig," **" ) else 
                  if( Pval <PvalGeneInfo)  sig = paste(sig," *" ) 
              
              freq <- freq[order( as.numeric(row.names(freq) )), ]
              freq[,1] <- freq[,1] *colSums(freq)[2]/colSums(freq)[1] # expected
              freq = freq[,c(2,1)] # reverse order
              barplot(t(freq), beside=TRUE,las=3,col=c("red","lightgrey"), ylab="Number of Genes",main= sig,
                      cex.lab=1.5, cex.axis= 2,cex.names=2, cex.main=1.5   )
              
              legend("topright", c("List","Expected"), pch=15, col=c("red","lightgrey"),bty="n", cex =2)
            }
          } else {             # Create empty plot
            plot(x = 0:1, y = 0:1, ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
            text(x = 0.5,y = 0.5,  # Add text to empty plot				 
                 "Chromosome plot not available.", 
                 cex = 1.8)
          }
          incProgress(1/8)
          
          
          # gene type
          if( sum(!is.na( x$gene_biotype) ) >= minGenes && length(unique(x$gene_biotype) ) > 2  && length(which(x$Set == "List") ) > minGenes ) {
            freq = table( x$gene_biotype,x$Set );
            freq <- as.matrix(freq[which( freq[,1]/colSums(freq)[1] >.01),])
            if(dim(freq)[2] >1 && dim(freq)[1]>1 ) {
              Pval = chisq.test(freq)$p.value
              sig=paste("Distribution by gene type \nChi-squared test P=",formatC(Pval, digits=2, format="G") )
              if( Pval <PvalGeneInfo2)  sig = paste(sig," ***" ) else 
                if( Pval <PvalGeneInfo1)  sig = paste(sig," **" ) else 
                  if( Pval <PvalGeneInfo)  sig = paste(sig," *" ) 
              freq <- freq[order(    freq[,1], decreasing=T), ]
              freq[,1] <- freq[,1] *colSums(freq)[2]/colSums(freq)[1]
              tem = gsub("protein_coding","Coding",rownames(freq));
              tem =gsub("pseudogene","pseudo",tem)
              tem =gsub("processed","proc",tem); 
              row.names(freq)= tem
              par(mar=c(20,6,4.1,2.1))
              freq = freq[,c(2,1)] # reverse order
              head(freq)
              
              barplot(t(freq), beside=TRUE, las=2, col=c("red","lightgrey"), ylab="Number of Genes",
                      main= sig, cex.lab=1.2, cex.axis= 1.2,cex.names=1.2, cex.main=1.2)
              legend("topright", c("List","Expected"), pch=15, col=c("red","lightgrey"),bty="n", cex=2)
              
            }
          } else {             # Create empty plot
            plot(x = 0:1, y = 0:1, ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
            text(x = 0.5,y = 0.5,  # Add text to empty plot				 
                 "Gene type plot not available.", 
                 cex = 1.8)
          }
          
          
          incProgress(1/8)
          par(mar=c(12,6,4.1,2.1))
          # N. exons
          
          if( sum(!is.na( x2$nExons) ) >= minGenes && length(unique(x2$nExons) ) > 2  && length(which(x2$Set == "List") ) > minGenes ) {
            freq = table( x2$nExons,x2$Set );
            freq <- as.matrix(freq[which( freq[,1]/colSums(freq)[1] >.02),])
            if(dim(freq)[2] >1 && dim(freq)[1]>1 ) {
              Pval = chisq.test(freq)$p.value
              sig=paste("Number of exons (coding genes only) \nChi-squared test P=",formatC(Pval, digits=2, format="G") )
              if( Pval <PvalGeneInfo2)  sig = paste(sig," ***" ) else 
                if( Pval <PvalGeneInfo1)  sig = paste(sig," **" ) else 
                  if( Pval <PvalGeneInfo)  sig = paste(sig," *" ) 
              #freq <- freq[order(    freq[,1], decreasing=T), ]
              freq[,1] <- freq[,1] *colSums(freq)[2]/colSums(freq)[1]
              freq = freq[,c(2,1)] # reverse order
              barplot(t(freq), beside=TRUE,las=2,col=c("red","lightgrey"), ylab="Number of Genes",
                      main= sig ,xlab =c("Number of exons"),cex.lab=1.5, cex.axis= 2,cex.names=1.5, cex.main=1.5)
              legend("topright", c("List","Expected"), pch=15, col=c("red","lightgrey"),bty="n",cex=2)
            }}  else {             # Create empty plot
              plot(x = 0:1, y = 0:1, ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
              text(x = 0.5,y = 0.5,  # Add text to empty plot				 
                   "Exon plot not available.", 
                   cex = 1.8)
            }
          incProgress(1/8)
          
          #Transcript count
          if( sum(!is.na( x2$transcript_count) ) >= minGenes && length(unique(x2$transcript_count) ) > 2  && length(which(x2$Set == "List") ) > minGenes ) {
            freq = table( x2$transcript_count,x2$Set );
            freq <- as.matrix(freq[which( freq[,1]/colSums(freq)[1] >.02),])
            if(dim(freq)[2] >1 && dim(freq)[1]>1 ) {
              Pval = chisq.test(freq)$p.value
              sig=paste("Number of transcript isoforms per coding gene \nChi-squared test P=",formatC(Pval, digits=2, format="G"))
              if( Pval <PvalGeneInfo2)  sig = paste(sig," ***" ) else 
                if( Pval <PvalGeneInfo1)  sig = paste(sig," **" ) else 
                  if( Pval <PvalGeneInfo)  sig = paste(sig," *" ) 
              freq <- freq[order(    freq[,1], decreasing=T), ]
              freq[,1] <- freq[,1] *colSums(freq)[2]/colSums(freq)[1]
              freq = freq[,c(2,1)] # reverse order
              barplot(t(freq), beside=TRUE,las=2,col=c("red","lightgrey"), ylab="Number of Genes",
                      main= sig,xlab =c("Number of transcripts per gene") ,cex.lab=1.5, cex.axis= 2,cex.names=1.5, cex.main=1.5 )
              legend("topright", c("List","Expected"), pch=15, col=c("red","lightgrey"),bty="n",cex=2)
            } } else {             # Create empty plot
              plot(x = 0:1, y = 0:1, ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
              text(x = 0.5,y = 0.5,  # Add text to empty plot				 
                   "Transcript plot not available.", 
                   cex = 1.8)
            }
          incProgress(1/8)
          
        } # if minGenes
        incProgress(1/8, detail = paste("Done"))	  })
    }) #isolate
  }, width=600,height = 1500)
  
  # density plots using ggplot2	
  output$genePlot2 <- renderPlot({
    if (input$goButton == 0  )    return()
    tem=input$selectOrg; 
    isolate( {
      withProgress(message="Ploting gene characteristics", {
        x = geneInfoLookup()
        x2 = x[which(x$gene_biotype == "protein_coding"),]  # only coding for some analyses
        
        # background genes	--------------------------------------------------------------   
        xB = geneInfoLookup_background()
        convertedB = converted_background()	   
        if(!is.null(xB) && 
           !is.null(convertedB) && 
           length( convertedB$IDs) < maxGenesBackground + 1) { # if more than 30k genes, ignore background genes.
          
          x <- x[ x$Set == "List", ] # remove background from selected genes
          xB <- xB[ xB$Set == "Genome", ] # remove Genome genes from background
          xB$Set <- "Background"
          x <- rbind(x, xB)
          x2 <- x[which(x$gene_biotype == "protein_coding"),]  # only coding for some analyses
        }
        # end background genes
        
        
        if(dim(x)[1]>=minGenes) # only making plots if more than 20 genes
        { # only plot when there 10 genes or more   # some columns have too many missing values
          # par(mfrow=c(10,1))
          # par(mar=c(8,6,8,2))
          
          # increase fonts
          theme_set(theme_gray(base_size = 20)) 
          
          #Coding Sequence length 
          if( sum(!is.na( x2$cds_length) ) >= minGenes && length(unique(x2$cds_length) ) > 2 
              && length(which(x2$Set == "List") ) > minGenes) {
            Pval = t.test(log(cds_length)~Set, data=x2 )$p.value
            sig = paste("P = ",formatC(Pval, digits=2, format="G"),sep="")
            if( Pval <PvalGeneInfo2)  sig = paste(sig," ***" ) else 
              if( Pval <PvalGeneInfo1)  sig = paste(sig," **" ) else 
                if( Pval <PvalGeneInfo)  sig = paste(sig," *" )  			
            
            
            
            p1 <- ggplot(x2, aes(cds_length, fill= Set, colour = Set) )+
              geom_density(alpha = 0.1) + 
              scale_x_log10() +
              labs(x = "Coding sequence length (bp)", y = "Density") +
              annotate("text",x= min(x2$cds_length)+50, y = .5, label=sig, size=8) +
              #annotate("text",x= max(x2$cds_length), y = densMode(x2$cds_length)$y, label=sig, size=8, hjust=1) +
              guides(color = guide_legend(nrow = 2)) +
              theme(legend.key = element_rect(color = NA, fill = NA), 
                    legend.key.size = unit(1.2, "line")) + 
              theme(plot.margin = unit(c(0,0,1,0), "cm"))
          } else p1 <- ggplot() + 
              geom_point() + 
              xlim(-10, 10) + 
              ylim(-10, 10) + 
              annotate("text", x = 0, y = 0, label = "Coding Sequence length plot not available.") +
              theme(legend.position = "none",
                    panel.grid = element_blank(),
                    axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank())
          
          incProgress(1/8)
          
          #Transcript length------------
          if( sum(!is.na( x2$transcript_length) ) >= minGenes && 
              length(unique(x2$transcript_length) ) > 2 && 
              length(which(x2$Set == "List") ) > minGenes ) {
            Pval = t.test(log(transcript_length)~Set, data=x2[which(!is.na(x2$transcript_length)),] )$p.value
            sig = paste("P = ",formatC(Pval, digits=2, format="G"),sep="")
            if( Pval <PvalGeneInfo2)  sig = paste(sig," ***" ) else 
              if( Pval <PvalGeneInfo1)  sig = paste(sig," **" ) else 
                if( Pval <PvalGeneInfo)  sig = paste(sig," *" ) 
            
            p2 <- ggplot(x2, aes(transcript_length, fill= Set, colour = Set) )+
              geom_density(alpha = 0.1) + 
              scale_x_log10() +
              annotate("text",x= min(x2$transcript_length)+100, y = .5, label=sig, size=8)+	
              #annotate("text",x= max(x2$transcript_length), y = densMode(x2$transcript_length)$y, label=sig, size=8, hjust=1) +
              labs(x = "Transcript length (bp)", y = "Density") +
              guides(color = guide_legend(nrow = 2)) +
              theme(legend.key = element_rect(color = NA, fill = NA), 
                    legend.key.size = unit(1.2, "line")) + 
              theme(plot.margin = unit(c(0,0,1,0), "cm"))		   
          }  else p2 <- ggplot() + 
            geom_point() + 
            xlim(-10, 10) + 
            ylim(-10, 10) + 
            annotate("text", x = 0, y = 0, label = "Transcript length plot not available.") +
            theme(legend.position = "none",
                  panel.grid = element_blank(),
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank())
          incProgress(1/8)
          
          #Genome span ------------		  
          
          if( sum(!is.na( x2$genomeSpan) ) >= minGenes && length(unique(x2$genomeSpan) ) > 2 && length(which(x2$Set == "List") ) > minGenes ) {
            Pval = t.test(log(genomeSpan)~Set, data=x2[which(!is.na(x2$genomeSpan)),] )$p.value
            sig = paste("P = ",formatC(Pval, digits=2, format="G"),sep="")
            if( Pval <PvalGeneInfo2)  sig = paste(sig," ***" ) else 
              if( Pval <PvalGeneInfo1)  sig = paste(sig," **" ) else 
                if( Pval <PvalGeneInfo)  sig = paste(sig," *" ) 
            p3 <- ggplot(x2, aes(genomeSpan, fill= Set, colour = Set) )+
              geom_density(alpha = 0.1) + 
              scale_x_log10() +
              annotate("text",x=  min(x2$genomeSpan)+200, y = .5, label=sig, size=8)+	
              #annotate("text",x= max(x2$genomeSpan), y = densMode(x2$genomeSpan)$y, label=sig, size=8, hjust=1) +
              labs(x = "Genome span (bp)", y = "Density") +
              guides(color = guide_legend(nrow = 2)) +
              theme(legend.key = element_rect(color = NA, fill = NA), 
                    legend.key.size = unit(1.2, "line")) + 
              theme(plot.margin = unit(c(0,0,1,0), "cm"))			   
          } else p3 <- ggplot() + 
            geom_point() + 
            xlim(-10, 10) + 
            ylim(-10, 10) + 
            annotate("text", x = 0, y = 0, label = "Genome span plot not available.") +
            theme(legend.position = "none",
                  panel.grid = element_blank(),
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank())		  
          
          incProgress(1/8)
          
          #5' UTR ------------
          
          if( sum(!is.na( x2$FiveUTR) ) >= minGenes && length(unique(x2$FiveUTR) ) > 2 && length(which(x2$Set == "List") ) > minGenes ) {
            Pval = t.test(log(FiveUTR)~Set, data=x2[which(!is.na(x2$FiveUTR) &x2$FiveUTR > 0 ),] )$p.value
            sig = paste("P = ",formatC(Pval, digits=2, format="G"),sep="")
            if( Pval <PvalGeneInfo2)  sig = paste(sig," ***" ) else 
              if( Pval <PvalGeneInfo1)  sig = paste(sig," **" ) else 
                if( Pval <PvalGeneInfo)  sig = paste(sig," *" ) 
            
            p4 <- ggplot(x2, aes(FiveUTR, fill= Set, colour = Set) )+
              geom_density(alpha = 0.1) + 
              scale_x_log10() +
              annotate("text",x= min(x2[ which(!is.na(x2$FiveUTR) &x2$FiveUTR > 0 ),'FiveUTR'])+5, 
                       y = .5, label=sig, size=8)+	
              #annotate("text",x= max(x2$FiveUTR), y = densMode(x2$FiveUTR)$y, label=sig, size=8, hjust=1) +
              labs(x = "5' UTR length (bp)", y = "Density")  +
              guides(color = guide_legend(nrow = 2)) +
              theme(legend.key = element_rect(color = NA, fill = NA), 
                    legend.key.size = unit(1.2, "line")) + 
              theme(plot.margin = unit(c(0,0,1,0), "cm"))	
          }	else p4 <- ggplot() + 
            geom_point() + 
            xlim(-10, 10) + 
            ylim(-10, 10) + 
            annotate("text", x = 0, y = 0, label = "5' UTR plot not available.") +
            theme(legend.position = "none",
                  panel.grid = element_blank(),
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank())
          
          incProgress(1/8)
          
          #3' UTR ------------	
          if( sum(!is.na( x2$ThreeUTR) ) >= minGenes && length(unique(x2$ThreeUTR) ) > 2 && length(which(x2$Set == "List") ) > minGenes ) {
            Pval = t.test(log(ThreeUTR)~Set, data=x2[which(!is.na(x2$ThreeUTR)&x2$ThreeUTR > 0 ),] )$p.value
            sig = paste("P = ",formatC(Pval, digits=2, format="G"),sep="")
            if( Pval <PvalGeneInfo2)  sig = paste(sig," ***" ) else 
              if( Pval <PvalGeneInfo1)  sig = paste(sig," **" ) else 
                if( Pval <PvalGeneInfo)  sig = paste(sig," *" ) 
            
            p5 <- ggplot(x2, aes(ThreeUTR, fill= Set, colour = Set) )+
              geom_density(alpha = 0.1) + 
              scale_x_log10() +
              annotate("text",x= min(x2[ which(!is.na(x2$ThreeUTR) &x2$ThreeUTR > 0 ),'ThreeUTR'])+5, y = .5, label=sig, size=8)+	
              #annotate("text",x= max(x2$ThreeUTR), y = densMode(x2$ThreeUTR)$y, label=sig, size=8, hjust=1) +
              labs(x = "3' UTR length (bp)", y = "Density")  +
              guides(color = guide_legend(nrow = 2)) +
              theme(legend.key = element_rect(color = NA, fill = NA), 
                    legend.key.size = unit(1.2, "line")) + 
              theme(plot.margin = unit(c(0,0,1,0), "cm"))	
          }	else p5 <- ggplot() + 
            geom_point() + 
            xlim(-10, 10) + 
            ylim(-10, 10) + 
            annotate("text", x = 0, y = 0, label = "3' UTR plot not available.") +
            theme(legend.position = "none",
                  panel.grid = element_blank(),
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank())
          incProgress(1/8)
          
          #GC content ------------		  
          if( sum(!is.na( x2$percentage_gc_content) ) >= minGenes && 
              length(unique(x2$percentage_gc_content) ) > 2 && 
              length(which(x2$Set == "List") ) > minGenes ) {
            Pval = t.test(percentage_gc_content~Set, 
                          data=x2[which(!is.na(x2$percentage_gc_content) &x2$percentage_gc_content > 0 ),] )$p.value
            sig = paste("P = ",formatC(Pval, digits=2, format="G"),sep="")
            if( Pval <PvalGeneInfo2)  sig = paste(sig," ***" ) else 
              if( Pval <PvalGeneInfo1)  sig = paste(sig," **" ) else 
                if( Pval <PvalGeneInfo)  sig = paste(sig," *" ) 
            
            p6 <- ggplot(x2, aes(percentage_gc_content, fill= Set, colour = Set) )+
              geom_density(alpha = 0.1) + 
              #annotate("text",x= min(x2$percentage_gc_content)+5, y = .02, label=sig, size=8)+	
              annotate("text",x= max(x2$percentage_gc_content), y = densMode(x2$percentage_gc_content)$y, label=sig, size=8, hjust=1) +
              labs(x = "GC content (%)", y = "Density")  +
              guides(color = guide_legend(nrow = 2)) +
              theme(legend.key = element_rect(color = NA, fill = NA), 
                    legend.key.size = unit(1.2, "line")) + 
              theme(plot.margin = unit(c(0,0,1,0), "cm"))		   
          }	else p6 <- ggplot() + 
            geom_point() + 
            xlim(-10, 10) + 
            ylim(-10, 10) + 
            annotate("text", x = 0, y = 0, label = "GC content plot not available.") +
            theme(legend.position = "none",
                  panel.grid = element_blank(),
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank())
          
          incProgress(1/8)	
          grid.arrange(p1,p2,p3,p4,p5,p6, ncol=1)
          
          
        }
        incProgress(1/8, detail = paste("Done"))	  })
    }) #isolate
  }, width=700,height = 3000)
  
  output$enrichChart <- renderPlot({
    if (input$goButton == 0  )    return()

    if(is.null(significantOverlaps() )  ) return(NULL)
    if(ncol(significantOverlaps()$x ) == 1 ) return(NULL) # no significant ones found.
    tem=input$selectOrg; tem = input$SortPathwaysPlot
    tem = input$SortPathwaysPlotX  
    tem = input$SortPathwaysPlotSize
    tem = input$SortPathwaysPlotColor
    tem = input$SortPathwaysPlotFontSize
    tem = input$SortPathwaysPlotMarkerSize
    tem = input$SortPathwaysPlotHighColor
    tem = input$SortPathwaysPlotLowColor
    tem = input$enrichChartType

    isolate( {
 
        goTable <- significantOverlaps()$x[, 1:5]

        # Remove spaces in col names
        colnames(goTable) <- gsub(" ","", colnames(goTable))



        x       = input$SortPathwaysPlotX  
        size    = input$SortPathwaysPlotSize
        colorBy = input$SortPathwaysPlotColor
        
        # convert to vector so that we can look up the readable names of columns 
        columns <- unlist(columnSelection)
        
        goTable$EnrichmentFDR = -log10(goTable$EnrichmentFDR)
        ix <- which(colnames(goTable) == input$SortPathwaysPlot)

        # sort the pathways
        if(ix >0 && ix < dim(goTable)[2])
            goTable <- goTable[order(goTable[, ix], decreasing = TRUE), ]
        # convert to factor so that the levels are not reordered by ggplot2
        goTable$Pathway <- factor(goTable$Pathway, levels = rev(goTable$Pathway) )                  

        p <- ggplot(goTable, aes_string(x=x, y="Pathway", size=size, color=colorBy)) +
               geom_point() +
               scale_color_continuous(low=input$SortPathwaysPlotLowColor, 
                                      high=input$SortPathwaysPlotHighColor, 
                                      name = names(columns)[columns == colorBy],
                                      guide=guide_colorbar(reverse=TRUE)) +
               scale_size(range=c(1, input$SortPathwaysPlotMarkerSize)) +
               xlab( names(columns)[columns == x]  ) +
               ylab(NULL) + 
               guides(size  = guide_legend(order = 2, title = names(columns)[columns == size]), 
                      color = guide_colorbar(order = 1)) +
               theme(axis.text=element_text( size = input$SortPathwaysPlotFontSize) ) 

        if(input$enrichChartType == "dotplot") {
          p <- p 
        } else if(input$enrichChartType == "lollipop") {
          p <- p + 
                 geom_segment(aes_string(x = 0, 
                                            xend = x, 
                                            y = "Pathway", 
                                            yend = "Pathway"),
                                  size=1) 
        } else if(input$enrichChartType == "barplot") {
          p <- ggplot(goTable, aes_string(x=x, y="Pathway", fill=colorBy)) +
               geom_col(width = 0.8, position = position_dodge(0.7)) +
               scale_fill_continuous(low=input$SortPathwaysPlotLowColor, 
                                      high=input$SortPathwaysPlotHighColor, 
                                      name = names(columns)[columns == colorBy],
                                      guide=guide_colorbar(reverse=TRUE)) +
               xlab( names(columns)[columns == x]  ) +
               ylab(NULL) + 
               theme(axis.text=element_text( size = input$SortPathwaysPlotFontSize) ) 

        }    
        


         p 
          


    }) #isolate
    
  } ) #,height = 800,height=400
    
  output$listSigPathways <- renderUI({
    tem = input$selectOrg
    if (input$goButton == 0 | is.null(significantOverlaps())) return(NULL)
    
    tem <- significantOverlaps();
    
    if(dim(tem$x)[2] ==1 ) return(NULL)  
    choices = tem$x[,5]		
    selectInput("sigPathways", label="Select a significant KEGG pathway to show diagram with your genes highlighted in red:"
                ,choices=choices)      
  })
  
  
  output$KeggImage <- renderImage({
    
    tem = input$sigPathways
    
    # First generate a blank image. Otherwise return(NULL) gives us errors.
    outfile <- tempfile(fileext='.png')
    png(outfile, width=400, height=300)
    
    frame()
    dev.off()
    blank <- list(src = outfile,
                  contentType = 'image/png',
                  width = 400,
                  height = 300,
                  alt = " ")	
    if (input$goButton == 0  )    return(blank)
    if(is.null( input$selectGO ) ) return(blank)
    if(input$selectGO != "KEGG") return(blank)
    if( is.null(significantOverlaps())) return(blank)	
    
    library(pathview,verbose=FALSE)
    
    # these two functions are from the pathview package, modified to write to a designated folder: temp.
    mypathview <- function (gene.data = NULL, cpd.data = NULL, pathway.id, species = "hsa", 
                            kegg.dir = ".", cpd.idtype = "kegg", gene.idtype = "entrez", 
                            gene.annotpkg = NULL, min.nnodes = 3, kegg.native = TRUE, 
                            map.null = TRUE, expand.node = FALSE, split.group = FALSE, 
                            map.symbol = TRUE, map.cpdname = TRUE, node.sum = "sum", 
                            discrete = list(gene = FALSE, cpd = FALSE), limit = list(gene = 1, 
                                                                                     cpd = 1), bins = list(gene = 10, cpd = 10), both.dirs = list(gene = T, 
                                                                                                                                                  cpd = T), trans.fun = list(gene = NULL, cpd = NULL), 
                            low = list(gene = "green", cpd = "blue"), mid = list(gene = "gray", 
                                                                                 cpd = "gray"), high = list(gene = "red", cpd = "yellow"), 
                            na.col = "transparent", ...) 
    {
      dtypes = !is.null(gene.data) + (!is.null(cpd.data))
      cond0 = dtypes == 1 & is.numeric(limit) & length(limit) > 
        1
      if (cond0) {
        if (limit[1] != limit[2] & is.null(names(limit))) 
          limit = list(gene = limit[1:2], cpd = limit[1:2])
      }
      if (is.null(trans.fun)) 
        trans.fun = list(gene = NULL, cpd = NULL)
      arg.len2 = c("discrete", "limit", "bins", "both.dirs", "trans.fun", 
                   "low", "mid", "high")
      for (arg in arg.len2) {
        obj1 = eval(as.name(arg))
        if (length(obj1) == 1) 
          obj1 = rep(obj1, 2)
        if (length(obj1) > 2) 
          obj1 = obj1[1:2]
        obj1 = as.list(obj1)
        ns = names(obj1)
        if (length(ns) == 0 | !all(c("gene", "cpd") %in% ns)) 
          names(obj1) = c("gene", "cpd")
        assign(arg, obj1)
      }
      if (is.character(gene.data)) {
        gd.names = gene.data
        gene.data = rep(1, length(gene.data))
        names(gene.data) = gd.names
        both.dirs$gene = FALSE
        ng = length(gene.data)
        nsamp.g = 1
      }
      else if (!is.null(gene.data)) {
        if (length(dim(gene.data)) == 2) {
          gd.names = rownames(gene.data)
          ng = nrow(gene.data)
          nsamp.g = 2
        }
        else if (is.numeric(gene.data) & is.null(dim(gene.data))) {
          gd.names = names(gene.data)
          ng = length(gene.data)
          nsamp.g = 1
        }
        else stop("wrong gene.data format!")
      }
      else if (is.null(cpd.data)) {
        stop("gene.data and cpd.data are both NULL!")
      }
      gene.idtype = toupper(gene.idtype)
      data(bods)
      if (species != "ko") {
        species.data = kegg.species.code(species, na.rm = T, 
                                         code.only = FALSE)
      }
      else {
        species.data = c(kegg.code = "ko", entrez.gnodes = "0", 
                         kegg.geneid = "K01488", ncbi.geneid = NA, ncbi.proteinid = NA, 
                         uniprot = NA)
        gene.idtype = "KEGG"
        msg.fmt = "Only KEGG ortholog gene ID is supported, make sure it looks like \"%s\"!"
        msg = sprintf(msg.fmt, species.data["kegg.geneid"])
        message("Note: ", msg)
      }
      if (length(dim(species.data)) == 2) {
        message("Note: ", "More than two valide species!")
        species.data = species.data[1, ]
      }
      species = species.data["kegg.code"]
      entrez.gnodes = species.data["entrez.gnodes"] == 1
      if (is.na(species.data["ncbi.geneid"])) {
        if (!is.na(species.data["kegg.geneid"])) {
          msg.fmt = "Mapping via KEGG gene ID (not Entrez) is supported for this species,\nit looks like \"%s\"!"
          msg = sprintf(msg.fmt, species.data["kegg.geneid"])
          message("Note: ", msg)
        }
        else {
          stop("This species is not annotated in KEGG!")
        }
      }
      if (is.null(gene.annotpkg)) 
        gene.annotpkg = bods[match(species, bods[, 3]), 1]
      if (length(grep("ENTREZ|KEGG|NCBIPROT|UNIPROT", gene.idtype)) < 
          1 & !is.null(gene.data)) {
        if (is.na(gene.annotpkg)) 
          stop("No proper gene annotation package available!")
        if (!gene.idtype %in% gene.idtype.bods[[species]]) 
          stop("Wrong input gene ID type!")
        gene.idmap = id2eg(gd.names, category = gene.idtype, 
                           pkg.name = gene.annotpkg, unique.map = F)
        gene.data = mol.sum(gene.data, gene.idmap)
        gene.idtype = "ENTREZ"
      }
      if (gene.idtype != "KEGG" & !entrez.gnodes & !is.null(gene.data)) {
        id.type = gene.idtype
        if (id.type == "ENTREZ") 
          id.type = "ENTREZID"
        kid.map = names(species.data)[-c(1:2)]
        kid.types = names(kid.map) = c("KEGG", "ENTREZID", "NCBIPROT", 
                                       "UNIPROT")
        kid.map2 = gsub("[.]", "-", kid.map)
        kid.map2["UNIPROT"] = "up"
        if (is.na(kid.map[id.type])) 
          stop("Wrong input gene ID type for the species!")
        message("Info: Getting gene ID data from KEGG...")
        gene.idmap = keggConv(kid.map2[id.type], species)
        message("Info: Done with data retrieval!")
        kegg.ids = gsub(paste(species, ":", sep = ""), "", names(gene.idmap))
        in.ids = gsub(paste0(kid.map2[id.type], ":"), "", gene.idmap)
        gene.idmap = cbind(in.ids, kegg.ids)
        gene.data = mol.sum(gene.data, gene.idmap)
        gene.idtype = "KEGG"
      }
      if (is.character(cpd.data)) {
        cpdd.names = cpd.data
        cpd.data = rep(1, length(cpd.data))
        names(cpd.data) = cpdd.names
        both.dirs$cpd = FALSE
        ncpd = length(cpd.data)
      }
      else if (!is.null(cpd.data)) {
        if (length(dim(cpd.data)) == 2) {
          cpdd.names = rownames(cpd.data)
          ncpd = nrow(cpd.data)
        }
        else if (is.numeric(cpd.data) & is.null(dim(cpd.data))) {
          cpdd.names = names(cpd.data)
          ncpd = length(cpd.data)
        }
        else stop("wrong cpd.data format!")
      }
      if (length(grep("kegg", cpd.idtype)) < 1 & !is.null(cpd.data)) {
        data(rn.list)
        cpd.types = c(names(rn.list), "name")
        cpd.types = tolower(cpd.types)
        cpd.types = cpd.types[-grep("kegg", cpd.types)]
        if (!tolower(cpd.idtype) %in% cpd.types) 
          stop("Wrong input cpd ID type!")
        cpd.idmap = cpd2kegg(cpdd.names, in.type = cpd.idtype)
        cpd.data = mol.sum(cpd.data, cpd.idmap)
      }
      warn.fmt = "Parsing %s file failed, please check the file!"
      if (length(grep(species, pathway.id)) > 0) {
        pathway.name = pathway.id
        pathway.id = gsub(species, "", pathway.id)
      }
      else pathway.name = paste(species, pathway.id, sep = "")
      kfiles = list.files(path = kegg.dir, pattern = "[.]xml|[.]png")
      npath = length(pathway.id)
      out.list = list()
      tfiles.xml = paste(pathway.name, "xml", sep = ".")
      tfiles.png = paste(pathway.name, "png", sep = ".")
      if (kegg.native) 
        ttype = c("xml", "png")
      else ttype = "xml"
      xml.file <- paste(kegg.dir, "/", tfiles.xml, sep = "")
      for (i in 1:npath) {
        if (kegg.native) 
          tfiles = c(tfiles.xml[i], tfiles.png[i])
        else tfiles = tfiles.xml[i]
        if (!all(tfiles %in% kfiles)) {
          dstatus = download.kegg(pathway.id = pathway.id[i], 
                                  species = species, kegg.dir = kegg.dir, file.type = ttype)
          if (dstatus == "failed") {
            warn.fmt = "Failed to download KEGG xml/png files, %s skipped!"
            warn.msg = sprintf(warn.fmt, pathway.name[i])
            message("Warning: ", warn.msg)
            return(invisible(0))
          }
        }
        if (kegg.native) {
          node.data = try(node.info(xml.file[i]), silent = T)
          if (class(node.data) == "try-error") {
            warn.msg = sprintf(warn.fmt, xml.file[i])
            message("Warning: ", warn.msg)
            return(invisible(0))
          }
          node.type = c("gene", "enzyme", "compound", "ortholog")
          sel.idx = node.data$type %in% node.type
          nna.idx = !is.na(node.data$x + node.data$y + node.data$width + 
                             node.data$height)
          sel.idx = sel.idx & nna.idx
          if (sum(sel.idx) < min.nnodes) {
            warn.fmt = "Number of mappable nodes is below %d, %s skipped!"
            warn.msg = sprintf(warn.fmt, min.nnodes, pathway.name[i])
            message("Warning: ", warn.msg)
            return(invisible(0))
          }
          node.data = lapply(node.data, "[", sel.idx)
        }
        else {
          gR1 = try(parseKGML2Graph2(xml.file[i], genes = F, 
                                     expand = expand.node, split.group = split.group), 
                    silent = T)
          node.data = try(node.info(gR1), silent = T)
          if (class(node.data) == "try-error") {
            warn.msg = sprintf(warn.fmt, xml.file[i])
            message("Warning: ", warn.msg)
            return(invisible(0))
          }
        }
        if (species == "ko") 
          gene.node.type = "ortholog"
        else gene.node.type = "gene"
        if ((!is.null(gene.data) | map.null) & sum(node.data$type == 
                                                   gene.node.type) > 1) {
          plot.data.gene = node.map(gene.data, node.data, node.types = gene.node.type, 
                                    node.sum = node.sum, entrez.gnodes = entrez.gnodes)
          kng = plot.data.gene$kegg.names
          kng.char = gsub("[0-9]", "", unlist(kng))
          if (any(kng.char > "")) 
            entrez.gnodes = FALSE
          if (map.symbol & species != "ko" & entrez.gnodes) {
            if (is.na(gene.annotpkg)) {
              warn.fmt = "No annotation package for the species %s, gene symbols not mapped!"
              warn.msg = sprintf(warn.fmt, species)
              message("Warning: ", warn.msg)
            }
            else {
              plot.data.gene$labels = NA # Try to fix this error: Error in $<-.data.frame: replacement has 97 rows, data has 103
              plot.data.gene$labels = eg2id(as.character(plot.data.gene$kegg.names), 
                                            category = "SYMBOL", pkg.name = gene.annotpkg)[, 
                                                                                           2]
              mapped.gnodes = rownames(plot.data.gene)
              node.data$labels[mapped.gnodes] = plot.data.gene$labels
            }
          }
          cols.ts.gene = node.color(plot.data.gene, limit$gene, 
                                    bins$gene, both.dirs = both.dirs$gene, trans.fun = trans.fun$gene, 
                                    discrete = discrete$gene, low = low$gene, mid = mid$gene, 
                                    high = high$gene, na.col = na.col)
        }
        else plot.data.gene = cols.ts.gene = NULL
        if ((!is.null(cpd.data) | map.null) & sum(node.data$type == 
                                                  "compound") > 1) {
          plot.data.cpd = node.map(cpd.data, node.data, node.types = "compound", 
                                   node.sum = node.sum)
          if (map.cpdname & !kegg.native) {
            plot.data.cpd$labels = cpdkegg2name(plot.data.cpd$labels)[, 
                                                                      2]
            mapped.cnodes = rownames(plot.data.cpd)
            node.data$labels[mapped.cnodes] = plot.data.cpd$labels
          }
          cols.ts.cpd = node.color(plot.data.cpd, limit$cpd, 
                                   bins$cpd, both.dirs = both.dirs$cpd, trans.fun = trans.fun$cpd, 
                                   discrete = discrete$cpd, low = low$cpd, mid = mid$cpd, 
                                   high = high$cpd, na.col = na.col)
        }
        else plot.data.cpd = cols.ts.cpd = NULL
        if (kegg.native) {
          pv.pars = my.keggview.native(plot.data.gene = plot.data.gene, 
                                       cols.ts.gene = cols.ts.gene, plot.data.cpd = plot.data.cpd, 
                                       cols.ts.cpd = cols.ts.cpd, node.data = node.data, 
                                       pathway.name = pathway.name[i], kegg.dir = kegg.dir, 
                                       limit = limit, bins = bins, both.dirs = both.dirs, 
                                       discrete = discrete, low = low, mid = mid, high = high, 
                                       na.col = na.col, ...)
        }
        else {
          pv.pars = keggview.graph(plot.data.gene = plot.data.gene, 
                                   cols.ts.gene = cols.ts.gene, plot.data.cpd = plot.data.cpd, 
                                   cols.ts.cpd = cols.ts.cpd, node.data = node.data, 
                                   path.graph = gR1, pathway.name = pathway.name[i], 
                                   map.cpdname = map.cpdname, split.group = split.group, 
                                   limit = limit, bins = bins, both.dirs = both.dirs, 
                                   discrete = discrete, low = low, mid = mid, high = high, 
                                   na.col = na.col, ...)
        }
        plot.data.gene = cbind(plot.data.gene, cols.ts.gene)
        if (!is.null(plot.data.gene)) {
          cnames = colnames(plot.data.gene)[-(1:8)]
          nsamp = length(cnames)/2
          if (nsamp > 1) {
            cnames[(nsamp + 1):(2 * nsamp)] = paste(cnames[(nsamp + 
                                                              1):(2 * nsamp)], "col", sep = ".")
          }
          else cnames[2] = "mol.col"
          colnames(plot.data.gene)[-(1:8)] = cnames
        }
        plot.data.cpd = cbind(plot.data.cpd, cols.ts.cpd)
        if (!is.null(plot.data.cpd)) {
          cnames = colnames(plot.data.cpd)[-(1:8)]
          nsamp = length(cnames)/2
          if (nsamp > 1) {
            cnames[(nsamp + 1):(2 * nsamp)] = paste(cnames[(nsamp + 
                                                              1):(2 * nsamp)], "col", sep = ".")
          }
          else cnames[2] = "mol.col"
          colnames(plot.data.cpd)[-(1:8)] = cnames
        }
        out.list[[i]] = list(plot.data.gene = plot.data.gene, 
                             plot.data.cpd = plot.data.cpd)
      }
      if (npath == 1) 
        out.list = out.list[[1]]
      else names(out.list) = pathway.name
      return(invisible(out.list))
    }# <environment: namespace:pathview>
    my.keggview.native <- function (plot.data.gene = NULL, plot.data.cpd = NULL, cols.ts.gene = NULL, 
                                    cols.ts.cpd = NULL, node.data, pathway.name, out.suffix = "pathview", 
                                    kegg.dir = ".", multi.state = TRUE, match.data = TRUE, same.layer = TRUE, 
                                    res = 400, cex = 0.25, discrete = list(gene = FALSE, cpd = FALSE), 
                                    limit = list(gene = 1, cpd = 1), bins = list(gene = 10, cpd = 10), 
                                    both.dirs = list(gene = T, cpd = T), low = list(gene = "green", 
                                                                                    cpd = "blue"), mid = list(gene = "gray", cpd = "gray"), 
                                    high = list(gene = "red", cpd = "yellow"), na.col = "transparent", 
                                    new.signature = TRUE, plot.col.key = FALSE, key.align = "x", 
                                    key.pos = "topright", ...) 
    {
      img <- readPNG(paste(kegg.dir, "/", pathway.name, ".png", 
                           sep = ""))
      width <- ncol(img)
      height <- nrow(img)
      cols.ts.gene = cbind(cols.ts.gene)
      cols.ts.cpd = cbind(cols.ts.cpd)
      nc.gene = max(ncol(cols.ts.gene), 0)
      nc.cpd = max(ncol(cols.ts.cpd), 0)
      nplots = max(nc.gene, nc.cpd)
      pn.suffix = colnames(cols.ts.gene)
      if (length(pn.suffix) < nc.cpd) 
        pn.suffix = colnames(cols.ts.cpd)
      if (length(pn.suffix) < nplots) 
        pn.suffix = 1:nplots
      if (length(pn.suffix) == 1) {
        pn.suffix = out.suffix
      }
      else pn.suffix = paste(out.suffix, pn.suffix, sep = ".")
      na.col = colorpanel2(1, low = na.col, high = na.col)
      if ((match.data | !multi.state) & nc.gene != nc.cpd) {
        if (nc.gene > nc.cpd & !is.null(cols.ts.cpd)) {
          na.mat = matrix(na.col, ncol = nplots - nc.cpd, nrow = nrow(cols.ts.cpd))
          cols.ts.cpd = cbind(cols.ts.cpd, na.mat)
        }
        if (nc.gene < nc.cpd & !is.null(cols.ts.gene)) {
          na.mat = matrix(na.col, ncol = nplots - nc.gene, 
                          nrow = nrow(cols.ts.gene))
          cols.ts.gene = cbind(cols.ts.gene, na.mat)
        }
        nc.gene = nc.cpd = nplots
      }
      out.fmt = "Working in directory %s"
      wdir = getwd()
      out.msg = sprintf(out.fmt, wdir)
      message("Info: ", out.msg)
      out.fmt = "Writing image file %s"
      multi.state = multi.state & nplots > 1
      if (multi.state) {
        nplots = 1
        pn.suffix = paste(out.suffix, "multi", sep = ".")
        if (nc.gene > 0) 
          cols.gene.plot = cols.ts.gene
        if (nc.cpd > 0) 
          cols.cpd.plot = cols.ts.cpd
      }
      for (np in 1:nplots) {
        # img.file = paste(pathway.name, pn.suffix[np], "png", 
        #    sep = ".")
        img.file = paste(kegg.dir,"/",pathway.name, ".",pn.suffix[np], ".png", 
                         sep = "")
        out.msg = sprintf(out.fmt, img.file)
        message("Info: ", out.msg)
        png(img.file, width = width, height = height, res = res)
        op = par(mar = c(0, 0, 0, 0))
        plot(c(0, width), c(0, height), type = "n", xlab = "", 
             ylab = "", xaxs = "i", yaxs = "i")
        if (new.signature) 
          img[height - 4:25, 17:137, 1:3] = 1
        if (same.layer != T) 
          rasterImage(img, 0, 0, width, height, interpolate = F)
        if (!is.null(cols.ts.gene) & nc.gene >= np) {
          if (!multi.state) 
            cols.gene.plot = cols.ts.gene[, np]
          if (same.layer != T) {
            render.kegg.node(plot.data.gene, cols.gene.plot, 
                             img, same.layer = same.layer, type = "gene", 
                             cex = cex)
          }
          else {
            img = render.kegg.node(plot.data.gene, cols.gene.plot, 
                                   img, same.layer = same.layer, type = "gene")
          }
        }
        if (!is.null(cols.ts.cpd) & nc.cpd >= np) {
          if (!multi.state) 
            cols.cpd.plot = cols.ts.cpd[, np]
          if (same.layer != T) {
            render.kegg.node(plot.data.cpd, cols.cpd.plot, 
                             img, same.layer = same.layer, type = "compound", 
                             cex = cex)
          }
          else {
            img = render.kegg.node(plot.data.cpd, cols.cpd.plot, 
                                   img, same.layer = same.layer, type = "compound")
          }
        }
        if (same.layer == T) 
          rasterImage(img, 0, 0, width, height, interpolate = F)
        pv.pars = list()
        pv.pars$gsizes = c(width = width, height = height)
        pv.pars$nsizes = c(46, 17)
        pv.pars$op = op
        pv.pars$key.cex = 2 * 72/res
        pv.pars$key.lwd = 1.2 * 72/res
        pv.pars$sign.cex = cex
        off.sets = c(x = 0, y = 0)
        align = "n"
        ucol.gene = unique(as.vector(cols.ts.gene))
        na.col.gene = ucol.gene %in% c(na.col, NA)
        if (plot.col.key & !is.null(cols.ts.gene) & !all(na.col.gene)) {
          off.sets = col.key(limit = limit$gene, bins = bins$gene, 
                             both.dirs = both.dirs$gene, discrete = discrete$gene, 
                             graph.size = pv.pars$gsizes, node.size = pv.pars$nsizes, 
                             key.pos = key.pos, cex = pv.pars$key.cex, lwd = pv.pars$key.lwd, 
                             low = low$gene, mid = mid$gene, high = high$gene, 
                             align = "n")
          align = key.align
        }
        ucol.cpd = unique(as.vector(cols.ts.cpd))
        na.col.cpd = ucol.cpd %in% c(na.col, NA)
        if (plot.col.key & !is.null(cols.ts.cpd) & !all(na.col.cpd)) {
          off.sets = col.key(limit = limit$cpd, bins = bins$cpd, 
                             both.dirs = both.dirs$cpd, discrete = discrete$cpd, 
                             graph.size = pv.pars$gsizes, node.size = pv.pars$nsizes, 
                             key.pos = key.pos, off.sets = off.sets, cex = pv.pars$key.cex, 
                             lwd = pv.pars$key.lwd, low = low$cpd, mid = mid$cpd, 
                             high = high$cpd, align = align)
        }
        if (new.signature) 
          pathview.stamp(x = 17, y = 20, on.kegg = T, cex = pv.pars$sign.cex)
        par(pv.pars$op)
        dev.off()
      }
      return(invisible(pv.pars))
    }
    
    # modify function in a package, change namespace
    # http://stackoverflow.com/questions/23279904/modifying-an-r-package-function-for-current-r-session-assigninnamespace-not-beh
    tmpfun <- get("keggview.native", envir = asNamespace("pathview"))
    environment(my.keggview.native) <- environment(tmpfun)
    attributes(my.keggview.native) <- attributes(tmpfun)  # don't know if this is really needed
    
    isolate({ 
      
      withProgress(message="Rendering KEGG pathway plot", {
        incProgress(1/5, "Loading the pathview package") 
        

        Species <- converted()$species[1,1]
        fold <- convertEnsembl2Entrez(converted()$IDs,  Species)

        fold <- fold$entrezgene_id  
        keggSpecies <- as.character( keggSpeciesID[which(keggSpeciesID[,1] == Species),3] )
        
        if(nchar( keggSpecies) <=2 ) return(blank) # not in KEGG
        
        # kegg pathway id
        incProgress(1/2, "Download pathway graph from KEGG.")
        pathID = keggPathwayID(input$sigPathways, Species, "KEGG",input$selectOrg)
        #cat("\nhere5  ",keggSpecies, " ",Species," ",input$sigPathways, "pathID:",pathID,"End", fold[1:5],names(fold)[1:5],"\n")
        #cat("\npathway:",is.na(input$sigPathways))
        #cat("\n",fold[1:5],"\n",keggSpecies,"\n",pathID)
        if(is.null(pathID) ) return(blank) # kegg pathway id not found.
        if(nchar(pathID)<3 ) return(blank)
        randomString <- gsub(".*file","",tempfile()) 
        tempFolder <- tempdir() # tempFolder = "temp";
        outfile <- paste( tempFolder,"/",pathID,".",randomString,".png",sep="")
        
        pv.out <- mypathview(gene.data = fold, pathway.id = pathID, kegg.dir = tempFolder,  out.suffix = randomString, species = keggSpecies, kegg.native=TRUE)
        
        
        # Return a list containing the filename
        list(src = outfile,
             contentType = 'image/png',
             width = "100%",
             height = "100%",
             alt = "KEGG pathway image.")
      }) 
    })
  }, deleteFile = TRUE)
  

  
}
