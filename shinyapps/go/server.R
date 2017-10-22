library(shiny)
library(RSQLite)

# relative path to data files
datapath = "../../data/"   # production server

Min_overlap <- 2
minSetSize = 3;
mappingCoverage = 0.60 # 60% percent genes has to be mapped for confident mapping
mappingEdge = 0.5  # Top species has 50% more genes mapped
maxTerms =30 # max number of enriched terms
PvalGeneInfo = 0.05; minGenes = 10 # min number of genes for ploting
# setwd("C:/Users/Xijin.Ge/Google Drive/research/Shiny/go")


cleanGeneSet <- function (x){
  # remove duplicate; upper case; remove special characters
  x <- unique( toupper( gsub("\n| ","",x) ) )
  x <- x[ which( nchar(x)>1) ]  # genes should have at least two characters
  return(x)
}
# read GMT files, does NO cleaning. Assumes the GMT files are created with cleanGeneSet()
readGMT <- function (fileName){
  x <- scan(fileName, what="", sep="\n")
  x <- strsplit(x, "\t")
  # Extract the first vector element and set it as the list element name
  names(x) <- sapply(x, `[[`, 1)
  x <- lapply(x, `[`, -c(1,2)) # 2nd element is comment, ignored
  x = x[which(sapply(x,length) > 1)]  # gene sets smaller than 1 is ignored!!!
  return(x)
}

sqlite  <- dbDriver("SQLite")
convert <- dbConnect(sqlite,paste0(datapath,"convertIDs.db"),flags=SQLITE_RO)  #read only mode
# keggSpeciesID = read.csv(paste0(datapath,"data_go/KEGG_Species_ID.csv"))
# List of GMT files in /gmt sub folder
gmtFiles = list.files(path = paste0(datapath,"pathwayDB"),pattern=".*\\.db")
gmtFiles = paste(datapath,"pathwayDB/",gmtFiles,sep="")
geneInfoFiles = list.files(path = paste0(datapath,"geneInfo"),pattern=".*GeneInfo\\.csv")
geneInfoFiles = paste(datapath,"geneInfo/",geneInfoFiles,sep="")
motifFiles = list.files(path = paste0(datapath,"motif"),pattern=".*\\.db")
motifFiles = paste(datapath,"motif/",motifFiles,sep="")

# # Create a list of GMT files in /gmt sub folder
# gmtFiles = list.files(path = "./pathwayDB",pattern=".*\\.db")
# gmtFiles = paste("pathwayDB/",gmtFiles,sep="")

# geneInfoFiles = list.files(path = "./geneInfo",pattern=".*GeneInfo\\.csv")
# geneInfoFiles = paste("geneInfo/",geneInfoFiles,sep="")

# motifFiles = list.files(path = "./motif",pattern=".*\\.db")
# motifFiles = paste("motif/",motifFiles,sep="")


# Create a list for Select Input options
orgInfo <- dbGetQuery(convert, paste("select distinct * from orgInfo " ))
orgInfo <- orgInfo[order(orgInfo$name),]
speciesChoice <- setNames(as.list( orgInfo$id ), orgInfo$name2 )
# add a defult element to list    # new element name       value
speciesChoice <- append( setNames( "BestMatch","Best matching species"), speciesChoice  )
# move one element to the 2nd place
move2 <- function(i) c(speciesChoice[1],speciesChoice[i],speciesChoice[-c(1,i)])
i= grep("Glycine max" ,names(speciesChoice)); speciesChoice <- move2(i)
i= grep("Zea mays" ,names(speciesChoice)); speciesChoice <- move2(i)
i= grep("Arabidopsis thaliana",names(speciesChoice)); speciesChoice <- move2(i)
i= grep("Saccharomyces cerevisiae" ,names(speciesChoice)); speciesChoice <- move2(i)
i= grep("Caenorhabditis elegans",names(speciesChoice)); speciesChoice <- move2(i)
i= grep("Danio rerio" ,names(speciesChoice)); speciesChoice <- move2(i)
i= grep("Bos taurus" ,names(speciesChoice)); speciesChoice <- move2(i)
i= grep("Rattus norvegicus" ,names(speciesChoice)); speciesChoice <- move2(i)
i= grep("Mus musculus",names(speciesChoice)); speciesChoice <- move2(i)
i= grep("Homo sapiens",names(speciesChoice)); speciesChoice <- move2(i)

GO_levels = dbGetQuery(convert, "select distinct id,level from GO
                                WHERE GO = 'biological_process'"  )
level2Terms = GO_levels[which(GO_levels$level %in% c(2,3))  ,1]  # level 2 and 3

idIndex <- dbGetQuery(convert, paste("select distinct * from idIndex " ))

quotes <- dbGetQuery(convert, " select * from quotes")
quotes = paste0("\"",quotes$quotes,"\"", " -- ",quotes$author,".       ")

# This function convert gene set names
# x="GOBP_mmu_mgi_GO:0000183_chromatin_silencing_at_rDNA"
# chromatin silencing at rDNA
proper=function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))
extract1 <- function (x) {
  words <- unlist ( strsplit(x,"_"))
  if(length( words )  <=4 ) return(gsub("_"," ",x)) else {
  words <- words[-c(1:4)]
  return( proper(paste(words,collapse = " ") ) )}
}
#find idType based on index
findIDtypeById <- function(x){ # find
  return( idIndex$idType[ as.numeric(x)] )
}
findSpeciesById <- function (speciesID){ # find species name use id
  return( orgInfo[which(orgInfo$id == speciesID),]  )
}
# just return name
findSpeciesByIdName <- function (speciesID){ # find species name use id
  return( orgInfo[which(orgInfo$id == speciesID),3]  )
}

# convert sorted species:idType combs into a list for repopulate species choice
matchedSpeciesInfo <- function (x) {
  a<- c()
  for( i in 1:length(x)) {
    a = c(a,paste( gsub("genes.*","",findSpeciesByIdName( as.numeric(gsub(" .*","",names(x[i])) ))), " (",
                   x[i]," mapped from ",findIDtypeById( gsub(".* ","",names(x[i]) ) ),")",sep="")
    ) }
  return(a )
}
# convert gene IDs to ensembl gene ids and find species
convertID <- function (query,selectOrg) {
  querySet <- cleanGeneSet( unlist( strsplit( toupper(query),'\t| |\n|\\,' )  ) )
  result <- dbGetQuery( convert,
                        paste( " select distinct id,ens,species from mapping where id IN ('", paste(querySet,collapse="', '"),   "')",sep="") )
  if( dim(result)[1] == 0  ) return(NULL)
  if(selectOrg == speciesChoice[[1]]) {
    comb = paste( result$species,result$idType)
    sortedCounts = sort( table(comb ),decreasing=T)
    recognized =names(sortedCounts[1]  )
    result <- result[which(comb == recognized )  , ]

	speciesMatched=sortedCounts
    names(speciesMatched )= sapply(as.numeric(gsub(" .*","",names(sortedCounts) ) ), findSpeciesByIdName  )
    speciesMatched <- as.data.frame( speciesMatched )
	if(length(sortedCounts) == 1) { # if only  one species matched
	  speciesMatched[1,1] <-paste( rownames(speciesMatched), "(",speciesMatched[1,1],")",sep="")
	 } else {# if more than one species matched
		speciesMatched[,1] <- as.character(speciesMatched[,1])
		speciesMatched[,1] <- paste( speciesMatched[,1]," (",speciesMatched[,2], ")", sep="")
		speciesMatched[1,1] <- paste( speciesMatched[1,1],"   ***Used in mapping***  To change, select from above and resubmit query.")
		speciesMatched <- as.data.frame(speciesMatched[,1])
	}

  } else { # if species is selected
    result <- result[which(result$species == selectOrg ) ,]
    if( dim(result)[1] == 0  ) return(NULL) #stop("ID not recognized!")
    speciesMatched <- as.data.frame(paste("Using selected species ", findSpeciesByIdName(selectOrg) )  )
  }
  result <- result[which(!duplicated(result[,2]) ),] # remove duplicates in ensembl_gene_id
  colnames(speciesMatched) = c("Matched Species (genes)" )
  conversionTable <- result[,1:2]; colnames(conversionTable) = c("User_input","ensembl_gene_id")
  conversionTable$Species = sapply(result[,3], findSpeciesByIdName )
  return(list(originalIDs = querySet,IDs=unique( result[,2]),
              species = findSpeciesById(result$species[1]),
              #idType = findIDtypeById(result$idType[1] ),
              speciesMatched = speciesMatched,
			  conversionTable = conversionTable
			  ) )
}

geneInfo <- function (converted,selectOrg){
  if(is.null(converted) ) return(as.data.frame("ID not recognized!") ) # no ID
  querySet <- converted$IDs
  if(length(querySet) == 0) return(as.data.frame("ID not recognized!") )
  ix = grep(converted$species[1,1],geneInfoFiles)
  if (length(ix) == 0 ) {return(as.data.frame("No matching gene info file found") )} else {
  # If selected species is not the default "bestMatch", use that species directly
  if(selectOrg != speciesChoice[[1]]) {
    ix = grep(findSpeciesById(selectOrg)[1,1], geneInfoFiles )
  }
  if(length(ix) == 1)  # if only one file           #WBGene0000001 some ensembl gene ids in lower case
  { x = read.csv(as.character(geneInfoFiles[ix]) ); x[,1]= toupper(x[,1]) }
  else # read in the chosen file
  { return(as.data.frame("Multiple geneInfo file found!") )   }

  Set = match(x$ensembl_gene_id, querySet)
  Set[which(is.na(Set))]="Genome"
  Set[which(Set!="Genome")] ="List"
  # x = cbind(x,Set) } # just for debuging
  return( cbind(x,Set) )}
 }

# Main function. Find a query set of genes enriched with functional category
FindOverlap <- function (converted,gInfo, GO,selectOrg,minFDR) {
  idNotRecognized = list(x=as.data.frame("ID not recognized!"),
                         groupings= as.data.frame("ID not recognized!")  )
  if(is.null(converted) ) return(idNotRecognized) # no ID
  querySet <- converted$IDs;
  if(length(querySet) == 0) return(idNotRecognized )

  ix = grep(converted$species[1,1],gmtFiles)
  totalGenes <- converted$species[1,7]

  if (length(ix) == 0 ) {return(idNotRecognized )}

  # If selected species is not the default "bestMatch", use that species directly
  if(selectOrg != speciesChoice[[1]]) {
    ix = grep(findSpeciesById(selectOrg)[1,1], gmtFiles )
	if (length(ix) == 0 ) {return(idNotRecognized )}
	totalGenes <- orgInfo[which(orgInfo$id == as.numeric(selectOrg)),7]
  }
  pathway <- dbConnect(sqlite,gmtFiles[ix])

  # Generate a list of geneset categories such as "GOBP", "KEGG" from file
  geneSetCategory <-  dbGetQuery(pathway, "select distinct * from categories " )
  geneSetCategory  <- geneSetCategory[,1]
  categoryChoices <- setNames(as.list( geneSetCategory ), geneSetCategory )
  categoryChoices <- append( setNames( "All","All available gene sets"), categoryChoices  )
  #change GOBO to the full description for display
  names(categoryChoices)[ match("GOBP",categoryChoices)  ] <- "GO Biological Process"
  names(categoryChoices)[ match("GOCC",categoryChoices)  ] <- "GO Cellular Component"
  names(categoryChoices)[ match("GOMF",categoryChoices)  ] <- "GO Molecular Function"

  sqlQuery = paste( " select distinct gene,pathwayID from pathway where gene IN ('", paste(querySet,collapse="', '"),"')" ,sep="")

  #cat(paste0("HH",GO,"HH") )

  if( GO != "All") sqlQuery = paste0(sqlQuery, " AND category ='",GO,"'")
  result <- dbGetQuery( pathway, sqlQuery  )
  if( dim(result)[1] ==0) {return(list( x=as.data.frame("No matching species or gene ID file!" )) )}

   # given a pathway id, it finds the overlapped genes, symbol preferred
  sharedGenesPrefered <- function(pathwayID) {
    tem <- result[which(result[,2]== pathwayID ),1]
	ix = match(tem, converted$conversionTable$ensembl_gene_id) # convert back to original
	tem2 <- unique( converted$conversionTable$User_input[ix] )
    if(length(unique(gInfo$symbol) )/dim(gInfo)[1] >.7  ) # if 70% genes has symbol in geneInfo
	{ ix = match(tem, gInfo$ensembl_gene_id);
	  tem2 <- unique( gInfo$symbol[ix] )      }
  return( paste( tem2 ,collapse=" ",sep="") )}



  x0 = table(result$pathwayID)
  x0 = as.data.frame( x0[which(x0>=Min_overlap)] )# remove low overlaps
  if(dim(x0)[1] <= 5 ) return(idNotRecognized) # no data
  colnames(x0)=c("pathwayID","overlap")
  pathwayInfo <- dbGetQuery( pathway, paste( " select distinct id,n,Description from pathwayInfo where id IN ('",
						paste(x0$pathwayID,collapse="', '"),   "') ",sep="") )
  x = merge(x0,pathwayInfo, by.x='pathwayID', by.y='id')

  x$Pval=phyper(x$overlap-1,length(querySet),totalGenes - length(querySet),as.numeric(x$n), lower.tail=FALSE );
  x$FDR = p.adjust(x$Pval,method="fdr")
  x <- x[ order( x$FDR)  ,]  # sort according to FDR


  # Gene groups for high level GOBP terms
  groups <- dbGetQuery( pathway,paste( " select distinct id, description from pathwayInfo
                       where golevel IN ( '2','3') ",sep="") )

  ix = match(groups$id, x0$pathwayID)
  if(length(groups)>0 && length(ix) >0)   groupings = as.data.frame("No grouping.")
  groups$ngenes <- x0$overlap[ix]
  groups <- groups[which(!is.na(ix) ),]
  groups <- groups[order(-groups$ngenes),]
  if(max(groups$ngenes)<=2) { groups = as.data.frame("Too few genes") } else {
		groupings = subset(groups,ngenes>2) # at least 10 genes
		if(dim(groups)[1] > 100) groups <- groups[1:100,]
	 groups = cbind(groups, sapply( groups$id, sharedGenesPrefered ) )
	 groups = groups[,-1]
	 groups = groups[,c(2,1,3)]
	 colnames(groups) = c("N","High level GO category", "Genes")
	}


  if(min(x$FDR) > minFDR) x=as.data.frame("No significant enrichment found!") else {
  x <- x[which(x$FDR < minFDR),]
  if(dim(x)[1] > maxTerms ) x = x[1:maxTerms,]
  x= cbind(x,sapply( x$pathwayID, sharedGenesPrefered ) )
  colnames(x)[7]= "Genes"
  x <- subset(x,select = c(FDR,overlap,n,description,Genes) )
  colnames(x) = c("Corrected P value (FDR)", "Genes in list", "Total genes in category","Functional Category","Genes"  )
  }

 dbDisconnect(pathway)
 return(list( x=x, groupings = groups, categoryChoices = categoryChoices ) )
}
                                     #, categoryChoices = categoryChoices
promoter <- function (converted,selectOrg, radio){
  idNotRecognized = as.data.frame("ID not recognized!")
  
  if(is.null(converted) ) 
    return(idNotRecognized) # no ID
  
  querySet <- converted$IDs;
  
  if(length(querySet) == 0) 
    return(idNotRecognized )
  ix = grep(converted$species[1,1],motifFiles)

  # If selected species is not the default "bestMatch", use that species directly
  if(selectOrg != speciesChoice[[1]]) {
    ix = grep(findSpeciesById(selectOrg)[1,1], motifFiles )
  }
    
  ix1 =grep(as.character(radio),motifFiles[ix]) # match 300bp or 600bp
  if(length(ix1) >0) ix = ix[ix1]   # if 600 is not found, use 300bp
  if (length(ix) == 0 ) {
    return(as.data.frame("No matching motif file found") )
  } 
  else {
    
    if(length(ix) > 1)  # if only one file
      return(as.data.frame("Multiple geneInfo file found!") )
  
    motifs <- dbConnect(sqlite,motifFiles[ix]) # makes a new file
  
    sqlQuery = paste( " select * from scores where row_names IN ('", paste(querySet,collapse="', '"),"')" ,sep="")
    result <- dbGetQuery( motifs, sqlQuery  )
    if( dim(result)[1] ==0) {
      return(list( x=as.data.frame("No matching species or gene ID file!" )) )
    }
    row.names(result) <- result$row_names; result <- result[,-1]
    TFstat <- as.data.frame( cbind(apply(result,2,mean),apply(result,2,sd) ) )
    colnames(TFstat) = c("scoreMean1","scoreSD1" )
    rownames(TFstat) = toupper( colnames(result) )
  
    TFs <- dbGetQuery(motifs, "select ID,TF_Name,Family_Name,DBID,Motif_ID,coreMotif,memo,nGenes,scoreSD,scoreMean from  TF_Information ")
    dbDisconnect(motifs)
    TFs$ID <- toupper(TFs$ID)
  
    TFs <- merge(TFs, TFstat, by.x = 'ID', by.y='row.names')
    TFs <- TFs[!is.na(TFs$scoreSD) ,]  #some TFs return NA -Inf
    n1 = dim(result)[1] # number of genes in query set
    TFs$scoreMean2 <- (TFs$scoreMean * TFs$nGenes - TFs$scoreMean1 *n1)/(TFs$nGenes - n1)
    #SD2 needs to be adjusted too, but ignored for now. use overall SD2
    # t test unequal variance statistic
    TFs$t <- (TFs$scoreMean1-TFs$scoreMean2)/ sqrt( TFs$scoreSD1^2/n1 + TFs$scoreSD^2/TFs$nGenes   )
    # degree of freedom
    TFs$df <- ( TFs$scoreSD1^2/n1 + TFs$scoreSD^2/TFs$nGenes)^2 / ((TFs$scoreSD1^2/n1)^2/(n1-1) +   (TFs$scoreSD^2/TFs$nGenes)^2/(TFs$nGenes-1))
    TFs$pVal =1-pt(TFs$t,df = TFs$df)  # t distribution
    TFs$FDR = p.adjust(TFs$pVal,method="fdr")
    TFs <- TFs[order(TFs$pVal) ,]
    TFs$scoreDiff = round(TFs$scoreMean1 - TFs$scoreMean2,0)
    #TFs <- TFs[order(-TFs$scoreDiff) ,]
  
    # does this transcription factor gene in this cluster?
    ix <- match(toupper( TFs$DBID), querySet) # assuming the DBID column in cisbp are ensembl gene ids
    TFs$note = ""
    if(sum(!is.na(ix)) >0) {
      TFs$note[which(!is.na(ix))] <- "* Query Gene"
    }
    TFs <- subset(TFs, FDR<0.25, select=c(coreMotif,TF_Name,Family_Name, pVal,FDR,scoreDiff, note ) )
    colnames(TFs) =c("Enriched motif in promoter", "TF","TF family","P val.","FDR","Score","Note"   )
    if(dim(TFs)[1] >30 ) 
      TFs <- TFs[1:30,]
    if(dim(TFs)[1] ==0) 
      return(as.data.frame("No significant TF binding motif detected.") ) else
    return( TFs )
   }
}


shinyServer(
  function(input, output,session){
    options(warn=-1)

    observe({  updateSelectInput(session, "selectOrg", choices = speciesChoice )      })
	#observe({  updateSelectInput(session, "selectGO", choices = list("All available gene sets" = "All", "GO Biological Process" = "GOBP","GO Molecular Function" = "GOMF","GO Cellular Component" = "GOCC",
     #                           "KEGG metabolic pathways" = "KEGG"), selected = "All")     })

	# this defines an reactive object that can be accessed from other rendering functions
	converted <- reactive({
	  if (input$goButton == 0)    return()

      isolate( {  convertID(input$input_text,input$selectOrg );})

	} )

	geneInfoLookup <- reactive({
	  if (input$goButton == 0)    return()

      isolate( { geneInfo(converted(),input$selectOrg )  }) # uses converted gene ids thru converted() call

	} )

	significantOverlaps <- reactive({
	  if (input$goButton == 0) return()
    isolate( {
  	  #gene info is passed to enable lookup of gene symbols
  	  tem = geneInfoLookup(); tem <- tem[which( tem$Set == "List"),]
  	  FindOverlap( converted(), tem, input$selectGO,input$selectOrg,input$minFDR )
    })
	})


  output$species <-renderTable({
    if (input$goButton == 0)    return()
    isolate( {  #tem <- convertID(input$input_text,input$selectOrg );
	  	  withProgress(message="Converting gene IDs", {
                  tem <- converted()
			incProgress(1, detail = paste("Done"))	  })

				  if( is.null(tem)) {as.data.frame("ID not recognized.")} else {
	              tem$speciesMatched }

      }) # avoid showing things initially
    }, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)


promoterData <-reactive({
  if (input$goButton == 0)    return()
	  tem = input$radio
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
	tem = input$radio
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
      isolate( {
          tem <- converted();
          tem2 <- geneInfoLookup()
				  if( is.null(tem)) {as.data.frame("ID not recognized.")} else {
  				  if(dim(tem2)[1] == 1) { tem$conversionTable }
            else { # if gene info is not available
              merged <- merge(tem$conversionTable,tem2,by='ensembl_gene_id')
    				  merged <- subset(merged,select=c(User_input,ensembl_gene_id,symbol,gene_biotype,Species,chromosome_name,start_position  ))

    				  tem3 <- as.data.frame(tem$originalIDs); colnames(tem3) = "User_input"
    				  merged <- merge(merged, tem3, all=T)
    				  merged$ensembl_gene_id[which(is.na(merged$ensembl_gene_id))] <- "Not mapped"
    				  merged <- merged[order(merged$ensembl_gene_id, decreasing =T),]
    				  merged <- merged[order(merged$gene_biotype),]
    				  merged <- merged[order( as.numeric(merged$chromosome_name)),]
    				  #merged <- merged[order( merged$start_position),]
    				  merged$start_position = merged$start_position/1e6
    				  colnames(merged) <- c("User ID", "Ensembl Gene ID", "Symbol",
    				     "Gene Type", "Species", "Chr", "Position (Mbp)" )
    				  i = 1:dim(merged)[1]
    				  merged = cbind(i,merged)
				  }
        }
      }) # avoid showing things initially
    })

output$conversionTable <-renderTable({
      if (input$goButton == 0)    return()   # still have problems when geneInfo is not found!!!!!
      isolate( {
	  conversionTableData()

      }) # avoid showing things initially
    }, digits = 4,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)


output$downloadGeneInfo <- downloadHandler(
	  filename = function() {"geneInfo.csv"},
		content = function(file) {
			write.csv(conversionTableData(), file, row.names=FALSE)
	    }
  )
output$table <-renderTable({
      if (input$goButton == 0  )    return()
	  myMessage = "Those genes seem interesting! Let me see what I can do.
	   I am comparing your query genes to all 150+ types of IDs across 111 species.
	  This can take up to 3 years. "
	  withProgress(message= sample(quotes,1),detail=myMessage, {
		tem <- significantOverlaps();
	  incProgress(1, detail = paste("Done"))	  })
	  if(dim(tem$x)[2] ==1 ) tem$x else tem$x[,1:4]  # If no significant enrichment found x only has 1 column.
    }, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)

output$downloadEnrichment <- downloadHandler(
	  filename = function() {"enrichment.csv"},
		content = function(file) {
			write.csv(significantOverlaps()$x, file, row.names=FALSE)
	    }
  )

    # this updates geneset categories based on species and file
    output$selectGO1 <- renderUI({
      if (input$goButton == 0 | class(significantOverlaps())=="data.frame"  )
       { selectInput("selectGO", label = NULL, # h6("Funtional Category"),
                  choices = list("All available gene sets" = "All", "GO Biological Process" = "GOBP","GO Molecular Function" = "GOMF","GO Cellular Component" = "GOCC",
                                "KEGG metabolic pathways" = "KEGG"), selected = "GOBP")  }	 else {

	  selectInput("selectGO", label=NULL,choices=significantOverlaps()$categoryChoices,selected="GOBP" )   }
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
    output$text1 <- renderText({
      #"Example mouse genes: Hus1 Rad1 Trp63 Trp73 Usp28 Rad9b Fanci Hus1b Cdk1 Cry1 D7Ertd443e Chek1 Foxo4 Zak Pea15a Mapkapk2 Brca1 E2f1"
	  "Example mouse genes: Hus1 Rad1 Trp63 Trp73 Usp28 Rad9b Fanci Hus1b Cdk1 Cry1 D7Ertd443e Chek1 Foxo4 Zak Pea15a Mapkapk2 Brca1 Taok1 Cdk5rap3 Ddx39b Mdm2 Fzr1 Rad17 Prkdc Cdkn1a Cdc5l Wac Thoc1 Prpf19 Rad9a Pidd1 Atrip Uimc1 Nek6 Atf2 E2f1 Nbn Rpa2 Rint1 Clock Chek2 Casp2 Blm Plk1 Brcc3 Hinfp Fem1b Tipin Atr Cdc14b Rfwd3 Ccar2 Foxn3 Atm Thoc5 Nek11 Fam175a Brsk1 Plk5 Rps27l Ints7 Dtl Tiprl Rbbp8 Clspn Cradd Rhno1 Bre Trp53 Taok2 Taok3 Ccnd1 Sox4 Msh2 Xpc Rad9a Rnaseh2b Fbxo4 Syf2 Cul4a Nek1 Mre11a Pml Ptpn11 Zfp830 Gigyf2 Mapk14 Bcat1 Fbxo31 Babam1 Cep63"
     # Ddx39b excluded as it has multiple ensembl gene ids?
	})

    output$genomePlot <- renderPlot({
	  if (input$goButton == 0  )    return()
	  isolate( {
       x = geneInfoLookup()
       converted1 = converted()

   	   #chromosomes
	   if((sum(!is.na( x$chromosome_name) ) >= minGenes && length(unique(x$chromosome_name) ) > 2 ) && length(which(x$Set == "List") ) > minGenes )
	   {
		   freq = table( x$chromosome_name,x$Set );
		   freq <- as.matrix(freq[which(nchar(row.names(freq))<3   ),])# remove unmapped chromosomes
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
		   text(1,y1+.2, paste(round(chrLength/1e6,0),"Mbp",sep=""),cex=2)
		   segments(0,y1+.01,1,y1+.01,col="blue")
           sapply(1:dim(a1)[1], function (i) text(a1$start_position[i], y1+.03, "|"))
		   if(dim(a1)[1] <100 && freq[i,2] >0)  # if more genes, do not show symbol
               sapply(1:dim(a1)[1], function (i) text(a1$start_position[i], y1, a1$preferedIDs[i], offset=0,srt=90, pos=2,cex=1.5))

    	}}

			}
		}


	   } )# isolate
	   }, height = 3000, width = 1000)


  output$genePlot <- renderPlot({
	   	  if (input$goButton == 0  )    return()
		    isolate( {
	  withProgress(message="Ploting gene characteristics", {
       x = geneInfoLookup()
	   x2 = x[which(x$gene_biotype == "protein_coding"),]  # only coding for some analyses
     if(dim(x)[1]>=minGenes) # only making plots if more than 20 genes
       { # only plot when there 10 genes or more   # some columns have too many missing values
	   par(mfrow=c(5,2))
   	   #chromosomes
	   if( sum(!is.na( x$chromosome_name) ) >= minGenes && length(unique(x$chromosome_name) ) > 2 && length(which(x$Set == "List") ) > minGenes )
	   {
		   freq = table( x$chromosome_name,x$Set );
		   freq <- as.matrix(freq[which(nchar(row.names(freq))<3   ),])# remove unmapped chromosomes
		   if(dim(freq)[2] >1 && dim(freq)[1]>1 ) { # some organisms do not have fully seuqence genome: chr. names: scaffold_99816
				Pval = chisq.test(freq)$p.value
				sig = paste("Distribution of query genes on chromosomes \nChi-squared test P=",formatC(Pval, digits=2, format="G") )
			   if( Pval <PvalGeneInfo)  sig = paste(sig," ****************" )
			   freq <- freq[order( as.numeric(row.names(freq) )), ]
				freq[,1] <- freq[,1] *colSums(freq)[2]/colSums(freq)[1] # expected
				freq = freq[,c(2,1)] # reverse order
			   barplot(t(freq), beside=TRUE,las=3,col=c("red","lightgrey"), ylab="Number of Genes",main= sig,xlab =c("Chromosome") )
		   legend("topright", c("List","Expected"), pch=15, col=c("red","lightgrey"),bty="n")
			}
		}
		incProgress(1/8)
	   # gene type
	    if( sum(!is.na( x$gene_biotype) ) >= minGenes && length(unique(x$gene_biotype) ) > 2  && length(which(x$Set == "List") ) > minGenes ) {
	 	freq = table( x$gene_biotype,x$Set );
		freq <- as.matrix(freq[which( freq[,1]/colSums(freq)[1] >.01),])
	   if(dim(freq)[2] >1 && dim(freq)[1]>1 ) {
	    Pval = chisq.test(freq)$p.value
		sig=paste("Distribution by gene type \nChi-squared test P=",formatC(Pval, digits=2, format="G") )
        if( Pval <PvalGeneInfo)  sig = paste(sig," ****************" )
		freq <- freq[order(    freq[,1], decreasing=T), ]
		freq[,1] <- freq[,1] *colSums(freq)[2]/colSums(freq)[1]
		tem = gsub("protein_coding","Coding",rownames(freq));
		tem =gsub("processed_pseudogene","proc_pseudo",tem)
	    tem =gsub("processed","proc",tem); #row.names(freq)= tem
		par(mar=c(12,4.1,4.1,2.1))
		freq = freq[,c(2,1)] # reverse order
        barplot(t(freq), beside=TRUE,las=2,col=c("red","lightgrey"), ylab="Number of Genes",
	      main= sig)
	    legend("topright", c("List","Expected"), pch=15, col=c("red","lightgrey"),bty="n")
		} }
		incProgress(1/8)

        # N. exons

		if( sum(!is.na( x2$nExons) ) >= minGenes && length(unique(x2$nExons) ) > 2  && length(which(x2$Set == "List") ) > minGenes ) {
		freq = table( x2$nExons,x2$Set );
		freq <- as.matrix(freq[which( freq[,1]/colSums(freq)[1] >.02),])
	    if(dim(freq)[2] >1 && dim(freq)[1]>1 ) {
	    Pval = chisq.test(freq)$p.value
		sig=paste("Number of exons (coding genes only) \nChi-squared test P=",formatC(Pval, digits=2, format="G") )
        if( Pval <PvalGeneInfo)  sig = paste(sig," ****************" )
		freq <- freq[order(    freq[,1], decreasing=T), ]
		freq[,1] <- freq[,1] *colSums(freq)[2]/colSums(freq)[1]
		freq = freq[,c(2,1)] # reverse order
        barplot(t(freq), beside=TRUE,las=2,col=c("red","lightgrey"), ylab="Number of Genes",
	      main= sig ,xlab =c("Number of exons"))
	    legend("topright", c("List","Expected"), pch=15, col=c("red","lightgrey"),bty="n")
		}}
		incProgress(1/8)

		#Transcript count
		if( sum(!is.na( x2$transcript_count) ) >= minGenes && length(unique(x2$transcript_count) ) > 2  && length(which(x2$Set == "List") ) > minGenes ) {
		freq = table( x2$transcript_count,x2$Set );
		freq <- as.matrix(freq[which( freq[,1]/colSums(freq)[1] >.02),])
		if(dim(freq)[2] >1 && dim(freq)[1]>1 ) {
	    Pval = chisq.test(freq)$p.value
		sig=paste("Number of transcript isoforms per coding gene \nChi-squared test P=",formatC(Pval, digits=2, format="G"))
        if( Pval <PvalGeneInfo)  sig = paste(sig," ****************" )
		freq <- freq[order(    freq[,1], decreasing=T), ]
		freq[,1] <- freq[,1] *colSums(freq)[2]/colSums(freq)[1]
		freq = freq[,c(2,1)] # reverse order
        barplot(t(freq), beside=TRUE,las=2,col=c("red","lightgrey"), ylab="Number of Genes",
	      main= sig,xlab =c("Number of transcripts per gene") )
	    legend("topright", c("List","Expected"), pch=15, col=c("red","lightgrey"),bty="n")
       } }
	   incProgress(1/8)

	  if( sum(!is.na( x2$cds_length) ) >= minGenes && length(unique(x2$cds_length) ) > 2 && length(which(x2$Set == "List") ) > minGenes) {
	   Pval = t.test(cds_length~Set, data=x2 )$p.value
	   sig = paste("Coding Sequence length \n T-test P=",formatC(Pval, digits=2, format="G"),sep="")
	   if( Pval <PvalGeneInfo)  sig = paste(sig," ****************" )
	   boxplot(cds_length~Set,ylim=c(0,4000), data=x,main=sig)
       }
	   incProgress(1/8)

	  if( sum(!is.na( x2$transcript_length) ) >= minGenes && length(unique(x2$transcript_length) ) > 2 && length(which(x2$Set == "List") ) > minGenes ) {
       Pval = t.test(transcript_length~Set, data=x2 )$p.value
	   sig = paste("Transcript length (coding genes only)\nT-test  P=",formatC(Pval, digits=2, format="G"),sep="")
	   if( Pval <PvalGeneInfo)  sig = paste(sig," ****************" )
	   boxplot(transcript_length~Set,ylim=c(0,6000), data=x,main=sig)
	  }
	  incProgress(1/8)

	  if( sum(!is.na( x2$genomeSpan) ) >= minGenes && length(unique(x2$genomeSpan) ) > 2 && length(which(x2$Set == "List") ) > minGenes ) {
       Pval = t.test(genomeSpan~Set, data=x2 )$p.value
	   sig = paste("Genome span (coding genes only) T-test \n P=",formatC(Pval, digits=2, format="G"),sep="")
	   if( Pval <PvalGeneInfo)  sig = paste(sig," ****************" )
	   boxplot(genomeSpan~Set,ylim=c(0,80000), data=x,main=sig)
	   }
	   incProgress(1/8)

	  if( sum(!is.na( x2$FiveUTR) ) >= minGenes && length(unique(x2$FiveUTR) ) > 2 && length(which(x2$Set == "List") ) > minGenes ) {
       Pval = t.test(FiveUTR~Set, data=x2 )$p.value
	   sig = paste("5' UTR length (coding genes only)\n T-test P=",formatC(Pval, digits=2, format="G"),sep="")
	   if( Pval <PvalGeneInfo)  sig = paste(sig," ****************" )
	   boxplot(FiveUTR~Set, data=x,ylim=c(0,400),main=sig)
	}
	incProgress(1/8)

	if( sum(!is.na( x2$ThreeUTR) ) >= minGenes && length(unique(x2$ThreeUTR) ) > 2 && length(which(x2$Set == "List") ) > minGenes ) {
       Pval = t.test(ThreeUTR~Set, data=x2 )$p.value
	   sig = paste("3' UTR length (coding genes only) \n T-test P=",formatC(Pval, digits=2, format="G"),sep="")
	   if( Pval <PvalGeneInfo)  sig = paste(sig," ****************" )
	   boxplot(ThreeUTR~Set, ylim=c(0,2000),data=x,main=sig)
	  }
	  incProgress(1/8)

	  if( sum(!is.na( x2$percentage_gc_content) ) >= minGenes && length(unique(x2$percentage_gc_content) ) > 2 && length(which(x2$Set == "List") ) > minGenes ) {
       Pval = t.test(percentage_gc_content~Set, data=x2 )$p.value
	   sig = paste("%GC content (coding genes only)\n T-test P=",formatC(Pval, digits=2, format="G"),sep="")
	   if( Pval <PvalGeneInfo)  sig = paste(sig," ****************" )
	   boxplot(percentage_gc_content~Set, ylim = c(20,80), data=x,main=sig)
	   }
     } # if minGenes
	 incProgress(1/8, detail = paste("Done"))	  })
	 }) #isolate
    }, height = 1500, width = 600)

    }
)
