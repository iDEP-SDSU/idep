library(shiny)
library(RSQLite)
library(ggplot2)
#library(grid)
library(gridExtra)
library(plotly)
library(reshape2)
library(visNetwork)

# relative path to data files
datapath = "../../data/data103/"   # production server

Min_overlap <- 2
minSetSize = 3;
mappingCoverage = 0.60 # 60% percent genes has to be mapped for confident mapping
mappingEdge = 0.5  # Top species has 50% more genes mapped
maxTerms =30 # max number of enriched terms; no longer used
PvalGeneInfo = 0.05; minGenes = 10 # min number of genes for ploting
PvalGeneInfo1 = 0.01
PvalGeneInfo2 = 0.001
pdf(NULL) # this prevents error Cannot open file 'Rplots.pdf'
ExampleGeneList=
"Hus1 Rad1 Tp63 Tp73 Usp28 Rad9b Fanci Hus1b 
Cdk1 Cry1 D7Ertd443e Chek1 Foxo4 Zak Pea15a 
Mapkapk2 Brca1 Taok1 Cdk5rap3 Ddx39b Mdm2 Fzr1 
Rad17 Prkdc Cdkn1a Cdc5l Wac Thoc1 Prpf19 Rad9a
Pidd1 Atrip Uimc1Nek6 Atf2 E2f1 Nbn Rpa2 Rint1 
Clock Chek2 Casp2 Blm Plk1 Brcc3 Hinfp Fem1b 
Tipin Atr Cdc14b Rfwd3 Ccar2 Foxn3 Atm Thoc5  
Rps27l Ints7 Dtl Tiprl Rbbp8 Clspn Cradd Rhno1  
Sox4 Msh2 Xpc Rad9a Rnaseh2b Fbxo4 Syf2 Cul4a 
Gigyf2 Mapk14 Bcat1 Fbxo31 Babam1 Cep63 Ccnd1
Nek11 Fam175a Brsk1 Plk5 Bre Tp53 Taok2 Taok3 
Nek1 Mre11a Pml Ptpn11 Zfp830 
"
ExampleGeneList= "ENSG00000078900
ENSG00000117614
ENSG00000117748
ENSG00000092853
ENSG00000143155
ENSG00000162889
ENSG00000143493
ENSG00000143476
ENSG00000095002
ENSG00000115966
ENSG00000204120
ENSG00000154767
ENSG00000164053
ENSG00000114670
ENSG00000182923
ENSG00000175054
ENSG00000073282
ENSG00000134852
ENSG00000137601
ENSG00000113456
ENSG00000151876
ENSG00000152942
ENSG00000188996
ENSG00000124766
ENSG00000198563
ENSG00000112062
ENSG00000124762
ENSG00000096401
ENSG00000136273
ENSG00000135249
ENSG00000106144
ENSG00000158941
ENSG00000253729
ENSG00000104320
ENSG00000081377
ENSG00000095787
ENSG00000170312
ENSG00000177595
ENSG00000110107
ENSG00000172613
ENSG00000110092
ENSG00000149311
ENSG00000048028
ENSG00000172273
ENSG00000149554
ENSG00000171792
ENSG00000060982
ENSG00000135679
ENSG00000169372
ENSG00000008405
ENSG00000151164
ENSG00000179295
ENSG00000135090
ENSG00000136104
ENSG00000139842
ENSG00000053254
ENSG00000185088
ENSG00000075131
ENSG00000169018
ENSG00000140464
ENSG00000140525
ENSG00000197299
ENSG00000166851
ENSG00000149930
ENSG00000168411
ENSG00000103264
ENSG00000141510
ENSG00000160551
ENSG00000012048
ENSG00000108465
ENSG00000079134
ENSG00000101773
ENSG00000185988
ENSG00000105325
ENSG00000105393
ENSG00000160469
ENSG00000101412
ENSG00000183765
ENSG00000100296
ENSG00000184481
ENSG00000185515"

# Wrapping long text by adding \n 
#  "Mitotic DNA damage checkpoint"  --> "Mitotic DNA damage\ncheckpoint"
# https://stackoverflow.com/questions/7367138/text-wrap-for-plot-titles
wrap_strings <- function( vector_of_strings, width = 30 ) { 
  as.character( sapply( vector_of_strings, FUN=function(x) 
  { paste(strwrap(x, width = width), collapse = "\n")}) )
}

# function to increase vertical spacing between legend keys
# @clauswilke https://stackoverflow.com/questions/11366964/is-there-a-way-to-change-the-spacing-between-legend-items-in-ggplot2
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)

  grid::rectGrob(
    width = grid::unit(0.6, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}
# register new key drawing function, 
# the effect is global & persistent throughout the R session
GeomBar$draw_key = draw_key_polygon3

# find peak values in density plots
# for adding annotation texts
# http://ianmadd.github.io/pages/PeakDensityDistribution.html
densMode <- function(x){
    td <- density(x, na.rm = TRUE)
    maxDens <- which.max(td$y)
    list(x=td$x[maxDens], y=td$y[maxDens])
}

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

STRING10_species = read.csv( paste0(datapath,"data_go/STRING11_species.csv") )

# Create a list for Select Input options
orgInfo <- dbGetQuery(convert, paste("select distinct * from orgInfo " ))
orgInfo <- orgInfo[order(orgInfo$name),]
annotatedSpeciesCounts <- sort( table(orgInfo$group) ) # total species, Ensembl, Plants, Metazoa, STRINGv10
speciesChoice <- setNames(as.list( orgInfo$id ), orgInfo$name2 )
# add a defult element to list    # new element name       value
speciesChoice <- append( setNames( "BestMatch","Best matching species"), speciesChoice  )
# move one element to the 2nd place
move2 <- function(i) c(speciesChoice[1],speciesChoice[i],speciesChoice[-c(1,i)])
i= which( names(speciesChoice) == "Vitis vinifera"); speciesChoice <- move2(i)
i= which( names(speciesChoice) == "Oryza sativa Japonica Group"); speciesChoice <- move2(i)
i= which( names(speciesChoice) == "Oryza sativa Indica Group"); speciesChoice <- move2(i)
i= which( names(speciesChoice) == "Glycine max"); speciesChoice <- move2(i)
i= which( names(speciesChoice) =="Zea mays"); speciesChoice <- move2(i)
i= which(names(speciesChoice) =="Arabidopsis thaliana"); speciesChoice <- move2(i)
i= which(names(speciesChoice) == "Saccharomyces cerevisiae"); speciesChoice <- move2(i)
i= which(names(speciesChoice)  == "Caenorhabditis elegans"); speciesChoice <- move2(i)
i= which(names(speciesChoice)  == "Drosophila melanogaster"); speciesChoice <- move2(i)
i= which(names(speciesChoice) =="Dog"); speciesChoice <- move2(i)
i= which(names(speciesChoice) =="Macaque"); speciesChoice <- move2(i)
i= which(names(speciesChoice) =="Chicken"); speciesChoice <- move2(i)
i= which(names(speciesChoice) =="Pig"); speciesChoice <- move2(i)
i= which(names(speciesChoice) =="Zebrafish" ); speciesChoice <- move2(i)
i= which(names(speciesChoice) == "Cow" ); speciesChoice <- move2(i)
i= which(names(speciesChoice) == "Rat" ); speciesChoice <- move2(i)
i= which(names(speciesChoice) == "Mouse"); speciesChoice <- move2(i)
i= which(names(speciesChoice) == "Human"); speciesChoice <- move2(i)

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

#Homo sapies --> hsapiens
shortSpeciesNames <- function(tem){
	 tem2 = strsplit(as.character(tem)," " ) 	   
	 return( tolower( paste0(substr(tem2[[1]][1],1,1), tem2[[1]][2]  ) ) )
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
 query <- gsub("\"|\'","",query) 
 	# remove " in gene ids, mess up SQL query				
	# remove ' in gene ids				
	# |\\.[0-9] remove anything after A35244.1 -> A35244  
	#  some gene ids are like Glyma.01G002100
	
  querySet <- cleanGeneSet( unlist( strsplit( toupper(query),'\t| |\n|\\,' )  ) )
 	
    if( selectOrg == "BestMatch") { # query all species
	  querySTMT <- paste( "select distinct id,ens,species from mapping where id IN ('", paste(querySet,collapse="', '"),"')",sep="")
    } else {  # organism has been selected query specific one
 	  querySTMT <- paste( "select distinct id,ens,species from mapping where species = '",selectOrg,
                          "' AND id IN ('", paste(querySet,collapse="', '"),"')",sep="")    
    }
	result <- dbGetQuery(convert, querySTMT)
  if( dim(result)[1] == 0  ) return(NULL)
  

  if(selectOrg == speciesChoice[[1]]) {
    comb = paste( result$species,result$idType)
    sortedCounts = sort( table(comb ),decreasing=T)
    # Try to use Ensembl instead of STRING-db genome annotation
    if( sortedCounts[1] <= sortedCounts[2] *1.1  # if the #1 species and #2 are close
         && as.numeric(names(sortedCounts[1])) > sum( annotatedSpeciesCounts[1:3])  # 1:3 are Ensembl species
         && as.numeric(names( sortedCounts[2] )) < sum( annotatedSpeciesCounts[1:3])    ) { # and #2 come earlier (ensembl) than #1
      tem <- sortedCounts[2]
      sortedCounts[2] <- sortedCounts[1]
      names(sortedCounts)[2] <- names(sortedCounts)[1]
       sortedCounts[1] <- tem
      names(sortedCounts)[1] <- names(tem)    
    } 
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

  result <- result[which(!duplicated(result[,1]) ),] # remove duplicates in query gene ids 
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
FindOverlap <- function (converted, gInfo, GO, selectOrg, minFDR, input_maxTerms, convertedB=NULL, gInfoB=NULL) {
  idNotRecognized = list(x=as.data.frame("ID not recognized!"),
                         groupings= as.data.frame("ID not recognized!")  )
  if(is.null(converted) ) return(idNotRecognized) # no ID
  querySet <- converted$IDs;
  if(length(querySet) == 0) return(idNotRecognized )

  ix = grep(converted$species[1,1],gmtFiles)
  totalGenes <- converted$species[1,7]

  errorMessage = list(x=as.data.frame("Annotation file cannot be found"),
                         groupings= as.data.frame("Annotation file cannot be found")  )
  if (length(ix) == 0 ) {return( errorMessage )}

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

  sqlQuery = paste( " select distinct gene,pathwayID from pathway where gene IN ('", paste(querySet, collapse="', '"),"')" ,sep="")

  if( GO != "All") sqlQuery = paste0(sqlQuery, " AND category ='",GO,"'")
  result <- dbGetQuery( pathway, sqlQuery  )
  if( dim(result)[1] ==0) {return(list( x=as.data.frame("No matching species or gene ID file!" )) )}

   # given a pathway id, it finds the overlapped genes, symbol preferred
  sharedGenesPrefered <- function(pathwayID) {
    tem <- result[which(result[,2]== pathwayID ),1]
    ix = match(tem, converted$conversionTable$ensembl_gene_id) # convert back to original
    tem2 <- unique( converted$conversionTable$User_input[ix] )
    #    if(length(unique(gInfo$symbol) )/dim(gInfo)[1] >.7  ) # if 70% genes has symbol in geneInfo
    #	{ ix = match(tem, gInfo$ensembl_gene_id);
    #	  tem2 <- unique( gInfo$symbol[ix] )      }
    return( paste( tem2 ,collapse=" ",sep="") )
	}
  
  x0 = table(result$pathwayID)

  x0 = as.data.frame( x0[which(x0>=Min_overlap)] )# remove low overlaps

  errorMessage = list(x=as.data.frame("Too few genes."),
                         groupings= as.data.frame("Too few genes.")  )
  if(dim(x0)[1] <= 5 ) return(errorMessage) # no data
  colnames(x0)=c("pathwayID","overlap")
  pathwayInfo <- dbGetQuery( pathway, paste( " select distinct id,n,Description from pathwayInfo where id IN ('",
						paste(x0$pathwayID,collapse="', '"),   "') ",sep="") )
  x = merge(x0,pathwayInfo, by.x='pathwayID', by.y='id')
  
  
  #Background genes----------------------------------------------------
  if(!is.null(convertedB) && 
     !is.null(gInfoB) && 
     length( convertedB$IDs) < 30000) { # if more than 30k genes, ignore background genes.
        querySetB <- convertedB$IDs;    
        # if background and selected genes matches to different organisms, error
        if( length( intersect( querySetB, querySet ) ) == 0 )    # if none of the selected genes are in background genes
          return(list( x=as.data.frame("None of the selected genes are in the background genes!" )) )
        
        querySetB <- unique( c( querySetB, querySet ) )  # just to make sure the background set includes the query set
        
        sqlQueryB = paste( " select distinct gene,pathwayID from pathway where gene IN ('", paste(querySetB, collapse="', '"),"')" ,sep="")    
        
        if( GO != "All") sqlQuery = paste0(sqlQuery, " AND category ='",GO,"'")
        resultB <- dbGetQuery( pathway, sqlQueryB  )
        if( dim(resultB)[1] ==0) {return(list( x=as.data.frame("No matching species or gene ID file!" )) )}    
        xB = table(resultB$pathwayID)
        rm(resultB)
        xB = as.data.frame( xB)
        colnames(xB)=c("pathwayID","overlapB")
        x2 = merge(x, xB, by='pathwayID', all.x = TRUE)       
        
        x$Pval=phyper(x2$overlap - 1,
                      length(querySet),
                      length(querySetB) - length(querySet),   
                      as.numeric(x2$overlapB), # use the number of genes in background set
                      lower.tail=FALSE ); 
        
      }   else { # original version without background genes
          x$Pval <- phyper(x$overlap - 1,
                        length(querySet),
                        totalGenes - length(querySet),   
                        as.numeric(x$n), 
                        lower.tail=FALSE );
          }
  
  x$FDR <- p.adjust(x$Pval, method="fdr")
  x <- x[order(x$FDR), ]  # sort according to FDR


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
  if(dim(x)[1] > as.integer(input_maxTerms) ) x = x[ 1:as.integer(input_maxTerms), ]
  x= cbind(x,sapply( x$pathwayID, sharedGenesPrefered ) )
  colnames(x)[7]= "Genes"
  x$n <- as.numeric(x$n) # convert total genes from character to numeric 10/21/19
  x <- subset(x,select = c(FDR,overlap,n,description,Genes) )
  colnames(x) = c("Enrichment FDR", "Genes in list", "Total genes","Functional Category","Genes"  )
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


mycolors = sort(rainbow(20))[c(1,20,10,11,2,19,3,12,4,13,5,14,6,15,7,16,8,17,9,18)] # 20 colors for kNN clusters

# a program for ploting enrichment results by highlighting the similarities among terms
# must have columns: Direction, adj.Pval   Pathways Genes
#  Direction	adj.Pval	nGenes	Pathways		Genes
#Down regulated	3.58E-59	131	Ribonucleoprotein complex biogenesis	36	Nsun5 Nhp2 Rrp15 
#Down regulated	2.55E-57	135	NcRNA metabolic process	23	Nsun5 Nhp2 Rrp15 Emg1 Ddx56 Rsl1d1 enrichmentPlot <- function( enrichedTerms){
# Up or down regulation is color-coded
# gene set size if represented by the size of marker
enrichmentPlot <- function( enrichedTerms, rightMargin=33) {
  if(class(enrichedTerms) != "data.frame") return(NULL)
  if(nrow(enrichedTerms) <=1 ) return(NULL)  # only one term or less
  library(dendextend) # customizing tree
  
  geneLists = lapply(enrichedTerms$Genes, function(x) unlist( strsplit(as.character(x)," " )   ) )
  names(geneLists)= enrichedTerms$Pathways

  # compute overlaps percentage--------------------

  n = length(geneLists)
  w <- matrix(NA, nrow = n, ncol = n)
# compute overlaps among all gene lists
    for (i in 1:n) {
        for (j in i:n) {
            u <- unlist(geneLists[i])
            v <- unlist(geneLists[j])
            w[i, j] = length(intersect(u, v))/length(unique(c(u,v)))
        }
    }
# the lower half of the matrix filled in based on symmetry
    for (i in 1:n) 
        for (j in 1:(i-1)) 
            w[i, j] = w[j,i] 
 

 # compute overlaps P value---------------------
  if(0) {
 total_elements = 30000
  n = length(geneLists)
  w <- matrix(rep(0,n*n), nrow = n, ncol = n)
# compute overlaps among all gene lists
    for (i in 1:n) {
        for (j in (i+1):n) {
            u <- unlist(geneLists[i])
            v <- unlist(geneLists[j])
            xx= length( intersect(u, v) )
			if(xx == 0)
				next;
			mm = length(u)
			nn <- total_elements - mm	
			kk = length(v)
			w[i,j] = -sqrt( -phyper(xx-1,mm,nn,kk, lower.tail=FALSE,log.p = TRUE ));
			
        }
    }
	

# the lower half of the matrix filled in based on symmetry
    for (i in 1:n) 
        for (j in 1:(i-1)) 
            w[i, j] = w[j,i] 
			
	# w =  w-min(w) 			
	# for( i in 1:n) 		w[i,i] = 0;
 
 }

  Terms = paste( sprintf("%-2.1e",as.numeric(enrichedTerms$adj.Pval)), 
				names(geneLists))
  rownames(w) = Terms
  colnames(w) = Terms
  par(mar=c(0,0,1,rightMargin)) # a large margin for showing 

  dend <- as.dist(1-w) %>%
	hclust (method="average") 
  ix = dend$order # permutated order of leaves

  leafType= as.factor( gsub(" .*","", enrichedTerms$Direction[ix] ) )
  #if(length(unique(enrichedTerms$Direction)  ) <=2 )
  if( max( nchar(enrichedTerms$Direction[ix] )) >= 1)   # if "Up regulated or Downregulated"; not "A", "B"
	#leafColors = c("green","red")  else  # mycolors # k-Means
	leafColors = mycolors[1:2] else
		{ 	# convert c("B","D","E") to c(2, 4, 5)
			#leafType= as.factor( gsub(" .*","", enrichedTerms$Direction[ix] ) )
			leafType= match(gsub(" .*","", enrichedTerms$Direction[ix] ), toupper(letters)   )   
			
			leafColors = mycolors 
		}
  #leafSize = unlist( lapply(geneLists,length) ) # leaf size represent number of genes
  #leafSize = sqrt( leafSize[ix] )  
  leafSize = -log10(as.numeric( enrichedTerms$adj.Pval[ix] ) ) # leaf size represent P values
  leafSize = .9*(leafSize-min(leafSize))/(max( leafSize )-min(leafSize)+1e-50) + .1   # scale more aggressively
  # leafSize = 1.*(leafSize)/max( leafSize ) + .1   # ratio scaling, less agressive
	dend %>% 
	as.dendrogram(hang=-1) %>%
	set("leaves_pch", 19) %>%   # type of marker
	set("leaves_cex", leafSize) %>% #Size
	set("leaves_col", leafColors[leafType]) %>% # up or down genes
	plot(horiz=TRUE)
	
  #legend("top",pch=19, col=leafColors[1:2],legend=levels(leafType),bty = "n",horiz =T  )
  # add legend using a second layer
  #	par(lend = 1)           # square line ends for the color legend
  # add_legend("top",pch=19, col=leafColors,legend=levels(leafType),bty = "n",horiz =T )

	
  
}


# numChar=100 maximum number of characters
# n=200  maximum number of nodes
# degree.cutoff = 0    Remove node if less connected
#from PPInfer
enrich.net2 <-  function (x, gene.set, node.id, node.name = node.id, pvalue, 
    n = 50, numChar = NULL, pvalue.cutoff = 0.05, edge.cutoff = 0.05, 
    degree.cutoff = 0, edge.width = function(x) {
        5 * x^2
    }, node.size = function(x) {
        2.5 * log10(x)
    }, group = FALSE, group.color = c("green","red" ), group.shape = c("circle", 
        "square"), legend.parameter = list("topright"), show.legend = TRUE, plotting=TRUE, 
    layoutButton = 0, ...) 
{
	library(igraph)
	set.seed(layoutButton)
    x <- data.frame(x, group)
    colnames(x)[length(colnames(x))] <- "Group"
    x <- x[as.numeric( x[, pvalue]) < pvalue.cutoff, ]
    x <- x[order(x[, pvalue]), ]
    n <- min(nrow(x), n)
    if (n == 0) {
        stop("no enriched term found...")
    }
    x <- x[1:n, ]
    index <- match(x[, node.id], names(gene.set))
    geneSets <- list()
    for (i in 1:n) {
        geneSets[[i]] <- gene.set[[index[i]]]
    }
    names(geneSets) <- x[, node.name]
    if (is.null(numChar)) {
        numChar <- max(nchar(as.character(x[, node.name])))
    }
    else {
        if (length(unique(substr(x[, node.name], 1, numChar))) < 
            nrow(x)) {
            numChar <- max(nchar(as.character(x[, node.name])))
            message("Note : numChar is too small.", "\n")
        }
    }
    x[, node.name] <- paste(substr(x[, node.name], 1, numChar), 
        ifelse(nchar(as.character(x[, node.name])) > numChar, 
            "...", ""), sep = "")
    w <- matrix(NA, nrow = n, ncol = n)

    for (i in 1:n) {
        for (j in i:n) {
            u <- unlist(geneSets[i])
            v <- unlist(geneSets[j])
            w[i, j] = length(intersect(u, v))/length(unique(c(u, 
                v)))
        }
    }
    list.edges <- stack(data.frame(w))
    list.edges <- cbind(list.edges[, 1], rep(x[, node.name], 
        n), rep(x[, node.name], each = n))
    list.edges <- list.edges[list.edges[, 2] != list.edges[,3], ]
    list.edges <- list.edges[!is.na(list.edges[, 1]), ]
    g <- graph.data.frame(list.edges[, -1], directed = FALSE)
    E(g)$width = edge.width(as.numeric(list.edges[, 1]))
    V(g)$size <- node.size(lengths(geneSets))
    g <- delete.edges(g, E(g)[as.numeric(list.edges[, 1]) < edge.cutoff])
    index.deg <- igraph::degree(g) >= degree.cutoff
    g <- delete.vertices(g, V(g)[!index.deg])
    x <- x[index.deg, ]
    index <- index[index.deg]
    if (length(V(g)) == 0) {
        stop("no categories greater than degree.cutoff...")
    }
    n <- min(nrow(x), n)
    x <- x[1:n, ]
    group.level <- sort(unique(group))
    pvalues <- log10( x[, pvalue] )
    for (i in 1:length(group.level)) {
        index <- x[, "Group"] == group.level[i]
        V(g)$shape[index] <- group.shape[i]
        group.pvalues <- pvalues[index]
        if (length(group.pvalues) > 0) {
            if (max(group.pvalues) == min(group.pvalues)) {
                V(g)$color[index] <- adjustcolor(group.color[i], 
                  alpha.f = 0.5)
            }
            else {
                V(g)$color[index] <- sapply(1 - .9* (group.pvalues - 
                  min(group.pvalues))/(max(group.pvalues) - min(group.pvalues)), 
                  function(x) {
                    adjustcolor(group.color[i], alpha.f =  .1 + x ) # change range?
                  })
            }
        }
    }
	if(plotting) { 
		plot(g, , vertex.label.dist = 1.2, ...)
		if (show.legend) {
			legend.parameter$legend <- group.level
			legend.parameter$text.col <- group.color
			legend.parameter$bty <- "n"	
			do.call(legend, legend.parameter)
		}}
    return(g)
}

enrichmentNetwork <- function(enrichedTerms, layoutButton=0, edge.cutoff = 5){
	geneLists = lapply(enrichedTerms$Genes, function(x) unlist( strsplit(as.character(x)," " )   ) )
	names(geneLists) = enrichedTerms$Pathways
	enrichedTerms$Direction = gsub(" .*","",enrichedTerms$Direction )

	g <- enrich.net2(enrichedTerms, geneLists, node.id = "Pathways", numChar = 100, 
	   pvalue = "adj.Pval",  pvalue.cutoff = 1, degree.cutoff = 0,
	   n = 200, group = enrichedTerms$Direction, vertex.label.cex = 1, 
       vertex.label.color = "black", show.legend = FALSE, 
       layoutButton = layoutButton, edge.cutoff = edge.cutoff) 

}

keggSpeciesID = read.csv(paste0(datapath,"data_go/KEGG_Species_ID.csv"))

# finds id index corresponding to entrez gene and KEGG for id conversion
idType_Entrez <- dbGetQuery(convert, paste("select distinct * from idIndex where idType = 'entrezgene_id'" ))
if(dim(idType_Entrez)[1] != 1) {cat("Warning! entrezgene ID not found!")}
idType_Entrez = as.numeric( idType_Entrez[1,1])
idType_KEGG <- dbGetQuery(convert, paste("select distinct * from idIndex where idType = 'kegg'" ))
if(dim(idType_KEGG)[1] != 1) {cat("Warning! KEGG ID not found!")}
idType_KEGG = as.numeric( idType_KEGG[1,1])

convertEnsembl2Entrez <- function (query,Species) { 
	querySet <- cleanGeneSet( unlist( strsplit( toupper(names( query)),'\t| |\n|\\,' )  ) )
	speciesID <- orgInfo$id[ which(orgInfo$ensembl_dataset == Species)]  # note uses species Identifying
	# idType 6 for entrez gene ID
	result <- dbGetQuery( convert,
						paste( " select  id,ens,species from mapping where ens IN ('", paste(querySet,collapse="', '"),
								"') AND  idType ='",idType_Entrez,"'",sep="") )	# slow
							
	if( dim(result)[1] == 0  ) return(NULL)
	result <- subset(result, species==speciesID, select = -species)

	ix = match(result$ens,names(query)  )

	tem <- query[ix];  names(tem) = result$id
	return(tem)
  
}

#Given a KEGG pathway description, found pathway ids
keggPathwayID <- function (pathwayDescription, Species, GO,selectOrg) {
	ix = grep(Species,gmtFiles)

	if (length(ix) == 0 ) {return(NULL)}
	
	# If selected species is not the default "bestMatch", use that species directly
	if(selectOrg != speciesChoice[[1]]) {  
		ix = grep(findSpeciesById(selectOrg)[1,1], gmtFiles )
		if (length(ix) == 0 ) {return(NULL )}
		totalGenes <- orgInfo[which(orgInfo$id == as.numeric(selectOrg)),7]
	}
	pathway <- dbConnect(sqlite,gmtFiles[ix],flags=SQLITE_RO)
	
	# change Parkinson's disease to Parkinson\'s disease    otherwise SQL 
	pathwayDescription <- gsub("\'","\'\'",pathwayDescription)
							
	pathwayInfo <- dbGetQuery( pathway, paste( " select * from pathwayInfo where description =  '", 
							pathwayDescription,   "' AND name LIKE '",GO,"%'",sep="") )
	dbDisconnect(pathway);
	if(dim(pathwayInfo)[1] != 1 ) {return(NULL) }
	tem = gsub(".*:","",pathwayInfo[1,2])  
	return( gsub("_.*","",tem) )
}

# not working, for updating GO cateory choices
gmtCategory <- function (converted,  selectOrg) {

	idNotRecognized = as.data.frame("ID not recognized!")
	if(is.null(converted) ) return(idNotRecognized) # no ID 
	querySet <- converted$IDs
	if(length(querySet) == 0) return(idNotRecognized )
	ix = grep(converted$species[1,1],gmtFiles)
	if (length(ix) == 0 ) {return(idNotRecognized )}
	
	# If selected species is not the default "bestMatch", use that species directly
	if(selectOrg != speciesChoice[[1]]) {  
		ix = grep(findSpeciesById(selectOrg)[1,1], gmtFiles )
		if (length(ix) == 0 ) {return(idNotRecognized )}
	}
	pathway <- dbConnect(sqlite,gmtFiles[ix],flags=SQLITE_RO)
	#cat(paste("selectOrg:",selectOrg) )
	# Generate a list of geneset categories such as "GOBP", "KEGG" from file
	geneSetCategory <-  dbGetQuery(pathway, "select distinct * from categories " ) 
	geneSetCategory  <- sort( geneSetCategory[,1] )
	categoryChoices <- setNames(as.list( geneSetCategory ), geneSetCategory )
	categoryChoices <- append( setNames( "All","All available gene sets"), categoryChoices  )
	
	# move one element to the 2nd place
	move1 <- function(i) c(categoryChoices[1],categoryChoices[i],categoryChoices[-c(1,i)])
	i = which( names(categoryChoices)  == "KEGG"); categoryChoices= move1(i);	
	i = which( names(categoryChoices)  == "GOMF"); categoryChoices= move1(i);	
	i = which( names(categoryChoices)  == "GOCC"); categoryChoices= move1(i);	
	i = which( names(categoryChoices)  == "GOBP"); categoryChoices= move1(i);
	#change GOBP to the full description for display
	names(categoryChoices)[ match("GOBP",categoryChoices)  ] <- "GO Biological Process"
	names(categoryChoices)[ match("GOCC",categoryChoices)  ] <- "GO Cellular Component"
	names(categoryChoices)[ match("GOMF",categoryChoices)  ] <- "GO Molecular Function"
	
	dbDisconnect(pathway)
	return(categoryChoices )
} 
 