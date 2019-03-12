library(shiny)
library(RSQLite)
library(gplots)
library(PGSEA)
Min_overlap <- 2
minSetSize = 3; 
mappingCoverage = 0.60 # 60% percent genes has to be mapped for confident mapping
mappingEdge = 0.5  # Top species has 50% more genes mapped
PvalGeneInfo = 0.05; minGenes = 10 # min number of genes for ploting
sqlite  <- dbDriver("SQLite")
convert <- dbConnect(sqlite,"../go/convertIDs.db")

hclust2 <- function(x, method="average", ...)
  hclust(x, method=method, ...)
dist2 <- function(x, ...)
  as.dist(1-cor(t(x), method="pearson"))
 detectGroups <- function (x){  # x are col names
 	groups=as.factor( gsub("[0-9]*$","",x) )
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

	geneChange <- function(x){
	  n = length(x)
	  if( n<4) return( max(x)-min(x)  ) else 
	  return(sort(x)[n-1] - sort(x)[2]   )
	}
myPGSEA  <- function (exprs, cl, range = c(25, 500), ref = NULL, center = TRUE, 
    p.value = 0.005, weighted = TRUE, nPermutation=100, enforceRange = TRUE, ...) 
{
    if (is(exprs, "ExpressionSet")) 
        exprs <- exprs(exprs)
    if (!is.list(cl)) 
        stop("cl need to be a list")
    if (!is.null(ref)) {
        if (!is.numeric(ref)) 
            stop("column index's required")
    }
    if (!is.null(ref)) {
        if (options()$verbose) 
            cat("Creating ratios...", "\n")
        ref_mean <- apply(exprs[, ref], 1, mean, na.rm = TRUE)
        exprs <- sweep(exprs, 1, ref_mean, "-")
    }
    if (center) 
        exprs <- scale(exprs, scale = FALSE)         # column centering is done
    results <- matrix(NA, length(cl), ncol(exprs))
    rownames(results) <- names(cl)
    colnames(results) <- colnames(exprs)
    mode(results) <- "numeric"
	Setsize = c(rep(0,length(cl)))     # gene set size vector
	mean2 = c(rep(0,length(cl)))     # mean of the range of means 
	meanSD = c(rep(0,length(cl)))     # SD of the range of means	
    if (is.logical(p.value)) 
        { p.results <- results; mean.results <- results;}
    for (i in 1:length(cl)) {              # for each gene list
		#cat("\nProcessing gene set",i);
        if (class(cl[[i]]) == "smc") {
            clids <- cl[[i]]@ids
        }
        else if (class(cl[[i]]) %in% c("GeneColorSet", "GeneSet")) {
            clids <- cl[[i]]@geneIds
        }
        else {
            clids <- cl[[i]]
        }
        if (options()$verbose) 
            cat("Testing region ", i, "\n")
        ix <- match(clids, rownames(exprs))
        ix <- unique(ix[!is.na(ix)])
        present <- sum(!is.na(ix))
		Setsize[i] <- present 
        if (present < range[1]) {
            if (options()$verbose) 
                cat("Skipping region ", i, " because too small-", 
                  present, ",\n")
            next
        }
        if (present > range[2]) {
            if (options()$verbose) 
                cat("Skipping region ", i, " because too large-", 
                  present, "\n")
            next
        }
        texprs <- exprs[ix, ]           # expression matrix for genes in gene set
        if (any(is.na(texprs))) 
            cat("Warning - 'NA' values within expression data, enrichment scores are estimates only.\n")
        if (!is.matrix(texprs)) 
            texprs <- as.matrix(texprs)
                            
        stat <- try(apply(texprs, 2, t.test, ...))
		means <- try(apply(texprs, 2, mean,trim=0.1))   # trim mean
		ps <- unlist(lapply(stat, function(x) x$p.value))
        stat <- unlist(lapply(stat, function(x) x$statistic))
        p.results[i, ] <- ps
		mean.results[i,] <- means
        results[i, ] <- as.numeric(stat)
		
		# permutation of gene sets of the same size
		if(nPermutation > 2 )  { # no permutation if <=2
			meansRanges = c(0, rep(nPermutation))
			for( k in 1:nPermutation ) {
				ix <- sample.int( dim(exprs)[1], length(ix) )
				texprs <- exprs[ix, ] 
				means <- try(apply(texprs, 2, mean,trim=0.1))   # trim mean
				meansRanges[k] = dynamicRange(means)
			}
			mean2[i] = mean(meansRanges)
			meanSD[i]= sd(meansRanges,na.rm=TRUE)   # NA are removed before calculating standard deviation
		}
    }
    return(list(results = results, p.results = p.results, means = mean.results, size=Setsize, mean2=mean2, meanSD=meanSD))
    
}

dynamicRange <- function( x ) {
y = sort(x)
   if(length(x)>=4)  k =2 else k =1;
   return( y[length(x)-k+1] - y[k]) 
}

# Create a list of GMT files in /gmt sub folder
gmtFiles = list.files(path = "../go/pathwayDB",pattern=".*\\.db")
gmtFiles = paste("../go/pathwayDB/",gmtFiles,sep="")

geneInfoFiles = list.files(path = "../go/geneInfo",pattern=".*GeneInfo\\.csv")
geneInfoFiles = paste("../go/geneInfo/",geneInfoFiles,sep="")

motifFiles = list.files(path = "../go/motif",pattern=".*\\.db")
motifFiles = paste("../go/motif/",motifFiles,sep="")


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
convertID <- function (query,selectOrg, selectGO) { 
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
  
  if(0){
  # generate a list of gene set categories
  ix = grep(findSpeciesById(result$species[1])[1,1],gmtFiles)
   
  if (length(ix) == 0 ) {categoryChoices = NULL}
  
  # If selected species is not the default "bestMatch", use that species directly
  if(selectOrg != speciesChoice[[1]]) {  
    ix = grep(findSpeciesById(selectOrg)[1,1], gmtFiles )
	if (length(ix) == 0 ) {categoryChoices = NULL}
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
  
  
  
  

  
  
  
  
  
  
  dbDisconnect(pathway)
  
  } #if (0)
  return(list(originalIDs = querySet,IDs=unique( result[,2]), 
              species = findSpeciesById(result$species[1]), 
              #idType = findIDtypeById(result$idType[1] ),
              speciesMatched = speciesMatched,
			  conversionTable = conversionTable
			  ) )
} 

geneInfo <- function (myData, converted,selectOrg){
  # query = scan("query_temp.txt",what=""); selectOrg ="BestMatch";
   # query = scan("zebrafish_test.gmt", what="" ); selectOrg ="BestMatch";
  # query = scan("Celegans_test.gmt", what="" ); selectOrg ="BestMatch";
  # query = scan("test_query_mouse_symbol.txt", what="" ); selectOrg ="BestMatch";
  #  query = scan("soy_test.txt", what="" );selectOrg ="BestMatch";
  # querySet <- cleanGeneSet( unlist( strsplit( toupper(query),'\t| |\n|\\,')))
  # converted = convertID( querySet,selectOrg)
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
  { x = read.csv(as.character(geneInfoFiles[ix]) ); x[,1]= toupper(x[,1]) } else # read in the chosen file 
  { return(as.data.frame("Multiple geneInfo file found!") )   }
  Set = match(x$ensembl_gene_id, querySet)
  Set[which(is.na(Set))]="Genome"
  Set[which(Set!="Genome")] ="List"
  # x = cbind(x,Set) } # just for debuging
  return( cbind(x,Set) )}
 }

# Main function. Find a query set of genes enriched with functional category
gmtCategory <- function (converted, convertedData, selectOrg) {
  idNotRecognized = as.data.frame("ID not recognized!")
  if(is.null(converted) ) return(idNotRecognized) # no ID 
  querySet <- rownames(convertedData)
  if(length(querySet) == 0) return(idNotRecognized )
  ix = grep(converted$species[1,1],gmtFiles)
  if (length(ix) == 0 ) {return(idNotRecognized )}
  
  # If selected species is not the default "bestMatch", use that species directly
  if(selectOrg != speciesChoice[[1]]) {  
    ix = grep(findSpeciesById(selectOrg)[1,1], gmtFiles )
	if (length(ix) == 0 ) {return(idNotRecognized )}
  }
  pathway <- dbConnect(sqlite,gmtFiles[ix])
 cat(paste("selectOrg:",selectOrg) )
  # Generate a list of geneset categories such as "GOBP", "KEGG" from file
  geneSetCategory <-  dbGetQuery(pathway, "select distinct * from categories " ) 
  geneSetCategory  <- geneSetCategory[,1]
  categoryChoices <- setNames(as.list( geneSetCategory ), geneSetCategory )
  categoryChoices <- append( setNames( "All","All available gene sets"), categoryChoices  )
  #change GOBO to the full description for display
  names(categoryChoices)[ match("GOBP",categoryChoices)  ] <- "GO Biological Process"
  names(categoryChoices)[ match("GOCC",categoryChoices)  ] <- "GO Cellular Component"
  names(categoryChoices)[ match("GOMF",categoryChoices)  ] <- "GO Molecular Function"
  
  dbDisconnect(pathway)
 return(categoryChoices )
} 
 
 
# Main function. Find a query set of genes enriched with functional category
readGeneSets <- function (converted, convertedData, GO,selectOrg, myrange) {
  idNotRecognized = as.data.frame("ID not recognized!")
  if(is.null(converted) ) return(idNotRecognized) # no ID 
  querySet <- rownames(convertedData)
  if(length(querySet) == 0) return(idNotRecognized )
  ix = grep(converted$species[1,1],gmtFiles)
  if (length(ix) == 0 ) {return(idNotRecognized )}
  
  # If selected species is not the default "bestMatch", use that species directly
  if(selectOrg != speciesChoice[[1]]) {  
    ix = grep(findSpeciesById(selectOrg)[1,1], gmtFiles )
	if (length(ix) == 0 ) {return(idNotRecognized )}
  }
  pathway <- dbConnect(sqlite,gmtFiles[ix])
 
  if(is.null(GO) ) GO <- "GOBP"   # initial value not properly set; enforcing  

 # get Gene sets
  querySet = rownames(convertedData)
  sqlQuery = paste( " select distinct gene,pathwayID from pathway where gene IN ('", paste(querySet,collapse="', '"),"')" ,sep="")
 # cat(paste0("\n\nhere:",GO,"There"))

  if( GO != "All") sqlQuery = paste0(sqlQuery, " AND category ='",GO,"'")
  result <- dbGetQuery( pathway, sqlQuery  )
  if( dim(result)[1] ==0) {return(list( x=as.data.frame("No matching species or gene ID file!" )) )}
  # list pathways and frequency of genes
  pathwayIDs = aggregate( result$pathwayID, by   = list(unique.values = result$pathwayID), FUN = length)
  pathwayIDs = pathwayIDs[which(pathwayIDs[,2]>= myrange[1] ),]
  pathwayIDs = pathwayIDs[which( pathwayIDs[,2] <= myrange[2] ),]
  if(dim(pathwayIDs)[1] ==0 ) geneSets = NULL;
  
   # convert pathways into lists like those generated by readGMT
   geneSets = lapply(pathwayIDs[,1], function(x)  result[which(result$pathwayID == x ),1]     )
   pathwayInfo <- dbGetQuery( pathway, paste( " select distinct id,Description from pathwayInfo where id IN ('", 
						paste(pathwayIDs[,1],collapse="', '"),   "') ",sep="") )
   ix = match( pathwayIDs[,1], pathwayInfo[,1])
   names(geneSets) <- pathwayInfo[ix,2]  

  dbDisconnect(pathway)
 return( geneSets )
} 

 
 
 PGSEApathway <- function (converted,convertedData, selectOrg,GO,gmt, myrange,Pval_pathway,top){
  	subtype = detectGroups(colnames(convertedData))
	Pvalue = 0.01  # cut off to report in PGSEA. Otherwise NA
	#Pval_pathway = 0.2   # cut off for P value of ANOVA test  to writ to file 
	# top = 30   # number of pathways to show
	if(length(gmt) ==0 ) return( list(pg3 = NULL, best = best ) )
    # centering by mean
	#pg = myPGSEA (convertedData - rowMeans(convertedData),
	 #            cl=gmt,range=myrange,p.value=TRUE, weighted=FALSE,nPermutation=100)
	pg = PGSEA (convertedData - rowMeans(convertedData),cl=gmt,range=myrange,p.value=TRUE, weighted=FALSE)
	
	pg2 = pg$results;
	pg2 = pg2[rowSums(is.na(pg2))<ncol(pg2) ,]  # remove se/wrts with all missing(non-signficant)
	if (dim(pg2)[1] < 2 ) return()
	best = max(abs(pg2))
	
	if(length(subtype) < 4 || length(unique(subtype)) <2 ||length(unique(subtype)) == dim(convertedData)[2] ) { 
	 pg2 = pg2[order(-apply(pg2,1,sd)     )   ,]
	 return( list(pg3 = pg2[1:top,], best = best ) )
	} 
    
	cat("\nComputing P values using ANOVA\n");
	pathPvalue <- function ( k){
	 return( summary(aov(pg2[k,]~subtype) )[[1]][["Pr(>F)"]][1] )
	}
	Pvalues = sapply(1:dim(pg2)[1], pathPvalue)
	Pvalues = p.adjust(Pvalues, "fdr")
	
	#if(min(Pvalues) > Pval_pathway ) return( list(pg3 = NULL, best = best ) )  else {  
    if(sort(Pvalues)[2] > Pval_pathway ) return( list(pg3 = NULL, best = best ) )  else {  

	NsigT = rowSums(pg$p.results<Pvalue)
	
	result=cbind( as.matrix(Pvalues),NsigT,pg2); 
	result = result[ order(result[,1])   ,]	
	
	result = result[which(result[,1] < Pval_pathway)    ,]
	#result = result[which(result[,2] >2)    ,]


	
	pg2 = result[,-2]

	# when there is only 1 left in the matrix pg2 becomes a vector
	if(sum( Pvalues<Pval_pathway) == 1) { pg3 = t( as.matrix(pg2));pg3 = rbind(pg3,pg3);} else
	{ if(dim(pg2)[1] > top ) {  pg3 = pg2[1:top,]; } else {  pg3 = pg2;  } }

	rownames(pg3) = sapply(rownames(pg3) , extract1)
	a=sprintf("%-1.0e",pg3[,1])
	rownames(pg3) = paste(a,rownames(pg3),sep=" ")
	pg3 =pg3[,-1]
 
    return( list(pg3 = pg3, best = best ) )
    }
 }
 
 if(0 ){
 x = read.csv("expression1.csv")
 x = read.csv("expression1_no_duplicate.csv")
 x = read.csv("mouse1.csv")
 x = read.csv("GSE40261.csv")
 x = x[order(x[,1]),]
    x = x[!duplicated(x[,1]),]
    rownames(x)= x[,1]
    x = x[,-1]
	
	tem = apply(x,1,max)
	x = x[which(tem> 1),] 
	
	x = log(x+abs( 10),2)
	tem = apply(x,1,sd)
	x = x[order(-tem),]
	
	selectOrg = "BestMatch"; GO="GOBP"; 
	myrange = c(15,1000)

	converted = convertID(rownames(x),selectOrg)
	
	head(converted$conversionTable)
	mapping = converted$conversionTable
	
	   rownames(x) = toupper(rownames(x))
	   x1 = merge(mapping[,1:2],x,  by.y = 'row.names', by.x = 'User_input')
	  	tem = apply(x1[,3:(dim(x1)[2]-2)],1,sd)
	   x1 = x1[order(x1[,2],-tem),]
	   x1 = x1[!duplicated(x1[,2]) ,]
	   rownames(x1) = x1[,2]
	   x1 = as.matrix(x1[,c(-1,-2)])
	
	convertedData = x1
    gmt = readGeneSets(converted, convertedData, GO,selectOrg, myrange)

subtype = detectGroups(colnames(convertedData))
	Pvalue = 1  # cut off to report in PGSEA. Otherwise NA
	Pval_pathway = 0.05   # cut off for P value of ANOVA test  to writ to file 
	top = 30   # number of pathways to show
    myrange = c(10,2000)
	
	pg = myPGSEA (x2,cl=gmt,range=myrange,p.value=TRUE, weighted=FALSE,nPermutation=1)
	
	tem = PGSEApathway (converted,convertedData, selectOrg,GO,gmt, myrange)
 }
 
 
 
 
 
 
shinyServer(

  function(input, output,session) {
   options(shiny.maxRequestSize = 50*1024^2) # 50MB file max 
   
   observe({  updateSelectInput(session, "selectOrg", choices = speciesChoice )      })
   
    # read data file
  readData <- reactive ({ 
    inFile <- input$file1
	inFile <- inFile$datapath
    if (is.null(input$file1) && input$goButton == 0)   return(NULL)
    if (is.null(input$file1) && input$goButton > 0 )   inFile = "expression1_no_duplicate.csv"
    x = read.csv(inFile)
	x[,1] = toupper(x[,1])
    x = x[order(x[,1]),]
    x = x[!duplicated(x[,1]),]
    rownames(x)= x[,1]
    x = x[,-1]
	
	tem = apply(x,1,max)
	x = x[which(tem> input$lowFilter),] 
	
	if (input$transform) x = log(x+abs( input$logStart),2)
	tem = apply(x,1,sd)
	x = x[order(-tem),]
	
    as.matrix(x)	
  })	
  	
  
	# this defines an reactive object that can be accessed from other rendering functions
	converted <- reactive({  
	  if (is.null(input$file1) && input$goButton == 0)    return()
       tem = input$selectOrg
      isolate( {  convertID(rownames(readData() ),input$selectOrg, input$selectGO );}) 

	} )

	convertedData <- reactive({  
	  if (is.null(input$file1) && input$goButton == 0)    return()
		tem = input$selectOrg
      isolate( {  
	   mapping <- converted()$conversionTable
	   x =readData()
	   rownames(x) = toupper(rownames(x))
	   x1 = merge(mapping[,1:2],x,  by.y = 'row.names', by.x = 'User_input')
	  	tem = apply(x1[,3:(dim(x1)[2]-2)],1,sd)
	   x1 = x1[order(x1[,2],-tem),]
	   x1 = x1[!duplicated(x1[,2]) ,]
	   rownames(x1) = x1[,2]
	   x1 = as.matrix(x1[,c(-1,-2)])
       return(x1)
	  }) 

	} )	
	
   GeneSets <- reactive({  
	  if (is.null(input$file1) && input$goButton == 0)    return()
	  tem = input$goButton2 
	  tem = input$selectOrg
      isolate( {
	  
	   readGeneSets( converted(), convertedData(), input$selectGO,input$selectOrg,c(input$minSetSize, input$maxSetSize)  )  }) 

	} )
     gmtCategory1 <- reactive({  
	  if (is.null(input$file1) && input$goButton == 0)    return()
	  tem = input$goButton2 
	  tem = input$selectOrg
      isolate( {
	  
	   gmtCategory(converted(), convertedData(), input$selectOrg) }) 

	} )
  output$contents <- renderTable({
    if (is.null(input$file1) && input$goButton == 0)   return(NULL)
	tem = input$selectOrg
    x <- readData()
    x[1:20,]
  },include.rownames=TRUE)
  
  output$species <-renderTable({   
      if (is.null(input$file1) && input$goButton == 0)    return()
      isolate( {  #tem <- convertID(input$input_text,input$selectOrg );
	  	  withProgress(message="Converting gene IDs", {
                  tem <- converted()
			incProgress(1, detail = paste("Done"))	  })
		  
				  if( is.null(tem)) {as.data.frame("ID not recognized.")} else {
	              tem$speciesMatched }

      }) # avoid showing things initially         
    }, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)
	
  
  output$EDA <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    x <- readData()
	 par(mfrow=c(2,1))
	par(mar=c(8,2,2,2))
	boxplot(x, las = 3, ylab="Expression values")
	myColors = rainbow(dim(x)[2])
	plot(density(x[,1]),col = myColors[1],xlab="Expresson values", ylab="Density", main= "Distribution of expression values")
	for( i in 2:dim(x)[2] )
	   lines(density(x[,i]),col=myColors[i] )
    legend("topright", colnames(x), lty=rep(1,dim(x)[2]), col=myColors )	
    
  }, height = 1000, width = 500)
  
    output$PCA <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    x <- readData()

	 pca.object <- prcomp(t(x))
	 
     plot( pca.object$x[,1], pca.object$x[,2], pch = 1,cex = 2,col = detectGroups(colnames(x)),
	     xlim=c(min(pca.object$x[,1]),max(pca.object$x[,1])*1.3   ),
	  xlab = "First principal component", ylab="Second Principal Component")
	   text( pca.object$x[,1], pca.object$x[,2],  pos=4, labels =colnames(x), offset=.5, cex=1)
	  
  }, height = 600, width = 600)
  
  output$heatmap <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    x <- readData()   # x = read.csv("expression1.csv")

	n=input$nGenes
	if(n>4000) n = 4000 # max
	# this will cutoff very large values, which could skew the color 
	x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
	x[x>cutoff] <- cutoff
	cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
	x[x< cutoff] <- cutoff
	
    groups = detectGroups(colnames(x) )
	groups.colors = rainbow(length(unique(groups) ) )

	#http://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
	lmat = rbind(c(5,4),c(0,1),c(3,2))
	lwid = c(1.5,6)
	lhei = c(1,.2,8)

	if( n>100) 
	heatmap.2(x, distfun = dist2,hclustfun=hclust2,
	 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
	,key=T, symkey=F
	,ColSideColors=groups.colors[ groups]
	,labRow=""
	,margins=c(10,0)
	,srtCol=45
	,lmat = lmat, lwid = lwid, lhei = lhei
	)

	if( n<=100) 
	heatmap.2(x, distfun = dist2,hclustfun=hclust2,
	 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
	,key=T, symkey=F,
	#,labRow=labRow
		,ColSideColors=groups.colors[ groups]
	,margins=c(10,12)
	,cexRow=1.5
	,srtCol=45
	,lmat = lmat, lwid = lwid, lhei = lhei
	)
	
  }, height = 800, width = 600)  
    
  heatmapData <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    x <- readData()
	geneSD = apply(x,1,sd)
	x = x[order(-geneSD),]
	n=input$nGenes
	if(n>4000) n = 4000 # max
	# this will cutoff very large values, which could skew the color 
	x1 = x[1:n,] 
	x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
	x[x>cutoff] <- cutoff
	cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
	x[x< cutoff] <- cutoff
	
	groups = detectGroups(colnames(x) )
	groups.colors = rainbow(length(unique(groups) ) )

	pdf(file=NULL,width =700, height =700)
	hy <- heatmap.2(x, distfun = dist2,hclustfun=hclust2,#labRow="",labCol="",
	density.info="none", trace="none", scale="none")
    dev.off()
	
	return(x1[ rev( hy$rowInd),hy$colInd])
	
  })  
  
  output$downloadData <- downloadHandler(
     filename = function() {"heatmap.csv"},
		content = function(file) {
      write.csv(heatmapData(), file)
	    }
  )

  
  
  
  
  
  output$debug <- renderTable({
      if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	  tem = input$selectOrg
	x <- convertedData()
	#tem = GeneSets()
	
	return( as.data.frame (x[1:20,] ))
	},include.rownames=TRUE)
	
	    # this updates geneset categories based on species and file
    output$selectGO1 <- renderUI({
	  tem = input$selectOrg
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectGO", label = NULL, # h6("Funtional Category"), 
                  choices = list("All available gene sets" = "All", "GO Biological Process" = "GOBP","GO Molecular Function" = "GOMF","GO Cellular Component" = "GOCC",
                                "KEGG metabolic pathways" = "KEGG"), selected = "GOBP")  }	 else { 
								
	  selectInput("selectGO", label=NULL,choices=gmtCategory(converted(), convertedData(), input$selectOrg)
	     ,selected = "GOBP" )   } 
	})
	
  
  output$PGSEAplot <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg
	tem = input$goButton2 # just to make it rerun the analysis
	isolate({ 
	withProgress(message="Running pathway analysis", {
	myrange = c(input$minSetSize, input$maxSetSize)
	subtype = detectGroups(colnames(convertedData() ))
	gmt = GeneSets()
	incProgress(2/8)
    if(length( GeneSets() )  == 0)  { plot.new(); text(0,1, "No gene sets!")} else {
	result = PGSEApathway(converted(),convertedData(), input$selectOrg,input$selectGO,
	             GeneSets(),  myrange, input$pathwayPvalCutoff, input$nPathwayShow 	)
					 
	if( is.null(result$pg3) ) { plot.new(); text(0.5,1, "No significant pathway found!")} else 
	smcPlot(result$pg3,factor(subtype),scale = c(-result$best, result$best), show.grid = T, margins = c(3,1, 13, 23), col = .rwb,cex.lab=0.5)
    }
	
	})
	})
    }, height = 800, width = 500)

})
