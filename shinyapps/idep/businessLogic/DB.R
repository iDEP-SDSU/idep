library('R6')
library(RSQLite,verbose=FALSE)	# for database connection

source('server.config')

DB.Manager <- R6Class("DB.Manager")

DB.Manager$set("public","sqlite", NULL)
DB.Manager$set("public","dbConvert", NULL)
DB.Manager$set("public","gmtFiles", NULL)	
DB.Manager$set("public","OrgInfo", NULL)	# this field does not allow direct access, use 'QueryOrgInfo' to get org info
DB.Manager$set("public","Quotes", NULL)

## Init function
DB.Manager$set("public","initialize",
	function(){
		self$sqlite <- dbDriver("SQLite")
		
		self$dbConvert <- dbConnect( self$sqlite, paste0(CONFIG_DATA_DATAPATH, "convertIDs.db"), flags=SQLITE_RO) 
		
		self$gmtFiles = list.files(path = paste0(CONFIG_DATA_DATAPATH,"pathwayDB"), pattern=".*\\.db")
		self$gmtFiles = paste(CONFIG_DATA_DATAPATH, "pathwayDB/", self$gmtFiles,sep="")
		dbQueryResult <- dbGetQuery(self$dbConvert, paste("select distinct * from orgInfo " ))
		self$OrgInfo <- dbQueryResult[order(dbQueryResult$name),]

		TmpQuotes <- dbGetQuery(self$dbConvert, " select * from quotes")
		self$Quotes <- paste0("\"",TmpQuotes$quotes,"\"", " -- ",TmpQuotes$author,".       ")

	}
)

DB.Manager$set("public", "GetSpeciesChoice",
    function(){
        orgInfo <- self$OrgInfo

		speciesChoice <- setNames(as.list( orgInfo$id ), orgInfo$name2 )
		# add a defult element to list    # new element name       value
		speciesChoice <- append( setNames( "NEW","**NEW SPECIES**"), speciesChoice  )
		speciesChoice <- append( setNames( "BestMatch","Best matching species"), speciesChoice  )

		# move one element to the 2nd place
		move2 <- function(i) c(speciesChoice[1:2],speciesChoice[i],speciesChoice[-c(1,2,i)])
		i= which( names(speciesChoice) == "Glycine max"); speciesChoice <- move2(i)
		i= which( names(speciesChoice) =="Zea mays"); speciesChoice <- move2(i)
		i= which(names(speciesChoice) =="Arabidopsis thaliana"); speciesChoice <- move2(i)
		i= which(names(speciesChoice) == "Saccharomyces cerevisiae"); speciesChoice <- move2(i)
		i= which(names(speciesChoice)  == "Caenorhabditis elegans"); speciesChoice <- move2(i)
		i= which(names(speciesChoice) =="Zebrafish" ); speciesChoice <- move2(i)
		i= which(names(speciesChoice) == "Cow" ); speciesChoice <- move2(i)
		i= which(names(speciesChoice) == "Rat" ); speciesChoice <- move2(i)
		i= which(names(speciesChoice) == "Mouse"); speciesChoice <- move2(i)
		i= which(names(speciesChoice) == "Human"); speciesChoice <- move2(i)
        return(speciesChoice)
    }
)

# find species information use id
DB.Manager$set("public", "findSpeciesById",
	function(speciesID){ 
		orgInfo <- self$OrgInfo
  		return( orgInfo[which(orgInfo$id == speciesID),]  )
	}
)


# find species name use id
DB.Manager$set("public", "findSpeciesNameById",
	function(speciesID){ 
		orgInfo <- self$OrgInfo
  		return( orgInfo[which(orgInfo$id == speciesID),3]  )
	}
)



# get id, ens, species from convertID mapping table
DB.Manager$set("public", "QuerySpeciesInfoFromConvertDBMapping",
	function(querySet){
		querySTMT <- paste( "select distinct id,ens,species from mapping where id IN ('", paste(querySet,collapse="', '"),"')",sep="")
		result <- dbGetQuery(self$dbConvert, querySTMT)
		return(result)
	}
)


#Find a query set of genes enriched with functional category
DB.Manager$set("public", "QueryGeneSetsFromPathway",
	function (converted, convertedData, GO, selectOrg, myrange) {
		idNotRecognized = as.data.frame("ID not recognized!")
		if(is.null(converted) ) return(idNotRecognized) # no ID 
		querySet <- rownames(convertedData)
		if(length(querySet) == 0) return(idNotRecognized )
		ix = grep(converted$species[1,1],self$gmtFiles)
		if (length(ix) == 0 ) {return(idNotRecognized )}
		
		# If selected species is not the default "bestMatch", use that species directly
		if(selectOrg != "BestMatch") {  
			ix = grep(self$findSpeciesById(selectOrg)[1,1], self$gmtFiles )
			if (length(ix) == 0 ) {return(idNotRecognized )}
		}
		pathway <- dbConnect(self$sqlite, self$gmtFiles[ix], flags=SQLITE_RO)
		
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
		#geneSets <- geneSets[ -which(duplicated(names(geneSets) ))] # remove geneSets with the same name
		dbDisconnect(pathway)
		return( geneSets )
	}
) 

DB.Manager$set("public", "cleanGeneSet",
    function(x){
        # remove duplicate; upper case; remove special characters
        x <- unique( toupper( gsub("\n| ","",x) ) )
        x <- x[which( nchar(x)>1) ]  # genes should have at least two characters
        return(x)
    }
)

DB.Manager$set("public", "FindOverlap",
    function(converted, gInfo, GO, selectOrg, minFDR, reduced = FALSE){
        maxTerms =15 # max number of enriched terms
        idNotRecognized = as.data.frame("ID not recognized!")
        
        if(is.null(converted) ) return(idNotRecognized) # no ID 
        if(is.null(GO)){
          return(NULL)
        }
        # only coding
        gInfo <- gInfo[which( gInfo$gene_biotype == "protein_coding"),]  
        querySet <- intersect( converted$IDs, gInfo[,1]);
        
        if(length(querySet) == 0) return(idNotRecognized )
        
        ix = grep(converted$species[1,1],self$gmtFiles )
        totalGenes <- converted$species[1,7]
        
        if (length(ix) == 0 ) {return(idNotRecognized )}
        
        # If selected species is not the default "bestMatch", use that species directly
        if(selectOrg != self$GetSpeciesChoice()[[1]]) {  
            ix = grep(self$findSpeciesById(selectOrg)[1,1], self$gmtFiles  )
            if (length(ix) == 0 ) {return(idNotRecognized )}
            totalGenes <- self$orgInfo[which(orgInfo$id == as.numeric(selectOrg)),7]
        }
        pathway <- dbConnect(self$sqlite,self$gmtFiles[ix],flags=SQLITE_RO)
        
            
        sqlQuery = paste( " select distinct gene,pathwayID from pathway where gene IN ('", paste(querySet,collapse="', '"),"')" ,sep="")
        #cat(paste0("HH",GO,"HH") )
        
        if( GO != "All") sqlQuery = paste0(sqlQuery, " AND category ='",GO,"'")
        result <- dbGetQuery( pathway, sqlQuery  )
        if( dim(result)[1] ==0) {return(as.data.frame("No matching species or gene ID file!" )) }

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
        
        if(dim(x)[1] > maxTerms ) x = x[1:maxTerms,]	
        
        if(min(x$FDR) > minFDR) x=as.data.frame("No significant enrichment found!") else {
            x <- x[which(x$FDR < minFDR),] 

            x= cbind(x,sapply( x$pathwayID, sharedGenesPrefered ) )
            colnames(x)[7]= "Genes"
            x <- subset(x,select = c(FDR,overlap,n,description,Genes) )
            colnames(x) = c("Corrected P value (FDR)", "Genes in list", "Total genes in category","Functional Category","Genes"  )
            
            # remove redudant gene sets
            if(reduced != FALSE ){  # reduced=FALSE no filtering,  reduced = 0.9 filter sets overlap with 90%
                n=  nrow(x)
                tem=rep(TRUE,n )
                geneLists = lapply(x$Genes, function(y) unlist( strsplit(as.character(y)," " )   ) )
                for( i in 2:n)
                    for( j in 1:(i-1) ) { 
                    if(tem[j]) { # skip if this one is already removed
                        commonGenes = length(intersect(geneLists[i] ,geneLists[j] ) )
                        if( commonGenes/ length(geneLists[j] ) > reduced )
                            tem[i] = FALSE	
                    }			
                    }								
                x <- x[which(tem),]		
            }
            

        }
                
        dbDisconnect(pathway)
        return(x)
    }
)

DB.Manager$set("public", "FindOverlapGMT",
# Given a gene set, finds significant overlaps with a gene set database  object 
    function( query, geneSet, minFDR=.2 ,minSize=2,maxSize=10000 ){
        total_elements = 30000
        Min_overlap <- 1
        maxTerms =10 # max number of enriched terms
        noSig <- as.data.frame("No significant enrichment found!")
        query <- self$cleanGeneSet(query)   # convert to upper case, unique()

        if(length(query) <=2) return(noSig)
        if(length(geneSet) <1) return(noSig)
        geneSet <- geneSet[which(sapply(geneSet,length) > minSize)]  # gene sets smaller than 1 is ignored!!!
        geneSet <- geneSet[which(sapply(geneSet,length) < maxSize)]  # gene sets smaller than 1 is ignored!!!
        result <- unlist( lapply(geneSet, function(x) length( intersect (query, x) ) ) )
        result <- cbind(unlist( lapply(geneSet, length) ), result )
        result <- result[ which(result[,2]>Min_overlap), ,drop=F]
        if(dim(result)[1] == 0) return( noSig)
        xx <- result[,2]
        mm <- length(query)
        nn <- total_elements - mm
        kk <- result[,1]
        Pval_enrich=phyper(xx-1,mm,nn,kk, lower.tail=FALSE );
        FDR <- p.adjust(Pval_enrich,method="fdr",n=length(geneSet) )
        result <- as.data.frame(cbind(FDR,result))
        result <- result[,c(1,3,2)]
        result$pathway = rownames(result)
        result$Genes = ""  # place holder just 
        colnames(result)= c("Corrected P value (FDR)", "Genes in list", "Total genes in category","Functional Category","Genes"  )
        result <- result[ which( result[,1] < minFDR),,drop=F]
        if( dim( result)[1] == 0) return(noSig) 
        if(min(FDR) > minFDR) return(noSig) 
        result <- result[order(result[,1] ),]
        if(dim(result)[1] > maxTerms ) result <- result[1:maxTerms,]

        return( result)
    }

)


# read gene set files in the GMT format, does NO cleaning. Assumes the GMT files are created with cleanGeneSet()
# See http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Gene_Set_Database_Formats
readGMT <- function (fileName){ 
  x <- scan(fileName, what="", sep="\n")
  x <- strsplit(x, "\t")
  # Extract the first vector element and set it as the list element name
  names(x) <- sapply(x, `[[`, 1)
  x <- lapply(x, `[`, -c(1,2)) # 2nd element is comment, ignored
  x = x[which(sapply(x,length) > 1)]  # gene sets smaller than 1 is ignored!!!
  return(x)
}



# Read gene sets GMT file
# This functions cleans and converts to upper case
readGMTRobust <- function (file1) {   # size restriction
	# Read in the first file 
	x <- scan(file1, what="", sep="\n")
	# x <- gsub("\t\t.","",x)     # GMT files saved by Excel has a lot of empty cells "\t\t\t\t"   "\t." means one or more tab
	x <- gsub(" ","",x)  # remove white space
	x <- toupper(x)    # convert to upper case

	#----Process the first file
	# Separate elements by one or more whitespace
	y <- strsplit(x, "\t")
	# Extract the first vector element and set it as the list element name
	names(y) <- sapply(y, `[[`, 1)
	#names(y) <- sapply(y, function(x) x[[1]]) # same as above
	# Remove the first vector element from each list element
	y <- lapply(y, `[`, -c(1,2))
	#y <- lapply(y, function(x) x[-1]) # same as above
	# remove duplicated elements
	for ( i in 1:length(y) )  y[[i]] <- cleanGeneSet(y[[i]])
	# check the distribution of the size of gene lists sapply(y, length) hold a vector of sizes
	if( max( sapply(y,length) ) <5) cat("Warning! Gene sets have very small number of genes!\n Please double check format.")
	y <- y[which(sapply(y,length) > 1)]  # gene sets smaller than 1 is ignored!!!

	return(y)
}

# Find correspond gmt category
DB.Manager$set("public", "gmtCategory",
    function(converted, convertedData, selectOrg, gmtFile) {
        if(selectOrg == "NEW" && !is.null(gmtFile) )
            return( list(Custom_GeneSet ="Custom" ) )
        idNotRecognized = as.data.frame("ID not recognized!")
        if(is.null(converted) ) return(idNotRecognized) # no ID 
        querySet <- rownames(convertedData)
        if(length(querySet) == 0) return(idNotRecognized )
        ix = grep(converted$species[1,1],self$gmtFiles)
        if (length(ix) == 0 ) {return(idNotRecognized )}
        
        speciesChoice <- self$GetSpeciesChoice()
        # If selected species is not the default "bestMatch", use that species directly
        if(selectOrg != speciesChoice[[1]]) {  
            ix = grep(findSpeciesById(selectOrg)[1,1], self$gmtFiles )
            if (length(ix) == 0 ) {return(idNotRecognized )}
        }
        pathway <- dbConnect(self$sqlite,self$gmtFiles[ix],flags=SQLITE_RO)
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
        return(categoryChoices)
    }
)
