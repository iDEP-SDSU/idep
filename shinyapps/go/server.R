library(shiny)
library(RSQLite)
library(ggplot2)
#library(grid)
library(gridExtra)

# relative path to data files
datapath = "../../data/data92/"   # production server

Min_overlap <- 2
minSetSize = 3;
mappingCoverage = 0.60 # 60% percent genes has to be mapped for confident mapping
mappingEdge = 0.5  # Top species has 50% more genes mapped
maxTerms =30 # max number of enriched terms
PvalGeneInfo = 0.05; minGenes = 10 # min number of genes for ploting
PvalGeneInfo1 = 0.01
PvalGeneInfo2 = 0.001
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

STRING10_species = read.csv("STRING10_species.csv")

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
        "square"), legend.parameter = list("topright"), show.legend = TRUE, plotting=TRUE, layoutButton = 0,
    ...) 
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
    pvalues <- x[, pvalue]
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
                V(g)$color[index] <- sapply(1 - (group.pvalues - 
                  min(group.pvalues))/(max(group.pvalues) - min(group.pvalues)), 
                  function(x) {
                    adjustcolor(group.color[i], alpha.f = x)
                  })
            }
        }
    }
	if(plotting) { 
		plot(g,, vertex.label.dist=0.8, ...)
		if (show.legend) {
			legend.parameter$legend <- group.level
			legend.parameter$text.col <- group.color
			legend.parameter$bty <- "n"	
			do.call(legend, legend.parameter)
		}}
    return(g)
}


enrichmentNetwork <- function(enrichedTerms,layoutButton=0){
	geneLists = lapply(enrichedTerms$Genes, function(x) unlist( strsplit(as.character(x)," " )   ) )
	names(geneLists)= enrichedTerms$Pathways
	enrichedTerms$Direction = gsub(" .*","",enrichedTerms$Direction )

	g <- enrich.net2(enrichedTerms, geneLists, node.id = "Pathways", numChar = 100, 
	   pvalue = "adj.Pval", edge.cutoff = 0.2, pvalue.cutoff = 1, degree.cutoff = 0,
	   n = 200, group = enrichedTerms$Direction, vertex.label.cex = 1, vertex.label.color = "black", show.legend=FALSE, layoutButton=layoutButton)

}

keggSpeciesID = read.csv(paste0(datapath,"data_go/KEGG_Species_ID.csv"))

# finds id index corresponding to entrez gene and KEGG for id conversion
idType_Entrez <- dbGetQuery(convert, paste("select distinct * from idIndex where idType = 'entrezgene'" ))
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
 
shinyServer(
  function(input, output,session){
    options(warn=-1)

    observe({  updateSelectInput(session, "selectOrg", choices = speciesChoice )      })

	 # update species for STRING-db related API access
	 
 # tried to solve the double reflashing problems	 
 #https://stackoverflow.com/questions/30991900/avoid-double-refresh-of-plot-in-shiny
 #observe({  	updateSelectInput(session, "speciesName", choices = sort(STRING10_species$official_name) ) 	})
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

significantOverlaps <- reactive({
	  if (input$goButton == 0 | is.null( input$selectGO) ) return()
	  isolate({ 
	  withProgress(message= sample(quotes,1),detail="enrichment analysis", {
  	  #gene info is passed to enable lookup of gene symbols
  	  tem = geneInfoLookup(); tem <- tem[which( tem$Set == "List"),]
  	  FindOverlap( converted(), tem, input$selectGO,input$selectOrg,input$minFDR )
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
  
output$EnrichmentTable <-renderTable({
      if (input$goButton == 0  )    return(NULL)

	  myMessage = "Those genes seem interesting! Let me see what I can do.
	   I am comparing your query genes to all 150+ types of IDs across 111 species.
	  This can take up to 3 years. "
	  if(is.null(significantOverlaps() ) ) return(NULL)
	  withProgress(message= sample(quotes,1),detail=myMessage, {
	  tem <- significantOverlaps();
	  incProgress(1, detail = paste("Done"))	  })

	  if(dim(tem$x)[2] >1 ) tem$x[,2] <- as.character(tem$x[,2])
	  if(dim(tem$x)[2] ==1 ) tem$x else tem$x[,1:4]  # If no significant enrichment found x only has 1 column.
    }, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)

significantOverlaps2 <- reactive({
    if (input$goButton == 0  )    return()

    tem <- significantOverlaps();
    if(dim(tem$x)[2] ==1 ) return(NULL)
	tem <- tem$x;
    colnames(tem)= c("adj.Pval","nGenesList","nGenesCategor","Pathways","Genes")
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
    if(is.null(significantOverlaps2())) return(NULL)

	enrichmentNetwork(significantOverlaps2(),layoutButton = input$layoutButton )

}, height=900)	  
 
output$enrichmentNetworkPlotDownload <- downloadHandler(
      filename = "enrichmentPlotNetworkPathway.tiff",
      content = function(file) {
	  tiff(file, width = 12, height = 12, units = 'in', res = 300, compression = 'lzw')
	  enrichmentNetwork(significantOverlaps2(),layoutButton = input$layoutButton )
        dev.off()
      })	  
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
		string_db <- STRINGdb$new( version="10", species=taxonomyID,
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
      if (input$goButton == 0  )    return(NULL)

	library(STRINGdb,verbose=FALSE)

	tem = input$STRINGdbGO
		taxonomyID = findTaxonomyID(  )
		if(is.null( taxonomyID ) ) return(NULL)		
		####################################

		#if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
		NoSig = as.data.frame("No significant enrichment found.")
		if(is.null(STRINGdb_geneList() ) ) return(NULL)
		
		isolate({
		withProgress(message=sample(quotes,1), detail ="Enrichment analysis", {
		#Intialization
		string_db <- STRINGdb$new( version="10", species=taxonomyID,
							   score_threshold=0, input_directory="" )
				
		# using expression data

		genes <- conversionTableData()
		#rownames(genes)= genes[,3]
		minGenesEnrichment=3
		if(is.null(genes) ) return(NULL) 

		if(dim(genes)[1] <= minGenesEnrichment ) return(NoSig) # if has only few genes
		
		fc = rep(1, dim(genes)[1] )		
		# GO
		results1 <- NULL; result <- NULL
		pp <- 0
		for( i in c(1) ) {
			#incProgress(1/2,detail = paste("Mapping gene ids")  )
			ids = STRINGdb_geneList()[[i]]
			if( length(ids) <= minGenesEnrichment) next; 			
			incProgress(1/3  )
			result <- string_db$get_enrichment( ids, category = input$STRINGdbGO, methodMT = "fdr", iea = TRUE )
			if(nrow(result) == 0 ) next; 
			if(nrow(result) > 30)  result <- result[1:30,]

			if( dim(result)[2] ==1) next;   # result could be NULL
			#if(i == 1) result$direction = "Up regulated"  else result$direction = "Down regulated"
			if (pp==0 ) { results1 <- result; pp = 1;} else  results1 = rbind(results1,result)
		}

		if ( pp == 0 ) return (NoSig)

		if ( is.null( results1) ) return (NoSig)
	
		if( dim(results1)[2] == 1 ) return(NoSig)  # Returns a data frame: "No significant results found!"
		
		results1= results1[,c(5,3, 6)]
		colnames(results1)= c("FDR","nGenes","GO terms or pathways")
		minFDR = 0.01

		if(min(results1$FDR) > minFDR ) results1 = as.data.frame("No signficant enrichment found.") else
		results1 = results1[which(results1$FDR < minFDR),]
		
		incProgress(1, detail = paste("Done")) 
		
		if(dim(results1)[2] != 3) return(NoSig)
		colnames(results1)= c("adj.Pval","nGenes","Pathways")
		results1$adj.Pval <- sprintf("%-2.1e",as.numeric(results1$adj.Pval) )	
		rownames(results1)=1:nrow(results1)
		#results1[ duplicated (results1[,1] ),1 ] <- ""  
		
		return( results1 )
		 })#progress
		}) #isolate						   

}) 

output$stringDB_GO_enrichment <- renderTable({
		if(is.null(stringDB_GO_enrichmentData() ) ) return(NULL)

		 stringDB_GO_enrichmentData()	   

}, digits = 0,spacing="s",include.rownames=F,striped=TRUE,bordered = TRUE, width = "auto",hover=T) 

output$STRING_enrichmentDownload <- downloadHandler(
		filename = function() {paste0("STRING_enrichment",input$STRINGdbGO,".csv")},
		content = function(file) {
			write.csv(stringDB_GO_enrichmentData(), file)
	    }
	) 


  
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
		string_db <- STRINGdb$new( version="10", species=taxonomyID,
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
		string_db <- STRINGdb$new( version="10", species=taxonomyID,
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
	  
	  selectInput("selectGO", label=NULL,
		choices=gmtCategory(converted(), input$selectOrg),
	    selected = "GOBP" )    	
		

		
		
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
	  tem=input$selectOrg; 
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

# barplots using R base graphics
output$genePlot <- renderPlot({
	   if (input$goButton == 0  )    return()
	  tem=input$selectOrg; 
	isolate( {
	  withProgress(message="Ploting gene characteristics", {
       x = geneInfoLookup()
	   x2 = x[which(x$gene_biotype == "protein_coding"),]  # only coding for some analyses
     if(dim(x)[1]>=minGenes) # only making plots if more than 20 genes
       { # only plot when there 10 genes or more   # some columns have too many missing values
	   par(mfrow=c(4,1))
	   par(mar=c(8,6,8,2))
   	   #chromosomes
	   if( sum(!is.na( x$chromosome_name) ) >= minGenes && length(unique(x$chromosome_name) ) > 2 && length(which(x$Set == "List") ) > minGenes )
	   {
		   freq = table( x$chromosome_name,x$Set );
		   freq <- as.matrix(freq[which(nchar(row.names(freq))<3   ),])# remove unmapped chromosomes
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
		tem =gsub("processed_pseudogene","proc_pseudo",tem)
	    tem =gsub("processed","proc",tem); #row.names(freq)= tem
		par(mar=c(20,6,4.1,2.1))
		freq = freq[,c(2,1)] # reverse order
		head(freq)
		
        barplot(t(freq), beside=TRUE,las=2,col=c("red","lightgrey"), ylab="Number of Genes",
	      main= sig,cex.lab=1.5, cex.axis= 1.5,cex.names=1.5, cex.main=1.5)

	    legend("topright", c("List","Expected"), pch=15, col=c("red","lightgrey"),bty="n", cex=2)
		if(0) { 
        plt = barplot(t(freq), beside=TRUE,las=2,col=c("red","lightgrey"), ylab="Number of Genes",
	      main= sig,cex.lab=1.5, cex.axis= 2, cex.names=1.5, cex.main=1.5, str=45, adj=1, xpd=TRUE,xaxt="n" )
		text( plt,par("usr")[3], labels = rownames(freq), srt = 45, adj=c(1.1,1.1), xpd=TRUE ,cex = 1.5  )
		}


		}
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
		}}
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
       } }
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
     if(dim(x)[1]>=minGenes) # only making plots if more than 20 genes
       { # only plot when there 10 genes or more   # some columns have too many missing values
	  # par(mfrow=c(10,1))
	   # par(mar=c(8,6,8,2))
	   
	  # increase fonts
	  theme_set(theme_gray(base_size = 25)) 
	   
      #Coding Sequence length 
	  if( sum(!is.na( x2$cds_length) ) >= minGenes && length(unique(x2$cds_length) ) > 2 && length(which(x2$Set == "List") ) > minGenes) {
	   Pval = t.test(log(cds_length)~Set, data=x2 )$p.value
	   sig = paste("P = ",formatC(Pval, digits=2, format="G"),sep="")
	   if( Pval <PvalGeneInfo2)  sig = paste(sig," ***" ) else 
			if( Pval <PvalGeneInfo1)  sig = paste(sig," **" ) else 
			if( Pval <PvalGeneInfo)  sig = paste(sig," *" )  			
			
	   p1 <- ggplot(x2, aes(cds_length, fill= Set, colour = Set) )+
			geom_density(alpha = 0.1) + 
			scale_x_log10() +
			labs(x = "Coding sequence length (bp)") +
			annotate("text",x= min(x2$cds_length)+50, y = .5, label=sig, size=8)				
       }

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
				annotate("text",x= min(x2$cds_length)+100, y = .5, label=sig, size=8)+	
				labs(x = "Transcript length (bp)")		   
		  }
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
			labs(x = "Genome span (bp)")		   
	   }			  

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
			labs(x = "5' UTR length(bp)")
	}	

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
			labs(x = "3' UTR length(bp)")
	  }	
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
				annotate("text",x= min(x2$percentage_gc_content)+5, y = .02, label=sig, size=8)+	
				labs(x = "GC content (%)")		   
	   }		  
		  
	   incProgress(1/8)	
	   grid.arrange(p1,p2,p3,p4,p5,p6, ncol=1)
	   
	   
		}
		 incProgress(1/8, detail = paste("Done"))	  })
	 }) #isolate
    }, width=700,height = 3000)

	
output$listSigPathways <- renderUI({
	tem = input$selectOrg
	if (input$goButton == 0 | is.null(significantOverlaps())) return(NULL)

	tem <- significantOverlaps();

	if(dim(tem$x)[2] ==1 ) return(NULL)  
	choices = tem$x[,4]		
	selectInput("sigPathways", label="Select a KEGG pathway to show diagram with query genes highlighted in red:"
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
    new.signature = TRUE, plot.col.key = TRUE, key.align = "x", 
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

	 fold = rep(1, length(converted()$IDs))
	 names(fold) <- converted()$IDs
	 Species <- converted()$species[1,1]
	 fold <- convertEnsembl2Entrez(fold,Species)
	 
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
	
})
