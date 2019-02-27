library(shiny)
library(RSQLite)
library(gplots)
library(ggplot2)
library(e1071) # computing kurtosis
library(limma) # Differential expression
library(DESeq2) # count data analysis
library(edgeR) # count data D.E.
library(gage) # pathway analysis
library(PGSEA) # pathway 
library(fgsea) # fast GSEA
library(ReactomePA) # pathway analysis
library(pathview)
library(reshape2) # for melt correlation matrix in heatmap
library(DT) # for renderDataTable

Min_overlap <- 2
minSetSize = 3; 
mappingCoverage = 0.60 # 60% percent genes has to be mapped for confident mapping
mappingEdge = 0.5  # Top species has 50% more genes mapped
PvalGeneInfo = 0.05; minGenes = 10 # min number of genes for ploting
kurtosis.log = 50  # log transform is enforced when kurtosis is big
kurtosis.warning = 10 # log transformation recommnded 
minGenesEnrichment = 2 # perform GO or promoter analysis only if more than this many genes

# this need to be removed. Also replace go to go for folder
#setwd("C:/Users/Xijin.Ge/Google Drive/research/Shiny/RNAseqer")
sqlite  <- dbDriver("SQLite")
convert <- dbConnect(sqlite,"../go/convertIDs.db")
set.seed(2)
mycolors = sort(rainbow(20))[c(1,20,10,11,2,19,3,12,4,13,5,14,6,15,7,16,8,17,9,18)] # 20 colors for kNN clusters

hclust2 <- function(x, method="average", ...)
  hclust(x, method=method, ...)
dist2 <- function(x, ...)
  as.dist(1-cor(t(x), method="pearson"))

keggSpeciesID = read.csv("KEGG_Species_ID.csv")  
 detectGroups <- function (x){  # x are col names
	tem <- gsub("[0-9]*$","",x) # Remove all numbers from end
	#tem = gsub("_Rep|_rep|_REP","",tem)
	tem <- gsub("_Rep$","",tem); # remove "_Rep" from end
	tem <- gsub("_rep$","",tem); # remove "_rep" from end
	tem <- gsub("_REP$","",tem)  # remove "_REP" from end
 	return( tem )
 }
 myheatmap <- function (x,n=-1) {
if(n == -1) n=dim(x)[1]
geneSD = apply(x,1,sd)
x = x[order(-geneSD),]
# this will cutoff very large values, which could skew the color 
x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
x[x>cutoff] <- cutoff
cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
x[x< cutoff] <- cutoff

hy <-  heatmap.2(x, distfun = dist2,hclustfun=hclust2,
 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
#,Colv=FALSE,
,key=F
,margins = c(6, 8)
)
}
myheatmap3 <- function (x,n=-1) {

if( n > 0 && n< dim(x)[1]   )
 { geneSD = apply(x,1,sd)
 x = x[order(-geneSD),]
# this will cutoff very large values, which could skew the color 
  x=as.matrix(x[1:n,])
  }
    x=as.matrix(x)-apply(x,1,mean)
cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
x[x>cutoff] <- cutoff
cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
x[x< cutoff] <- cutoff
groups = detectGroups(colnames(x))
colnames(x) = detectGroups(colnames(x))
hy <-  heatmap.2(x, dendrogram ="row",distfun = dist2,hclustfun=hclust2,
 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
#,Colv=FALSE,
,labRow=""
#,labCol=""
,ColSideColors=mycolors[ groups]
,key=F
,margins = c(6, 20)
)
if(0) {
	lmat = rbind(c(5,4),c(0,1),c(3,2))
	lwid = c(1.5,6)
	lhei = c(1,.2,8)



if( dim(x)[1]>100) 
	heatmap.2(x, distfun = dist2,hclustfun=hclust2,
	 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
	,key=T, symkey=F
	,ColSideColors=mycolors[ groups]
	,labRow=""
	,margins=c(6,8)
	,srtCol=45
	#,lmat = lmat, lwid = lwid, lhei = lhei
	)

	if( dim(x)[1] <=100) 
	heatmap.2(x, distfun = dist2,hclustfun=hclust2,
	 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
	,key=T, symkey=F,
	#,labRow=labRow
		,ColSideColors=mycolors[ groups]
	,margins=c(6,8)
	,cexRow=1.5
	,srtCol=45
	#,lmat = lmat, lwid = lwid, lhei = lhei
	)
  }

}

# randomly samples genes
myheatmap4 <- function (x,n=-1) {

if( n > 0 && n< dim(x)[1]   )
 { 
  ix = sample(1:dim(x)[1], n)
# this will cutoff very large values, which could skew the color 
  x=as.matrix(x[ix,])
  }
    x=as.matrix(x)-apply(x,1,mean)
cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
x[x>cutoff] <- cutoff
cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
x[x< cutoff] <- cutoff
groups = detectGroups(colnames(x))
colnames(x) = detectGroups(colnames(x))
hy <-  heatmap.2(x, dendrogram ="row",distfun = dist2,hclustfun=hclust2,
 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
,Colv=FALSE,
,labRow=""
#,labCol=""
,ColSideColors=mycolors[ groups]
,key=F
,margins = c(6, 20)
)

}


myheatmap2 <- function (x,bar,n=-1 ) {
# number of genes to show
ngenes = as.character( table(bar))
if(length(bar) >n && n != -1) {ix = sort( sample(1:length(bar),n) ); bar = bar[ix]; x = x[ix,]  }

# this will cutoff very large values, which could skew the color 
x=as.matrix(x)-apply(x,1,mean)
cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
x[x>cutoff] <- cutoff
cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
x[x< cutoff] <- cutoff
colnames(x)= detectGroups(colnames(x))
 heatmap.2(x,  Rowv =F,Colv=F, dendrogram ="none",
 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
,key=F, labRow = F,
,RowSideColors = mycolors[bar]
,margins = c(8, 24)
,srtCol=45
)

legend.text = paste("Cluster ", toupper(letters)[unique(bar)], " (N=", ngenes,")", sep="")
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
    legend = legend.text, # category labels
    col = mycolors,  # color key
    lty= 1,             # line style
    lwd = 10            # line width
)

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

convertEnsembl2Entrez <- function (query,Species) { 
  querySet <- cleanGeneSet( unlist( strsplit( toupper(names( query)),'\t| |\n|\\,' )  ) )
  speciesID <- orgInfo$id[ which(orgInfo$ensembl_dataset == Species)]  # note uses species Identifying
  # idType 6 for entrez gene ID
  result <- dbGetQuery( convert,
                       paste( " select  id,ens,species from mapping where ens IN ('", paste(querySet,collapse="', '"),
                			"') AND  idType ='6'",sep="") )	# slow
						
  if( dim(result)[1] == 0  ) return(NULL)
  result <- subset(result, species==speciesID, select = -species)

  ix = match(result$ens,names(query)  )

  tem <- query[ix];  names(tem) = result$id
  return(tem)
  
}

convertEnsembl2KEGG <- function (query,Species) {  # not working
  querySet <- cleanGeneSet( unlist( strsplit( toupper(names( query)),'\t| |\n|\\,' )  ) )
  speciesID <- orgInfo$id[ which(orgInfo$ensembl_dataset == Species)]  # note uses species Identifying
  # idType 6 for entrez gene ID
  result <- dbGetQuery( convert,
                       paste( " select  id,ens,species from mapping where ens IN ('", paste(querySet,collapse="', '"),
                			"') AND  idType ='107'",sep="") )	# slow
						
  if( dim(result)[1] == 0  ) return(NULL)
  result <- subset(result, species==speciesID, select = -species)

  ix = match(result$ens,names(query)  )

  tem <- query[ix];  names(tem) = result$id
  return(tem)  
}

geneInfo <- function (converted,selectOrg){
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
FindOverlap <- function (converted,gInfo, GO,selectOrg,minFDR) {
  maxTerms =10 # max number of enriched terms
  idNotRecognized = as.data.frame("ID not recognized!")
  
  if(is.null(converted) ) return(idNotRecognized) # no ID 
  
  # only coding
  gInfo <- gInfo[which( gInfo$gene_biotype == "protein_coding"),]  
  querySet <- intersect( converted$IDs, gInfo[,1]);
  
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
 
 
  if(min(x$FDR) > minFDR) x=as.data.frame("No significant enrichment found!") else {
  x <- x[which(x$FDR < minFDR),] 
  if(dim(x)[1] > maxTerms ) x = x[1:maxTerms,]
  x= cbind(x,sapply( x$pathwayID, sharedGenesPrefered ) )
  colnames(x)[7]= "Genes"
  x <- subset(x,select = c(FDR,overlap,n,description,Genes) )
  colnames(x) = c("Corrected P value (FDR)", "Genes in list", "Total genes in category","Functional Category","Genes"  )
  }
        
 dbDisconnect(pathway)
 return(x )
} 
                                     #, categoryChoices = categoryChoices 
#Given a pathway description, found gene ids
# slow, needs to add index for pathway ID when creating.
pathwayGenes <- function (pathwayDescription, converted, GO,selectOrg) {
  ix = grep(converted$species[1,1],gmtFiles)

  if (length(ix) == 0 ) {return(idNotRecognized )}
  
  # If selected species is not the default "bestMatch", use that species directly
  if(selectOrg != speciesChoice[[1]]) {  
    ix = grep(findSpeciesById(selectOrg)[1,1], gmtFiles )
	if (length(ix) == 0 ) {return(idNotRecognized )}
	totalGenes <- orgInfo[which(orgInfo$id == as.numeric(selectOrg)),7]
  }
  pathway <- dbConnect(sqlite,gmtFiles[ix])
 
  pathwayInfo <- dbGetQuery( pathway, paste( " select * from pathwayInfo where description =  '", 
						pathwayDescription,   "' AND name LIKE '",GO,"%'",sep="") )
  if(dim(pathwayInfo)[1] != 1 ) { cat("\nPathway not mapped, or not uniquely mapped"); dbDisconnect(pathway); return(NULL) }
  
  #pthwayInfo table: id name description, n, memo, golevel
  # id                                  name   description   n                                                                      memo golevel
 #  1 GOBP_gmx_ens_GO:0000003_reproduction  Reproduction  366  http://amigo.geneontology.org/cgi-bin/amigo/term_details?term=GO:0000003       2

  sqlQuery = paste( " select distinct gene from pathway where pathwayID = '", pathwayInfo[1,1],"'" ,sep="")
  result <- dbGetQuery( pathway, sqlQuery  )
	
  if( dim(result)[1] ==0) {dbDisconnect(pathway); return(NULL) }
   dbDisconnect(pathway)
  return( result[,1])
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
  pathway <- dbConnect(sqlite,gmtFiles[ix])
 
  pathwayInfo <- dbGetQuery( pathway, paste( " select * from pathwayInfo where description =  '", 
						pathwayDescription,   "' AND name LIKE '",GO,"%'",sep="") )
  dbDisconnect(pathway);
  if(dim(pathwayInfo)[1] != 1 ) {return(NULL) }
  tem = gsub(".*:","",pathwayInfo[1,2])  
  return( gsub("_.*","",tem) )
  } 
  
  
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
   #geneSets <- geneSets[ -which(duplicated(names(geneSets) ))] # remove geneSets with the same name
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
	
	result = result[which(result[,1] < Pval_pathway),,drop=F]
	#result = result[which(result[,2] >2)    ,]
	pg2 = result[,-2]

	# when there is only 1 left in the matrix pg2 becomes a vector
	if(sum( Pvalues<Pval_pathway) == 1) { pg3 = t( as.matrix(pg2));pg3 = rbind(pg3,pg3);} else
	{ if(dim(pg2)[1] > top ) {  pg3 = pg2[1:top,]; } else {  pg3 = pg2;  } }

	rownames(pg3) = sapply(rownames(pg3) , extract1)
	a=sprintf("%-1.0e",pg3[,1])
	rownames(pg3) = paste(a,rownames(pg3),sep=" ")
	pg3 =pg3[,-1]
	
	pg3 <- pg3[order( -apply(pg3,1,sd)    ),] # sort by SD
 
    return( list(pg3 = pg3, best = best ) )
    }
 }

DEG.limma <- function (x, maxP_limma=.1, minFC_limma=2, rawCounts,countsDEGMethods,priorCounts, dataFormat){
 topGenes = list();  limmaTrend = FALSE
 if( dataFormat == 2) {   # if normalized data
	eset = new("ExpressionSet", exprs=as.matrix(x)) } else { # counts data
		if (countsDEGMethods == 1 ) { # limma-trend method selected for counts data
			dge <- DGEList(counts=rawCounts);
			dge <- calcNormFactors(dge)
			eset <- cpm(dge, log=TRUE, prior.count=priorCounts)
			limmaTrend = TRUE
		}
 }

 groups = colnames(x)
 groups = detectGroups( groups)
 g =  unique(groups)  
 
 # check for replicates, removes samples without replicates
 reps = as.matrix(table(groups)) # number of replicates per biological sample
 if ( sum( reps[,1] >= 2) <2 ) # if less than 2 samples with replicates
 return( list(results= NULL, comparisons = NULL, Exp.type="Failed to parse sample names to define groups. 
       Cannot perform DEGs and pathway analysis. Please double check column names! Use WT_Rep1, WT_Rep2 etc. ", topGenes=NULL)) 
 # remove samples without replicates
 g <- rownames(reps)[which(reps[,1] >1)]
 ix <- which( groups %in% g)  
 groups <- groups[ix]   
 x<- x[,ix]; rawCounts <- rawCounts[,ix] 
 
 
 if(length(g) ==2 ) { 
  g= unique(groups)
 comparisons <-  paste(g[2],"-",g[1],sep="")
 design <- model.matrix(~0+groups)
 colnames(design) <- g
 
 if( !is.null(rawCounts) && countsDEGMethods == 2) {  # voom
     v <- voom(rawCounts, design); fit <- lmFit(v, design) } else 
	fit <- lmFit(eset, design)      # regular limma
	
cont.matrix <- makeContrasts(contrasts=comparisons, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, trend=limmaTrend)

 # calls differential gene expression 1 for up, -1 for down
 results <- decideTests(fit2, p.value=maxP_limma, lfc=log2(minFC_limma) )
 #vennDiagram(results,circle.col=rainbow(5))
 topGenes1 =topTable(fit2, number = 1e12,sort.by="M" )
 if (dim(topGenes1)[1] != 0) {
 topGenes1 = topGenes1[,c('logFC','adj.P.Val')] 
 # topGenes1[,1] <-  -1* topGenes1[,1] # reverse direction
 topGenes[[1]] <- topGenes1 }
  # log fold change is actually substract of means. So if the data is natral log transformed, it shoudl be natral log.
 Exp.type = "2 biological samples."
 }
  
  if(length(g) > 2 ) { 
 design <- model.matrix(~ 0+factor(groups))
 colnames(design) <- gsub(".*)","",colnames(design))
 
 if( !is.null(rawCounts) && countsDEGMethods == 2) {  # voom
     v <- voom(rawCounts, design); fit <- lmFit(v, design) } else 
		fit <- lmFit(eset, design)
 
 fit <- eBayes(fit, trend=limmaTrend)
 
 comparisons = ""
  for( i in 1:(length(g)-1) )
    for (j in (i+1):length(g)) 
	 comparisons = c(comparisons,paste(g[j],"-",g[i],sep="" ) )
  comparisons <- comparisons[-1]

contrast1 <- makeContrasts(contrasts=comparisons[1], levels=design)
for( kk in 2:length(comparisons) )
     contrast1<-  cbind(contrast1,makeContrasts(contrasts=comparisons[kk], levels=design)   )
 Exp.type = paste(length(g)," samples detected.")	 
 
# if factorial design 2x2, 2x3, 3x5 etc.
	# all samples must be something like WT_control_rep1
if ( sum(sapply(strsplit(g,"_"),length) == 2 ) == length(g) ) {
	comparisons = ""
	  for( i in 1:(length(g)-1) )
		for (j in (i+1):length(g)) 
		 if( strsplit(g[i],"_")[[1]][1] == strsplit(g[j],"_")[[1]][1]| strsplit(g[i],"_")[[1]][2] == strsplit(g[j],"_")[[1]][2]) # only compare WT_control vs. WT_treatment
			comparisons = c(comparisons,paste(g[j],"-",g[i],sep="" ) )
	  comparisons <- comparisons[-1]
	 
	extract_treatment <- function (x) paste( gsub( ".*_","",unlist( strsplit(x,"-")) ), collapse="-")
	extract_genotype <- function (x) gsub( "_.*","",unlist( strsplit(x,"-")) )[1]
	extract_treatment_counting <- unique( gsub( ".*_","",unlist( strsplit(g,"-")) ))
	treatments = sapply(comparisons, extract_treatment)
	genotypes = sapply(comparisons, extract_genotype)
	Exp.type = paste( Exp.type, "\nFactorial design:",length(unique(genotypes)),"X", length( extract_treatment_counting ), sep="" )
	contrast1 <- makeContrasts(contrasts=comparisons[1], levels=design)
	for( kk in 2:length(comparisons) )
		 contrast1<-  cbind(contrast1,makeContrasts(contrasts=comparisons[kk], levels=design)   )
	contrast.names = colnames(contrast1)
	for ( kk in 1:(length(comparisons)-1) ) {
	   for( kp in (kk+1):length(comparisons)) 
		  if( treatments[kp]== treatments[kk] ) 
		   {  
			  contrast1 = cbind(contrast1, contrast1[,kp]- contrast1[,kk] )
			  contrast.names = c(contrast.names, paste("Diff:",  genotypes[kp], "-", genotypes[kk],"(",gsub("-",".vs.",treatments[kp]),")",sep="" ) )
		   }   
	}
	colnames(contrast1)=contrast.names
	comparisons = contrast.names
	}	
	 
fit2 <- contrasts.fit(fit, contrast1)
fit2 <- eBayes(fit2, trend=limmaTrend)
#topTable(fit2, coef=1, adjust="BH")
results <- decideTests(fit2, p.value=maxP_limma, lfc= log2(minFC_limma ))
#vennDiagram(results[,1:5],circle.col=rainbow(5))
 top <- function (comp) {
   tem <- topTable(fit2, number = 1e12,coef=comp,sort.by="M" ) 
   if(dim(tem)[1] == 0) return (1) else   return( tem[,c(1,5)])   # no reversing
                                           #{ tem <-tem[,c(1,5)]; tem[,1] <- -1*tem[,1]; return(tem)  }    # FC is reversed ???
                                           
 }  # no significant gene returns 1, otherwise a data frame
 topGenes <- lapply(comparisons, top)
 topGenes <- setNames(topGenes, comparisons )
 ix <- which( unlist( lapply(topGenes, class) ) == "numeric")
 if( length(ix)>0) topGenes <- topGenes[ - ix ]
  # if (length(topGenes) == 0) topGenes = NULL;
}
 return( list(results= results, comparisons = comparisons, Exp.type=Exp.type, topGenes=topGenes)) 
}

DEG.DESeq2 <- function (  rawCounts,maxP_limma=.05, minFC_limma=2){
 groups = as.character ( detectGroups( colnames( rawCounts ) ) )
 g = unique(groups)# order is reversed
 
 # check for replicates, removes samples without replicates
 reps = as.matrix(table(groups)) # number of replicates per biological sample
 if ( sum( reps[,1] >= 2) <2 ) # if less than 2 samples with replicates
 return( list(results= NULL, comparisons = NULL, Exp.type="Failed to parse sample names to define groups. 
       Cannot perform DEGs and pathway analysis. Please double check column names! Use WT_Rep1, WT_Rep2 etc. ", topGenes=NULL)) 
 # remove samples without replicates
 g <- rownames(reps)[which(reps[,1] >1)]
 ix <- which( groups %in% g)  
 groups <- groups[ix]   
 rawCounts <- rawCounts[,ix] 
  
 
 Exp.type = paste(length(g)," samples detected.")
  comparisons = ""
  for( i in 1:(length(g)-1) )
    for (j in (i+1):length(g)) 
	 comparisons = c(comparisons,paste(g[j],"-",g[i],sep="" ) )
  comparisons <- comparisons[-1]
  
colData = cbind(colnames(rawCounts), groups )

# Set up the DESeqDataSet Object and run the DESeq pipeline
dds = DESeqDataSetFromMatrix(countData=rawCounts,
                              colData=colData,
                              design=~groups)
dds = DESeq(dds)  # main function

result1 = NULL; allCalls = NULL;
topGenes = list(); pk = 1 # counter
pp=0 # first results?
for( kk in 1:length(comparisons) ) {
	tem = unlist( strsplit(comparisons[kk],"-") )
	selected = results(dds, contrast=c("groups", tem[1], tem[2]) ) #, lfcThreshold=log2(minFC_limma))
	# selected = subset(res, padj < maxP_limma )
	if(dim(selected)[1] == 0 ) next; # no significant genes
	selected = selected[order(-abs(selected$log2FoldChange)),] 

	selected$calls =0   
	selected$calls [which( selected$log2FoldChange > log2(minFC_limma) & selected$padj < maxP_limma ) ]  <-  1
	selected$calls [ which( selected$log2FoldChange <  -log2(minFC_limma) & selected$padj < maxP_limma ) ] <-  -1
	colnames(selected)= paste( as.character(comparisons[kk]), "___",colnames(selected),sep="" )
	selected = as.data.frame(selected)
	if (pp==0){  # if first one with significant genes, collect gene list and Pval+ fold
	  result1 = selected; pp = 1; 
	  # selected[,2] <- -1 * selected[,2] # reverse fold change direction
	  topGenes[[1]] = selected[,c(2,6)]; 
	  names(topGenes)[1] = comparisons[kk]; } else 
		  { result1 = merge(result1,selected,by="row.names"); 
		    rownames(result1) = result1[,1]; 
			result1 <- result1[,-1]
			pk= pk+1; 
			# selected[,2] <- -1 * selected[,2] # reverse fold change direction
			topGenes[[pk]] = selected[,c(2,6)]; 
			names(topGenes)[pk] = comparisons[kk]; 
		  }
}
#if( length(comparisons) == 1) topGenes <- topGenes[[1]] # if only one comparison, topGenes is not a list, just a data frame itself.
if(! is.null(result1)) { 
# note that when you only select 1 column from a data frame it automatically converts to a vector. drop =FALSE prevents that.
allCalls = as.matrix( result1[,grep("calls",colnames(result1)), drop = FALSE  ] )
colnames(allCalls)= gsub("___.*","", colnames(allCalls))
colnames(allCalls)= gsub("\\.","-", colnames(allCalls)) # note that samples names should have no "."
}
 return( list(results= allCalls, comparisons = comparisons, Exp.type=Exp.type, topGenes=topGenes)) 
}

promoter <- function (converted,selectOrg, radio){
  idNotRecognized = as.data.frame("ID not recognized!") 
  if(is.null(converted) ) return(idNotRecognized) # no ID 
  querySet <- converted$IDs;
  if(length(querySet) == 0) return(idNotRecognized )
    ix = grep(converted$species[1,1],motifFiles)

  # If selected species is not the default "bestMatch", use that species directly
  if(selectOrg != speciesChoice[[1]]) {  
    ix = grep(findSpeciesById(selectOrg)[1,1], motifFiles )
  }
  ix1 =grep(as.character(radio),motifFiles[ix]) # match 300bp or 600bp
  if(length(ix1) >0) ix = ix[ix1]   # if 600 is not found, use 300bp
  
   if (length(ix) == 0 ) {return(as.data.frame("No matching motif file found") )} else { 
  if(length(ix) > 1)  # if only one file          
    return(as.data.frame("Multiple geneInfo file found!") )   
 
 motifs <- dbConnect(sqlite,motifFiles[ix]) # makes a new file
    
  sqlQuery = paste( " select * from scores where row_names IN ('", paste(querySet,collapse="', '"),"')" ,sep="")
   result <- dbGetQuery( motifs, sqlQuery  )
  if( dim(result)[1] ==0) {return(as.data.frame("No matching species or gene ID file!" ) )}
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
   TFs$df <- ( TFs$scoreSD1^2/n1 + TFs$scoreSD^2/TFs$nGenes)^2 /
    (   (TFs$scoreSD1^2/n1)^2/(n1-1) +   (TFs$scoreSD^2/TFs$nGenes)^2/(TFs$nGenes-1)   )
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
colnames(TFs) =c("Motif", "TF","TF family","List","FDR","Score","Note"   )
if(dim(TFs)[1] >20 ) TFs <- TFs[1:20,]
if(dim(TFs)[1] ==0) return(as.data.frame("No significant TF binding motif detected.") ) else
return( TFs )
 }
}

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
            kegg.geneid = "K01488", ncbi.geneid = "")
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
            msg.fmt = "Only native KEGG gene ID is supported for this species,\nmake sure it looks like \"%s\"!"
            msg = sprintf(msg.fmt, species.data["kegg.geneid"])
            message("Note: ", msg)
        }
        else {
            stop("This species is not annotated in KEGG!")
        }
    }
    if (is.null(gene.annotpkg)) 
        gene.annotpkg = bods[match(species, bods[, 3]), 1]
    if (length(grep("ENTREZ|KEGG", gene.idtype)) < 1 & !is.null(gene.data)) {
        if (is.na(gene.annotpkg)) 
            stop("No proper gene annotation package available!")
        if (!gene.idtype %in% gene.idtype.bods[[species]]) 
            stop("Wrong input gene ID type!")
        gene.idmap = id2eg(gd.names, category = gene.idtype, 
            pkg.name = gene.annotpkg, unique.map = F)
        gene.data = mol.sum(gene.data, gene.idmap)
        gene.idtype = "ENTREZ"
    }
    if (gene.idtype == "ENTREZ" & !entrez.gnodes & !is.null(gene.data)) {
        message("Info: Getting gene ID data from KEGG...")
        gene.idmap = keggConv("ncbi-geneid", species)
        message("Info: Done with data retrieval!")
        kegg.ids = gsub(paste(species, ":", sep = ""), "", names(gene.idmap))
        ncbi.ids = gsub("ncbi-geneid:", "", gene.idmap)
        gene.idmap = cbind(ncbi.ids, kegg.ids)
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
    tfiles.xml = paste(pathway.name, ".xml", sep = "")
    tfiles.png = paste(pathway.name, ".png", sep = "")
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
            pv.pars = my.keggview.native( plot.data.gene = plot.data.gene, 
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
}
# <environment: namespace:pathview>
my.keggview.native <- function ( plot.data.gene = NULL, plot.data.cpd = NULL, cols.ts.gene = NULL, 
    cols.ts.cpd = NULL, node.data, pathway.name, out.suffix = "pathview", 
    kegg.dir = ".", multi.state = TRUE, match.data = TRUE, same.layer = TRUE, 
    res = 300, cex = 0.25, discrete = list(gene = FALSE, cpd = FALSE), 
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
        img.file = paste(kegg.dir,"/",pathway.name, ".",pn.suffix[np], ".png", 
            sep = "")
		message("here:",img.file)
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

 if(0 ){ # pathway testing
 x = read.csv("expression.csv")
 x = read.csv("expression1_no_duplicate.csv")
 x = read.csv("mouse1.csv")
 x = read.csv("GSE40261.csv")
 x = read.csv("GSE52778_All_Sample_FPKM_Matrix.csv")
 x = read.csv("exampleData/GSE87194.csv")
  x = read.csv("exampleData/expression_3groups.csv")
 x = x[order(x[,1]),]
    x = x[!duplicated(x[,1]),]
    rownames(x)= x[,1]
    x = x[,-1]
	
	tem = apply(x,1,max)
	x = x[which(tem> 1),] 
	
	x = log(x+abs( 1),2)
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
	
	pg = myPGSEA (x,cl=gmt,range=myrange,p.value=TRUE, weighted=FALSE,nPermutation=1)
	result = PGSEApathway (converted,convertedData, selectOrg,GO,gmt, myrange,.05,30)
	smcPlot(result$pg3,factor(subtype),scale = c(-result$best, result$best), show.grid = T, margins = c(3,1, 13, 23), col = .rwb,cex.lab=0.5)
	smcPlot(result$pg3,factor(subtype),scale = c(-max(result$pg3), max(result$pg3)), show.grid = T, margins = c(3,1, 13, 23), col = .rwb,cex.lab=0.5)
	
	pca = 100*prcomp(t(x))$rotation 
	Npca = 10
	if (Npca > dim(pca)[2]) { Npca = dim(pca)[2] } else pca <-  pca[,1:Npca]
	#pca = pca[,1:5]
	pg = myPGSEA (pca,cl=gmt,range=myrange,p.value=TRUE, weighted=FALSE,nPermutation=1)
	
	# correcting for multiple testing
	p.matrix = pg$p.result
	tem = p.adjust(as.numeric(p.matrix),"fdr")
	p.matrix = matrix(tem, nrow=dim(p.matrix)[1], ncol = dim(p.matrix)[2] )
	rownames(p.matrix) = rownames(pg$p.result); colnames(p.matrix) = colnames(pg$p.result)

	# using absolute value to rank 
 	#selected = unlist( apply(pg$result, 2, function(y) which( rank(y) >= length(y)-3.1)   ) )
    
	# using p value to rank #
	#selected = unlist( apply(p.matrix, 2, function(x) which( rank(x,ties.method='first') <= 5)   ) )
	selected =c()
	for( i in 1:dim(p.matrix)[2]) {
	 tem = which( rank(p.matrix[,i],ties.method='first') <= 3) 
	 #tem = which( rank(pg$result[,i],ties.method='first') >= dim(p.matrix)[1]-3.1)
	 names(tem) = paste("PC",i," ", rownames(p.matrix)[tem], sep="" )
	 selected = c(selected, tem)
	}
	rowids = gsub(" .*","",names(selected))
	rowids = as.numeric( gsub("PC","",rowids) )
	pvals = p.matrix[ cbind(selected,rowids) ]
	a=sprintf("%-1.0e",pvals)
	tem = pg$result[selected,]
	rownames(tem) = paste(a,names(selected)); #colnames(tem)= paste("PC",colnames(tem),sep="")
	
	tem = tem[!duplicated(selected),] 
	#tem = t(tem); tem = t( (tem - apply(tem,1,mean)) ) #/apply(tem,1,sd) )
	smcPlot(tem,scale =  c(-max(tem), max(tem)), show.grid = T, margins = c(3,1, 13, 23), col = .rwb,cex.lab=0.5)
	
	############  testing D.E.G.
	limma = DEG.limma(convertedData, .1, .5,rawCounts=NULL,countsDEGMethods=1,priorCounts=3, dataFormat=2)
	genes = limma$results
	  if( is.null(genes) ) return(NULL)
	  ix = match(limma$comparisons, colnames(genes)) 
	  query = rownames(genes)[which(genes[,ix] != 0)]
      iy = match(query, rownames(convertedData  ) )
	  convertedData[iy,]
      iz= match( detectGroups(colnames(convertedData)), unlist(strsplit( limma$comparisons, "-"))	  )
      iz = which(!is.na(iz))
      myheatmap( convertedData[iy,iz] )
	 # convertedData()[iy,iz]
 
  	  # rawCounts = read.csv("exampleData/airway_GSE52778.csv", row.names=1)
	   x = read.csv("exampleData/GSE87194.csv") ; x[,1] = toupper(x[,1]);   x = x[order(x[,1]),];    x = x[!duplicated(x[,1]),] #rownames(x)= x[,1]; rawCounts= x[,-1]
	   rawCounts = rawCounts[which(apply(rawCounts,1,max )>10 ) ,]
	 # res =DEG.DESeq2(rawCounts, .05, 2)
	 # res = 
	 # tem = res$topGenes
	 # head( res$results )
	 
	 # res2= DEG.limma(rawCounts, .05, 2,rawCounts, 1 ,3) 
 
 
 
 
 
  ######### testing GAGE
   fc = apply(x1[,4:6],1,mean)- apply(x1[,1:3],1,mean)
   paths <- gage(fc, gsets = gmt, ref = NULL, samp = NULL)
   
   paths <- as.data.frame(paths)
   path1 <- rownames(paths)[1]
   
  
   
	  x = read.csv("exampleData/airway_GSE52778.csv", row.names=1)
	 #x = read.csv("exampleData/GSE87194.csv") ; 
	 x=read.csv("GSE37704_sailfish_genecounts.csv");
	 #x = read.csv("exampleData/counts_test_data_3groups.csv") ;
	 x = read.csv("exampleData/hoppe 2 samples.csv")
	 
	 x[,1] = toupper(x[,1]);  
	 colnames(x)[1]= "User_input"

	
	selectOrg = "BestMatch"; GO="KEGG"; 
	myrange = c(15,2000)

	converted = convertID(x[,1],selectOrg)
	mapping = converted$conversionTable
	x = merge(mapping[,1:2],x,   by = 'User_input')
	tem = apply(x[,3:(dim(x)[2]-2)],1,sum)
	x = x[order(x[,2],-tem),]
	x = x[!duplicated(x[,2]) ,]
	rownames(x) = x[,2]
	x = as.matrix(x[,c(-1,-2)])
	
	convertedData = x
    gmt = readGeneSets(converted, convertedData, GO,selectOrg, myrange)

    res =DEG.DESeq2(x, .25, 1)
	
	res <- DEG.limma (x, maxP_limma=.2, minFC_limma=2, x,countsDEGMethods=2,priorCounts=3, dataFormat=1)

    top1 <- res$topGenes[[1]]
	
    head(top1)	
	
	 paths <- gage(top1[,1,drop=F], gsets = gmt, ref = NULL, samp = NULL)
	  paths <-  rbind(paths$greater,paths$less)
	  if(dim(paths)[1] < 1 | dim(paths)[2]< 6 ) return( noSig )
	  top1 <- paths[,c('stat.mean','set.size','q.val')]
	  colnames(top1)= c("stat.mean","Set Size","FDR")
	  top1 <- top1[order(top1[,3]) ,]  
	  if ( length( which( top1[,3] <=  .9   ) ) == 0 )
	    return( noSig)
	  top1 <- top1[which(top1[,3] <=  .9 ) ,]
	  if(dim(top1)[1] > 30 ) 
	     top1 <- top1[1:30,]
	  top1
	
 
   	res = DEG.limma(x, .05, 2,NULL, 1,3 )

 ## fgsea
     top1 <- res$topGenes[[1]]
    head(top1)	
	colnames(top1)= c("Fold","FDR")
	fold = top1[,1]; names(fold) <- rownames(top1)
 	 paths <- fgsea(pathways = gmt, 
                  stats = fold,
                  minSize=15,
                  maxSize=2000,
                  nperm=10000)
 if(dim(paths)[1] < 1  ) return( noSig )
      paths <- as.data.frame(paths)
	  top1 <- paths[,c(4,5,7,3)]
	  rownames(top1) <- paths[,1]
	  colnames(top1)= c("ES","NES","Set Size","FDR")
	  top1 <- top1[order(top1[,4]) ,]  
	  if ( length( which( top1[,4] <=  input$pathwayPvalCutoff   ) ) == 0 )
	    return( noSig)
	  top1 <- top1[which(top1[,4] <=  input$pathwayPvalCutoff ) ,]
	  if(dim(top1)[1] > input$nPathwayShow ) 
	     top1 <- top1[1:input$nPathwayShow,]
		 
	  top1
	  
 # testing visualize KEGG pathway
  
    query = x[1:500,1]
	 fc= convertEnsembl2Entrez (query, Species)  
	 
	 fc = log2(fc/mean(fc))
	 
	      top1 <- res$topGenes[[1]]
	 top1 <- top1[,1]; names(top1)= rownames(res$topGenes[[1]] )
	  Species = converted$species[1,1] 	
	  
	system.time(  fc <- convertEnsembl2Entrez (top1, Species)  )
		 fc = sort(fc,decreasing =T)
	 
       head(fc)
	 
	  system.time ( y<- gsePathway(fc, nPerm=1000,
                minGSSize=15, pvalueCutoff=0.5,
                pAdjustMethod="BH", verbose=FALSE) )
      res <- as.data.frame(y)
     head(res)
	  
	  
	  # testing mouse 
	  top1 = limma$topGenes[[1]]
	   top1 <- top1[,1]; names(top1)= rownames(limma$topGenes[[1]] )
	    Species = converted$species[1,1] 	
	   	system.time(  fc <- convertEnsembl2Entrez (top1, Species)  )
			 fc = sort(fc,decreasing =T)
			  system.time ( y<- gsePathway(fc, nPerm=1000,organism = "mouse",
                minGSSize=15, pvalueCutoff=0.5,
                pAdjustMethod="BH", verbose=FALSE) )
      res <- as.data.frame(y)
     head(res)
		ensemblSpecies <- c("hsapiens_gene_ensembl","rnorvegicus_gene_ensembl", "mmusculus_gene_ensembl",
		"celegans_gene_ensembl","scerevisiae_gene_ensembl", "drerio_gene_ensembl", "dmelanogaster_gene_ensembl")
         ReactomePASpecies= c("human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly" )

	
  #### testing KEGG pathway graph
 x=read.csv("GSE37704_sailfish_genecounts.csv");
	 #x = read.csv("exampleData/counts_test_data_3groups.csv") ;
	 x[,1] = toupper(x[,1]);  
	 colnames(x)[1]= "User_input"

	
	selectOrg = "BestMatch"; GO="KEGG"; 
	myrange = c(15,2000)

	converted = convertID(x[,1],selectOrg)
	mapping = converted$conversionTable
	x = merge(mapping[,1:2],x,   by = 'User_input')
	tem = apply(x[,3:(dim(x)[2]-2)],1,sum)
	x = x[order(x[,2],-tem),]
	x = x[!duplicated(x[,2]) ,]
	rownames(x) = x[,2]
	x = as.matrix(x[,c(-1,-2)])
	
	convertedData = x
    gmt = readGeneSets(converted, convertedData, GO,selectOrg, myrange)

    res =DEG.DESeq2(x, .25, 1)
	
    top1 <- res$topGenes[[1]]
	
    head(top1)	
	
	 paths <- gage(top1[,1,drop=F], gsets = gmt, ref = NULL, samp = NULL)
	  paths <-  rbind(paths$greater,paths$less)
	  
	selectedPathway = rownames(paths)[1]
   # [1] "Cytokine-cytokine receptor interaction"

	    Species <- converted$species[1,1]
	  
	 fold = top1[,1]; names(fold) <- rownames(top1)
	 fold <- convertEnsembl2Entrez(fold,Species)
	 
     keggSpecies <- as.character( keggSpeciesID[which(keggSpeciesID[,1] == Species),3] )
	 
     if(nchar( keggSpecies) <=2 ) return(blank) # not in KEGG
	cat("here5  ",keggSpecies, " ",Species," ",input$sigPathways)
	 # kegg pathway id
	pathID = keggPathwayID(selectedPathway, Species, "KEGG",selectOrg)

	cat("\n",fold[1:5],"\n",keggSpecies,"\n",pathID)
    if(is.null(pathID) ) return(blank) # kegg pathway id not found.	
	pv.out <- pathview(gene.data = fold, pathway.id = pathID, species = keggSpecies, kegg.native=TRUE)
	
	
	
	
	#######################################
	# testing for species not recognized 
	x=read.csv("exampleData/Wu_wet_vs_control - new species.csv");
	 #x = read.csv("exampleData/counts_test_data_3groups.csv") ;
	 x[,1] = toupper(x[,1]);  
	 colnames(x)[1]= "User_input"

	
	selectOrg = "BestMatch"; GO="KEGG"; 
	myrange = c(15,2000)

	converted = convertID(x[,1],selectOrg)
	mapping = converted$conversionTable
	
	
	
	  
 }
 
 debug_prep <- function () {
 
  inFile = "C:/Users/Xijin.Ge/Google Drive/research/Shiny/RNAseqer/expression1_no_duplicate.csv"
 # inFile = "C:/Users/Xijin.Ge/Google Drive/research/Shiny/RNAseqer/GSE52778_All_Sample_FPKM_Matrix.csv"
 lowFilter = 1; logStart = 1
    x = read.csv(inFile)
	x[,1] = toupper(x[,1])
    x = x[order(x[,1]),]
    x = x[!duplicated(x[,1]),]
    rownames(x)= x[,1]
    x = x[,-1]
	
	tem = apply(x,1,max)
	x = x[which(tem> lowFilter),] 
	
	x = log(x+abs( logStart),2)
	tem = apply(x,1,sd)
	x = x[order(-tem),]
 
    ###########Converted data
	convertedID = convertID(rownames(x ),selectOrg="BestMatch", selectGO = "GOBP" );#"gmax_eg_gene"

	mapping <- convertedID$conversionTable

	   rownames(x) = toupper(rownames(x))
	   x1 = merge(mapping[,1:2],x,  by.y = 'row.names', by.x = 'User_input')
	  	tem = apply(x1[,3:(dim(x1)[2]-2)],1,sd)
	   x1 = x1[order(x1[,2],-tem),]
	   x1 = x1[!duplicated(x1[,2]) ,]
	   rownames(x1) = x1[,2]
	   x1 = as.matrix(x1[,c(-1,-2)])

	   
     	tem = apply(x1,1,sd)
	x1 = x1[order(-tem),]
	x=x1
	 head(x)
	 
    #################################
	#testing Kmeans
	
	x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	#x = 100* x / apply(x,1,sum)  # this is causing problem??????
	#x = x - apply(x,1,mean)  # this is causing problem??????
	#colnames(x) = gsub("_.*","",colnames(x))
	set.seed(2)
	# determining number of clusters
	k=6
	
	cl = kmeans(x,k,iter.max = 50)
	#myheatmap(cl$centers)	

	hc <- hclust2(dist2(cl$centers-apply(cl$centers,1,mean) )  )# perform cluster for the reordering of samples
	tem = match(cl$cluster,hc$order) #  new order 

	x = x[order(tem),]

	bar = sort(tem)
	#myheatmap2(x, bar)
    # GO
	pp=0
	for( i in 1:k) {
	#incProgress(1/k, , detail = paste("Cluster",toupper(letters)[i]) )
	query = rownames(x)[which(bar == i)]
	convertedID = convertID(query,"BestMatch", selectGO = "GOBP" );#"gmax_eg_gene"
	tem = geneInfo(convertedID,"BestMatch") #input$selectOrg ) ;
	tem <- tem[which( tem$Set == "List"),] 
	
	
	#selectOrg = input$selectOrg
	selectOrg ="BestMatch"
	
	result = FindOverlap (convertedID,tem, "GOBP",selectOrg,1) 
	if( dim(result)[2] ==1) next;   # result could be NULL
	result$Genes = toupper(letters)[i] 
	if (pp==0 ) { results = result; pp = 1;} else  results = rbind(results,result)
	
	
	}
	results= results[,c(5,1,2,4)]
	colnames(results)= c("Cluster","FDR","Genes","GO BP Terms")
	minFDR = 0.05
	if(min(results$FDR) > minFDR ) results = as.matrix("No signficant enrichment found.") else
	results = results[which(results$FDR < minFDR),]
	 
	 
	 
	 
	 
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
    # if (is.null(input$file1) && input$goButton > 0 )   inFile = "expression1_no_duplicate.csv"
	if (is.null(input$file1) && input$goButton > 0 )   inFile = "GSE37704_sailfish_genecounts.csv"
	tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) # these are needed to make it responsive to changes
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	isolate({ 

	withProgress(message="Reading and pre-processing ", {
	if (is.null( input$dataFileFormat )) return(NULL)
	dataTypeWarning =0
    if (input$dataFileFormat == 2 ) {  # if FPKM
	    incProgress(1/3,"Pre-processing data")
		x <- read.csv(inFile)	# try CSV
        if(dim(x)[2] <= 2 )   # if less than 3 columns, try tab-deliminated
			x <- read.table(inFile, sep="\t",header=TRUE)		
		x[,1] <- toupper(x[,1])
		x = x[order(- apply(x[,2:dim(x)[2]],1,sd) ),]
		x <- x[!duplicated(x[,1]) ,]
		rownames(x) <- x[,1]
		x <- as.matrix(x[,c(-1)])
		
		if ( is.integer(x) ) dataTypeWarning = 1;  # Data appears to be read counts
		
		tem <- apply(x,1,max)
		x = x[which(tem> input$lowFilter),]  # max by row is at least 		
		x = x[which(apply(x,1, function(y) max(y)- min(y) ) > 0  ),]  # remove rows with all the same 
		# Takes log if log is selected OR kurtosis is big than 100
		#class(input$transform)
		mean.kurtosis = mean( apply(x,2, kurtosis) )
		if ( (input$transform == TRUE) | (mean.kurtosis > kurtosis.log ) ) x = log(x+abs( input$logStart),2)
		tem = apply(x,1,sd)
		x = x[order(-tem),]
		rawCounts=NULL
		
   } else {  # counts data
		incProgress(1/3, "Pre-processing counts data")
		
		tem = input$CountsDEGMethod; tem = input$countsTransform
		x = read.csv(inFile)
		if(dim(x)[2] <= 2 ) x <- read.table(inFile, sep="\t",header=TRUE)	
		# x = read.csv("exampleData/GSE87194.csv")
		# x = read.csv("exampleData/airway_GSE52778.csv")
		x[,1] = toupper(x[,1])
		x = x[order(x[,1],-rowSums(x[,2:dim(x)[2]]) ),]  # sort by gene id, then sum of counts
		x = x[!duplicated(x[,1]),]
		rownames(x)= x[,1]
		x = as.matrix( x[,-1] )
		
		if(!is.integer(x) ) dataTypeWarning = -1  # data not seems to be read counts
		x <- round(x,0) # enforce the read counts to be integers. Sailfish outputs has decimal points.
		x <- x[ which( apply(x,1,max) >= input$minCounts ) , ] # remove all zero counts
		
		# x = x[which(apply(x,1, function(y) max(y)- min(y) ) > 0  ),]  # remove rows with all the same 
		if(0){
			# remove genes with low expression by counts per million (CPM)
			dge <- DGEList(counts=x); dge <- calcNormFactors(dge)
			myCPM <- cpm(dge, prior.counts = 3 )
			x <- x[which(rowSums(  myCPM > input$minCounts)  > 1 ),]  # at least two samples above this level
			rm(dge); rm(myCPM)
		}
		rawCounts = x;
		
		# construct DESeqExpression Object
		# colData = cbind(colnames(x), as.character(detectGroups( colnames(x) )) )
		tem = rep("A",dim(x)[2]); tem[1] <- "B"   # making a fake design matrix to allow process, even when there is no replicates
		colData = cbind(colnames(x), tem )
		colnames(colData)  = c("sample", "groups")
		dds <- DESeqDataSetFromMatrix(countData = x,
                                  colData = colData,
                                  design = ~ groups)
		dds <- estimateSizeFactors(dds) # estimate size factor for use in normalization later for started log method
		incProgress(1/2,"transforming raw counts")
		# regularized log  or VST transformation
		if( input$CountsTransform == 3 && dim(counts(dds))[2] <= 10 )  # rlog is slow, only do it with 10 samples
		    { x <- rlog(dds, blind=FALSE); x <- assay(x) } else {
			if ( input$CountsTransform == 2 )    # vst is fast but aggressive
			   { x <- vst(dds, blind=FALSE); x <- assay(x)  } else{  # normalized by library sizes and add a constant.
				    x <- log2( counts(dds, normalized=TRUE) + input$countsLogStart )   # log(x+c)
				}
			}
			
		mean.kurtosis = mean( apply(x,2, kurtosis) )
								  
		
   }
       incProgress(1, "Done.")
     
    return( list(data = as.matrix(x), mean.kurtosis = mean.kurtosis, rawCounts = rawCounts, dataTypeWarning=dataTypeWarning) )
	 })
	})
  })	
  
output$text.transform <- renderText({
      if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	  inFile <- input$file1
	  
	  k.value =  readData()$mean.kurtosis	  
      tem = paste( "Mean Kurtosis =  ", round(k.value,2), ".\n",sep = "")
	  if( k.value > kurtosis.log) tem = paste(tem, " Detected extremely large numbers with kurtosis >", kurtosis.log,
   	  ". Log transformation is automatically applied, even user selects otherwise.", sep="") else if (k.value>kurtosis.warning)
	  {tem = paste(tem, " Detected  large numbers with kurtosis >",
   	  kurtosis.warning,". Log transformation is recommended.", sep="") }
	  
	  if(readData()$dataTypeWarning == 1 ) {
	      tem = paste (tem, " ------Warning!!! Data matrix contains all integers. It seems to be read counts!!! Please select appropriate data type in the previous page and reload file.")}
	  if(readData()$dataTypeWarning == -1 ) {
	       tem = paste (tem, "  ------Warning!!! Data matrix is not all integers. Data does not look like read counts data. Please select appropriate data type in the previous page and reload file.")}
	
	  tem
  
	})	

output$fileFormat <- renderText({
		" Upload CSV (comma separated value) files containing gene expression data derived from RNA-seq or microarray.
		For RNA-seq data, gene-level read count data is recommended, but FPKM or RPKM data also acceptable.
	  iDEP can convert most types of common gene ids to Ensembl gene IDs. Name samples in each column carefully as iDEP parses them to define sample groups. 
	  Replicates should be denotated by \"_Rep1,2,3...\" at the end. For example, 
	  Control_Rep1, Control_Rep2, Treatment_Rep1, Treatment_Rep2. For factorial experiment design, use underscore \"_\" to
	  separate factors such as genetic background and treatment. For example, 
      WT_control_Rep1, WT_control_Rep2, WT_Treatment_Rep1, WT_treatment_Rep2,
		Mu_control_Rep1, Mu_control_Rep2, Mu_Treatment_Rep1, Mu_treatment_Rep2.
		Also, avoid using hyphen \"-\" or dot \".\" in sample names."
	})	

	# this defines an reactive object that can be accessed from other rendering functions
converted <- reactive({  
	  if (is.null(input$file1) && input$goButton == 0)    return(NULL)
       tem = input$selectOrg
      isolate( {
   #   withProgress(message="Converting gene ids", {
	  convertID(rownames(readData()$data ),input$selectOrg, input$selectGO );
	         
	 # })
	  }) 

	} )

		# this defines an reactive object that can be accessed from other rendering functions
allGeneInfo <- reactive({  
	  if (is.null(input$file1) && input$goButton == 0)    return(NULL)
       tem = input$selectOrg
      isolate( {
      withProgress(message="Looking up gene annotation", {
	  geneInfo(converted(),input$selectOrg); 
			 })
	  })
	} )

convertedData <- reactive({  
	  if (is.null(input$file1) && input$goButton == 0)    return()
		tem = input$selectOrg;  tem = input$lowFilter ; tem = input$transform
      if( is.null(converted() ) ) return( readData()$data) # if id or species is not recognized use original data.
	  isolate( {  
	   mapping <- converted()$conversionTable
	   
	   x =readData()$data
	   rownames(x) = toupper(rownames(x))
	  
	  # any gene not recognized by the database is disregarded
	  # x1 = merge(mapping[,1:2],x,  by.y = 'row.names', by.x = 'User_input')

	   # the 3 lines keeps the unrecogized genes using original IDs
	   x1 = merge(mapping[,1:2],x,  by.y = 'row.names', by.x = 'User_input', all.y=TRUE)
	   ix = which(is.na(x1[,2]) )
	   x1[ix,2] = x1[ix,1] # original IDs used
	   #write.csv(x1,"tem.csv")
	   
	   tem = apply(x1[,3:(dim(x1)[2])],1,sd)
	   x1 = x1[order(x1[,2],-tem),]
	   x1 = x1[!duplicated(x1[,2]) ,]
	   rownames(x1) = x1[,2]
	   x1 = as.matrix(x1[,c(-1,-2)])
	   tem = apply(x1,1,sd)
	   x1 = x1[order(-tem),]  # sort again by SD
       return(x1)
	  }) 

	} )	
	
GeneSets <- reactive({  
	  if (is.null(input$file1) && input$goButton == 0)    return()

	  tem = input$selectOrg
	  tem = input$selectGO
	  tem =input$minSetSize; tem= input$maxSetSize
      isolate( {	  
	   readGeneSets( converted(), convertedData(), input$selectGO,input$selectOrg,c(input$minSetSize, input$maxSetSize)  )  }) 
	} )

gmtCategory1 <- reactive({  
	  if (is.null(input$file1) && input$goButton == 0)    return()
	  tem = input$selectOrg
      isolate( {
	  
	   gmtCategory(converted(), convertedData(), input$selectOrg) }) 

	} )
 
output$contents <- renderTable({
   inFile <- input$file1
	inFile <- inFile$datapath
    if (is.null(input$file1) && input$goButton == 0)   return(NULL)
#    if (is.null(input$file1) && input$goButton > 0 )   inFile = "expression1_no_duplicate.csv"
	if (is.null(input$file1) && input$goButton > 0 )   inFile = "GSE37704_sailfish_genecounts.csv"

	tem = input$selectOrg
	isolate({
	x <- read.csv(inFile)
	if(dim(x)[2] <= 2 ) x <- read.table(inFile, sep="\t",header=TRUE)	# not CSV
    #x <- readData()$data
     x[1:20,]
	})
  },include.rownames=FALSE)

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
    x <- readData()$data
	 par(mfrow=c(3,1))
	par(mar=c(8,4,4,4))
	myColors = rainbow(dim(x)[2])
	plot(density(x[,1]),col = myColors[1],xlab="Expresson values", ylab="Density", main= "Distribution of transformed data")
	for( i in 2:dim(x)[2] )
	lines(density(x[,i]),col=myColors[i] )
    legend("topright", colnames(x), lty=rep(1,dim(x)[2]), col=myColors )	
   # boxplot of first two samples, often technical replicates
	boxplot(x, las = 2, ylab="Transformed expression levels", main="Distribution of transformed data")
	plot(x[,1:2],xlab=colnames(x)[1],ylab=colnames(x)[2], main="Scatter plot of first two samples")
	
   }, height = 1500, width = 500)

output$correlationMatrix <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    # heatmap of correlation matrix
	x <- readData()$data
	maxGene <- apply(x,1,max)
	x <- x[which(maxGene > quantile(maxGene)[1] ) ,] # remove bottom 25% lowly expressed genes, which inflate the PPC
	
   melted_cormat <- melt(round(cor(x),2), na.rm = TRUE)
# melted_cormat <- melted_cormat[which(melted_cormat[,1] != melted_cormat[,2] ) , ]
	# Create a ggheatmap
	ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
	 geom_tile(color = "white")+
	 scale_fill_gradient2(low = "green", high = "red",  mid = "white", 
		space = "Lab",  limit = c(min(melted_cormat[,3]) ,max(melted_cormat[,3])), midpoint = median(melted_cormat[,3]),
		name="Pearson\nCorrelation"
	  ) +
	  theme_minimal()+ # minimal theme
	 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
		size = 12, hjust = 1))+
	 coord_fixed()
	# print(ggheatmap)
	 ggheatmap + 
	geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
	theme(
	  axis.title.x = element_blank(),
	  axis.title.y = element_blank(),
	  panel.grid.major = element_blank(),
	  panel.border = element_blank(),
	  panel.background = element_blank(),
	  axis.ticks = element_blank(),
	 legend.justification = c(1, 0),
	  legend.position = c(0.6, 0.7),
	 legend.direction = "horizontal")+
	 guides(fill = FALSE) + ggtitle("Pearson's Correlation Coefficient (all genes)")

 
  } , height = 600, width = 600)
   
output$downloadProcessedData <- downloadHandler(
     filename = function() {"processed_Data.csv"},
		content = function(file) {
      write.csv( 
	  #convertedData()
	merge(allGeneInfo()[,c('ensembl_gene_id','symbol')], round(convertedData(),4),by.x="ensembl_gene_id", by.y ="row.names", all.y=T )
	  , file, row.names=FALSE) 
	  
	  # write.csv(readData(), file)
	    }
 )
output$examineData <- DT::renderDataTable({
   inFile <- input$file1
	inFile <- inFile$datapath
    if (is.null(input$file1) && input$goButton == 0)   return(NULL)

	tem = input$selectOrg
	isolate({
	merge(allGeneInfo()[,c('ensembl_gene_id','symbol')], round(convertedData(),2),by.x="ensembl_gene_id", by.y ="row.names", all.y=T )
	})
  })

output$totalCounts <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    if (is.null(readData()$rawCounts))   return(NULL)
    par(mar=c(10,2,2,2))
    x <- readData()$rawCounts
	barplot( colSums(x)/1e6, col="green",las=3, main="Total read counts (millions)")

})
output$PCA <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    
	x <- convertedData();
     if(input$PCA_MDS ==1) {   #PCA
	 pca.object <- prcomp(t(x))
	 # par(mfrow=c(2,1))
	if(0){
     plot( pca.object$x[,1], pca.object$x[,2], pch = 1,cex = 2,col = detectGroups(colnames(x)),
	     xlim=c(min(pca.object$x[,1]),max(pca.object$x[,1])*1.5   ),
		xlab = "First principal component", ylab="Second Principal Component")
		text( pca.object$x[,1], pca.object$x[,2],  pos=4, labels =colnames(x), offset=.5, cex=.8)
		}
		
	pcaData = as.data.frame(pca.object$x[,1:2]); pcaData = cbind(pcaData,detectGroups(colnames(x)) )
	colnames(pcaData) = c("PC1", "PC2", "Type")
	percentVar=round(100*summary(pca.object)$importance[2,1:2],0)
	p=ggplot(pcaData, aes(PC1, PC2, color=Type, shape = Type)) + geom_point(size=5) 
	p=p+xlab(paste0("PC1: ",percentVar[1],"% variance")) 
	p=p+ylab(paste0("PC2: ",percentVar[2],"% variance")) 
	p=p+ggtitle("Principal component analysis (PCA)")+coord_fixed(ratio=1.0)+ 
     theme(plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio=1)
	   print(p)
	   }
	# variance chart
	# plot(pca.object,type="bar", xlab="Principal Components", main ="Variances explained")
	
	# pathways
	if(input$PCA_MDS ==2) {  
	withProgress(message="Running pathway analysis", {
	pca.object <- prcomp(t(x))
	pca = 100*pca.object$rotation 
	Npca = 5
	if (Npca > dim(pca)[2]) { Npca = dim(pca)[2] } else pca <-  pca[,1:Npca]
	#pca = pca[,1:5]
	if(is.null(GeneSets() ) ) return(NULL)  # no species recognized
	if(length(GeneSets() ) <= 1 ) return(NULL)
	#cat("\n\nGene Sets:",length( GeneSets()))
	pg = myPGSEA (pca,cl=GeneSets(),range=c(15,2000),p.value=TRUE, weighted=FALSE,nPermutation=1)
	incProgress(2/8)
	# correcting for multiple testing
	p.matrix = pg$p.result
	tem = p.adjust(as.numeric(p.matrix),"fdr")
	p.matrix = matrix(tem, nrow=dim(p.matrix)[1], ncol = dim(p.matrix)[2] )
	rownames(p.matrix) = rownames(pg$p.result); colnames(p.matrix) = colnames(pg$p.result)


	selected =c()
	for( i in 1:dim(p.matrix)[2]) {
	  tem = which( rank(p.matrix[,i],ties.method='first') <= 5)  # rank by P value
	 #tem = which( rank(pg$result[,i],ties.method='first') >= dim(p.matrix)[1]-3.1) # rank by mean
	 names(tem) = paste("PC",i," ", rownames(p.matrix)[tem], sep="" )
	 selected = c(selected, tem)
	}
	rowids = gsub(" .*","",names(selected))
	rowids = as.numeric( gsub("PC","",rowids) )
	pvals = p.matrix[ cbind(selected,rowids) ]
	a=sprintf("%-1.0e",pvals)
	tem = pg$result[selected,]
	rownames(tem) = paste(a,names(selected)); #colnames(tem)= paste("PC",colnames(tem),sep="")
	
	tem = tem[!duplicated(selected),] 
	incProgress(3/8)
	#tem = t(tem); tem = t( (tem - apply(tem,1,mean)) ) #/apply(tem,1,sd) )

	smcPlot(tem,scale =  c(-max(tem), max(tem)), show.grid = T, margins = c(3,1, 13, 23), col = .rwb,cex.lab=0.5, main="Pathways analysis on PCA")
	 } )
	 }
	 
	if(input$PCA_MDS ==3) {  # MDS
	 fit = cmdscale( dist2(t(x) ), eig=T, k=2)
	 
	# par(pin=c(5,5))
	if(0) {
	plot( fit$points[,1],fit$points[,2],pch = 1,cex = 2,col = detectGroups(colnames(x)),
	     xlim=c(min(fit$points[,1]),max(fit$points[,1])*1.5   ),
	  xlab = "First dimension", ylab="Second dimension"  )
	 text( fit$points[,1], fit$points[,2],  pos=4, labels =colnames(x), offset=.5, cex=1)
	}
	pcaData = as.data.frame(fit$points[,1:2]); pcaData = cbind(pcaData,detectGroups(colnames(x)) )
	colnames(pcaData) = c("x1", "x2", "Type")
	
	p=ggplot(pcaData, aes(x1, x2, color=Type, shape = Type)) + geom_point(size=5) 
	#p=p+xlab(paste0("PC1: ",percentVar[1],"% variance")) 
	#p=p+ylab(paste0("PC2: ",percentVar[2],"% variance")) 
	p=p+ggtitle("Multidimensional scaling (MDS)")+ coord_fixed(ratio=1.)+ 
     theme(plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio=1)
	   print(p)
	
	 }
	  
  }, height = 500, width = 500)

PCAdata <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    x <- readData()$data
	
	 result = prcomp(t(x))$x[,1:2]
	 fit = cmdscale( dist2(t(x) ), eig=T, k=2)
     result = cbind( result, fit$points[,1:2] )
	 colnames(result) = c("PCA.x","PCA.y","MDS.x", "MDS.y")
	 return( result)		  
  })
  
output$downloadPCAData <- downloadHandler(
     filename = function() {"PCA_and_MDS.csv"},
		content = function(file) {
          write.csv(PCAdata(), file) 
	    }
  )
  
output$KNN <- renderPlot({ # Kmeans clustering
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	withProgress(message="k-means clustering", {
    x <- convertedData()
	#x <- readData()
	#par(mfrow=c(1,2))
	n=input$nGenesKNN
	if(n>6000) n = 6000 # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
	#x1 <- x;
	#x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	x = 100* x[1:n,] / apply(x[1:n,],1,sum) 
	#x = x - apply(x,1,mean)  # this is causing problem??????
	#colnames(x) = gsub("_.*","",colnames(x))
	set.seed(2)
	k=input$nClusters
	
	cl = kmeans(x,k,iter.max = 50)
	#myheatmap(cl$centers)	
 
   incProgress(.3, detail = paste("Heatmap..."))
	hc <- hclust2(dist2(cl$centers-apply(cl$centers,1,mean) )  )# perform cluster for the reordering of samples
	tem = match(cl$cluster,hc$order) #  new order 
	x = x[order(tem),] ; 	bar = sort(tem)
	
	myheatmap2(x-apply(x,1,mean), bar,1000)
	
	incProgress(1, detail = paste("Done")) }) #progress 
  } , height = 500, width = 500)
output$KmeansNclusters <- renderPlot({ # Kmeans clustering
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	withProgress(message="k-means clustering", {
    x <- convertedData()
	#x <- readData()
	#par(mfrow=c(1,2))
	n=input$nGenesKNN
	if(n>6000) n = 6000 # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
	#x1 <- x;
	#x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	x = 100* x[1:n,] / apply(x[1:n,],1,sum)  # this is causing problem??????
	#x = x - apply(x,1,mean)  # this is causing problem??????
	#colnames(x) = gsub("_.*","",colnames(x))
	set.seed(2)
	# determining number of clusters
	incProgress(.3, detail = paste("Performing k-means..."))
	
	k = 30
	wss <- (nrow(x)-1)*sum(apply(x,2,var))
	  for (i in 2:k) wss[i] <- sum(kmeans(x,centers=i,iter.max = 30)$withinss)
	plot(1:k, wss, type="b", xlab="Number of Clusters (k)",
		 ylab="Within groups sum of squares")

	
	incProgress(1, detail = paste("Done")) }) #progress 
  } , height = 500, width = 500)
  
KNNdata <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	withProgress(message="Generating data", {
    x <- convertedData()
	#x <- readData()
	#par(mfrow=c(2,1))
	n=input$nGenesKNN
	if(n>6000) n = 6000 # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
	
	#x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	x = 100* x[1:n,] / apply(x[1:n,],1,sum) 
	#x = x - apply(x,1,mean)  # this is causing problem??????
	#colnames(x) = gsub("_.*","",colnames(x))
	set.seed(2)
	# determining number of clusters

	k=input$nClusters
	
	cl = kmeans(x,k,iter.max = 50)
	#myheatmap(cl$centers)	


	hc <- hclust2(dist2(cl$centers-apply(cl$centers,1,mean) )  )# perform cluster for the reordering of samples
	tem = match(cl$cluster,hc$order) #  new order 
	x = x[order(tem),] ; 	bar = sort(tem)
	
	
	})
	#myheatmap2(x, bar)
	Cluster = toupper(letters)[bar]
	x = cbind(Cluster,x)
	
		# add gene symbol
	ix = match( rownames(x), allGeneInfo()[,1])
	return( cbind(as.character( allGeneInfo()$symbol)[ix],x))
	
	 #progress 
  })
  
output$downloadDataKNN <- downloadHandler(
     filename = function() {"Kmeans.csv"},
		content = function(file) {
      write.csv(KNNdata(), file)
	    }
  )
output$KNN_GO <- renderTable({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectGO3
	if( is.null(input$selectGO3 ) ) return (NULL)
	if( input$selectGO3 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.") )#No matching species
   	withProgress(message="GO Enrichment", {
	x <- convertedData()
	#x <- readData()
	#par(mfrow=c(2,1))
	n=input$nGenesKNN
	# if(n>6000) n = 6000 # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
	
	x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	#x = 100* x / apply(x,1,sum)  # this is causing problem??????
	#x = x - apply(x,1,mean)  # this is causing problem??????
	#colnames(x) = gsub("_.*","",colnames(x))
	set.seed(2)
	# determining number of clusters
	k=input$nClusters
	
	cl = kmeans(x,k,iter.max = 50)
	#myheatmap(cl$centers)	

	hc <- hclust2(dist2(cl$centers-apply(cl$centers,1,mean) )  )# perform cluster for the reordering of samples
	tem = match(cl$cluster,hc$order) #  new order 

	x = x[order(tem),]

	bar = sort(tem)
	#myheatmap2(x, bar)
    # GO
	pp=0
	for( i in 1:k) {
	incProgress(1/k, , detail = paste("Cluster",toupper(letters)[i]) )
	query = rownames(x)[which(bar == i)]
	convertedID = convertID(query,input$selectOrg, selectGO = input$selectGO3 );#"gmax_eg_gene"
	tem = geneInfo(convertedID,input$selectOrg) #input$selectOrg ) ;
	tem <- tem[which( tem$Set == "List"),]
	tem <- tem[which( tem$gene_biotype == "protein_coding"),]
	selectOrg = input$selectOrg
	#selectOrg ="BestMatch"
	
	result = FindOverlap (convertedID,tem, input$selectGO3,selectOrg,1) 
	if( dim(result)[2] ==1) next;   # result could be NULL
	result$Genes = toupper(letters)[i] 
	if (pp==0 ) { results <- result; pp <- 1;} else  results <- rbind(results,result)
	}
	if( pp == 0) return( as.data.frame("No enrichment found.") )
	results= results[,c(5,1,2,4)]
	colnames(results)= c("Cluster","FDR","Genes","Pathways")
	minFDR = 0.01
	if(min(results$FDR) > minFDR ) results = as.data.frame("No signficant enrichment found.") else
	results = results[which(results$FDR < minFDR),]
	incProgress(1, detail = paste("Done")) 
	}) #progress
	colnames(results)[2] = "adj.Pval"
	results$Genes <- as.character(results$Genes)
	results$Cluster[which( duplicated(results$Cluster) ) ] <- ""
	results
  }, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)

output$KmeansPromoter <- renderTable({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectGO3; tem = input$radioPromoterKmeans; tem=input$nGenesKNN; tem=input$nClusters
	if( is.null(input$selectGO3 ) ) return (NULL)
	if( is.null(limma()$results) ) return(NULL)

	isolate({ 
   	withProgress(message="Promoter analysis", {
	x <- convertedData()
	#x <- readData()
	#par(mfrow=c(2,1))
	n=input$nGenesKNN
	# if(n>6000) n = 6000 # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
	
	x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	#x = 100* x / apply(x,1,sum)  # this is causing problem??????
	#x = x - apply(x,1,mean)  # this is causing problem??????
	#colnames(x) = gsub("_.*","",colnames(x))
	set.seed(2)
	# determining number of clusters
	k=input$nClusters
	cl = kmeans(x,k,iter.max = 50)

	hc <- hclust2(dist2(cl$centers-apply(cl$centers,1,mean) )  )# perform cluster for the reordering of samples
	tem = match(cl$cluster,hc$order) #  new order 
	x = x[order(tem),]
	bar = sort(tem)
	
	results1 <- NULL; result <- NULL 
	pp<- 0
	for( i in 1:k ) {
	incProgress(1/k, , detail = paste("Cluster",toupper(letters)[i]) )
	query = rownames(x)[which(bar == i)]
	
	convertedID = convertID(query,input$selectOrg, input$selectGO2 );#"gmax_eg_gene"
	result <- promoter( convertedID,input$selectOrg,input$radioPromoterKmeans )
	
	if( is.null(result)  ) next;   # result could be NULL
	if(  dim(result)[2] ==1) next;
	result$List = toupper(letters)[i]    
	if (pp==0 ) { results1 <- result; pp <- 1 } else  { results1 = rbind(results1,result) }
	}

	incProgress(1, detail = paste("Done")) 
	}) #progress
	
	if( is.null(results1) ) {as.data.frame("No significant motif enrichment found.")} else {
		results1[ duplicated (results1[,4] ),4 ] <- ""
	  results1[,c(4,1:3,5)]
	  }
	})
  }, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)
  
output$heatmap <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    x <- readData()$data   # x = read.csv("expression1.csv")
	withProgress(message="Reading and pre-processing ", {
	n=input$nGenes
	if(n>6000) n = 6000 # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
	# this will cutoff very large values, which could skew the color 
	x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	cutoff = median(unlist(x)) + 4*sd (unlist(x)) 
	x[x>cutoff] <- cutoff
	cutoff = median(unlist(x)) - 4*sd (unlist(x)) 
	x[x< cutoff] <- cutoff
	
    groups = detectGroups(colnames(x) )
	groups.colors = rainbow(length(unique(groups) ) )

	#http://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
	lmat = rbind(c(5,4),c(0,1),c(3,2))
	lwid = c(1.5,6)
	lhei = c(1,.2,8)

	if( n>110) 
	heatmap.2(x, distfun = dist2,hclustfun=hclust2,
	 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
	,key=T, symkey=F
	,ColSideColors=groups.colors[ as.factor(groups)]
	,labRow=""
	,margins=c(10,0)
	,srtCol=45
	,lmat = lmat, lwid = lwid, lhei = lhei
	)

	if( n<=110) 
	heatmap.2(x, distfun = dist2,hclustfun=hclust2,
	 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
	,key=T, symkey=F,
	#,labRow=labRow
	,ColSideColors=groups.colors[ groups]
	,margins=c(10,12)
	,cexRow=1
	,srtCol=45
	,lmat = lmat, lwid = lwid, lhei = lhei
	)
	incProgress(1,"Done")
	})
  }, height = 800, width = 400)  
    
heatmapData <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    x <- readData()$data
	geneSD = apply(x,1,sd)
	x = x[order(-geneSD),]
	
	n=input$nGenes
	if(n>6000) n = 6000 # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
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
	
	x1 <- x1[ rev( hy$rowInd),hy$colInd]
	# add gene symbol
	ix = match( rownames(x1), allGeneInfo()[,1])
	return( cbind(as.character( allGeneInfo()$symbol)[ix],x1) )

	
  })  
  
output$downloadData <- downloadHandler(
     filename = function() {"heatmap.csv"},
		content = function(file) {
			write.csv(heatmapData(), file)
	    }
  )
  
output$debug <- renderTable({
      if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	  tem = input$selectOrg; tem = input$lowFilter ; tem = input$transform
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

output$selectGO2 <- renderUI({
	  tem = input$selectOrg
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectGO2", label = NULL, # h6("Funtional Category"), 
                  choices = list("All available gene sets" = "All", "GO Biological Process" = "GOBP","GO Molecular Function" = "GOMF","GO Cellular Component" = "GOCC",
                                "KEGG metabolic pathways" = "KEGG"), selected = "GOBP")  }	 else { 
								
	  selectInput("selectGO2", label=NULL,choices=gmtCategory(converted(), convertedData(), input$selectOrg)
	     ,selected = "GOBP" )   } 
	})
	
output$selectGO3 <- renderUI({
	  tem = input$selectOrg
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectGO3", label = NULL, # h6("Funtional Category"), 
                  choices = list("All available gene sets" = "All", "GO Biological Process" = "GOBP","GO Molecular Function" = "GOMF","GO Cellular Component" = "GOCC",
                                "KEGG metabolic pathways" = "KEGG"), selected = "GOBP")  }	 else { 
								
	  selectInput("selectGO3", label=NULL,choices=gmtCategory(converted(), convertedData(), input$selectOrg)
	     ,selected = "GOBP" )   } 
	})

output$PGSEAplot <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO;		tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	
	if(is.null(input$selectGO ) ) return (NULL)
	if(input$selectGO == "ID not recognized!" ) return( NULL)
	isolate({ 
	withProgress(message="Running pathway analysis", {
	myrange = c(input$minSetSize, input$maxSetSize)
	genes = convertedData()
if (is.null(input$selectContrast1 ) ) return(NULL)
	incProgress(1/4,"Retrieing gene sets")
	gmt = GeneSets()
	incProgress(2/4,"Runing PGSEA.")
	#if(0){
	iz= match( detectGroups(colnames(genes)), unlist(strsplit( input$selectContrast1, "-"))	  )
    iz = which(!is.na(iz))
	if (grepl("Diff:",input$selectContrast1) == 1) iz=1:(dim(genes)[2]) 
	genes = genes[,iz]
	#}
	subtype = detectGroups(colnames(genes )) 
    if(length( GeneSets() )  == 0)  { plot.new(); text(0.5,0.5, "No gene sets!")} else {
	result = PGSEApathway(converted(),genes, input$selectOrg,input$selectGO,
	             GeneSets(),  myrange, input$pathwayPvalCutoff, input$nPathwayShow 	)
					 
	if( is.null(result$pg3) ) { plot.new(); text(0.5,1, "No significant pathway found!")} else 
	smcPlot(result$pg3,factor(subtype),scale = c(-max(result$pg3), max(result$pg3)), 
	show.grid = T, margins = c(3,1, 13, 23), col = .rwb,cex.lab=0.5, main="Pathway Analysis:PGSEA")
    }
	
	})
	})
    }, height = 800, width = 500)

output$PGSEAplotAllSamples <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	
	if(is.null(input$selectGO ) ) return (NULL)
	if(input$selectGO == "ID not recognized!" ) return( NULL)
	isolate({ 
	withProgress(message="Running pathway analysis", {
	myrange = c(input$minSetSize, input$maxSetSize)
	genes = convertedData()
if (is.null(input$selectContrast1 ) ) return(NULL)
	incProgress(1/4,"Retrieing gene sets")
	gmt = GeneSets()
	incProgress(2/8, "Runing PGSEA")
	subtype = detectGroups(colnames(genes )) 
    if(length( GeneSets() )  == 0)  { plot.new(); text(0,1, "No gene sets!")} else {
	result = PGSEApathway(converted(),genes, input$selectOrg,input$selectGO,
	             GeneSets(),  myrange, input$pathwayPvalCutoff, input$nPathwayShow 	)
					 
	if( is.null(result$pg3) ) { plot.new(); text(0.5,1, "No significant pathway found!")} else 
	smcPlot(result$pg3,factor(subtype),scale = c(-max(result$pg3), max(result$pg3)), 
	show.grid = T, margins = c(3,1, 13, 23), col = .rwb,cex.lab=0.5, main="Pathway Analysis:PGSEA")
    }
	
	})
	})
    }, height = 800, width = 500)

output$gagePathway <- renderTable({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO;	tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	
	isolate({ 
	gagePathwayData()
	})
  },digits=0,align="l",include.rownames=FALSE,striped=TRUE,bordered = TRUE, width = "auto",hover=T)
 
gagePathwayData <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)

	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO
	tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	
	if(is.null(input$selectGO ) ) return (NULL)
	if(input$selectGO == "ID not recognized!" ) return( as.data.frame("Gene ID not recognized." ))
	isolate({ 
	withProgress(message="Running pathway analysis using GAGE", {
	if (is.null(input$selectContrast1 ) ) return(NULL)
	myrange = c(input$minSetSize, input$maxSetSize)
	noSig = as.data.frame("No significant pathway found.")
	if( length(limma()$topGenes) == 0 ) return(noSig)
	if(length( limma()$comparisons)  ==1 )  
    { top1=limma()$topGenes[[1]]  
	} else {
	  top = limma()$topGenes
	  ix = match(input$selectContrast1, names(top))
	  if( is.na(ix)) return (noSig)
	  top1 <- top[[ix]]; 
	  }
	  if(dim(top1)[1] == 0 ) return (noSig)
	  colnames(top1)= c("Fold","FDR")
	  incProgress(1/4,"Retrieing gene sets")
	  gmt = GeneSets() 
      if(length( GeneSets() )  == 0)  { return(as.data.frame("No gene set found!"))}
	  #converted = convertID(rownames(top1),input$selectOrg)
	  #
	   #gmt = readGeneSets(converted, top1, input$selectGO, input$selectOrg, myrange )
     # cat("Sets",length(gmt))
	 incProgress(2/4,"Runing GAGE")
	 fold = top1[,1]; names(fold) <- rownames(top1)
	 
	 if(input$absoluteFold) fold <- abs(fold)
	 paths <- gage(fold, gsets = gmt, ref = NULL, samp = NULL)

	  paths <-  rbind(paths$greater,paths$less)
	 # write.csv(paths,"tem.csv")	  
	 # cat( dim(paths) )
	  if(dim(paths)[1] < 1 | dim(paths)[2]< 6 ) return( noSig )
	  top1 <- paths[,c('stat.mean','set.size','q.val')]
	  colnames(top1)= c("statistic","Genes","adj.Pval")
	  top1 <- top1[order(top1[,3]) ,]  
	  if ( length( which( top1[,3] <=  input$pathwayPvalCutoff   ) ) == 0 )
	    return( noSig)
	  top1 <- top1[which(top1[,3] <=  input$pathwayPvalCutoff ) ,,drop=FALSE]
	  if(dim(top1)[1] > input$nPathwayShow ) 
	     top1 <- top1[1:input$nPathwayShow, ,drop=FALSE]
		 
		top1 <- as.data.frame(top1)
		top1 <- cbind(rep( input$selectContrast1, dim(top1)[1]),row.names(top1), top1); 
		top1$statistic <- as.character( round(as.numeric(top1$statistic),4)); 
		top1$adj.Pval <- sprintf("%-2.1e",as.numeric(top1$adj.Pval) )
		top1[,2] <- as.character(top1[,2]);top1[,1] <- as.character(top1[,1])
		colnames(top1)[1] <- "Direction"
		if(input$pathwayMethod == 1 ) p.m <- "GAGE"
		else if(input$pathwayMethod == 2 ) p.m <- "PGSEA"
		else if(input$pathwayMethod == 3 ) p.m <- "GSEA"
		else if(input$pathwayMethod == 4 ) p.m <- "PGSEA_All"
		else if(input$pathwayMethod == 5 ) p.m <- "ReactomePA"
		colnames(top1)[2] <- paste(p.m," analysis:", gsub("-"," vs ",input$selectContrast1 ) )
		top1[ which( top1[,3] >0),1 ] <- gsub("-"," > ",input$selectContrast1 )
		top1[ which( top1[,3] <0),1 ] <- gsub("-"," < ",input$selectContrast1 )
		top1 <- top1[order( top1[,1], -abs(as.numeric( top1[,3]) ) ) ,]
		top1[ duplicated (top1[,1] ),1 ] <- ""
		#write.csv(top1,"tem.csv")
	  return( top1)
	}) })
  })

output$fgseaPathway <- renderTable({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO; tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	
	isolate({ 
	fgseaPathwayData()
	})
  },digits=0,align="l",include.rownames=FALSE,striped=TRUE,bordered = TRUE, width = "auto",hover=T)
 
fgseaPathwayData <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO; tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold
	if(is.null(input$selectGO ) ) return (NULL)
	if(input$selectGO == "ID not recognized!" ) return( as.data.frame("Gene ID not recognized." ))
	isolate({ 
	withProgress(message="Running pathway analysis using GAGE", {
	if (is.null(input$selectContrast1 ) ) return(NULL)
	myrange = c(input$minSetSize, input$maxSetSize)
	noSig = as.data.frame("No significant pathway found.")
	if( length(limma()$topGenes) == 0 ) return(noSig)
	if(length( limma()$comparisons)  ==1 )  
    { top1=limma()$topGenes[[1]]  
	} else {
	  top = limma()$topGenes
	  ix = match(input$selectContrast1, names(top))
	  if( is.na(ix)) return (noSig)
	  top1 <- top[[ix]]; 
	  }
	  if(dim(top1)[1] == 0 ) return (noSig)
	  colnames(top1)= c("Fold","FDR")
	  incProgress(1/4,"Retrieing gene sets")
	  gmt = GeneSets() 
     if(length( GeneSets() )  == 0)  { return(as.data.frame("No gene set found!"))}

	  #converted = convertID(rownames(top1),input$selectOrg)
	  #
	   #gmt = readGeneSets(converted, top1, input$selectGO, input$selectOrg, myrange )
     # cat("Sets",length(gmt))
	 incProgress(2/4,"Runing GSEA using the fgsea package.")
	 fold = top1[,1]; names(fold) <- rownames(top1)
	 if(input$absoluteFold) fold <- abs(fold) # use absolute value of fold change, disregard direction
	 paths <- fgsea(pathways = gmt, 
                  stats = fold,
                  minSize=input$minSetSize,
                  maxSize=input$maxSetSize,
                  nperm=5000)
	 # paths <-  rbind(paths$greater,paths$less)
	  if(dim(paths)[1] < 1  ) return( noSig )
	       paths <- as.data.frame(paths)
       paths <- paths[order(-abs( paths[,5])) ,]  # sort by NES
	  # paths <- paths[order( paths[,3]) ,]  # sort by FDR
	  top1 <- paths[,c(1,5,7,3)]
	  # rownames(top1) <- paths[,1] #paste(1:dim(paths)[1],": ",paths[,1],sep="" )
	  colnames(top1)= c("Pathway", "NES","Genes","adj.Pval")
	  
	  if ( length( which( top1[,4] <=  input$pathwayPvalCutoff   ) ) == 0 )
	    return( noSig)
	  top1 <- top1[which(top1[,4] <=  input$pathwayPvalCutoff ) , ,drop=FALSE]
	  if(dim(top1)[1] > input$nPathwayShow ) 
	     top1 <- top1[1:input$nPathwayShow,,drop=FALSE]
	 #top1 <- cbind(row.names(top1), top1); colnames(top1)[1] <-input$selectContrast1 	
		top1 <- as.data.frame(top1)
		top1 <- cbind(rep( input$selectContrast1, dim(top1)[1]), top1); 
		top1[,4] = as.character( round(as.numeric(top1[,4]),4)); 
		top1$adj.Pval <- sprintf("%-2.1e",as.numeric(top1$adj.Pval) )
		top1[,1] <- as.character(top1[,1])
		colnames(top1)[1] <- "Direction"
		if(input$pathwayMethod == 1 ) p.m <- "GAGE"
		else if(input$pathwayMethod == 2 ) p.m <- "PGSEA"
		else if(input$pathwayMethod == 3 ) p.m <- "GSEA"
		else if(input$pathwayMethod == 4 ) p.m <- "PGSEA_All"
		else if(input$pathwayMethod == 5 ) p.m <- "ReactomePA"
		colnames(top1)[2] <- paste(p.m," analysis:", gsub("-"," vs ",input$selectContrast1 ) )
		top1[ which( as.numeric( top1[,3]) >0),1 ] <- gsub("-"," > ",input$selectContrast1 )
		top1[ which( as.numeric( top1[,3]) <0),1 ] <- gsub("-"," < ",input$selectContrast1 )
		top1 <- top1[order( top1[,1], -abs(as.numeric( top1[,3]) ) ) ,]
		top1[ duplicated (top1[,1] ),1 ] <- ""	 
	    top1[,3] = as.character( round(as.numeric(top1[,3]),4));
	 return( top1)
	}) })
  })

ReactomePAPathwayData <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO; tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	
	if(is.null(input$selectGO ) ) return (NULL)
	if(input$selectGO == "ID not recognized!" ) return( as.data.frame("Gene ID not recognized." ))
	isolate({ 
	withProgress(message="Running pathway analysis using ReactomePA", {
	if (is.null(input$selectContrast1 ) ) return(NULL)
	
	ensemblSpecies <- c("hsapiens_gene_ensembl","rnorvegicus_gene_ensembl", "mmusculus_gene_ensembl",
		               "celegans_gene_ensembl","scerevisiae_gene_ensembl", "drerio_gene_ensembl", "dmelanogaster_gene_ensembl")
    ReactomePASpecies= c("human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly" )
    # cbind(ensemblSpecies,ReactomePASpecies)  # double check mapping
    
	
	myrange = c(input$minSetSize, input$maxSetSize)
	noSig = as.data.frame("No significant pathway found.")
	if( length(limma()$topGenes) == 0 ) return(noSig)
	if(length( limma()$comparisons)  ==1 )  
    { top1=limma()$topGenes[[1]]  
	} else {
	  top = limma()$topGenes
	  ix = match(input$selectContrast1, names(top))
	  if( is.na(ix)) return (noSig)
	  top1 <- top[[ix]]; 
	  }
	  if(dim(top1)[1] == 0 ) return (noSig)
	  colnames(top1)= c("Fold","FDR")
	  incProgress(1/4,"Retrieing gene sets")
	  # gmt = GeneSets() 

	  #converted = convertID(rownames(top1),input$selectOrg)
	  #
	   #gmt = readGeneSets(converted, top1, input$selectGO, input$selectOrg, myrange )
     # cat("Sets",length(gmt))
	 
	 fold = top1[,1]; names(fold) <- rownames(top1)
	 if(input$absoluteFold) fold <- abs(fold) # use absolute value of fold change, disregard direction

      incProgress(2/4,"Convert to Entrez Gene ID")
      Species <- converted()$species[1,1]
      ix <- match( Species,ensemblSpecies )	
	  if(is.na(ix) ) return(as.data.frame("Species not coverted by ReactomePA package!"))
	  
	  fold <- convertEnsembl2Entrez (fold, Species)  
		 fold <- sort(fold,decreasing =T)
	  incProgress(3/4,"Runing enrichment analysis using fgsea")
	  paths <- gsePathway(fold, nPerm=5000, organism = ReactomePASpecies[ix],
                minGSSize= input$minSetSize, 
				maxGSSize= input$maxSetSize,
				pvalueCutoff=0.5,
                pAdjustMethod="BH", verbose=FALSE)
      paths <- as.data.frame(paths)
	  
	  if(is.null(paths) ) return( noSig)
	  if(dim(paths)[1] ==0 ) return( noSig)

	 # paths <-  rbind(paths$greater,paths$less)
	  if(dim(paths)[1] < 1  ) return( noSig )
	       paths <- as.data.frame(paths)
       paths <- paths[order(-abs( paths[,5])) ,]  # sort by NES
	  # paths <- paths[order( paths[,3]) ,]  # sort by FDR
	  top1 <- paths[,c(2,5,3,7)]
	  # rownames(top1) <- paths[,1] #paste(1:dim(paths)[1],": ",paths[,1],sep="" )
	  colnames(top1)= c("Pathway", "NES","Genes","adj.Pval")
	  
	  if ( length( which( top1[,4] <=  input$pathwayPvalCutoff   ) ) == 0 )
	    return( noSig)
	  top1 <- top1[which(top1[,4] <=  input$pathwayPvalCutoff ) ,,drop=FALSE]
	  if(dim(top1)[1] > input$nPathwayShow ) 
	     top1 <- top1[1:input$nPathwayShow,,drop=FALSE]
	 #top1 <- cbind(row.names(top1), top1); colnames(top1)[1] <-input$selectContrast1 	
		top1 <- as.data.frame(top1)
		top1 <- cbind(rep( input$selectContrast1, dim(top1)[1]), top1); 
		top1[,4] = as.character( round(as.numeric(top1[,4]),4)); 
		top1$adj.Pval <- sprintf("%-2.1e",as.numeric(top1$adj.Pval) )
		top1[,1] <- as.character(top1[,1])
		colnames(top1)[1] <- "Direction"
		if(input$pathwayMethod == 1 ) p.m <- "GAGE"
		else if(input$pathwayMethod == 2 ) p.m <- "PGSEA"
		else if(input$pathwayMethod == 3 ) p.m <- "GSEA"
		else if(input$pathwayMethod == 4 ) p.m <- "PGSEA_All"
		else if(input$pathwayMethod == 5 ) p.m <- "ReactomePA"
		colnames(top1)[2] <- paste(p.m," analysis:", gsub("-"," vs ",input$selectContrast1 ) )
		top1[ which( as.numeric( top1[,3]) >0),1 ] <- gsub("-"," > ",input$selectContrast1 )
		top1[ which( as.numeric( top1[,3]) <0),1 ] <- gsub("-"," < ",input$selectContrast1 )
		top1 <- top1[order( top1[,1], -abs(as.numeric( top1[,3]) ) ) ,]
		top1[ duplicated (top1[,1] ),1 ] <- ""	 
	  top1[,3] = as.character( round(as.numeric(top1[,3]),4));
	 return( top1)
	}) })
  })

output$ReactomePAPathway <- renderTable({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO; tem = input$selectContrast1
    tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	
	isolate({ 
	ReactomePAPathwayData()
	})
  },digits=0,align="l",include.rownames=FALSE,striped=TRUE,bordered = TRUE, width = "auto",hover=T)
  
PGSEAplot.data <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg
	
	isolate({ 
	withProgress(message="Running pathway analysis", {
	myrange = c(input$minSetSize, input$maxSetSize)
	genes = convertedData()
if (is.null(input$selectContrast1 ) ) return(NULL)
	gmt = GeneSets()
	incProgress(2/8)
	#if(0){
	iz= match( detectGroups(colnames(genes)), unlist(strsplit( input$selectContrast1, "-"))	  )
    iz = which(!is.na(iz))
	if (grepl("Diff:",input$selectContrast1) == 1) iz=1:(dim(genes)[2]) 
	genes = genes[,iz]
	#}
	subtype = detectGroups(colnames(genes )) 
    if(length( GeneSets() )  == 0)  { plot.new(); text(0,1, "No gene sets!")} else {
	result = PGSEApathway(converted(),genes, input$selectOrg,input$selectGO,
	             GeneSets(),  myrange, input$pathwayPvalCutoff, input$nPathwayShow 	)
					 
	if( is.null(result$pg3) ) { return(as.matrix("No significant pathway!"))} else 
	   return( result$pg3)
    }
	
	})
	})
    })

output$download.PGSEAplot.data <- downloadHandler(
     filename = function() {"PGSEA_pathway_anova.csv"},
		content = function(file) {
			write.csv(PGSEAplot.data(), file)
	    }
  )
  
output$listSigPathways <- renderUI({
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	
      if (is.null(input$file1)&& input$goButton == 0 | is.null(gagePathwayData()))
	  
       { selectInput("sigPathways", label = NULL, # h6("Funtional Category"), 
                  choices = list("All" = "All"), selected = "All")  }	 else { 
		choices <- "All"  # default, sometimes these methods returns "No significant pathway found"
		if( input$pathwayMethod == 1) { if(!is.null(gagePathwayData())) if(dim(gagePathwayData())[2] >1) choices <- gagePathwayData()[,2] }
		else if( input$pathwayMethod == 3) { if(!is.null(fgseaPathwayData())) if(dim(fgseaPathwayData())[2] >1) choices <- fgseaPathwayData()[,2] }
		else if( input$pathwayMethod == 5) { if(!is.null(ReactomePAPathwayData())) if(dim(ReactomePAPathwayData())[2] >1) choices <- ReactomePAPathwayData()[,2] }
		selectInput("sigPathways", label="Show gene expression of a selected pathway:",choices=choices)
	        } 
	})
	
selectedPathwayData <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg
	tem = input$sigPathways; 
	if(is.null(gagePathwayData() ) ) return(NULL)
	if(is.null( input$sigPathways))  return (NULL) 
	#cat(input$sigPathways)
    isolate({ 
	withProgress (message="Retrieving genes",{ 
    if(input$sigPathways == "All") return (NULL) 
	#as.data.frame(input$sigPathways )
	genes <- pathwayGenes(input$sigPathways,converted(), input$selectGO1,input$selectOrg )
	incProgress(1/2,"Merging data")
	
    x <-  convertedData()[which(rownames(convertedData()) %in% genes) ,]
	ix = match( rownames(x), allGeneInfo()[,1])
	rownames(x) <- paste(rownames(x),":", as.character( allGeneInfo()$symbol)[ix]) 
	return( x )
	
     }) })
})

output$downloadSelectedPathwayData <- downloadHandler(
    # filename = function() {"Selected_Pathway_detail.csv"},
	  filename = function() {paste(input$selectContrast1,"(",input$sigPathways,")",".csv",sep="")},
		content = function(file) {
			write.csv(selectedPathwayData(), file)
	    }
  )
  
output$selectedPathwayData1 <- renderTable({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg
	tem = input$sigPathways; 
	if(is.null(gagePathwayData() ) ) return(NULL)
	if(is.null( input$sigPathways))  return (NULL) 
	# cat(input$sigPathways)

    isolate({ 
	selectedPathwayData()
    }) 
})

output$selectedPathwayHeatmap <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg
	tem = input$sigPathways; 
	if(is.null(gagePathwayData() ) ) return(NULL)
	if(is.null( input$sigPathways))  return (NULL) 
	if( is.null(selectedPathwayData()) ) return(NULL)
	# cat(input$sigPathways)

    isolate({ 
	x = selectedPathwayData()
	if(dim(x)[1]<=2 | dim(x)[2]<=2 ) return(NULL)
	groups = detectGroups(colnames(x))
	withProgress(message="Generating heatmap",{ 
	# this will cutoff very large values, which could skew the color 
	x=as.matrix(x)-apply(as.matrix(x),1,mean)
	cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
	x[x>cutoff] <- cutoff
	cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
	x[x< cutoff] <- cutoff
	
	lmat = rbind(c(5,4),c(0,1),c(3,2))
	lwid = c(1.5,6)
	lhei = c(1,.2,8)

	if( dim(x)[1]>100) 
	heatmap.2(x, distfun = dist2,hclustfun=hclust2,
	 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
	,key=T, symkey=F
	,ColSideColors=mycolors[ groups]
	,labRow=""
	,margins=c(10,0)
	,srtCol=45
	,lmat = lmat, lwid = lwid, lhei = lhei
	#,main ="Title"
	)

	if( dim(x)[1]<=100) 
	heatmap.2(x, distfun = dist2,hclustfun=hclust2,
	 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
	,key=T, symkey=F,
	#,labRow=labRow
		,ColSideColors=mycolors[ groups]
	,margins=c(10,12)
	,cexRow=1.5
	,srtCol=45
	,lmat = lmat, lwid = lwid, lhei = lhei
	#,main ="Title"
	)
	incProgress(1,"Done")
    }) 
	})
}, height = 800, width = 400)

output$KeggImage <- renderImage({
   # First generate a blank image. Otherse return(NULL) gives us errors.
    outfile <- tempfile(fileext='.png')
    png(outfile, width=400, height=300)
    frame()
	dev.off()
    blank <- list(src = outfile,
         contentType = 'image/png',
         width = 400,
         height = 300,
         alt = " ")
    if (is.null(input$file1)&& input$goButton == 0)   return(blank)

	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO
	tem = input$selectContrast
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	
	tem = input$sigPathways; 
	if(is.null( input$selectGO ) ) return(blank)
	if(input$selectGO != "KEGG") return(blank)
	if(is.null(gagePathwayData() ) ) return(blank)
	if(is.null( input$sigPathways))  return (blank) 
	# if( is.null(selectedPathwayData()) ) return(blank)
	
	isolate({ 
	withProgress(message="Rendering KEGG pathway plot", {
	
	if (is.null(input$selectContrast1 ) ) return(blank)
	
	if(input$sigPathways == "All") return (blank) 
	
	if( length(limma()$topGenes) == 0 ) return(blank)

	# get fold change
	if(length( limma()$comparisons)  ==1 )  
    { top1=limma()$topGenes[[1]]  
	} else {
	  top = limma()$topGenes
	  ix = match(input$selectContrast1, names(top))
	  if( is.na(ix)) return (blank)
	  top1 <- top[[ix]]; 
	  }
	  if(dim(top1)[1] == 0 ) return (blank)
	 # cat("here5")
	  colnames(top1)= c("Fold","FDR")
      Species <- converted()$species[1,1]
	  
	 fold = top1[,1]; names(fold) <- rownames(top1)
	 fold <- convertEnsembl2Entrez(fold,Species)
	 
     keggSpecies <- as.character( keggSpeciesID[which(keggSpeciesID[,1] == Species),3] )
	 
     if(nchar( keggSpecies) <=2 ) return(blank) # not in KEGG

	 # kegg pathway id
	incProgress(1/2, "Download pathway graph from KEGG.")
	pathID = keggPathwayID(input$sigPathways, Species, "KEGG",input$selectOrg)
	#cat("here5  ",keggSpecies, " ",Species," ",input$sigPathways, "pathID:",pathID,"End", fold[1:5],names(fold)[1:5],"\n")
	#cat("pathway:",is.na(input$sigPathways))
	#cat("\n",fold[1:5],"\n",keggSpecies,"\n",pathID)
    if(is.null(pathID) ) return(blank) # kegg pathway id not found.
	randomString <- gsub(".*file","",tempfile()) 
	pv.out <- mypathview(gene.data = fold, pathway.id = pathID, kegg.dir = "temp",  out.suffix = randomString, species = keggSpecies, kegg.native=TRUE)
	#
	outfile <- paste( "temp/",pathID,".",randomString,".png",sep="")
   
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
       width = "200%",
        height = "200%",
         alt = "KEGG pathway image.")
		}) 
	})
  }, deleteFile = TRUE)
	
limma <- reactive({  
  if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$CountsDEGMethod; tem = input$countsLogStart
	isolate({ 
	withProgress(message="Identifying Differentially expressed genes", {
	if(input$dataFileFormat == 1 ) {  # if count data
		 if(input$CountsDEGMethod == 3 )   # if DESeq2 method
		  # rawCounts = read.csv("exampleData/airway_GSE52778.csv", row.names=1)
		 # res =DEG.DESeq2(rawCounts, .05, 2) 
		  # res1 =DEG.limma(rawCounts, .1, 1.5,rawCounts, 2,3) 
			return( DEG.DESeq2(readData()$rawCounts,input$limmaPval, input$limmaFC)  )
		if(input$CountsDEGMethod < 3 )    # voom or limma-trend
			return( DEG.limma(convertedData(), input$limmaPval, input$limmaFC,readData()$rawCounts, input$CountsDEGMethod,priorCounts=input$countsLogStart,input$dataFileFormat) )
	} else { # normalized data
	 return( DEG.limma(convertedData(), input$limmaPval, input$limmaFC,readData()$rawCounts, input$CountsDEGMethod,priorCounts=input$countsLogStart,input$dataFileFormat) )
	}
	
	
	})
	})
	})	

output$text.limma <- renderText({
      if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
 	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	  limma()$Exp.type
  
	})	
  
output$vennPlot <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	isolate({ 
	
	results = limma()$results
	if(dim(results)[2] >5) results <- results[,1:5]
	vennDiagram(results,circle.col=rainbow(5))
 
	
	})
    }, height = 600, width = 600)

output$listComparisons <- renderUI({
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectContrast", label = NULL, # h6("Funtional Category"), 
                  choices = list("All" = "All"), selected = "All")  }	 else { 
	  selectInput("selectContrast", label="Select a comparison to examine. A-B means A vs. B.",choices=limma()$comparisons
	     )   } 
	})

output$listComparisonsPathway <- renderUI({
	tem = input$selectOrg

      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectContrast1", label = NULL, # h6("Funtional Category"), 
                  choices = list("All" = "All"), selected = "All")  }	 else { 
	  selectInput("selectContrast1", label="Select a comparison to analyze:",choices=limma()$comparisons
	     )   } 
	})
		
DEG.data <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; 
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts; tem = input$lowFilter; tem=input$transform; tem = input$logStart
	isolate({ 
	  genes = limma()$results
	  genes = as.data.frame( genes[which( rowSums(genes) != 0 ),] )
	  colnames(genes) = colnames( limma()$results )
	  genes = merge(genes,convertedData(), by='row.names')
	  colnames(genes)[1] = "1: upregulation, -1: downregulation"
	  	# add gene symbol
	ix = match( genes[,1], allGeneInfo()[,1])
	genes <- cbind(as.character( allGeneInfo()$symbol)[ix],genes) 
	colnames(genes)[1] = "Symbol"
	genes <- genes[,c(2,1,3:dim(genes)[2]) ]
	return(genes)
	})
    })

output$selectedHeatmap <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast;
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts; tem = input$lowFilter; tem=input$transform; tem = input$logStart
	withProgress(message="Generating heatmap", {
	if( is.null(input$selectContrast)) return(NULL)
	if( is.null(limma()$results)) return(NULL)
	isolate({ 
	  genes = limma()$results
	  if( is.null(genes) ) return(NULL)
	  ix = match(input$selectContrast, colnames(genes)) 
	  if(is.null(ix)) reutn(NULL); if(is.na(ix)) return (NULL)
	  if( sum(abs(genes[,ix] )  ) <= 1 ) return(NULL) # no significant genes for this comparison
	  query = rownames(genes)[which(genes[,ix] != 0)]
      iy = match(query, rownames(convertedData()  ) )
	 iz= match( detectGroups(colnames(convertedData())), unlist(strsplit( input$selectContrast, "-"))	  )
     iz = which(!is.na(iz))
	 if (grepl("Diff:",input$selectContrast) == 1) iz=1:(dim(convertedData())[2]) # if it is factor design use all samples
	 genes = convertedData()[iy,iz]
	 groups = detectGroups(colnames(genes) )
	 N1 = sum(groups == groups[1] ) # number of samples in type above
	 fc = rowMeans(genes[,1:N1])-rowMeans(genes[,(N1+1):length(groups)])
	 genes = genes[order(fc),]
	 fc = sort(fc)
	 bar = (fc>0 )+1
	 incProgress(1/2 )
 	 myheatmap2( genes,bar,200 )
	 incProgress(1, detail = paste("Done")) 
	
	 })
	})
	
    }, height = 400, width = 500)

selectedHeatmap.data <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast;
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts; tem = input$lowFilter; tem=input$transform; tem = input$logStart
	if( is.null(input$selectContrast)) return(NULL)
	if( is.null(limma()$results) ) return(NULL)
	isolate({ 
	  genes <- limma()$results
	  if( is.null(genes) ) return(NULL)
	  ix = match(input$selectContrast, colnames(genes)) 
	  if(is.null(ix)) return(NULL)
	  if(is.na(ix)) return(NULL)
	  if( sum(abs(genes[,ix] )  ) <= 1 ) return(NULL) # no significant genes for this comparison
	  if(dim(genes)[2] < ix ) return(NULL)
	  query = rownames(genes)[which(genes[,ix] != 0)]
	  if(length(query) == 0) return(NULL)
      iy = match(query, rownames(convertedData()  ) )
	  
     iz= match( detectGroups(colnames(convertedData())), unlist(strsplit( input$selectContrast, "-"))	  )
     iz = which(!is.na(iz))
	 if (grepl("Diff:",input$selectContrast) == 1) iz=1:(dim(convertedData())[2]) # if it is factor design use all samples
	 genes = convertedData()[iy,iz]
	 groups = detectGroups(colnames(genes) )
	 N1 = sum(groups == groups[1] ) # number of samples in type above
	 fc = rowMeans(genes[,1:N1])-rowMeans(genes[,(N1+1):length(groups)])
	 genes = genes[order(fc),]
 	 return(genes)
	})
    })
output$download.selectedHeatmap.data <- downloadHandler(
     filename = function() {paste("Diff_genes_heatmap_",input$selectContrast,".csv",sep="")},
		content = function(file) {
			write.csv(selectedHeatmap.data(), file)
	    }
  )
	
output$download.DEG.data <- downloadHandler(
     filename = function() {"Diff_expression_all_comparisons.csv"},
		content = function(file) {
			write.csv(DEG.data(), file,row.names=FALSE)
	    }
  )
  
output$geneList <- renderTable({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts; tem = input$lowFilter; tem=input$transform; tem = input$logStart
	noSig = as.data.frame("No significant genes find!")
	if( is.null(input$selectContrast) ) return(NULL)
	if( is.null( limma()$comparisons ) ) return(NULL) # if no significant genes found
	if( length(limma()$topGenes) == 0 ) return(noSig)
	if(length( limma()$comparisons)  ==1 )  
    { top1=limma()$topGenes[[1]]  
	} else {
	  top = limma()$topGenes
	  ix = match(input$selectContrast, names(top))
	  if( is.na(ix)) return (noSig)
	  top1 <- top[[ix]]; 
	  }
	  if(dim(top1)[1] == 0 ) return (noSig)
	  colnames(top1)= c("Fold","FDR")
	  #top1 = merge(top1,convertedData(), by='row.names')
	  #colnames(top1)[1] = "Genes"
	  top1 = top1[order(-abs(top1$Fold)) ,]
	  if ( length( which( top1$FDR <=  input$limmaPval  &  abs(top1$Fold)  >= log2(input$limmaFC) ) ) == 0 )
	    return( noSig)
	  top1 <- top1[which(top1$FDR <=  input$limmaPval ) ,]
	  top1 <- top1[which(abs(top1$Fold)  >= log2( input$limmaFC)) ,]
	  top1$Top_Genes <- rownames(top1)
	  top1 <- top1[,c(3,1,2)]
	  
	  # if new species
	  if( input$selectGO2 == "ID not recognized!" ) return (top1); 
	  
	  #convertedID = convertID(top1[,1],input$selectOrg, "GOBP" );#"gmax_eg_gene"
	  # tem <- geneInfo(convertedID,input$selectOrg) #input$selectOrg ) ;
	 #  tem <- geneInfo(converted(),input$selectOrg)
	  top1 <- merge(top1, allGeneInfo(), by.x ="Top_Genes", by.y="ensembl_gene_id",all.x=T )
	  
      if ( sum( is.na(top1$band)) == dim(top1)[1] ) top1$chr = top1$chromosome_name else
	  	top1$chr = paste( top1$chromosome_name, top1$band,sep="")
	  
	  top1 <- top1[,c('Top_Genes','Fold','FDR','symbol','chr','gene_biotype')]
	  
	 
	#  ix = match(top1[,1], tem$ensembl_gene_id)
	 # if( sum(is.na( tem$Symbol[ix]) ) != length(ix) ) 
	  # { top1 <- cbind(top1, tem$Symbol[ix]); colnames(top1)[4]= "Symbol" }
	  top1 = top1[order(-abs(as.numeric( top1$Fold))) ,]
	  top1$FDR <- sprintf("%-2.1e",top1$FDR )
	  colnames(top1) <- c("Ensembl ID", "log2 Fold Change", "Adj.Pval", "Symbol","Chr","Type")
	  if ( sum( is.na(top1$Symbol)) == dim(top1)[1] ) top1 <- top1[,-4] 

	  top1
	
  },digits=2,align="l",include.rownames=F,striped=TRUE,bordered = TRUE, width = "auto",hover=T)
  #, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T,include.rownames=TRUE)
output$volcanoPlot <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts; tem = input$lowFilter; tem=input$transform; tem = input$logStart
	if( is.null(input$selectContrast) ) return(NULL)
	if( is.null( limma()$comparisons ) ) return(NULL) # if no significant genes found
	if( length(limma()$topGenes) == 0 ) return(NULL)
	isolate({ 
	withProgress(message="Generating volcano plot",{ 
	
	if(length( limma()$comparisons)  ==1 )  
    { top1=limma()$topGenes[[1]]  
	} else {
	  top = limma()$topGenes
	  ix = match(input$selectContrast, names(top))
	  if( is.na(ix)) return (NULL)
	  top1 <- top[[ix]]; 
	  }
	  incProgress(1/2)
	  if(dim(top1)[1] == 0 ) return (NULL)
	  colnames(top1)= c("Fold","FDR")
	 top1 <- as.data.frame(top1) # convert to data frame
     top1 <- top1[which(!(is.na(top1$Fold)|is.na(top1$FDR)    )),] # remove NA's 
	 top1$upOrDown <- 1
	 #write.csv(top1,"tem.csv")
	 top1$upOrDown[ which(top1$FDR <=  input$limmaPval& top1$Fold  >= log2( input$limmaFC)) ]  <- 2
	 top1$upOrDown[ which(top1$FDR <=  input$limmaPval & top1$Fold  <= -log2( input$limmaFC)) ]  <- 3
	 plot(top1$Fold,-log10(top1$FDR),col = c("grey30", "red","blue")[top1$upOrDown],
	 pch =16, cex = .3, xlab= "log2 fold change", ylab = "- log10 (FDR)")    
     legend("bottomright",c("Upregulated","Downregulated"),fill = c("red","blue") )
	 
	}) })

  },height=450, width=500)

output$scatterPlot <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts; tem = input$lowFilter; tem=input$transform; tem = input$logStart
	if( is.null(input$selectContrast) ) return(NULL)
	if( is.null( limma()$comparisons ) ) return(NULL) # if no significant genes found
	if( length(limma()$topGenes) == 0 ) return(NULL)
	isolate({ 
	withProgress(message="Generating scatter plot with all genes",{ 
	if(length( limma()$comparisons)  ==1 )  
    { top1=limma()$topGenes[[1]]  
	} else {
	  top = limma()$topGenes
	  ix = match(input$selectContrast, names(top))
	  if( is.na(ix)) return (NULL)
	  top1 <- top[[ix]]; 
	  }
	  if(dim(top1)[1] == 0 ) return (NULL)
	  colnames(top1)= c("Fold","FDR")
	 top1 <- as.data.frame(top1) # convert to data frame
     top1 <- top1[which(!(is.na(top1$Fold)|is.na(top1$FDR)    )),] # remove NA's 
	 top1$upOrDown <- 1
	 #write.csv(top1,"tem.csv")
	 top1$upOrDown[ which(top1$FDR <=  input$limmaPval& top1$Fold  >= log2( input$limmaFC)) ]  <- 2
	 top1$upOrDown[ which(top1$FDR <=  input$limmaPval & top1$Fold  <= -log2( input$limmaFC)) ]  <- 3

     incProgress(1/2)
	 if (grepl("Diff:",input$selectContrast) == 1) return(NULL) # iz=1:(dim(convertedData())[2]) # if it is factor design use all samples
  	 # average expression
     samples = unlist(strsplit( input$selectContrast, "-"))
     iz= match( detectGroups(colnames(convertedData())), samples[1]	  )
     iz = which(!is.na(iz))
	 genes <- convertedData()[,iz]
	 genes <- as.data.frame(genes)
	 genes$average1 <- apply( genes,1,mean)

     iz= match( detectGroups(colnames(convertedData())), samples[2]	  )
     iz = which(!is.na(iz))
	 # genes <- cbind( genes,convertedData()[,iz] )
	 genes$average2 <- apply(convertedData()[,iz] ,1,mean)	 
	genes <- merge(genes,top1,by="row.names")
 # write.csv(genes, "tem.csv")
	plot(genes$average2,genes$average1,col = c("grey45","red","blue")[genes$upOrDown],
	 pch =16, cex = .3, xlab= paste("Average expression in", samples[2] ), ylab = paste("Average expression in", samples[1] ))    
     legend("bottomright",c("Upregulated","Downregulated"),fill = c("red","blue") )

		})

	})
	 
  },height=450, width=500)
  
output$geneListGO <- renderTable({
    
	if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	if( is.null(input$selectContrast)) return(NULL)
	if( is.null( input$selectGO2) ) return (NULL)
	if( input$selectGO2 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.")) #No matching species

	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO2
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts; tem = input$lowFilter; tem=input$transform; tem = input$logStart
	if( is.null(limma()$results) ) return(NULL)
	if( is.null(selectedHeatmap.data()) ) return(NULL) # this has to be outside of isolate() !!!
	NoSig = as.data.frame("No significant enrichment found.")
	isolate({
 	withProgress(message="GO Enrichment", {
		
 
	genes <- selectedHeatmap.data()
	if(is.null(genes) ) return(NULL) 
	if(dim(genes)[1] <= minGenesEnrichment ) return(NoSig) # if has only few genes
	
	groups = detectGroups(colnames(genes) )
	N1 = sum(groups == groups[1] ) # number of samples in type above
	if( dim(genes)[2] < N1) return(NoSig)
	fc = rowMeans(genes[,1:N1])-rowMeans(genes[,(N1+1):length(groups)])
	
	# GO
	results1 <- NULL; result <- NULL
	pp <- 0
	for( i in c(1,-1) ) {
		incProgress(1/2 )
		
		if( length(which(fc*i<0)) <= minGenesEnrichment) next; 
		query = rownames(genes)[which(fc*i<0)]
		convertedID = convertID(query,input$selectOrg, input$selectGO2 );#"gmax_eg_gene"
		if( length(convertedID) <= minGenesEnrichment) next; 
		tem = geneInfo(convertedID,input$selectOrg) #input$selectOrg ) ;
		tem <- tem[which( tem$Set == "List"),]
		tem <- tem[which( tem$gene_biotype == "protein_coding"),]
		result <- FindOverlap (convertedID,tem, input$selectGO2,input$selectOrg,1) 
	# cat("\n",i," ",class(result),"\ncontent",result)
	if( dim(result)[2] ==1) next;   # result could be NULL
		if(i == 1) result$Genes = "A"  else result$Genes = "B"
		if (pp==0 ) { results1 <- result; pp = 1;} else  results1 = rbind(results1,result)
	}

	if ( pp == 0 ) return (NoSig)
	if ( is.null( results1) ) return (NoSig)
	if( dim(results1)[2] == 1 ) return(NoSig)  # Returns a data frame: "No significant results found!"
	
	results1= results1[,c(5,1,2,4)]
	colnames(results1)= c("List","FDR","Genes","GO terms or pathways")
	minFDR = 0.01
	if(min(results1$FDR) > minFDR ) results1 = as.data.frame("No signficant enrichment found.") else
	results1 = results1[which(results1$FDR < minFDR),]
	
	incProgress(1, detail = paste("Done")) 
	
	if(dim(results1)[2] != 4) return(NoSig)
	colnames(results1)= c("Direction","adj.Pval","Genes","Pathways")
	
	results1$adj.Pval <- sprintf("%-2.1e",as.numeric(results1$adj.Pval) )
	results1[,1] <- as.character(results1[,1])
	tem <- results1[,1]
	results1[ which(tem == "A"),1] <- gsub("-"," > ",input$selectContrast )
	results1[ which(tem == "B"),1] <- gsub("-"," < ",input$selectContrast )
	results1[ duplicated (results1[,1] ),1 ] <- ""
	
	results1
	 })#progress
	}) #isolate
  }, digits = 0,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)
 #   output$selectedHeatmap <- renderPlot({       hist(rnorm(100))    })

output$DEG.Promoter <- renderTable({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	if( is.null(input$selectContrast)) return(NULL)

	tem = input$selectOrg; tem = input$radio.promoter
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO2
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts; tem = input$lowFilter; tem=input$transform; tem = input$logStart
	if( is.null(limma()$results) ) return(NULL)
	if( is.null(selectedHeatmap.data()  ) ) return(NULL)
	#if( !is.data.frame(selectedHeatmap.data()  ) ) return(NULL)
	
  isolate({ 
   	withProgress(message="Promoter analysis", {
	genes <- selectedHeatmap.data()
	if(is.null(genes)) return(NULL)
	# if( !is.data.frame(genes)) return (NULL)
	# cat("\nHere",paste(genes[1,],collapse=" ") )
	if( dim(genes)[1] < minGenesEnrichment ) return (NULL) # skip if less than 5 genes total
	groups = detectGroups(colnames(genes) )
	N1 = sum(groups == groups[1] ) # number of samples in type above
	fc = rowMeans(genes[,1:N1])-rowMeans(genes[,(N1+1):length(groups)])
	
	# GO
	pp <- 0; results1 <- NULL; result <- NULL
	for( i in c(1, -1) ) {
	incProgress(1/2 )	
	query = rownames(genes)[which(fc*i<0)]
	if(length(query) < minGenesEnrichment) next; 
	convertedID = convertID(query,input$selectOrg, input$selectGO2 );#"gmax_eg_gene"
	if(length(convertedID) < minGenesEnrichment) next; 
	result <- promoter( convertedID,input$selectOrg,input$radio.promoter )
	
	if( is.null(result)  ) next;   # result could be NULL
	if(  dim(result)[2] ==1) next;
	
	if(i == 1) result$List ="A"  else result$List ="B" 
	if (pp==0 ) { results1 <- result; pp <- 1 } else  { results1 = rbind(results1,result) }
	}

	incProgress(1, detail = paste("Done")) 
	}) #progress
	
	if( is.null(results1)) {as.data.frame("No significant motif enrichment found.")} else {
	results1 <- results1[,c(4,1:3,5)]
	tem <- results1[,1]
	results1[ which( tem == "A"),1] <- gsub("-"," > ",input$selectContrast )
	results1[ which( tem == "B"),1] <- gsub("-"," < ",input$selectContrast )
	results1[ duplicated (results1[,1] ),1 ] <- ""
	results1
	   }
  })
  }, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)
	 
})
