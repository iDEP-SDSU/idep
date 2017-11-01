## PLAN dplyr should be used for all filter and mutate process

iDEPversion = "iDEP 0.43"
################################################################
# R packages
################################################################
# R packages, installed by:
# install.packages(c("shiny","RSQLite","gplots","ggplot2","e1071","shinyAce","shinyBS","reshape2","DT","plotly" ) )
# auto install
Rlibs = c("shiny","RSQLite","gplots","ggplot2","e1071","shinyAce","shinyBS","reshape2","DT","plotly","statmod" )
notInstalled = setdiff(Rlibs, rownames(installed.packages()))
if(length(notInstalled)>0)
	install.packages(notInstalled)

# To test these packages, start an R session and paste these lines below.
#library(shiny)   	# for Shiny interface
library(RSQLite)	# for database connection
library(gplots)		# for hierarchical clustering
library(ggplot2)	# graphics
library(e1071) 		# computing kurtosis
library(reshape2) 	# for melt correlation matrix in heatmap
library(DT) 		# for renderDataTable
library(plotly) 	# for interactive heatmap

# Bioconductor packages
#source("https://bioconductor.org/biocLite.R")
#biocLite(c( "limma", "DESeq2","edgeR","gage", "PGSEA", "fgsea", "ReactomePA", "pathview","PREDA","PREDAsampledata","sfsmisc","lokern","multtest" ))
# annotation packages needed by pathview; will be installed automatically if runing on Windows
#biocLite( c( "org.Ag.eg.db","org.At.tair.db","org.Bt.eg.db","org.Ce.eg.db","org.Cf.eg.db","org.Dm.eg.db","org.Dr.eg.db","org.EcK12.eg.db","org.EcSakai.eg.db","org.Gg.eg.db","org.Hs.eg.db","org.Hs.ipi.db","org.Mm.eg.db","org.Mmu.eg.db","org.Pf.plasmo.db","org.Pt.eg.db","org.Rn.eg.db","org.Sc.sgd.db","org.Sco.eg.db","org.Ss.eg.db","org.Tgondii.eg.db","org.Xl.eg.db")  )
# auto install 
biocLibs = c( "limma", "DESeq2","edgeR","gage", "PGSEA", "fgsea", "ReactomePA", "pathview","PREDA","PREDAsampledata","sfsmisc","lokern","multtest","hgu133plus2.db" )
notInstalled = setdiff(biocLibs, rownames(installed.packages()))
if(length(notInstalled)>0) { 
	source("https://bioconductor.org/biocLite.R")
	biocLite(notInstalled)
}

#library(DESeq2) # count data analysis
#library(edgeR) # count data D.E.

#library(PGSEA) # pathway 

library(sfsmisc)
library(lokern)
library(multtest)


################################################################
# Global variables
################################################################

pdf(NULL) # this prevents error Cannot open file 'Rplots.pdf'
Min_overlap <- 2
minSetSize = 3; 
mappingCoverage = 0.60 # 60% percent genes has to be mapped for confident mapping
mappingEdge = 0.5  # Top species has 50% more genes mapped
PvalGeneInfo = 0.05; minGenes = 10 # min number of genes for ploting
kurtosis.log = 50  # log transform is enforced when kurtosis is big
kurtosis.warning = 10 # log transformation recommnded 
minGenesEnrichment = 2 # perform GO or promoter analysis only if more than this many genes
PREDA_Permutations =1000
set.seed(2) # seed for random number generator
mycolors = sort(rainbow(20))[c(1,20,10,11,2,19,3,12,4,13,5,14,6,15,7,16,8,17,9,18)] # 20 colors for kNN clusters
#Each row of this matrix represents a color scheme;
heatColors = rbind(      greenred(75),     bluered(75),     colorpanel(75,"green","black","magenta"),colorpanel(75,"blue","yellow","red") )
rownames(heatColors) = c("Green-Black-Red","Blue-White-Red","Green-Black-Magenta","Blue-Yellow-Red")
colorChoices = setNames(1:dim(heatColors)[1],rownames(heatColors)) # for pull down menu

################################################################
#   Input files
################################################################

# this need to be removed. Also replace go to go for folder
#  setwd("C:/Users/Xijin.Ge/Google Drive/research/Shiny/RNAseqer")

# relative path to data files
datapath = "../../data/"   # production server
#datapath = "../../../go/"  # windows
#datapath = "../go/" # digital ocean

sqlite  <- dbDriver("SQLite")
convert <- dbConnect(sqlite,paste0(datapath,"convertIDs.db"),flags=SQLITE_RO)  #read only mode
keggSpeciesID = read.csv(paste0(datapath,"data_go/KEGG_Species_ID.csv"))
# List of GMT files in /gmt sub folder
gmtFiles = list.files(path = paste0(datapath,"pathwayDB"),pattern=".*\\.db")
gmtFiles = paste(datapath,"pathwayDB/",gmtFiles,sep="")
geneInfoFiles = list.files(path = paste0(datapath,"geneInfo"),pattern=".*GeneInfo\\.csv")
geneInfoFiles = paste(datapath,"geneInfo/",geneInfoFiles,sep="")
motifFiles = list.files(path = paste0(datapath,"motif"),pattern=".*\\.db")
motifFiles = paste(datapath,"motif/",motifFiles,sep="")
demoDataFile = paste0(datapath,"data_go/GSE37704_sailfish_genecounts.csv") #"expression1_no_duplicate.csv"

################################################################
#   Utility functions
################################################################

# Functions for hierarchical clustering
hclust2 <- function(x, method="average", ...)  # average linkage
  hclust(x, method=method, ...)
hclust.ward.D <- function(x, method="ward.D", ...)  # average linkage
  hclust(x, method=method, ...)
hclust.ward.D2 <- function(x, method="ward.D2", ...)  # average linkage
  hclust(x, method=method, ...)
hclust.single <- function(x, method="single", ...)  # average linkage
  hclust(x, method=method, ...)
hclust.mcquitty <- function(x, method="mcquitty", ...)  # average linkage
  hclust(x, method=method, ...)
hclust.median <- function(x, method="median", ...)  # average linkage
  hclust(x, method=method, ...)
hclust.centroid <- function(x, method="centroid", ...)  # average linkage
  hclust(x, method=method, ...)
  
hclustFuns <- list( averge = hclust2, complete=hclust, single=hclust.single,
					median=hclust.median, centroid=hclust.centroid, mcquitty=hclust.mcquitty)
hclustChoices = setNames(1:length(hclustFuns),names(hclustFuns)) # for pull down menu

dist2 <- function(x, ...)   # distance function = 1-PCC (Pearson's correlation coefficient)
  as.dist(1-cor(t(x), method="pearson"))
  
dist3 <- function(x, ...)   # distance function = 1-abs(PCC) (Pearson's correlation coefficient)
  as.dist(1-abs(cor(t(x), method="pearson")))   
  
# List of distance functions 
distFuns <- list(Correlation=dist2, Euclidean=dist,AbsolutePCC=dist3)
distChoices = setNames(1:length(distFuns),names(distFuns)) # for pull down menu

# Given a set of numbers, find the difference between 2nd largest and 2nd smallest
# 2,3,5,6,1   --> 5-2 = 3
geneChange <- function(x){
	n = length(x)
	if( n<4) return( max(x)-min(x)  ) else 
	return(sort(x)[n-1] - sort(x)[2]   )
}

dynamicRange <- function( x ) {
	y = sort(x)
	if(length(x)>=4)  k =2 else k =1;
	return( y[length(x)-k+1] - y[k]) 
}  

# Define sample groups based on column names
 detectGroups <- function (x){  # x are col names
	tem <- gsub("[0-9]*$","",x) # Remove all numbers from end
	#tem = gsub("_Rep|_rep|_REP","",tem)
	tem <- gsub("_$","",tem); # remove "_" from end
	tem <- gsub("_Rep$","",tem); # remove "_Rep" from end
	tem <- gsub("_rep$","",tem); # remove "_rep" from end
	tem <- gsub("_REP$","",tem)  # remove "_REP" from end
 	return( tem )
 }

# heatmap with color bar define gene groups
myheatmap2 <- function (x,bar=NULL,n=-1,mycolor=1,clusterNames=NULL ) {
	# number of genes to show
	ngenes = as.character( table(bar))
	if(length(bar) >n && n != -1) {ix = sort( sample(1:length(bar),n) ); bar = bar[ix]; x = x[ix,]  }

	# this will cutoff very large values, which could skew the color 
	x=as.matrix(x)-apply(x,1,mean)
	cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
	x[x>cutoff] <- cutoff
	cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
	x[x< cutoff] <- cutoff
	#colnames(x)= detectGroups(colnames(x))
	if(is.null(bar)) # no side colors
		heatmap.2(x,  Rowv =F,Colv=F, dendrogram ="none",
		col=heatColors[as.integer(mycolor),], density.info="none", trace="none", scale="none", keysize=.3
		,key=F, labRow = F,
		#,RowSideColors = mycolors[bar]
		,margins = c(8, 24)
		,srtCol=45
		) else
		heatmap.2(x,  Rowv =F,Colv=F, dendrogram ="none",
		col=heatColors[as.integer(mycolor),], density.info="none", trace="none", scale="none", keysize=.3
		,key=F, labRow = F,
		,RowSideColors = mycolors[bar]
		,margins = c(8, 24)
		,srtCol=45
		)
	if(!is.null(bar)) { 

		legend.text = paste("Cluster ", toupper(letters)[unique(bar)], " (N=", ngenes,")", sep="") 
		if( !is.null( clusterNames ) && length(clusterNames)>= length( unique(bar) ) )  
			legend.text = paste(clusterNames[ 1:length( unique(bar) )  ], " (N=", ngenes,")", sep="") 
		
		par(lend = 1)           # square line ends for the color legend
		legend("topright",      # location of the legend on the heatmap plot
		legend = legend.text, # category labels
		col = mycolors,  # color key
		lty= 1,             # line style
		lwd = 10 )           # line width
		}
}



# adding sample legends to heatmap; this is for the main heatmap
# https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 6), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

# Clean up gene sets. Remove spaces and other control characters from gene names  
cleanGeneSet <- function (x){
  # remove duplicate; upper case; remove special characters
  x <- unique( toupper( gsub("\n| ","",x) ) )
  x <- x[which( nchar(x)>1) ]  # genes should have at least two characters
  return(x)
}

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

################################################################
#   Main functions
################################################################

# Given a gene set, finds significant overlaps with a gene set database  object 
findOverlapGMT <- function ( query, geneSet, minFDR=.2 ,minSize=2,maxSize=10000 ){ 
	#geneSets <- readGMT("exampleData/MousePath_TF_gmt.gmt")
	#query <-  geneSets[['TF_MM_FRIARD_C-REL']] 
	#query <- query[1:60]
	total_elements = 30000
	Min_overlap <- 1
	maxTerms =10 # max number of enriched terms
	noSig <- as.data.frame("No significant enrichment found!")
	query <- cleanGeneSet(query)   # convert to upper case, unique()

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

# Runs pathway analysis using PGSEA; this is copied and revised from PGSEA package
myPGSEA  <- function (exprs, cl, range = c(25, 500), ref = NULL, center = TRUE, 
    p.value = 0.005, weighted = TRUE, nPermutation=100, enforceRange = TRUE, ...) {
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


# [ConvertDB Class START]
# prepare species list

# Create a list for Select Input options
orgInfo <- dbGetQuery(convert, paste("select distinct * from orgInfo " ))
orgInfo <- orgInfo[order(orgInfo$name),]
speciesChoice <- setNames(as.list( orgInfo$id ), orgInfo$name2 )
# add a defult element to list    # new element name       value
speciesChoice <- append( setNames( "NEW","**NEW SPECIES**"), speciesChoice  )
speciesChoice <- append( setNames( "BestMatch","Best matching species"), speciesChoice  )

# move one element to the 2nd place
move2 <- function(i) c(speciesChoice[1:2],speciesChoice[i],speciesChoice[-c(1,i)])
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
# [ConvertDB Class END]

# This function convert gene set names
# x="GOBP_mmu_mgi_GO:0000183_chromatin_silencing_at_rDNA"
# chromatin silencing at rDNA
proper <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))

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
	querySet <- cleanGeneSet( unlist( strsplit( toupper(query),'\t| |\n|\\,')))
	# querySet is ensgene data for example, ENSG00000198888, ENSG00000198763, ENSG00000198804
	querySTMT <- paste( "select distinct id,ens,species from mapping where id IN ('", paste(querySet,collapse="', '"),"')",sep="")
	result <- dbGetQuery(convert, querySTMT)
	if( dim(result)[1] == 0  ) return(NULL)
	if(selectOrg == speciesChoice[[1]]) {
		comb = paste( result$species,result$idType)
		sortedCounts = sort(table(comb),decreasing=T)
		recognized =names(sortedCounts[1])
		result <- result[which(comb == recognized),]
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
	result <- result[which(!duplicated(result[,1]) ),] # remove duplicates in user ID
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

convertEnsembl2KEGG <- function (query,Species) {  # not working
	querySet <- cleanGeneSet( unlist( strsplit( toupper(names( query)),'\t| |\n|\\,' )  ) )
	speciesID <- orgInfo$id[ which(orgInfo$ensembl_dataset == Species)]  # note uses species Identifying
	# idType 6 for entrez gene ID
	result <- dbGetQuery( convert,
						paste( " select  id,ens,species from mapping where ens IN ('", paste(querySet,collapse="', '"),
								"') AND  idType ='",idType_KEGG,"'",sep="") )	# slow
							
	if( dim(result)[1] == 0  ) return(NULL)
	result <- subset(result, species==speciesID, select = -species)

	ix = match(result$ens,names(query)  )

	tem <- query[ix];  names(tem) = result$id
	return(tem)  
}

# retrieve detailed info on genes
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

gmtCategory <- function (converted, convertedData, selectOrg,gmtFile) {
	if(selectOrg == "NEW" && !is.null(gmtFile) )
		return( list(Custom_GeneSet ="Custom" ) )
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
	library(PGSEA)
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

if(0){ # for testing LIMMA
	x = read.csv("C:/Users/Xijin.Ge/Google Drive/research/Shiny/RNAseqer/doc/Hoxa1-1/GSE50813_reduced.csv")
	rownames(x) = x[,1]
	x = x[,-1]
	maxP_limma=.1; minFC_limma=2; rawCounts=NULL; countsDEGMethods=2;priorCounts=4; dataFormat=2;
} 

# Differential expression using LIMMA 
DEG.limma <- function (x, maxP_limma=.1, minFC_limma=2, rawCounts,countsDEGMethods,priorCounts, dataFormat, selectedComparisons=NULL, sampleInfo = NULL,modelFactors=NULL, blockFactor = NULL){
	library(limma) # Differential expression
	library(statmod)
	
	# many different situations: 1. just use sample names 2. just one factor  3. two factors no interaction
	# 4. two factors with interaction   5. block factor 
	
	
	topGenes = list();  limmaTrend = FALSE
	if( dataFormat == 2) {   # if normalized data
		eset = new("ExpressionSet", exprs=as.matrix(x)) } else { # counts data
			if (countsDEGMethods == 1 ) { # limma-trend method selected for counts data
				#dge <- DGEList(counts=rawCounts);
				#dge <- calcNormFactors(dge, method = "TMM")
				#eset <- cpm(dge, log=TRUE, prior.count=priorCounts)
				eset = new("ExpressionSet", exprs=as.matrix(x))  # use transformed data for limma 
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
		comparisons <-  paste(g[2],"-",g[1],sep="")  # "Mutant-WT"
		
		# no sample file, but user selected comparisons using column names
		if( is.null(modelFactors) & length( selectedComparisons) >0  ) 	
			comparisons = selectedComparisons
		
		design <- model.matrix(~0+groups)
		colnames(design) <- g
		
		if( !is.null(rawCounts) && countsDEGMethods == 2) {  # voom
			dge <- DGEList(counts=rawCounts);
			dge <- calcNormFactors(dge, method = "TMM")  # normalization
			v <- voom(dge, design); fit <- lmFit(v, design) } else 
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
		Exp.type = "2 sample groups."
	}
	
	if(length(g) > 2 ) { # more than two sample groups
	
		design <- model.matrix(~ 0+factor(groups))
		colnames(design) <- gsub(".*)","",colnames(design))
		
		if( !is.null(rawCounts) && countsDEGMethods == 2) {  # voom
			v <- voom(rawCounts, design); 
			fit <- lmFit(v, design) 
		} else 
			fit <- lmFit(eset, design)
		
		fit <- eBayes(fit, trend=limmaTrend)
		
		comparisons = ""
		for( i in 1:(length(g)-1) )
			for (j in (i+1):length(g)) 
			comparisons = c(comparisons,paste(g[j],"-",g[i],sep="" ) )
		comparisons <- comparisons[-1]

		# no sample file, but user selected comparisons using column names
		if( is.null(modelFactors) & length( selectedComparisons) >0  ) 	
			comparisons = selectedComparisons
		
		contrast1 <- makeContrasts(contrasts=comparisons[1], levels=design)
		if(length(comparisons)>1 )
			for( kk in 2:length(comparisons) )
				contrast1<-  cbind(contrast1,makeContrasts(contrasts=comparisons[kk], levels=design)   )
		Exp.type = paste(length(g)," sample groups detected.")	 
		
		# if factorial design 2x2, 2x3, 3x5 etc.
			# all samples must be something like WT_control_rep1
		if ( sum(sapply(strsplit(g,"_"),length) == 2 ) == length(g) ) {
		
			#comparisons
			comparisons = ""
			for( i in 1:(length(g)-1) )
				for (j in (i+1):length(g)) 
				if( strsplit(g[i],"_")[[1]][1] == strsplit(g[j],"_")[[1]][1]| strsplit(g[i],"_")[[1]][2] == strsplit(g[j],"_")[[1]][2]) # only compare WT_control vs. WT_treatment
					comparisons = c(comparisons,paste(g[j],"-",g[i],sep="" ) )
			comparisons <- comparisons[-1]

			#factors genotype treatment levels
			extract_treatment <- function (x) paste( gsub( ".*_","",unlist( strsplit(x,"-")) ), collapse="-")
			extract_genotype <- function (x) gsub( "_.*","",unlist( strsplit(x,"-")) )[1]
			extract_treatment_counting <- unique( gsub( ".*_","",unlist( strsplit(g,"-")) ))
			treatments = sapply(comparisons, extract_treatment)
			genotypes = sapply(comparisons, extract_genotype)
			Exp.type = paste( Exp.type, "\nFactorial design:",length(unique(genotypes)),"X", length( extract_treatment_counting ), sep="" )

			# pairwise contrasts
			contrast1 <- makeContrasts(contrasts=comparisons[1], levels=design)
			for( kk in 2:length(comparisons) )
				contrast1<-  cbind(contrast1,makeContrasts(contrasts=comparisons[kk], levels=design)   )
			contrast.names = colnames(contrast1)

			# interaction contrasts
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

		} # if interaction terms	
		

		# if sample information is uploaded and user selected factors and comparisons
		if( !is.null(modelFactors) & length( selectedComparisons) >0  ) {
			Exp.type = paste("Model: ~", paste(modelFactors,collapse="+") )
			interactionTerm = FALSE # default value to be re-write if needed
			
			# model factors that does not contain interaction terms
			# modelFactors "genotype", "condition", "genotype:condition"
			keyModelFactors = modelFactors[ !grepl(":",modelFactors) ]
			
			# "genotype: control vs. mutant"
			factorsVector= gsub(":.*","",selectedComparisons) # corresponding factors for each comparison
			# remove factors not used in comparison, these are batch effects/pairs/blocks, 
			# keyModelFactors = keyModelFactors[ !is.na(match(keyModelFactors, factorsVector)) ]	
			
			# if a factor is selected both in block and main factors, then use it as block factor
			keyModelFactors = keyModelFactors[ is.na(match(keyModelFactors, blockFactor)) ]	
		
            #------Design matrix 			
			sampleInfo2 = sampleInfo[,keyModelFactors,drop=F] # remove factors not used.
			groups = apply(sampleInfo2,1, function(x) paste(x,collapse="_")  )
			g =  unique(groups)  
			
			design <- model.matrix(~ 0+factor(groups))
			colnames(design) <- gsub(".*)","",colnames(design))
		
			if( !is.null(rawCounts) && countsDEGMethods == 2) {  # voom
				v <- voom(rawCounts, design); 
				fit <- lmFit(v, design) 
			} else 
				fit <- lmFit(eset, design)
		
			fit <- eBayes(fit, trend=limmaTrend)	
			
			#-----------Making comaprisons
			if( length(keyModelFactors) != 2 | length(blockFactor) >1 )  { # if only one factor, or more than two then use all pairwise comparisons
				comparisons = gsub(".*: ","",selectedComparisons)
				comparisons = gsub(" vs\\. ","-",comparisons)
			} else if( length(keyModelFactors) == 2 ){ # two key factors

				if( sum(grepl(":",modelFactors)	>0) ) {  # interaction? 
					interactionTerm =TRUE
					# all pairwise comparisons
					comparisons = ""
					for( i in 1:(length(g)-1) )
						for (j in (i+1):length(g)) 
						if( strsplit(g[i],"_")[[1]][1] == strsplit(g[j],"_")[[1]][1]| strsplit(g[i],"_")[[1]][2] == strsplit(g[j],"_")[[1]][2]) # only compare WT_control vs. WT_treatment
							comparisons = c(comparisons,paste(g[j],"-",g[i],sep="" ) )
					comparisons <- comparisons[-1]
					
					# pairwise contrasts
					contrast1 <- makeContrasts(contrasts=comparisons[1], levels=design)
					if(length(comparisons)>1 )
					for( kk in 2:length(comparisons) )
						contrast1<-  cbind(contrast1,makeContrasts(contrasts=comparisons[kk], levels=design)   )
					contrast.names = colnames(contrast1)		
				
					# all possible interactions
					# interaction contrasts
					
					contrast2 = NULL
					contrast.names =""
					for ( kk in 1:(dim(contrast1)[2]-1) ) {
					for( kp in (kk+1):dim(contrast1)[2]) 
						#if( treatments[kp]== treatments[kk] ) 
						{  if(is.null(contrast2))  
							contrast2 = contrast1[,kp]- contrast1[,kk]  else					
							contrast2 = cbind(contrast2, contrast1[,kp]- contrast1[,kk] )
							
							contrast.names = c(contrast.names, paste0("Diff:",  colnames(contrast1)[kp], ".vs.", colnames(contrast1)[kk] ) )
						}   
					}
					colnames(contrast2)=contrast.names[-1]
					
					# remove nonsense contrasts from interactions
					 contrast2 = contrast2[,which(apply(abs(contrast2),2,max)==1),drop=F]
					 contrast2 = contrast2[,which(apply(abs(contrast2),2,sum)==4),drop=F]	
					 contrast2 = t( unique(t(contrast2)) ) # remove duplicate columns		
					
					# remove unwanted contrasts involving more than three levels in either factor
					keep= c()
					for( i in 1:dim(contrast2)[2]) {
						tem = rownames(contrast2)[ contrast2[ ,i ] != 0   ]
						tem1 = unique (  unlist(gsub("_.*","", tem) ) )
						tem2 = unique (  unlist(gsub(".*_","", tem) ) )
						if( length(tem1) == 2 & length(tem2) ==2 )
						keep = c(keep, colnames(contrast2) [i] )
					}
					contrast2 = contrast2[,keep,drop=F]
					comparisons2 = colnames(contrast2) 				 
				}

				# "stage: MN vs. EN"  -->  c("MN_AB-EN_AB", "EN_Nodule-EN_AB") 
				#  comparisons in all levels of the other factor 
				transformComparisons <- function (comparison1){
					tem = gsub(".*: ","",comparison1)
					tem = unlist(strsplit(tem, " vs\\. ") ) # control  mutant							
					factor1= gsub(":.*","",comparison1)

					ix = match(factor1, keyModelFactors) # 1: first factor, 2: 2nd factor
					otherFactor = keyModelFactors[3-ix]   # 3-1 = 2; 3-1=1
					otherFactorLevels = unique( sampleInfo2[,otherFactor] )				
					comparisons = c( )
					
					for (factorLevels in otherFactorLevels) {
						if( ix == 1){
							comparisons = c( comparisons, paste(paste0( tem, "_",factorLevels),collapse="-") )
						} else {
							comparisons = c( comparisons,  paste(paste0(factorLevels, "_", tem),collapse="-") )
						}
					}
					return(comparisons)		
				}	
				
				comparisons = unlist( sapply(selectedComparisons, transformComparisons ))
				comparisons = as.vector(comparisons)
			
			} # two factors
			
			# make contrasts
			contrast1 <- makeContrasts(contrasts=comparisons[1], levels=design)
			if(length(comparisons) >1 )
			for( kk in 2:length(comparisons) )
				contrast1<-  cbind(contrast1,makeContrasts(contrasts=comparisons[kk], levels=design)   )

			if( interactionTerm ) {  # if interaction terms
				contrast1 = cbind(contrast1,contrast2)
				contrast.names = c(colnames(contrast1), colnames(contrast2) )
				comparisons = c(comparisons,comparisons2)
			}
			
			# block design to remove batch effect or paired samples
			# corfit <- duplicateCorrelation(eset,design,block=targets$Subject)
			# corfit$consensus
			#Then this inter-subject correlation is input into the linear model fit:
			# fit <- lmFit(eset,design,block=targets$Subject,correlation=corfit$consensus)

			if(length(blockFactor) == 1 ) { # if a factor is selected as block
				#blockFactor = modelFactors[ which( is.na(ix) ) ]
				block = sampleInfo[, blockFactor]  # the column not used
			
				if( !is.null(rawCounts) && countsDEGMethods == 2) {  # voom
					v <- voom(rawCounts, design);
					corfit <- duplicateCorrelation(v,design,block= block)			
					fit <- lmFit(v, design,block=block,correlation=corfit$consensus) 
				} else {
					corfit <- duplicateCorrelation(eset,design,block= block)
					fit <- lmFit(eset, design,block=block,correlation=corfit$consensus)
				}
				fit <- eBayes(fit, trend=limmaTrend)			
	
			} # block factors


		} # use selected factors

		
		fit2 <- contrasts.fit(fit, contrast1)
		fit2 <- eBayes(fit2, trend=limmaTrend)
		#topTable(fit2, coef=1, adjust="BH")
		results <- decideTests(fit2, p.value=maxP_limma, lfc= log2(minFC_limma ))
		#vennDiagram(results[,1:5],circle.col=rainbow(5))

		# extract fold change for each comparison
		# there is issues with direction of foldchange. Sometimes opposite
		top <- function (comp) {
			tem <- topTable(fit2, number = 1e12,coef=comp,sort.by="M" ) 
			if(dim(tem)[1] == 0) { return (1) 
			} else	{ 			
				# compute fold change for the first gene (ranked by absolute value)
				tem2 = as.numeric( x[ which(rownames(x)== rownames(tem)[1]) , ] )
				names(tem2) = colnames(x) 
				
				if( 0){ # if not interaction term
					# compute real fold change in A-B comparison
					realFC =  mean ( tem2[which( groups == gsub("-.*","",comp))] ) - # average in A
								mean ( tem2[which( groups == gsub(".*-","",comp))] ) # average in B
					if( realFC * tem[1,1] <0 ) # if reversed
						tem[,1] <- -1*tem[,1]; # reverse direction if needed
				}				
				return( tem[,c(1,5)]) 
			}  
													
		}  # no significant gene returns 1, otherwise a data frame
		
			
		topGenes <- lapply(comparisons, top)
		topGenes <- setNames(topGenes, comparisons )
		ix <- which( unlist( lapply(topGenes, class) ) == "numeric")
		if( length(ix)>0) topGenes <- topGenes[ - ix ]
		# if (length(topGenes) == 0) topGenes = NULL;
	}
	return( list(results= results, comparisons = comparisons, Exp.type=Exp.type, topGenes=topGenes)) 
}

# Differential expression using DESeq2
DEG.DESeq2 <- function (  rawCounts,maxP_limma=.05, minFC_limma=2, selectedComparisons=NULL, sampleInfo = NULL,modelFactors=NULL, blockFactor = NULL){
	library(DESeq2) # count data analysis
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
	
	
	Exp.type = paste(length(g)," sample groups detected.")
	comparisons = ""
	for( i in 1:(length(g)-1) )
		for (j in (i+1):length(g)) 
		comparisons = c(comparisons,paste(g[j],"-",g[i],sep="" ) )
	comparisons <- comparisons[-1]

	colData = cbind(colnames(rawCounts), groups )

	# no sample file, but user selected comparisons using column names
	if( is.null(modelFactors) & length( selectedComparisons) >0  ) 	
		comparisons = selectedComparisons
		
	comparisons2 = comparisons	 # this is for showing comparison names, which might be different from internally	
	# Set up the DESeqDataSet Object and run the DESeq pipeline
	dds = DESeqDataSetFromMatrix(countData=rawCounts,
								colData=colData,
								design=~groups)
								
	
	if( is.null(modelFactors)  ) 
		dds = DESeq(dds)  	else  
	{    # using selected factors and comparisons
		# build model
		modelFactors = c(modelFactors,blockFactor) # block factor is just added in. 
		
		factors = modelFactors   # selected facters and interactions: c( "strain", "treatment",  "strain:treatment")
		factors = factors[ !grepl(":",factors )]   # non-interaction terms
		
		colData = sampleInfo  
		factorsCoded = toupper(letters )[1: dim(colData)[2] ]   # Factors are encoded as "A", "B", "C"; this avoid illigal letters
		names(factorsCoded) =  colnames(colData)  # this is for look up; each column of sampleInfo  
		colnames(colData) = factorsCoded # all columns named A B C D 
		# base model
        DESeq2.Object= paste("dds = DESeqDataSetFromMatrix(countData=rawCounts, colData=colData, design=~ ", 
							 paste( factorsCoded[factors],collapse="+")) # only use selected factors		
		Exp.type = paste("Model: ~", paste(modelFactors,collapse="+") )
		# interaction terms like strain:treatment
		Interactions = modelFactors[ grepl(":",modelFactors )]
		# create model
		if( length(Interactions)>0 ) { # if there is interaction
			for( interactionTerms in Interactions) {
				interactingFactors = unlist(strsplit(interactionTerms,":" ) )  # split strain:treatment as "strain" and "mutant"
				tem = paste(factorsCoded [ interactingFactors ],collapse=":")   # convert "strain:mutant" to "A:B"
				DESeq2.Object = paste(DESeq2.Object, "+",tem)
			}			
		}
		DESeq2.Object= paste( DESeq2.Object, ")") # ends the model

		eval(parse(text = DESeq2.Object) )
		dds = DESeq(dds)  # main function		
		
		# comparisons 
		# "group: control vs. mutant"
		comparisons = gsub(".*: ","",selectedComparisons)
		comparisons = gsub(" vs\\. ","-",comparisons)
		factorsVector= gsub(":.*","",selectedComparisons) # corresponding factors for each comparison
		comparisons2 = comparisons

		# Note that with interaction terms, not all meaningful comparisons is listed for selection. 
		# this is complex. Only under reference level.
		
		# comparisons due to interaction terms
		if( length(Interactions)>0 ) { # if there is interaction
			interactionComparisons = resultsNames(dds)
			interactionComparisons = interactionComparisons[ grepl("\\.",interactionComparisons )   ]
	
			comparisons = c(comparisons,interactionComparisons )
			
			# translate comparisons generated in interaction terms back to real factor names
			interactionComparisons2 = interactionComparisons
			for ( i in 1:length(interactionComparisons2 ) ) {
				tem = unlist(strsplit(interactionComparisons2[i],"\\." ) )
				tem_factors = substr(tem,1,1) 
				tem_factors = names(factorsCoded)[factorsCoded == tem_factors]  # get the first letters and translate into real factor names
				interactionComparisons2[i] <- paste0( "Diff:",tem_factors[1], "_",substr(tem[1],2,nchar(tem[1]) ),".",
													          tem_factors[2], "_",substr(tem[2],2,nchar(tem[2]) ) 
													)				
			}
			comparisons2 = c(comparisons2,interactionComparisons2 )
		}
	}	
	
	result1 = NULL; allCalls = NULL;
	topGenes = list(); pk = 1 # counter
	pp=0 # first results?
	for( kk in 1:length(comparisons) ) {
		tem = unlist( strsplit(comparisons[kk],"-") )
		
		if(is.null(modelFactors)) 
			selected = results(dds, contrast=c("groups", tem[1], tem[2]) )   else {
			if(!grepl("\\.", comparisons[kk] ) )    # if not interaction term: they contain .  interaction term
				selected = results(dds, contrast=c( factorsCoded[ factorsVector[kk]  ],tem[1], tem[2]) ) else # either A, B, C ...
				selected = results(dds, name=comparisons[kk] ) # interaction term
			}

		selected$calls =0   
		selected$calls [which( selected$log2FoldChange > log2(minFC_limma) & selected$padj < maxP_limma ) ]  <-  1
		selected$calls [ which( selected$log2FoldChange <  -log2(minFC_limma) & selected$padj < maxP_limma ) ] <-  -1
		colnames(selected)= paste( as.character(comparisons2[kk]), "___",colnames(selected),sep="" )
		selected = as.data.frame(selected)
		if (pp==0){  # if first one with significant genes, collect gene list and Pval+ fold
		result1 = selected; pp = 1; 
		# selected[,2] <- -1 * selected[,2] # reverse fold change direction
		topGenes[[1]] = selected[,c(2,6)]; 
		names(topGenes)[1] = comparisons2[kk]; } else 
			{ result1 = merge(result1,selected,by="row.names"); 
				rownames(result1) = result1[,1]; 
				result1 <- result1[,-1]
				pk= pk+1; 
				# selected[,2] <- -1 * selected[,2] # reverse fold change direction
				topGenes[[pk]] = selected[,c(2,6)]; 
				names(topGenes)[pk] = comparisons2[kk]; 
			}
	}
	#if( length(comparisons) == 1) topGenes <- topGenes[[1]] # if only one comparison, topGenes is not a list, just a data frame itself.
	if(! is.null(result1)) { 
	# note that when you only select 1 column from a data frame it automatically converts to a vector. drop =FALSE prevents that.
	allCalls = as.matrix( result1[,grep("calls",colnames(result1)), drop = FALSE  ] )
	colnames(allCalls)= gsub("___.*","", colnames(allCalls))
	colnames(allCalls)= gsub("\\.","-", colnames(allCalls)) # note that samples names should have no "."
	}
	return( list(results= allCalls, comparisons = comparisons2, Exp.type=Exp.type, topGenes=topGenes)) 
}

# Find enriched TF binding motifs in promoters
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
 
if(0) {  # testing

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

	tem = apply(x,1,function(y) sum(y>2) )
	
	
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



################################################################
#   Server function
################################################################
	
shinyServer(
function(input, output,session) {
  
 #----------------------------------------------- 
 # Available datasets 
 # readData()$data: transformed data, readData()$rawCounts: Counts data. NULL if non-count data.
 # converted():id conversion results with many components such as: converted()$originalIDs,  converted()$IDs: converted IDs
               #,converted()$species,   converted()$speciesMatched,  converted()$conversionTable
 # allGeneInfo(): returns all information in the geneInfo file for each gene
 # geneSets(): gene set as a list for pathway analysis
  
	options(shiny.maxRequestSize = 50*1024^2) # 50MB file max for upload

	observe({  updateSelectInput(session, "selectOrg", choices = speciesChoice )      })
	observe({  updateSelectInput(session, "heatColors1", choices = colorChoices )      })
	observe({  updateSelectInput(session, "distFunctions", choices = distChoices )      })
	observe({  updateSelectInput(session, "hclustFunctions", choices = hclustChoices )      })

	################################################################
	#   Read data
	################################################################
 
	# read data file and do filtering and transforming
	readData <- reactive ({
		library(edgeR) # count data D.E.
		library(DESeq2) # count data analysis
		inFile <- input$file1
		inFile <- inFile$datapath

		if(is.null(input$file1) && input$goButton == 0)   return(NULL)
		if(is.null(input$file1) && input$goButton > 0 )   inFile = demoDataFile
		tem = input$dataFileFormat
		if(!is.null(input$dataFileFormat)) # these are needed to make it responsive to changes
			if(input$dataFileFormat== 1){  
				tem = input$minCounts 
				tem= input$NminSamples
				tem = input$countsLogStart
				tem = input$CountsTransform 
			}
		if(!is.null(input$dataFileFormat))
			if(input$dataFileFormat== 2){ 
				tem = input$transform; 
				tem = input$logStart; 
				tem= input$lowFilter 
		}

		isolate({
			withProgress(message="Reading and pre-processing ", {
				if (is.null( input$dataFileFormat )) return(NULL)
				dataTypeWarning =0
				dataType =c(TRUE)

				#---------------Read file
				x <- read.csv(inFile)	# try CSV
				if(dim(x)[2] <= 2 )   # if less than 3 columns, try tab-deliminated
					x <- read.table(inFile, sep="\t",header=TRUE)	
				#-------Remove non-numeric columns, except the first column
				
				for(i in 2:dim(x)[2])
					dataType = c( dataType, is.numeric(x[,i]) )
				if(sum(dataType) <=2) return (NULL)  # only less than 2 columns are numbers
				x <- x[,dataType]  # only keep numeric columns

				dataSizeOriginal = dim(x); dataSizeOriginal[2] = dataSizeOriginal[2] -1
				
				x[,1] <- toupper(x[,1])
				x[,1] <- gsub(" ","",x[,1]) # remove spaces in gene ids
				x = x[order(- apply(x[,2:dim(x)[2]],1,sd) ),]  # sort by SD
				x <- x[!duplicated(x[,1]) ,]  # remove duplicated genes
				rownames(x) <- x[,1]
				x <- as.matrix(x[,c(-1)])

				# missng value for median value
				if(sum(is.na(x))>0) {# if there is missing values
					rowMeans <- apply(x,1, function (y)  median(y,na.rm=T))
					for( i in 1:dim(x)[1] )
					x[i, which( is.na(x[i,]) )  ]  <- rowMeans[i]
				}

				# Compute kurtosis
				mean.kurtosis = mean(apply(x,2, kurtosis))

				if (input$dataFileFormat == 2 ) {  # if FPKM, microarray
					incProgress(1/3,"Pre-processing data")

					if ( is.integer(x) ) dataTypeWarning = 1;  # Data appears to be read counts

					#-------------filtering
					tem <- apply(x,1,max)
					x <- x[which(tem > input$lowFilter),]  # max by row is at least 		
					x <- x[which(apply(x,1, function(y) max(y)- min(y) ) > 0  ),]  # remove rows with all the same levels

					#--------------Log transform
					# Takes log if log is selected OR kurtosis is big than 100
					if ( (input$transform == TRUE) | (mean.kurtosis > kurtosis.log ) ) 
						x = log(x+abs( input$logStart),2)

					tem <- apply(x,1,sd) 
					x <- x[order(-tem),]  # sort by SD
					rawCounts = NULL
				} else {  # counts data
					incProgress(1/3, "Pre-processing counts data")
					tem = input$CountsDEGMethod; tem = input$countsTransform
					# data not seems to be read counts
					if(!is.integer(x) & mean.kurtosis < kurtosis.log ) {
						dataTypeWarning = -1
					}
					# not used as some counts data like those from CRISPR screen
					#validate(   # if Kurtosis is less than a threshold, it is not read-count
					#	need(mean.kurtosis > kurtosis.log, "Data does not seem to be read count based on distribution. Please double check.")
					# )
					x <- round(x,0) # enforce the read counts to be integers. Sailfish outputs has decimal points.
					#x <- x[ which( apply(x,1,max) >= input$minCounts ) , ] # remove all zero counts
					

					# remove genes if it does not at least have minCounts in at least NminSamples
					#x <- x[ which( apply(x,1,function(y) sum(y>=input$minCounts)) >= input$NminSamples ) , ]  # filtering on raw counts
					# using counts per million (CPM) for filtering out genes.
                                             # CPM matrix                  #N samples > minCounts
					x <- x[ which( apply( cpm(DGEList(counts = x)), 1,  function(y) sum(y>=input$minCounts)) >= input$NminSamples ) , ] 
					

					if(0){  # disabled
						# remove genes with low expression by counts per million (CPM)
						dge <- DGEList(counts=x); dge <- calcNormFactors(dge)
						myCPM <- cpm(dge, prior.counts = 3 )
						x <- x[which(rowSums(  myCPM > input$minCounts)  > 1 ),]  # at least two samples above this level
						rm(dge); rm(myCPM)
					}
					rawCounts = x; # ??? 
					# construct DESeqExpression Object
					# colData = cbind(colnames(x), as.character(detectGroups( colnames(x) )) )
					tem = rep("A",dim(x)[2]); tem[1] <- "B"   # making a fake design matrix to allow process, even when there is no replicates
					colData = cbind(colnames(x), tem )
					colnames(colData)  = c("sample", "groups")
					dds <- DESeqDataSetFromMatrix(countData = x, colData = colData, design = ~ groups)
					dds <- estimateSizeFactors(dds) # estimate size factor for use in normalization later for started log method

					incProgress(1/2,"transforming raw counts")
					# regularized log  or VST transformation
					if( input$CountsTransform == 3 ) { # rlog is slow, only do it with 10 samples
						x <- rlog(dds, blind=TRUE); x <- assay(x) } 
					else {
						if ( input$CountsTransform == 2 ) {    # vst is fast but aggressive
							x <- vst(dds, blind=TRUE)
							x <- assay(x)  
						} else{  # normalized by library sizes and add a constant.
							x <- log2( counts(dds, normalized=TRUE) + input$countsLogStart )   # log(x+c) 
							#x <- cpm(DGEList(counts = x),log=TRUE, prior.count=input$countsLogStart )  #log CPM from edgeR
							#x <- x-min(x)  # shift values to avoid negative numbers
						}
					}
				}
				dataSize = dim(x);
				incProgress(1, "Done.")
				finalResult <- list(data = as.matrix(x), mean.kurtosis = mean.kurtosis, rawCounts = rawCounts, dataTypeWarning=dataTypeWarning, dataSize=c(dataSizeOriginal,dataSize) )
				return(finalResult)
			})
		})
	})

	readSampleInfo <- reactive ({
		inFile <- input$file2
		inFile <- inFile$datapath

		if(is.null(input$file2) && input$goButton == 0)   return(NULL)
		if(is.null(readData() ) ) return(NULL)

		isolate({
				if (is.null( input$dataFileFormat )) return(NULL)
				dataTypeWarning =0
				dataType =c(TRUE)

				#---------------Read file
				x <- read.csv(inFile,row.names=1,header=T,colClasses="character")	# try CSV
				if(dim(x)[2] <= 2 )   # if less than 3 columns, try tab-deliminated
					x <- read.table(inFile, row.names=1,sep="\t",header=TRUE,colClasses="character")	
				#----------------Matching with column names of expression file
				ix = match(toupper(colnames(readData()$data)), toupper(colnames(x)) ) 
				ix = ix[which(!is.na(ix))] # remove NA

				validate(
				  need(length(unique(ix) ) == dim(readData()$data)[2] 
				       & dim(x)[1]>=1  # at least one row
					 ,"Error!!! Sample information file not recognized. Sample names must be exactly the same. Each row is a factor. Each column represent a sample.  Please see documentation on format.")
				)
				if( length(unique(ix) ) == dim(readData()$data)[2]) { # matches exactly
				x = x[,ix]
				# if the levels of different factors are the same, it may cause problems
				if( sum( apply(x, 1, function(y) length(unique(y)))) > length(unique(unlist(x) ) ) ) {
					tem2 =apply(x,2, function(y) paste0( names(y),y)) # factor names are added to levels
					rownames(tem2) = rownames(x)
					x <- tem2				
				}
				return(t( x ) )			
				} else retrun(NULL)
							
				
		})
	})
	
	output$sampleInfoTable <- renderTable({

		if (is.null(input$file2) )   return(NULL)
		isolate({
			t(readSampleInfo() )
		})
	},include.rownames=TRUE) 	
	
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
			tem = paste (tem, "  ------Warning!!! Data does not look like read counts data. Please select appropriate data type in the previous page and reload file.")}
		tem
	})	

	# provide info on number of genes passed filter
	output$nGenesFilter <- renderText({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		inFile <- input$file1
		if( is.null(readData()) ) return(NULL)
		tem = readData()$dataSize
		return( paste(tem[4], "samples.", tem[3], "of the ", tem[1], "genes pased filter."  ) )
	})	

	# Show info on file format	
	output$fileFormat <- renderUI({
		i = "<h3>Done. Ready to load data files.</h3>"
		i = c(i,"Users can upload a CSV or tab-delimited text file with the first column as gene IDs. 
		For RNA-seq data, read count per gene is recommended.
		Also accepted are normalized expression data based on FPKM, RPKM, or DNA microarray data. iDEP can convert most types of common gene IDs to Ensembl gene IDs, which is used 
			internally for enrichment and pathway analyses. iDEP parses column names to define sample groups. To define 3 biological samples (Control,
		TreatmentA, TreatmentB) with 2 replicates each, column names should be:")
		i = c(i," <strong> Ctrl_1, Ctrl_2, TrtA_1, TrtA_2, TrtB_1, TrtB_2</strong>.") 
		i = c(i,"For factorial design, use underscore \"_\" to separate factors such as genetic background 
		(wide type vs. mutant:WT vs. Mu) and experimental condition (Ctrl vs. Trt). 
		Currently, only two factors are allowed. To define an 2x2 factorial design, use column names like:")
		i = c(i,"<strong>WT_Ctrl_1, WT_Ctrl_2, WT_Trt_1, WT_Trt_2, Mu_Ctrl_1, Mu_Ctrl_2, Mu_Trt_1, Mu_Trt_2</strong>") 
		
		HTML(paste(i, collapse='<br/>') )
	})
	# this defines an reactive object that can be accessed from other rendering functions

	converted <- reactive({
		if (is.null(input$file1) && input$goButton == 0)    return(NULL)
		tem = input$selectOrg;
		isolate( {
		#   withProgress(message="Converting gene ids", {
		# cat (paste("\nID:",input$selectOrg) )
		convertID(rownames(readData()$data ),input$selectOrg, input$selectGO );
		# })
		}) 
	})

	# this defines an reactive object that can be accessed from other rendering functions
	allGeneInfo <- reactive({
		if (is.null(input$file1) && input$goButton == 0)    return(NULL)
		tem = input$selectOrg;
		isolate( {
		withProgress(message="Looking up gene annotation", {
		geneInfo(converted(),input$selectOrg); 
				})
		})
	})

	convertedData <- reactive({
		if (is.null(input$file1) && input$goButton == 0) return()  
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
		####################################
		if( is.null(converted() ) ) return( readData()$data) # if id or species is not recognized use original data.
		isolate( {  
			withProgress(message="Converting data ... ", {
				mapping <- converted()$conversionTable
				# cat (paste( "\nData:",input$selectOrg) )
				x =readData()$data
				rownames(x) = toupper(rownames(x))
				# any gene not recognized by the database is disregarded
				# x1 = merge(mapping[,1:2],x,  by.y = 'row.names', by.x = 'User_input')
				# the 3 lines keeps the unrecogized genes using original IDs
				x1 = merge(mapping[,1:2],x,  by.y = 'row.names', by.x = 'User_input', all.y=TRUE)
				ix = which(is.na(x1[,2]) )
				x1[ix,2] = x1[ix,1] # original IDs used
				#write.csv(x1,"tem.csv")
				#multiple matched IDs, use the one with highest SD
				tem = apply(x1[,3:(dim(x1)[2])],1,sd)
				x1 = x1[order(x1[,2],-tem),]
				x1 = x1[!duplicated(x1[,2]) ,]
				rownames(x1) = x1[,2]
				x1 = as.matrix(x1[,c(-1,-2)])
				tem = apply(x1,1,sd)
				x1 = x1[order(-tem),]  # sort again by SD
				incProgress(1, "Done.")
			})
		return(x1)
		})
	})
	
	GeneSets <- reactive({
		if (is.null(input$file1) && input$goButton == 0)	return()
		tem = input$selectOrg
		tem = input$selectGO
		tem =input$minSetSize; tem= input$maxSetSize; 
		isolate( {
			if(input$selectOrg == "NEW" && !is.null(input$gmtFile) ) # new species 
			{     inFile <- input$gmtFile; inFile <- inFile$datapath
				return( readGMTRobust(inFile) ) }else
		return( readGeneSets( converted(), convertedData(), input$selectGO,input$selectOrg,c(input$minSetSize, input$maxSetSize)  ) ) }) 
	})

####### [TODO] Kevin Indentation Work 10/5 #######

	output$contents <- renderTable({
	   inFile <- input$file1
		inFile <- inFile$datapath
		if (is.null(input$file1) && input$goButton == 0)   return(NULL)
	#    if (is.null(input$file1) && input$goButton > 0 )   inFile = "expression1_no_duplicate.csv"
		if (is.null(input$file1) && input$goButton > 0 )   inFile = demoDataFile

		tem = input$selectOrg
		isolate({
		x <- read.csv(inFile)
		if(dim(x)[2] <= 2 ) x <- read.table(inFile, sep="\t",header=TRUE)	# not CSV
		#x <- readData()$data
		 x[1:20,]
		})
	  },include.rownames=FALSE)

# show first 20 rows of data
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

# show first 20 rows of processed data; not used
	output$debug <- renderTable({
      if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	  tem = input$selectOrg; tem = input$lowFilter ; tem = input$transform
	x <- convertedData()
	#tem = GeneSets()
	
	return( as.data.frame (x[1:20,] ))
	},include.rownames=TRUE)


################################################################
#   Pre-process
################################################################
	
	output$EDA <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	####################################

	
    x <- readData()$data
	 par(mfrow=c(3,1))
	par(mar=c(14,6,4,4))
	myColors = rainbow(dim(x)[2])
	plot(density(x[,1]),col = myColors[1], lwd=2,
	  xlab="Expresson values", ylab="Density", main= "Distribution of transformed data",
	  cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, ylim=c(0, max(density(x[,1])$y)+.02 ) )
	  
	for( i in 2:dim(x)[2] )
	lines(density(x[,i]),col=myColors[i], lwd=2 )
	if(dim(x)[2]< 31 ) # if too many samples do not show legends
		legend("topright", cex=1.2,colnames(x), lty=rep(1,dim(x)[2]), col=myColors )	
   # boxplot of first two samples, often technical replicates
   
	boxplot(x, las = 2, ylab="Transformed expression levels", main="Distribution of transformed data"
		,cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=2)
	plot(x[,1:2],xlab=colnames(x)[1],ylab=colnames(x)[2], main="Scatter plot of first two samples",cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
	
   }, height = 1600, width = 800)
   
	output$genePlot <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }

	
	tem = input$geneSearch ; tem = input$genePlotBox; tem = input$useSD
	isolate({
    x <- convertedData()
	
	Symbols <- rownames(x)
	if( input$selectOrg != "NEW") {
		ix = match( rownames(x), allGeneInfo()[,1])
		if( sum( is.na(allGeneInfo()$symbol )) != dim(allGeneInfo() )[1] ) {  # symbol really exists? 
			Symbols = as.character( allGeneInfo()$symbol[ix] )
			Symbols[which( nchar(Symbols) <= 2 ) ] <- rownames(x) [which( nchar(Symbols) <= 2 ) ] 
			}
	   }
	x = as.data.frame(x)
	x$Genes = Symbols
    #write.csv(x,"tem.csv")
	# Search for genes
	#ix = grep("HOXA",toupper(x$Genes) )
	# ix = grep(toupper(input$geneSearch),toupper(x$Genes))  # sox --> Tsox  
	# matching from the beginning of symbol
	ix = which(regexpr(  paste("^" , toupper(input$geneSearch),sep="")   ,toupper(x$Genes)) > 0)
	
	if(grepl(" ", input$geneSearch)  )  # if there is space character, do exact match
		ix = match(gsub(" ","", toupper(input$geneSearch)),x$Genes)
	
	# too few or too many genes found
	if(length(ix) == 0 | length(ix) > 50 ) return(NULL)
	  # no genes found
	 	 
	mdf = melt(x[ix,],id.vars="Genes", value.name="value", variable.name="samples")
	# bar plot of individual samples
	p1 <- ggplot(data=mdf, aes(x=samples, y=value, group = Genes, shape=Genes, colour = Genes)) +
		geom_line() +
		geom_point( size=5,  fill="white")+ #shape=21  circle
		#theme(axis.text.x = element_text(size=16,angle = 45, hjust = 1)) +
		labs(y="Transformed expression level") +
		coord_cartesian(ylim = c(0, max(mdf$value)))
	p1 <- p1 + theme(plot.title = element_text(size = 16,hjust = 0.5)) + # theme(aspect.ratio=1) +
	 theme(axis.text.x = element_text(angle=45, size = 16, hjust=1),
	       axis.text.y = element_text( size = 16),
		   axis.title.x = element_blank(),
		   axis.title.y = element_text( size = 16) ) +
	theme(legend.text=element_text(size=12))	
		
		
	#ggplotly(p) %>% layout(margin = list(b = 250,l=100))  # prevent cutoff of sample names

	# Barplot with error bars
	mdf$count = 1
	g = detectGroups(mdf$samples)
	Means = aggregate(mdf$value,by=list( g, mdf$Genes ), FUN = mean, na.rm=TRUE  )
	SDs = aggregate(mdf$value,by=list( g, mdf$Genes ), FUN = sd, na.rm=TRUE  )
	Ns = aggregate(mdf$count, by= list(g, mdf$Genes) , FUN = sum  )
	summarized = cbind(Means,SDs[,3],Ns[,3])
	colnames(summarized)= c("Samples","Genes","Mean","SD","N")
	summarized$SE = summarized$SD / sqrt(summarized$N)	
		
	#http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization
	p2 <- ggplot(summarized, aes(x=Genes, y=Mean,fill=Samples) ) + # data & aesthetic mapping
		geom_bar(stat="identity", position=position_dodge()) + # bars represent average
		geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.2,position=position_dodge(.9)) +
		labs(y="Expression Level") 
	if(input$useSD == 1) { 
	p2 <- ggplot(summarized, aes(x=Genes, y=Mean,fill=Samples) ) + # data & aesthetic mapping
		geom_bar(stat="identity", position=position_dodge()) + # bars represent average
		geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.2,position=position_dodge(.9)) +
		labs(y="Expression Level") 
	}
	
	p2 <- p2 +  theme(plot.title = element_text(size = 16,hjust = 0.5)) + # theme(aspect.ratio=1) +
	 theme(axis.text.x = element_text(angle=45, size = 16, hjust=1),
	       axis.text.y = element_text( size = 16),
		   axis.title.x = element_blank(),
		   axis.title.y = element_text( size = 16) ) +
	theme(legend.text=element_text(size=16))
	
	if( input$genePlotBox == 1)  p1 else p2
	
	})
   })
   

	processedData <- reactive({
      if (is.null(input$file1) && input$goButton == 0)    return()
	  
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	####################################
	 
		if(input$selectOrg == "NEW") return(  convertedData() ) else { 

			withProgress(message="Preparing data for download ", {
			
			tem <-  merge(allGeneInfo()[,c('ensembl_gene_id','symbol')], round(convertedData(),4),by.x="ensembl_gene_id", by.y ="row.names", all.y=T )   
			tem[,2] = paste(" ",tem[,2]) # add space to gene symbol to avoid auto convertion of symbols to dates by Excel 
			incProgress(1/2, "mapping")
			#tem[,1] = paste(" ",tem[,1]) # add space
			 tem2 = merge( converted()$conversionTable, tem, by.x = "ensembl_gene_id",by.y="ensembl_gene_id")
			 tem2 <- tem2[,-3] # remove species column
			 ix = which( tem2[,3] == "  NA") # remove NA's 
			 tem2[ix,3] <- ""
			 incProgress(1, "Done.")
			 # add original data
			# tem3 <- merge( readData()$data, tem2, by.x = "row.names", by.y = "User_input", all.x=TRUE )
			# return(converted()$conversionTable) #return(readData()$data)
			# colnames(tem)[1] = "Ensembl_or_Original_ID"
			})
			return(tem2)
			
			
		}
	  })
  
	output$downloadProcessedData <- downloadHandler(
		filename = function() {"processed_Data.csv"},
		content = function(file) {
      write.csv( processedData(), file, row.names=FALSE )	    
	})
 
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
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	####################################
	
    par(mar=c(16,2,2,2))
    x <- readData()$rawCounts
	barplot( colSums(x)/1e6, col="green",las=3, cex.axis=1.3, cex =1.5, main="Total read counts (millions)")

},width=400) # height is automatic, this enables the display of other plots below.


################################################################
#   Heatmaps
################################################################

	output$listFactorsHeatmap <- renderUI({
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	
    if (is.null(input$file2) ) # if sample info is uploaded and correctly parsed.
       { return(NULL) }	 else { 
	  selectInput("selectFactorsHeatmap", label="Sample color bar:",choices=colnames(readSampleInfo())
	     )   } 
	})  

	# conventional heatmap.2 plot
	output$heatmap1 <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    x <- readData()$data   # x = read.csv("expression1.csv")
	withProgress(message="Reading and pre-processing ", {
	n=input$nGenes
	#if(n>6000) n = 6000 # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
	# this will cutoff very large values, which could skew the color 
	x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	
	# standardize by gene
	if(input$geneNormalize) 
		x <- x / apply(x,1,sd)
	# standardize by sample
	if(input$sampleNormalize){
		x <- scale(x, center = TRUE, scale = apply(x,2,sd)) 
	}
	
	
	cutoff = median(unlist(x)) + input$heatmapCutoff * sd (unlist(x)) 
	x[x>cutoff] <- cutoff
	cutoff = median(unlist(x)) - input$heatmapCutoff *sd (unlist(x)) 
	x[x< cutoff] <- cutoff
	
    groups = detectGroups(colnames(x) )
	# if sample info file is uploaded us that info:

	if(!is.null(input$file2) &&  !is.null(input$selectFactorsHeatmap) ) { 
		ix = match(input$selectFactorsHeatmap, colnames(readSampleInfo() ) ) 
		groups = readSampleInfo()[,ix]
	}	
	
	groups.colors = rainbow(length(unique(groups) ) )

	#http://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
	lmat = rbind(c(0,4),c(0,1),c(3,2),c(5,0))
	lwid = c(2,6) # width of gene tree; width of heatmap
	lhei = c(1.5,.2,8,1.1)
	#layout matrix
	#		 [,1] [,2]
	#	[1,]    0    4
	#	[2,]    0    1
	#	[3,]    3    2
	#	[4,]    5    0
	# 4--> column tree; 1--> column color bar; 2--> heatmap; 3-> row tree; 5--> color key.
	# height of 4 rows is specified by lhei; width of columns is given by lwid


	par(mar = c(5, 4, 1.4, 0.2))
	

	if( n>110) 
	heatmap.2(x, distfun = distFuns[[as.integer(input$distFunctions)]]
		,hclustfun=hclustFuns[[as.integer(input$hclustFunctions)]]
		,Colv=!input$noSampleClustering
		#col=colorpanel(75,"green","black","magenta")  ,
		#col=bluered(75),
		#col=greenred(75), 
		,col= heatColors[as.integer(input$heatColors1),]
		,density.info="none", trace="none", scale="none", keysize=.5
		,key=T, symkey=F
		,ColSideColors=groups.colors[ as.factor(groups)]
		,labRow=""
		,margins=c(10,0)
		,srtCol=45
		,cexCol=2  # size of font for sample names
		,lmat = lmat, lwid = lwid, lhei = lhei
		)

	if( n<=110) 
	heatmap.2(x, distfun =  distFuns[[as.integer(input$distFunctions)]]
		,hclustfun=hclustFuns[[as.integer(input$hclustFunctions)]]
		,Colv=!input$noSampleClustering		
		,col= heatColors[as.integer(input$heatColors1),], density.info="none", trace="none", scale="none", keysize=.5
		,key=T, symkey=F,
		#,labRow=labRow
		,ColSideColors=groups.colors[ as.factor(groups)]
		,margins=c(18,12)
		,cexRow=1
		,srtCol=45
		,cexCol=1.8  # size of font for sample names
		,lmat = lmat, lwid = lwid, lhei = lhei
	)
	
	
	par(lend = 1)           # square line ends for the color legend
	add_legend("topleft",
		legend = unique(groups), # category labels
		col = groups.colors[ unique(as.factor(groups))],  # color key
		lty= 1,             # line style
		lwd = 10            # line width
	)
	
	incProgress(1,"Done")
	})

} , height = 800, width = 400 )  

# interactive heatmap with plotly
	output$heatmap <- renderPlotly({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
   # x <- readData()$data   # x = read.csv("expression1.csv")
    x <- convertedData()
	withProgress(message="Rendering heatmap ", {
	n=input$nGenesPlotly
	#if(n>6000) n = 6000 # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
	# this will cutoff very large values, which could skew the color 
	x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
	x[x>cutoff] <- cutoff
	cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
	x[x< cutoff] <- cutoff

	# adding gene symbol
	ix <- match( rownames(x), allGeneInfo()[,1])
	geneSymbols <- as.character( allGeneInfo()$symbol)[ix]
	# if missing or duplicated, use Ensembl ID
	ix <- which(nchar( geneSymbols) <=1 | duplicated(geneSymbols ) );	geneSymbols[ ix ] <- rownames(x)[ix]
	rownames( x) = geneSymbols;

	incProgress(1/2, "Clustering of genes")	
	clust <- x %>% 
	  dist2() %>% 
	  hclust2()

	# Get order
	ord <- clust$order

	# Re-arrange based on order
	df <- t( x[ord,] )%>%
	   melt()
	   
	colnames(df)[1:2] <- c("X","Y")
    colorNames = unlist(strsplit(tolower(rownames(heatColors)[ as.integer(input$heatColors1)   ]),"-" ) )
	p <- df %>%
	  ggplot(aes(X, Y, fill = value)) + 
		   geom_tile()+ scale_fill_gradient2(low = colorNames[1], mid = colorNames[2],high = colorNames[3]) +
		   theme(axis.title.y=element_blank(),   # remove y labels
		   # axis.text.y=element_blank(),  # keep gene names for zooming
			axis.ticks.y=element_blank(),
			axis.title.x=element_blank()) +
			theme(axis.text.x = element_text(size=10,angle = 45, hjust = 1))

		incProgress(1,"Done")
	ggplotly(p) %>% 
		layout(margin = list(b = 150,l=200))  # prevent cutoff of sample names

	})
  })  
  
	heatmapData <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    x <- readData()$data

	
	n=input$nGenes
	#if(n>6000) n = 6000 # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
	# this will cutoff very large values, which could skew the color 
	x1 = x[1:n,] 
	x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	cutoff = median(unlist(x)) + 4*sd (unlist(x)) 
	x[x>cutoff] <- cutoff
	cutoff = median(unlist(x)) - 4*sd (unlist(x)) 
	x[x< cutoff] <- cutoff
	
	groups = detectGroups(colnames(x) )
	groups.colors = rainbow(length(unique(groups) ) )

	#pdf(file=NULL,width =700, height =700)
	hy <- heatmap.2(x, distfun = distFuns[[as.integer(input$distFunctions)]]
		,hclustfun=hclustFuns[[as.integer(input$hclustFunctions)]]
		,density.info="none", trace="none", scale="none")
    #dev.off()
	
	# if not new species, add gene symbol
	if( input$selectOrg == "NEW") return(NULL) else { 
		x1 <- x1[ rev( hy$rowInd),hy$colInd]
		# add gene symbol
		ix = match( rownames(x1), allGeneInfo()[,1])
		x1 <- cbind(as.character( allGeneInfo()$symbol)[ix],x1)
		return( x1 )
	}
	

	
  })  
  
	output$downloadData <- downloadHandler(
		filename = function() {"heatmap.csv"},
		content = function(file) {
			write.csv(heatmapData(), file)
	    }
	)

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
		size = 15, hjust = 1))+
	 theme(axis.text.y = element_text( 
		size = 15))+
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
	 guides(fill = FALSE) # + ggtitle("Pearson's Correlation Coefficient (all genes)")

 
  }  )#, height = 500, width = 500)

################################################################
#   PCA
################################################################
	output$listFactors <- renderUI({
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	
      if (is.null(input$file2) )
       { return(HTML("Upload a sample info file to customize this plot.") ) }	 else { 
	  selectInput("selectFactors", label="Color:",choices=colnames(readSampleInfo())
	     )   } 
	})
	output$listFactors2 <- renderUI({
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	
      if (is.null(input$file2) )
       { return(NULL) }	 else { 
	   tem <- colnames(readSampleInfo() )
	   if(length(tem)>1) { tem2 = tem[1]; tem[1] <- tem[2]; tem[1] = tem2; } # swap 2nd factor with first
	  selectInput("selectFactors2", label="Shape:",choices=tem)
	        } 
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
	if(is.null(input$file2)) { 
		p=ggplot(pcaData, aes(PC1, PC2, color=Type, shape = Type)) + geom_point(size=5) 
		} else {
		pcaData = cbind(pcaData,readSampleInfo() )
		p=ggplot(pcaData, aes_string("PC1", "PC2", color=input$selectFactors,shape=input$selectFactors2)) + geom_point(size=5) 
		
		}
	p=p+xlab(paste0("PC1: ",percentVar[1],"% variance")) 
	p=p+ylab(paste0("PC2: ",percentVar[2],"% variance")) 
	p=p+ggtitle("Principal component analysis (PCA)")+coord_fixed(ratio=1.0)+ 
     theme(plot.title = element_text(size = 16,hjust = 0.5)) + theme(aspect.ratio=1) +
	 theme(axis.text.x = element_text( size = 16),
	       axis.text.y = element_text( size = 16),
		   axis.title.x = element_text( size = 16),
		   axis.title.y = element_text( size = 16) ) +
	theme(legend.text=element_text(size=16))
	   print(p)
	   }
	# variance chart
	# plot(pca.object,type="bar", xlab="Principal Components", main ="Variances explained")
	
	# pathways
	if(input$PCA_MDS ==2) {  
	withProgress(message="Running pathway analysis", {
	library(PGSEA)
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
	

	if(is.null(input$file2)) { 
	p=ggplot(pcaData, aes(x1, x2, color=Type, shape = Type)) + geom_point(size=5) 
	} else {
		pcaData = cbind(pcaData,readSampleInfo() )
		p=ggplot(pcaData, aes_string("x1", "x2", color=input$selectFactors,shape=input$selectFactors2)) + geom_point(size=5) 
		}
	p=p+xlab("Dimension 1") 
	p=p+ylab("Dimension 2") 
	p=p+ggtitle("Multidimensional scaling (MDS)")+ coord_fixed(ratio=1.)+ 
     theme(plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio=1) +
	 	 theme(axis.text.x = element_text( size = 16),
	       axis.text.y = element_text( size = 16),
		   axis.title.x = element_text( size = 16),
		   axis.title.y = element_text( size = 16) ) +
	theme(legend.text=element_text(size=16))
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

################################################################
#   K-means
################################################################
  
	Kmeans <- reactive({ # Kmeans clustering
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  {  
			tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform 
		}
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) { 
			tem = input$transform; tem = input$logStart; tem= input$lowFilter 
		}
	####################################
	
	withProgress(message="k-means clustering", {
    x <- convertedData()
	#x <- readData()
	#par(mfrow=c(1,2))
	n=input$nGenesKNN
	#if(n>6000) n = 6000 # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
	#x1 <- x;
	#x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	#x = 100* x[1:n,] / apply(x[1:n,],1,sum) 
	x = 100* x[1:n,] / apply(x[1:n,],1,function(y) sum(abs(y))) # L1 norm
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
		incProgress(1, detail = paste("Done")) }) #progress 
	#myheatmap2(x-apply(x,1,mean), bar,1000)
	return(list( x = x, bar = bar)) 

  } )
  
	output$KmeansHeatmap <- renderPlot({ # Kmeans clustering
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$heatColors1
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	####################################
	
	if( is.null(Kmeans()) ) return(NULL)
	withProgress(message="Creating heatmap", {
   
	myheatmap2(Kmeans()$x-apply(Kmeans()$x,1,mean), Kmeans()$bar,1000,mycolor=input$heatColors1)
	
	incProgress(1, detail = paste("Done")) }) #progress 
  } , height = 500)
  
	output$KmeansNclusters <- renderPlot({ # Kmeans clustering
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	withProgress(message="k-means clustering", {
    x <- convertedData()
	#x <- readData()
	#par(mfrow=c(1,2))
	n=input$nGenesKNN
	#if(n>6000) n = 6000 # max
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
		par(mar=c(4,5,4,4))
	plot(1:k, wss, type="b", xlab="Number of Clusters (k)",
		 ylab="Within groups sum of squares",
		 cex=2,cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2	,xaxt="n"	 )
	axis(1, at = seq(1, 30, by = 2),cex.axis=1.5,cex=1.5)
	
	incProgress(1, detail = paste("Done")) }) #progress 
  } , height = 500, width = 550)
  
	KmeansData <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	if( is.null(Kmeans()) ) return(NULL)

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	####################################
	
	#myheatmap2(x, bar)
	Cluster <- toupper(letters)[Kmeans()$bar]
	x <- cbind(Cluster,Kmeans()$x)
	
		# add gene symbol
	if( input$selectOrg != "NEW") 
	{ ix <- match( rownames(x), allGeneInfo()[,1])
	  x <- cbind(as.character( allGeneInfo()$symbol)[ix],x) }
	return(x)
	
	 #progress 
  })
  
	output$downloadDataKmeans <- downloadHandler(
		filename = function() {"Kmeans.csv"},
			content = function(file) {
      write.csv(KmeansData(), file)
	    }
	)
	output$KmeansGO <- renderTable({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectGO3
	if( is.null(input$selectGO3 ) ) return (NULL)
	if( input$selectGO3 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.") )#No matching species
   	if( is.null(Kmeans()) ) return(NULL)
	if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	####################################
	withProgress(message="GO Enrichment", {
		# GO
		pp=0
		minFDR = 0.01

		for( i in 1:input$nClusters) {
			incProgress(1/input$nClusters, , detail = paste("Cluster",toupper(letters)[i]) )

			query = rownames(Kmeans()$x)[which(Kmeans()$bar == i)]
			if(input$selectOrg == "NEW" && !is.null( input$gmtFile) ){ 
				result <- findOverlapGMT( query, GeneSets(),1) 
			} else {
				convertedID <- converted()
				convertedID$IDs <- query
				result = FindOverlap (convertedID,allGeneInfo(),input$selectGO3,input$selectOrg,1) 
			}
			if( dim(result)[2] ==1) next;   # result could be NULL
			result$Genes = toupper(letters)[i] 
			if (pp==0 ) { results <- result; pp <- 1;
			} else {
				results <- rbind(results,result)
			}
		}

		if(pp == 0) return( as.data.frame("No enrichment found."))
		results= results[,c(5,1,2,4)]
		colnames(results)= c("Cluster","FDR","Genes","Pathways")
		if(min(results$FDR) > minFDR ) results = as.data.frame("No signficant enrichment found.") else
		results = results[which(results$FDR < minFDR),]
		incProgress(1, detail = paste("Done")) 
	}) #progress

	if( dim(results)[2] ==1)  return ( as.matrix("No significant enrichment.") )
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
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	####################################

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



################################################################
#   Differential gene expression
################################################################
	output$listFactorsDE <- renderUI({
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	 # Note that " in HTML needs to be denoted with \"    with  the escape character \.
    if (is.null(input$file2) ) {# if sample info is uploaded and correctly parsed.
        return(HTML("<font size = \"3\">A <a href=\"https://idepsite.wordpress.com/data-format/\">sample information file</a> 
	     can be uploaded to build a linear model according to experiment design. </font>"))          
	} else { 	   
		factors = colnames(readSampleInfo())
		choices = setNames(factors, factors  )
		interactions = apply(t(combn(factors,2)),1, function(x) paste(x,collapse=":"))
		choices = append( choices,setNames( interactions, paste(interactions,"interaction") ))
							  
		checkboxGroupInput("selectFactorsModel", 
                              h4("1. Select one or two main factors for the model"), 
                              choices = choices,
                              selected = NULL)	   
	  
	        } 
	})  

	output$listBlockFactorsDE <- renderUI({
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	 # Note that " in HTML needs to be denoted with \"    with  the escape character \.
    if (is.null(input$file2) ) {# if sample info is uploaded and correctly parsed.
		return(NULL)		   
		} else { 
	   
		factors = colnames(readSampleInfo())
		choices = setNames(factors, factors  )
		checkboxGroupInput("selectBlockFactorsModel", 
                              h4("Select a factor for batch effect or paired samples, if needed. (Slow)"), 
                              choices = choices,
                              selected = NULL)  
	  
	        } 
	})	
	
	output$listModelComparisons <- renderUI({
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		
		if (is.null(input$file2) | is.null(input$selectFactorsModel) ) # if sample info is uploaded and correctly parsed.
		   { # if using sample names
		   
		   	factors = as.character ( detectGroups( colnames( readData()$data ) ) )
			factors = unique(factors)# order is reversed
			comparisons = apply(t(combn(factors,2)),1, function(x) paste(x,collapse=" vs. "))
			comparisons = c(comparisons, apply(t(combn(rev(factors),2)),1, function(x) paste(x,collapse=" vs. ")) )	
			comparisons = sort(comparisons)
			choices =  setNames(gsub(" vs\\. ","-",comparisons), comparisons )
			checkboxGroupInput("selectModelComprions", 
									  h4("Select comparisons among sample groups:"), 
									  choices = choices,
									  selected = choices[[1]])	
	   
		   }	 else { 
				choices = list()
				for( selectedFactors in input$selectFactorsModel) { 
					ix = match(selectedFactors, colnames(readSampleInfo() ) )
					if(is.na(ix) ) next;   # if column not found, skip
					factors = unique( readSampleInfo()[,ix])
					comparisons = apply(t(combn(factors,2)),1, function(x) paste(x,collapse=" vs. "))
					comparisons = c(comparisons, apply(t(combn(rev(factors),2)),1, function(x) paste(x,collapse=" vs. ")) )	
					comparisons = sort(comparisons)
					comparisons = paste0(selectedFactors,": ",comparisons)
					choices = append( choices, setNames(comparisons, comparisons ))
				}
				if(length(choices)==0 ) return(NULL) else
				checkboxGroupInput("selectModelComprions", 
									  h4("2. Select one or more comparisons:"), 
									  choices = choices,
									  selected = choices[[1]])	   
				} 
	}) 

	limma <- reactive({  
	if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$CountsDEGMethod; tem = input$countsLogStart
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem=input$CountsDEGMethod
	tem = input$submitModelButton  # this is used to make it responsive only when model and comparisons are completed
	####################################
	
	isolate({ 
	withProgress(message="Identifying Differentially expressed genes", {
	if(input$dataFileFormat == 1 ) {  # if count data
		 if(input$CountsDEGMethod == 3 )   # if DESeq2 method
		  # rawCounts = read.csv("exampleData/airway_GSE52778.csv", row.names=1)
		 # res =DEG.DESeq2(rawCounts, .05, 2) 
		  # res1 =DEG.limma(rawCounts, .1, 1.5,rawCounts, 2,3) 
			return( DEG.DESeq2(readData()$rawCounts,input$limmaPval, input$limmaFC, input$selectModelComprions, readSampleInfo(),input$selectFactorsModel, input$selectBlockFactorsModel)  )
		if(input$CountsDEGMethod < 3 )    # voom or limma-trend
			return( DEG.limma(convertedData(), input$limmaPval, input$limmaFC,readData()$rawCounts, input$CountsDEGMethod,priorCounts=input$countsLogStart,input$dataFileFormat, input$selectModelComprions, readSampleInfo(),input$selectFactorsModel, input$selectBlockFactorsModel) )
	} else { # normalized data
	 return( DEG.limma(convertedData(), input$limmaPval, input$limmaFC,readData()$rawCounts, input$CountsDEGMethod,priorCounts=input$countsLogStart,input$dataFileFormat, input$selectModelComprions, readSampleInfo(),input$selectFactorsModel, input$selectBlockFactorsModel) )
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
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem=input$CountsDEGMethod
	tem = input$selectFactorsModel # responsive to changes in model and comparisons
	tem = input$selectModelComprions
	if(is.null(input$selectComparisonsVenn) ) return(NULL)
	####################################
	
	isolate({ 
	
		results = limma()$results
		results = results[,input$selectComparisonsVenn,drop=FALSE] # only use selected comparisons
		if(dim(results)[2] >5) results <- results[,1:5]
		vennDiagram(results,circle.col=rainbow(5))

	})
    }, height = 600, width = 600)

	output$listComparisonsVenn <- renderUI({
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectContrast", label = NULL, # h6("Funtional Category"), 
                  choices = list("All" = "All"), selected = "All")  
		}	 else { 
				choices = setNames(limma()$comparisons, limma()$comparisons  )
				checkboxGroupInput("selectComparisonsVenn", 
									  h4("Select up to 5 comparisons"), 
									  choices = choices,
									  selected = choices)	

	     } 
	})

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

	output$listComparisonsGenome <- renderUI({
	tem = input$selectOrg

      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectContrast1", label = NULL, # h6("Funtional Category"), 
                  choices = list("All" = "All"), selected = "All")  }	 else { 
	  selectInput("selectContrast2", label="Select a comparison to analyze:",choices=limma()$comparisons
	     )   } 
	})
	
	DEG.data <- reactive({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC; 
		
			tem = input$CountsDEGMethod; 	
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
		####################################
		
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
	tem = input$heatColors1
	tem = input$CountsDEGMethod; 	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	####################################
	
	withProgress(message="Generating heatmap", {
	if( is.null(input$selectContrast)) return(NULL)
	if( is.null(limma()$results)) return(NULL)
	if( is.null(selectedHeatmap.data() )) return(NULL) 
	isolate({ 
	
		 bar = selectedHeatmap.data()$bar +2;
		 bar[bar==3] =2
	  	 myheatmap2( selectedHeatmap.data()$genes,bar,200,mycolor=input$heatColors1,c("Down","Up") )
	 
	 incProgress(1, detail = paste("Done")) 
	
	 })
	})
	
    }, height = 400, width = 500)

	selectedHeatmap.data <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast;
	tem = input$heatColors1
	tem = input$CountsDEGMethod; 	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	####################################
	

	if( is.null(input$selectContrast)) return(NULL)
	if( is.null(limma()$results)) return(NULL)
	#if( is.null(selectedHeatmap.data() )) return(NULL) 
	isolate({ 
		  genes <- limma()$results
		  if( is.null(genes) ) return(NULL)
		  if(!grepl("Diff:", input$selectContrast) ) {  # if not interaction term
			ix = match(input$selectContrast, colnames(genes)) 
		  } else {
			#mismatch in comparison names for interaction terms for DESeq2
			#Diff:water_Wet.genetic_Hy 	 in the selected Contrast
			#Diff-water_Wet-genetic_Hy   in column names
			tem = gsub("Diff-","Diff:" ,colnames(genes))
			tem = gsub("-","\\.",tem)
			ix = match(input$selectContrast, tem) 
			
			if(is.na(ix) ) # this is for limma package
				ix = match(input$selectContrast, colnames(genes)) 			
			
		  }
	  
		  if(is.null(ix)) return(NULL)
		  if(is.na(ix)) return(NULL)
		  if( sum(abs(genes[,ix] )  ) <= 1 ) return(NULL) # no significant genes for this comparison
		  if(dim(genes)[2] < ix ) return(NULL)
		  query = rownames(genes)[which(genes[,ix] != 0)]
		  if(length(query) == 0) return(NULL)
		  iy = match(query, rownames(convertedData()  ) )
		  
		  # find sample related to the comparison
		 iz= match( detectGroups(colnames(convertedData())), unlist(strsplit( input$selectContrast, "-"))	  )
		 iz = which(!is.na(iz))		 
		 if ( !is.null(input$file2) & !is.null(input$selectFactorsModel) & length(input$selectModelComprions)>0 ) {
			comparisons = gsub(".*: ","",input$selectModelComprions)   # strings like: "groups: mutant vs. control"
			comparisons = gsub(" vs\\. ","-",comparisons)		
			factorsVector= gsub(":.*","",input$selectModelComprions) # corresponding factors
			ik = match( input$selectContrast, comparisons )   # selected contrast lookes like: "mutant-control"
			if (is.na(ik)) iz=1:(dim(convertedData())[2])  else {  # interaction term, use all samples		
				selectedfactor= factorsVector[ ik ] # corresponding factors
				iz= match( readSampleInfo()[,selectedfactor], unlist(strsplit( input$selectContrast, "-"))	  )
				iz = which(!is.na(iz))
			}		 
		 }

		 if (grepl("Diff:",input$selectContrast)) iz=1:(dim(convertedData())[2]) # if it is factor design use all samples
		 if( is.na(iz)[1] | length(iz)<=1 ) iz=1:(dim(convertedData())[2]) 

		# color bar
		 bar = genes[,ix]
		 bar = bar[bar!=0]

		 # retreive related data		 
		 genes = convertedData()[iy,iz]
		 
		 genes = genes[order(bar),] # needs to be sorted because myheatmap2 does not reorder genes
		 bar = sort(bar)


		 return(list(genes=genes, bar=bar ))

	
	 })
	})
	
	output$download.selectedHeatmap.data <- downloadHandler(
		filename = function() {paste("Diff_genes_heatmap_",input$selectContrast,".csv",sep="")},
			content = function(file) {
			write.csv(geneListDataExport(), file, row.names=FALSE)
	    }
	)
	
	output$download.DEG.data <- downloadHandler(
		filename = function() {"Diff_expression_all_comparisons.csv"},
		content = function(file) {
			write.csv(DEG.data(), file,row.names=FALSE)
	    }
	)

	# Top DEGs  
	output$geneList <- renderTable({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	
	tem = input$CountsDEGMethod; 	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	####################################
	
	noSig = as.data.frame("No significant genes find!")
	if( is.null(input$selectContrast) ) return(NULL)
	if( is.null( limma()$comparisons ) ) return(NULL) # if no significant genes found
	if( length(limma()$topGenes) == 0 ) return(noSig)
	if( is.null( geneListData() ) ) return(NULL)
	if( dim(geneListData() )[1]> 50 )
		return( geneListData()[1:50,] )
	else return( geneListData() )
	
	
  },digits=2,align="l",include.rownames=F,striped=TRUE,bordered = TRUE, width = "auto",hover=T)
	#, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T,include.rownames=TRUE)
  
	geneListDataExport <- reactive({
		if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
			tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts; tem= input$NminSamples; tem = input$lowFilter; tem=input$transform; tem = input$logStart
		noSig = as.data.frame("No significant genes find!")
		if( is.null(input$selectContrast) ) return(NULL)
		if( is.null( limma()$comparisons ) ) return(NULL) # if no significant genes found
		if( length(limma()$topGenes) == 0 ) return(noSig)
		if( is.null( geneListData() ) ) return(NULL)
		if( input$selectOrg == "NEW" )
			tem <- merge(geneListData(), convertedData(), by.x = 'Top_Genes',by.y = 'row.names') else
		tem <- merge(geneListData(), convertedData(), by.x = 'Ensembl ID',by.y = 'row.names') 
		tem <- tem[order( -sign(tem[,2] ), -abs(tem[,2])),]
		tem$Regulation = "Up"
		tem$Regulation[which(tem[,2]<0 )] <- "Down"
		tem <- tem[,c(dim(tem)[2],1:( dim(tem)[2]-1) )  ]
		return( tem )
	
  })
	#, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T,include.rownames=TRUE)

	geneListData <- reactive({
		if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
			tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem=input$transform; tem = input$logStart
		tem = input$selectFactorsModel # responsive to changes in model and comparisons
		tem = input$selectModelComprions
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
		  if( input$selectGO2 == "ID not recognized!" | input$selectOrg == "NEW") return (top1); 
		  
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
		  top1$FDR <- sprintf("%-3.2e",top1$FDR )
		  colnames(top1) <- c("Ensembl ID", "log2 Fold Change", "Adj.Pval", "Symbol","Chr","Type")
		  if ( sum( is.na(top1$Symbol)) == dim(top1)[1] ) top1 <- top1[,-4] 

		  return(top1)
	
  })

	output$volcanoPlot <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts; tem= input$NminSamples;tem = input$lowFilter; tem=input$transform; tem = input$logStart
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
	 par(mar=c(5,5,1,1))
	 plot(top1$Fold,-log10(top1$FDR),col = c("grey30", "red","blue")[top1$upOrDown],
	 pch =16, cex = .3, xlab= "log2 fold change", ylab = "- log10 (FDR)",
	 cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2	)    
     legend("bottomright",c("Upregulated","Downregulated"),fill = c("red","blue"),cex=1.1 )
	 
	}) })

  },height=450, width=500)
	output$volcanoPlotly <- renderPlotly({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts; tem= input$NminSamples;tem = input$lowFilter; tem=input$transform; tem = input$logStart
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
	 
	 top1$FDR = -log10(top1$FDR)
	 
	 if(0) { 
	 par(mar=c(5,5,1,1))
	 plot(top1$Fold,-log10(top1$FDR),col = c("grey30", "red","blue")[top1$upOrDown],
	 pch =16, cex = .3, xlab= "log2 fold change", ylab = "- log10 (FDR)",
	 cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2	)    
     legend("bottomright",c("Upregulated","Downregulated"),fill = c("red","blue"),cex=1.1 )
	 }
	 
	 # remove non-DEGs
	top1 = top1[which(top1$upOrDown != 1),]  
	top1$type ="Upregulated";
	top1$type[which(top1$upOrDown ==3)] <- "Downregulated"	
	
	# adding gene symbol
	ix <- match( rownames(top1), allGeneInfo()[,1])
	geneSymbols <- as.character( allGeneInfo()$symbol)[ix]
	# if missing or duplicated, use Ensembl ID
	ix <- which(nchar( geneSymbols) <=1 | duplicated(geneSymbols ) );	geneSymbols[ ix ] <- rownames(top1)[ix]
	top1$Symbols = geneSymbols;
	
    p <- plot_ly(data= top1, x = ~Fold, y= ~FDR, color = ~type, text = ~Symbols)%>% 
	layout(xaxis = list(size =35,title = "log2 fold change"), 
           yaxis = list(size =35, title = "- log10 (FDR)"),
		   showlegend = FALSE)
		   
	ggplotly(p)
	}) 
	})

  })

	output$scatterPlot <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem=input$transform; tem = input$logStart
	if( is.null(input$selectContrast) ) return(NULL)
	if( is.null( limma()$comparisons ) ) return(NULL) # if no significant genes found
	if( length(limma()$topGenes) == 0 ) return(NULL)
	if(grepl("Diff:", input$selectContrast) ) 	return(NULL) # if interaction term related comparison
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
	# find sample using sample info file 
	if ( !is.null(input$file2) & !is.null(input$selectFactorsModel) & length(input$selectModelComprions)>0 ) {
		comparisons = gsub(".*: ","",input$selectModelComprions)   # strings like: "groups: mutant vs. control"
		comparisons = gsub(" vs\\. ","-",comparisons)		
		factorsVector= gsub(":.*","",input$selectModelComprions) # corresponding factors
		ik = match( input$selectContrast, comparisons )   # selected contrast lookes like: "mutant-control"
		selectedFactor= factorsVector[ ik ] # corresponding factors
		cat(selectedFactor)
		iz= match( readSampleInfo()[,selectedFactor], unlist(strsplit( input$selectContrast, "-"))[1]	  )
		iz = which(!is.na(iz))
	}
	 
	 
	 genes <- convertedData()[,iz]
	 genes <- as.data.frame(genes)
	 genes$average1 <- apply( genes,1,mean)

     iz= match( detectGroups(colnames(convertedData())), samples[2]	  )
     iz = which(!is.na(iz))
	# find sample using sample info file 
	if ( !is.null(input$file2) & !is.null(input$selectFactorsModel) & length(input$selectModelComprions)>0 ) {
		comparisons = gsub(".*: ","",input$selectModelComprions)   # strings like: "groups: mutant vs. control"
		comparisons = gsub(" vs\\. ","-",comparisons)		
		factorsVector= gsub(":.*","",input$selectModelComprions) # corresponding factors
		ik = match( input$selectContrast, comparisons )   # selected contrast lookes like: "mutant-control"
		selectedfactor= factorsVector[ ik ] # corresponding factors
		iz= match( readSampleInfo()[,selectedfactor], unlist(strsplit( input$selectContrast, "-"))[2]	  )
		iz = which(!is.na(iz))
	}
	 # genes <- cbind( genes,convertedData()[,iz] )
	 genes$average2 <- apply(convertedData()[,iz] ,1,mean)	 
	genes <- merge(genes,top1,by="row.names")
 # write.csv(genes, "tem.csv")
 	 par(mar=c(5,5,1,1))
	plot(genes$average2,genes$average1,col = c("grey45","red","blue")[genes$upOrDown],
	 pch =16, cex = .3, xlab= paste("Average expression in", samples[2] ), 
	 ylab = paste("Average expression in", samples[1] ),
	 cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2	)    
     legend("bottomright",c("Upregulated","Downregulated"),fill = c("red","blue"),cex=1.3 )

		})

	})
	 
  },height=450, width=500)

	output$scatterPlotly <- renderPlotly({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem=input$transform; tem = input$logStart
	if( is.null(input$selectContrast) ) return(NULL)
	if( is.null( limma()$comparisons ) ) return(NULL) # if no significant genes found
	if( length(limma()$topGenes) == 0 ) return(NULL)
	if(grepl("Diff:", input$selectContrast) ) 	return(NULL) # if interaction term related comparison
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
	# find sample using sample info file 
	if ( !is.null(input$file2) & !is.null(input$selectFactorsModel) ) {
		comparisons = gsub(".*: ","",input$selectModelComprions)   # strings like: "groups: mutant vs. control"
		comparisons = gsub(" vs\\. ","-",comparisons)		
		factorsVector= gsub(":.*","",input$selectModelComprions) # corresponding factors
		ik = match( input$selectContrast, comparisons )   # selected contrast lookes like: "mutant-control"
		selectedFactor= factorsVector[ ik ] # corresponding factors
		cat(selectedFactor)
		iz= match( readSampleInfo()[,selectedFactor], unlist(strsplit( input$selectContrast, "-"))[1]	  )
		iz = which(!is.na(iz))
	}
	 
	 genes <- convertedData()[,iz]
	 genes <- as.data.frame(genes)
	 genes$average1 <- apply( genes,1,mean)

     iz= match( detectGroups(colnames(convertedData())), samples[2]	  )
     iz = which(!is.na(iz))

	# find sample using sample info file 
	if ( !is.null(input$file2) & !is.null(input$selectFactorsModel) ) {
		comparisons = gsub(".*: ","",input$selectModelComprions)   # strings like: "groups: mutant vs. control"
		comparisons = gsub(" vs\\. ","-",comparisons)		
		factorsVector= gsub(":.*","",input$selectModelComprions) # corresponding factors
		ik = match( input$selectContrast, comparisons )   # selected contrast lookes like: "mutant-control"
		selectedfactor= factorsVector[ ik ] # corresponding factors
		iz= match( readSampleInfo()[,selectedfactor], unlist(strsplit( input$selectContrast, "-"))[2]	  )
		iz = which(!is.na(iz))
	}
	 # genes <- cbind( genes,convertedData()[,iz] )
	 genes$average2 <- apply(convertedData()[,iz] ,1,mean)	 
	genes <- merge(genes,top1,by="row.names")
    
	# remove non-DEGs
	genes = genes[which(genes$upOrDown != 1),]  
	genes$type ="Upregulated";
	genes$type[which(genes$upOrDown ==3)] <- "Downregulated"	
	
	# adding gene symbol
	ix <- match( genes[,1], allGeneInfo()[,1])
	geneSymbols <- as.character( allGeneInfo()$symbol)[ix]
	# if missing or duplicated, use Ensembl ID
	ix <- which(nchar( geneSymbols) <=1 | duplicated(geneSymbols ) );	geneSymbols[ ix ] <- genes[ix,1]
	genes$Symbols = geneSymbols;
	
    p <- plot_ly(data= genes, x = ~average2, y= ~average1, color = ~type, text = ~Symbols)%>% 
	layout(xaxis = list(size =35,title = paste("Average expression in", samples[2] )), 
           yaxis = list(size =35, title = paste("Average expression in", samples[1] )),
		   showlegend = FALSE)
		   
	ggplotly(p)

		})

	})
	 
  })

	output$geneListGO <- renderTable({		
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if( is.null(input$selectContrast)) return(NULL)
		if( is.null( input$selectGO2) ) return (NULL)
		if( input$selectGO2 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.")) #No matching species

		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO2
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem=input$transform; tem = input$logStart
		if( is.null(limma()$results) ) return(NULL)
		if( is.null(selectedHeatmap.data()) ) return(NULL) # this has to be outside of isolate() !!!
		if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
		NoSig = as.data.frame("No significant enrichment found.")
		isolate({
		withProgress(message="GO Enrichment", {
			
		# using expression data
		genes <- selectedHeatmap.data()$genes
		if(is.null(genes) ) return(NULL) 
		if(dim(genes)[1] <= minGenesEnrichment ) return(NoSig) # if has only few genes
		
		fc = selectedHeatmap.data()$bar		
		# GO
		results1 <- NULL; result <- NULL
		pp <- 0
		for( i in c(1,-1) ) {
			incProgress(1/2 )
			
			if( length(which(fc*i<0)) <= minGenesEnrichment) next; 
			query = rownames(genes)[which(fc*i<0)]
			if( length(query) <= minGenesEnrichment) next; 	
			
			if(input$selectOrg == "NEW" && !is.null( input$gmtFile) ){
				result <- findOverlapGMT( query, GeneSets(),1) 
			} else  { 
				convertedID <- converted()
				convertedID$IDs <- query
				result = FindOverlap (convertedID,allGeneInfo(), input$selectGO2,input$selectOrg,1) }

			if( dim(result)[2] ==1) next;   # result could be NULL
			if(i == -1) result$Genes = "Up regulated"  else result$Genes = "Down regulated"
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
	tem = input$minCounts; tem= input$NminSamples;tem = input$lowFilter; tem=input$transform; tem = input$logStart
	if( is.null(limma()$results) ) return(NULL)
	if( is.null(selectedHeatmap.data()  ) ) return(NULL)
	#if( !is.data.frame(selectedHeatmap.data()  ) ) return(NULL)
	
  isolate({ 
   	withProgress(message="Promoter analysis", {
	genes <- selectedHeatmap.data()$genes
	if(is.null(genes)) return(NULL)
	# if( !is.data.frame(genes)) return (NULL)
	# cat("\nHere",paste(genes[1,],collapse=" ") )
	if( dim(genes)[1] < minGenesEnrichment ) return (NULL) # skip if less than 5 genes total

	fc = selectedHeatmap.data()$bar
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
		
		if(i == -1) result$List ="Up regulated"  else result$List ="Down regulated" 
		if (pp==0 ) { results1 <- result; pp <- 1 } else  { results1 = rbind(results1,result) }
	}

	incProgress(1, detail = paste("Done")) 
	}) #progress
	
	if( is.null(results1)) {
		as.data.frame("No significant motif enrichment found.")
	} else {
		results1 <- results1[,c(4,1:3,5)]
		tem <- results1[,1]
		results1[ duplicated (results1[,1] ),1 ] <- ""
		results1
	}
  })
  }, digits = -1,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)
  
  
################################################################
#   Pathway analysis
################################################################
 
# this updates geneset categories based on species and file
	output$selectGO1 <- renderUI({
	  tem = input$selectOrg
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectGO", label = NULL, # h6("Funtional Category"), 
                  choices = list("All available gene sets" = "All", "GO Biological Process" = "GOBP","GO Molecular Function" = "GOMF","GO Cellular Component" = "GOCC",
                                "KEGG metabolic pathways" = "KEGG"), selected = "GOBP")  }	 else { 
								
	  selectInput("selectGO", label=NULL,choices=gmtCategory(converted(), convertedData(), input$selectOrg,input$gmtFile)
	     ,selected = "GOBP" )   } 
	})

	output$selectGO2 <- renderUI({
	  tem = input$selectOrg
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectGO2", label = NULL, # h6("Funtional Category"), 
                  choices = list("All available gene sets" = "All", "GO Biological Process" = "GOBP","GO Molecular Function" = "GOMF","GO Cellular Component" = "GOCC",
                                "KEGG metabolic pathways" = "KEGG"), selected = "GOBP")  }	 else { 
								
	  selectInput("selectGO2", label=NULL,choices=gmtCategory(converted(), convertedData(), input$selectOrg,input$gmtFile)
	     ,selected = "GOBP" )   } 
	})
	
	output$selectGO3 <- renderUI({
	  tem = input$selectOrg
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectGO3", label = NULL, # h6("Funtional Category"), 
                  choices = list("All available gene sets" = "All", "GO Biological Process" = "GOBP","GO Molecular Function" = "GOMF","GO Cellular Component" = "GOCC",
                                "KEGG metabolic pathways" = "KEGG"), selected = "GOBP")  }	 else { 
								
	  selectInput("selectGO3", label=NULL,choices=gmtCategory(converted(), convertedData(), input$selectOrg,input$gmtFile)
	     ,selected = "GOBP" )   } 
	})

	output$selectGO4 <- renderUI({
	  tem = input$selectOrg
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectGO4", label = NULL, # h6("Funtional Category"), 
                  choices = list("All available gene sets" = "All", "GO Biological Process" = "GOBP","GO Molecular Function" = "GOMF","GO Cellular Component" = "GOCC",
                                "KEGG metabolic pathways" = "KEGG"), selected = "GOBP")  }	 else { 
								
	  selectInput("selectGO4", label=NULL,choices=gmtCategory(converted(), convertedData(), input$selectOrg,input$gmtFile)
	     ,selected = "GOBP" )   } 
	})

	output$PGSEAplot <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	library(PGSEA)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO;		tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
	
	if(is.null(input$selectGO ) ) return (NULL)
	if(input$selectGO == "ID not recognized!" ) return( NULL)
	isolate({ 
	withProgress(message="Running pathway analysis", {
	myrange = c(input$minSetSize, input$maxSetSize)
	genes = convertedData()
if (is.null(input$selectContrast1 ) ) return(NULL)
	incProgress(1/4,"Retrieving gene sets")
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
	library(PGSEA)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
	
	if(is.null(input$selectGO ) ) return (NULL)
	if(input$selectGO == "ID not recognized!" ) return( NULL)
	isolate({ 
	withProgress(message="Running pathway analysis", {
	myrange = c(input$minSetSize, input$maxSetSize)
	genes = convertedData()
if (is.null(input$selectContrast1 ) ) return(NULL)
	incProgress(1/4,"Retrieving gene sets")
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

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
			
	isolate({ 
	gagePathwayData()
	})
  },digits=0,align="l",include.rownames=FALSE,striped=TRUE,bordered = TRUE, width = "auto",hover=T)
 
	gagePathwayData <- reactive({
	library(gage) # pathway analysis	
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)

	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO
	tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
	
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
	  incProgress(1/4,"Retrieving gene sets")
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
		top1[ which( top1[,3] >0),1 ] <- "Up" #gsub("-"," > ",input$selectContrast1 )
		top1[ which( top1[,3] <0),1 ] <- "Down" # gsub("-"," < ",input$selectContrast1 )
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

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
	
	isolate({ 
	fgseaPathwayData()
	})
  },digits=0,align="l",include.rownames=FALSE,striped=TRUE,bordered = TRUE, width = "auto",hover=T)
 
	fgseaPathwayData <- reactive({
	library(fgsea) # fast GSEA
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO; tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold
	if(is.null(input$selectGO ) ) return (NULL)
	if(input$selectGO == "ID not recognized!" ) return( as.data.frame("Gene ID not recognized." ))

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
	
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
	  incProgress(1/4,"Retrieving gene sets")
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
		top1[ which( as.numeric( top1[,3]) >0),1 ] <- "Up" #gsub("-"," > ",input$selectContrast1 )
		top1[ which( as.numeric( top1[,3]) <0),1 ] <- "Down" #gsub("-"," < ",input$selectContrast1 )
		top1 <- top1[order( top1[,1], -abs(as.numeric( top1[,3]) ) ) ,]
		top1[ duplicated (top1[,1] ),1 ] <- ""	 
	    top1[,3] = as.character( round(as.numeric(top1[,3]),4));
	 return( top1)
	}) })
  })

	ReactomePAPathwayData <- reactive({
	library(ReactomePA) # pathway analysis
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO; tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	
	if(is.null(input$selectGO ) ) return (NULL)
	if(input$selectGO == "ID not recognized!" ) return( as.data.frame("Gene ID not recognized." ))

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
	
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
	  incProgress(1/4,"Retrieving gene sets")
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
	  incProgress(3/4,"Runing enrichment analysis using ReactomePA")
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
		top1[ which( as.numeric( top1[,3]) >0),1 ] <- "Up" #gsub("-"," > ",input$selectContrast1 )
		top1[ which( as.numeric( top1[,3]) <0),1 ] <- "Down" #gsub("-"," < ",input$selectContrast1 )
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

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
	
	isolate({ 
	ReactomePAPathwayData()
	})
  },digits=0,align="l",include.rownames=FALSE,striped=TRUE,bordered = TRUE, width = "auto",hover=T)
  
	PGSEAplot.data <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################

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

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
			

      if (is.null(input$file1)&& input$goButton == 0 | is.null(gagePathwayData()))
	  
       { selectInput("sigPathways", label = NULL, # h6("Funtional Category"), 
                  choices = list("All" = "All"), selected = "All")  }	 else { 
		choices <- "All"  # default, sometimes these methods returns "No significant pathway found"
		if( input$pathwayMethod == 1) { if(!is.null(gagePathwayData())) if(dim(gagePathwayData())[2] >1) choices <- gagePathwayData()[,2] }
		else if( input$pathwayMethod == 3) { if(!is.null(fgseaPathwayData())) if(dim(fgseaPathwayData())[2] >1) choices <- fgseaPathwayData()[,2] }
		else if( input$pathwayMethod == 5) { if(!is.null(ReactomePAPathwayData())) if(dim(ReactomePAPathwayData())[2] >1) choices <- ReactomePAPathwayData()[,2] }
		selectInput("sigPathways", label="Select a pathway to show expression pattern of related genes:",choices=choices)
	        } 
	})
	
	selectedPathwayData <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg
	tem = input$sigPathways; 
	if(is.null(gagePathwayData() ) ) return(NULL)
	if(is.null( input$sigPathways))  return (NULL) 
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
	
	#cat(input$sigPathways)
    isolate({ 
	withProgress (message="Retrieving genes",{ 
    if(input$sigPathways == "All") return (NULL) 
	ix <- which(names(GeneSets() ) == input$sigPathways   ) # find the gene set
	if(length(ix) == 0 ) return(NULL)
    genes <- GeneSets()[[ix]] # retrieve genes
	# genes <- pathwayGenes(input$sigPathways,converted(), input$selectGO1,input$selectOrg )
	incProgress(1/2,"Merging data")
	
    x <-  convertedData()[which(rownames(convertedData()) %in% genes) ,]
	if( input$selectOrg != "NEW") {
	ix = match( rownames(x), allGeneInfo()[,1])
	if( sum( is.na(allGeneInfo()$symbol )) != dim(allGeneInfo() )[1] )  # symbol really exists? 
	   rownames(x) <- paste(rownames(x),":", as.character( allGeneInfo()$symbol)[ix]) }
	
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
  
	output$selectedPathwayHeatmap <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg
	tem = input$sigPathways; 
	if(is.null(gagePathwayData() ) ) return(NULL)
	if(is.null( input$sigPathways))  return (NULL) 
	if( is.null(selectedPathwayData()) ) return(NULL)
	# cat(input$sigPathways)
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$heatColors1
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
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
	lhei = c(.5,.2,8)
	
	
	if( dim(x)[1]>200) 
	heatmap.2(x, distfun = dist2,hclustfun=hclust2,
	 col=heatColors[as.integer(input$heatColors1),], density.info="none", trace="none", scale="none", keysize=.5
	,key=T, symkey=F
	,dendrogram = "row"
	,ColSideColors=mycolors[ groups]
	,labRow=""
	,margins=c(10,22)
	,srtCol=45
	,lmat = lmat, lwid = lwid, lhei = lhei
	#,main ="Title"
	)

	if( dim(x)[1]<=200) 
	heatmap.2(x, distfun = dist2,hclustfun=hclust2,
	 col=heatColors[as.integer(input$heatColors1),], density.info="none", trace="none", scale="none", keysize=.5
	,key=T, symkey=F, 
	dendrogram = "row",
	,labRow=gsub(".*:","",rownames(x))
	,ColSideColors=mycolors[ groups]
	,margins=c(10,22)
	,cexRow=1.2
	,srtCol=45
	,lmat = lmat, lwid = lwid, lhei = lhei
	#,main ="Title"
	)
	incProgress(1,"Done")
    }) 
	})
}, height = 1800, width = 600)

	output$KeggImage <- renderImage({
	library(pathview)

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
	#incProgress(1/2, outfile)
	pathID = keggPathwayID(input$sigPathways, Species, "KEGG",input$selectOrg)
	#cat("here5  ",keggSpecies, " ",Species," ",input$sigPathways, "pathID:",pathID,"End", fold[1:5],names(fold)[1:5],"\n")
	#cat("pathway:",is.na(input$sigPathways))
	#cat("\n",fold[1:5],"\n",keggSpecies,"\n",pathID)
    if(is.null(pathID) ) return(blank) # kegg pathway id not found.
	randomString <- gsub(".*file","",tempfile()) 
	 tempFolder <- tempdir() # tempFolder = "temp";
	outfile <- paste( tempFolder,"/",pathID,".",randomString,".png",sep="")
	wd = getwd()
	setwd(tempFolder)

	#pv.out <- mypathview(gene.data = fold, pathway.id = pathID, kegg.dir = tempFolder,  out.suffix = randomString, species = keggSpecies, kegg.native=TRUE)
	pv.out <- pathview(gene.data = fold, pathway.id = pathID, out.suffix = randomString, species = keggSpecies, kegg.native=TRUE)
    setwd(wd)
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
       width = "200%",
        height = "200%",
         alt = "KEGG pathway image.")
		}) 
	})
  }, deleteFile = TRUE)



################################################################
#   Chromosome
################################################################
  

# visualizing fold change on chrs. 
	output$genomePlotly <- renderPlotly({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		#if(is.null(genomePlotDataPre() ) ) return(NULL)
		
		tem = input$selectOrg ; 
		tem = input$selectContrast2
		if (is.null(input$selectContrast2 ) ) return(NULL)
		if( input$selectOrg == "NEW") return(NULL)
		if( length(limma()$topGenes) == 0 ) return(NULL)
		
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
		tem = input$CountsDEGMethod;
		tem=input$limmaPvalViz; 
		tem=input$limmaFCViz
		####################################
		
	  isolate({ 
		withProgress(message="Visualzing expression on the genome", {

		
		if(length( limma()$comparisons)  ==1 )  
		{ top1=limma()$topGenes[[1]]  
		} else {
		  top = limma()$topGenes
		  ix = match(input$selectContrast2, names(top))
		  if( is.na(ix)) return (NULL)
		  top1 <- top[[ix]]; 
		  }
		  if(dim(top1)[1] == 0 ) return (NULL)
		  colnames(top1)= c("Fold","FDR")
		  
		 # write.csv(merge(top1,allGeneInfo(), by.x="row.names",by.y="ensembl_gene_id"  ),"tem.csv"  )
		 x <- merge(top1,allGeneInfo(), by.x="row.names",by.y="ensembl_gene_id"  )
		 
		 x <- x[order(x$chromosome_name,x$start_position),]
		 tem = sort( table( x$chromosome_name), decreasing=T) 
		 chromosomes <- names( tem[tem >= 5 ] )  # chromosomes with less than 100 genes are excluded
		 if(length(chromosomes) > 50) chromosomes <- chromosomes[1:50]  # at most 50 chromosomes
		 chromosomes <- chromosomes[ nchar(chromosomes)<=12] # chr. name less than 10 characters
		 chromosomes = chromosomes[order(as.numeric(chromosomes) ) ]
		 # chromosomes = chromosomes[!is.na(as.numeric(chromosomes) ) ]
		 chromosomesNumbers = as.numeric(chromosomes)
		 # convert chr.x to numbers
		  j = max( chromosomesNumbers,na.rm=T) 
		  for( i in 1:length( chromosomes)) {
		   if ( is.na(chromosomesNumbers[i]) ) 
		   { chromosomesNumbers[i] <- j+1; j <- j+1; }
		 }
		  
		 x <- x[which(x$chromosome_name %in% chromosomes   ),]
		 x <- droplevels(x)
		 
		# find the number coding for chromosome 
		 getChrNumber <- function (chrName){
		 return( chromosomesNumbers[ which( chromosomes == chrName)] )
		 }
		  x$chrNum = 1 # numeric coding
		  x$chrNum <- unlist( lapply( x$chromosome_name, getChrNumber) )
		 
		 x$Row.names <- as.character( x$Row.names)
	 
		 # if symbol is missing use Ensembl id
		 x$symbol = as.character(x$symbol)	 
		 ix = which(is.na(x$symbol))
		 ix2 = which(nchar(as.character(x$symbol))<= 2 )
		 ix3 = which( duplicated(x$symbol))
		 ix = unique( c(ix,ix2,ix3))
		 x$symbol[ix] <- x$Row.names[ix] 
					 
	 
		 ##################################
		 # plotting
		# x = read.csv("tem_genome.csv")
		x = x[!is.na(x$chromosome_name),]
		x = x[!is.na(x$start_position),]	
		# only keep significant genes
		#x = x[which(x$FDR<input$limmaPval),]
		# x = x[which(abs(x$Fold) > log2( input$limmaFC)),]
		
		ix = which( (x$FDR< as.numeric(input$limmaPvalViz)) & (abs(x$Fold) > as.numeric(input$limmaFCViz) ) )
		
		if (length(ix) == 0) { 
				p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
							 geom_blank() + ggtitle("No genes found.") +
							 theme(axis.title.x=element_blank(),axis.title.y=element_blank())
		} else { 
			x = x[ix,]		
			
			x$start_position = x$start_position/1000000 # Mbp
			chrD = 20 # distance between chrs.
			foldCutoff = 3   # max log2 fold 
			
			x$Fold = x$Fold / sd(x$Fold)  # standardize fold change
			
			x$Fold[which(x$Fold > foldCutoff )] = foldCutoff   # log2fold within -5 to 5
			x$Fold[which(x$Fold <   -1*foldCutoff )] = -1*foldCutoff 
			x$Fold = 4* x$Fold
			

			
			x$y = x$chrNum*chrD + x$Fold
			chrLengthTable = aggregate(start_position~chrNum, data=x,max )
			chrTotal = dim(chrLengthTable)[1]	
			x$R = as.factor(sign(x$Fold))
			
			colnames(x)[ which(colnames(x) == "start_position")] = "x"

			p= ggplot(x, aes(x = x, y = y, colour = R, text = symbol ) ) + geom_point(shape = 20, size = .2)
			
				#label y with chr names
			p <- p +  scale_y_continuous(labels = paste("chr",chromosomes[chrLengthTable$chrNum],sep=""), breaks = chrD* (1:chrTotal), limits = c(0, 
				chrD*chrTotal + 5) )
			# draw horizontal lines for each chr.
			for( i in 1:dim(chrLengthTable)[1] )
				p = p+ annotate( "segment",x = 0, xend = chrLengthTable$start_position[i],
					y = chrLengthTable$chrNum[i]*chrD, yend = chrLengthTable$chrNum[i]*chrD)
			# change legend		http://ggplot2.tidyverse.org/reference/scale_manual.html
			p=p+scale_colour_manual(name="",   # customize legend text
				values=c("blue", "red"),
				breaks=c("1","-1"),
				labels=c("Up", "Dn"))	
			p = p + xlab("Position on chrs. (Mbp)") +	 theme(axis.title.y=element_blank())					 
			p= p + theme(legend.position="none")
			# p <- ggplot(mtcars, aes(x=hp, y=mpg)) + geom_point(shape=20, size=17)
		   # p=p+ geom_smooth(method = "lm",se=FALSE)
			#p+ geom_line(aes(y=rollmean(y,7, na.pad=TRUE)))
		  
		   # Customize hover text https://cran.r-project.org/web/packages/plotly/plotly.pdf
		   # style(ggplotly(p),hoverinfo="text")  # not working
		  }  # if no genes else
			ggplotly(p)
		
			}) # progress
		}) # isloate
	  })

# pre-calculating PREDA, so that changing FDR cutoffs does not trigger entire calculation
	genomePlotDataPre <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)

	tem = input$selectOrg ; 
	tem = input$selectContrast2
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
	
  isolate({ 
	withProgress(message="Identifying differentially expressed genomic regions using PREDA", {
	if (is.null(input$selectContrast2 ) ) return(NULL)
	
	if( length(limma()$topGenes) == 0 ) return(NULL)
	
	if(length( limma()$comparisons)  ==1 )  
    { top1=limma()$topGenes[[1]]  
	} else {
	  top = limma()$topGenes
	  ix = match(input$selectContrast2, names(top))
	  if( is.na(ix)) return (NULL)
	  top1 <- top[[ix]]; 
	  }
	  if(dim(top1)[1] == 0 ) return (NULL)
	  colnames(top1)= c("Fold","FDR")
	  
	 # write.csv(merge(top1,allGeneInfo(), by.x="row.names",by.y="ensembl_gene_id"  ),"tem.csv"  )
	 x <- merge(top1,allGeneInfo(), by.x="row.names",by.y="ensembl_gene_id"  )
	 
	 #hist(top1[,1])
	 
	 ########## PREDA
	infofile <- system.file("sampledata", "GeneExpression", "sampleinfoGE_PREDA.txt", package = "PREDAsampledata")
	sampleinfo<-read.table(infofile, sep="\t", header=TRUE)
	head(sampleinfo)

	data(ExpressionSetRCC)
	GEstatisticsForPREDA<-statisticsForPREDAfromEset(ExpressionSetRCC, statisticType="tstatistic", referenceGroupLabel="normal", classVector=sampleinfo[,"Class"])
	analysesNames(GEstatisticsForPREDA)
	
	GEGenomicAnnotations<-eset2GenomicAnnotations(ExpressionSetRCC, retain.chrs=1:22) # needs hgu133plus2.db package
	GEGenomicAnnotationsForPREDA<-GenomicAnnotations2GenomicAnnotationsForPREDA(GEGenomicAnnotations, reference_position_type="median")
	GEDataForPREDA<-MergeStatisticAnnotations2DataForPREDA(GEstatisticsForPREDA, GEGenomicAnnotationsForPREDA, sortAndCleanNA=TRUE)

	#x = read.csv("PREDA_test.csv")

	#x <- x[,-1] # remove the row index, not needed in Shiny

	 x <- x[order(x$chromosome_name,x$start_position),]
	 tem = sort( table( x$chromosome_name), decreasing=T) 
	 chromosomes <- names( tem[tem > 100 ] )  # chromosomes with less than 100 genes are excluded
	 if(length(chromosomes) > 50) chromosomes <- chromosomes[1:50]  # at most 50 chromosomes

	 chromosomes = chromosomes[order(as.numeric(chromosomes) ) ]
	 # chromosomes = chromosomes[!is.na(as.numeric(chromosomes) ) ]
	 chromosomesNumbers = as.numeric(chromosomes)
	 # convert chr.x to numbers
	  j = max( chromosomesNumbers,na.rm=T) 
	  for( i in 1:length( chromosomes)) {
	   if ( is.na(chromosomesNumbers[i]) ) 
	   { chromosomesNumbers[i] <- j+1; j <- j+1; }
	 }
	  
	 x <- x[which(x$chromosome_name %in% chromosomes   ),]
	 x <- droplevels(x)
	 
	# find the number coding for chromosome 
	 getChrNumber <- function (chrName){
	 return( chromosomesNumbers[ which( chromosomes == chrName)] )
	 }
	  x$chrNum = 1 # numeric coding
	  x$chrNum <- unlist( lapply( x$chromosome_name, getChrNumber) )
	 
	 x$Row.names <- as.character( x$Row.names)
	 fold = x$Fold; names(fold) = x$Row.names
	 
	# write.csv(x,"tem.csv")
	 x <- x[!duplicated(x$Row.names),]
	 #rownames(x) = x$Row.names
	 
	 myData <- GEDataForPREDA
	  
	 myData@position <- as.integer( x$start_position+ x$genomeSpan/2 )
	 myData@ids <- as.character( x$Row.names)
	 myData@chr <- as.integer( x$chrNum )
	 myData@start <- x$start_position
	 myData@end <- x$start_position + x$genomeSpan
	 myData@strand <- rep(1,dim(x)[1])
	 myData@chromosomesLabels <- chromosomes
	 myData@chromosomesNumbers <- as.integer( chromosomesNumbers )
	 myData@statistic <- as.matrix( x[,2,drop=FALSE] )
	 myData@analysesNames <- "Test"
	 myData@testedTail <- "both"
	 myData@optionalAnnotations <- as.matrix( x[,1:2,drop=FALSE] )
	 myData@optionalAnnotationsHeaders <- c("a","b")
	 
	incProgress(1/2, "Runing 1000x permutations using PREDA. This may take up to 5 minutes.")
	set.seed(2)
	GEanalysisResults<-PREDA_main(myData,nperms = PREDA_Permutations)

	return( GEanalysisResults )
	}) # progress
	}) # isloate
  })

# results from PREDA
	genomePlotData <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    if(is.null(genomePlotDataPre() ) ) return(NULL)
	
	tem = input$selectOrg ; 
	tem = input$selectContrast2
	tem = input$StatisticCutoff
	tem = input$RegionsPvalCutoff
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
	
  isolate({ 
	withProgress(message="Identifying differentially expressed genomic regions using PREDA", {
	if (is.null(input$selectContrast2 ) ) return(NULL)
	
	if( length(limma()$topGenes) == 0 ) return(NULL)
	
	if(length( limma()$comparisons)  ==1 )  
    { top1=limma()$topGenes[[1]]  
	} else {
	  top = limma()$topGenes
	  ix = match(input$selectContrast2, names(top))
	  if( is.na(ix)) return (NULL)
	  top1 <- top[[ix]]; 
	  }
	  if(dim(top1)[1] == 0 ) return (NULL)
	  colnames(top1)= c("Fold","FDR")
	  
	 # write.csv(merge(top1,allGeneInfo(), by.x="row.names",by.y="ensembl_gene_id"  ),"tem.csv"  )
	 x <- merge(top1,allGeneInfo(), by.x="row.names",by.y="ensembl_gene_id"  )
	 
	 x <- x[order(x$chromosome_name,x$start_position),]
	 tem = sort( table( x$chromosome_name), decreasing=T) 
	 chromosomes <- names( tem[tem > 100 ] )  # chromosomes with less than 100 genes are excluded
	 if(length(chromosomes) > 50) chromosomes <- chromosomes[1:50]  # at most 50 chromosomes

	 chromosomes = chromosomes[order(as.numeric(chromosomes) ) ]
	 # chromosomes = chromosomes[!is.na(as.numeric(chromosomes) ) ]
	 chromosomesNumbers = as.numeric(chromosomes)
	 # convert chr.x to numbers
	  j = max( chromosomesNumbers,na.rm=T) 
	  for( i in 1:length( chromosomes)) {
	   if ( is.na(chromosomesNumbers[i]) ) 
	   { chromosomesNumbers[i] <- j+1; j <- j+1; }
	 }
	  
	 x <- x[which(x$chromosome_name %in% chromosomes   ),]
	 x <- droplevels(x)
	 
	# find the number coding for chromosome 
	 getChrNumber <- function (chrName){
	 return( chromosomesNumbers[ which( chromosomes == chrName)] )
	 }
	  x$chrNum = 1 # numeric coding
	  x$chrNum <- unlist( lapply( x$chromosome_name, getChrNumber) )
	 
	 x$Row.names <- as.character( x$Row.names)
	 fold = x$Fold; names(fold) = x$Row.names
	 
	# write.csv(x,"tem.csv")
	 x <- x[!duplicated(x$Row.names),]
	 #rownames(x) = x$Row.names
	 
	incProgress(1/2, "Summarizing statistics")
	GEanalysisResults <- genomePlotDataPre(); 
	
	genomic_regions_UP<-PREDAResults2GenomicRegions(GEanalysisResults, qval.threshold=input$RegionsPvalCutoff, smoothStatistic.tail="upper", smoothStatistic.threshold= input$StatisticCutoff)
	genomic_regions_DOWN<-PREDAResults2GenomicRegions(GEanalysisResults, qval.threshold=input$RegionsPvalCutoff, smoothStatistic.tail="lower", smoothStatistic.threshold= -1 *input$StatisticCutoff)

	if( is.null(genomic_regions_UP$Test) &&  is.null(genomic_regions_UP$Test)  ) return(-1) # no significant regions
	
	regions <- 0 
	if( !is.null(genomic_regions_UP$Test)) { 
	dataframe_UPregions<-GenomicRegions2dataframe(genomic_regions_UP[[1]])
	dataframe_UPregions$Regulation <- "Up"
	regions = dataframe_UPregions
	}

	if( !is.null(genomic_regions_DOWN$Test) ){ 
	dataframe_DOWNregions<-GenomicRegions2dataframe(genomic_regions_DOWN[[1]])
	dataframe_DOWNregions$Regulation <- "Down"
	if ( class(regions) != "data.frame" ) regions <- dataframe_DOWNregions else # if UP regions is NULL
	 regions = rbind(regions,dataframe_DOWNregions  )
	}

	if( class(regions) != "data.frame") return(-1)

	Regions <- regions[,c(4,1:3)] 
	Regions <- Regions[which(Regions$end != Regions$start),]
	Regions$size = round((Regions$end - Regions$start )/1000000,3)
	Regions$chr <- chromosomes[ Regions$chr ] # convert from chr. number to names

	# find the gene indices in the up or down regulated regions
	regulatedGenes <- function (i) {
	ix =  which( Regions$chr == x$chromosome_name[i] &
		 Regions$start < x$start_position[i] #&(Regions$end > x$start_position[i] + x$genomeSpan[i])
		 &(Regions$end > x$start_position[i] )) # if the start position is within the region
	if( length(ix) == 0 | length(ix) >1  ) return(NA) else return( ix )
	}

	regionID = unlist( lapply(1:dim(x)[1], regulatedGenes) )
	x1 = x[which(!is.na(regionID)),]
	regionID = regionID[!is.na(regionID)]
	x1 = cbind(regionID,Regions[regionID, ], x1[,c('symbol', 'Row.names', 'Fold','FDR', 'band','start_position')])
	x1 = x1[order(x1$regionID,x1$start_position ), ]	
    colnames(x1)[8] = "Ensembl"; colnames(x1)[4] = "Region.Start";colnames(x1)[5] = "Region.End"; colnames(x1)[12] = "Gene.Start"

	# number of genes
	tem = table(x1$regionID)
	Regions$Ngenes = 0
	Regions$Ngenes[as.integer(names(tem) ) ] <- tem

	# cytoband per region
	tem = unique(x1[,c('regionID','band')])
	tem$band = gsub("\\..*","",tem$band)
	tem = unique(tem)
	Regions$band =""
	Regions$ID <- 1:dim(Regions)[1]
	for ( i in 1:dim(Regions)[1] )
	  Regions$band[i] <- paste( tem[ which( tem[,1] == Regions$ID[i]),2],  collapse=";" )

	# Genes
	Regions$Genes <- ""
	for ( i in 1:dim(Regions)[1] )
	  Regions$Genes[i] <- paste( x1$symbol[ which( x1[,1] == Regions$ID[i])],  collapse=" " )
	return( list(mainResult = GEanalysisResults, Regions = Regions, Genes = x1,legend.x=max(x$start_position)*.6, legend.y=max(chromosomesNumbers)-3  ) )
	}) # progress
	}) # isloate
  })

# Using PREDA to identify significant genomic regions 
	output$genomePlot <- renderPlot({
	library(PREDA)  # showing expression on genome
	library(PREDAsampledata) 
	library(hgu133plus2.db)
	
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    if( is.null( genomePlotData() ) ) return(NULL)
	tem = input$selectOrg ; 
	tem = input$selectContrast2
	tem = input$StatisticCutoff
	tem = input$RegionsPvalCutoff
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################

isolate({ 
  if( class(genomePlotData()) != "list" ) { plot.new(); text(0.2,1, "No significant regions found!")} else {
   GEanalysisResults <- genomePlotData()$mainResult;
     if( is.null( GEanalysisResults ) ) return(NULL)	
	genomic_regions_UP<-PREDAResults2GenomicRegions(GEanalysisResults, qval.threshold=input$RegionsPvalCutoff, smoothStatistic.tail="upper", smoothStatistic.threshold= -1 *input$StatisticCutoff)
	genomic_regions_DOWN<-PREDAResults2GenomicRegions(GEanalysisResults, qval.threshold=input$RegionsPvalCutoff, smoothStatistic.tail="lower", smoothStatistic.threshold= -1 *input$StatisticCutoff)
	 
	checkplot<-genomePlot(GEanalysisResults, genomicRegions=c(genomic_regions_UP, genomic_regions_DOWN), grouping=c(1, 1), scale.positions="Mb", region.colors=c("red","blue"))
	legend(x=genomePlotData()$legend.x, y=genomePlotData()$legend.y , legend=c("UP", "DOWN"), fill=c("red","blue"))
	} #if  else
 	 })
  }, height = 800, width = 1000)
 
	output$downloadRegions <- downloadHandler(
		filename = function() {paste("Diff_Chr_Regions_",input$selectContrast2,".csv",sep="")},
		content = function(file) {
			write.csv(genomePlotData()$Regions, file, row.names=FALSE)
	    }
	)
  
	output$downloadGenesInRegions <- downloadHandler(
		filename = function() {paste("Genes_in_Diff_Chr_Regions_",input$selectContrast2,".csv",sep="")},
		content = function(file) {
		write.csv(genomePlotData()$Genes, file, row.names=FALSE)
	}
	)

	output$chrRegionsList <- renderTable({
  if (is.null(input$file1) && input$goButton == 0)   return(NULL)

  	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
  
  tem = input$selectContrast2
  if( is.null( genomePlotData() ) |class(genomePlotData()) != "list"   ) return(NULL)
      tem = genomePlotData()$Regions
	  tem = tem[,c(8,1:7,9)]
	  colnames(tem)[1] ="RegionID"
	  tem = tem[,-7]
	  
	  tem <- tem[,c(2,3,6,7)]
	  colnames(tem) = c("Dir.","Chr","Mbp","Band")
	  tem$Mbp <- round(tem$Mbp,2)
	  tem

  },rownames= FALSE)

	output$chrRegions <- DT::renderDataTable({
  if (is.null(input$file1) && input$goButton == 0)   return(NULL)
  
 	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
  
  tem = input$selectContrast2
      tem = genomePlotData()$Regions
	  tem = tem[,c(8,1:7,9)]
	  colnames(tem)[1] ="RegionID"
	  tem = tem[,-7]
	 
	  tem

  },rownames= FALSE)

	output$genesInChrRegions <- DT::renderDataTable({
  if (is.null(input$file1) && input$goButton == 0)   return(NULL)
  
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
	tem = input$CountsDEGMethod;
	####################################
  
  tem = input$selectContrast2

      tem = genomePlotData()$Genes
	   tem <- tem[,-c(5,10,12)]
	  tem$Fold <- round(tem$Fold,3)
	  colnames(tem)[2]="Dir"
	  tem


  },rownames= FALSE)  

################################################################
#   Biclustering
################################################################
   biclustering <- reactive({
	  if (is.null(input$file1) && input$goButton == 0)   return(NULL)

		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
		tem = input$CountsDEGMethod;
		
		tem = input$nGenesBiclust
		tem = input$biclustMethod
		####################################   
   
		isolate({
			library(biclust)
			# library(QUBIC) # have trouble installing on Linux
			# library(runibic) # 
			
		
			x <- convertedData()
			n=input$nGenesBiclust	
			if(n>dim(x)[1]) n = dim(x)[1] # max	as data
			if(n<10) n = 10 
			if(n> 2000 ) n = 2000 			
			x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
			write.csv(x,"tem.csv")
			#res <- biclust::biclust(as.matrix( x), method = BCQU()) 
			runR = paste( "res <- biclust::biclust(as.matrix( x), method =", input$biclustMethod ,")" )
			eval(parse(text = runR ) )
			
			
			return(list( x=x, res=res)	)	 
		
		
		} )
   
   
   
   } )
   
  	output$listBiclusters <- renderUI({
		tem = input$selectOrg
		tem = input$biclustMethod
		tem = input$nGenesBiclust
		if (is.null(biclustering() ) ){ # if sample info is uploaded and correctly parsed.
		   return(NULL)	   
		} else { 
			if( biclustering()$res@Number == 0 ) { 
				return(NULL) 
			}	else  
					selectInput(	"selectBicluster", 
									label="Select a cluster",
									choices= 1:biclustering()$res@Number  
								)

				

		}			
	}) 

  	output$biclusterInfo <- renderText({
		tem = input$nGenesBiclust
		tem = input$selectBicluster
		tem = input$biclustMethod
		if(  is.null(input$selectBicluster)) return("No cluster found! Perhaps sample size is too small." )
		if (is.null(biclustering()  ) ){ # if sample info is uploaded and correctly parsed.
		   return("No cluster found! Perhaps sample size is too small." )	   
		} else {
			res = biclustering()$res
			if( res@Number == 0 ) { 
				return( "No cluster found! Perhaps sample size is too small."  ) 
			}	else  { 
					x = biclust::bicluster(biclustering()$x, res, as.numeric( input$selectBicluster)  )[[1]]
					paste(res@Number,"clusters found. Cluster ",input$selectBicluster, " has", dim(x)[1], "genes correlated across", dim(x)[2], "samples." )	
				}
				

		}			
	}) 
  
	output$biclustHeatmap <- renderPlot ({
	  if (is.null(input$file1) && input$goButton == 0)   return(NULL)

		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
		tem = input$CountsDEGMethod;
		tem = input$nGenesBiclust
		tem = input$heatColors1
		tem = input$selectBicluster
		tem = input$biclustMethod
		####################################  
		if( is.null(biclustering() ) ) return(NULL)
		if( is.null(input$selectBicluster ) ) return(NULL)
		tem = input$biclustMethod
		isolate({ 
			res = biclustering()$res
			if( res@Number == 0 ) { plot.new(); text(0.5,0.5, "No cluster found!")} else {
		
				x = biclust::bicluster(biclustering()$x, res, as.numeric( input$selectBicluster)  )[[1]]

					par(mar = c(5, 4, 1.4, 0.2))

				if( dim(x)[1] <=30 )
					heatmap.2(x,  Rowv =T,Colv=F, dendrogram ="none",
					col=heatColors[as.integer(input$heatColors1),], density.info="none", trace="none", scale="none", keysize=.2
					,key=F, #labRow = T,
					,margins = c(8, 24)
					,cexRow=1
					#,srtCol=45
					,cexCol=1.  # size of font for sample names
					) else 		
					heatmap.2(x,  Rowv =T,Colv=F, dendrogram ="none",
					col=heatColors[as.integer(input$heatColors1),], density.info="none", trace="none", scale="none", keysize=.2
					,key=F, labRow = F,
					,margins = c(8, 4)
					,cexRow=1
					#,srtCol=45
					,cexCol=1.  # size of font for sample names

					)				

				
			
			}


		})
   
  }, height = 400, width = 400)

	output$geneListBclustGO <- renderTable({		
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if( is.null( input$selectGO4) ) return (NULL)
		if( input$selectGO4 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.")) #No matching species

		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO4
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem=input$transform; tem = input$logStart

		if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
		NoSig = as.data.frame("No significant enrichment found.")
		
		if( is.null(biclustering() ) ) return(NULL)
		if( is.null(input$selectBicluster ) ) return(NULL)
		tem = input$nGenesBiclust	
		tem = input$biclustMethod
		
		isolate({
			withProgress(message="GO Enrichment", {
			
			res = biclustering()$res
			if( res@Number == 0 ) return(as.data.frame("No clusters found!") ) 
			x = biclust::bicluster(biclustering()$x, res, as.numeric( input$selectBicluster)  )[[1]]
			
			if( dim(x)[1]  <= minGenesEnrichment) return(NoSig) 
			query = rownames(x)
			
			if(input$selectOrg == "NEW" && !is.null( input$gmtFile) ){
				result <- findOverlapGMT( query, GeneSets(),1) 
			} else  { 
				convertedID <- converted()
				convertedID$IDs <- query
				result = FindOverlap (convertedID,allGeneInfo(), input$selectGO4,input$selectOrg,1) }
				
			result$Genes = "Up regulated"
			

			results1 = result
			if ( is.null( results1) ) return (NoSig)
			if( dim(results1)[2] <5 ) return(NoSig)  # Returns a data frame: "No significant results found!"

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

			results1[ duplicated (results1[,1] ),1 ] <- ""
			
			results1[,-1]			


		 })#progress
		}) #isolate
  }, digits = 0,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)
	#   output$selectedHeatmap <- renderPlot({       hist(rnorm(100))    })

	output$geneListBicluster <- renderTable({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if( is.null( input$selectGO4) ) return (NULL)
		if( input$selectGO4 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.")) #No matching species

		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO4
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem=input$transform; tem = input$logStart

		if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
		NoSig = as.data.frame("No significant enrichment found.")
		
		if( is.null(biclustering() ) ) return(NULL)
		if( is.null(input$selectBicluster ) ) return(NULL)
		tem = input$nGenesBiclust
		tem = input$biclustMethod
	
		res = biclustering()$res
		if( res@Number == 0 ) return(as.data.frame("No clusters found!") ) 
		top1 = biclust::bicluster(biclustering()$x, res, as.numeric( input$selectBicluster)  )[[1]]	
		top2 = top1  
		  # if new species
		  if( input$selectGO4 == "ID not recognized!" | input$selectOrg == "NEW") return (top1); 
		  
		  top1 <- merge(top1, allGeneInfo(), by.x ="row.names", by.y="ensembl_gene_id",all.x=T )
		  
		  if ( sum( is.na(top1$band)) == dim(top1)[1] ) top1$chr = top1$chromosome_name else
			top1$chr = paste( top1$chromosome_name, top1$band,sep="")

		  geneMean = apply(top1[, colnames(top2)  ],1,mean )
		  top1 = cbind(top1, geneMean)
		  
		  top1 <- top1[order(geneMean,decreasing=T),] # sort by the 2nd column
		  
		  top1 <- top1[,c('Row.names','geneMean','symbol','chr','gene_biotype')]
		  
			top1$geneMean <- sprintf("%-4.2f",as.numeric(top1$geneMean) )

		  colnames(top1) <- c("Ensembl ID", "Mean","Symbol","Chr","Type")
		  if ( sum( is.na(top1$Symbol)) == dim(top1)[1] ) top1 <- top1[,-3] 
		  
		  if(dim(top1)[1] > 1000 ) top1 = top1[1:1000,] # at most 1000 genes are shown
		  return(top1)
	
  }, digits = 0,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)

  	biclustData <- reactive({
  		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if( is.null( input$selectGO4) ) return (NULL)
		if( input$selectGO4 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.")) #No matching species

		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO4
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem=input$transform; tem = input$logStart

		if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
		
		if( is.null(biclustering() ) ) return(NULL)
		if( is.null(input$selectBicluster ) ) return(NULL)
		tem = input$nGenesBiclust
		tem = input$biclustMethod
		if(biclustering()$res@Number == 0) return(NULL)
	 isolate({ 
		res = biclustering()$res
		
		# given a dataframe, cover it to a string object.
		dataframe2String <- function(x){
			sampleNames = paste0("\t",paste(colnames(x),collapse="\t"))
			data = paste(rownames(x),apply(x,1, paste, collapse="\t"),sep="\t")
			return( paste(c(sampleNames, data),collapse="\n" ) )
		}
		

		x = biclust::bicluster(biclustering()$x, res, 1:res@Number)
		a=paste("Total of ", res@Number, "clusters\n")
		for ( i in 1:res@Number) {
		   x1 = x[[i]]
		   if( input$selectGO4 != "ID not recognized!" & input$selectOrg != "NEW") {
				  x2 <- merge(x1, allGeneInfo(), by.x ="row.names", by.y="ensembl_gene_id",all.x=T )
				 rownames(x2)= paste(x2$symbol,  x2$Row.names )
				 x1 = x2[, 2:(dim(x1)[2]+1) ]
	
		   }
	  
		   a = paste(a, "\n\nCluster", i,":", dim( x1 )[1], "genes", dim(x1)[2], "samples\n")
	       a = paste0(a, dataframe2String( round(x1,4)  ) )
	
		}
		return(a)	
		})
	})

	output$download.biclust.data <- downloadHandler(
		filename = function() {"biclustering_result.txt"},
			content = function(file) {
			write(biclustData(), file)
	    }
	)
  
################################################################
#   Session Info
################################################################
  
	# output user settings and session info
	output$RsessionInfo <- renderUI({
	if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		
	i = "<h4> <a href=\"mailto:Xijin.Ge@SDSTATE.EDU?Subject=iDEP\" target=\"_top\">Email us</a> for questions, 
	suggestions, or data contributions. Stay connected via <a href=\"https://groups.google.com/d/forum/idep\"> user group</a>, or
	   <a href=\"https://twitter.com/useIDEP\"> Twitter. </a></h4>"  
	i = c(i,"<h2>R as in Reproducibility</h2>")
	i = c(i,"We recommend users to save the following details about their analyses.")
	i = c(i,paste("Analysis were conducted using the awesome ",iDEPversion, 
		", hosted at http://ge-lab.org on ", date(),".",sep="") )
		
	
	i = c(i,"If you really, really need it, here is my embarrassingly  
	messy <a href=\"https://github.com/gexijin/iDEP\"> source code.</a>")
	
	#inFile <- input$file1
	#inFile <- inFile$datapath
	#if (is.null(input$file1) && input$goButton == 0)   return(NULL)
	#if (is.null(input$file1) && input$goButton > 0 )   inFile = paste(demoDataFile,"(Demo file)")	  
	#i = c(i,paste("Input file:", inFile ))
	i = c(i, paste("<br><strong>Data</strong><br>Species:",converted()$species[2]))	
	i = c(i, paste("Number of samples:",dim(convertedData() )[2]))	
	i = c(i, paste("Number of genes converted and filtered:",dim(convertedData() )[1]))	
	i = c(i, limma()$Exp.type)	
	
	if(input$dataFileFormat == 2 )  {  # if FPKM, microarray
		i = c(i,"Input file type: normalized experssiopn file")
		i = c(i,"<br><strong>Pre-processing and exploratory data analysis settings:</strong>", paste("Low filter: lowFilter=", input$lowFilter) ) 
		if(input$transform == TRUE | readData()$mean.kurtosis > kurtosis.log) 
			i = c(i, paste("Log-transformed: x' = log2( x +",input$logStart,")")) else 
			i = c(i, "No log-transformation.")
			
		}else { # read counts
			i = c(i,"Input file type: RNA-seq read count file")    
			i = c(i,"<br><strong>Pre-processing and exploratory data analysis settings:</strong>", paste("Min. counts: minCounts=", input$minCounts), paste("Min. counts samples: NminSamples=", input$NminSamples) ) 
			transforms1 = c("Started log: log2(x+c)" ,"VST: variance stabilizing transform", "rlog: regularized log")
			i = c(i, paste("Counts data transformation method: ", transforms1[as.numeric(input$CountsTransform)] ) )
			if ( input$CountsTransform == 1)
				i = c(i,paste("Pseudo count: c=",input$countsLogStart ) )
			DEchoices = c("limma-trend", "limma-voom","DESeq2")
			i = c(i, paste("Method for differential expression:  CountsDEGMethod=", 
				input$CountsDEGMethod, "(",DEchoices[ as.numeric(input$CountsDEGMethod)],")"))
		}

	i = c(i, paste("number of genes in heatmap: nGenes=", input$nGenes))
	i = c(i, paste("number of genes in k-means clustering: nGenesKNN=", input$nGenesKNN))
	i = c(i, paste("number of clusters in k-means clustering: nClusters=", input$nClusters))
	i = c(i, paste("Promoter analysis for k-means clustering: radioPromoterKmeans=", input$radioPromoterKmeans,"bp"))
		
	i = c( i, paste("<br><strong> Differential expression settings:</strong>") )
	i = c(i, paste("FDR cutoff: limmaPval=", input$limmaPval) )
	i = c(i, paste("Fold-change cutoff: limmaFC=", input$limmaFC) )	
	i = c(i, paste("Promoter analysis for DEGs: radio.promoter=", input$radio.promoter,"bp"))

	i = c( i, paste("<br> <strong>Pathway analysis settings:</strong>") )
	PAMethods =c("GAGE", "PGSEA", "GSEA (preranked fgsea)", "PGSEA w/ all samples", "ReactomePA")
	i = c(i, paste("Pathway analysis methods: pathwayMethod=", PAMethods[as.numeric(input$pathwayMethod)]))
	i = c(i, paste("FDR cutoff: pathwayPvalCutoff=", input$pathwayPvalCutoff) )
	i = c(i, paste("Min size for gene set: minSetSize=", input$minSetSize) )
	i = c(i, paste("Max size for gene set: maxSetSize=", input$maxSetSize) )
	if(input$absoluteFold == 1) 
		i = c(i, paste("Using absolute values for fold change: absoluteFold=", input$absoluteFold) )
	
	i = c( i, paste("<br> <strong>PREDA settings:</strong>") )
	i = c(i, paste("FDR cutoff: RegionsPvalCutoff=", input$RegionsPvalCutoff) )
	i = c(i, paste("FDR cutoff: StatisticCutoff=", input$StatisticCutoff) )
	
	i =c(i,"<br><h3>R session info: </h3>")
	i = c(i,capture.output(sessionInfo()) )
	HTML(paste(i, collapse='<br/>') )
  })
  
})  # shiny Server
