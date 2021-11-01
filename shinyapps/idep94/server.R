# iDEP server logic, By Steven Ge Xijin.Ge@sdstate.edu
# Integrated Differential Gene Expression and Pathway analysis
# hosted at http://ge-lab.org/idep/
# manuscript: https://www.biorxiv.org/content/early/2018/04/20/148411 

iDEPversion = "iDEP 0.94"     

################################################################
# R packages
################################################################
# R packages, installed by:
#auto install
#Rlibs = c("shiny","RSQLite","gplots","ggplot2","e1071","shinyAce","shinyBS","reshape2","DT","plotly","statmod","biclust","WGCNA","Rtsne","feather","shinyjs","reactable")
# notInstalled = setdiff(Rlibs, rownames(installed.packages()))
# if(length(notInstalled)>0)
# 	install.packages(notInstalled)

library(RSQLite,verbose=FALSE)	# for database connection
library(gplots,verbose=FALSE)		# for hierarchical clustering
library(ggplot2,verbose=FALSE)	# graphics
library(e1071,verbose=FALSE) 		# computing kurtosis
library(DT,verbose=FALSE) 		# for renderDataTable
library(plotly,verbose=FALSE) 	# for interactive heatmap
library(reshape2,verbose=FALSE) 	# for melt correlation matrix in heatmap
library(visNetwork) # interative network graphs
# Bioconductor packages
#source("https://bioconductor.org/biocLite.R")
#biocLite(c( "limma", "DESeq2","edgeR","gage", "PGSEA", "fgsea", "ReactomePA", "pathview","PREDA","PREDAsampledata","sfsmisc","lokern","multtest","dplyr"))
# annotation packages needed by pathview; will be installed automatically if runing on Windows
#biocLite( c( "org.Ag.eg.db","org.At.tair.db","org.Bt.eg.db","org.Ce.eg.db",
#"org.Cf.eg.db","org.Dm.eg.db","org.Dr.eg.db","org.EcK12.eg.db","org.EcSakai.eg.db",
#"org.Gg.eg.db","org.Hs.eg.db","org.Hs.ipi.db","org.Mm.eg.db","org.Mmu.eg.db",
#"org.Pf.plasmo.db","org.Pt.eg.db","org.Rn.eg.db","org.Sc.sgd.db","org.Sco.eg.db",
#"org.Ss.eg.db","org.Tgondii.eg.db","org.Xl.eg.db")  )
# auto install 
# biocLibs = c( "limma", "DESeq2","edgeR","gage", "PGSEA", "fgsea", "ReactomePA", "pathview","PREDA","PREDAsampledata","sfsmisc","lokern","multtest","hgu133plus2.db","impute" )
# notInstalled = setdiff(biocLibs, rownames(installed.packages()))
# if(length(notInstalled)>0) { 
# 	source("https://bioconductor.org/biocLite.R")
# 	biocLite(notInstalled)
# }

#library(DESeq2) # count data analysis
#library(edgeR) # count data D.E.

#library(PGSEA) # pathway 

#library(sfsmisc,verbose=FALSE)   #required by PREDA
#library(lokern,verbose=FALSE)	#required by PREDA
#library(multtest,verbose=FALSE)	#required by PREDA

# KEGG and WGCNA generate temporary files. Needs to be deleted regularily. 

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
redudantGeneSetsRatio = 0.9
maxGeneClustering = 12000  # max genes for hierarchical clustering and k-Means clustering. Slow if larger
maxGeneWGCNA = 3000 # max genes for co-expression network
maxFactors =6  # max number of factors in DESeq2 models
set.seed(2) # seed for random number generator
mycolors = sort(rainbow(20))[c(1,20,10,11,2,19,3,12,4,13,5,14,6,15,7,16,8,17,9,18)] # 20 colors for kNN clusters
#Each row of this matrix represents a color scheme;
maxSamples = 100   # DESeq2 gets really slow when more than 50 samples
maxSamplesDefault = 30   # change default from DESeq2 to limma
maxComparisons = 20 # max number of pair wise comparisons in DESeq2
hmcols <- colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF",
"#E0F3F8", "#91BFDB", "#4575B4")))(75)
heatColors = rbind(  greenred(75),     bluered(75),     
                     colorpanel(75,"green", "black","magenta"),
                     colorpanel(75,"blue", "yellow","red"),hmcols )
rownames(heatColors) = c("Green-Black-Red", "Blue-White-Red", "Green-Black-Magenta",
                         "Blue-Yellow-Red", "Blue-white-brown")
colorChoices = setNames(1:dim(heatColors)[1], rownames(heatColors)) # for pull down menu
maxSamplesEDAplot = 100  # max number of samples for EDA plots
STRING_DB_VERSION <- "11.0" # what version of STRINGdb needs to be used 

################################################################
#   Input files
################################################################

# relative path to data files
datapath = "../../data/data104/"   # production server

sqlite  <- dbDriver("SQLite")
convert <- dbConnect( sqlite, paste0(datapath, "convertIDs.db"), flags=SQLITE_RO)  #read only mode
keggSpeciesID = read.csv(paste0(datapath, "data_go/KEGG_Species_ID.csv"))
# List of GMT files in /gmt sub folder
gmtFiles = list.files(path = paste0(datapath,"pathwayDB"), pattern=".*\\.db")
gmtFiles = paste(datapath, "pathwayDB/", gmtFiles,sep="")
geneInfoFiles = list.files(path = paste0(datapath, "geneInfo"), pattern=".*GeneInfo\\.csv")
geneInfoFiles = paste(datapath, "geneInfo/", geneInfoFiles,sep="")
motifFiles = list.files(path = paste0(datapath,"motif"), pattern=".*\\.db")
motifFiles = paste(datapath, "motif/", motifFiles,sep="")
#demoDataFile = paste0(datapath,"data_go/GSE37704_sailfish_genecounts.csv") #"expression1_no_duplicate.csv"
#demoDataFile = paste0(datapath,"data_go/BcellGSE71176_p53.csv") # GSE71176
#demoDataFile2 = paste0(datapath,"data_go/BcellGSE71176_p53_sampleInfo.csv") # sample Info file
demoDataFile = paste0(datapath, "data_go/BcellGSE71176_p53.csv") # GSE71176
demoDataFile2 = paste0(datapath, "data_go/BcellGSE71176_p53_sampleInfo.csv") # sample Info file
quotes <- dbGetQuery(convert, " select * from quotes")
quotes = paste0("\"",quotes$quotes,"\"", " -- ",quotes$author,".       ")

STRING10_species = read.csv(paste0(datapath, "data_go/STRING11_species.csv"))
# File needs to be updated when STRING updates, using the following commands
# library(STRINGdb)
# species = get_STRING_species(version="10", species_name=NULL)
# write.csv(species,"STRING10_species.csv")
# Also this STRINGdb package downloads a lot of file from the website. Needs to clean the temp folder from time to time. 

# prepare species list

# Create a list for Select Input options
orgInfo <- dbGetQuery(convert, paste("select distinct * from orgInfo " ))
orgInfo <- orgInfo[order(orgInfo$name),]
annotatedSpeciesCounts <- sort( table(orgInfo$group) ) # total species, Ensembl, Plants, Metazoa, STRINGv10
speciesChoice <- setNames(as.list( orgInfo$id ), orgInfo$name2 )
# add a defult element to list    # new element name       value

speciesChoice <- append( setNames( "NEW","**NEW SPECIES**"), speciesChoice  )
speciesChoice <- append( setNames( "BestMatch","Best matching species"), speciesChoice  )

# move one element to the 2nd place
move2 <- function(i) c(speciesChoice[1:2],speciesChoice[i],speciesChoice[-c(1,2,i)])
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
################################################################
#   Utility functions
################################################################

# Functions for hierarchical clustering
hclust2 <- function(x, method="average", ...)  # average linkage
  hclust(x, method=method, ...)
hclust.ward.D <- function(x, method="ward.D", ...)  # ward.D linkage
  hclust(x, method=method, ...)
hclust.ward.D2 <- function(x, method="ward.D2", ...)  # ward.D2 linkage
  hclust(x, method=method, ...)
hclust.single <- function(x, method="single", ...)  # single linkage
  hclust(x, method=method, ...)
hclust.mcquitty <- function(x, method="mcquitty", ...)  # mcquitty linkage
  hclust(x, method=method, ...)
hclust.median <- function(x, method="median", ...)  # median linkage
  hclust(x, method=method, ...)
hclust.centroid <- function(x, method="centroid", ...)  # centroid linkage
  hclust(x, method=method, ...)
  
hclustFuns <- list( averge   = hclust2, 
                    complete = hclust, 
                    single   = hclust.single,
					median   = hclust.median, 
                    centroid = hclust.centroid, 
                    mcquitty = hclust.mcquitty)

hclustChoices = setNames(1:length(hclustFuns),names(hclustFuns)) # for pull down menu

dist2 <- function(x, ...)   # distance function = 1-PCC (Pearson's correlation coefficient)
  as.dist(1-cor(t(x), method="pearson"))
  
dist3 <- function(x, ...)   # distance function = 1-abs(PCC) (Pearson's correlation coefficient)
  as.dist(1-abs(cor(t(x), method="pearson")))   
  
# List of distance functions 
distFuns <- list(Correlation=dist2, Euclidean=dist,AbsolutePCC=dist3)
distChoices = setNames(1:length(distFuns),names(distFuns)) # for pull down menu


geneChange <- function(x){
  # Given a set of numbers, find the difference between 2nd largest and 2nd smallest
  # 2,3,5,6,1   --> 5-2 = 3
  n <- length(x)
  if( n<4) { 
    return( max(x)-min(x) )
  } else { 
	return( sort(x)[n-1] - sort(x)[2] )
  }
}

dynamicRange <- function( x ) {
  # Given a set of numbers, find the difference between 2nd largest and 2nd smallest 
  y <- sort(x)
  if(length(x)>=4) {
    k <- 2 
  } else {
    k <- 1
  }
  return( y[length(x)-k+1] - y[k] ) 
}  


 detectGroups <- function (x, sampleInfo = NULL){  # x are col names
# parsing samples by either the name or using a data frame of sample infos. 
# Note that each row of the sampleInfo data frame represents a sample.
# Revised 4-19-2020
  if(is.null(sampleInfo)) {
  #if(1) {
    # Define sample groups based on column names
    # Args:
    #   x are vector of characters, column names in data file
    # Returns: 
    #   a character vector, representing sample groups.
    g <- gsub("[0-9]*$","",x) # Remove all numbers from end
    #g = gsub("_Rep|_rep|_REP","", g)
    g <- gsub("_$", "", g); # remove "_" from end
    g <- gsub("_Rep$", "", g); # remove "_Rep" from end
    g <- gsub("_rep$", "", g); # remove "_rep" from end
    g <- gsub("_REP$", "", g)  # remove "_REP" from end
    return( g ) 
  } else {
    
    # the orders of samples might not be the same.The total number of samples might also differ
    iy = match(x, row.names(sampleInfo))
    sampleInfo2 = sampleInfo[iy, , drop = FALSE]
    
   if(ncol(sampleInfo2) == 1) {  # if there's only one factor
     g = sampleInfo2[, 1] 
     
     } else {   # multiple columns/factors

      #old comments: This does not work as sometimes there is no replicates wt vs mt, 1, 2, 3 paired replicates as factors 
      g = unlist( apply(sampleInfo2, 1, function (y) paste(y, collapse = "_")) )
      names(g) = row.names(sampleInfo2)
      
      if( min( table(g) ) ==  1 ) # no replicates? 
         g = sampleInfo2[, 1]        

     }
   }
   
  return( as.character(g) )
 }

 plotGenes <- function(convertedData, allGeneInfo, readSampleInfo, geneSearch, genePlotBox, useSD, selectOrg){
   # plot the expression of one or more genes in the preprocess tab
   x <- convertedData
   
   Symbols <- rownames(x)
   
   if( selectOrg != "NEW" &&  ncol(allGeneInfo) != 1 ) {
     ix = match( rownames(x), allGeneInfo[,1])
     if( sum( is.na(allGeneInfo$symbol )) != dim(allGeneInfo )[1] ) {  # symbol really exists? 
       Symbols = as.character( allGeneInfo$symbol[ix] )
       Symbols[which( nchar(Symbols) <= 2 ) ] <- rownames(x) [which( nchar(Symbols) <= 2 ) ] 
     }
   }
   x = as.data.frame(x)
   x$Genes = Symbols
   
   # matching from the beginning of symbol
   searchWord = gsub("^ ", "", geneSearch )
   ix = which(regexpr(  paste("^" , toupper(searchWord),sep="")   ,toupper(x$Genes)) > 0)
   if(grepl(" $", searchWord)  )  # if there is space character at the end, do exact match
     ix = match(gsub(" ","", toupper(searchWord)), toupper(x$Genes) )
   
   if(grepl(",|;", searchWord)  ) { # if there is comma or semicolon, split into multiple words
     Words <- unlist( strsplit(searchWord,",|;") ) # split words
     Words <- gsub(" ", "", Words)
     ix = match( toupper(Words), toupper(x$Genes) )
   }
   ix = ix[!is.na(ix)] # remove NAs
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
   g = detectGroups(mdf$samples, readSampleInfo)
   mdf$g = g	
   
   options(dplyr.summarise.inform = FALSE)
   #calculate mean, SD, N, per gene per condition
   summarized <- mdf %>% 
     group_by(g, Genes) %>%  
     summarise(Mean = mean(value), SD = sd(value), N = sum(count))
   colnames(summarized)= c("Samples","Genes","Mean","SD","N")
   summarized$SE = summarized$SD / sqrt(summarized$N)	
   
   if(grepl(",|;", searchWord)  ) { # re-order according to user input, not alphabetically
     levels <- unique(summarized$Genes)
     iy <- match(toupper(Words), toupper(levels) )
     levels <- levels[iy]
     summarized$Genes <- factor( summarized$Genes, levels = levels) 
   }
   
   #http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization
   p2 <- ggplot(summarized, aes(x=Genes, y=Mean,fill=Samples) ) + # data & aesthetic mapping
     geom_bar(stat="identity", position=position_dodge()) + # bars represent average
     geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.2,position=position_dodge(.9)) +
     labs(y="Expression Level") 
   if(useSD == 1) { 
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
   
   if( genePlotBox == 1)  
     return(p1) else 
     return(p2)
   
 }
 
# heatmap with color bar define gene groups
myheatmap2 <- function (x,bar=NULL,n=-1,mycolor=1,clusterNames=NULL,sideColors=NULL ) {
	# number of genes to show
	ngenes = as.character( table(bar))
	if(length(bar) >n && n != -1) {ix = sort( sample(1:length(bar),n) ); bar = bar[ix]; x = x[ix,]  }
	if(! is.null(bar) )
		if(is.null(sideColors) ) 
			sideColors = mycolors

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
			,key=F, labRow = F
			#,RowSideColors = mycolors[bar]
			,margins = c(8, 24)
			,srtCol=45
		) else
		heatmap.2(x,  Rowv =F,Colv=F, dendrogram ="none",
			col=heatColors[as.integer(mycolor),], density.info="none", trace="none", scale="none", keysize=.3
			,key=F, labRow = F
			,RowSideColors = sideColors[bar]
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
		col = sideColors,  # color key
		lty= 1,             # line style
		lwd = 10 )           # line width
		}
}

 #Change comparison names in limma from "KO_ko-WT_ko" to    "KO-WT_for_ko" 
changeNames <- function (comp) {
 	 if( !grepl(".*_.*-.*_.*",comp )) return(comp)
	 levels4 = unlist( strsplit( unlist( strsplit(comp,"-") ), "_") ) #double split!
	 if(length(levels4)!=4) return(comp)
	 ix = which(duplicated(levels4))
	 return( paste0( paste0( levels4[-c(ix, ix-2)] ,collapse ="-"), "_for_",levels4[ix] ) )
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
# updated 10/15; some changes not included in Gavin's new version
convertID <- function (query,selectOrg) {
	querySet <- cleanGeneSet( unlist( strsplit( toupper(query),'\t| |\n|\\,')))
	# querySet is ensgene data for example, ENSG00000198888, ENSG00000198763, ENSG00000198804
    querSetString <- paste0("('", paste(querySet,collapse="', '"),"')")
	# ('ENSG00000198888', 'ENSG00000198763', 'ENSG00000198804')

	if(selectOrg == speciesChoice[[1]]) {# if best match

	  #First send a query to determine the species
	  query_species <- paste0( "select species, idType, COUNT(species) as freq from mapping where id IN ", 
	                      querSetString," GROUP by species,idType")
	  species_ranked <- dbGetQuery(convert, query_species)

	  if( dim(species_ranked)[1] == 0  ) return(NULL)	  	
	  sortedCounts <- species_ranked$freq 
	  names(sortedCounts) <- paste(species_ranked$species, species_ranked$idType)
	  sortedCounts <- sort(sortedCounts, decreasing = TRUE)

		# Try to use Ensembl instead of STRING-db genome annotation
		if(length(sortedCounts) > 1) # if more than 1 species matched
        if( sortedCounts[1] <= sortedCounts[2] *1.1  # if the #1 species and #2 are close
             && as.numeric( gsub(" .*", "", names(sortedCounts[1]))) > sum( annotatedSpeciesCounts[1:3])  # 1:3 are Ensembl species
             && as.numeric( gsub(" .*", "", names(sortedCounts[2]))) < sum( annotatedSpeciesCounts[1:3])    ) {
		  tem <- sortedCounts[2]
		  sortedCounts[2] <- sortedCounts[1]
		  names(sortedCounts)[2] <- names(sortedCounts)[1]
		   sortedCounts[1] <- tem
		  names(sortedCounts)[1] <- names(tem)    
		} 
		recognized =names(sortedCounts[1])
		
		speciesMatched=sortedCounts
		speciesMatched <- as.data.frame( speciesMatched )	
		orgName <- sapply(as.numeric(gsub(" .*","",names(sortedCounts) ) ), findSpeciesByIdName  )
		speciesMatched <- cbind( orgName,  speciesMatched)

		if(length(sortedCounts) == 1) { # if only  one species matched
		   speciesMatched[1,1] <-paste( speciesMatched[1,1], "(",speciesMatched[1,2],")",sep="")
		   speciesMatched <- speciesMatched[, 1, drop = FALSE]
		} else {# if more than one species matched
            speciesMatched <- speciesMatched[!duplicated(speciesMatched[, 1]), ] # same species different mapping (ensembl, arayexpress, hpa)
			speciesMatched[,1] <- as.character(speciesMatched[,1])
			speciesMatched[,1] <- paste( speciesMatched[,1]," (",speciesMatched[,2], ")", sep="") 
			speciesMatched[1,1] <- paste( speciesMatched[1,1],"   ***Used in mapping***  To change, select from above and resubmit query.") 	
			speciesMatched <- as.data.frame(speciesMatched[,1])
		}

	
		querySTMT <- paste0("select distinct id,ens,species,idType from mapping where ",  
		                    " species = '", gsub(" .*","", recognized), "'",
		                    " AND idType = '", gsub(".* ","", recognized ), "'",
		                    " AND id IN ", querSetString)
		
		result <- dbGetQuery(convert, querySTMT)
		if( dim(result)[1] == 0  ) return(NULL)		
		
		
	} else { # if species is selected

	  querySTMT <- paste0( "select distinct id,ens,species,idType from mapping where species = '", selectOrg,
	                      "' AND id IN ", querSetString) 
	  result <- dbGetQuery(convert, querySTMT)

	  if( dim(result)[1] == 0  ) return(NULL)
		result <- result[which(result$species == selectOrg ) ,]
		if( dim(result)[1] == 0  ) return(NULL) #stop("ID not recognized!")
		speciesMatched <- as.data.frame(paste("Using selected species ", findSpeciesByIdName(selectOrg) )  )
	}
	result <- result[which(!duplicated(result[,2]) ),] # remove duplicates in ensembl_gene_id
	result <- result[which(!duplicated(result[,1]) ),] # remove duplicates in user ID
	colnames(speciesMatched) = c("Matched Species (genes)" ) 
	conversionTable <- result[,1:2]; colnames(conversionTable) = c("User_input","ensembl_gene_id")
	conversionTable$Species = sapply(result[,3], findSpeciesByIdName )

	return(list(originalIDs = querySet,
                IDs=unique( result[,2]), 
				species = findSpeciesById(result$species[1]), 
				#idType = findIDtypeById(result$idType[1] ),
				speciesMatched = speciesMatched,
				conversionTable = conversionTable
				) )
}

# finds id index corresponding to entrez gene and KEGG for id conversion ; from ensembl 100 changed from 'entrezgene' to 'entrezgene_id'
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
	{ x = read.csv(as.character(geneInfoFiles[ix]) ); 
      x[,1]= toupper(x[,1]) 
      # if symbol is missing use Ensembl IDs
      x$symbol[ is.na( x$symbol) ] <- x[, 1]
      # if duplicated symbol, paste Ensembl id to the end
      n_occur <- data.frame(table(x$symbol))
      ix_duplicated <- which(n_occur$Freq > 1) # rows with duplicated symbols
      x$symbol[ix_duplicated] <- paste(x$symbol[ix_duplicated], x[ix_duplicated, 1] )

    } else # read in the chosen file 
	{ return(as.data.frame("Multiple geneInfo file found!") )   }
	Set = match(x$ensembl_gene_id, querySet)
	Set[which(is.na(Set))]="Genome"
	Set[which(Set!="Genome")] ="List"
	# x = cbind(x,Set) } # just for debuging
	return( cbind(x,Set) )}
 }

hyperText <- function (textVector, urlVector){
  # for generating pathway lists that can be clicked.
  # Function that takes a vector of strings and a vector of URLs
  # and generate hyper text 
  # add URL to Description 
  # see https://stackoverflow.com/questions/30901027/convert-a-column-of-text-urls-into-active-hyperlinks-in-shiny
  # see https://stackoverflow.com/questions/21909826/r-shiny-open-the-urls-from-rendertable-in-a-new-tab
  if( sum(is.null(urlVector) ) == length(urlVector) )
     return(textVector)
 
  if(length(textVector) != length(urlVector))
    return(textVector)

  #------------------URL correction
  # URL changed from http://amigo.geneontology.org/cgi-bin/amigo/term_details?term=GO:0000077 
  #                  http://amigo.geneontology.org/amigo/term/GO:0000077
  urlVector <- gsub("cgi-bin/amigo/term_details\\?term=", "amigo/term/", urlVector )
  urlVector <- gsub(" ", "", urlVector )


  # first see if URL is contained in memo
  ix <- grepl("http:", urlVector, ignore.case = TRUE)  
  if(sum(ix) > 0) { # at least one has http?
    tem <- paste0("<a href='", 
      urlVector, "' target='_blank'>",
      textVector, 
      "</a>" )
    # only change the ones with URL
    textVector[ix] <- tem[ix]
  }
  return(textVector)
}

removeHypertext <- function( df ){
# Given data frame, remove hypertext in all columns
  if(class(df) == "data.frame") {
    for( i in 1:dim(df)[2])
      if(is.character(df[, i]))
        df[, i] <- gsub(".*'_blank'>|</a>", "", df[, i]) 
    }
  return(df)
}

# Main function. Find a query set of genes enriched with functional category
FindOverlap <- function (converted,gInfo, GO,selectOrg,minFDR, reduced = FALSE, convertedDataBackground = NULL) {
    maxGenesBackground <- 30000
    minGenesBackground <- 2000
	maxTerms =15 # max number of enriched terms
    maxPvalFilter = 0.3
	idNotRecognized = as.data.frame("ID not recognized!")

	if(is.null(converted) ) return(idNotRecognized) # no ID 
	
    querySet <- converted$IDs

  if(!is.null(gInfo) )
     if( class(gInfo) == "data.frame" )
       if(dim(gInfo)[1] > 1) {  # some species does not have geneInfo. STRING
	     # only coding
	     querySet <- intersect( querySet, 
                                gInfo[which( gInfo$gene_biotype == "protein_coding"), 1] )
	}

	if(length(querySet) == 0) return(idNotRecognized )

	ix = grep(converted$species[1,1],gmtFiles)
	totalGenes <- converted$species[1,7]
    errorMessage = as.data.frame("Annotation file cannot be found")
    
	if (length(ix) == 0 ) {return(errorMessage )}
	
	# If selected species is not the default "bestMatch", use that species directly
	if(selectOrg != speciesChoice[[1]]) {  
		ix = grep(findSpeciesById(selectOrg)[1,1], gmtFiles )
		if (length(ix) == 0 ) {return(idNotRecognized )}
		totalGenes <- orgInfo[which(orgInfo$id == as.numeric(selectOrg)),7]
	}

	pathway <- dbConnect(sqlite,gmtFiles[ix],flags=SQLITE_RO)	
		
  if( GO != "All") {
     sqlQuery = paste( " select distinct gene,pathwayID from pathway where category='", GO, "'",
                          " AND gene IN ('", paste(querySet, collapse="', '"),"')" ,sep="")
   } else {
     sqlQuery = paste( " select distinct gene,pathwayID from pathway where gene IN ('", 
                        paste(querySet, collapse="', '"),"')" ,sep="")
   }

	result <- dbGetQuery( pathway, sqlQuery  )

	if( dim(result)[1] ==0) {return(as.data.frame("No matching species or gene ID file!" )) }

	# given a pathway id, it finds the overlapped genes, symbol preferred
	sharedGenesPrefered <- function(pathwayID) {
		tem <- result[which(result[,2]== pathwayID ),1]
		ix = match(tem, converted$conversionTable$ensembl_gene_id) # convert back to original
		tem2 <- unique( converted$conversionTable$User_input[ix] )
        if(!is.null(gInfo) )
          if( class(gInfo) == "data.frame")
            if( dim(gInfo)[1] > 1)
              if(length(unique(gInfo$symbol) )/dim(gInfo)[1] >.7  ) { # if 70% genes has symbol in geneInfo
    	        ix = match(tem, gInfo$ensembl_gene_id);
    	        tem2 <- unique( gInfo$symbol[ix] )      }
	     return( paste( tem2 ,collapse=" ",sep="") )
    }
	
	x0 = table(result$pathwayID)					
	x0 = as.data.frame( x0[which(x0>=Min_overlap)] )# remove low overlaps
    errorMessage = as.data.frame("Too few genes.")
	if(dim(x0)[1] <= 2 ) return(errorMessage) # no data
	colnames(x0)=c("pathwayID","overlap")
	pathwayInfo <- dbGetQuery( pathway, paste( " select distinct id,n,description,memo from pathwayInfo where id IN ('", 
							paste(x0$pathwayID,collapse="', '"),   "') ",sep="") )
	
    # create hypertext
    pathwayInfo$description <- hyperText(pathwayInfo$description, pathwayInfo$memo )
   
    pathwayInfo <- subset( pathwayInfo, select = -c(memo))
	x = merge(x0,pathwayInfo, by.x='pathwayID', by.y='id')
    
    # filtered pathways with enrichment ratio less than one
    x <- x[ which( x$overlap/ length(querySet) / (as.numeric(x$n) / totalGenes ) > 1)  ,]
    x$Pval <- phyper(x$overlap - 1,
				length(querySet),
				totalGenes - length(querySet),   
				as.numeric(x$n), 
				lower.tail=FALSE );
    # further filter by P value; if Pval is big, we assume that using the background genes will not change that.
    #x <- subset(x, Pval < maxPvalFilter)

	  #Background genes----------------------------------------------------

  if(!is.null(convertedDataBackground))
  if( length( intersect( convertedDataBackground$IDs, querySet ) ) > 0.5 * length(querySet) &&  # Species matching? At least half of the query in background
     length(convertedDataBackground$IDs) > minGenesBackground &&  # if too few genes, use all
     length(convertedDataBackground$IDs) < maxGenesBackground + 1) { # if more than 30k genes, ignore background genes.
     
     querySetB <- convertedDataBackground$IDs # all genes in the converted GEnes  
     if(!is.null(gInfo) )
         if(dim(gInfo)[1] > 1) {  # some species does not have geneInfo. STRING
	          # only coding
	          querySetB <- intersect( querySetB, 
                             gInfo[which( gInfo$gene_biotype == "protein_coding"),1] )
	      }  

        if( length( intersect( querySetB, querySet ) ) == 0 )    # if none of the selected genes are in background genes
          return(list( x=as.data.frame("None of the selected genes are in the background genes!" )) )
        
        querySetB <- unique( c( querySetB, querySet ) )  # just to make sure the background set includes the query set
        
        sqlQueryB = paste( " select distinct gene,pathwayID from pathway where gene IN ('", 
                           paste(querySetB, collapse="', '"),"')" ,sep="")    
        # restrict to pathways with genes matching querySet. This improves the speed drastically
        sqlQueryB = paste0(sqlQueryB, " AND pathwayID IN ('", paste(x$pathwayID, collapse="', '"),"')"  )
        
        if( GO != "All") sqlQueryB = paste0(sqlQueryB, " AND category ='",GO,"'")
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
        
      }
  
  # end background genes------------------------------------------------------------
  


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
		if(reduced != FALSE && dim(x)[1] > 5){  # reduced=FALSE no filtering,  reduced = 0.9 filter sets overlap with 90%
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
	pathway <- dbConnect(sqlite,gmtFiles[ix],flags=SQLITE_RO)
	
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
	library(PGSEA,verbose=FALSE)
	Pvalue = 0.01  # cut off to report in PGSEA. Otherwise NA
	#Pval_pathway = 0.2   # cut off for P value of ANOVA test  to writ to file 
	# top = 30   # number of pathways to show
	if(length(gmt) ==0 ) return( list(pg3 = NULL, best = 1 ) )
    # centering by mean
	#pg = myPGSEA (convertedData - rowMeans(convertedData),
	 #            cl=gmt,range=myrange,p.value=TRUE, weighted=FALSE,nPermutation=100)

	#if( class(convertedData) != "data.frame" | class(convertedData) != "matarix") return( list(pg3 = NULL, best = 1 ) )
	#if( dim(convertedData)[2] <2 ) return( list(pg3 = NULL, best = 1 ) )
	
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
	a=sprintf("%-3.2e",pg3[,1])
	rownames(pg3) = paste(a,rownames(pg3),sep=" ")
	pg3 =pg3[,-1]
	
	pg3 <- pg3[order( -apply(pg3,1,sd)    ),] # sort by SD
 
    return( list(pg3 = pg3, best = best ) )
    }
 }


# Differential expression using LIMMA 
DEG.limma <- function (x, maxP_limma=.1, minFC_limma=2, rawCounts,countsDEGMethods,priorCounts, dataFormat, selectedComparisons=NULL, sampleInfo = NULL,modelFactors=NULL, blockFactor = NULL){
	library(limma,verbose=FALSE) # Differential expression
	library(statmod,verbose=FALSE)
	
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
	groups =  detectGroups( groups, sampleInfo) 
	g =  unique(groups)  
	
	# check for replicates, removes samples without replicates
	reps = as.matrix(table(groups)) # number of replicates per biological sample
	if ( sum( reps[,1] >= 2) <2    ) # if less than 2 samples with replicates
	return( list(results= NULL, comparisons = NULL, Exp.type="Failed to parse sample names to define groups. 
		Cannot perform DEGs and pathway analysis. Please double check column names! Use WT_Rep1, WT_Rep2 etc. ", topGenes=NULL)) 
	# remove samples without replicates
	g <- rownames(reps)[which(reps[,1] >1)]
	ix <- which( groups %in% g)  
	groups <- groups[ix]   
	x<- x[,ix]; rawCounts <- rawCounts[,ix] 
	
	if(length(g) ==2 ) {  # just two groups
		g= unique(groups)
		comparisons <-  paste(g[2],"-",g[1],sep="")  # "Mutant-WT"
		
		# no sample file, but user selected comparisons using column names
		if( is.null(modelFactors) & length( selectedComparisons) >0  ) 	
			comparisons <- selectedComparisons

        # set reference level based on the order in which the levels appear
        # the first appearing level is set as reference; otherwise, we get up and down-regulation reversed.
        groups <- factor(groups, levels = g) 

		design <- model.matrix( ~ 0 + groups )
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

		topGenes1 =topTable(fit2, number = 1e12,sort.by="M" )
		if (dim(topGenes1)[1] != 0) {
		topGenes1 = topGenes1[,c('logFC','adj.P.Val')] 
		# topGenes1[,1] <-  -1* topGenes1[,1] # reverse direction
		topGenes[[1]] <- topGenes1 }
		# log fold change is actually substract of means. So if the data is natral log transformed, it shoudl be natral log.
		Exp.type = "2 sample groups."

       # Up and down regulation seems to be reversed? 3-2-2019
       #results <- -1 * results # holds up or down regulated genes, marked as 1 or -1.
       #for(i in 1:length(topGenes))
         #topGenes[[i]][ ,1] <- -1 * topGenes[[i]][ ,1]
	}
	
	if(length(g) > 2 ) { # more than two sample groups
	    # set reference level based on the order in which the levels appear
        # the first appearing level is set as reference; otherwise, we get up and down-regulation reversed.
        groups <- factor(groups, levels = g)

		design <- model.matrix(~ 0 + groups)
		#colnames(design) <- gsub(".*)","",colnames(design))  
		colnames(design) <- gsub("^groups","",colnames(design))

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
					contrast.names = c(contrast.names, paste("I:",  genotypes[kp], "-", genotypes[kk],"(",gsub("-",".vs.",treatments[kp]),")",sep="" ) )
				}   
			}
			colnames(contrast1)=contrast.names
			comparisons = contrast.names

		} # if interaction terms	
		

		# if sample information is uploaded and user selected factors and comparisons
		if( !is.null(modelFactors) & length( selectedComparisons) >0  ) {
			Exp.type = paste("Model: ~", paste(modelFactors,collapse=" + ") )
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

            groups <- factor(groups, levels = g)

		    design <- model.matrix(~ 0 + groups)
		    #colnames(design) <- gsub(".*)","",colnames(design))  
		    colnames(design) <- gsub("^groups","",colnames(design))
            	
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
							
							contrast.names = c(contrast.names, paste0("I:",  colnames(contrast1)[kp], ".vs.", colnames(contrast1)[kk] ) )
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

			if(length(blockFactor) >= 1 ) { # if a factor is selected as block
				if(length(blockFactor) >= 1 ) 
					blockFactor = blockFactor[1] # if multiple use the first one

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

		#colnames(results) = unlist(sapply( colnames(results), changeNames  ) )
		#comparisons3 <-  unlist(sapply( comparisons, changeNames  ) ) 

		# extract fold change for each comparison
		# there is issues with direction of foldchange. Sometimes opposite
		top <- function (comp) {
			tem <- topTable(fit2, number = 1e12,coef=comp,sort.by="M" ) 
			if(dim(tem)[1] == 0) { return (1) 
			} else	{ 			
				# compute fold change for the first gene (ranked by absolute value)
				tem2 = as.numeric( x[ which(rownames(x)== rownames(tem)[1]) , ] )
				names(tem2) = colnames(x) 
					
				return( tem[,c(1,5)]) 
			}  
													
		}  # no significant gene returns 1, otherwise a data frame
		
		
		topGenes <- lapply(comparisons, top)
		topGenes <- setNames(topGenes, comparisons )

		ix <- which( unlist( lapply(topGenes, class) ) == "numeric")
		if( length(ix)>0) topGenes <- topGenes[ - ix ]
		# if (length(topGenes) == 0) topGenes = NULL;
	}
	
	#cat("\n", names(topGenes) )
	#cat("\n", colnames(results))
	#cat("\n", comparisons3)
	#comparisons <- comparisons3
	# it does not make any sense! comparisons can not be changed!
	#comparisons =comparisons3	
	#cat("\n", comparisons)	
    #write.csv(results,"results.csv"); write.csv(topGenes,"topGenes.csv")	
    
	return( list(results= results, comparisons=comparisons, Exp.type=Exp.type, topGenes=topGenes)) 
}

# Differential expression using DESeq2
DEG.DESeq2 <- function (  rawCounts,maxP_limma=.05, minFC_limma=2, selectedComparisons=NULL, sampleInfo = NULL,modelFactors=NULL, blockFactor = NULL, referenceLevels=NULL){
	library(DESeq2,verbose=FALSE) # count data analysis
    #library("BiocParallel")
	groups = as.character ( detectGroups( colnames( rawCounts ), sampleInfo) )
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
	
	# if too many samples 
	if(ncol(rawCounts)  > maxSamples) { 
		#comparisons <- comparisons[1:maxComparisons]
		return( list(results= NULL, comparisons = NULL, 
			Exp.type = paste(Exp.type," Too many samples for DESeq2. Please choose limma-voom or limma-trend." ),
			topGenes=NULL))
		}	
		
	# all pair-wise comparisons

	comparisons = ""
	for( i in 1:(length(g)-1) )
		for (j in (i+1):length(g)) 
		comparisons = c(comparisons,paste(g[j],"-",g[i],sep="" ) )
	comparisons <- comparisons[-1]

   # if too many comparisons 
	if(length(comparisons)  > maxComparisons) { 
		Exp.type = paste(Exp.type," Too many comparisons. Only first",maxComparisons, "of the ", length(comparisons), 
			"comparisons calculated. Please choose comparisons. " )
		comparisons <- comparisons[1:maxComparisons]
		}	
	
	
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

		factors = modelFactors   # selected factors and interactions: c( "strain", "treatment",  "strain:treatment")
		factors = factors[ !grepl(":",factors )]   # non-interaction terms
		# interaction terms like strain:treatment
		Interactions = modelFactors[ grepl(":",modelFactors )]
		
		colData = sampleInfo  
		factorsCoded = toupper(letters )[1: dim(colData)[2] ]   # Factors are encoded as "A", "B", "C"; this avoid illigal letters
		names(factorsCoded) =  colnames(colData)  # this is for look up; each column of sampleInfo  
		colnames(colData) = factorsCoded # all columns named A B C D 

		colData = as.data.frame(colData)
		
		# set reference levels for factors
		if(! is.null( referenceLevels) ) {   # c("genotype:wt", "treatment:control" )
			# first factor
			for ( refs in referenceLevels)
				if(! is.null( refs) ) {
					ix = match(gsub(":.*","",refs), colnames(sampleInfo) ) # corresponding column id for factor
					colData[,ix] = as.factor( colData[,ix] )
					colData[,ix] = relevel(colData[,ix],gsub(".*:","",refs)  )
				}
		}
				
		# base model
        DESeq2.Object= paste("dds = DESeqDataSetFromMatrix(countData=rawCounts, colData=colData, design=~ ", 
							 paste( factorsCoded[factors],collapse="+")) # only use selected factors		
		Exp.type = paste("Model: ~", paste(modelFactors,collapse=" + ") )

		
		# create model
		if( length(Interactions)>0 ) { # if there is interaction
			for( interactionTerms in Interactions) {
				interactingFactors = unlist(strsplit(interactionTerms,":" ) )  # split strain:treatment as "strain" and "mutant"
				tem = paste(factorsCoded [ interactingFactors ],collapse=":")   # convert "strain:mutant" to "A:B"
				DESeq2.Object = paste(DESeq2.Object, " + ",tem)
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
		
		# comparison2 holds names for display with real factor names
		# comparison  is used in calculation it is A, B, C for factors
		comparisons2 = comparisons    
		#comparisons2 = gsub(" vs\\. ","-",selectedComparisons) 
		#comparisons2 = gsub(":","_",comparisons2)		
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
				
				tem_factors[1] = names(factorsCoded)[factorsCoded == tem_factors[1]]  # get the first letter and translate into real factor names
				tem_factors[2] = names(factorsCoded)[factorsCoded == tem_factors[2]]  # get the 2nd letters and translate into real factor names

				interactionComparisons2[i] <- paste0( "I:",tem_factors[1], "_",substr(tem[1],2,nchar(tem[1]) ),".",
													          tem_factors[2], "_",substr(tem[2],2,nchar(tem[2]) ) 
													)				
			}
			comparisons2 = c(comparisons2,interactionComparisons2 )
		}
	} # if selected factors	
	
	# extract contrasts according to comparisons defined above
	result1 = NULL; allCalls = NULL;
	topGenes = list(); pk = 1 # counter
	pp=0 # first results?

	for( kk in 1:length(comparisons) ) {
		tem = unlist( strsplit(comparisons[kk],"-") )
		
		if(is.null(modelFactors)) # if just group comparison using sample names
			selected = results(dds, contrast=c("groups", tem[1], tem[2]) )   else {
			if(!grepl("\\.", comparisons[kk] ) )    # if not interaction term: they contain .  interaction term
				selected = results(dds, contrast=c( factorsCoded[ factorsVector[kk] ],tem[1], tem[2]) ) else # either A, B, C ...
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
				names(topGenes)[pk] = comparisons2[kk];  # assign name to comprison
			}
	}

	Interactions = c()
	if( !is.null(modelFactors) )
		Interactions = modelFactors[ grepl(":",modelFactors )]
		
#---  add comprisons for non-reference levels. It adds to the results1 object.	
	if( length(Interactions)>0 ) { # if there is interaction
		factorLookup=c() # a factor whose values are factors and names are factor and level combination conditionTreated, genotypeWT
		levelLookup = c()
		
		for( i in 1:dim(sampleInfo)[2]) {
			sampleInfo2 = unique(sampleInfo)
			tem = rep(toupper(letters)[i],dim(sampleInfo2)[1]  )
			names(tem) = paste0(toupper(letters)[i],sampleInfo2[,i])
			factorLookup = c(factorLookup,tem)  
			
			tem = as.character( sampleInfo2[,i] )
			names(tem) = paste0(toupper(letters)[i],sampleInfo2[,i])
			levelLookup = c(levelLookup, tem)
		}
		
		# split  genotypeI.conditionTrt --> c("genotype","I","conditoin","Trt")
		splitInteractionTerms <- function (term) {
			if(!grepl("\\.",term) ) return(NULL)
			terms2 = unlist(strsplit(term,"\\.") )
					 # factor1, level1, factor2, level2
			return(c(factorLookup[terms2[1]], levelLookup[terms2[1]],factorLookup[terms2[2]], levelLookup[terms2[2]]   ) )
		}
		# none interaction terms 
		NoneInterTerms = resultsNames(dds)[ !grepl( "\\.", resultsNames(dds)) ]
		NoneInterTerms=NoneInterTerms[-1]
		allInteractionTerms = resultsNames(dds)[ grepl( "\\.", resultsNames(dds)) ]


		for( kk in 1:length(NoneInterTerms) ) { # for each none interaction term
			if(!is.null(modelFactors) ) {# if not just group comparison using sample names
					#current factor
					cFactor = gsub("_.*","",NoneInterTerms[kk] )
				
					for(interactionTerm in allInteractionTerms ) {
					
						splited = splitInteractionTerms (interactionTerm)  # 4 components
						if (cFactor != splited[1] & cFactor != splited[3]  ) 
							next;						
						
						selected = results(dds, list(c( NoneInterTerms[kk],interactionTerm ) ) ) 
						comparisonName = paste0( NoneInterTerms[kk],"__", gsub("\\.","",interactionTerm) )
						
						if( cFactor == splited[1] )
							otherLevel = splited[4] else otherLevel = splited[2]
							
						comparisonName = paste0(#names(factorsCoded)[which(factorsCoded==cFactor)], # real factor name
												gsub("_vs_","-", substr(NoneInterTerms[kk], 3, nchar(NoneInterTerms[kk]  )  )), # the comparison
												"_for_",otherLevel)
						comparisons2 = c(comparisons2, comparisonName)
						selected$calls =0   
						selected$calls [which( selected$log2FoldChange > log2(minFC_limma) & selected$padj < maxP_limma ) ]  <-  1
						selected$calls [ which( selected$log2FoldChange <  -log2(minFC_limma) & selected$padj < maxP_limma ) ] <-  -1
						colnames(selected)= paste( comparisonName, "___",colnames(selected),sep="" )
						selected = as.data.frame(selected)
						if (pp==0){  # if first one with significant genes, collect gene list and Pval+ fold
							result1 = selected; pp = 1; 
							# selected[,2] <- -1 * selected[,2] # reverse fold change direction
							topGenes[[1]] = selected[,c(2,6)]; 
							names(topGenes)[1] = comparisonName; } else 
							{ result1 = merge(result1,selected,by="row.names"); 
								rownames(result1) = result1[,1]; 
								result1 <- result1[,-1]
								pk= pk+1; 
								# selected[,2] <- -1 * selected[,2] # reverse fold change direction
								topGenes[[pk]] = selected[,c(2,6)]; 
								names(topGenes)[pk] = comparisonName;  # assign name to comprison
							}
					} #for	
						
			} #if
		} #for
	
	
	} #if



#---
	#if( length(comparisons) == 1) topGenes <- topGenes[[1]] # if only one comparison, topGenes is not a list, just a data frame itself.
	if(! is.null(result1)) { 
		# note that when you only select 1 column from a data frame it automatically converts to a vector. drop =FALSE prevents that.
		allCalls = as.matrix( result1[,grep("calls",colnames(result1)), drop = FALSE  ] )
		colnames(allCalls)= gsub("___.*","", colnames(allCalls))
		colnames(allCalls)= gsub("\\.","-", colnames(allCalls)) # note that samples names should have no "."
		colnames(allCalls)= gsub("^I-","I:", colnames(allCalls))
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
	
	motifs <- dbConnect(sqlite,motifFiles[ix],flags=SQLITE_RO) # makes a new file

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

# find sample index for selected comparisons
findContrastSamples <- function(selectContrast, allSampleNames,sampleInfo=NULL, selectFactorsModel=NULL,selectModelComprions =NULL , referenceLevels=NULL, countsDEGMethod=NULL, dataFileFormat=NULL ){
	iz= match( detectGroups(allSampleNames), unlist(strsplit( selectContrast, "-"))	  )
	iz = which(!is.na(iz))
	
	# has design file, but didn't select factors
	if ( !is.null(sampleInfo) & is.null(selectFactorsModel) & length(selectModelComprions)==0 ) {
	  
	  findSamples <- function( factorLevel ) { 
	    # given a factor level such as "wt", return a vector indicating the samples with TRUE FALST
	    #  p53_mock_1  p53_mock_2  p53_mock_3  p53_mock_4    p53_IR_1    p53_IR_2    p53_IR_3    p53_IR_4 null_mock_1 null_mock_2 
	    #  TRUE        TRUE        TRUE        TRUE        TRUE        TRUE        TRUE        TRUE       FALSE       FALSE 
  	  tem = apply(sampleInfo, 2, function(y) y == factorLevel )
  	  colSums(tem) > 0
  	  tem <- tem[, colSums(tem) > 0]
  	  return(tem)
	  }
	  
	  
	  sample1 <- gsub("-.*","",selectContrast)
	  level1 <- gsub("_.*","",sample1)
	  level2 <- gsub(".*_","",sample1)
	  iz <- which( findSamples( level1 ) & findSamples( level2 ) )
	  
	  sample2 <- gsub(".*-","",selectContrast)
	  level1 <- gsub("_.*","",sample2)
	  level2 <- gsub(".*_","",sample2)
	  iz <- c(iz, which( findSamples( level1 ) & findSamples( level2 ) ))	
	  
	  
	}
	
	#Has design file and chose factors
	if ( !is.null(sampleInfo) & !is.null(selectFactorsModel) & length(selectModelComprions)>0 ) {

		comparisons = gsub(".*: ","",selectModelComprions)   # strings like: "groups: mutant vs. control"
		comparisons = gsub(" vs\\. ","-",comparisons)		
		factorsVector= gsub(":.*","",selectModelComprions) # corresponding factors

		  # if read counts data and DESeq2
		if(dataFileFormat==1 & countsDEGMethod == 3) { # if DESeq2
			contrast = gsub("_for_.*","",selectContrast) # could be "wt-mu"   or "wt-mu_for_conditionB"
			ik = match( contrast, comparisons )   # selected contrast lookes like: "mutant-control"

			otherFactorLevel = gsub(".*_for_","",selectContrast)
			# find the corresponding factor for the other factor
			otherFactor=" "
			if(nchar( otherFactorLevel ) >0){
				for( eachFactor in colnames(sampleInfo) )
					if ( otherFactorLevel %in%  sampleInfo[,eachFactor ] )
						otherFactor = eachFactor		
			}
			
			if (is.na(ik)) iz=1:(length(allSampleNames))  else {  # interaction term, use all samples		
				selectedfactor= factorsVector[ ik ] # corresponding factors

				iz = which(sampleInfo[,selectedfactor] %in%  unlist(strsplit( contrast, "-")) )

				#filter by other factors: reference level
				if(! is.null( referenceLevels) ) {   # c("genotype:wt", "treatment:control" )
					for ( refs in referenceLevels)
						if(! is.null( refs) & gsub(":.*","",refs) != selectedfactor ) {
							currentFactor = gsub(":.*","",refs)
							if(nchar( otherFactorLevel ) >0 & currentFactor == otherFactor ) { # if not reference level
								iz = intersect( iz, which(sampleInfo[,currentFactor] == otherFactorLevel  ) )
							} else
								iz = intersect( iz, which(sampleInfo[,currentFactor] == gsub(".*:","",refs)  ) )
							
						}
				}			
				iz = iz[which(!is.na(iz))]
				
			# switching from limma to DESeq2 causes problem, as reference level is not defined.

			} 
			
		} else {  # not DESeq2
			
					# given level find corresponding sample ids
				    findIDsFromLevel <- function (aLevel){
						# find factor
						currentFactor=""
						for( eachFactor in colnames(sampleInfo) )
							if ( aLevel %in%  sampleInfo[,eachFactor ] )
								currentFactor = eachFactor			
						if(nchar(currentFactor) >0 ) 
							return( which(sampleInfo[,currentFactor ] %in% aLevel   )  ) else return(NULL)
					}		
					
					if( !grepl(".*_.*-.*_.*",selectContrast )) iz = c()
					levels4 = unlist( strsplit( unlist( strsplit(selectContrast,"-") ), "_") ) #double split!
					if(length(levels4)!=4) { 
						iz = c() 
					} else {
						iz = intersect( findIDsFromLevel(levels4[1]),  findIDsFromLevel(levels4[2])  ) # first sample
						iz = c(iz, intersect( findIDsFromLevel(levels4[3]),  findIDsFromLevel(levels4[4])  )   ) # 2nd sample
						}
				} #else
			
			
			
		}	
	 
	 if (grepl("I:",selectContrast)) iz=1:length(allSampleNames) # if it is factor design use all samples
	 if( is.na(iz)[1] | length(iz)<=1 )    iz=1:length(allSampleNames)

	return(iz)
}

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

  Terms = paste( sprintf("%-1.0e",as.numeric(enrichedTerms$adj.Pval)), 
				names(geneLists))
  rownames(w) = Terms
  colnames(w) = Terms
  par(mar=c(0,0,1,rightMargin)) # a large margin for showing 

  dend <- as.dist(1-w) %>%
	hclust (method="average") 
  ix = dend$order # permutated order of leaves

  leafType= as.factor( gsub(" .*","", enrichedTerms$Direction[ix] ) )
  if(length(unique(enrichedTerms$Direction)  ) ==2 )
	leafColors = c("green","red")  else  # mycolors
	leafColors = mycolors
	
  #leafSize = unlist( lapply(geneLists,length) ) # leaf size represent number of genes
  #leafSize = sqrt( leafSize[ix] )  
  leafSize = -log10(as.numeric( enrichedTerms$adj.Pval[ix] ) ) # leaf size represent P values
  leafSize = 1.5*leafSize/max( leafSize ) + .2
  
	dend %>% 
	as.dendrogram(hang=-1) %>%
	set("leaves_pch", 19) %>%   # type of marker
	set("leaves_cex", leafSize) %>% #Size
	set("leaves_col", leafColors[leafType]) %>% # up or down genes
	plot(horiz=TRUE)
	
  #legend("top",pch=19, col=leafColors[1:2],legend=levels(leafType),bty = "n",horiz =T  )
  # add legend using a second layer
  	par(lend = 1)           # square line ends for the color legend
	add_legend("top",pch=19, col=leafColors,legend=levels(leafType),bty = "n",horiz =T 

	)
  
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

showGeneIDs <- function(species, nGenes = 10){
# Given a species ID, this function returns 10 gene ids for each idType
    if(species == "BestMatch")
      return(as.data.frame("Select a species above.") )

	idTypes <- dbGetQuery( convert,
						paste0( " select DISTINCT idType from mapping where species = '", species,"'") )	# slow
    idTypes <- idTypes[,1, drop = TRUE]
    
    if(nGenes > 100) nGenes <- 100; # upper limit
    
    # for each id Type
    for(k in 1:length(idTypes)){
        # retrieve 500 gene ids and then random choose 10
		result <- dbGetQuery( convert,
                       paste0( " select  id,idType from mapping where species = '", species,"' 
                                 AND idType ='", idTypes[k], "' 
                                 LIMIT ", 50 * nGenes) )
       result <- result[sample(1:(50 * nGenes), nGenes), ]
       if(k == 1) { 
          resultAll <- result 
       } else { 
         resultAll <- rbind(resultAll, result)
       }
     }

     # Names of idTypes
     idNames <- dbGetQuery( convert,
                            paste0( " SELECT id,idType from idIndex where id IN ('",
                                    paste(idTypes,collapse="', '"),  "') "))
     
     resultAll <- merge(resultAll, idNames, by.x = "idType", by.y = "id")
     
     

     #library(dplyr)
     resultAll <- resultAll %>% 
       select(id, idType.y) %>%
       group_by(idType.y) %>%
       summarise(Examples = paste0(id, collapse = "; "))

       colnames(resultAll)[1] <- "ID Type"
        # put symbols first, refseq next, followed by ensembls. Descriptions (long gnee names) last
        resultAll <- resultAll[ order( grepl("ensembl", resultAll$'ID Type'), decreasing = TRUE), ]    
        resultAll <- resultAll[ order( grepl("refseq", resultAll$'ID Type'), decreasing = TRUE), ]      
        resultAll <- resultAll[ order( grepl("symbol", resultAll$'ID Type'), decreasing = TRUE), ]
        resultAll <- resultAll[ order( grepl("description", resultAll$'ID Type'), decreasing = FALSE), ]
    
    return(resultAll)

}


# Wrapping long text by adding \n 
#  "Mitotic DNA damage checkpoint"  --> "Mitotic DNA damage\ncheckpoint"
# https://stackoverflow.com/questions/7367138/text-wrap-for-plot-titles
wrap_strings <- function( vector_of_strings, width = 30 ) { 
  as.character( sapply( vector_of_strings, FUN=function(x) 
  { paste(strwrap(x, width = width), collapse = "\n")}) )
}

	# output Parameter: This function outputs values in the input$variable to input_variable 
	op <- function (i,para, annot=NULL, textInfo = FALSE) {
		if(is.null(para)  )
			return( paste0(i,"\n ", gsub("\\$","_",
					deparse(substitute(para)))," <- NULL \t#",annot)  )	else
		if(length(para) == 1  ) { # only 1 elements
			if(! textInfo )
				return( paste0(i,"\n ", gsub("\\$","_",
						deparse(substitute(para)))," <- ", para, "\t#",annot)  ) else		
				return( paste0(i,"\n ", gsub("\\$","_",
					deparse(substitute(para)))," <- \'", para, "'\t#",annot) )
		} else # a vector, note must be string vector   c("p53","treatment")
				return( paste0(i,"\n ", gsub("\\$","_",
					deparse(substitute(para)))," <- c('", paste0(para, collapse="','" ), "')\t#",annot) )

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
options(shiny.maxRequestSize = 200*1024^2) # 200MB file max for upload
observe({  updateSelectizeInput(session, "selectOrg", choices = speciesChoice, selected = speciesChoice[1], server = TRUE )      })

  # for gene ID example
  observe({  updateSelectizeInput(session, "userSpecieIDexample", choices = speciesChoice, selected = speciesChoice[1] )      })  

observe({  updateSelectInput(session, "heatColors1", choices = colorChoices )      })
observe({  updateSelectInput(session, "distFunctions", choices = distChoices )      })
observe({  updateSelectInput(session, "hclustFunctions", choices = hclustChoices )      })
# update species for STRING-db related API access
observe({  	updateSelectizeInput(session, "speciesName", choices = sort(STRING10_species$official_name) ) 	})

	################################################################
	#   Read data
	################################################################
 
	# read data file and do filtering and transforming
readData <- reactive ({
		inFile <- input$file1
		inFile <- inFile$datapath

		if(is.null(input$file1) && input$goButton == 0)   return(NULL)
		if(is.null(input$file1) && input$goButton > 0 )   inFile = demoDataFile
		tem = input$dataFileFormat; tem=input$missingValue
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
				tem = input$NminSamples2
		}

		isolate({
			withProgress(message="Reading and pre-processing ", {
				# these packages moved here to reduce loading time
				library(edgeR,verbose=FALSE) # count data D.E.
				library(DESeq2,verbose=FALSE) # count data analysis

				if (is.null( input$dataFileFormat )) return(NULL)
				dataTypeWarning =0
				dataType =c(TRUE)

				#---------------Read file
				x <- read.csv(inFile,quote = "",comment.char="")	# try CSV
				if(dim(x)[2] <= 2 )   # if less than 3 columns, try tab-deliminated
					x <- read.table(inFile, sep="\t",header=TRUE,quote = "",comment.char="")	
				#-------Remove non-numeric columns, except the first column
				
				for(i in 2:dim(x)[2])
					dataType = c( dataType, is.numeric(x[,i]) )
				if(sum(dataType) <=2) return (NULL)  # only less than 2 columns are numbers
				x <- x[,dataType]  # only keep numeric columns
				
				# rows with all missing values
				ix = which( apply(x[,-1],1, function(y) sum( is.na(y) ) ) != dim(x)[2]-1 )
				x <- x[ix,]
				
				dataSizeOriginal = dim(x); dataSizeOriginal[2] = dataSizeOriginal[2] -1
				
				x[,1] <- toupper(x[,1])
				#x[,1] <- gsub(" |\"|\'|\\.[0-9]{1,2}$", "", x[ , 1]) 
                 x[,1] <- gsub(" |\"|\'", "", x[ , 1]) 
				             # remove spaces in gene ids
				                 # remove " in gene ids, mess up SQL query				
				                      # remove ' in gene ids		
				                        # remove one or two digits after "." at the end.
                                          #A35244.1 -> A35244  or A35244.23 -> A35244, but not more than two.  GLYMA.18G52160 stays the same.	
				
				x = x[order(- apply(x[,2:dim(x)[2]],1,sd) ),]  # sort by SD
				x <- x[!duplicated(x[,1]) ,]  # remove duplicated genes
				rownames(x) <- x[,1]
				x <- as.matrix(x[,c(-1)])
				
				# remove "-" or "." from sample names
				colnames(x) = gsub("-","",colnames(x))
				colnames(x) = gsub("\\.","",colnames(x))

				
				#cat("\nhere",dim(x))
				# missng value for median value
				if(sum(is.na(x))>0) {# if there is missing values
					if(input$missingValue =="geneMedian") { 
						rowMedians <- apply(x,1, function (y)  median(y,na.rm=T))
						for( i in 1:dim(x)[2] ) {
							ix = which(is.na(x[,i]) )
							x[ix,i] <- rowMedians[ix]						
						}
							
					} else if(input$missingValue =="treatAsZero") {
						x[is.na(x) ] <- 0					
					} else if (input$missingValue =="geneMedianInGroup") {
						sampleGroups = detectGroups( colnames(x))
						for (group in unique( sampleGroups) ){		
							samples = which( sampleGroups == group )
							rowMedians <- apply(x[,samples, drop=F],1, function (y)  median(y,na.rm=T))
							for( i in  samples ) { 
								ix = which(is.na(x[ ,i] ) )	
								if(length(ix) >0 )
									x[ix, i  ]  <- rowMedians[ix]
							}										
						}
						
						# missing for entire sample group, use median for all samples
						if(sum(is.na(x) )>0 ) { 
							rowMedians <- apply(x,1, function (y)  median(y,na.rm=T))
							for( i in 1:dim(x)[2] ) {
								ix = which(is.na(x[,i]) )
								x[ix,i] <- rowMedians[ix]						
							}						
						}
					}
				}

				# Compute kurtosis
				mean.kurtosis = mean(apply(x,2, kurtosis),na.rm=T)
				rawCounts = NULL
				pvals= NULL
				if (input$dataFileFormat == 2 ) {  # if FPKM, microarray
					incProgress(1/3,"Pre-processing data")

					if ( is.integer(x) ) dataTypeWarning = 1;  # Data appears to be read counts

					#-------------filtering
					#tem <- apply(x,1,max)
					#x <- x[which(tem > input$lowFilter),]  # max by row is at least 
					x <- x[ which( apply( x, 1,  function(y) sum(y >= input$lowFilter)) >= input$NminSamples2 ) , ] 

					
					x <- x[which(apply(x,1, function(y) max(y)- min(y) ) > 0  ),]  # remove rows with all the same levels

					#--------------Log transform
					# Takes log if log is selected OR kurtosis is big than 100
					if ( (input$transform == TRUE) | (mean.kurtosis > kurtosis.log ) ) 
						x = log(x+abs( input$logStart),2)

					tem <- apply(x,1,sd) 
					x <- x[order(-tem),]  # sort by SD

				} else 
					if( input$dataFileFormat == 1) {  # counts data
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
					x <- x[ which( apply( cpm(DGEList(counts = x)), 1,  
							function(y) sum(y>=input$minCounts)) >= input$NminSamples ) , ] 
					

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
					if( input$CountsTransform == 3 ) { # rlog is slow
						#if(dim(x)[2]<=50 )   # disable the menu when large sample size	?				 
						{ x <- rlog(dds, blind=TRUE); x <- assay(x)  } # else 
						# x <- log2( counts(dds, normalized=TRUE) + input$countsLogStart ) 
						 }  

						else {
						if ( input$CountsTransform == 2 ) {    # vst is fast but aggressive
							x <- vst(dds, blind=TRUE)
							x <- assay(x)  
						} else{  # normalized by library sizes and add a constant.
							x <- log2( counts(dds, normalized=TRUE) + input$countsLogStart )   # log(x+c) 
							# This is equivalent to below. But the prior counts is more important
							#x <- cpm(DGEList(counts = x),log=TRUE, prior.count=input$countsLogStart )  #log CPM from edgeR
							#x <- x-min(x)  # shift values to avoid negative numbers
						}
					}
				} else 
					if( input$dataFileFormat == 3)	{  # other data type

						#neg_lfc neg_fdr pos_lfc pos_fdr 
						#11       1      11       1 
						
						n2 = ( dim(x)[2] %/% 2) # 5 --> 2

						if(!input$noFDR) { 	 # if we have FDR
							pvals = x [,2*(1:n2 ),drop=FALSE ]  # 2, 4, 6
							x = x[, 2*(1:n2 )-1,drop=FALSE]   # 1, 3, 5		
                            if(dim(x)[2] == 1) { # if only 1 column of FDR and P value, add another column
                              placeholder <- rep(1, dim(x)[1])
                              pvals <- cbind(pvals, placeholder)
                              zero_placeholder <- rep(0, dim(x)[1])
                              x <- cbind(x, zero_placeholder)
                            }

		
						}	
						
						if(0) {  
						ix =  which(apply(x,1, function(y) max(y)- min(y) ) > 0  )
						x <- x[ix,]  # remove rows with all the same levels
						if(!is.null(pvals) )
							pvals = pvals[ix,]
						 }
					}
					
					
				dataSize = dim(x);

                # this updates sample sizes for scatter plot in the pre-process tab
                sampleChoice <- setNames(as.list( 1:( dim(x)[2] ) ), colnames( x ) )
                observe({  updateSelectInput(session, "scatterX", choices = sampleChoice, selected = sampleChoice[1]) })
                observe({  updateSelectInput(session, "scatterY", choices = sampleChoice, selected = sampleChoice[2]) })

				validate( need(dim(x)[1]>5 & dim(x)[2]>=1 , 
					"Data file not recognized. Please double check."))

				incProgress(1, "Done.")
				
				sampleInfoDemo=NULL
				if( input$goButton >0)
					sampleInfoDemo <- t( read.csv(demoDataFile2,row.names=1,header=T,colClasses="character") )

					finalResult <- list(data = as.matrix(x), mean.kurtosis = mean.kurtosis, rawCounts = rawCounts, dataTypeWarning=dataTypeWarning, dataSize=c(dataSizeOriginal,dataSize),sampleInfoDemo=sampleInfoDemo, pvals =pvals )
				return(finalResult)
			})
		})
	})

readSampleInfo <- reactive ({
		if( is.null(input$file2) && !is.null( readData()$sampleInfoDemo ) ) return( readData()$sampleInfoDemo   )
		inFile <- input$file2
		inFile <- inFile$datapath

		if(is.null(input$file2) && input$goButton == 0)   return(NULL)
		if(is.null(readData() ) ) return(NULL)
		#if(input$goButton2 > 0 )   inFile = demoDataFile2
		
		isolate({
				if (is.null( input$dataFileFormat )) return(NULL)
				dataTypeWarning =0
				dataType =c(TRUE)

				#---------------Read file
				x <- read.csv(inFile,row.names=1,header=T,colClasses="character")	# try CSV
				if(dim(x)[2] <= 2 )   # if less than 3 columns, try tab-deliminated
					x <- read.table(inFile, row.names=1,sep="\t",header=TRUE,colClasses="character")
				# remove "-" or "." from sample names
				colnames(x) = gsub("-","",colnames(x))
				colnames(x) = gsub("\\.","",colnames(x))	

				#----------------Matching with column names of expression file
				ix = match(toupper(colnames(readData()$data)), toupper(colnames(x)) ) 
				ix = ix[which(!is.na(ix))] # remove NA

				validate(
				  need(length(unique(ix) ) == dim(readData()$data)[2] 
				       & dim(x)[1]>=1 & dim(x)[1] <500 # at least one row, it can not be more than 500 rows
					 ,"Error!!! Sample information file not recognized. Sample names must be exactly the same. Each row is a factor. Each column represent a sample.  Please see documentation on format.")
				)
				
				#-----------Double check factor levels, change if needed
				# remove "-" or "." from factor levels
				for( i in 1:dim(x)[1]) {
				   x[i,] = gsub("-","",x[i,])
				   x[i,] = gsub("\\.","",x[i,])
				   x[i,] = gsub(" ","",x[i,])				   
				}
				# if levels from different factors match
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

		if (is.null(readSampleInfo() ) )   return(NULL)
		isolate({
			tem = t(readSampleInfo() )
			tem = cbind(rownames(tem),tem)
			colnames(tem)[1] <- "Study_design"
			return(tem)
		})
	},include.rownames=FALSE,striped=TRUE,bordered = TRUE, width = "auto",hover=T) 	
	
output$textTransform <- renderText({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		inFile <- input$file1
		k.value =  readData()$mean.kurtosis	  
		tem = paste( "Mean Kurtosis =  ", round(k.value,2), ".\n",sep = "")
		if( k.value > kurtosis.log) tem = paste(tem, " Detected extremely large values.  When kurtosis >", kurtosis.log,
		", log transformation is automatically applied.", sep="") else if (k.value>kurtosis.warning)
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
		if( is.null(convertedData() ) ) return(NULL)
		tem = input$noIDConversion; tem=input$missingValue
		tem = readData()$dataSize
		ix = match( toupper( rownames(convertedData())), toupper(converted()$conversionTable$ensembl_gene_id  ) )
		nMatched = sum( !is.na(ix) )
		if( input$noIDConversion) 
		return( paste(tem[1], "genes in", tem[4], "samples.", tem[3], " genes passed filter (see above). Original gene IDs used." ) )  else
		return( paste(tem[1], "genes in", tem[4], "samples. Of the", tem[3], " genes passed filter (see above), ", 
			nMatched, " were converted to Ensembl gene IDs in our database. 
			  The remaining ", tem[3]-nMatched, "genes were kept in the data using original IDs." ) ) 	  
			  
			  
			  
	})	

	# Show info on file format	
output$fileFormat <- renderUI({
  shinyjs::hideElement(id = 'loadMessage')
		i = "<h4>Ready to load data files.</h3>"
#		i = c(i,"Users can upload a CSV or tab-delimited text file with the first column as gene IDs. 
#		For RNA-seq data, read count per gene is recommended.
#		Also accepted are normalized expression data based on FPKM, RPKM, or DNA microarray data. iDEP can convert most types of common gene IDs to Ensembl gene IDs, which is used 
#			internally for enrichment and pathway analyses. iDEP parses column names to define sample groups. To define 3 biological samples (Control,
#		TreatmentA, TreatmentB) with 2 replicates each, column names should be:")
#		i = c(i," <strong> Ctrl_1, Ctrl_2, TrtA_1, TrtA_2, TrtB_1, TrtB_2</strong>.") 
#		i = c(i,"For more complex experimental design, users can upload a <a href=\"https://idepsite.wordpress.com/data-format/\" target=\"_blank\">sample information file</a>  with samples in columns and factors (genotypes and conditions) in rows. 
#		       With such a file, users can define a statistic model according to study design, which enables them to control the effect for batch effects or paired samples. 
#	         or detect interactions between factors (how mutant responds differently to treatment than wild-type).") 
		
		HTML(paste(i, collapse='<br/>') )
	})

converted <- reactive({
		if (is.null(input$file1) && input$goButton == 0)    return(NULL)
		tem = input$selectOrg;
		isolate( {
		
		convertID(rownames(readData()$data ),input$selectOrg);

		# converted()$conversionTable: Not matched is skipped
		#User_input	ensembl_gene_id	Species
		#MTURN	ENSMUSG00000038065	Mouse
		#MTUS1	ENSMUSG00000045636	Mouse


		}) 
	})


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
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
		####################################
		if( is.null(converted() ) ) return( readData()$data) # if id or species is not recognized use original data.
		isolate( {  
			withProgress(message="Converting data ... ", {
			
				if(input$noIDConversion) return( readData()$data )
				
				mapping <- converted()$conversionTable
				# cat (paste( "\nData:",input$selectOrg) )
				x =readData()$data

				rownames(x) = toupper(rownames(x))
				# any gene not recognized by the database is disregarded
				# x1 = merge(mapping[,1:2],x,  by.y = 'row.names', by.x = 'User_input')
				# the 3 lines keeps the unrecogized genes using original IDs
				x1 = merge(mapping[,1:2],x,  by.y = 'row.names', by.x = 'User_input', all.y=TRUE)

				# original IDs used if ID is not matched in database
				ix = which(is.na(x1[,2]) )
				x1[ix,2] = x1[ix,1] 
				
				#multiple matched IDs, use the one with highest SD
				tem = apply(x1[,3:(dim(x1)[2])],1,sd)
				x1 = x1[order(x1[,2],-tem),]
				x1 = x1[!duplicated(x1[,2]) ,]
				rownames(x1) = x1[,2]
				x1 = as.matrix(x1[,c(-1,-2)])
				tem = apply(x1,1,sd)
				x1 = x1[order(-tem),]  # sort again by SD
				incProgress(1, "Done.")
			
				return(x1)
		})
		})
	})
	
convertedCounts <- reactive({
		if (is.null(input$file1) && input$goButton == 0) return()  
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2 }
		####################################
		if( is.null(converted() ) ) return( readData()$rawCounts) # if id or species is not recognized use original data.
		isolate( {  
			if(input$noIDConversion) return( readData()$rawCounts )
			withProgress(message="Converting data ... ", {

				
				mapping <- converted()$conversionTable
				# cat (paste( "\nData:",input$selectOrg) )
				x =readData()$rawCounts
				if(is.null(x)) return(NULL) else 
				{ 
					rownames(x) = toupper(rownames(x))
					# any gene not recognized by the database is disregarded
					# x1 = merge(mapping[,1:2],x,  by.y = 'row.names', by.x = 'User_input')
					# the 3 lines keeps the unrecogized genes using original IDs
					x1 = merge(mapping[,1:2],x,  by.y = 'row.names', by.x = 'User_input', all.y=TRUE)
					ix = which(is.na(x1[,2]) )
					x1[ix,2] = x1[ix,1] # original IDs used
					
					#multiple matched IDs, use the one with highest SD
					tem = apply(x1[,3:(dim(x1)[2])],1,sd)
					x1 = x1[order(x1[,2],-tem),]
					x1 = x1[!duplicated(x1[,2]) ,]
					rownames(x1) = x1[,2]
					x1 = as.matrix(x1[,c(-1,-2)])
					tem = apply(x1,1,sd)
					x1 = x1[order(-tem),]  # sort again by SD
					incProgress(1, "Done.")
					return(x1)
				} # here
			})
		
		})
	})

convertedPvals <- reactive({
		if (is.null(input$file1) && input$goButton == 0) return()  
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
		####################################
		if( is.null(converted() ) ) return( readData()$pvals) # if id or species is not recognized use original data.
		if( is.null(readData()$pvals) ) return(NULL)
		isolate( {  
			withProgress(message="Converting data ... ", {
			
				if(input$noIDConversion) return( readData()$pvals )
				
				mapping <- converted()$conversionTable
				# cat (paste( "\nData:",input$selectOrg) )
				x =readData()$pvals

				rownames(x) = toupper(rownames(x))
				# any gene not recognized by the database is disregarded
				# x1 = merge(mapping[,1:2],x,  by.y = 'row.names', by.x = 'User_input')
				# the 3 lines keeps the unrecogized genes using original IDs
				x1 = merge(mapping[,1:2],x,  by.y = 'row.names', by.x = 'User_input', all.y=TRUE)

				# original IDs used if ID is not matched in database
				ix = which(is.na(x1[,2]) )
				x1[ix,2] = x1[ix,1] 
				
				#multiple matched IDs, use the one with highest SD
				tem = apply(x1[,3:(dim(x1)[2])],1,sd)
				x1 = x1[order(x1[,2],-tem),]
				x1 = x1[!duplicated(x1[,2]) ,]
				rownames(x1) = x1[,2]
				x1 = as.matrix(x1[,c(-1,-2)])
				tem = apply(x1,1,sd)
				x1 = x1[order(-tem),]  # sort again by SD
				incProgress(1, "Done.")
			
				return(x1)
		})
		})
	})

GeneSetsPCA <- reactive({
		if (is.null(input$file1) && input$goButton == 0)	return()
		tem = input$selectOrg
		tem = input$selectGO6
		tem =input$minSetSize; tem= input$maxSetSize; 
		isolate( {
			if(input$selectOrg == "NEW" && !is.null(input$gmtFile) ) # new species 
			{     inFile <- input$gmtFile; inFile <- inFile$datapath
				return( readGMTRobust(inFile) ) }else
		return( readGeneSets( converted(), convertedData(), input$selectGO6,input$selectOrg,c(input$minSetSize, input$maxSetSize)  ) ) }) 
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
# show first 20 rows of data
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
	  },include.rownames=FALSE,striped=TRUE,bordered = TRUE, width = "auto",hover=T)


output$species <-renderTable({   
      if (is.null(input$file1) && input$goButton == 0)    return()
      isolate( {  #tem <- convertID(input$input_text,input$selectOrg );
	  	  withProgress(message=sample(quotes,1), detail="Converting gene IDs", {
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
  

# show first 20 rows of processed data; not used
output$debug <- renderTable({
      if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	  tem = input$selectOrg; tem = input$lowFilter ; tem =input$NminSamples2; tem = input$transform
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
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2 }
    sample1 <- as.integer( input$scatterX); sample2 <- as.integer( input$scatterY );
	####################################

	
    x <- readData()$data
	memo =""
	if( ncol(x) > maxSamplesEDAplot ) { # if too many samples, just show the first 40 or 60
		#part= sample(1:ncol(x), maxSamplesEDAplot)
		part = 1:maxSamplesEDAplot
		x <- x[, part ]
		memo =paste("(only showing", maxSamplesEDAplot, "samples)")
		}
	groups = as.factor( detectGroups(colnames(x ), readSampleInfo() ) )
	if(nlevels(groups)<=1 | nlevels(groups) >20 )  
	   col1 = "green"  else
	   col1 = rainbow(nlevels(groups))[ groups ]	
	   
	par(mfrow=c(3,1))
	par(mar=c(18,8,4,4))
	myColors = rainbow(dim(x)[2])
	
	if(ncol(x) < 31 )  
	   cexFactor = 2.2 else
	   cexFactor =1.6
    
	plot(x[,c(sample1,sample2)],xlab=colnames(x)[sample1],ylab=colnames(x)[sample2], 
         main="Scatter plot of transformed expression in two samples",
	     cex.lab=2.2, cex.axis=2.2, cex.main=2, cex.sub=2)
    text( x = min(x[, sample1]) + .75 * ( max( x[, sample1] ) - min( x[, sample1]) ),
          y = min(x[, sample2]) + .25 * ( max( x[, sample2] ) - min( x[, sample2]) ),
          paste0("R = ", round(cor(x[, sample1], x[, sample2]),3 )), 
          cex = 3
        )
              
	#------------------------boxplot
	boxplot(x, las = 2,  main=paste("Distribution of transformed data",memo)
		,cex.lab=2.2,  cex.axis=cexFactor, cex.main=2, cex.sub=2,col=col1)	
	
	#----------------------- density plot
     maxDensity = max( apply(x,2, function(y) max(density(y)$y ) ) )	
	plot(density(x[,1]),col = col1[1], lwd=2,
	  xlab="Expression values", main= paste("Density plot of transformed data",memo),
	  cex.lab=2.2, cex.axis=2.2, cex.main=2, cex.sub=2, ylim=c(0, maxDensity+0.01 )  )  #ylim=c(0,1)
	  
	for( i in 2:dim(x)[2] )
	lines(density(x[,i]),col=col1[i],  lwd=1 )
	if(nlevels(groups)< 31 ) # if too many samples do not show legends
		legend("topright", levels(groups), lty=rep(1,nlevels(groups)), col=rainbow(nlevels(groups)), cex = 2 )	

   


	
   }, height = 1600)
EDA4download <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2 }
    sample1 <- as.integer( input$scatterX); sample2 <- as.integer( input$scatterY );
	####################################
    
    # total read counts plots
	par(mfrow=c(4,1))
#	par(mar=c(18,8,4,4))
    par(mar=c(20,4,2,2))
    x <- readData()$rawCounts
	memo =""
	if( ncol(x) > maxSamplesEDAplot ) { 
		#part= sample(1:ncol(x), maxSamplesEDAplot)
		part=1:maxSamplesEDAplot
		x <- x[,part]
		memo =paste(" (only showing", maxSamplesEDAplot, "samples)")
	}
	groups = as.factor( detectGroups(colnames(x ), readSampleInfo() ) )
	if(nlevels(groups)<=1 | nlevels(groups) >20)  
	   col1 = "green"  else
	   col1 = rainbow(nlevels(groups))[ groups ]	
	   
	if(ncol(x) < 31 )  
	   cexFactor = 1.5 else
	   cexFactor =1


	   
	barplot( colSums(x)/1e6, col=col1,las=3, 
		cex.axis=1.5,    # expansion factor for numeric axis labels.
		cex.names=cexFactor,  # expansion factor for axis names (bar labels).
		main=paste("Total read counts (millions)", memo) ) 

    # boxplot, density plot and scatter plot
	
    x <- readData()$data
	memo =""
	if( ncol(x) > maxSamplesEDAplot ) { # if too many samples, just show the first 40 or 60
		#part= sample(1:ncol(x), maxSamplesEDAplot)
		part = 1:maxSamplesEDAplot
		x <- x[, part ]
		memo =paste("(only showing", maxSamplesEDAplot, "samples)")
		}
	groups = as.factor( detectGroups(colnames(x ), readSampleInfo() ) )
	if(nlevels(groups)<=1 | nlevels(groups) >20 )  
	   col1 = "green"  else
	   col1 = rainbow(nlevels(groups))[ groups ]	
	   

	myColors = rainbow(dim(x)[2])
	
	if(ncol(x) < 31 )  
	   cexFactor = 2.2 else
	   cexFactor =1.6
	#------------------------boxplot
	boxplot(x, las = 2,  main=paste("Distribution of transformed data",memo)
		,cex.lab=2.2,  cex.axis=cexFactor, cex.main=2, cex.sub=2,col=col1)	
	
	#----------------------- density plot
     maxDensity = max( apply(x,2, function(y) max(density(y)$y ) ) )	
	plot(density(x[,1]),col = col1[1], lwd=2,
	  xlab="Expression values", main= paste("Density plot of transformed data",memo),
	  cex.lab=2.2, cex.axis=2.2, cex.main=2, cex.sub=2, ylim=c(0, maxDensity+0.01 )  )  #ylim=c(0,1)
	  
	for( i in 2:dim(x)[2] )
	lines(density(x[,i]),col=col1[i],  lwd=1 )
	if(nlevels(groups)< 31 ) # if too many samples do not show legends
		legend("topright", levels(groups), lty=rep(1,nlevels(groups)), col=rainbow(nlevels(groups)), cex = 2 )	

	plot(x[,c(sample1,sample2)],xlab=colnames(x)[sample1],ylab=colnames(x)[sample2], 
         main="Scatter plot of transformed expression in two samples",
	     cex.lab=2.2, cex.axis=2.2, cex.main=2, cex.sub=2)
    text( x = min(x[, sample1]) + .75 * ( max( x[, sample1] ) - min( x[, sample1]) ),
          y = min(x[, sample2]) + .25 * ( max( x[, sample2] ) - min( x[, sample2]) ),
          paste0("R = ", round(cor(x[, sample1], x[, sample2]),3 )), 
          cex = 3
        )
	
   })
output$downloadEDAplot <- downloadHandler(
      filename = "EDA.eps",
      content = function(file) {
	  cairo_ps(file, width = 6, height = 24, pointsize = 9)
        EDA4download()
      dev.off()
      })     
output$genePlot <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat ; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}

	
	tem = input$geneSearch ; tem = input$genePlotBox; tem = input$useSD
	isolate({
	  
	  plotGenes( convertedData(), allGeneInfo(), readSampleInfo(), input$geneSearch, input$genePlotBox, input$useSD, input$selectOrg)

	})
   })

genePlot4Download <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat ; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}

	
	tem = input$geneSearch ; tem = input$genePlotBox; tem = input$useSD
	isolate({
	  
	  plotGenes( convertedData(), allGeneInfo(), readSampleInfo(), input$geneSearch, input$genePlotBox, input$useSD, input$selectOrg)
	  
	})
   })

output$downloadGenePlot <- downloadHandler(
      filename = "genePlot.eps",
      content = function(file) {
	cairo_ps(file, width = 8, height = 6, points = 8 )
        genePlot4Download()
        dev.off()
      })     

processedData <- reactive({
		  if (is.null(input$file1) && input$goButton == 0)    return()
		  
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
		####################################
            # sometimes users upload unknow species but not choosing "NEW". 
			if(input$selectOrg == "NEW" | ncol(allGeneInfo()) == 1 ) return(  convertedData() ) else { 

				withProgress(message="Preparing data for download ", {
				
				tem <-  merge(allGeneInfo()[,c('ensembl_gene_id','symbol')], round(convertedData(),4),by.x="ensembl_gene_id", by.y ="row.names", all.y=T )   
				tem[,2] = paste(" ",tem[,2]) # add space to gene symbol to avoid auto convertion of symbols to dates by Excel 
				incProgress(1/2, "mapping")
				#tem[,1] = paste(" ",tem[,1]) # add space
				 tem2 = merge( converted()$conversionTable, tem, by.x = "ensembl_gene_id",by.y="ensembl_gene_id", all.y=TRUE)
				 tem2 <- tem2[,-3] # remove species column
				 ix = which( tem2[,3] == "  NA") # remove NA's 
				 tem2[ix,3] <- ""
				 
				 ix = which(is.na(tem2[,2]) ) # not in mapping table
				 userIDs = tem2[,2]; userIDs[ix]=tem2[ix,1]	; userIDs = paste0(" ",userIDs)	# prevents Excel auto conversion		 
				 tem2[,2] = tem2[,1]; tem2[ix,2]= ""
				 tem2[,1]= userIDs;
				 colnames(tem2)[1:2]= c("User_ID sorted by SD","Ensembl_gene_id")
				 
				 # sort by sd
				 tem2 = tem2[ order( -apply(tem2[,-3:-1],1,sd )   )  ,]
				 
				 rownames(tem2)=1:nrow(tem2)
				 
				 incProgress(1, "Done.")
				 # add original data
				# tem3 <- merge( readData()$data, tem2, by.x = "row.names", by.y = "User_input", all.x=TRUE )
				# return(converted()$conversionTable) #return(readData()$data)
				# colnames(tem)[1] = "Ensembl_or_Original_ID"
				})
				return(tem2)
				
			
		}
	  })

processedCountsData <- reactive({
		  if (is.null(input$file1) && input$goButton == 0)    return()
		  
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
		####################################
			
		if(input$selectOrg == "NEW"| ncol(allGeneInfo()) == 1) return(  convertedData() ) else { 

			withProgress(message="Preparing data for download ", {

			if(is.null(convertedCounts() ) ) return(NULL) else 
			{

				withProgress(message="Preparing data for download ", {
				
				tem <-  merge(allGeneInfo()[,c('ensembl_gene_id','symbol')], round(convertedCounts(),4),by.x="ensembl_gene_id", by.y ="row.names", all.y=T )   
				tem[,2] = paste(" ",tem[,2]) # add space to gene symbol to avoid auto convertion of symbols to dates by Excel 
				incProgress(1/2, "mapping")
				#tem[,1] = paste(" ",tem[,1]) # add space
				 tem2 = merge( converted()$conversionTable, tem, by.x = "ensembl_gene_id",by.y="ensembl_gene_id", all.y=TRUE)
				 tem2 <- tem2[,-3] # remove species column
				 ix = which( tem2[,3] == "  NA") # remove NA's 
				 tem2[ix,3] <- ""
				 
				 ix = which(is.na(tem2[,2]) ) # not in mapping table
				 userIDs = tem2[,2]; userIDs[ix]=tem2[ix,1]	; userIDs = paste0(" ",userIDs)	# prevents Excel auto conversion		 
				 tem2[,2] = tem2[,1]; tem2[ix,2]= ""
				 tem2[,1]= userIDs;
				 colnames(tem2)[1:2]= c("User_ID sorted by SD","Ensembl_gene_id")
				 # sort by sd
				 tem2 = tem2[ order( -apply(log2(10+ tem2[,-3:-1]),1,sd )   )  ,]
				 rownames(tem2)=1:nrow(tem2)
				 incProgress(1, "Done.")
				 # add original data
				# tem3 <- merge( readData()$data, tem2, by.x = "row.names", by.y = "User_input", all.x=TRUE )
				# return(converted()$conversionTable) #return(readData()$data)
				# colnames(tem)[1] = "Ensembl_or_Original_ID"
				})
				return(tem2)
			}

			}) # progress
			
		}
	  })  

output$downloadProcessedData <- downloadHandler(
		filename = function() {"Processed_Data.csv"},
		content = function(file) {
      write.csv( processedData(), file)	    
	})

	
output$downloadConvertedCounts <- downloadHandler(
		filename = function() {"Converted_Counts_Data.csv"},
		content = function(file) {
      write.csv( processedCountsData(), file, row.names=FALSE)	    
	})
 

output$examineData <- DT::renderDataTable({
   inFile <- input$file1
	inFile <- inFile$datapath
    if (is.null(input$file1) && input$goButton == 0)   return(NULL)

	tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
	isolate({
	
	if(input$selectOrg == "NEW"| ncol(allGeneInfo()) == 1) return(  round(convertedData(),2) ) else
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
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	####################################
	
    par(mar=c(20,4,2,2))
    x <- readData()$rawCounts
	memo =""
	if( ncol(x) > maxSamplesEDAplot ) { 
		#part= sample(1:ncol(x), maxSamplesEDAplot)
		part=1:maxSamplesEDAplot
		x <- x[,part]
		memo =paste(" (only showing", maxSamplesEDAplot, "samples)")
	}
	groups = as.factor( detectGroups(colnames(x ), readSampleInfo() ) )
	if(nlevels(groups)<=1 | nlevels(groups) >20)  
	   col1 = "green"  else
	   col1 = rainbow(nlevels(groups))[ groups ]	
	   
	if(ncol(x) < 31 )  
	   cexFactor = 1.5 else
	   cexFactor =1
	   
	barplot( colSums(x)/1e6, col=col1,las=3, 
		cex.axis=1.5,    # expansion factor for numeric axis labels.
		cex.names=cexFactor,  # expansion factor for axis names (bar labels).
		main=paste("Total read counts (millions)", memo) ) 

},height =500) 

# detecting sequencing depth bias
output$readCountsBias <- renderText({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    if (is.null(readData()$rawCounts))   return(NULL)
	
	totalCounts = colSums(readData()$rawCounts) 
	groups = as.factor( detectGroups(colnames(readData()$rawCounts ), readSampleInfo() ) )

	tem = NULL
	# ANOVA of total read counts vs sample groups parsed by sample name
	pval = summary( aov(totalCounts ~ groups ))[[1]][["Pr(>F)"]][1]

	means = aggregate(totalCounts, by=list(groups), mean)
	maxMinRatio = max(means[,2])/min(means[,2])
	if(is.null(pval)) {
	  tem <- NULL
	} else if(pval <0.05)
	  tem = paste("Warning! Sequencing depth bias detected. Total read counts are significantly different among sample groups (p=",
				sprintf("%-3.2e",pval),") based on ANOVA. Total read counts max/min =",round(maxMinRatio,2))

	# ANOVA of total read counts vs factors in experiment design
	if(!is.null(readSampleInfo()   )  ) {
	  y <- readSampleInfo()
		for (j in 1:ncol(y) ) { 
		pval = summary( aov(totalCounts ~ as.factor(y[,j])))[[1]][["Pr(>F)"]][1]

		if(pval <0.01)
		tem = paste(tem, " Total read counts seem to be correlated with factor",colnames(y)[j], 
					"(p=",  sprintf("%-3.2e",pval),").  ")
	  }
	 }
	return( tem )
		  
}) 

################################################################
#   Heatmaps
################################################################

output$listFactorsHeatmap <- renderUI({
	tem = input$selectOrg; 
	tem=input$limmaPval; tem=input$limmaFC
	
    if (is.null(readSampleInfo() ) ) # if sample info is uploaded and correctly parsed.
       { return(NULL) }	 else { 
	  selectInput("selectFactorsHeatmap", label="Sample color bar:",choices= c(colnames(readSampleInfo()), "Sample_Name"))   } 
	})  

	# conventional heatmap.2 plot
output$heatmap1 <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  {  
			tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform 
		}
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) { 
			tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2
		}


	####################################

    x <- readData()$data   # x = read.csv("expression1.csv")
	withProgress(message=sample(quotes,1), detail ="Runing hierarchical clustering ", {
	n=input$nGenes
	#if(n>6000) n = 6000 # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
	if(n<10) n = 10 # min
	# this will cutoff very large values, which could skew the color 
	if(input$geneCentering)
		x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	
	# standardize by gene
	if(input$geneNormalize) 
		x <- x / apply(x,1,sd)
		
	# row centering and normalize
	x <- scale(x, center = input$sampleCentering, scale = input$sampleNormalize) 

	if(ncol(x) < 20 )  
	   cexFactor = 2 else
	if(ncol(x) < 31 ) 
	   cexFactor = 1.5 else
	cexFactor =1	
	
	cutoff = median(unlist(x)) + input$heatmapCutoff * sd (unlist(x)) 
	x[x>cutoff] <- cutoff
	cutoff = median(unlist(x)) - input$heatmapCutoff *sd (unlist(x)) 
	x[x< cutoff] <- cutoff
	
    groups = detectGroups(colnames(x) )
	# if sample info file is uploaded us that info:

	if(!is.null(readSampleInfo()) &&  !is.null(input$selectFactorsHeatmap) ) { 
		if(input$selectFactorsHeatmap == "Sample_Name" )
			groups = detectGroups(colnames(x) ) else 
			{ 	ix = match(input$selectFactorsHeatmap, colnames(readSampleInfo() ) ) 
				groups = readSampleInfo()[,ix]
			}
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
	

	if( n>110) # no labels when too many genes
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
		,cexCol=cexFactor  # size of font for sample names
		,lmat = lmat, lwid = lwid, lhei = lhei
		)

	if( n<=110) 
	heatmap.2(x, distfun =  distFuns[[as.integer(input$distFunctions)]]
		,hclustfun=hclustFuns[[as.integer(input$hclustFunctions)]]
		,Colv=!input$noSampleClustering		
		,col= heatColors[as.integer(input$heatColors1),], density.info="none", trace="none", scale="none", keysize=.5
		,key=T, symkey=F
		#,labRow=labRow
		,ColSideColors=groups.colors[ as.factor(groups)]
		,margins=c(18,12)
		,cexRow=1
		,srtCol=45
		,cexCol=cexFactor  # size of font for sample names
		,lmat = lmat, lwid = lwid, lhei = lhei
	)
	
	if(length(unique(groups) ) <= 30 ) {  # only add legend when there is less categories
		par(lend = 1)           # square line ends for the color legend
		add_legend("topleft",
			legend = unique(groups), # category labels
			col = groups.colors[ unique(as.factor(groups))],  # color key
			lty= 1,             # line style
			lwd = 10            # line width
		)
	}
	
	incProgress(1,"Done")
	})


} , height = 900 )  #, width = 600

#heatmap for download
plotHeatmap1 <- reactive ({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  {  
			tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform 
		}
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) { 
			tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2
		}
	####################################
		
    x <- readData()$data   # x = read.csv("expression1.csv")
	withProgress(message="Runing hierarchical clustering ", {
	n=input$nGenes
	#if(n>6000) n = 6000 # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
	# this will cutoff very large values, which could skew the color 
	if(input$geneCentering)
		x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	
	# standardize by gene
	if(input$geneNormalize) 
		x <- x / apply(x,1,sd)
		
	# row centering and normalize
	x <- scale(x, center = input$sampleCentering, scale = input$sampleNormalize) 
	
	cutoff = median(unlist(x)) + input$heatmapCutoff * sd (unlist(x)) 
	x[x>cutoff] <- cutoff
	cutoff = median(unlist(x)) - input$heatmapCutoff *sd (unlist(x)) 
	x[x< cutoff] <- cutoff
	
    groups = detectGroups(colnames(x) )
	# if sample info file is uploaded us that info:

	if(!is.null(readSampleInfo()) &&  !is.null(input$selectFactorsHeatmap) ) { 
		if(input$selectFactorsHeatmap == "Sample_Name" )
			groups = detectGroups(colnames(x) ) else 
			{ 	ix = match(input$selectFactorsHeatmap, colnames(readSampleInfo() ) ) 
				groups = readSampleInfo()[,ix]
			}
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
	
	if(ncol(x) < 20 )  
	   cexFactor = 2 else
	if(ncol(x) < 31 ) 
	   cexFactor = 1.5 else
	cexFactor =1	

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
		,margins=c(8,0)
		,srtCol=45
		,cexCol=cexFactor  # size of font for sample names
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
		,cexCol=cexFactor  # size of font for sample names
		,lmat = lmat, lwid = lwid, lhei = lhei
	)
	
	if(length(unique(groups) ) <= 30 ) {  # only add legend when there is less categories
	par(lend = 1)           # square line ends for the color legend
	add_legend("topleft",
		legend = unique(groups), # category labels
		col = groups.colors[ unique(as.factor(groups))],  # color key
		lty= 1,             # line style
		lwd = 10            # line width
	)
	}
	incProgress(1,"Done")
	})

}  ) 
	# conventional heatmap.2 plot

output$downloadHeatmap1 <- downloadHandler(
      filename = "heatmap.eps",
      content = function(file) {
	cairo_ps(file, width = 10, height = 15)
        plotHeatmap1()
        dev.off()
      })  

# interactive heatmap with plotly
output$heatmapPlotly <- renderPlotly({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  {  
			tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform 
		}
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) { 
			tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2
		}
	####################################
	   x <- convertedData()
	withProgress(message="Rendering heatmap ", {
	n=input$nGenesPlotly
	#if(n>6000) n = 6000 # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
	if(n<10) n = 10 # min
	
	if(input$geneCentering)
		x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	# standardize by gene
	if(input$geneNormalize) 
		x <- x / apply(x,1,sd)
	# row centering and normalize
	x <- scale(x, center = input$sampleCentering, scale = input$sampleNormalize) 

	cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
	x[x>cutoff] <- cutoff
	cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
	x[x< cutoff] <- cutoff

	# adding gene symbol
	ix <- match( rownames(x), allGeneInfo()[,1])
	geneSymbols <- as.character( allGeneInfo()$symbol)[ix]
	# if missing or duplicated, use Ensembl ID
	ix <- which(nchar( geneSymbols) <=1 | duplicated(geneSymbols ) | is.null(geneSymbols)  | is.na(geneSymbols)    );	geneSymbols[ ix ] <- rownames(x)[ix]
	rownames( x) = geneSymbols;

	incProgress(1/2, "Clustering of genes")	
	
    # clustering genes------
	clust <- x %>% 
	  dist2() %>% 
	  hclust2()
	# Get order
	ord_row <- clust$order

	#clustering samples --------
	if( input$noSampleClustering )
	  ord_column = 1:ncol(x) else { 
		clust <- t(x) %>% 
		  dist2() %>% 
		  hclust2()
		# Get order
		ord_column <- clust$order	
	}
	
	# Re-arrange based on order
	df <- t( x[ord_row,ord_column] )%>%
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

correlationMatrixData <- reactive({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		# heatmap of correlation matrix
		x <- readData()$data
		maxGene <- apply(x,1,max)
		x <- x[which(maxGene > quantile(maxGene)[1] ) ,] # remove bottom 25% lowly expressed genes, which inflate the PPC
		
	    round(cor(x),3)
})
output$downloadCorrelationMatrix <- downloadHandler(
		filename = function() {"correlationMatrix.csv"},
		content = function(file) {
			write.csv(correlationMatrixData(), file)
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
			name="Pearson's \nCorrelation") +
			theme_minimal()+ # minimal theme
			theme(axis.text.x = element_text(angle = 45, vjust = 1, size=14,hjust = 1))+
			theme(axis.text.y = element_text( size = 14 ))+
			coord_fixed()
		# print(ggheatmap)
		 if(input$labelPCC && ncol(x)<20)
				ggheatmap <- ggheatmap +  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4)
				
		ggheatmap + 
		  theme(axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				panel.grid.major = element_blank(),
				panel.border = element_blank(),
				panel.background = element_blank(),
				axis.ticks = element_blank(),
				legend.justification = c(1, 0),
				legend.position = c(0.6, 0.7),
				legend.direction = "horizontal")+
				guides(fill = FALSE) # + ggtitle("Pearson's Correlation Coefficient (all genes)")

			# why legend does not show up??????	
		 

 
  }, height = 600, width = 700  )#)

correlationMatrix4Download <- reactive({
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
			name="Pearson's \nCorrelation") +
			theme_minimal()+ # minimal theme
			theme(axis.text.x = element_text(angle = 45, vjust = 1, size=14,hjust = 1))+
			theme(axis.text.y = element_text( size = 14 ))+
			coord_fixed()
		# print(ggheatmap)
		 if(input$labelPCC && ncol(x)<20)
				ggheatmap <- ggheatmap +  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4)
				
		ggheatmap <- ggheatmap + 
		  theme(axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				panel.grid.major = element_blank(),
				panel.border = element_blank(),
				panel.background = element_blank(),
				axis.ticks = element_blank(),
				legend.justification = c(1, 0),
				legend.position = c(0.6, 0.7),
				legend.direction = "horizontal")+
				guides(fill = FALSE) # + ggtitle("Pearson's Correlation Coefficient (all genes)")
        print(ggheatmap)
			# why legend does not show up??????	
 
  })#)
output$downloadCorrelationMatrixPlot <- downloadHandler(
      filename = "correlation_matrix.eps",
      content = function(file) {
	  cairo_ps(file, width = 8, height = 8, points = 7)
        correlationMatrix4Download()
        dev.off()
      }) 
  
output$sampleTree <- renderPlot({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		# heatmap of correlation matrix
		x <- readData()$data
		maxGene <- apply(x,1,max)
		x <- x[which(maxGene > quantile(maxGene)[1] ) ,] # remove bottom 25% lowly expressed genes, which inflate the PPC
		if(input$geneCentering)
			x=as.matrix(x)-apply(x,1,mean)
		
		# standardize by gene
		if(input$geneNormalize) 
			x <- x / apply(x,1,sd)
			
		# row centering and normalize
		x <- scale(x, center = input$sampleCentering, scale = input$sampleNormalize) 
		
		#plot(as.dendrogram(hclust2( dist2(t(x)))), xlab="", ylab="1 - Pearson C.C.", type = "rectangle")
		plot(as.dendrogram(  hclustFuns[[as.integer(input$hclustFunctions)]] ( 
				distFuns[[as.integer(input$distFunctions)]](t(x)))) 
				,xlab="", ylab=paste( names(distFuns)[as.integer(input$distFunctions)],"(", 
				names(hclustFuns)[as.integer(input$hclustFunctions)],"linkage",")"   ), type = "rectangle")


# distFuns[[as.integer(input$distFunctions)]]
  } )#, height = 500, width = 500)

sampleTree4download <- reactive({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		# heatmap of correlation matrix
		x <- readData()$data
		maxGene <- apply(x,1,max)
		x <- x[which(maxGene > quantile(maxGene)[1] ) ,] # remove bottom 25% lowly expressed genes, which inflate the PPC
		if(input$geneCentering)
			x=as.matrix(x)-apply(x,1,mean)
		
		# standardize by gene
		if(input$geneNormalize) 
			x <- x / apply(x,1,sd)
			
		# row centering and normalize
		x <- scale(x, center = input$sampleCentering, scale = input$sampleNormalize) 
		
		#plot(as.dendrogram(hclust2( dist2(t(x)))), xlab="", ylab="1 - Pearson C.C.", type = "rectangle")
		plot(as.dendrogram(  hclustFuns[[as.integer(input$hclustFunctions)]] ( 
				distFuns[[as.integer(input$distFunctions)]](t(x)))) 
				,xlab="", ylab=paste( names(distFuns)[as.integer(input$distFunctions)],"(", 
				names(hclustFuns)[as.integer(input$hclustFunctions)],"linkage",")"   ), type = "rectangle")


    # distFuns[[as.integer(input$distFunctions)]]
  } )#, height = 500, width = 500)

output$downloadSampleTree <- downloadHandler(
      filename = "sample_tree.eps",
      content = function(file) {
	  cairo_ps(file, width = 8, height = 6)
        sampleTree4download()
        dev.off()
      }) 
output$distributionSD_heatmap <- renderPlot({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if( is.null(Kmeans()) ) return(NULL)

		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
				
		tem = input$nGenes
		####################################
		withProgress(message="Calculating SD distribution", {
		isolate({ 

		
		SDs=apply(convertedData(),1,sd)
		maxSD = mean(SDs)+ 4*sd(SDs)
		SDs[ SDs > maxSD] = maxSD
		
		top = input$nGenes
		if(top > length(SDs)) top = length(SDs)
		Cutoff=sort(SDs,decreasing=TRUE)[top] 

		SDs = as.data.frame(SDs)

		p <- ggplot(SDs, aes(x=SDs)) + 
		  geom_density(color="darkblue", fill="lightblue") +
		  labs(x = "Standard deviations of all genes", y="Density")+
		  geom_vline(aes(xintercept=Cutoff),
					color="red", linetype="dashed", size=1) +
					annotate("text", x = Cutoff + 0.4*sd(SDs[,1]), 
						y = 1,colour = "red", label = paste0("Top ", top)) +
					theme(axis.text=element_text(size=14),
					axis.title=element_text(size=16,face="bold"))
			incProgress(1)
		p
		})
	
	})
	
	
	 #progress 
  }, height = 600, width = 800,res=120 )

################################################################
#   PCA
################################################################
output$listFactors <- renderUI({
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	
      if (is.null(readSampleInfo()) )
       { return(HTML("Upload a sample info file to customize this plot.") ) }	 else { 
	  selectInput("selectFactors", label="Color:",choices=c( colnames(readSampleInfo()), "Sample_Name")
	      , selected = "Sample_Name")   } 
	})

	
output$listFactors2 <- renderUI({
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	
      if (is.null(readSampleInfo()) )
       { return(NULL) }	 else { 
	   tem <- c( colnames(readSampleInfo()), "Sample_Name")
	   if(length(tem)>1) { tem2 = tem[1]; tem[1] <- tem[2]; tem[1] = tem2; } # swap 2nd factor with first
	  selectInput("selectFactors2", label="Shape:",choices=tem, selected = "Sample_Name")
	        } 
	})

# note this function is cloned below for download. Any changes made to it, a similar change need to be made to PCAplots4Download	
output$PCA <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	#tem = input$selectGO6
    PCAxy <- c(as.integer( input$PCAx ),as.integer( input$PCAy) ) # selected principal components
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  {  
			tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform 
		}
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) { 
			tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2
		}
	####################################
	# for showing shapes on ggplot2. The first 6 are default. Default mapping can only show 6 types.
	shapes = c(16,17,15,3,7,8,   1,2,4:6,9:15,18:25  )

	withProgress(message=sample(quotes,1), detail ="Running ", {

	x <- convertedData();
	
	#---PCA----------------------------------------------- 
     if(input$PCA_MDS ==1) {   #PCA
		 incProgress(1/3,detail="PCA")
		 pca.object <- prcomp(t(x))
		 # par(mfrow=c(2,1))
		if(0){
		 plot( pca.object$x[,1], pca.object$x[,2], pch = 1,cex = 2,col = detectGroups(colnames(x)),
			 xlim=c(min(pca.object$x[,1]),max(pca.object$x[,1])*1.5   ),
			xlab = "First principal component", ylab="Second Principal Component")
			text( pca.object$x[,1], pca.object$x[,2],  pos=4, labels =colnames(x), offset=.5, cex=.8)
			}
		pcaData = as.data.frame(pca.object$x[, PCAxy]); 
        pcaData = cbind(pcaData,detectGroups(colnames(x), readSampleInfo()) )
		colnames(pcaData) = c("PC1", "PC2", "Sample_Name")
		percentVar=round(100*summary(pca.object)$importance[2, PCAxy],0)
		if(is.null(readSampleInfo())) { 
			p=ggplot(pcaData, aes(PC1, PC2, color=Sample_Name, shape = Sample_Name)) 
			} else {
			pcaData = cbind(pcaData,readSampleInfo() )
			p=ggplot(pcaData, 
                     aes_string( "PC1", 
                                 "PC2", 
                                 color=input$selectFactors,
                                 shape=input$selectFactors2))  
			}
		 if(ncol(x)<20) # change size depending of # samples
			p <- p + geom_point(size=5)  else if(ncol(x)<50)
			 p <- p + geom_point(size=3)  else 
			 p <- p + geom_point(size=2)
			 

		p <- p+	 scale_shape_manual(values= shapes)	 
		
		p=p+xlab(paste0("PC", input$PCAx ,": ", percentVar[1],"% variance")) 
		p=p+ylab(paste0("PC", input$PCAy ,": ",percentVar[2],"% variance")) 
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

	#---PGSEA----------------------------------------------- 	
	if(input$PCA_MDS ==2) {  # pathway
		incProgress(1/8, detail="PGSEA")
		library(PGSEA,verbose=FALSE)
		pca.object <- prcomp(t(x))
		pca = 100*pca.object$rotation 
		Npca = 5
		Nterms = 6
		if (Npca > dim(pca)[2]) { Npca = dim(pca)[2] } else pca <-  pca[,1:Npca]
		#pca = pca[,1:5]
		if(is.null(GeneSetsPCA() ) ) return(NULL)  # no species recognized
		if(length(GeneSetsPCA() ) <= 1 ) return(NULL)
		#cat("\n\nGene Sets:",length( GeneSets()))
		pg = myPGSEA (pca,cl=GeneSetsPCA(),range=c(15,2000),p.value=TRUE, weighted=FALSE,nPermutation=1)
		incProgress(2/8)
		# correcting for multiple testing
		p.matrix = pg$p.result
		tem = p.adjust(as.numeric(p.matrix),"fdr")
		p.matrix = matrix(tem, nrow=dim(p.matrix)[1], ncol = dim(p.matrix)[2] )
		rownames(p.matrix) = rownames(pg$p.result); colnames(p.matrix) = colnames(pg$p.result)


		selected =c()
		for( i in 1:dim(p.matrix)[2]) {
		  tem = which( rank(p.matrix[,i],ties.method='first') <= Nterms)  # rank by P value
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
		incProgress(3/8, detail ="Running PGSEA")
		#tem = t(tem); tem = t( (tem - apply(tem,1,mean)) ) #/apply(tem,1,sd) )

		smcPlot(tem,scale =  c(-max(tem), max(tem)), show.grid = T, margins = c(1,1, 3, 35),
			col = .rwb,cex.lab=0.8)

	}
	 
	#---MDS----------------------------------------------- 
	if(input$PCA_MDS ==3) {  # MDS
		 fit = cmdscale( dist2(t(x) ), eig=T, k=2)
		 incProgress(1/3,detail = " MDS")
		# par(pin=c(5,5))
		if(0) {
		plot( fit$points[,1],fit$points[,2],pch = 1,cex = 2,col = detectGroups(colnames(x), readSampleInfo()),
			 xlim=c(min(fit$points[,1]),max(fit$points[,1])*1.5   ),
		  xlab = "First dimension", ylab="Second dimension"  )
		 text( fit$points[,1], fit$points[,2],  pos=4, labels =colnames(x), offset=.5, cex=1)
		}
		pcaData = as.data.frame(fit$points[,1:2]); pcaData = cbind(pcaData,detectGroups(colnames(x), readSampleInfo()) )
		colnames(pcaData) = c("x1", "x2", "Sample_Name")
		

		if(is.null(readSampleInfo())) { 
		p=ggplot(pcaData, aes(x1, x2, color=Sample_Name, shape = Sample_Name))  
		} else {
			pcaData = cbind(pcaData,readSampleInfo() )
			p=ggplot(pcaData, aes_string("x1", "x2", color=input$selectFactors,shape=input$selectFactors2))  
			}
			
		if(ncol(x)<20) # change size depending of # samples
			p <- p + geom_point(size=5)  else if(ncol(x)<50)
			 p <- p + geom_point(size=3)  else 
			 p <- p + geom_point(size=2)
		p <- p+	 scale_shape_manual(values= shapes)	 
		
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

	#--t-SNE----------------------------------------------------------------
	if(input$PCA_MDS ==4) {  # t-SNE
		incProgress(1/4, detail=" t-SNE")
		library(Rtsne,verbose=FALSE)
		 set.seed(input$tsneSeed2)
		 tsne <- Rtsne(t(x), dims = 2, perplexity=1, verbose=FALSE, max_iter = 400)

		pcaData = as.data.frame(tsne$Y); pcaData = cbind(pcaData,detectGroups(colnames(x), readSampleInfo()) )
		colnames(pcaData) = c("x1", "x2", "Sample_Name")
		
		#pcaData$Sample_Name = as.factor( pcaData$Sample_Name)

		if(is.null(readSampleInfo())) { 
			p=ggplot(pcaData, aes(x1, x2, color=Sample_Name, shape = Sample_Name)) 
		} else {
			pcaData = cbind(pcaData,readSampleInfo() )
			p=ggplot(pcaData, aes_string("x1", "x2", color=input$selectFactors,shape=input$selectFactors2)) 
			}
			
		if(ncol(x)<20) # change size depending of # samples
			p <- p + geom_point(size=5)  else if(ncol(x)<50)
			 p <- p + geom_point(size=3)  else 
			 p <- p + geom_point(size=2)
		p <- p+	 scale_shape_manual(values= shapes)	  
		p=p+xlab("Dimension 1") 
		p=p+ylab("Dimension 2") 
		p=p+ggtitle("t-SNE plot")+ coord_fixed(ratio=1.)+ 
		 theme(plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio=1) +
			 theme(axis.text.x = element_text( size = 16),
			   axis.text.y = element_text( size = 16),
			   axis.title.x = element_text( size = 16),
			   axis.title.y = element_text( size = 16) ) +
		theme(legend.text=element_text(size=16))
		   print(p)
		
	 }
	 incProgress(1)
		 }) # progress

  }, height = 700, width = 700)

PCAplots4Download <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	#tem = input$selectGO6
    PCAxy <- c(as.integer( input$PCAx ),as.integer( input$PCAy) ) # selected principal componen	
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  {  
			tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform 
		}
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) { 
			tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2
		}
	####################################
	# for showing shapes on ggplot2. The first 6 are default. Default mapping can only show 6 types.
	shapes = c(16,17,15,3,7,8,   1,2,4:6,9:15,18:25  )

	withProgress(message="Running ", {

	x <- convertedData();
	
	#---PCA----------------------------------------------- 
     if(input$PCA_MDS ==1) {   #PCA
		 incProgress(1/3,detail="PCA")
		 pca.object <- prcomp(t(x))
		 # par(mfrow=c(2,1))
		if(0){
		 plot( pca.object$x[,1], pca.object$x[,2], pch = 1,cex = 2,col = detectGroups(colnames(x), readSampleInfo()),
			 xlim=c(min(pca.object$x[,1]),max(pca.object$x[,1])*1.5   ),
			xlab = "First principal component", ylab="Second Principal Component")
			text( pca.object$x[,1], pca.object$x[,2],  pos=4, labels =colnames(x), offset=.5, cex=.8)
			}
			
		pcaData = as.data.frame(pca.object$x[, PCAxy]); pcaData = cbind(pcaData,detectGroups(colnames(x), readSampleInfo()) )
		colnames(pcaData) = c("PC1", "PC2", "Sample_Name")
		percentVar=round(100*summary(pca.object)$importance[2, PCAxy],0)
		if(is.null(readSampleInfo())) { 
			p=ggplot(pcaData, aes(PC1, PC2, color=Sample_Name, shape = Sample_Name)) 
			} else {
			pcaData = cbind(pcaData,readSampleInfo() )
			p=ggplot(pcaData, aes_string("PC1", "PC2", color=input$selectFactors,shape=input$selectFactors2))  

			}
		 if(ncol(x)<20) # change size depending of # samples
			p <- p + geom_point(size=5)  else if(ncol(x)<50)
			 p <- p + geom_point(size=3)  else 
			 p <- p + geom_point(size=2)
			 

		p <- p+	 scale_shape_manual(values= shapes)	 
		
		p=p+xlab(paste0("PC", input$PCAx ,": ", percentVar[1], "% variance")) 
		p=p+ylab(paste0("PC", input$PCAy ,": ", percentVar[2], "% variance")) 
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

	#---PGSEA----------------------------------------------- 	
	if(input$PCA_MDS ==2) {  # pathway
		incProgress(1/8, detail="PGSEA")
		library(PGSEA,verbose=FALSE)
		pca.object <- prcomp(t(x))
		pca = 100*pca.object$rotation 
		Npca = 5
		Nterms = 6
		if (Npca > dim(pca)[2]) { Npca = dim(pca)[2] } else pca <-  pca[,1:Npca]
		#pca = pca[,1:5]
		if(is.null(GeneSetsPCA() ) ) return(NULL)  # no species recognized
		if(length(GeneSetsPCA() ) <= 1 ) return(NULL)
		#cat("\n\nGene Sets:",length( GeneSets()))
		pg = myPGSEA (pca,cl=GeneSetsPCA(),range=c(15,2000),p.value=TRUE, weighted=FALSE,nPermutation=1)
		incProgress(2/8)
		# correcting for multiple testing
		p.matrix = pg$p.result
		tem = p.adjust(as.numeric(p.matrix),"fdr")
		p.matrix = matrix(tem, nrow=dim(p.matrix)[1], ncol = dim(p.matrix)[2] )
		rownames(p.matrix) = rownames(pg$p.result); colnames(p.matrix) = colnames(pg$p.result)


		selected =c()
		for( i in 1:dim(p.matrix)[2]) {
		  tem = which( rank(p.matrix[,i],ties.method='first') <= Nterms)  # rank by P value
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
		incProgress(3/8, detail ="Running PGSEA")
		#tem = t(tem); tem = t( (tem - apply(tem,1,mean)) ) #/apply(tem,1,sd) )

		smcPlot(tem,scale =  c(-max(tem), max(tem)), show.grid = T, margins = c(1,1, 3, 30),
			col = .rwb,cex.lab=0.8)

	}
	 
	#---MDS----------------------------------------------- 
	if(input$PCA_MDS ==3) {  # MDS
		 fit = cmdscale( dist2(t(x) ), eig=T, k=2)
		 incProgress(1/3,detail = " MDS")
		# par(pin=c(5,5))
		if(0) {
		plot( fit$points[,1],fit$points[,2],pch = 1,cex = 2,col = detectGroups(colnames(x)),
			 xlim=c(min(fit$points[,1]),max(fit$points[,1])*1.5   ),
		  xlab = "First dimension", ylab="Second dimension"  )
		 text( fit$points[,1], fit$points[,2],  pos=4, labels =colnames(x), offset=.5, cex=1)
		}
		pcaData = as.data.frame(fit$points[,1:2]); pcaData = cbind(pcaData,detectGroups(colnames(x), readSampleInfo()) )
		colnames(pcaData) = c("x1", "x2", "Sample_Name")
		

		if(is.null(readSampleInfo())) { 
		p=ggplot(pcaData, aes(x1, x2, color=Sample_Name, shape = Sample_Name))  
		} else {
			pcaData = cbind(pcaData,readSampleInfo() )
			p=ggplot(pcaData, aes_string("x1", "x2", color=input$selectFactors,shape=input$selectFactors2))  
			}
			
		if(ncol(x)<20) # change size depending of # samples
			p <- p + geom_point(size=5)  else if(ncol(x)<50)
			 p <- p + geom_point(size=3)  else 
			 p <- p + geom_point(size=2)
		p <- p+	 scale_shape_manual(values= shapes)	 
		
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

	#--t-SNE----------------------------------------------------------------
	if(input$PCA_MDS ==4) {  # t-SNE
		incProgress(1/4, detail=" t-SNE")
		library(Rtsne,verbose=FALSE)
		 set.seed(input$tsneSeed2)
		 tsne <- Rtsne(t(x), dims = 2, perplexity=1, verbose=FALSE, max_iter = 400)

		pcaData = as.data.frame(tsne$Y); pcaData = cbind(pcaData,detectGroups(colnames(x), readSampleInfo()) )
		colnames(pcaData) = c("x1", "x2", "Sample_Name")
		
		#pcaData$Sample_Name = as.factor( pcaData$Sample_Name)

		if(is.null(readSampleInfo())) { 
			p=ggplot(pcaData, aes(x1, x2, color=Sample_Name, shape = Sample_Name)) 
		} else {
			pcaData = cbind(pcaData,readSampleInfo() )
			p=ggplot(pcaData, aes_string("x1", "x2", color=input$selectFactors,shape=input$selectFactors2)) 
			}
			
		if(ncol(x)<20) # change size depending of # samples
			p <- p + geom_point(size=5)  else if(ncol(x)<50)
			 p <- p + geom_point(size=3)  else 
			 p <- p + geom_point(size=2)
		p <- p+	 scale_shape_manual(values= shapes)	  
		p=p+xlab("Dimension 1") 
		p=p+ylab("Dimension 2") 
		p=p+ggtitle("t-SNE plot")+ coord_fixed(ratio=1.)+ 
		 theme(plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio=1) +
			 theme(axis.text.x = element_text( size = 16),
			   axis.text.y = element_text( size = 16),
			   axis.title.x = element_text( size = 16),
			   axis.title.y = element_text( size = 16) ) +
		theme(legend.text=element_text(size=16))
		   print(p)
		
	 }
	 incProgress(1)
		 }) # progress

  })

 

output$downloadPCA <- downloadHandler(
      filename = "PCA_MDS_tSNE.eps",
      content = function(file) {
	  cairo_ps(file, width = 6, height = 6)
	  PCAplots4Download()
        dev.off()
      })    
  
PCAdata <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    x <- readData()$data
	
	 result = prcomp(t(x))$x[,1:5]
	 fit = cmdscale( dist2(t(x) ), eig=T, k=2)
     result = cbind( result, fit$points[,1:2] )
	 library(Rtsne,verbose=FALSE)
	 set.seed(input$tsneSeed2)
	 tsne <- Rtsne(t(x), dims = 2, perplexity=1, verbose=FALSE, max_iter = 400)
	
	result = cbind( result, tsne$Y)
	 
	 
	 colnames(result) = c("PCA.1","PCA.2","PCA.3","PCA.4","PCA.5","MDS.x", "MDS.y", "tSNE.x", "tSNE.y")
	 return( result)		  
  })

# correlation PCA with factors  
output$PCA2factor <- renderUI({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	if(is.null(readSampleInfo())) return(NULL)
	npc = 5 # number of Principal components
    x <- readData()$data
	y <- readSampleInfo()

	pca.object <- prcomp(t(x))
	pcaData = as.data.frame(pca.object$x[,1:npc]); 
	pvals = matrix(1,nrow=npc,ncol=ncol(y))
	for (i in 1:npc ){
		for (j in 1:ncol(y) )
			pvals[i,j] =summary( aov(pcaData[,i] ~ as.factor(y[,j])))[[1]][["Pr(>F)"]][1]
	}
	
	pvals = pvals * npc* ncol(y)   # correcting for multiple testing
	pvals[pvals>1] = 1

	colnames(pvals) = colnames(y)
	rownames(pvals) = paste0("PC",1:npc)
	a="<h4>Correlation between Principal Components (PCs) with factors </h4>"
	nchar0 = nchar(a)
	for ( i in 1:npc){
		j = which.min(pvals[i,])
		if(pvals[i,j]< 0.05) a=paste0(a,rownames(pvals)[i], 
					" is correlated with ", colnames(pvals)[j],
					" (p=",   sprintf("%-3.2e",pvals[i,j]),").<br>")
	}
	if(nchar(a) == nchar0 ) return(NULL) else 
	  return( HTML(a) )
		  
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
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  {  
			tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform 
		}
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) { 
			tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2
		}
	tem = input$KmeansReRun
	####################################
	
	withProgress(message=sample(quotes,1), detail ="k-means clustering", {
    x <- convertedData()
	#x <- readData()
	#par(mfrow=c(1,2))
	n=input$nGenesKNN
	if(n>maxGeneClustering) n = maxGeneClustering # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
	if(n<10) n = 10 # min
	#x1 <- x;
	#x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
	#x = 100* x[1:n,] / apply(x[1:n,],1,sum) 
	x = x[1:n,]
	if( input$kmeansNormalization == 'L1Norm')
		x = 100* x / apply(x,1,function(y) sum(abs(y))) else # L1 norm
	if( input$kmeansNormalization == 'geneMean')
		x = x - apply(x,1,mean)  else # this is causing problem??????
	if( input$kmeansNormalization == 'geneStandardization')	
		x = (x - apply(x,1,mean) ) / apply(x,1,sd)
	#colnames(x) = gsub("_.*","",colnames(x))
	set.seed(input$KmeansReRun)
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
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$heatColors1; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
	tem = input$KmeansReRun
	####################################
	
	if( is.null(Kmeans()) ) return(NULL)
	withProgress(message="Creating heatmap", {
   
	myheatmap2(Kmeans()$x-apply(Kmeans()$x,1,mean), Kmeans()$bar,1000,mycolor=input$heatColors1)
	
	incProgress(1, detail = paste("Done")) }) #progress 
  }, width=600, height = 500)

KmeansHeatmap4Download <- reactive({ # Kmeans clustering
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$heatColors1; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
	tem = input$KmeansReRun
	####################################
	
	if( is.null(Kmeans()) ) return(NULL)
	withProgress(message="Creating heatmap", {
   
	myheatmap2(Kmeans()$x-apply(Kmeans()$x,1,mean), Kmeans()$bar,1000,mycolor=input$heatColors1)
	
	incProgress(1, detail = paste("Done")) }) #progress 
  })

output$downloadKmeansHeatmap <- downloadHandler(
      filename = "Kmeans_heatmap.eps",
      content = function(file) {
	  cairo_ps(file, width = 8, height = 12)
	  KmeansHeatmap4Download()
        dev.off()
      })
  
output$KmeansNclusters <- renderPlot({ # Kmeans clustering
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	withProgress(message="k-means clustering", {
    x <- convertedData()
	#x <- readData()
	#par(mfrow=c(1,2))
	n=input$nGenesKNN
	#if(n>6000) n = 6000 # max
	if(n>dim(x)[1]) n = dim(x)[1] # max	as data
	if(n<10) n = 10 # min
	x = x[1:n,]
	if( input$kmeansNormalization == 'L1Norm')
		x = 100* x / apply(x,1,function(y) sum(abs(y))) else # L1 norm
	if( input$kmeansNormalization == 'geneMean')
		x = x - apply(x,1,mean)  else # this is causing problem??????
	if( input$kmeansNormalization == 'geneStandardization')	
		x = (x - apply(x,1,mean) ) / apply(x,1,sd)
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
	tem = input$KmeansReRun
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
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
  
output$tSNEgenePlot <- renderPlot({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if( is.null(Kmeans()) ) return(NULL)

		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
				
		tem = input$colorGenes; tem = input$seedTSNE
			tem = input$KmeansReRun
		####################################
		withProgress(message="Runing t-SNE algorithm", {
		isolate({ 
			Cluster <- Kmeans()$bar
			train <- as.data.frame( cbind(Cluster,Kmeans()$x) )

			library(Rtsne,verbose=FALSE)
			incProgress(1/3)
			train = unique(train)
			Cluster = train$Cluster	
			set.seed(input$seedTSNE)
			
			## Executing the algorithm on curated data
			tsne <- Rtsne(train[,-1], dims = 2, perplexity=30, verbose=FALSE, max_iter = 400)
			incProgress(2/3)

			nClusters = length(unique(Cluster) )
			if(input$colorGenes) {			
				plot(tsne$Y[,1], tsne$Y[,2], pch = (0:(nClusters-1))[Cluster], cex = 1.,col = mycolors[Cluster], xlab="X",ylab="Y")
				legend("topright",toupper(letters)[1:nClusters], pch = 0:(nClusters-1), col=mycolors, title="Cluster"  )
			} else
				plot(tsne$Y[,1], tsne$Y[,2],  cex = 1., xlab="X",ylab="Y")

			incProgress(1)
		})
	
	})
	
	
	 #progress 
  }, height = 800, width = 800,res=120 )

output$distributionSD <- renderPlot({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if( is.null(Kmeans()) ) return(NULL)

		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
				
		tem = input$colorGenes; tem = input$seedTSNE
			tem = input$KmeansReRun
		####################################
		withProgress(message="Calculating SD distribution", {
		isolate({ 

		
		SDs=apply(convertedData(),1,sd)
		maxSD = mean(SDs)+ 4*sd(SDs)
		SDs[ SDs > maxSD] = maxSD
		
		top = input$nGenesKNN
		if(top > length(SDs)) top = length(SDs)
		Cutoff=sort(SDs,decreasing=TRUE)[top] 

		SDs = as.data.frame(SDs)

		p <- ggplot(SDs, aes(x=SDs)) + 
		  geom_density(color="darkblue", fill="lightblue") +
		  labs(x = "Standard deviations of all genes", y="Density")+
		  geom_vline(aes(xintercept=Cutoff),
					color="red", linetype="dashed", size=1) +
					annotate("text", x = Cutoff + 0.4*sd(SDs[,1]), 
					y = 1,colour = "red", label = paste0("Top ", top))+				
					theme(axis.text=element_text(size=14),
						axis.title=element_text(size=16,face="bold"))
			incProgress(1)
		p
		})
	
	})
	
	
	 #progress 
  }, height = 600, width = 800,res=120 )

distributionSD4Download <- reactive({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if( is.null(Kmeans()) ) return(NULL)

		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
				
		tem = input$colorGenes; tem = input$seedTSNE
			tem = input$KmeansReRun
		####################################
		withProgress(message="Calculating SD distribution", {
		isolate({ 

		
		SDs=apply(convertedData(),1,sd)
		maxSD = mean(SDs)+ 4*sd(SDs)
		SDs[ SDs > maxSD] = maxSD
		
		top = input$nGenesKNN
		if(top > length(SDs)) top = length(SDs)
		Cutoff=sort(SDs,decreasing=TRUE)[top] 

		SDs = as.data.frame(SDs)

		p <- ggplot(SDs, aes(x=SDs)) + 
		  geom_density(color="darkblue", fill="lightblue") +
		  labs(x = "Standard deviations of all genes", y="Density")+
		  geom_vline(aes(xintercept=Cutoff),
					color="red", linetype="dashed", size=1) +
					annotate("text", x = Cutoff + 0.4*sd(SDs[,1]), 
					y = 1,colour = "red", label = paste0("Top ", top))+				
					theme(axis.text=element_text(size=14),
						axis.title=element_text(size=16,face="bold"))
			incProgress(1)
		print(p)
		})
	
	})
	
	
	 #progress 
  } )

output$downloadDistributionSD <- downloadHandler(
      filename = "gene_SD_distribution.eps",
      content = function(file) {
	  cairo_ps(file, width = 6, height = 4)
      distributionSD4Download()
        dev.off()
      })
output$downloadDistributionSD1 <- downloadHandler(
      filename = "gene_SD_distribution.eps",
      content = function(file) {
	  cairo_ps(file, width = 6, height = 4)
      distributionSD4Download()
      dev.off()
      })	
output$downloadDataKmeans <- downloadHandler(
		filename = function() {"Kmeans.csv"},
			content = function(file) {
      write.csv(KmeansData(), file)
	    }
	)	
	
KmeansGOdata <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectGO3
	if( is.null(input$selectGO3 ) ) return (NULL)
	if( input$selectGO3 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.") )#No matching species
   	if( is.null(Kmeans()) ) return(NULL)
	if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
	tem = input$KmeansReRun
	tem = input$nGenesKNN;
	tem = input$kmeansNormalization
	tem = input$nClusters
	tem = input$removeRedudantSets
    tem = input$useFilteredAsBackground
	####################################
	withProgress(message=sample(quotes,1), detail ="GO Enrichment", {
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
				if(input$removeRedudantSets) reduced = redudantGeneSetsRatio else reduced = FALSE

                if(input$useFilteredAsBackground) {
                   convertedDataBackground <- converted() 
                } else {
                  convertedDataBackground <- NULL
                 }


				result <- FindOverlap( converted = convertedID,
				                       gInfo = allGeneInfo(),
				                       GO = input$selectGO3,
				                       selectOrg = input$selectOrg,
				                       minFDR = minFDR,
				                       reduced = reduced,
				                       convertedDataBackground = convertedDataBackground
				                       )
			}
			if( is.null(result)) next;   # result could be NULL

			if( dim(result)[2] ==1) next;   # result could be NULL
			result$direction = toupper(letters)[i] 
			if (pp==0 ) { results <- result; pp <- 1;
			} else {
				results <- rbind(results,result)
			}
		}

		if(pp == 0) return( as.data.frame("No enrichment found."))
		results= results[,c(6,1,2,4,5)]
		colnames(results)= c("Cluster","FDR","nGenes","Pathways","Genes")
		if(min(results$FDR) > minFDR ) results = as.data.frame("No signficant enrichment found.") else
		results = results[which(results$FDR < minFDR),]
		incProgress(1, detail = paste("Done")) 
	}) #progress
	if( is.null(results) )  return ( as.matrix("No significant enrichment.") )	
	if( class(results) != "data.frame")  return ( as.matrix("No significant enrichment.") )
	if( dim(results)[2] ==1)  return ( as.matrix("No significant enrichment.") )
	colnames(results)[2] = "adj.Pval"
	#results$Genes <- as.character(results$Genes)
	#results$Cluster[which( duplicated(results$Cluster) ) ] <- ""
	results
  })

output$KmeansGO <- renderTable({	
  if(is.null(KmeansGOdata())) return(NULL)
  tem = input$removeRedudantSets
  results1 = KmeansGOdata()
  if(dim(results1)[2] == 1) return( results1) else{
	results1$adj.Pval <- sprintf("%-2.1e",as.numeric(results1$adj.Pval) )
	results1[,1] <- as.character(results1[,1])
	results1[ duplicated (results1[,1] ),1 ] <- ""  
	
    return( results1[,-5])
	 }
  }, digits = 0,spacing="s",striped=TRUE,bordered = TRUE, width = "auto", hover=T, 
    sanitize.text.function = function(x) x) # hyperText

output$downloadKmeansGO <- downloadHandler(
		filename = function() {"KmeansEnrichment.csv"},
		content = function(file) {
			write.csv( removeHypertext( KmeansGOdata() ), file, row.names=FALSE)
	    }
	) 
	
output$enrichmentPlotKmeans <- renderPlot({
    if(is.null(KmeansGOdata())) return(NULL)
	
	tem = input$selectGO3
	if( is.null(input$selectGO3 ) ) return (NULL)
	if( input$selectGO3 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.") )#No matching species
   	if( is.null(Kmeans()) ) return(NULL)
	if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
	tem = input$KmeansReRun; 
	tem = input$nGenesKNN;
	tem = input$kmeansNormalization
	tem = input$nClusters
	tem = input$removeRedudantSets
	####################################
	
	
	tem1 = removeHypertext( KmeansGOdata() )
	colnames(tem1)[1]="Direction"
	enrichmentPlot(tem1, 46  )

}, width = 800, height = 1600)


output$enrichmentPlotKmeans4Download <- downloadHandler(
      filename = "enrichmentPlotKmeans.eps",
      content = function(file) {
	  cairo_ps(file, width = 10, height = 16)
	  tem = removeHypertext( KmeansGOdata() )

	  colnames(tem)[1]="Direction"
	  enrichmentPlot(tem,41  )
        dev.off()
      })
	
output$KmeansPromoter <- renderTable({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectGO3; tem = input$radioPromoterKmeans; tem=input$nGenesKNN; tem=input$nClusters
	if( is.null(input$selectGO3 ) ) return (NULL)
	#if( is.null(limma()$results) ) return(NULL)
   	if( is.null(Kmeans()) ) return(NULL)	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$KmeansReRun
	####################################

	isolate({ 
   	withProgress(message="Promoter analysis", {


	results1 <- NULL; result <- NULL 
	pp<- 0
	for( i in 1:input$nClusters ) {
	incProgress(1/input$nClusters, , detail = paste("Cluster",toupper(letters)[i]) )
	#query = rownames(x)[which(bar == i)]
	query = rownames(Kmeans()$x)[which(Kmeans()$bar == i)]	
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
#   Differential gene expression  1
################################################################
output$listFactorsDE <- renderUI({
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	 # Note that " in HTML needs to be denoted with \"    with  the escape character \.
    if (is.null(readSampleInfo()) ) {# if sample info is uploaded and correctly parsed.
        return(HTML("<font size = \"2\">A <a href=\"https://idepsite.wordpress.com/data-format/\">sample information file</a> 
	     can be uploaded to build a linear model according to experiment design. </font>")) 
	 } else {
	#} else if ( !(input$dataFileFormat==1& input$CountsDEGMethod==3 )  ){  # disable factor choosing for DESeq2	   
		factors = colnames(readSampleInfo())
		choices = setNames(factors, factors  )
		#interactions = apply(t(combn(factors,2)),1, function(x) paste(x,collapse=":"))
		#choices = append( choices,setNames( interactions, paste(interactions,"interaction") ))
		title1 = "1. Select 1 or 2 main factors. Or leave it blank and just choose pairs of sample groups below."	
		if ( input$dataFileFormat==1& input$CountsDEGMethod==3  )
					title1 = "1. Select 6 or less main factors. Or skip this step and just choose pairs of sample groups below."	
		checkboxGroupInput("selectFactorsModel", 
                              h5(title1), 
                              choices = choices,
                              selected = NULL)	   	  
	        } 
	})  

	
output$listBlockFactorsDE <- renderUI({
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	 # Note that " in HTML needs to be denoted with \"    with  the escape character \.
    if (is.null(readSampleInfo()) ) {# if sample info is uploaded and correctly parsed.
		return(NULL)		   
		 } else { 
	  # } else if ( !(input$dataFileFormat==1& input$CountsDEGMethod==3 )  ){ 
		factors = colnames(readSampleInfo())
		
		#only show factors not in main modelFactors
		
		factors = setdiff(factors, input$selectFactorsModel  )
		if(length(factors) == 0 ) return(NULL)
		choices = setNames(factors, factors  )
		
		title1 = "Select a factor for batch effect or paired samples, if needed."	
		if ( input$dataFileFormat==1& input$CountsDEGMethod==3  )
			title1 = "Select factors for batch effects or paired samples, if needed."	
		
		checkboxGroupInput("selectBlockFactorsModel", 
                              h5(title1), 
                              choices = choices,
                              selected = NULL)  
	  
	        } 
	})	
	
	
output$listModelComparisons <- renderUI({
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if (is.null(readSampleInfo()) | is.null(input$selectFactorsModel) ) # if sample info is uploaded and correctly parsed.
		   { # if using sample names
		   
		   	factors = as.character ( detectGroups( colnames( readData()$data ) ) )
			factors = unique(factors)# order is reversed
			comparisons = apply(t(combn(factors,2)),1, function(x) paste(x,collapse=" vs. "))
			comparisons = c(comparisons, apply(t(combn(rev(factors),2)),1, function(x) paste(x,collapse=" vs. ")) )	
			comparisons = sort(comparisons)
			choices =  setNames(gsub(" vs\\. ","-",comparisons), comparisons )
			checkboxGroupInput("selectModelComprions", 
									  h5("Select comparisons among sample groups:"), 
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
					
				} # for each factor
				if(length(choices)==0 ) return(NULL) else
				checkboxGroupInput("selectModelComprions", 
									  h5("2. Select one or more comparisons:"), 
									  choices = choices,
									  selected = choices[[1]])	   
				} 
	}) 

	
output$listInteractionTerms <- renderUI({
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if (is.null(readSampleInfo()) | is.null(input$selectFactorsModel) ) # if sample info is uploaded and correctly parsed.
		   {  return(NULL)
	   
		   }	 else { 

				selectedFactors = input$selectFactorsModel
				selectedFactors = selectedFactors[ !grepl(":",selectedFactors  ) ]
				if(length(selectedFactors)<=1 ) return(NULL) 
				#if ( !(input$dataFileFormat==1 & input$CountsDEGMethod==3)   )  return(NULL) # if not using DESeq2
				# for ( i in 1:length(selectedFactors) ) {
				interactions = apply(t(combn(selectedFactors,2)),1, function(x) paste(x,collapse=":"))
				
				choices = setNames(interactions,interactions)

				 checkboxGroupInput(  'selectInteractions', 
									  label=h5("Interaction terms between factors(e.g. genotypes repond differently to treatment?):"),
									  choices= choices,selected = NULL)
									  
					
				
				
				} #else 
	}) 

	
	# set limits for selections of factors. 
observe({
		if(length(input$selectFactorsModel) > maxFactors) # less than 4 factors
			updateCheckboxGroupInput(session, "selectFactorsModel", selected= tail(input$selectFactorsModel,maxFactors))
		
		if( input$CountsDEGMethod !=3 ) { # if using the limma package
			if(length(input$selectFactorsModel) > 2) # less than 2 factors
				updateCheckboxGroupInput(session, "selectFactorsModel", selected= tail(input$selectFactorsModel,2))

			if(length(input$selectBlockFactorsModel) > 1) # less than 1 factors
				updateCheckboxGroupInput(session, "selectBlockFactorsModel", selected= tail(input$selectBlockFactorsModel,2))
				
		}
		if(length(input$selectComparisonsVenn) >5 )
					updateCheckboxGroupInput(session, "selectComparisonsVenn", selected= tail(input$selectComparisonsVenn,5))
		
	})
	
	
output$experimentDesign <- renderText({
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if (is.null(readSampleInfo()) | is.null(input$selectFactorsModel) ) # if sample info is uploaded and correctly parsed.
		   {  return(NULL)
	   
		   }	 else { 
					model = paste("Model: expression ~ ",  paste(input$selectFactorsModel, collapse=" + "  ))
					if(!is.null(input$selectBlockFactorsModel )  )
						model = paste0(model," + ", paste(input$selectBlockFactorsModel, collapse=" + " )    )					
					if(!is.null(input$selectInteractions )  )
						model = paste0(model," + ", paste(input$selectInteractions, collapse=" + " )    )
					return( model  )									
				} #else 
	})

	
output$selectReferenceLevels1 <- renderUI({
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if (is.null(readSampleInfo()) | is.null(input$selectFactorsModel) ) # if sample info is uploaded and correctly parsed.
		   {  return(NULL)
	   
		   }	 else { 

				selectedFactors = input$selectFactorsModel
				selectedFactors = selectedFactors[ !grepl(":",selectedFactors  ) ]
				if(length(selectedFactors)==0 ) return(NULL) 
				if ( !(input$dataFileFormat==1 & input$CountsDEGMethod==3)   )  return(NULL) # if not using DESeq2
				# for ( i in 1:length(selectedFactors) ) {
				 i = 1; 
				 if( is.na( match(selectedFactors[i], colnames( readSampleInfo() )   )   )   )  return(NULL)
				 factorLevels = unique(readSampleInfo()[, selectedFactors[i] ] )

				 selectInput(  paste0("referenceLevelFactor",i), 
									  label=h5(paste ("Reference/baseline level for",selectedFactors[i])),
									  choices= setNames(as.list( paste0(selectedFactors[i],":", factorLevels )), factorLevels ))
									  
					
				
				
				} #else 
	}) 

	
output$selectReferenceLevels2 <- renderUI({
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if (is.null(readSampleInfo()) | is.null(input$selectFactorsModel) ) # if sample info is uploaded and correctly parsed.
		   {  return(NULL)
	   
		   }	 else { 

				selectedFactors = input$selectFactorsModel
				selectedFactors = selectedFactors[ !grepl(":",selectedFactors  ) ]
				if(length(selectedFactors) < 2 ) return(NULL) 
				if ( !(input$dataFileFormat==1 & input$CountsDEGMethod==3)   )  return(NULL) # if not using DESeq2
				 i = 2; 
				 if( is.na( match(selectedFactors[i], colnames( readSampleInfo() )   )   )   )  return(NULL)
				 factorLevels = unique(readSampleInfo()[, selectedFactors[i] ] )

				 selectInput(  paste0("referenceLevelFactor",i), 
									  label=h5( paste ("Reference/baseline level for",selectedFactors[i])),
									  choices= setNames(as.list( paste0(selectedFactors[i],":", factorLevels )), factorLevels ))
									                             # "genotype:wt"
					
				
				
				} #else 
	})	

	
output$selectReferenceLevels3 <- renderUI({
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if (is.null(readSampleInfo()) | is.null(input$selectFactorsModel) ) # if sample info is uploaded and correctly parsed.
		   {  return(NULL)
	   
		   }	 else { 

				selectedFactors = input$selectFactorsModel
				selectedFactors = selectedFactors[ !grepl(":",selectedFactors  ) ]
				if(length(selectedFactors) < 2 ) return(NULL) 
				if ( !(input$dataFileFormat==1 & input$CountsDEGMethod==3)   )  return(NULL) # if not using DESeq2
				 i = 3; 
				 if( is.na( match(selectedFactors[i], colnames( readSampleInfo() )   )   )   )  return(NULL)
				 factorLevels = unique(readSampleInfo()[, selectedFactors[i] ] )

				 selectInput(  paste0("referenceLevelFactor",i), 
									  label=h5( paste ("Reference/baseline level for",selectedFactors[i])),
									  choices= setNames(as.list( paste0(selectedFactors[i],":", factorLevels )), factorLevels ))
									                             # "genotype:wt"
					
				
				
				} #else 
	})

	
output$selectReferenceLevels4 <- renderUI({
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if (is.null(readSampleInfo()) | is.null(input$selectFactorsModel) ) # if sample info is uploaded and correctly parsed.
		   {  return(NULL)
	   
		   }	 else { 

				selectedFactors = input$selectFactorsModel
				selectedFactors = selectedFactors[ !grepl(":",selectedFactors  ) ]
				if(length(selectedFactors) < 2 ) return(NULL) 
				if ( !(input$dataFileFormat==1 & input$CountsDEGMethod==3)   )  return(NULL) # if not using DESeq2
				 i = 4; 
				 if( is.na( match(selectedFactors[i], colnames( readSampleInfo() )   )   )   )  return(NULL)
				 factorLevels = unique(readSampleInfo()[, selectedFactors[i] ] )

				 selectInput(  paste0("referenceLevelFactor",i), 
									  label=h5( paste ("Reference/baseline level for",selectedFactors[i])),
									  choices= setNames(as.list( paste0(selectedFactors[i],":", factorLevels )), factorLevels ))
									                             # "genotype:wt"
					
				
				
				} #else 
	})	

	
output$selectReferenceLevels5 <- renderUI({
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if (is.null(readSampleInfo()) | is.null(input$selectFactorsModel) ) # if sample info is uploaded and correctly parsed.
		   {  return(NULL)
	   
		   }	 else { 

				selectedFactors = input$selectFactorsModel
				selectedFactors = selectedFactors[ !grepl(":",selectedFactors  ) ]
				if(length(selectedFactors) < 2 ) return(NULL) 
				if ( !(input$dataFileFormat==1 & input$CountsDEGMethod==3)   )  return(NULL) # if not using DESeq2
				 i = 5; 
				 if( is.na( match(selectedFactors[i], colnames( readSampleInfo() )   )   )   )  return(NULL)
				 factorLevels = unique(readSampleInfo()[, selectedFactors[i] ] )

				 selectInput(  paste0("referenceLevelFactor",i), 
									  label=h5( paste ("Reference/baseline level for",selectedFactors[i])),
									  choices= setNames(as.list( paste0(selectedFactors[i],":", factorLevels )), factorLevels ))
									                             # "genotype:wt"
					
				
				
				} #else 
	})	

	
output$selectReferenceLevels6 <- renderUI({
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if (is.null(readSampleInfo()) | is.null(input$selectFactorsModel) ) # if sample info is uploaded and correctly parsed.
		   {  return(NULL)
	   
		   }	 else { 

				selectedFactors = input$selectFactorsModel
				selectedFactors = selectedFactors[ !grepl(":",selectedFactors  ) ]
				if(length(selectedFactors) < 2 ) return(NULL) 
				if ( !(input$dataFileFormat==1 & input$CountsDEGMethod==3)   )  return(NULL) # if not using DESeq2
				 i = 6; 
				 if( is.na( match(selectedFactors[i], colnames( readSampleInfo() )   )   )   )  return(NULL)
				 factorLevels = unique(readSampleInfo()[, selectedFactors[i] ] )

				 selectInput(  paste0("referenceLevelFactor",i), 
									  label=h5( paste ("Reference/baseline level for",selectedFactors[i])),
									  choices= setNames(as.list( paste0(selectedFactors[i],":", factorLevels )), factorLevels ))
									                             # "genotype:wt"
					
				
				
				} #else 
	})	

factorReferenceLevels <- reactive({
	if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$submitModelButton  # this is used to make it responsive only when model and comparisons are completed
	
	return( c(  input$referenceLevelFactor1,
				input$referenceLevelFactor2,
				input$referenceLevelFactor3,
				input$referenceLevelFactor4,
				input$referenceLevelFactor5,
				input$referenceLevelFactor6 )
	)

})
	
limma <- reactive({  
	if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$CountsDEGMethod; tem = input$countsLogStart
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem=input$CountsDEGMethod
	tem = input$submitModelButton  # this is used to make it responsive only when model and comparisons are completed
	####################################
	
	isolate({ 
	withProgress(message=sample(quotes,1), detail ="Identifying differentially expressed genes", {
	if(input$dataFileFormat == 1 ) {  # if count data
		 if(input$CountsDEGMethod == 3 ) {    # if DESeq2 method
				# rawCounts = read.csv("exampleData/airway_GSE52778.csv", row.names=1)
				# res =DEG.DESeq2(rawCounts, .05, 2) 
				# res1 =DEG.limma(rawCounts, .1, 1.5,rawCounts, 2,3) 
					
			return(   DEG.DESeq2( convertedCounts(),input$limmaPval, input$limmaFC,
									input$selectModelComprions, readSampleInfo(),
									c(input$selectFactorsModel,input$selectInteractions), 
									input$selectBlockFactorsModel, factorReferenceLevels() )  )
			}
			if(input$CountsDEGMethod < 3 )    # voom or limma-trend
				return( DEG.limma(convertedData(), input$limmaPval, input$limmaFC,
									convertedCounts(), input$CountsDEGMethod,
									priorCounts=input$countsLogStart,input$dataFileFormat,
									input$selectModelComprions, readSampleInfo(),
									c(input$selectFactorsModel,input$selectInteractions),
									input$selectBlockFactorsModel) )
	} else if (input$dataFileFormat == 2 ){ # normalized data
	 return( DEG.limma(convertedData(), input$limmaPval, input$limmaFC,
						convertedCounts(), input$CountsDEGMethod,
						priorCounts=input$countsLogStart,input$dataFileFormat,
						input$selectModelComprions, readSampleInfo(),
						c(input$selectFactorsModel,input$selectInteractions),
						input$selectBlockFactorsModel) )
	} else {   # dataFileFormat == 3 user just uploaded fold change matrix
	
		x = convertedData()
		
		pvals = convertedPvals()
		if(!is.null(pvals) ) {
		  ix = match(rownames(x), rownames(pvals))
		  pvals = pvals[ix,]
		}


		# looks like ratio data, take log2
		if( sum(round(apply(x,2, median) + .2) == 1 ) == dim(x)[2] & min(x) > 0) 
			x = log2(x)
		
		Exp.type = "None standard data without replicates."
		all.Calls = x # fake calls
		for( i in 1: dim(all.Calls)[2]) { 
			tem <- all.Calls[,i]
			all.Calls[which( tem <= log2(input$limmaFC) & tem >=  -log2(input$limmaFC) ) ,i] = 0			
			all.Calls[which( tem >  log2(input$limmaFC) ) ,i] = 1
			all.Calls[which( tem < -log2(input$limmaFC) ) ,i] = -1		
			if(!is.null(pvals) ) 
				all.Calls[ which( pvals[,i] > input$limmaPval),i] = 0
		}
		comparisons = colnames(all.Calls)
		extractColumn <- function (i) {
			topGenes = as.data.frame( convertedData()[,i,drop=FALSE])
			if(is.null(pvals) ) topGenes$FDR = 0 else 
				topGenes$FDR = pvals[,i]# fake fdr
				
			colnames(topGenes) = c("Fold","FDR")
			return(topGenes)	
		} 
		topGenes = lapply( 1:dim( x )[2], extractColumn )
		topGenes <- setNames(topGenes, colnames(x ) )
		
		return( list(results= all.Calls, comparisons = colnames(x ), Exp.type=Exp.type, topGenes=topGenes) )
	}
		
	})
	})
	})	

	
output$textLimma <- renderText({
      if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
 	tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
	tem=input$limmaPval; tem=input$limmaFC
	tem = input$submitModelButton 
	limma()$Exp.type
	#paste( limma()$comparisons,collapse=" "     )  
	})	

    
output$vennPlot <- renderPlot({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
		tem=input$CountsDEGMethod
		tem = input$selectFactorsModel # responsive to changes in model and comparisons
		tem = input$selectModelComprions
		tem = input$submitModelButton 
		tem = input$UpDownRegulated
		if(is.null(input$selectComparisonsVenn) ) return(NULL)
		####################################
		
		isolate({ 
		
			results = limma()$results

			# split by up or down regulation
			if(input$UpDownRegulated) { 	
			    resultUp = results; 
				resultUp[resultUp < 0 ] <- 0;
				colnames(resultUp) = paste0("Up_", colnames(resultUp))
			    resultDown = results; 
				resultDown[resultDown > 0] <- 0;
				colnames(resultDown) = paste0("Down_", colnames(resultDown))				
				results <- cbind(resultUp, resultDown)
			}			
			
			ixa = c()
			for (comps in  input$selectComparisonsVenn) { 
				 if(!grepl("^I:|^I-|^Up_I:|^Up_I-|^Down_I:|^Down_I-", comps) ) {  # if not interaction term
					ix = match(comps, colnames(results)) 
				  } else {
						#mismatch in comparison names for interaction terms for DESeq2
						#I:water_Wet.genetic_Hy 	 in the selected Contrast
						#Diff-water_Wet-genetic_Hy   in column names
						tem = gsub("^I-","I:" ,colnames(results))
						tem = gsub("-","\\.",tem)
						ix = match(comps, tem) 

						if(is.na(ix) ) # this is for limma package
							ix = match(comps, colnames(results)) 						
					  }
					ixa = c(ixa,ix)
				  }
					
			results = results[,ixa,drop=FALSE] # only use selected comparisons
			if(dim(results)[2] >5) results <- results[,1:5]
			colnames(results) = gsub("^I-","I:" ,colnames(results))	
			
			vennDiagram(results,circle.col=rainbow(5), cex=c(1.,1, 0.7) ) # part of limma package

		})
    }, height = 600, width = 600)

vennPlot4Download <- reactive({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
		tem=input$CountsDEGMethod
		tem = input$selectFactorsModel # responsive to changes in model and comparisons
		tem = input$selectModelComprions
		tem = input$submitModelButton 
		tem = input$UpDownRegulated
		if(is.null(input$selectComparisonsVenn) ) return(NULL)
		####################################
		
		isolate({ 
		
			results = limma()$results

			# split by up or down regulation
			if(input$UpDownRegulated) { 	
			    resultUp = results; 
				resultUp[resultUp < 0 ] <- 0;
				colnames(resultUp) = paste0("Up_", colnames(resultUp))
			    resultDown = results; 
				resultDown[resultDown > 0] <- 0;
				colnames(resultDown) = paste0("Down_", colnames(resultDown))				
				results <- cbind(resultUp, resultDown)
			}			
			
			ixa = c()
			for (comps in  input$selectComparisonsVenn) { 
				 if(!grepl("^I:|^I-|^Up_I:|^Up_I-|^Down_I:|^Down_I-", comps) ) {  # if not interaction term
					ix = match(comps, colnames(results)) 
				  } else {
						#mismatch in comparison names for interaction terms for DESeq2
						#I:water_Wet.genetic_Hy 	 in the selected Contrast
						#Diff-water_Wet-genetic_Hy   in column names
						tem = gsub("^I-","I:" ,colnames(results))
						tem = gsub("-","\\.",tem)
						ix = match(comps, tem) 

						if(is.na(ix) ) # this is for limma package
							ix = match(comps, colnames(results)) 						
					  }
					ixa = c(ixa,ix)
				  }
					
			results = results[,ixa,drop=FALSE] # only use selected comparisons
			if(dim(results)[2] >5) results <- results[,1:5]
			colnames(results) = gsub("^I-","I:" ,colnames(results))	
			
			vennDiagram(results,circle.col=rainbow(5), cex=c(1.,1, 0.7) ) # part of limma package

		})
    })

output$DownloadVenn <- downloadHandler(
      filename = "VennDiagram.eps",
      content = function(file) {
	  cairo_ps(file, width = 6, height = 6)
	  vennPlot4Download()
      dev.off()
      })

output$sigGeneStats <- renderPlot({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		if(is.null(limma()$results) ) return(NULL)
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
		tem=input$CountsDEGMethod
		tem = input$selectFactorsModel # responsive to changes in model and comparisons
		tem = input$selectModelComprions
		tem = input$submitModelButton 

		####################################
		
		isolate({ 
		
		results = limma()$results

		 library(reshape2)
		 Up =  apply(results, 2, function(x) sum(x == 1) )
		 Down = apply(results, 2, function(x) sum(x == -1) ) 
		 stats = rbind(Up, Down)
				 
		 gg <- melt(stats)

		 colnames(gg) = c("Regulation","Comparisons","Genes")
		 
		# gg <- within(gg, Regulation <- factor(
		#			Regulation, levels=names(sort(table(Regulation), decreasing=FALSE   ) ) )  )
		 
		# gg$Regulation <- factor( gg$Regulation, levels=c("Up","Down")  )
		# gg$Regulation = as.factor(gg$Regulation)
		 #gg$Regulation = relevel(gg$Regulation,"Up")
		 
		 p= ggplot(gg, aes(x=Comparisons, y=Genes, fill=  Regulation )  )+
			 geom_bar(position="dodge", stat="identity") + coord_flip() +
			 theme(legend.position = "top") + 
			 scale_fill_manual(values=c("red", "blue")) +
			 ylab("Number of differntially expressed genes") +
			 theme(axis.title.y=element_blank(),
				axis.text=element_text(size=14)) +  
			 geom_text(aes(label=Genes), position=position_dodge(width=0.9), vjust=0.5, hjust =0)
# updated 2020 
		p
			

		})
    })

sigGeneStats4Download <- reactive({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		if(is.null(limma()$results) ) return(NULL)
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
		tem=input$CountsDEGMethod
		tem = input$selectFactorsModel # responsive to changes in model and comparisons
		tem = input$selectModelComprions
		tem = input$submitModelButton 

		####################################
		
		isolate({ 
		
		results = limma()$results

		 library(reshape2)
		 Up =  apply(results, 2, function(x) sum(x == 1) )
		 Down = apply(results, 2, function(x) sum(x == -1) ) 
		 stats = rbind(Up, Down)
				 
		 gg <- melt(stats)

		 colnames(gg) = c("Regulation","Comparisons","Genes")
		 
		# gg <- within(gg, Regulation <- factor(
		#			Regulation, levels=names(sort(table(Regulation), decreasing=FALSE   ) ) )  )
		 
		# gg$Regulation <- factor( gg$Regulation, levels=c("Up","Down")  )
		# gg$Regulation = as.factor(gg$Regulation)
		 #gg$Regulation = relevel(gg$Regulation,"Up")
		 
		 p= ggplot(gg, aes(x=Comparisons, y=Genes, fill=  Regulation )  )+
			 geom_bar(position="dodge", stat="identity") + coord_flip() +
			 theme(legend.position = "top") + 
			 scale_fill_manual(values=c("red", "blue")) +
			 ylab("Number of differntially expressed genes") +
			 theme(axis.title.y=element_blank(),
				axis.text=element_text(size=14)) 
			
		print(p)
		})
    })

output$downloadSigGeneStats <- downloadHandler(
      filename = "DEG_stats.eps",
      content = function(file) {
	  cairo_ps(file, width = 8, height = 6, points = 8)
	  sigGeneStats4Download()
        dev.off()
      })	
output$sigGeneStatsTable <- renderTable({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		if(is.null(limma()$results) ) return(NULL)		
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
		tem=input$CountsDEGMethod
		tem = input$selectFactorsModel # responsive to changes in model and comparisons
		tem = input$selectModelComprions
		tem = input$submitModelButton 

		####################################
		
		isolate({ 
		
		results = limma()$results
        if(dim(results)[2] == 1) { # if only one comparison
           Up = sum( results == 1)
           Down = sum( results == -1 )
           stats = c( colnames(results), Up, Down)
           stats = t( as.data.frame( stats ) )
           row.names( stats ) = colnames(results)
           colnames( stats ) = c("Comparison","Up", "Down")

        } else {  #More than one comparisons

		   Up =  apply(results, 2, function(x) sum(x == 1) )
		   Down = apply(results, 2, function(x) sum(x == -1) ) 
		   stats = rbind(Up, Down)
		   stats = t(stats)
		   stats=cbind(rownames(stats), stats)
		   colnames(stats)[1]="Comparisons"
		   stats <- stats[dim(stats)[1]:1,] # reverse row order, to be the same with plot
        }
		 return(as.data.frame(stats))

		})
    }, digits = 0,spacing="s",include.rownames=F,striped=TRUE,bordered = TRUE, width = "auto",hover=T)
	
output$listComparisonsVenn <- renderUI({
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	tem = input$submitModelButton 
	tem = input$UpDownRegulated	
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectComparisonsVenn", label = NULL, # h6("Funtional Category"), 
                  choices = list("All" = "All"), selected = "All")  
		}	 else { 
				choices = setNames(limma()$comparisons, limma()$comparisons  )
				if(input$UpDownRegulated) {
				  tem = c( paste0("Up_", limma()$comparisons),paste0("Down_", limma()$comparisons) )
				  choices = setNames(tem, tem)
				
				}
				
				choices3 = choices;
				if(length(choices3)>3) choices3 = choices[1:3]  # by default only 3 are selected
				checkboxGroupInput("selectComparisonsVenn", 
									  h4("Select up to 5 comparisons"), 
									  choices = choices,
									  selected = choices3)	

	     } 
	})

	
output$listComparisons <- renderUI({
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC
	tem = input$submitModelButton ; tem = input$CountsDEGMethod	
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectContrast", label = NULL,  
                  choices = list("All" = "All"), selected = "All")  
				  }	 else { 
			selectInput("selectContrast", label="Select a comparison to examine. \"A-B\" means A vs. B (See heatmap). Interaction terms start with \"I:\"",choices=limma()$comparisons
	     )   } 
	})

	
output$listComparisonsPathway <- renderUI({
	tem = input$selectOrg
	tem = input$submitModelButton ; tem = input$CountsDEGMethod	
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectContrast1", label = NULL, # h6("Funtional Category"), 
                  choices = list("All" = "All"), selected = "All")  }	 else { 

	  selectInput("selectContrast1", label="Select a comparison to analyze:",choices=limma()$comparisons
	     )   } 
	})

	
output$listComparisonsGenome <- renderUI({
	tem = input$selectOrg
	tem = input$submitModelButton ; tem = input$CountsDEGMethod	
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectContrast1", label = NULL, # h6("Funtional Category"), 
                  choices = list("All" = "All"), selected = "All")  }	 else { 
	  selectInput("selectContrast2", label="Select a comparison to visualize:",choices=limma()$comparisons
	     )   } 
	})
	
	
DEG.data <- reactive({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC; 
		
			tem = input$CountsDEGMethod; 	
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
		tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
		tem= input$selectModelComprions;  tem= input$selectInteractions
		tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
		tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
		tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
    	####################################
		
		isolate({ 
			#genes = limma()$results
			genes = limma()$topGenes[[1]]
			if(length(limma()$topGenes ) >1 )
			  for( i in 2:length(limma()$topGenes) )
				genes <- cbind(genes, limma()$topGenes[[i]])
			genes$data = " " # add an empty column
			genes = merge(genes,convertedData(), by='row.names')

			colnames(genes) = gsub("\\.","-",colnames(genes))
			# add gene symbol
			ix = match( genes[,1], allGeneInfo()[,1])
			genes <- cbind(as.character( allGeneInfo()$symbol)[ix],genes) 
			colnames(genes)[1] = "Symbol"
			genes <- genes[,c(2,1,3:dim(genes)[2]) ]	
			
			return(genes)
		})
		})

output$download.DEG.data <- downloadHandler(
		filename = function() {"Diff_expression_all_comparisons.csv"},
		content = function(file) {
			write.csv(DEG.data(), file,row.names=FALSE)
	    }
	)	
	
AllGeneListsGMT <- reactive({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		tem = input$selectOrg
		tem=input$limmaPval; tem=input$limmaFC
		
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
		tem=input$CountsDEGMethod
		tem = input$selectFactorsModel # responsive to changes in model and comparisons
		tem = input$selectModelComprions
		tem = input$submitModelButton 

		####################################
		
		isolate({ 
			# upregulated genes
			resultsUp <- limma()$results
			resultsUp[which( resultsUp < 0 )] <- 0
			colnames(resultsUp) <- paste0( "UP_", colnames(resultsUp) )

			# down-regulated genes
			resultsDown <- limma()$results
			resultsDown[which( resultsDown > 0 )] <- 0
			colnames(resultsDown) <- paste0( "Down_", colnames(resultsDown) )

			results2 <- cbind(resultsUp, resultsDown)
            # reorder columns
            columnOrder = c()
            for( i in 1:dim(resultsUp)[2] )
               columnOrder = c(columnOrder, i, i + dim(resultsUp)[2] )
            results2 <- results2[,columnOrder]    

			geneList1 <- function (i) {
				ix = which(results2[,i] !=0 )
				return(paste0( gsub("^UP_|^Down_","", colnames(results2)[i]),"\t", 
                               gsub("_.*","", colnames(results2)[i]),", ", length(ix)," genes\t",
						paste(rownames(results2 )[ix], collapse="\t" ) ) )			
			}
			tem = sapply(1:dim(results2)[2],geneList1 )
				
			return( paste(tem, collapse="\n") )
		})
    })

output$downloadGeneListsGMT <- downloadHandler(
		filename = function() {"All_gene_lists_GMT.txt"},
		content = function(file) {
			write(AllGeneListsGMT(), file)
	    }
	)
################################################################
#   Differential gene expression  2
################################################################		
				
output$selectedHeatmap <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast;
	tem = input$heatColors1
	tem = input$CountsDEGMethod; 	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
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

selectedHeatmap4Download <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast;
	tem = input$heatColors1
	tem = input$CountsDEGMethod; 	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
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
	
    })

output$downloadSelectedHeatmap <- downloadHandler(
      filename = "Heatmap_comparison.eps",
      content = function(file) {
	  cairo_ps(file, width = 6, height = 8, points = 8)
	  selectedHeatmap4Download()
        dev.off()
      })	

selectedHeatmap.data <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	
	tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast;
	tem = input$heatColors1
	tem = input$CountsDEGMethod; 	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ;tem= input$NminSamples; tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	####################################	

	if( is.null(input$selectContrast)) return(NULL)
	if( is.null(limma()$results)) return(NULL)
	#if( is.null(selectedHeatmap.data() )) return(NULL) 
	isolate({ 
		  genes <- limma()$results
		  if( is.null(genes) ) return(NULL)
		  if(!grepl("I:", input$selectContrast) ) {  # if not interaction term
			ix = match(input$selectContrast, colnames(genes)) 
		  } else {
			#mismatch in comparison names for interaction terms for DESeq2
			#I:water_Wet.genetic_Hy 	 in the selected Contrast
			#Diff-water_Wet-genetic_Hy   in column names
			tem = gsub("I-","I:" ,colnames(genes))
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
		  

		 iz = findContrastSamples( input$selectContrast, 
									colnames(convertedData()),
									readSampleInfo(),
									input$selectFactorsModel,
									input$selectModelComprions, 
									factorReferenceLevels(),
									input$CountsDEGMethod,
									input$dataFileFormat
								)
	
		# color bar
		 bar = as.vector( genes[,ix]  ); # new R versions stopped autoconvert single column data frames to vectors.
		 names(bar) = row.names( genes[,ix] )
		 bar = bar[bar!=0]

		 # retreive related data		 
		 genes = convertedData()[iy,iz,drop=FALSE]
		 
		 genes = genes[order(bar),,drop=FALSE] # needs to be sorted because myheatmap2 does not reorder genes
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
	
	# Top DEGs  
output$geneList <- renderTable({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	
	tem = input$CountsDEGMethod; 	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
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
			tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts; tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
		tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
		tem= input$selectModelComprions;  tem= input$selectInteractions
		tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
		tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
		tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 

		noSig = as.data.frame("No significant genes find!")
		if( is.null(input$selectContrast) ) return(NULL)
		if( is.null( limma()$comparisons ) ) return(NULL) # if no significant genes found
		if( length(limma()$topGenes) == 0 ) return(noSig)
		if( is.null( geneListData() ) ) return(NULL)
		if( input$selectOrg == "NEW" | input$selectGO2 == "ID not recognized!" )
			tem <- merge(geneListData(), convertedData(), by.x = 'Top_Genes',by.y = 'row.names') else
		tem <- merge(geneListData(), convertedData(), by.x = 'Ensembl ID',by.y = 'row.names') 
		tem <- tem[order( -sign(tem[,2] ), -abs(tem[,2])),]
		tem$Regulation = "Up"
		tem$Regulation[which(tem[,2]<0 )] <- "Down"
		tem <- tem[,c(dim(tem)[2],1:( dim(tem)[2]-1) )  ]
		return( tem )
	
  })

# list of DEGs ranked by absolute values of lfc
# Ensembl ID	 log2 Fold Change	 Adj.Pval	 Symbol	 Chr	 Type	
geneListData <- reactive({
		if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
			tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
		tem = input$selectFactorsModel # responsive to changes in model and comparisons
		tem = input$selectModelComprions
		tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
		tem= input$selectModelComprions;  tem= input$selectInteractions
		tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
		tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
		tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
		####################################
		
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
		  if( input$selectGO2 == "ID not recognized!" 
		      | input$selectOrg == "NEW" 
		      | dim(allGeneInfo())[1] == 1)  # stringDB species
		    return (top1); 
		  
		  #convertedID = convertID(top1[,1],input$selectOrg, "GOBP" );#"gmax_eg_gene"
		  # tem <- geneInfo(convertedID,input$selectOrg) #input$selectOrg ) ;
		 #  tem <- geneInfo(converted(),input$selectOrg)
		  if(dim(allGeneInfo())[1] >1) { # Has geneInfo? 
		  top1 <- merge(top1, allGeneInfo(), by.x ="Top_Genes", by.y="ensembl_gene_id",all.x=T )
		  } else {  #StringDB species does not have gene info.
		   top1$symbol = top1$Top_Genes
		   top1$chr = NA		   
		   top1$gene_biotype = NA	
		   top1$band = NA
		   top1$chromosome_name = NA
		  }
		  
		  if ( sum( is.na(top1$band)) == dim(top1)[1] ) top1$chr = top1$chromosome_name else
			top1$chr = paste( top1$chromosome_name, top1$band, sep="")
		  
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
		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts; tem= input$NminSamples;tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 

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
  
volcanoPlot4Download <- reactive({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts; tem= input$NminSamples;tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 

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

  })

output$downloadVolcanoPlot <- downloadHandler(
      filename = "volcanoPlot.eps",
      content = function(file) {
	  cairo_ps(file, width = 8, height = 8)
	  volcanoPlot4Download()
        dev.off()
      })
  
output$volcanoPlotly <- renderPlotly({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts; tem= input$NminSamples;tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################
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
		tem = input$selectOrg; tem = input$noIDConversion
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################
	if( is.null(input$selectContrast) ) return(NULL)
	if( is.null( limma()$comparisons ) ) return(NULL) # if no significant genes found
	if( length(limma()$topGenes) == 0 ) return(NULL)
	if(grepl("I:", input$selectContrast) ) 	return(NULL) # if interaction term related comparison
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
	 
		 iz = findContrastSamples(	input$selectContrast, 
									colnames(convertedData()),
									readSampleInfo(),
									input$selectFactorsModel,
									input$selectModelComprions, 
									factorReferenceLevels(),
									input$CountsDEGMethod,
									input$dataFileFormat
								)
	 
	 genes <- convertedData()[,iz]
	 
	 g = detectGroups(colnames(genes), readSampleInfo())
	 
	 if(length(unique(g))  > 2) { plot.new(); text(0.5,0.5, "Not available.") } else{
		average1 <- apply( genes[, which( g == unique(g)[1] ) ],1,mean)

		average2 <- apply(  genes[, which( g == unique(g)[2] ) ],1,mean)

		genes2 <- cbind(average1,average2)
		rownames(genes2) = rownames(genes)
		genes2 <-  merge(genes2,top1,by="row.names")

		
		par(mar=c(5,5,1,1))
		plot(genes2$average2,genes2$average1,col = c("grey45","red","blue")[genes2$upOrDown],
		pch =16, cex = .3, xlab= paste("Average expression in", unique(g)[2] ), 
		ylab = paste("Average expression in", unique(g)[1] ),
		cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2	)    
		legend("bottomright",c("Upregulated","Downregulated"),fill = c("red","blue"),cex=1.3 ) 
	 }
	 
		})

	})
	 
  },height=450, width=500)

scatterPlot4Download <- reactive({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg; tem = input$noIDConversion
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################
	if( is.null(input$selectContrast) ) return(NULL)
	if( is.null( limma()$comparisons ) ) return(NULL) # if no significant genes found
	if( length(limma()$topGenes) == 0 ) return(NULL)
	if(grepl("I:", input$selectContrast) ) 	return(NULL) # if interaction term related comparison
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
	 
		 iz = findContrastSamples(	input$selectContrast, 
									colnames(convertedData()),
									readSampleInfo(),
									input$selectFactorsModel,
									input$selectModelComprions, 
									factorReferenceLevels(),
									input$CountsDEGMethod,
									input$dataFileFormat
								)
	 
	 genes <- convertedData()[,iz]
	 
	 g = detectGroups(colnames(genes))
	 
	 if(length(unique(g))  > 2) { plot.new(); text(0.5,0.5, "Not available.") } else{
		average1 <- apply( genes[, which( g == unique(g)[1] ) ],1,mean)

		average2 <- apply(  genes[, which( g == unique(g)[2] ) ],1,mean)

		genes2 <- cbind(average1,average2)
		rownames(genes2) = rownames(genes)
		genes2 <-  merge(genes2,top1,by="row.names")

		
		par(mar=c(5,5,1,1))
		plot(genes2$average2,genes2$average1,col = c("grey45","red","blue")[genes2$upOrDown],
		pch =16, cex = .3, xlab= paste("Average expression in", unique(g)[2] ), 
		ylab = paste("Average expression in", unique(g)[1] ),
		cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2	)    
		legend("bottomright",c("Upregulated","Downregulated"),fill = c("red","blue"),cex=1.3 )


	 
	 }
	 
		})

	})
	 
  })

output$downloadScatterPlot <- downloadHandler(
      filename = "scatterPlot.eps",
      content = function(file) {
	  cairo_ps(file, width = 8, height = 8)
	  scatterPlot4Download()
        dev.off()
      })
	  
output$scatterPlotly <- renderPlotly({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################
	if( is.null(input$selectContrast) ) return(NULL)
	if( is.null( limma()$comparisons ) ) return(NULL) # if no significant genes found
	if( length(limma()$topGenes) == 0 ) return(NULL)
	if(grepl("I:", input$selectContrast) ) 	return(NULL) # if interaction term related comparison
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


	 iz = findContrastSamples(	input$selectContrast, 
									colnames(convertedData()),
									readSampleInfo(),
									input$selectFactorsModel,
									input$selectModelComprions, 
									factorReferenceLevels(),
									input$CountsDEGMethod,
									input$dataFileFormat
								)
	 
	 genes <- convertedData()[,iz]
	 
	 g = detectGroups(colnames(genes))
	 
	 if(length(unique(g))  > 2) { plot.new(); text(0.5,0.5, "Not available.") } else{
		average1 <- apply( genes[, which( g == unique(g)[1] ) ],1,mean)

		average2 <- apply(  genes[, which( g == unique(g)[2] ) ],1,mean)

		genes2 <- cbind(average1,average2)
		rownames(genes2) = rownames(genes)
	
		genes <-  merge(genes2,top1,by="row.names")
	 
	 
	 
	 
 
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
	layout(xaxis = list(size =35,title = paste("Average expression in", unique(g)[2]  )), 
           yaxis = list(size =35, title = paste("Average expression in", unique(g)[1]  )),
		   showlegend = FALSE)
		   
	ggplotly(p)
	}
		})

	})
	 
  })

  
output$MAplot <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg; tem = input$noIDConversion
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################
	if( is.null(input$selectContrast) ) return(NULL)
	if( is.null( limma()$comparisons ) ) return(NULL) # if no significant genes found
	if( length(limma()$topGenes) == 0 ) return(NULL)
	if(grepl("I:", input$selectContrast) ) 	return(NULL) # if interaction term related comparison
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
	 
		 iz = findContrastSamples(	input$selectContrast, 
									colnames(convertedData()),
									readSampleInfo(),
									input$selectFactorsModel,
									input$selectModelComprions, 
									factorReferenceLevels(),
									input$CountsDEGMethod,
									input$dataFileFormat
								)
	 


		average1 <- as.data.frame( apply( convertedData()[,iz],1,mean) )
		colnames(average1) = "Average"
		rownames(average1) = rownames(convertedData())
		
		genes2 <-  merge(average1,top1,by="row.names")

		par(mar=c(5,5,1,1))
		plot(genes2$Average,genes2$Fold,col = c("grey45","red","blue")[genes2$upOrDown],
		pch =16, cex = .3, xlab= "Average expression", 
		ylab = "Log2 fold change",
		cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2	)   
			abline(h=0)
		legend("bottomright",c("Upregulated","Downregulated"),fill = c("red","blue"),cex=1.3 )


	 

	 
		})

	})
	 
  },height=450, width=500)

MAplot4Download <- reactive({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg; tem = input$noIDConversion
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################
	if( is.null(input$selectContrast) ) return(NULL)
	if( is.null( limma()$comparisons ) ) return(NULL) # if no significant genes found
	if( length(limma()$topGenes) == 0 ) return(NULL)
	if(grepl("I:", input$selectContrast) ) 	return(NULL) # if interaction term related comparison
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
	 
		 iz = findContrastSamples(	input$selectContrast, 
									colnames(convertedData()),
									readSampleInfo(),
									input$selectFactorsModel,
									input$selectModelComprions, 
									factorReferenceLevels(),
									input$CountsDEGMethod,
									input$dataFileFormat
								)
	 


		average1 <- as.data.frame( apply( convertedData()[,iz],1,mean) )
		colnames(average1) = "Average"
		rownames(average1) = rownames(convertedData())
		
		genes2 <-  merge(average1,top1,by="row.names")

		par(mar=c(5,5,1,1))
		plot(genes2$Average,genes2$Fold,col = c("grey45","red","blue")[genes2$upOrDown],
		pch =16, cex = .3, xlab= "Average expression", 
		ylab = "Log2 fold change",
		cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2	)   
			abline(h=0)
		legend("bottomright",c("Upregulated","Downregulated"),fill = c("red","blue"),cex=1.3 )
		})

	})
	 
  })

output$downloadMAPlot <- downloadHandler(
      filename = "MA_Plot.eps",
      content = function(file) {
	  cairo_ps(file, width = 8, height = 8)
	  MAplot4Download()
        dev.off()
      })

output$MAplotly <- renderPlotly({
    if (is.null(input$file1)&& input$goButton == 0  )   return(NULL)
		tem = input$selectOrg; tem = input$noIDConversion
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################
	if( is.null(input$selectContrast) ) return(NULL)
	if( is.null( limma()$comparisons ) ) return(NULL) # if no significant genes found
	if( length(limma()$topGenes) == 0 ) return(NULL)
	if(grepl("I:", input$selectContrast) ) 	return(NULL) # if interaction term related comparison
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
	 
		 iz = findContrastSamples(	input$selectContrast, 
									colnames(convertedData()),
									readSampleInfo(),
									input$selectFactorsModel,
									input$selectModelComprions, 
									factorReferenceLevels(),
									input$CountsDEGMethod,
									input$dataFileFormat
								)
	 


		average1 <- as.data.frame( apply( convertedData()[,iz],1,mean) )
		colnames(average1) = "Average"
		rownames(average1) = rownames(convertedData())
		
		genes <-  merge(average1,top1,by="row.names")

		
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
		
		p <- plot_ly(data= genes, x = ~Average, y= ~Fold, color = ~type, text = ~Symbols, type = 'scatter')%>% 
		layout(xaxis = list(size =35,title = "Average expression"), 
			   yaxis = list(size =35, title = "Log2 fold change"),
			   showlegend = FALSE)
			   
		ggplotly(p)
			
			})

	})
	 
  })
  
geneListGOTable <- reactive({		
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if( is.null(input$selectContrast)) return(NULL)
		if( is.null( input$selectGO2) ) return (NULL)
		if( input$selectGO2 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.")) #No matching species

		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO2
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
		tem = input$removeRedudantSets
        tem = input$UseFilteredGenesEnrich
		####################################
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
				if(input$removeRedudantSets) reduced = redudantGeneSetsRatio else reduced = FALSE
                if(input$UseFilteredGenesEnrich) {
                   convertedDataBackground <- converted() 
                } else {
                  convertedDataBackground <- NULL
                 }

				result = FindOverlap (convertedID,allGeneInfo(), input$selectGO2,input$selectOrg,1, reduced, convertedDataBackground) }

			if( dim(result)[2] ==1) next;   # result could be NULL
			if(i == -1) result$direction = "Up regulated"  else result$direction = "Down regulated"
			if (pp==0 ) { results1 <- result; pp = 1;} else  results1 = rbind(results1,result)
		}

		if ( pp == 0 ) return (NoSig)
		if ( is.null( results1) ) return (NoSig)
		if( dim(results1)[2] == 1 ) return(NoSig)  # Returns a data frame: "No significant results found!"
		
		results1= results1[,c(6,1,2,4,5)]
		colnames(results1)= c("List","FDR","Genes","GO terms or pathways","Genes")
		minFDR = 0.01
		if(min(results1$FDR) > minFDR ) results1 = as.data.frame("No signficant enrichment found.") else
		results1 = results1[which(results1$FDR < minFDR),]
		
		incProgress(1, detail = paste("Done")) 
		
		if(dim(results1)[2] != 5) return(NoSig)
		colnames(results1)= c("Direction","adj.Pval","nGenes","Pathways","Genes")
		rownames(results1)=1:nrow(results1)

		
		return( results1 )
		 })#progress
		}) #isolate
})		


output$geneListGO <- renderTable({	
  if(is.null(geneListGOTable())) return(NULL)
  results1 = geneListGOTable()
  tem = input$removeRedudantSets
  if(dim(results1)[2] ==1) return(results1) else { 
	results1$adj.Pval <- sprintf("%-2.1e",as.numeric(results1$adj.Pval) )
	results1[,1] <- as.character(results1[,1])
	results1[ duplicated (results1[,1] ),1 ] <- ""  
	
	return( results1[,-5] )
	}	
  }, digits = 0,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T, 
  sanitize.text.function = function(x) x ) # for hypertext rendering
	#   output$selectedHeatmap <- renderPlot({       hist(rnorm(100))    })	

output$downloadGOTerms <- downloadHandler(
		filename = function() {"Enriched.csv"},
		content = function(file) {
			write.csv( removeHypertext( geneListGOTable() ), file)
	    }
	)

output$enrichmentPlotDEG2 <- renderPlot({
    if(is.null(geneListGOTable())) return(NULL)
	tem = input$removeRedudantSets
    
	enrichmentPlot(removeHypertext( geneListGOTable() ), 45  )

}, height=600, width=800)

output$enrichmentPlotDEG24Download <- downloadHandler(
      filename = "enrichmentPlotDEG2.eps",
      content = function(file) {
	  cairo_ps(file, width = 10, height = 6)
	  enrichmentPlot(removeHypertext( geneListGOTable() ),41  )
        dev.off()
      })
	  
output$enrichmentNetworkPlot <- renderPlot({
    if(is.null(geneListGOTable())) return(NULL)
	tem = input$removeRedudantSets

	enrichmentNetwork_old_remove_later(removeHypertext( geneListGOTable() ),
                       layout_change = input$layoutButton2 )

}, height=900, width=900)	  

# define a network
networkDEG <- reactive({
    if(is.null(geneListGOTable())) return(NULL)
    if(is.null( input$wrapTextNetworkDEG )) return(NULL)
	tem = input$removeRedudantSets

    network <- removeHypertext( geneListGOTable() )
    
    if(is.null( input$upORdownRegDEG )) return(NULL)
    if(input$upORdownRegDEG != "Both")
       network <- network[ grepl(input$upORdownRegDEG, network$Direction), ]
    if(dim(network)[1] == 0) return(NULL)

    if(input$wrapTextNetworkDEG)
      network$Pathways <- wrap_strings( network$Pathways ) # wrap long pathway names using default width of 30 10/21/19

    g <- enrichmentNetwork(network,layoutButton = input$layoutVisDEG, edge.cutoff = input$edgeCutoffDEG )

    data1 <- toVisNetworkData(g)
    
    # Color codes: https://www.rapidtables.com/web/color/RGB_Color.html
    data1$nodes$shape <- "dot"
    # remove the color change of nodes
    #data1$nodes <- subset(data1$nodes, select = -color)
    
    data1$nodes$size <- 5 + data1$nodes$size^2 
    
    return(data1)
})

  # note the same code is used twice as above. They need to be updated together!!!	  
output$visNetworkDEG <- renderVisNetwork({
    if(is.null(geneListGOTable())) return(NULL)
    if(dim(geneListGOTable())[1] == 1) return(NULL)
    if(is.null( input$wrapTextNetworkDEG )) return(NULL)
	tem = input$removeRedudantSets
    if(is.null(networkDEG() )) return(NULL)

    data1 <- networkDEG()
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
  
output$visNetworkDEGDownload <- downloadHandler(
    filename = "enrichmentPlotNetwork_DEG.html",
    content = function(file) {

    # same code as above to generate the network
	tem = input$removeRedudantSets

    data1 <- networkDEG()
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
      )  %>% 
        visSave(file = file, background = "white")
      
    })  


  output$downloadNodesDEG <- downloadHandler(
    filename = function() {"network_nodes.csv"},
    content = function(file) {      
      write.csv(networkDEG()$nodes, file, row.names=FALSE)
    }
  )
  output$downloadEdgesDEG <- downloadHandler(
    filename = function() {"network_edges.csv"},
    content = function(file) {    
      write.csv(networkDEG()$edges, file, row.names=FALSE)
    }
  )  
	  
output$DEG.Promoter <- renderTable({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	if( is.null(input$selectContrast)) return(NULL)

	tem = input$selectOrg; tem = input$radio.promoter; tem = input$noIDConversion; tem=input$missingValue
	tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO2
	tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
	tem = input$minCounts; tem= input$NminSamples;tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################
	
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

#----------------------------------------------------
# STRING-db functionality
# find Taxonomy ID from species official name 
findTaxonomyID <- reactive({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		#if (input$submit2STRINGdb == 0)   return(NULL)
		if( is.null(input$selectContrast)) return(NULL)
		if( is.null( input$selectGO2) ) return (NULL)
		 
		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO2
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
		tem = input$removeRedudantSets
		
    if(!is.null(input$speciesName) ) { # if species name is entered
	   ix = match(input$speciesName, STRING10_species$official_name)
	   } else if( input$selectGO2 != "ID not recognized!" )
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
	library(STRINGdb,verbose=FALSE)
						   
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		#if (input$submit2STRINGdb == 0)   return(NULL)
		if( is.null(input$selectContrast)) return(NULL)

		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO2
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
		tem = input$removeRedudantSets


		####################################
		if( is.null(limma()$results) ) return(NULL)
		if( is.null(geneListData()) ) return(NULL) # this has to be outside of isolate() !!!
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
		genes <- geneListData()
		colnames(genes)[1:2]=c("gene","lfc")
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
    # tem=table(STRING10_species$kingdom)
    # tem=paste(tem, names(tem), sep=" ", collapse=", ")
	# tem =paste0("Total species in STRING:",tem)
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
						   
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		#if (input$submit2STRINGdb == 0)   return(NULL)
		if( is.null(input$selectContrast)) return(NULL)
		if( is.null( input$selectGO2) ) return (NULL)

		if( is.null(STRINGdb_geneList() ) ) return("No genes mapped by STRINGdb. Please enter or double-check species name above.")
		if(! is.null(STRINGdb_geneList() ) ) { 
			tem=paste0( 100*round(STRINGdb_geneList()$ratio,3), "% genes mapped by STRING web server.")
			if(STRINGdb_geneList()$ratio <0.3 ) tem = paste(tem, "Warning!!! Very few gene mapped. Double check if the correct species is selected.")
			return( tem  )
		}
}) 

stringDB_GO_enrichmentData <- function(input, output) {
	library(STRINGdb,verbose=FALSE)
						   
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		#if (input$submit2STRINGdb == 0)   return(NULL)
		if( is.null(input$selectContrast)) return(NULL)
		if( is.null( input$selectGO2) ) return (NULL)
		#if( input$selectGO2 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.")) #No matching species

		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO2
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
		tem = input$removeRedudantSets

		tem = input$STRINGdbGO
		taxonomyID = findTaxonomyID(  )
		if(is.null( taxonomyID ) ) return(NULL)		
		####################################
		if( is.null(limma()$results) ) return(NULL)
		if( is.null(geneListData()) ) return(NULL) # this has to be outside of isolate() !!!
		#if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
		if(is.null(STRINGdb_geneList() ) ) return(NULL)
		
		stringDB_GO_enrichmentDataR <- reactive({
		  withProgress(message=sample(quotes,1), detail ="Enrichment analysis", {
		    #Intialization
		    string_db <- STRINGdb$new( version=STRING_DB_VERSION, species=taxonomyID,
		                               score_threshold=0, input_directory="" )
		    
		    # using expression data
		    genes <- selectedHeatmap.data()$genes
		    fc = selectedHeatmap.data()$bar 
		    resultFilter <- NULL; result <- NULL
		    minGenesEnrichment <- 1
		    if(is.null(genes)) {
		      return(NULL) 
		    } else if (dim(genes)[1] <= minGenesEnrichment ) {
		      return(NULL) # if has only few genes
		    } else {
		      pp <- 0
		      for( i in c(1:2) ) {
		        if( length(which(fc*i<0)) <= minGenesEnrichment) {
		          next
		        }
		        query = rownames(genes)[which(fc*i<0)]
		        if( length(query) <= minGenesEnrichment) {
		          next 
		        }
		        
		        ids = STRINGdb_geneList()[[i]]
		        
		        if( length(ids) <= minGenesEnrichment) {
		          next
		        } 
		        incProgress(1/3  )
		        result <- string_db$get_enrichment( ids, category = input$STRINGdbGO, methodMT = "fdr", iea = TRUE )
		        if(nrow(result) == 0 || is.null(result)) {
		          next
		        } 
		        
		        if (i == 1) {
		          result$direction = "Up regulated"
		        } else {
		          result$direction = "Down regulated"
		        }
		        
		        if (pp == 0) {
		          resultFilter <- result
		          pp = 1
		        } else {
		          resultFilter <- rbind(resultFilter, result)
		        }
		      } #end of for
		      
		      if(nrow(resultFilter) == 0 || is.null(resultFilter) || pp == 0) {
		        return(NULL)
		      } else {
		        incProgress(1/3)
		        if(min(resultFilter$fdr) > input$STRINGFDR) {
		          return (NULL)
		        } else {
		          resultFilter <- resultFilter[which(resultFilter$fdr < input$STRINGFDR),]
		          incProgress(1, detail = paste("Done")) 
		          return(resultFilter)
		        } #end of check minFDR
		      }# check results 
		    } # end of check genes if 
		  })#progress
		}) #reactive
		
		result <- stringDB_GO_enrichmentDataR()
		if (is.null(result)) {
		  result <- as.data.frame("No significant enrichment found.")
		} else {
		  resultDownload <- result
		  result <- dplyr::select(result,
		                          c('fdr','number_of_genes','term',
		                            'description'))
		  colnames(result) <- c('FDR','nGenes','GO terms or pathways',
		                        'Description')
		  if(nrow(result) > 30) {
		    result <- result[1:30,] 
		  }
		} #end of if else 
		output$stringDB_GO_enrichment <- renderTable(result,
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
		    write.csv(resultDownload, file)
		  }
		) #downloadHandler
		
} #end of stringDB_GO_enrichmentData

observeEvent(input$submit2STRINGdb, {
    stringDB_GO_enrichmentData(input = input, 
                               output = output) 
}) # end of observeEvent
   
output$stringDB_network1 <- renderPlot({
	library(STRINGdb)
						   
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		#if (input$submit2STRINGdb == 0)   return(NULL)
		if( is.null(input$selectContrast)) return(NULL)
		if( is.null( input$selectGO2) ) return (NULL)
		#if( input$selectGO2 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.")) #No matching species

		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO2
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
		tem = input$removeRedudantSets

		tem = input$STRINGdbGO
		tem = input$nGenesPPI
		taxonomyID = findTaxonomyID( )
		if(is.null( taxonomyID ) ) return(NULL)		
		####################################
		if( is.null(limma()$results) ) return(NULL)
		if( is.null(geneListData()) ) return(NULL) # this has to be outside of isolate() !!!
		# if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
		NoSig = as.data.frame("No significant enrichment found.")
		if(is.null(STRINGdb_geneList() ) ) return(NULL)
		
		isolate({
		withProgress(message=sample(quotes,1), detail ="Enrichment analysis", {
		#Intialization
		string_db <- STRINGdb$new( version=STRING_DB_VERSION, species=taxonomyID,
							   score_threshold=0, input_directory="" )
		# only up regulated is ploted		
		for( i in c(1:1) ) {
			incProgress(1/2,detail = paste("Plotting network")  )
			
		
			ids = STRINGdb_geneList()[[i]]
			if(length(ids)> input$nGenesPPI )  # n of genes cannot be more than 400
				ids <- ids[1:input$nGenesPPI]
			incProgress(1/3  )
			string_db$plot_network( ids,add_link=FALSE)

		}

		 })#progress
		}) #isolate						   
}, width = 1000, height=600)

output$stringDB_network_link <- renderUI({
		library(STRINGdb,verbose=FALSE)
						   
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		#if (input$submit2STRINGdb == 0)   return(NULL)
		if( is.null(input$selectContrast)) return(NULL)
		if( is.null( input$selectGO2) ) return (NULL)
		#if( input$selectGO2 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.")) #No matching species

		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO2
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart
		tem = input$removeRedudantSets

		tem = input$STRINGdbGO
		tem = input$nGenesPPI
		taxonomyID = findTaxonomyID( )
		if(is.null( taxonomyID ) ) return(NULL)		
		
		####################################
		if( is.null(limma()$results) ) return(NULL)
		if( is.null(geneListData()) ) return(NULL) # this has to be outside of isolate() !!!
		#if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
		NoSig = as.data.frame("No significant enrichment found.")
		if(is.null(STRINGdb_geneList() ) ) return(NULL)
		
		isolate({
		withProgress(message=sample(quotes,1), detail ="PPI Enrichment and link", {
		#Intialization
		string_db <- STRINGdb$new( version=STRING_DB_VERSION, species=taxonomyID,
							   score_threshold=0, input_directory="" )
			# upregulated
		   ids = STRINGdb_geneList()[[1]]
			if(length(ids)> input$nGenesPPI )  # n of genes cannot be more than 400
				ids <- ids[1:input$nGenesPPI]
			incProgress(1/4  )
			link1 = string_db$get_link( ids)
			Pval1 = string_db$get_ppi_enrichment( ids)
			tem = "<h5> Interactive and annotated PPI networks among DEGs: <br/>  "
			tem = paste(tem, "<a href=\"", link1, "\" target=\"_blank\"> Up-regulated; </a>"  )

			# downregulated
			ids = STRINGdb_geneList()[[2]]
			if(length(ids)> input$nGenesPPI )
				ids <- ids[1:input$nGenesPPI]
			incProgress(2/4  )
			link2 = string_db$get_link( ids)
			Pval2 = string_db$get_ppi_enrichment( ids)				
			tem = paste(tem, " &nbsp  <a href=\"", link2, "\"target=\"_blank\"> Down-regulated; </a>"  )
			
			# both up and down with color code
			incProgress(3/4  )
			geneTable = STRINGdb_geneList()$geneTable
			if(nrow(geneTable)> input$nGenesPPI ) 
				geneTable <- geneTable[1:input$nGenesPPI,]
			geneTable =  string_db$add_diff_exp_color( geneTable, logFcColStr="lfc" ) 
			payload_id <- string_db$post_payload( geneTable$STRING_id,colors=geneTable$color )
			link3 = string_db$get_link(geneTable$STRING_id, payload_id = payload_id)			
			tem = paste(tem, " &nbsp  <a href=\"", link3, "\"target=\"_blank\"> Both with fold-changes color coded.</a></h5>"  )

            tem2 = paste("<h5> PPI enrichment P values: ")  
			tem2 = paste0(tem2,"Up-regulated: ", sprintf("%-3.2e",Pval1[1]), " &nbsp Down-regulated: ", sprintf("%-3.2e",Pval2[1]),".")
			tem2 = paste(tem2, " Small P values indicate more PPIs among DEGs than background. </h5>" )
			tem = paste(tem2,tem )
			return(HTML(tem))	
		
			incProgress(1  )

		 })#progress
		}) #isolate			

}) 

################################################################
#   Pathway analysis
################################################################
 
# this updates geneset categories based on species and file
output$selectGO1 <- renderUI({   # gene set for pathway analysis
	  tem = input$selectOrg;
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectGO", label = h5("Select genesets (Choose KEGG to show pathway diagrams):"), 
                  choices = list("All available gene sets" = "All", "GO Biological Process" = "GOBP","GO Molecular Function" = "GOMF","GO Cellular Component" = "GOCC",
                                "KEGG metabolic pathways" = "KEGG"), selected = "GOBP")  }	 else { 
								
	  selectInput("selectGO", label=h5("Select genesets (Choose KEGG to show pathway diagrams):"),
		choices=gmtCategory(converted(), convertedData(), input$selectOrg,input$gmtFile)
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
	
output$selectGO5 <- renderUI({
	  tem = input$selectOrg
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectGO5", label = NULL, # h6("Funtional Category"), 
                  choices = list("All available gene sets" = "All", "GO Biological Process" = "GOBP","GO Molecular Function" = "GOMF","GO Cellular Component" = "GOCC",
                                "KEGG metabolic pathways" = "KEGG"), selected = "GOBP")  }	 else { 
								
	  selectInput("selectGO5", label=NULL,choices=gmtCategory(converted(), convertedData(), input$selectOrg,input$gmtFile)
	     ,selected = "GOBP" )   } 
	})

output$selectGO6 <- renderUI({  # for GSEA using PCA loadings
	  tem = input$selectOrg
      if (is.null(input$file1)&& input$goButton == 0 )
       { selectInput("selectGO6", label = h5("Select genesets"), # h6("Funtional Category"), 
                  choices = list("All available gene sets" = "All", "GO Biological Process" = "GOBP","GO Molecular Function" = "GOMF","GO Cellular Component" = "GOCC",
                                "KEGG metabolic pathways" = "KEGG"), selected = "GOBP")  }	 else { 
								
	  selectInput("selectGO6", label=h5("Select genesets"),choices=gmtCategory(converted(), convertedData(), input$selectOrg,input$gmtFile)
	     ,selected = "GOBP" )   } 
	})

output$PGSEAplot <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	library(PGSEA,verbose=FALSE)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO;		tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################	
	if(is.null(input$selectGO ) ) return (NULL)
	if(input$selectGO == "ID not recognized!" ) return( NULL)
	isolate({ 
	withProgress(message=sample(quotes,1), detail ="Running pathway analysis", {
	myrange = c(input$minSetSize, input$maxSetSize)
	genes = convertedData()
if (is.null(input$selectContrast1 ) ) return(NULL)
	incProgress(1/4,"Retrieving gene sets")
	gmt = GeneSets()
	incProgress(2/4,"Runing PGSEA.")

	# find related samples
	iz = findContrastSamples(input$selectContrast1, colnames(convertedData()),readSampleInfo(),
										input$selectFactorsModel,input$selectModelComprions, 
										factorReferenceLevels(),input$CountsDEGMethod,
										input$dataFileFormat  )

	
	
	genes = genes[,iz]	
	
	subtype = detectGroups(colnames(genes )) 
    if(length( GeneSets() )  == 0)  { plot.new(); text(0.5,0.5, "No gene sets!")} else {
	result = PGSEApathway(converted(),genes, input$selectOrg,input$selectGO,
	             GeneSets(),  myrange, input$pathwayPvalCutoff, input$nPathwayShow 	)
					 
	if( is.null(result$pg3) ) { plot.new(); text(0.5,1, "No significant pathway found!")} else 
	smcPlot(result$pg3,factor(subtype),scale = c(-max(result$pg3), max(result$pg3)), 
	show.grid = T, margins = c(3,1, 13, 38), col = .rwb,cex.lab=0.5)
    }
	
	})
	})
    }, height = 800, width = 800)

PGSEAplot4Download <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	library(PGSEA,verbose=FALSE)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO;		tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
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

	# find related samples
	iz = findContrastSamples(input$selectContrast1, colnames(convertedData()),readSampleInfo(),
										input$selectFactorsModel,input$selectModelComprions, 
										factorReferenceLevels(),input$CountsDEGMethod,
										input$dataFileFormat  )

	
	
	genes = genes[,iz]	
	
	subtype = detectGroups(colnames(genes )) 
    if(length( GeneSets() )  == 0)  { plot.new(); text(0.5,0.5, "No gene sets!")} else {
	result = PGSEApathway(converted(),genes, input$selectOrg,input$selectGO,
	             GeneSets(),  myrange, input$pathwayPvalCutoff, input$nPathwayShow 	)
					 
	if( is.null(result$pg3) ) { plot.new(); text(0.5,1, "No significant pathway found!")} else 
	smcPlot(result$pg3,factor(subtype),scale = c(-max(result$pg3), max(result$pg3)), 
	show.grid = T, margins = c(3,1, 13, 38), col = .rwb,cex.lab=0.5)
    }
	
	})
	})
    })

output$PGSEAplot.Download <- downloadHandler(
      filename = function() {paste0("PGSEA",input$selectContrast1,"(",input$selectGO,")",".eps")},
      content = function(file) {
	  cairo_ps(file, width = 10, height = 8)
	  PGSEAplot4Download()
        dev.off()
      })
	  
PGSEAplot.data <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO;		tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 

	isolate({ 
	withProgress(message="Running pathway analysis", {
	myrange = c(input$minSetSize, input$maxSetSize)
	genes = convertedData()
if (is.null(input$selectContrast1 ) ) return(NULL)
	gmt = GeneSets()
	incProgress(2/8)

	if(0) { 
	iz= match( detectGroups(colnames(genes)), unlist(strsplit( input$selectContrast1, "-"))	  )
    iz = which(!is.na(iz))
	if (grepl("I:",input$selectContrast1) == 1) iz=1:(dim(genes)[2]) 
	if (length(iz) == 0) iz=1:(dim(genes)[2]) 
	}
	
			  # find sample related to the comparison
		 iz= match( detectGroups(colnames(convertedData())), unlist(strsplit( input$selectContrast1, "-"))	  )
		 iz = which(!is.na(iz))		 
		 if ( !is.null(readSampleInfo()) & !is.null(input$selectFactorsModel) & length(input$selectModelComprions)>0 ) {
			comparisons = gsub(".*: ","",input$selectModelComprions)   # strings like: "groups: mutant vs. control"
			comparisons = gsub(" vs\\. ","-",comparisons)		
			factorsVector= gsub(":.*","",input$selectModelComprions) # corresponding factors
			ik = match( input$selectContrast1, comparisons )   # selected contrast lookes like: "mutant-control"
			if (is.na(ik)) iz=1:(dim(convertedData())[2])  else {  # interaction term, use all samples		
				selectedfactor= factorsVector[ ik ] # corresponding factors
				iz= match( readSampleInfo()[,selectedfactor], unlist(strsplit( input$selectContrast1, "-"))	  )
				iz = which(!is.na(iz))				
			}
		 }

		 if (grepl("I:",input$selectContrast1)) iz=1:(dim(convertedData())[2]) # if it is factor design use all samples
		 if( is.na(iz)[1] | length(iz)<=1 )    iz=1:(dim(convertedData())[2]) 
		#cat("\n IZ:",iz)
	
	genes = genes[,iz]

	subtype = detectGroups(colnames(genes )) 
    if(length( GeneSets() )  == 0)  { return(as.data.frame("No significant pathway!"))} else {
	result = PGSEApathway(converted(),genes, input$selectOrg,input$selectGO,
	             GeneSets(),  myrange, input$pathwayPvalCutoff, input$nPathwayShow 	)
					 
	if( is.null(result$pg3) ) { return(as.data.frame("No significant pathway!"))} else 
	   return( as.data.frame(result$pg3) )
    }
	
	})
	})
    })
	
output$PGSEAplotAllSamples <- renderPlot({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	library(PGSEA,verbose=FALSE)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO ; tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################
	
	if(is.null(input$selectGO ) ) return (NULL)
	if(input$selectGO == "ID not recognized!" ) return( NULL)
	isolate({ 
	withProgress(message=sample(quotes,1), detail ="Running pathway analysis", {
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
	show.grid = T, margins = c(3,1, 13, 38), col = .rwb,cex.lab=0.5)
    }
	
	})
	})
    }, height = 800, width = 800)

PGSEAplotAllSamples4download <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	library(PGSEA,verbose=FALSE)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO ;  tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
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
	show.grid = T, margins = c(3,1, 13, 38), col = .rwb,cex.lab=0.5)
    }
	
	})
	})
    })

output$PGSEAplotAllSamples.Download <- downloadHandler(
      filename = function() {paste0("PGSEA_all_samples_",input$selectContrast1,"(",input$selectGO,")",".eps")},
      content = function(file) {
	  cairo_ps(file, width = 10, height = 8)
	  PGSEAplotAllSamples4download()
        dev.off()
      })
	  
PGSEAplotAllSamples.data <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	library(PGSEA,verbose=FALSE)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO  ; tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
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
					 
	if( is.null(result$pg3) ) { return(as.data.frame("No significant pathway!"))} else 
	result = PGSEApathway(converted(),genes, input$selectOrg,input$selectGO,
	             GeneSets(),  myrange, input$pathwayPvalCutoff, input$nPathwayShow 	)
					 
	if( is.null(result$pg3) ) { return(as.data.frame("No significant pathway!"))} else 
	   return( as.data.frame(result$pg3) )
    }
	
	})
	})
    })

output$gagePathway <- renderTable({

    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO;	tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter ; tem =input$NminSamples2}
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################
			
	isolate({ 
	gagePathwayData()
	})
  },digits=0,align="l",include.rownames=FALSE,striped=TRUE,bordered = TRUE, width = "auto",hover=T)
 
 
gagePathwayData <- reactive({
	library(gage,verbose=FALSE) # pathway analysis	
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)

	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO
	tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold; tem =input$pathwayMethod

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	tem=input$GenePvalCutoff
	####################################
	
	if(is.null(input$selectGO ) ) return (NULL)
	if(input$selectGO == "ID not recognized!" ) return( as.data.frame("Gene ID not recognized." ))
	isolate({ 
		withProgress(message=sample(quotes,1), detail ="Running pathway analysis using GAGE", {
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

		  top1 = top1[which(top1$FDR <input$GenePvalCutoff) , ]
 
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

		  return( top1)
		}) # progress
	}) #isloate
  })

  
output$fgseaPathway <- renderTable({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO; tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################
	
	isolate({ 
	fgseaPathwayData()
	})
  },digits=0,align="l",include.rownames=FALSE,striped=TRUE,bordered = TRUE, width = "auto",hover=T)
 
 
fgseaPathwayData <- reactive({
	library(fgsea,verbose=FALSE) # fast GSEA
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO; tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold; tem =input$pathwayMethod
	if(is.null(input$selectGO ) ) return (NULL)
	if(input$selectGO == "ID not recognized!" ) return( as.data.frame("Gene ID not recognized." ))

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	tem=input$GenePvalCutoff
	####################################
	
	isolate({ 
	withProgress(message=sample(quotes,1), detail ="Running pathway analysis using fgsea", {
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
	  
	  # remove some genes
	  top1 = top1[which(top1$FDR <input$GenePvalCutoff) , ]
	  
	  incProgress(1/4,"Retrieving gene sets")
	  gmt = GeneSets() 
     if(length( GeneSets() )  == 0)  { return(as.data.frame("No gene set found!"))}

	  #converted = convertID(rownames(top1),input$selectOrg)
	  #
	   #gmt = readGeneSets(converted, top1, input$selectGO, input$selectOrg, myrange )
     # cat("Sets",length(gmt))
	 incProgress(2/4,message=sample(quotes,1), detail = "Runing GSEA using the fgsea package.")
	 fold = top1[,1]; names(fold) <- rownames(top1)
	 if(input$absoluteFold) fold <- abs(fold) # use absolute value of fold change, disregard direction
	 paths <- fgsea(pathways = gmt, 
                  stats = fold,
                  minSize=input$minSetSize,
                  maxSize=input$maxSetSize,
				  nproc = 6, # cpu cores
                  nperm=100000)
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
	library(ReactomePA,verbose=FALSE) # pathway analysis
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO; tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	; tem =input$pathwayMethod
	if(is.null(input$selectGO ) ) return (NULL)
	if(input$selectGO == "ID not recognized!" ) return( as.data.frame("Gene ID not recognized." ))

	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	tem=input$GenePvalCutoff
	####################################
	
	isolate({ 
	withProgress(message=sample(quotes,1), detail ="Running pathway analysis using ReactomePA", {
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

	  # remove some genes
	  top1 = top1[which(top1$FDR <input$GenePvalCutoff) , ]
	  
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
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################
	
	isolate({ 
	ReactomePAPathwayData()
	})
  },digits=0,align="l",include.rownames=FALSE,striped=TRUE,bordered = TRUE, width = "auto",hover=T)

  
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
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################
			

      if (is.null(input$file1)&& input$goButton == 0 | is.null(gagePathwayData()))
	  
       { selectInput("sigPathways", label = NULL, # h6("Funtional Category"), 
                  choices = list("All" = "All"), selected = "All")  }	 else { 
		choices <- "All"  # default, sometimes these methods returns "No significant pathway found"
		if( input$pathwayMethod == 1) { 
			if(!is.null(gagePathwayData())) 
				if(dim(gagePathwayData())[2] >1) 
				choices <- gagePathwayData()[,2] 
		} 
		if( input$pathwayMethod == 2) {
			if(!is.null(PGSEAplot.data()))  
				if(dim(PGSEAplot.data())[2] >1) 
					{ 	pathways <- as.data.frame( PGSEAplot.data())
						choices <- substr(rownames(pathways),10, nchar( rownames(pathways)) )
					}
					
		}
		if( input$pathwayMethod == 3) 
		{ 	if(!is.null(fgseaPathwayData())) 
			if(dim(fgseaPathwayData())[2] >1) 
				choices <- fgseaPathwayData()[,2] 
		} 
		if( input$pathwayMethod == 4) 
			if(!is.null(PGSEAplotAllSamples.data())) 
				if(dim(PGSEAplotAllSamples.data())[2] >1) {
					pathways <- as.data.frame( PGSEAplotAllSamples.data())
					choices <-  substr(rownames(pathways),10, nchar( rownames(pathways)) )
				
				}
				
		if( input$pathwayMethod == 5) {
			if(!is.null(ReactomePAPathwayData())) 
				if(dim(ReactomePAPathwayData())[2] >1) 
					choices <- ReactomePAPathwayData()[,2] 
		}
		
		selectInput("sigPathways", label="Select a pathway to show expression pattern of related genes on a heatmap or a KEGG pathway diagram:"
						,choices=choices)
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
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
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
	
	# find related samples	
	iz = findContrastSamples(input$selectContrast1, colnames(convertedData()),readSampleInfo(),
										input$selectFactorsModel,input$selectModelComprions, 
										factorReferenceLevels(),input$CountsDEGMethod ,
										input$dataFileFormat )
	x <-  convertedData()[which(rownames(convertedData()) %in% genes), iz ]
	if( input$selectOrg != "NEW") {
		ix = match( rownames(x), allGeneInfo()[,1])
		if( sum( is.na(allGeneInfo()$symbol )) != dim(allGeneInfo() )[1] )  # symbol really exists? 
		   rownames(x) <- paste(rownames(x),":", as.character( allGeneInfo()$symbol)[ix])
   }
	
	return( x )
	
     }) })
})


output$downloadSelectedPathwayData <- downloadHandler(
		# filename = function() {"Selected_Pathway_detail.csv"},
		filename = function() {paste(input$selectContrast1,"_method_",input$sigPathways,".csv",sep="")},
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
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$heatColors1; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
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
	
	# sometimes, a gene can be all zero in selected samples.
	x <- x[which(apply(x,1,sd)>0) ,]
	
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
    if (is.null(input$file1)&& input$goButton == 0)   return(blank)

	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO; tem = input$noIDConversion; tem=input$missingValue
	tem = input$selectContrast
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold	
	tem = input$sigPathways; 
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################
	
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

	if(is.null( input$selectGO ) ) return(blank)
	if(input$selectGO != "KEGG") return(blank)
	if(is.null(gagePathwayData() ) ) return(blank)
	if(is.null( input$sigPathways))  return (blank) 
	# if( is.null(selectedPathwayData()) ) return(blank)

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
	#cat("\nhere5  ",keggSpecies, " ",Species," ",input$sigPathways, "pathID:",pathID,"End", fold[1:5],names(fold)[1:5],"\n")
	#cat("\npathway:",is.na(input$sigPathways))
	#cat("\n",fold[1:5],"\n",keggSpecies,"\n",pathID)
    if(is.null(pathID) ) return(blank) # kegg pathway id not found.
	if(nchar(pathID)<3 ) return(blank)
	randomString <- gsub(".*file","",tempfile()) 
	tempFolder <- tempdir() # tempFolder = "temp";
	outfile <- paste( tempFolder,"/",pathID,".",randomString,".png",sep="")
	
	pv.out <- mypathview(gene.data = fold, pathway.id = pathID, kegg.dir = tempFolder,  out.suffix = randomString, species = keggSpecies, kegg.native=TRUE)
	if(0) {  # this works but sometimes messes up with current folder, shiny cannot find database files.
		wd = getwd()
		setwd(tempFolder)
		#pv.out <- mypathview(gene.data = fold, pathway.id = pathID, kegg.dir = tempFolder,  out.suffix = randomString, species = keggSpecies, kegg.native=TRUE)
		try( 
		 pv.out <- pathview(gene.data = fold, pathway.id = pathID, out.suffix = randomString, species = keggSpecies, kegg.native=TRUE)
		)
		setwd(wd)
	}
	
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
       width = "100%",
        height = "100%",
         alt = "KEGG pathway image.")
		}) 
	})
  }, deleteFile = TRUE)

# list of pathways with details
pathwayListData  <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem=input$limmaPval; tem=input$limmaFC
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold; tem =input$pathwayMethod
	if(is.null(input$selectGO ) ) return (NULL)
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem= input$selectFactorsModel;    tem= input$selectBlockFactorsModel; 
	tem= input$selectModelComprions;  tem= input$selectInteractions
	tem= input$referenceLevelFactor1; tem= input$referenceLevelFactor2;
	tem= input$referenceLevelFactor3; tem= input$referenceLevelFactor4; 
	tem= input$referenceLevelFactor5; tem= input$referenceLevelFactor6; 
	####################################
	pathways = NULL
	if( input$pathwayMethod == 1)  
		if(!is.null(gagePathwayData())) 
			if(dim(gagePathwayData())[2] >1) { 
				pathways <- gagePathwayData()
				colnames(pathways)[2] ="Pathways"; 	
				colnames(pathways)[4] ="nGenes"; 
			}
	if( input$pathwayMethod == 3) 
		if(!is.null(fgseaPathwayData())) 
			if(dim(fgseaPathwayData())[2] >1) {
				pathways <- fgseaPathwayData()
				colnames(pathways)[2] ="Pathways"; 	
				colnames(pathways)[4] ="nGenes"; 
			}
	
	if( input$pathwayMethod == 2) 
		if(!is.null(PGSEAplot.data())) 
			if(dim(PGSEAplot.data())[2] >1) {
				pathways <- as.data.frame( PGSEAplot.data())
				pathways$Pathways = substr(rownames(pathways),10, nchar( rownames(pathways)) )
				pathways$adj.Pval = gsub(" .*","", rownames(pathways))
				pathways$Direction ="Diff"
				
				}
	if( input$pathwayMethod == 4) 
		if(!is.null(PGSEAplotAllSamples.data())) 
			if(dim(PGSEAplotAllSamples.data())[2] >1) {
				pathways <- as.data.frame( PGSEAplotAllSamples.data())
				pathways$Pathways = substr(rownames(pathways),10, nchar( rownames(pathways)) )
				pathways$adj.Pval = gsub(" .*","", rownames(pathways))
				pathways$Direction ="Diff"
				
				}			
	if( is.null( pathways) ) return(NULL)
	
	# if no gene set data, return pathway list
	if(is.null(GeneSets() ) ) return(pathways) 
	

	pathways$adj.Pval = as.numeric(pathways$adj.Pval)

	if(nrow(pathways)>1)  # sometimes only one pathway is in the table
	for( i in 2:nrow(pathways) )
		if(nchar(pathways$Direction[i]) <=1)
			pathways$Direction[i] = pathways$Direction[i-1]
			
	# gene symbol matching symbols 
	probeToGene = NULL
	if( input$selectGO != "ID not recognized!" & input$selectOrg != "NEW")
	if(sum(is.na( allGeneInfo()$symbol ) )/ dim( allGeneInfo() )[1] <.5 ) { # if more than 50% genes has symbol
		probeToGene = allGeneInfo()[,c("ensembl_gene_id","symbol")]
		probeToGene$symbol = gsub(" ","",probeToGene$symbol)

		ix = which( is.na(probeToGene$symbol) |
					nchar(probeToGene$symbol)<2 | 
					toupper(probeToGene$symbol)=="NA" |  
					toupper(probeToGene$symbol)=="0"  ) 			
		probeToGene[ix,2] = probeToGene[ix,1]  # use gene ID

	}		

	
	pathways$Genes =""
	# looking up genes for each pathway
	for(i in 1:nrow(pathways) ){ 
		ix <- which(names(GeneSets() ) == pathways$Pathways[i]   ) # find the gene set
		if(length(ix) != 0 ) { 
			genes <- GeneSets()[[ix]] # retrieve genes
			
			if(!is.null(probeToGene) ) { 
				iy = match(genes,probeToGene[,1])
				genes = probeToGene[iy,2]
			}
			
			pathways$Genes[i] = paste(genes, collapse=" ")
		}
	}
	return(pathways)


})

output$downloadPathwayListData <- downloadHandler(
		# filename = function() {"Selected_Pathway_detail.csv"},
		filename = function() {paste(input$selectContrast1,"(",input$pathwayMethod,")",".csv",sep="")},
			content = function(file) {
			write.csv(pathwayListData(), file)
	    }
	)
	
output$enrichmentPlotPathway <- renderPlot({
    if(is.null(pathwayListData())) return(NULL)
	tem = input$removeRedudantSets
	enrichmentPlot(pathwayListData(), 45  )

}, height=600, width=800)

output$enrichmentPlotPathway4Download <- downloadHandler(
      filename = "enrichmentPlotPathway.eps",
      content = function(file) {
	  cairo_ps(file, width = 10, height = 6, pointsize = 10)
	  enrichmentPlot(pathwayListData(),41  )
        dev.off()
      })

output$enrichmentNetworkPlotPathway <- renderPlot({
    if(is.null(pathwayListData())) return(NULL)
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO; tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold; tem =input$pathwayMethod
	if(is.null(input$selectGO ) ) return (NULL)
	
	enrichmentNetwork(pathwayListData(),layout_change = input$layoutButton3 )
}, height=900, width=900)	  


# define a network. This code block is almost exactly the same as that for the network in DEG2
networkPA <- reactive({
    if(is.null(pathwayListData())) return(NULL)
    if(is.null( input$wrapTextNetworkPA )) return(NULL)

	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO; tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold; tem =input$pathwayMethod

    network <- pathwayListData()

    if(is.null( input$upORdownRegPA )) return(NULL)
    if(input$upORdownRegPA != "Both")
       network <- network[ grepl(input$upORdownRegPA, network$Direction), ]
    if(dim(network)[1] == 0) return(NULL)

    if(input$wrapTextNetworkPA)
      network$Pathways <- wrap_strings( network$Pathways ) # wrap long pathway names using default width of 30 10/21/19

    g <- enrichmentNetwork(network,layoutButton = input$layoutVisPA, edge.cutoff = input$edgeCutoffPA )

    data1 <- toVisNetworkData(g)
    
    # Color codes: https://www.rapidtables.com/web/color/RGB_Color.html
    data1$nodes$shape <- "dot"
    # remove the color change of nodes
    #data1$nodes <- subset(data1$nodes, select = -color)
    
    data1$nodes$size <- 5 + data1$nodes$size^2 
    
    return(data1)
})

  # note the same code is used twice as above. They need to be updated together!!!	  
output$visNetworkPA <- renderVisNetwork({
    if(is.null(pathwayListData())) return(NULL)
    if(dim(pathwayListData())[1] == 1) return(NULL)
    if(is.null( input$wrapTextNetworkPA )) return(NULL)
    if(is.null(networkPA() )) return(NULL)

	tem = input$selectOrg ; #tem = input$listComparisonsPathway
	tem = input$selectGO; tem = input$selectContrast1
	tem = input$minSetSize; tem = input$maxSetSize; tem=input$pathwayPvalCutoff; 
	tem=input$nPathwayShow; tem=input$absoluteFold; tem =input$pathwayMethod


    data1 <- networkPA()

    # render network of pathways
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
 
output$visNetworkPADownload <- downloadHandler(
    filename = "Pathway_Overlaping_Network_DEG.html",
    content = function(file) {

    data1 <- networkPA()
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
      )  %>% 
        visSave(file = file, background = "white")
      
    })  


output$downloadNodesPA <- downloadHandler(
    filename = function() {"Pathway_network_nodes.csv"},
    content = function(file) {      
      write.csv(networkPA()$nodes, file, row.names=FALSE)
    }
  )

output$downloadEdgesDEG <- downloadHandler(
    filename = function() {"Pathway_network_edges.csv"},
    content = function(file) {    
      write.csv(networkPA()$edges, file, row.names=FALSE)
    }
  )  

################################################################
#   Chromosome
################################################################

# visualizing fold change on chrs. 
output$genomePlotly <- renderPlotly({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)		
		tem = input$selectOrg ; 
		tem = input$selectContrast2
		if (is.null(input$selectContrast2 ) ) return(NULL)
		if( input$selectOrg == "NEW" | ncol(allGeneInfo() )==1 ) return(NULL)
		if( length(limma()$topGenes) == 0 ) return(NULL)
		library(dplyr)
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
		tem = input$CountsDEGMethod;
		tem=input$limmaPvalViz; 
		tem=input$limmaFCViz
		tem = input$submitModelButton 
		####################################

		tem = input$MAwindowSize
        tem = input$MAwindowSteps
        tem = input$MAwindowCutoff
        tem = input$ignoreNonCoding
        tem = input$chRegionPval
        tem = input$labelGeneSymbol

		####################################
		
	  isolate({ 
		withProgress(message=sample(quotes,1), detail ="Visualzing expression on the genome", {
		# default plot
		fake = data.frame(a=1:3,b=1:3)
		p <- ggplot(fake, aes(x = a, y = b)) +
							 geom_blank() + ggtitle("No genes with position info.") +
							 theme(axis.title.x=element_blank(),axis.title.y=element_blank())

		if(length( limma()$comparisons)  ==1 )  {
			top1=limma()$topGenes[[1]]  
		} else {
		  top = limma()$topGenes
		  ix = match(input$selectContrast2, names(top))
		  if( is.na(ix)) return (ggplotly(p))
		  top1 <- top[[ix]]; 
		}
		  if(dim(top1)[1] == 0 ) return (ggplotly(p))
		  colnames(top1)= c("Fold","FDR")

		 x <- merge(top1,allGeneInfo(), by.x="row.names",by.y="ensembl_gene_id"  )

         colnames(x)[which(colnames(x) == "Row.names")] <- "ensembl_gene_id"

      
        # only coding genes? 
        if(input$ignoreNonCoding) {
          x <- subset(x, gene_biotype == "protein_coding")
        }


        incProgress(0.1)
		 # if no chromosomes found. For example if user do not convert gene IDs.
		 if( dim(x)[1] >5  ) { 

           x <- x[order(x$chromosome_name,x$start_position),]
  
           x$ensembl_gene_id <- as.character( x$ensembl_gene_id)
   
           # if symbol is missing use Ensembl id
           x$symbol = as.character(x$symbol)  
           ix = which(is.na(x$symbol))
           ix2 = which(nchar(as.character(x$symbol))<= 2 )
           ix3 = which( duplicated(x$symbol))
           ix = unique( c(ix,ix2,ix3))
           x$symbol[ix] <- x$ensembl_gene_id[ix] 

           x = x[!is.na(x$chromosome_name),]
           x = x[!is.na(x$start_position),]
            


             tem = sort( table( x$chromosome_name), decreasing=T)
             ch <- names( tem[tem >= 1 ] )  # ch with less than 100 genes are excluded
             if(length(ch) > 50) ch <- ch[1:50]  # at most 50 ch
             ch <- ch[ nchar(ch)<=12] # ch. name less than 10 characters
             ch = ch[order(as.numeric(ch) ) ]
             tem <- ch
             ch <- 1:(length(ch))  # the numbers are continous from 1 to length(ch)
             names(ch) <- tem  # the names are real chr. names


             x <- x[which(x$chromosome_name %in% names(ch)),]
             x <- droplevels(x)

             x$chNum <- 1 # numeric encoding
             x$chNum <- ch[ x$chromosome_name ]

            # use max position as chr. length   before filtering
           chLengthTable = aggregate(start_position~chromosome_name, data=x,max )
              # add chr. numer 
             chLengthTable$chNum <-  ch[ chLengthTable$chromosome_name ]
             chLengthTable <- chLengthTable[!is.na( chLengthTable$chNum ), ]
             chLengthTable <- chLengthTable[order(chLengthTable$chNum), c(3,2)]
             chLengthTable <- chLengthTable[order(chLengthTable$chNum), ]
             chLengthTable$start_position <- chLengthTable$start_position/1e6

  
           # only keep significant genes
          ix = which( (x$FDR< as.numeric(input$limmaPvalViz)) &
                        (abs(x$Fold) > log2( as.numeric(input$limmaFCViz) ) ) )

           if (length(ix) > 5) { 

             # remove nonsignificant / not selected genes
             x0 <- x   # keep a copy
             x = x[ix, ]  



              # prepare coordinates
             x$start_position = x$start_position/1000000 # Mbp
             chD = 30 # distance between chs.
             foldCutoff = 4   # max log2 fold 
    
             # 
             x$Fold[which(x$Fold > foldCutoff )] = foldCutoff   # log2fold within -5 to 5
             x$Fold[which(x$Fold <   -1*foldCutoff )] = -1*foldCutoff 
             x$Fold = 4 * x$Fold
    
             x$y = x$chNum*chD + x$Fold
             chTotal = dim(chLengthTable)[1] 
             x$R = as.factor(sign(x$Fold))
    
             colnames(x)[ which(colnames(x) == "start_position")] = "x"

             incProgress(0.3)
             # plotting ----------------------------------

             p <- ggplot() +  # don't define x and y, so that we could plot use two datasets
                  geom_point(data = x, aes(x = x, y = y, colour = R, text = symbol), shape = 20, size = 0.2 ) 

            if(input$labelGeneSymbol)
                p <- p + geom_text(data = x, aes(x = x, y = y, label = symbol),
                                    check_overlap = FALSE, angle = 45, size = 2, vjust = 0, nudge_y = 4 )
             #label y with ch names
             p <- p +  scale_y_continuous(labels = paste("chr", names(ch[chLengthTable$chNum]),sep=""), 
                                          breaks = chD* (1:chTotal), 
                                          limits = c(0, chD*(chTotal + 1) + 5) )
             # draw horizontal lines for each ch.
             for( i in 1:dim(chLengthTable)[1] )
               p = p+ annotate( "segment",x = 0, xend = chLengthTable$start_position[i],
                                y = chLengthTable$chNum[i]*chD, yend = chLengthTable$chNum[i]*chD)
             # change legend  http://ggplot2.tidyverse.org/reference/scale_manual.html
             p <- p + scale_colour_manual(name="",   # customize legend text
                                          values=c("red", "blue"),
                                          breaks=c("1","-1"),
                                          labels=c("Up", "Dn")) 
             p <- p + xlab("Position on chrs. (Mbp)") +  theme(axis.title.y=element_blank())      
             p <- p + theme(legend.position="none")

             incProgress(0.5)
             # add trend lines------------------------------------------
             x0 <- x0[x0$chromosome_name %in% unique(x$chromosome_name), ]
             x0$chNum <- 1 # numeric encoding
             x0$chNum <- ch[ x0$chromosome_name ]
             x0$start_position = x0$start_position/1e6 # Mbp             

             windowSize = as.numeric( input$MAwindowSize )#Mb            
             steps = as.numeric( input$MAwindowSteps ) # step size is then windowSize / steps       
             cutoff <- as.numeric(input$MAwindowCutoff) 

             x0$Fold <-  x0$Fold - mean(x0$Fold) # centering             

             incProgress(0.6)
                             
             for(i in 0:(steps-1)) {
               #step size is  windowSize/steps   
               # If windowSize=10 and steps = 2; then step size is 5Mb
               # 1.3 becomes 5, 11.2 -> 15 for step 1
               # 1.3 -> -5
               x0$x <- ( floor((x0$start_position - i * windowSize / steps)/ windowSize )  
                        + 0.5 + i / steps ) * windowSize

               movingAverage1 <- x0 %>%
                 select(chNum, x, Fold) %>%
                 filter( x >= 0) %>%   # beginning bin can be negative for first bin in the 2nd step
                 group_by(chNum, x) %>%
                 summarize( ma = mean(Fold),
                            n = n(),
                            pval = ifelse( n() >= 3 && sd(Fold) > 0, t.test(Fold)$p.value, 0 ) ) %>%
                 filter(!is.na(pval)) # na when only 1 data point?

               if(i == 0) {
                 movingAverage <- movingAverage1
               } else {
                 movingAverage <- rbind(movingAverage, movingAverage1)  
               }        
             }
 
            
             # translate fold to y coordinates
             movingAverage <- movingAverage %>%
                filter(n >= 3) %>%
#                mutate( pval = p.adjust(pval, method = "fdr", n = length(pval[pval < 0.9]) ) )
                mutate( pval = p.adjust(pval, method = "fdr" ) ) %>%
                filter( pval < as.numeric(input$chRegionPval) ) %>%
                mutate( y = ifelse(ma > 0, 1, -1)) %>% # upper bound
                mutate(y = chNum * chD + 3 * y) %>%
                mutate( ma = ifelse(ma > 0, 1, -1)) %>%
                mutate( ma = as.factor(ma))

              # significant regions are marked as horizontal error bars 
             if(dim(movingAverage)[1] > 0) {
               p <- p +
                 geom_errorbarh(data = movingAverage, aes(x = x, 
                                                          y = y, 
                                                          xmin = x -windowSize/2, 
                                                          xmax = x + windowSize/2,
                                                          colour = ma), 
                                 size = 2, 
                                 height = 15 )

                 # label significant regions
                 sigCh <- sort(table(movingAverage$chNum), decreasing = TRUE)
                 sigCh <- names(ch)[ as.numeric(names(sigCh)) ]
                 if(length(sigCh) <= 5) { # more than 5 just show 5
                   sigCh <- paste0("chr", sigCh, collapse = ", ")
                 } else {
                   sigCh <- sigCh[1:5]
                   sigCh <- paste0("chr", sigCh, collapse = ", ")                  
                   sigCh <- paste0(sigCh,", ...")
                 }

                 sigCh <- paste(dim(movingAverage)[1], 
                                " enriched regions \n(",
                                round( sum(chLengthTable$start_position)/ windowSize * steps * as.numeric(input$chRegionPval), 2),
                                          " expected)  detected on:\n ", sigCh)
                 
               p <- p + annotate(geom = "text", 
                          x = max(x$x) * 0.70,
                          y = max(x$y) * 0.90,
                          label = sigCh)


             }

         } # have genes after filter
			
      }  # have 5+ genes to begin with
              incProgress(1)
	  ggplotly(p)
    }) # progress
  }) # isolate
})

 
	  
# pre-calculating PREDA, so that changing FDR cutoffs does not trigger entire calculation
genomePlotDataPre <- reactive({
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	if( input$selectOrg == "NEW" | ncol(allGeneInfo() )==1 ) return(NULL)
	tem = input$selectOrg ; 
	tem = input$selectContrast2
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem = input$submitModelButton 
	####################################
	
  isolate({ 
	withProgress(message=sample(quotes,1), detail ="Identifying differentially expressed genomic regions using PREDA", {
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
	if( input$selectOrg == "NEW" | ncol(allGeneInfo() )==1 ) return(NULL)	
	tem = input$selectOrg ; 
	tem = input$selectContrast2
	tem = input$StatisticCutoff
	tem = input$RegionsPvalCutoff
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem = input$submitModelButton
	####################################
	
  isolate({ 
	withProgress(message=sample(quotes,1), detail ="Identifying differentially expressed genomic regions using PREDA", {
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
	library(PREDA,verbose=FALSE)  # showing expression on genome
	library(PREDAsampledata,verbose=FALSE) 
	library(hgu133plus2.db,verbose=FALSE)
	if( input$selectOrg == "NEW" | ncol(allGeneInfo() )==1 ) return(NULL)	
    if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
    if( is.null( genomePlotData() ) ) return(NULL)
	tem = input$selectOrg ; 
	tem = input$selectContrast2
	tem = input$StatisticCutoff
	tem = input$RegionsPvalCutoff
	
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem = input$submitModelButton
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
	if( input$selectOrg == "NEW" | ncol(allGeneInfo() )==1 ) return(NULL)
  	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem = input$submitModelButton
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
	if( input$selectOrg == "NEW" | ncol(allGeneInfo() )==1 ) return(NULL)  
 	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem = input$submitModelButton
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
	if( input$selectOrg == "NEW" | ncol(allGeneInfo() )==1 ) return(NULL)  
	##################################  
	# these are needed to make it responsive to changes in parameters
	tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
	if( !is.null(input$dataFileFormat) ) 
    	if(input$dataFileFormat== 1)  
    		{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
	if( !is.null(input$dataFileFormat) )
    	if(input$dataFileFormat== 2) 
    		{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
	tem = input$CountsDEGMethod;
	tem = input$submitModelButton
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
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
		tem = input$CountsDEGMethod;
		
		tem = input$nGenesBiclust
		tem = input$biclustMethod
		####################################   
   
		isolate({
			withProgress(message="Runing biclustering", {
			library(biclust,verbose=FALSE)
			

			if(input$biclustMethod == "BCQU()" )
				library(QUBIC,verbose=FALSE) # have trouble installing on Linux
			if(input$biclustMethod == "BCUnibic()" )
				library(runibic,verbose=FALSE) # 
			
		
			x <- convertedData()
			n=input$nGenesBiclust	
			if(n>dim(x)[1]) n = dim(x)[1] # max	as data
			if(n<10) n = 10 
			if(n> 2000 ) n = 2000 			
			x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
			
			if( input$biclustMethod == "BCXmotifs()"  )
				x<-discretize(x)
			incProgress(1/2)
			#res <- biclust::biclust(as.matrix( x), method = BCQU()) 
			runR = paste( "res <- biclust::biclust(as.matrix( x), method =", input$biclustMethod ,")" )
			eval(parse(text = runR ) )
			
			incProgress(1)
			return(list( x=x, res=res)	)	 
			})
		
		} )
   
   
   
   } )
   
   
output$listBiclusters <- renderUI({
		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
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
		if (is.null(input$file1) && input$goButton == 0)   return(NULL)		
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
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
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

		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO4
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart

		if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
		NoSig = as.data.frame("No significant enrichment found.")
		
		if( is.null(biclustering() ) ) return(NULL)
		if( is.null(input$selectBicluster ) ) return(NULL)
		tem = input$nGenesBiclust	
		tem = input$biclustMethod
		tem = input$removeRedudantSets
		
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
				if(input$removeRedudantSets) reduced = redudantGeneSetsRatio else reduced = FALSE
				result = FindOverlap (convertedID,allGeneInfo(), input$selectGO4,input$selectOrg,1, reduced) }
				
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
  }, digits = 0,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T, 
     sanitize.text.function = function(x) x)
	#   output$selectedHeatmap <- renderPlot({       hist(rnorm(100))    })

	
output$geneListBicluster <- renderTable({
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if( is.null( input$selectGO4) ) return (NULL)
		if( input$selectGO4 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.")) #No matching species

		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO4
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart

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
		  # if  new species
		if(  input$selectGO4 == "ID not recognized!" | 
		     input$selectOrg == "NEW" | 
		     dim(allGeneInfo())[1] == 1) {
			top1 = as.data.frame(rownames(top1))
			colnames(top1)="Genes"
		} else {	# add gene info	
		  
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
		}
		  return(top1)
	
  }, digits = 0,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T)

  
biclustData <- reactive({
  		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if( is.null( input$selectGO4) ) return (NULL)
		if( input$selectGO4 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.")) #No matching species

		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO4
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart

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
#   Co-expression network by WGCNA 
################################################################
 
wgcna <- reactive ({
	  if (is.null(input$file1) && input$goButton == 0)   return(NULL)

		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
		tem = input$CountsDEGMethod;
		
		tem = input$mySoftPower;
		tem = input$nGenesNetwork
		tem = input$minModuleSize
		
		####################################   
   
		isolate({
			withProgress(message="Constructing co-expression network using WGCNA", {
			
			#http://pklab.med.harvard.edu/scw2014/WGCNA.html
			library(WGCNA,verbose=FALSE)
			incProgress(1/3)
			x <- convertedData()
			n=input$nGenesNetwork	
			if(n>dim(x)[1]) n = dim(x)[1] # max	as data
			if(n<50) return(NULL)
			if(dim(x)[2] <4) return(NULL)
			if(n> maxGeneWGCNA ) n = maxGeneWGCNA 			
			#x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
			
			datExpr=t(x[1:n,])

			subGeneNames=colnames(datExpr)
			
			#Choosing a soft-threshold to fit a scale-free topology to the network
			powers = c(c(1:10), seq(from = 12, to=20, by=2));
			sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")

			incProgress(1/2)
			softPower = input$mySoftPower;

			#calclute the adjacency matrix
			adj= adjacency(datExpr,type = "unsigned", power = softPower);

			#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
			TOM=TOMsimilarityFromExpr(datExpr,networkType = "unsigned", TOMType = "unsigned", power = softPower);
			colnames(TOM) =subGeneNames
			rownames(TOM) =subGeneNames
			
			##########################################
			# module detection

			library(flashClust,verbose=FALSE)

			geneTree = flashClust(as.dist(1-TOM),method="average");

			# Set the minimum module size
			minModuleSize = input$minModuleSize;

			# Module identification using dynamic tree cut
			dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
			#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
			# table(dynamicMods)
			
			dynamicColors = labels2colors(dynamicMods)
			#table(dynamicColors)

			#discard the unassigned genes, and focus on the rest  # this causes problems
			#restGenes= (dynamicColors != "grey")
			#diss1=1-TOMsimilarityFromExpr(datExpr[,restGenes], power = softPower)
			#colnames(diss1) =rownames(diss1) =subGeneNames[restGenes]
			
			#colorCode = as.character(dynamicColors[restGenes])

			moduleInfo = cbind( subGeneNames, dynamicColors, dynamicMods)
			moduleInfo = moduleInfo[which(moduleInfo[,2] != "grey") ,] # remove genes not in any modules
			moduleInfo = moduleInfo[order(moduleInfo[,3]),] # sort
			
			n.modules = length(unique(dynamicColors) ) -1 ; nGenes = dim(moduleInfo)[1]	
			
			return(list(x = t(datExpr),powers=powers,sft=sft, TOM = TOM, dynamicColors = dynamicColors, moduleInfo = moduleInfo,n.modules=n.modules, nGenes =nGenes) )

			
			}) #progress
		}) # isolate
	})

output$moduleStatistics <- renderText({
		if(is.null(wgcna() ) ) return(NULL)
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2; tem =input$NminSamples2; tem =input$NminSamples2; tem =input$NminSamples2 }
		tem = input$CountsDEGMethod;
		tem = input$mySoftPower;
		
		####################################  
		paste( "A network of", wgcna()$nGenes,"genes was divided into ",wgcna()$n.modules, "modules." )
	})	

	
output$softPower <- renderPlot({
		if(is.null(wgcna() ) ) return(NULL)
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter }
		tem = input$CountsDEGMethod;
		tem = input$mySoftPower;
		
		####################################  
		
		sft = wgcna()$sft;
		powers=wgcna()$powers;
		# Plot the results
		#sizeGrWindow(9, 5)
		par(mfrow = c(1,2));
		cex1 = 0.9;

		# Scale-free topology fit index as a function of the soft-thresholding power
		plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
		text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

		# Red line corresponds to using an R^2 cut-off
		abline(h=0.80,col="red")

		# Mean connectivity as a function of the soft-thresholding power
		plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
		text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
	})

	
output$modulePlot <- renderPlot({
		if(is.null(wgcna() ) ) return(NULL)
		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
		tem = input$CountsDEGMethod;
		tem = input$mySoftPower;
		tem = input$nGenesNetwork		
		tem = input$minModuleSize
		####################################  
		
		diss1 = 1-wgcna()$TOM;
		dynamicColors = wgcna()$dynamicColors
		
		hier1=flashClust(as.dist(diss1), method="average" )
		#set the diagonal of the dissimilarity to NA 
		diag(diss1) = NA;
		plotDendroAndColors(hier1, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

		#Visualize the Tom plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
		#TOMplot(diss1, hier1, colorCode)			


		})
		
		
output$networkHeatmap <- renderPlot({ 
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)

		##################################  
		# these are needed to make it responsive to changes in parameters
		tem = input$selectOrg;  tem = input$dataFileFormat; tem = input$heatColors1; tem = input$noIDConversion; tem=input$missingValue
		if( !is.null(input$dataFileFormat) ) 
			if(input$dataFileFormat== 1)  
				{  tem = input$minCounts ; tem= input$NminSamples;tem = input$countsLogStart; tem=input$CountsTransform }
		if( !is.null(input$dataFileFormat) )
			if(input$dataFileFormat== 2) 
				{ tem = input$transform; tem = input$logStart; tem= input$lowFilter; tem =input$NminSamples2 }
		####################################
		tem = input$minModuleSize
		tem = input$mySoftPower;
		tem = input$nGenesNetwork
		if( is.null(wgcna()) ) return(NULL)
		withProgress(message="Creating heatmap", {
			x <- convertedData()
			n=input$nGenesNetwork	
			if(n>dim(x)[1]) n = dim(x)[1] # max	as data
			if(n<50) return(NULL)
			if(n> 2000 ) n = 2000
			
			x = x[1:n,]
			
			x2 = merge(x,wgcna()$moduleInfo, by.x="row.names",by.y="subGeneNames",all.y=TRUE )
			#x2 = merge(x,moduleInfo, by.x="row.names",by.y="subGeneNames",all.y=TRUE )		
			x2 = x2[order(x2$dynamicMods),]
			bar = as.numeric( x2$dynamicMods )

			x2 = x2[,colnames(x)]
			
			clusterNames = paste("Module",1:max(bar))
			sideColors= standardColors( max(bar) )# WGCNA function for module colors
	   
		myheatmap2(  x2-apply(x2,1,mean), bar,1000,mycolor=input$heatColors1, clusterNames=clusterNames, sideColors= sideColors)
		
		incProgress(1, detail = paste("Done")) }) #progress 
  } ,height = 500,width = 500)

  
moduleData <- reactive({
  		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)

		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO4
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart

		if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
		
		if( is.null(wgcna()) ) return(NULL)
		tem = input$mySoftPower;
		tem = input$nGenesNetwork		
		tem = input$minModuleSize
		
	 isolate({ 
			
			x2 = merge(wgcna()$moduleInfo, wgcna()$x, by.y="row.names",by.x="subGeneNames",all.x=TRUE )
			x2 = x2[order(x2$dynamicMods),]
			
			if( input$selectGO5 != "ID not recognized!" & input$selectOrg != "NEW")
			if( !is.null(allGeneInfo() )   ) {
				x2 <- merge(x2, allGeneInfo(), by.x ="subGeneNames", by.y="ensembl_gene_id",all.x=T )
				rownames(x2)= paste0(x2$symbol,"__",  x2$subGeneNames )
				#x1 = x2[, 2:(dim(x1)[2]+1) ]
			}
			x2 = x2[order(x2$dynamicMods),]
			colnames(x2)[1:3] = c("gene id","Module color","Module#")
		   
			
			
			return(x2)	
		})
	})

	
output$download.WGCNA.Module.data <- downloadHandler(
		filename = function() {"WGCNA_modules.csv"},
			content = function(file) {
			write.csv(moduleData(), file)
	    }
	)	
	
	
output$networkModuleGO <- renderTable({		
		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		if( is.null( input$selectGO5) ) return (NULL)
		if( input$selectGO5 == "ID not recognized!" ) return ( as.matrix("Gene ID not recognized.")) #No matching species

		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; tem = input$selectGO5
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart

		if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
		NoSig = as.data.frame("No significant enrichment found.")
		

		tem = input$mySoftPower;
		tem = input$nGenesNetwork		
		tem = input$minModuleSize
		tem = input$selectWGCNA.Module
		tem = input$removeRedudantSets
		if( is.null(wgcna()) ) return(NULL)
		
		isolate({
			withProgress(message="GO Enrichment", {
	
			#module = gsub(".* ","",input$selectWGCNA.Module)
			module = unlist(strsplit(input$selectWGCNA.Module," " ) )[2]	
			moduleColors = wgcna()$dynamicColors 
			inModule = (moduleColors==module);
			
			if( input$selectWGCNA.Module == "Entire network")
			  inModule = rep(TRUE, length(inModule))

			probes = rownames(wgcna()$x   )
			query  = probes[inModule];
		
			if( length(query)  <= minGenesEnrichment) return(NoSig) 

			
			if(input$selectOrg == "NEW" && !is.null( input$gmtFile) ){
				result <- findOverlapGMT( query, GeneSets(),1) 
			} else  { 
				convertedID <- converted()
				convertedID$IDs <- query
				if(input$removeRedudantSets) reduced = redudantGeneSetsRatio else reduced = FALSE
				result = FindOverlap (convertedID,allGeneInfo(), input$selectGO5,input$selectOrg,1, reduced) }
				
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
  }, digits = 0,spacing="s",striped=TRUE,bordered = TRUE, width = "auto",hover=T,
sanitize.text.function = function(x) x)
	#   output$selectedHeatmap <- renderPlot({       hist(rnorm(100))    })


output$listWGCNA.Modules <- renderUI({
		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem = input$mySoftPower;
		tem = input$nGenesNetwork		
		tem = input$minModuleSize
		if (is.null(wgcna() ) ){ # if sample info is uploaded and correctly parsed.
		   return(NULL)	   
		} else { 
			if( dim(wgcna()$moduleInfo)[1] == 0  ) {  # if no module
				return(NULL) 
			}	else  { 
					modules = unique(wgcna()$moduleInfo[, c("dynamicMods","dynamicColors")] )
					moduleList = apply(modules,1,paste,collapse=". ")
					moduleList  = paste0( moduleList, " (", table(wgcna()$moduleInfo[,"dynamicMods"] )," genes)"  )
					moduleList = c(moduleList,"Entire network")
					selectInput(	"selectWGCNA.Module", 
									label="Select a module",
									choices= moduleList
								)

				}

		}			
	}) 

	
exportModuleNetwork <- reactive({
  		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)

		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; 
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart

		if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
		
		if( is.null(wgcna()) ) return(NULL)
		tem = input$mySoftPower;
		tem = input$nGenesNetwork		
		tem = input$minModuleSize
		tem = input$selectWGCNA.Module
		tem=input$edgeThreshold 
		tem=input$topGenesNetwork;
		
	 isolate({ 
		outfile <- tempfile(fileext='.txt')
		withProgress(message="Extracting module", {

		#module = gsub(".* ","",input$selectWGCNA.Module)
		module = unlist(strsplit(input$selectWGCNA.Module," " ) )[2]		
		
		moduleColors = wgcna()$dynamicColors 
		inModule = (moduleColors==module);
		
		if( input$selectWGCNA.Module == "Entire network")
		  inModule = rep(TRUE, length(inModule))
		
		datExpr = t(wgcna()$x )
		probes = colnames(datExpr)
		modProbes = probes[inModule];
		
		modTOM = wgcna()$TOM[inModule, inModule];
		dimnames(modTOM) = list(modProbes, modProbes)

		nTop = input$topGenesNetwork;
		
		if( nTop > 1000) nTop = 1000; 
		
		IMConn = softConnectivity(datExpr[, modProbes]);
		top = (rank(-IMConn) <= nTop)
		
		# adding symbols 
		probeToGene = NULL
		if( input$selectGO5 != "ID not recognized!" & input$selectOrg != "NEW")
	    if(sum(is.na( allGeneInfo()$symbol ) )/ dim( allGeneInfo() )[1] <.5 ) { # if more than 50% genes has symbol
			probeToGene = allGeneInfo()[,c("ensembl_gene_id","symbol")]
			probeToGene$symbol = gsub(" ","",probeToGene$symbol)
			ix = which( is.na(probeToGene$symbol) |
						nchar(probeToGene$symbol)<2 | 
						toupper(probeToGene$symbol)=="NA" |  
						toupper(probeToGene$symbol)=="0"  ) 			
			probeToGene[ix,2] = probeToGene[ix,1]  # use gene ID

		}		
	
		incProgress(1/2, "Writing to file")
		if(input$noIDConversion) { 
			vis = exportNetworkToVisANT(modTOM[top,top],
				file = outfile,
				weighted = TRUE,
				threshold = input$edgeThreshold )	
		} else 				
			vis = exportNetworkToVisANT(modTOM[top,top],
				file = outfile,
				weighted = TRUE,
				threshold = input$edgeThreshold,
				probeToGene = probeToGene )
		
		return(outfile)	
		incProgress(1,"Done")
		})
		})
	})
	
	
output$downloadSelectedModule <- downloadHandler(
	  filename <- function() {
	  paste0("Module",input$selectWGCNA.Module,".txt")
		
	  },
	  content <- function(file) {
		file.copy(exportModuleNetwork() , file)
	  },
	  contentType = "text file"
	)
	  	
output$moduleNetwork <- renderPlot({
  		if (is.null(input$file1)&& input$goButton == 0)   return(NULL)

		tem = input$selectOrg; tem = input$noIDConversion; tem=input$missingValue
		tem=input$limmaPval; tem=input$limmaFC; tem = input$selectContrast; 
		tem = input$CountsDEGMethod; tem = input$countsLogStart; tem = input$CountsTransform
		tem = input$minCounts;tem= input$NminSamples; tem = input$lowFilter; tem =input$NminSamples2; tem=input$transform; tem = input$logStart

		if(input$selectOrg == "NEW" && is.null( input$gmtFile) ) return(NULL) # new but without gmtFile
		
		if( is.null(wgcna()) ) return(NULL)
		tem = input$mySoftPower;
		tem = input$nGenesNetwork		
		tem = input$minModuleSize
		tem = input$selectWGCNA.Module
		tem=input$edgeThreshold 
		tem=input$topGenesNetwork;
		tem = input$networkLayout
		if(is.null(input$selectWGCNA.Module ) ) return(NULL)
		
	 isolate({ 
		outfile <- tempfile(fileext='.txt')
		withProgress(message="Extracting module", {

		#module = gsub(".* ","",input$selectWGCNA.Module)
		module = unlist(strsplit(input$selectWGCNA.Module," " ) )[2]
		
		moduleColors = wgcna()$dynamicColors 
		inModule = (moduleColors==module);
		
		if( input$selectWGCNA.Module == "Entire network")
		  inModule = rep(TRUE, length(inModule))
		
		datExpr = t(wgcna()$x )
		probes = colnames(datExpr)
		modProbes = probes[inModule];
		
		modTOM = wgcna()$TOM[inModule, inModule];
		dimnames(modTOM) = list(modProbes, modProbes)

		nTop = input$topGenesNetwork;
		
		if( nTop > 1000) nTop = 1000; 
		
		IMConn = softConnectivity(datExpr[, modProbes]);
		top = (rank(-IMConn) <= nTop)
		
		# adding symbols 
		probeToGene = NULL
		if( input$selectGO5 != "ID not recognized!" 
		    & input$selectOrg != "NEW"
		    & dim(allGeneInfo())[1] > 1)
	    if(sum(is.na( allGeneInfo()$symbol ) )/ dim( allGeneInfo() )[1] <.5 ) { # if more than 50% genes has symbol
			probeToGene = allGeneInfo()[,c("ensembl_gene_id","symbol")]
			probeToGene$symbol = gsub(" ","",probeToGene$symbol)

			ix = which( is.na(probeToGene$symbol) |
						nchar(probeToGene$symbol)<2 | 
						toupper(probeToGene$symbol)=="NA" |  
						toupper(probeToGene$symbol)=="0"  ) 			
			probeToGene[ix,2] = probeToGene[ix,1]  # use gene ID

		}
		
		net <- modTOM[top,top] >input$edgeThreshold	

		
		for( i in 1:dim(net)[1])  # remove self connection
			net[i,i] = FALSE
		if(!is.null(probeToGene) & !input$noIDConversion) { # if gene symbol exist
			ix = match( colnames(net), probeToGene[,1])		
			colnames(net) = probeToGene[ix,2]
			ix = match( rownames(net), probeToGene[,1])		
			rownames(net) = probeToGene[ix,2]		
		}
		library(igraph,verbose=FALSE)
		#plot(graph_from_data_frame(d=data.frame(1:10,ncol=2)  ,directed=F) )
		# http://www.kateto.net/wp-content/uploads/2016/01/NetSciX_2016_Workshop.pdf
		plot( graph_from_adjacency_matrix( net, mod ="undirected" ), 
				vertex.label.color="black", vertex.label.dist=3,vertex.size=7)
		
		if(0){ 
		incProgress(1/2, "Writing to file")
		vis = exportNetworkToVisANT(,
			file = outfile,
			weighted = TRUE,
			threshold = input$edgeThreshold,
			probeToGene = probeToGene )
		
		return(outfile)
		}
		incProgress(1,"Done")
		})
		})
	})
	
	
################################################################
#   Session Info
################################################################
  
	# output user settings and session info
output$RsessionInfo <- renderUI({
	if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
		
	i =""
	i = c(i,"We recommend users to save the following details about their analyses.")
	i = c(i,paste("Analyses of read count data were conducted using  ",iDEPversion, 
		", hosted at http://ge-lab.org/idep/ on ", date(),".",sep="") )	
	i = c(i,"<a href=\"https://github.com/iDEP-SDSU/idep\"target=\"_blank\"> Source code on Github.</a>")
	
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
	# output user settings and session info

output$downloadConvertedCountsRtab <- downloadHandler(
		filename = function() {"Downloaded_Converted_Data.csv"},
		content = function(file) {
      write.csv( processedCountsData(), file, row.names=FALSE)	    
	})
output$downloadProcessedDataRtab <- downloadHandler(
		filename = function() {"Downloaded_Converted_Data.csv"},
		content = function(file) {
      write.csv( convertedData(), file)	    
	}) 
output$downloadRcode <- downloadHandler(
		filename = function() {"iDEP_Rcode.R"},
		content = function(file) {
      write( Rcode(), file)	    
	})
	
# find the corresponding file names for gene info, for download	
geneInfoFileName <- reactive({
  	if (is.null(input$file1)&& input$goButton == 0)   return(NULL)

	ix = grep(converted()$species[1,1],geneInfoFiles)
	if (length(ix) == 0 ) {return(NULL)} else {
		# If selected species is not the default "bestMatch", use that species directly
		if(input$selectOrg != speciesChoice[[1]]) {  
			ix = grep(findSpeciesById(input$selectOrg)[1,1], geneInfoFiles )
		}
		if(length(ix) == 1)  # if only one file           #WBGene0000001 some ensembl gene ids in lower case
			{x = as.character(geneInfoFiles[ix]); 	return(x) }else # read in the chosen file 
			return(NULL )   
	}

})	
	
output$downloadGeneInfo <- downloadHandler(
  filename <- function() {
    gsub(".*/","",geneInfoFileName())
  },
  content <- function(file) {
    file.copy(as.character(geneInfoFileName()), file)
  }
)
# sampleInfo file name
sampleInfoFileName <- reactive({
  	if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	if( is.null(input$file2) && !is.null( readData()$sampleInfoDemo ) ) # if using demo data
		return(demoDataFile2) else 
	if( is.null(input$file2) )
		return(NULL)  else
	{ 	
	  inFile <- input$file2    # user data
	  inFile <- inFile$datapath
	  return(inFile)		
	}

})	
output$downloadSampleInfoFile <- downloadHandler(
  filename <- function() {
    "Downloaded_sampleInfoFile.csv"
  },
  content <- function(file) {
    file.copy(sampleInfoFileName(), file)
  }
)
# find the corresponding file names for gene info, for download	
pathwayFileName <- reactive({
  	if (is.null(input$file1)&& input$goButton == 0)   return(NULL)

	ix = grep(converted()$species[1,1],gmtFiles)
	if (length(ix) == 0 ) {return(NULL)} else {
		# If selected species is not the default "bestMatch", use that species directly
		if(input$selectOrg != speciesChoice[[1]]) {  
			ix = grep(findSpeciesById(input$selectOrg)[1,1], gmtFiles )
		}
		if(length(ix) == 1)  # if only one file           #WBGene0000001 some ensembl gene ids in lower case
			{x = as.character(gmtFiles[ix]); 	return(x) }else # read in the chosen file 
			return(NULL )   
	}

})	
	
output$downloadPathwayFile <- downloadHandler(
  filename <- function() {
    gsub(".*/","",pathwayFileName())
  },
  content <- function(file) {
    file.copy(pathwayFileName(), file)
  }
)
output$downloadRfunctions <- downloadHandler(
  filename <- function() {
    "iDEP_core_functions.R"
  },
  content <- function(file) {
    file.copy("iDEP_core_functions.R", file)
  }
)
output$downloadRcode <- downloadHandler(
		filename = function() {"iDEP_customized_R_code.R"},
		content = function(file) {
      write( Rcode(), file)	    
	})	

# generate R code based on user input	
Rcode <- reactive({
	if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	
	i = "# R code for stand-alone iDEP analysis"
	i = paste(i,"\n# by Steven Xijin Ge, South Dakota State University,  gexijin@gmail.com ")
	i = paste(i, "\n# Generated by ",iDEPversion, 
		" hosted at http://ge-lab.org/idep/ ", date() )	
	i = paste(i, "\n\n##########################\n# 1. Read data\n##########################" )
	i= paste(i,"\n setwd('C:/Users/Xijin.Ge/Downloads')   # Needs to be changed ")
    i=paste(i, "\n source('iDEP_core_functions.R')")

	i= paste(i,"\n\n # Input files")
	i = paste(i, "\n # Expression file has to use Ensembl for gene ID. Otherwise, use custom pathway database with matching IDs.")
	i = paste(i,"\n inputFile <- 'Downloaded_Converted_Data.csv'  # Expression matrix")  
	if(is.null( sampleInfoFileName() ) ) # if no experiment design file
		i = paste0(i,"\n sampleInfoFile <- NULL ")	 else	
		i = paste0(i,"\n # Experiment design file\n sampleInfoFile <- 'Downloaded_sampleInfoFile.csv'  ")
	i = paste0(i,"\n  #Gene symbols, location etc. \n geneInfoFile <- '",gsub(".*/","",geneInfoFileName()), "' ")
	i = paste0(i,"\n # pathway database in SQL; can be GMT format ")
	i = paste0(i, "\n geneSetFile <- '",gsub(".*/","",pathwayFileName()), "'   ")
	i = paste0(i,"\n STRING10_speciesFile <- 'https://raw.githubusercontent.com/iDEP-SDSU/idep/master/shinyapps/idep/STRING10_species.csv'")

	
	i= paste(i,"\n\n # Parameters")	
	i = op(i,input$missingValue, "Missing values imputation method", TRUE)
	i = op(i, input$dataFileFormat, "1- read counts, 2 FKPM/RPKM or DNA microarray")
	
	if(input$dataFileFormat == 1) { # counts data
		i = op(i, input$minCounts, "Min counts")		
		i = op(i, input$NminSamples, "Minimum number of samples ") 
		i = op(i, input$countsLogStart, "Pseudo count for log CPM")
		i = op(i, input$CountsTransform, "Methods for data transformation of counts. 1-EdgeR's logCPM; 2-VST; 3-rlog")		
	}
	
	if(input$dataFileFormat == 2) {  # FPKM or microrray data
		i = op(i, input$lowFilter, "Min expression score for FPKM")
		i = op(i, input$NminSamples2, "Minimum number of samples ")
		i = op(i, input$logStart, "Pseudo count for FPKM ")
		i = op(i, input$transform, "Log transformation")	
	}
	
	i = paste(i,"\n\n #Read data files\n readData.out <- readData(inputFile)")	
	if(input$dataFileFormat== 2)
		i = paste(i,"\n textTransform() ")
	
	if(!is.null( sampleInfoFileName() ) )
		i = paste(i,"\n readSampleInfo.out <- readSampleInfo(sampleInfoFile)") else 
		i = paste(i,"\n readSampleInfo.out <- NULL") 	
	i = paste(i, "\n input_selectOrg =\"NEW\" ")
	i = op(i, input$selectGO, "Gene set category", TRUE)	# selectGO is for pathway analysis
	i = paste(i, "\n input_noIDConversion = TRUE ")
	
	i = paste(i,"\n allGeneInfo.out <- geneInfo(geneInfoFile)")	
	i = paste(i, 
"\n converted.out = NULL 
 convertedData.out <- convertedData()	 
 nGenesFilter() ")
 
 	if(input$dataFileFormat == 1) { 
	 i = paste(i, "\n convertedCounts.out <- convertedCounts()  # converted counts, just for compatibility") 
	 i = paste(i, "\n readCountsBias()  # detecting bias in sequencing depth")
	 }			
	# 2. Pre-process
	i = paste(i, "\n\n##########################\n# 2. Pre-Process \n##########################" )
	i = paste(i,"\n parDefault = par() \n par(mar=c(12,4,2,2))")	
	if(input$dataFileFormat == 1 ) { # if read counts data
		i = paste(i,"\n # barplot of total read counts		
 x <- readData.out$rawCounts
 groups = as.factor( detectGroups(colnames(x ) ) )
 if(nlevels(groups)<=1 | nlevels(groups) >20 )  
  col1 = 'green'  else
  col1 = rainbow(nlevels(groups))[ groups ]")	
		i = paste(i, "\n barplot( colSums(readData.out$rawCounts)/1e6, 
		col=col1,las=3, main=\"Total read counts (millions)\") ")
	}
	i = paste(i,"\n\n # Box plot")
	i = paste(i,"\n x = readData.out$data")
	i = paste(i, "\n boxplot(x, las = 2, col=col1,
    ylab='Transformed expression levels',
    main='Distribution of transformed data')")

	i = paste(i,"\n\n # Density plot")
	i = paste(i, "\n par(parDefault) \n densityPlot()      ")

	i = paste(i,"\n\n # Scatter plot of the first two samples")
	i = paste(i, "\n plot(x[,1:2],xlab=colnames(x)[1],ylab=colnames(x)[2], 
    main='Scatter plot of first two samples')")
	
	i = paste(i, "\n\n #plot gene or gene family\n input_selectOrg =\"BestMatch\" ")	
	i = op(i, input$geneSearch, "Gene ID for searching", TRUE)	# selectGO is for pathway analysis	 
	i = paste(i, "\n genePlot() ")	 
	i = op(i, input$useSD, "Use standard deviation instead of standard error in error bar?", TRUE)	# selectGO is for pathway analysis	 
	i = paste(i, "\n geneBarPlotError()      ")

	i = paste(i, "\n\n##########################\n# 3. Heatmap \n##########################" )
	i = paste(i,
"\n # hierarchical clustering tree
 x <- readData.out$data
 maxGene <- apply(x,1,max)
 # remove bottom 25% lowly expressed genes, which inflate the PPC
 x <- x[which(maxGene > quantile(maxGene)[1] ) ,] 
 plot(as.dendrogram(hclust2( dist2(t(x)))), ylab=\"1 - Pearson C.C.\", type = \"rectangle\")" )

	i = paste(i,"\n #Correlation matrix")
	i = op(i, input$labelPCC, "Show correlation coefficient?", FALSE)	# selectGO is for pathway analysis	 
	i = paste(i,"\n correlationMatrix()")
	i = paste(i,"\n\n # Parameters for heatmap")
	i = op(i, input$nGenes, "Top genes for heatmap")
	i = op(i, input$geneCentering, "centering genes ?")
	i = op(i, input$sampleCentering, "Center by sample?")
	i = op(i, input$geneNormalize, "Normalize by gene?")
	i = op(i, input$sampleNormalize, "Normalize by sample?")
	i = op(i, input$noSampleClustering, "Use original sample order")
	i = op(i, input$heatmapCutoff, "Remove outliers beyond number of SDs ")
	i = op(i, input$distFunctions, "which distant funciton to use")
	i = op(i, input$hclustFunctions, "Linkage type")
	i = op(i, input$heatColors1, "Colors")
	i = op(i, input$selectFactorsHeatmap, "Sample coloring factors", TRUE)
	i = paste(i, "\n\n staticHeatmap() #Legends not showing due to margin.\n #For a better figure, plot to a file")
	i = paste(i, "\n dev.off() # close plot")
	i = paste(i, "\n tiff('heatmap.tiff', width = 10, height = 15, units = 'in', res = 300, compression = 'lzw')")
	i = paste(i, "\n staticHeatmap() \n dev.off() \n  browseURL('heatmap.tiff') # show heatmap in browser")
	i = paste(i, "\n heatmapPlotly() # interactive heatmap using Plotly")	
 
	i = paste(i, "\n\n##########################\n# 4. k-Means clustering \n##########################" )

	i = op(i, input$nGenesKNN, "Number of genes fro k-Means")	
	i = op(i, input$nClusters, "Number of clusters")
	i = paste(i, "\n maxGeneClustering = 12000")
	i = op(i, input$kmeansNormalization, "Normalization", TRUE)	
	i = op(i, input$KmeansReRun, "Random seed")
	
	i = paste(i, "\n\n distributionSD()  #Distribution of standard deviations")
	i = paste(i, "\n KmeansNclusters()  #Number of clusters")

	i = paste(i, "\n\n Kmeans.out = Kmeans()   #Running K-means")	
	i = paste(i, "\n KmeansHeatmap()   #Heatmap for k-Means \n")
	
	i = paste(i, "\n\n #Read gene sets for enrichment analysis \n sqlite  <- dbDriver('SQLite')")	
	i = op(i, input$selectGO3, "Gene set category", TRUE)	
	i = op(i, input$minSetSize, "Min gene set size")
	i = op(i, input$maxSetSize, "Max gene set size")
	i = paste(i, "\n GeneSets.out <-readGeneSets( geneSetFile,
    convertedData.out, input_selectGO3,input_selectOrg,
    c(input_minSetSize, input_maxSetSize)  ) ")	
	i = paste(i, "\n # Alternatively, users can use their own GMT files by
 #GeneSets.out <- readGMTRobust('somefile.GMT') ")	
	
	i = paste(i, "\n KmeansGO()  #Enrichment analysis for k-Means clusters\n")	
	
	i = op(i, input$seedTSNE, "Random seed for t-SNE")
	i = op(i, input$colorGenes, "Color genes in t-SNE plot?")	
	i = paste(i, "\n tSNEgenePlot()  #Plot genes using t-SNE")		

	i = paste(i, "\n\n##########################\n# 5. PCA and beyond \n##########################" )
	i = op(i, input$selectFactors, "Factor coded by color", TRUE) 
	i = op(i, input$selectFactors2, "Factor coded by shape", TRUE) 

 	i = op(i, input$tsneSeed2, "Random seed for t-SNE") 
	i = paste(i, "\n #PCA, MDS and t-SNE plots\n PCAplot()	\n MDSplot()\n tSNEplot() ")

	i = paste(i, "\n\n #Read gene sets for pathway analysis using PGSEA on principal components")
	if(is.null(input$selectGO6 ) ) 
		i = paste(i, "\n input_selectGO6 <- 'GOBP'") else 
		i = op(i, input$selectGO6, "Gene set category", TRUE) 
	i = paste(i, "\n GeneSets.out <-readGeneSets( geneSetFile,
    convertedData.out, input_selectGO6,input_selectOrg,
    c(input_minSetSize, input_maxSetSize)  ) ")

	i = paste(i,"\n PCApathway() # Run PGSEA analysis" );
	i = paste(i,"\n cat( PCA2factor() )   #The correlation between PCs with factors");

	i = paste(i, "\n\n##########################\n# 6. DEG1 \n##########################" )
	i = op(i, input$CountsDEGMethod, "DESeq2= 3,limma-voom=2,limma-trend=1 ")
	i = op(i, input$limmaPval, "FDR cutoff")
	i = op(i, input$limmaFC, "Fold-change cutoff")
	i = op(i, input$selectModelComprions, "Selected comparisons", TRUE)
	i = op(i, input$selectFactorsModel, "Selected comparisons", TRUE)
	
	i = op(i, input$selectInteractions, "Selected comparisons", TRUE)
	i = op(i, input$selectBlockFactorsModel, "Selected comparisons", TRUE)
	if( is.null(factorReferenceLevels() ) )
		i = paste(i, "\n factorReferenceLevels.out <- NULL") else
		i = paste0(i, "\n factorReferenceLevels.out <- c('", paste0(factorReferenceLevels(), collapse="','" ),"')")
	i = paste(i,"\n\n limma.out <- limma()
 limma.out$comparisons
 DEG.data.out <- DEG.data()")
	if(is.null(input$selectComparisonsVenn) ) 
		i = paste(i, "\n input_selectComparisonsVenn = limma.out$comparisons[1:3] # use first three comparisons") else
		i = op(i, input$selectComparisonsVenn, "Selected comparisons for Venn diagram", TRUE)
 	i = op(i, input$UpDownRegulated, "Split up and down regulated genes") 
	i = paste(i, "\n vennPlot() # Venn diagram
 sigGeneStats() # number of DEGs as figure
 sigGeneStatsTable() # number of DEGs as table")

 
	i = paste(i, "\n\n##########################\n# 7. DEG2 \n##########################" ) 
	if(is.null(input$selectContrast) ) 
		i = paste(i, "\n input_selectComparisonsVenn = limma.out$comparisons[1] # use first  comparisons") else
		i = op(i, input$selectContrast, "Selected comparisons", TRUE)
	i = paste(i, "\n selectedHeatmap.data.out <- selectedHeatmap.data()
 selectedHeatmap()   # heatmap for DEGs in selected comparison
\n # Save gene lists and data into files
 write.csv( selectedHeatmap.data()$genes, 'heatmap.data.csv') 
 write.csv(DEG.data(),'DEG.data.csv' )
 write(AllGeneListsGMT() ,'AllGeneListsGMT.gmt')\n")
	
	if( is.null( input$selectGO2 ))
		i = paste(i,"\n input_selectGO2 = 'GOBP' ") else
		i = op(i, input$selectGO2, "Gene set category",TRUE) 
	i = paste(i,"\n geneListData.out <- geneListData()
 volcanoPlot() 
 scatterPlot()
 MAplot() 
 geneListGOTable.out <- geneListGOTable() ")
	i = paste(i, "\n # Read pathway data again \n GeneSets.out <-readGeneSets( geneSetFile,
    convertedData.out, input_selectGO2,input_selectOrg,
    c(input_minSetSize, input_maxSetSize)  ) ")	
 	i = op(i, input$removeRedudantSets, "Remove highly redundant gene sets?") 
	i = paste(i,"\n geneListGO()");
	
#enrichmentPlot(geneListGOTable.out, rightMargin=25  )
#enrichmentNetwork(geneListGOTable.out )	
	i = paste(i,"\n\n # STRING-db API access" )
	i = paste(i,"\n STRING10_species = read.csv(STRING10_speciesFile) ")
	i = paste(i,"\n ix = grep('Mus musculus', STRING10_species$official_name )
 findTaxonomyID.out <- STRING10_species[ix,1] # find taxonomyID
 findTaxonomyID.out  
 # users can also skip the above and assign NCBI taxonomy id directly by
 # findTaxonomyID.out = 10090 # mouse 10090, human 9606 etc.
 STRINGdb_geneList.out <- STRINGdb_geneList() #convert gene lists")
	
 	i = op(i, input$STRINGdbGO, "'Process', 'Component', 'Function', 'KEGG', 'Pfam', 'InterPro'", TRUE) 
	i = paste(i,"\n stringDB_GO_enrichmentData()")
	
	i = paste(i,"\n\n # PPI network retrieval and analysis")	
 	i = op(i, input$nGenesPPI, "Number of top genes for PPI retrieval and analysis") 
	i = paste(i,"\n stringDB_network1(1) #Show PPI network")
	i = paste(i,"\n write(stringDB_network_link(), 'PPI_results.html') # write results to html file")
	i = paste(i,"\n browseURL('PPI_results.html') # open in browser" )
	

	i = paste(i, "\n\n##########################\n# 8. Pathway analysis \n##########################" )  
   	i = op(i, input$selectContrast1, "select Comparison", TRUE) 
	i = paste(i,"\n #input_selectContrast1 = limma.out$comparisons[3] # manually set")
	i = op(i, input$selectGO, "Gene set category", TRUE) 
	i = paste(i,"\n #input_selectGO='custom' # if custom gmt file" ) 
  	i = op(i, input$minSetSize, "Min size for gene set")  
  	i = op(i, input$maxSetSize, "Max size for gene set")  
	i = paste(i, "\n # Read pathway data again \n GeneSets.out <-readGeneSets( geneSetFile,
    convertedData.out, input_selectGO,input_selectOrg,
    c(input_minSetSize, input_maxSetSize)  ) ")	
  	i = op(i, input$pathwayPvalCutoff, "FDR cutoff") 
  	i = op(i, input$nPathwayShow, "Top pathways to show") 	
  	i = op(i, input$absoluteFold, "Use absolute values of fold-change?") 	
  	i = op(i, input$GenePvalCutoff, "FDR to remove genes") 

	i = paste(i,"\n\n input_pathwayMethod = 1  # 1  GAGE
 gagePathwayData.out <- gagePathwayData()  # pathway analysis using GAGE  
 gagePathwayData.out
  pathwayListData.out = pathwayListData()
 enrichmentPlot(pathwayListData.out, 25  )
 enrichmentNetwork(pathwayListData.out )
 enrichmentNetworkPlotly(pathwayListData.out)

 input_pathwayMethod = 3  # 1  fgsea 
 fgseaPathwayData.out <- fgseaPathwayData() #Pathway analysis using fgsea
 fgseaPathwayData.out
 pathwayListData.out = pathwayListData()
 enrichmentPlot(pathwayListData.out, 25  )
 enrichmentNetwork(pathwayListData.out )
 enrichmentNetworkPlotly(pathwayListData.out) 
  
 PGSEAplot() # pathway analysis using PGSEA" )
	
# PGSEA for selected contrast not showing up

	i = paste(i, "\n\n##########################\n# 9. Chromosome \n##########################" )  
   	i = op(i, input$selectContrast2, "select Comparison", TRUE) 
	i = paste(i,"\n #input_selectContrast2 = limma.out$comparisons[3] # manually set")
  	i = op(i, input$limmaPvalViz, "FDR to filter genes") 
  	i = op(i, input$limmaFCViz, "FDR to filter genes") 
	i = paste(i, "\n genomePlotly() # shows fold-changes on the genome")
	
	
	i = paste(i, "\n\n##########################\n# 10. Bicluster \n##########################" )  
  	i = op(i, input$nGenesBiclust, "Top genes for biclustering") 
  	i = op(i, input$biclustMethod, "Method: 'BCCC', 'QUBIC', 'runibic' ...", TRUE) 
	i = paste(i, "\n biclustering.out = biclustering()  # run analysis\n")
  	i = op(i, input$selectBicluster, "select a cluster") 	
	i = paste(i, "\n biclustHeatmap()   # heatmap for selected cluster" )
  	i = op(i, input$selectGO4, "Gene set", TRUE)
	i = paste(i, "\n # Read pathway data again \n GeneSets.out <-readGeneSets( geneSetFile,
    convertedData.out, input_selectGO4,input_selectOrg,
    c(input_minSetSize, input_maxSetSize)  ) ")		
	i = paste(i, "\n geneListBclustGO()  # enrichment analysis")
	
	i = paste(i, "\n\n##########################\n# 11. Co-expression network \n##########################" )  
  	i = op(i, input$mySoftPower, "SoftPower to cutoff") 
  	i = op(i, input$nGenesNetwork, "Number of top genes") 
  	i = op(i, input$minModuleSize, "Module size minimum") 
	i = paste(i, "\n wgcna.out = wgcna()   # run WGCNA
 softPower()  # soft power curve
 modulePlot()  # plot modules
 listWGCNA.Modules.out = listWGCNA.Modules() #modules\n")
  	i = op(i, input$selectGO5, "Gene set", TRUE) 
	i = paste(i, "\n # Read pathway data again \n GeneSets.out <-readGeneSets( geneSetFile,
    convertedData.out, input_selectGO5,input_selectOrg,
    c(input_minSetSize, input_maxSetSize)  ) ")	
  	i = op(i, input$selectWGCNA.Module, "Select a module", TRUE) 
  	i = op(i, input$topGenesNetwork, "SoftPower to cutoff") 
  	i = op(i, input$edgeThreshold, "Number of top genes") 
	i = paste(i, "\n moduleNetwork()	# show network of top genes in selected module\n")
	

  	i = op(i, input$removeRedudantSets, "Remove redundant gene sets")
	i = paste(i, "\n networkModuleGO()	# Enrichment analysis of selected module")
 
	return(i)
 })

output$downloadRcodeMarkdown <- downloadHandler(
		filename = function() {"iDEP_R_Markdown.Rmd"},
		content = function(file) {
      write( RcodeMarkdown(), file)	    
	})
	
# generate R code based on user input	
RcodeMarkdown <- reactive({
	if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
	

	i = "---
title: 'Reproducing iDEP analyses with auto-generated R Markdown'
author: "
	i = paste(i,iDEPversion, " http://ge-lab.org/idep/, originally by Steven Xijin.Ge@sdstate.edu ")
	i = paste(i,"\ndate:", date() )
	i = paste(i," 
output: html_document
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width=6, fig.height=5, fig.align = 'center') 
```")
		
	i = paste(i, "\n This R markdown file was auto-generated by the [iDEP website](http://ge-lab.org/idep/). It is assumed that users have analyzed their data with iDEP by clicking through all the tabs and have downloaded the related files to a folder. 
	\n\n## 1. Read data ")	
		
#	i = paste(i, "\n\n##########################\n# 1. Read data\n##########################" )
	i = paste(i,"\nFirst we set up the working directory to where the files are saved.   ")
	i= paste(i,"\n```{r, message=FALSE} \n setwd('C:/Users/Xijin.Ge/Downloads')   # Needs to be changed ")
	
  	i = paste(i,"\n```\nR packages and iDEP core Functions. 
	Users can also download the iDEP_core_functions.R file. 
	Many R packages needs to be installed first. This may take hours. 
	Each of these packages took years to develop.So be a patient thief. Sometimes dependencies needs to be installed manually. 
	If you are using an older version of R, and having trouble with package installation,
	try un-install the current version of R, delete all folders and files 
	(C:/Program Files/R/R-3.4.3), and 
	reinstall from scratch.  \n```{r, message=FALSE  } ") # R Markdown 

	i=paste(i, "\n if(file.exists('iDEP_core_functions.R'))
	source('iDEP_core_functions.R') else 
    source('https://raw.githubusercontent.com/iDEP-SDSU/idep/master/shinyapps/idep/iDEP_core_functions.R')")
 	i = paste(i,"\n```  \nWe are using the downloaded gene expression file where gene IDs has 
	been converted to Ensembl gene IDs. This is because the ID conversion database is too large
	to download. You can use your original file if your file uses Ensembl ID, or you do not want 
	to use the pathway files available in iDEP (or it is not available).  \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	
#	i= paste(i,"\n\n # Input files")
#	i = paste(i, "\n ### Expression file has to use Ensembl for gene ID. Otherwise, use custom pathway database with matching IDs.")
	i = paste(i,"\n inputFile <- 'Downloaded_Converted_Data.csv'  # Expression matrix")  
	if(is.null( sampleInfoFileName() ) ) # if no experiment design file
		i = paste0(i,"\n sampleInfoFile <- NULL ")	 else	
		i = paste0(i,"\n sampleInfoFile <- 'Downloaded_sampleInfoFile.csv' # Experiment design file ")
	i = paste0(i,"\n geneInfoFile <- '",gsub(".*/","",geneInfoFileName()), "' #Gene symbols, location etc. ")
	i = paste0(i, "\n geneSetFile <- '",gsub(".*/","",pathwayFileName()), "'  # pathway database in SQL; can be GMT format ")
	i = paste0(i,"\n STRING10_speciesFile <- 'https://raw.githubusercontent.com/iDEP-SDSU/idep/master/shinyapps/idep/STRING10_species.csv'")
 	i = paste(i,"\n``` \nParameters for reading data  \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	

	i = op(i,input$missingValue, "Missing values imputation method", TRUE)
	i = op(i, input$dataFileFormat, "1- read counts, 2 FKPM/RPKM or DNA microarray")
	
	if(input$dataFileFormat == 1) { # counts data
		i = op(i, input$minCounts, "Min counts")		
		i = op(i, input$NminSamples, "Minimum number of samples ") 
		i = op(i, input$countsLogStart, "Pseudo count for log CPM")
		i = op(i, input$CountsTransform, "Methods for data transformation of counts. 1-EdgeR's logCPM 2-VST, 3-rlog")		
	}
	
	if(input$dataFileFormat == 2) {  # FPKM or microrray data
		i = op(i, input$lowFilter, "Min expression score for FPKM")
		i = op(i, input$NminSamples2, "Minimum number of samples ")
		i = op(i, input$logStart, "Pseudo count for FPKM ")
		i = op(i, input$transform, "Log transformation")	
	}
	
 	i = paste(i,"\n```\n```{r, message=FALSE  } ") # R Markdown block-------------------------------------			
	i = paste(i,"\n\n readData.out <- readData(inputFile)")	
	i = paste(i,"\n library(knitr)   #  install if needed. for showing tables with kable")	
	i = paste(i,"\n kable( head(readData.out$data) )    # show the first few rows of data")
	if(input$dataFileFormat== 2)
		i = paste(i,"\n textTransform() ")	
	i = paste(i,"\n```\n```{r, message=FALSE  } ") # R Markdown block

	if(!is.null( sampleInfoFileName() ) )
		i = paste(i,"\n readSampleInfo.out <- readSampleInfo(sampleInfoFile) \n kable( readSampleInfo.out ) ") else 
		i = paste(i,"\n readSampleInfo.out <- NULL") 	
 	i = paste(i,"\n```\n```{r, message=FALSE  } ") # R Markdown block-------------------------------------			
	i = paste(i, "\n input_selectOrg =\"NEW\" ")
	i = op(i, input$selectGO, "Gene set category", TRUE)	# selectGO is for pathway analysis
	i = paste(i, "\n input_noIDConversion = TRUE ")
	
	i = paste(i,"\n allGeneInfo.out <- geneInfo(geneInfoFile)")	
	i = paste(i, 
"\n converted.out = NULL 
 convertedData.out <- convertedData()	 
 nGenesFilter() ")

 	if(input$dataFileFormat == 1) { 
	 i = paste(i, "\n convertedCounts.out <- convertedCounts()  # converted counts, just for compatibility") 

	 }
 	i = paste(i,"\n```\n\n## 2. Pre-process \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------			
	# 2. Pre-process
#	i = paste(i, "\n\n##########################\n# 2. Pre-Process \n##########################" )
	i = paste(i,"\n# Read counts per library \n parDefault = par() \n par(mar=c(12,4,2,2))")	
	if(input$dataFileFormat == 1 ) { # if read counts data
		i = paste(i,"\n # barplot of total read counts
 x <- readData.out$rawCounts
 groups = as.factor( detectGroups(colnames(x ) ) )
 if(nlevels(groups)<=1 | nlevels(groups) >20 )  
  col1 = 'green'  else
  col1 = rainbow(nlevels(groups))[ groups ]				
		")
		i = paste(i, "\n barplot( colSums(x)/1e6, 
		col=col1,las=3, main=\"Total read counts (millions)\") ")
		
	 i = paste(i, "\n readCountsBias()  # detecting bias in sequencing depth")		
		
	}
 	i = paste(i,"\n```\n```{r, message=FALSE  } ") # R Markdown block-------------------------------------		
	i = paste(i,"\n\n # Box plot")
	i = paste(i,"\n x = readData.out$data")
	i = paste(i, "\n boxplot(x, las = 2, col=col1,
    ylab='Transformed expression levels',
    main='Distribution of transformed data')")
 	i = paste(i,"\n```\n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	i = paste(i,"\n\n #Density plot")
	i = paste(i, "\n par(parDefault) \n densityPlot()      ")
 	i = paste(i,"\n```\n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	i = paste(i,"\n\n # Scatter plot of the first two samples")
	i = paste(i, "\n plot(x[,1:2],xlab=colnames(x)[1],ylab=colnames(x)[2], 
    main='Scatter plot of first two samples')")
 	i = paste(i,"\n```\n```{r, message=FALSE  } ") # R Markdown block-------------------------------------		
	i = paste(i, "\n\n ####plot gene or gene family\n input_selectOrg =\"BestMatch\" ")	
	i = op(i, input$geneSearch, "Gene ID for searching", TRUE)	# selectGO is for pathway analysis	 
	i = paste(i, "\n genePlot() ")	 
 	i = paste(i,"\n```\n```{r, message=FALSE  } ") # R Markdown block-------------------------------------		
	i = op(i, input$useSD, "Use standard deviation instead of standard error in error bar?", TRUE)	# selectGO is for pathway analysis	 
	i = paste(i, "\n geneBarPlotError()      ")
 	i = paste(i,"\n```\n\n## 3. Heatmap  \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	
#	i = paste(i, "\n\n##########################\n# 3. Heatmap \n##########################" )
	i = paste(i,
"\n # hierarchical clustering tree
 x <- readData.out$data
 maxGene <- apply(x,1,max)
 # remove bottom 25% lowly expressed genes, which inflate the PPC
 x <- x[which(maxGene > quantile(maxGene)[1] ) ,] 
 plot(as.dendrogram(hclust2( dist2(t(x)))), ylab=\"1 - Pearson C.C.\", type = \"rectangle\")" )
	i = paste(i,"\n```\n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	

	i = paste(i,"\n #Correlation matrix")
	i = op(i, input$labelPCC, "Show correlation coefficient?", FALSE)	# selectGO is for pathway analysis	 
	i = paste(i,"\n correlationMatrix()")
	i = paste(i,"\n```\n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	

	i = paste(i,"\n\n # Parameters for heatmap")
	i = op(i, input$nGenes, "Top genes for heatmap")
	i = op(i, input$geneCentering, "centering genes ?")
	i = op(i, input$sampleCentering, "Center by sample?")
	i = op(i, input$geneNormalize, "Normalize by gene?")
	i = op(i, input$sampleNormalize, "Normalize by sample?")
	i = op(i, input$noSampleClustering, "Use original sample order")
	i = op(i, input$heatmapCutoff, "Remove outliers beyond number of SDs ")
	i = op(i, input$distFunctions, "which distant funciton to use")
	i = op(i, input$hclustFunctions, "Linkage type")
	i = op(i, input$heatColors1, "Colors")
	i = op(i, input$selectFactorsHeatmap, "Sample coloring factors", TRUE)
	i = paste(i, "\n png('heatmap.png', width = 10, height = 15, units = 'in', res = 300)")
	i = paste(i, "\n staticHeatmap() \n dev.off() ")
	i = paste(i,"\n```\n  ![heatmap] (heatmap.png)   \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	

	i = paste(i, "\n heatmapPlotly() # interactive heatmap using Plotly")	
 	i = paste(i,"\n```\n\n## 4. K-means clustering   \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	

#	i = paste(i, "\n\n##########################\n# 4. k-Means clustering \n##########################" )

	i = op(i, input$nGenesKNN, "Number of genes fro k-Means")	
	i = op(i, input$nClusters, "Number of clusters")
	i = paste(i, "\n maxGeneClustering = 12000")
	i = op(i, input$kmeansNormalization, "Normalization", TRUE)	
	i = op(i, input$KmeansReRun, "Random seed")
	i = paste(i, "\n\n distributionSD()  #Distribution of standard deviations")
 	i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------		
	i = paste(i, "\n KmeansNclusters()  #Number of clusters")
 	i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	i = paste(i, "\n\n Kmeans.out = Kmeans()   #Running K-means")	
	i = paste(i, "\n KmeansHeatmap()   #Heatmap for k-Means \n")
 	i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------		
	i = paste(i, "\n\n #Read gene sets for enrichment analysis \n sqlite  <- dbDriver('SQLite')")	
	i = op(i, input$selectGO3, "Gene set category", TRUE)	
	i = op(i, input$minSetSize, "Min gene set size")
	i = op(i, input$maxSetSize, "Max gene set size")
	i = paste(i, "\n GeneSets.out <-readGeneSets( geneSetFile,
    convertedData.out, input_selectGO3,input_selectOrg,
    c(input_minSetSize, input_maxSetSize)  ) ")	
	i = paste(i, "\n # Alternatively, users can use their own GMT files by
 #GeneSets.out <- readGMTRobust('somefile.GMT') ")	
	
	i = paste(i, "\n results <- KmeansGO()  #Enrichment analysis for k-Means clusters	
 results$adj.Pval <- format( results$adj.Pval,digits=3 )
 kable( results, row.names=FALSE)")
 	i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	
	i = op(i, input$seedTSNE, "Random seed for t-SNE")
	i = op(i, input$colorGenes, "Color genes in t-SNE plot?")	
	i = paste(i, "\n tSNEgenePlot()  #Plot genes using t-SNE")		
 	i = paste(i,"\n```\n\n## 5. PCA and beyond   \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	
#	i = paste(i, "\n\n##########################\n# 5. PCA and beyond \n##########################" )
	if(is.null(input$selectFactors ) )
		i = paste(i,"\n input_selectFactors <- 'Sample_Name' ") else 
		i = op(i, input$selectFactors, "Factor coded by color", TRUE) 
	if(is.null(input$selectFactors2 ) )
		i = paste(i,"\n input_selectFactors2 <- 'Sample_Name' ") else 		
		i = op(i, input$selectFactors2, "Factor coded by shape", TRUE) 

 	i = op(i, input$tsneSeed2, "Random seed for t-SNE") 
	i = paste(i, "\n #PCA, MDS and t-SNE plots\n PCAplot() ")
 	i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------		
	i = paste(i, "\n MDSplot()")
 	i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------		
	i = paste(i,"\n tSNEplot() ")
 	i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	i = paste(i, "\n\n #Read gene sets for pathway analysis using PGSEA on principal components")
	if(is.null(input$selectGO6 ) ) 
		i = paste(i, "\n input_selectGO6 <- 'GOBP'") else 
		i = op(i, input$selectGO6, "Gene set category", TRUE) 
	i = paste(i, "\n GeneSets.out <-readGeneSets( geneSetFile,
    convertedData.out, input_selectGO6,input_selectOrg,
    c(input_minSetSize, input_maxSetSize)  ) ")

	i = paste(i,"\n PCApathway() # Run PGSEA analysis" );
 	i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------		
	i = paste(i,"\n cat( PCA2factor() )   #The correlation between PCs with factors");
 	i = paste(i,"\n```\n\n## 6. DEG1   \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	
#	i = paste(i, "\n\n##########################\n# 6. DEG1 \n##########################" )
	i = op(i, input$CountsDEGMethod, "DESeq2= 3,limma-voom=2,limma-trend=1 ")
	i = op(i, input$limmaPval, "FDR cutoff")
	i = op(i, input$limmaFC, "Fold-change cutoff")
	i = op(i, input$selectModelComprions, "Selected comparisons", TRUE)
	i = op(i, input$selectFactorsModel, "Selected comparisons", TRUE)
	
	i = op(i, input$selectInteractions, "Selected comparisons", TRUE)
	i = op(i, input$selectBlockFactorsModel, "Selected comparisons", TRUE)
	if( is.null(factorReferenceLevels() ) )
		i = paste(i, "\n factorReferenceLevels.out <- NULL") else
		i = paste0(i, "\n factorReferenceLevels.out <- c('", paste0(factorReferenceLevels(), collapse="','" ),"')")
	i = paste(i,"\n\n limma.out <- limma()
 DEG.data.out <- DEG.data()
 limma.out$comparisons")
 
  	i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------		
	if(is.null(input$selectComparisonsVenn) ) 
		i = paste(i, "\n input_selectComparisonsVenn = limma.out$comparisons[1:3] # use first three comparisons") else
		i = op(i, input$selectComparisonsVenn, "Selected comparisons for Venn diagram", TRUE)
 	i = op(i, input$UpDownRegulated, "Split up and down regulated genes") 
	i = paste(i, "\n vennPlot() # Venn diagram")
  	i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------		
	i = paste(i, "\n  sigGeneStats() # number of DEGs as figure")
  	i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------			
	i = paste(i, "\n  sigGeneStatsTable() # number of DEGs as table")
 	i = paste(i,"\n```\n\n## 7. DEG2   \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
 
#	i = paste(i, "\n\n##########################\n# 7. DEG2 \n##########################" ) 
	if(is.null(input$selectContrast) ) 
		i = paste(i, "\n input_selectContrast = limma.out$comparisons[1] # use first  comparisons") else
		i = op(i, input$selectContrast, "Selected comparisons", TRUE)
	i = paste(i, "\n selectedHeatmap.data.out <- selectedHeatmap.data()
 selectedHeatmap()   # heatmap for DEGs in selected comparison
\n # Save gene lists and data into files
 write.csv( selectedHeatmap.data()$genes, 'heatmap.data.csv') 
 write.csv(DEG.data(),'DEG.data.csv' )
 write(AllGeneListsGMT() ,'AllGeneListsGMT.gmt')\n")
  	i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	
	if( is.null( input$selectGO2 ))
		i = paste(i,"\n input_selectGO2 = 'GOBP'  # gene set category") else
		i = op(i, input$selectGO2, "Gene set category",TRUE) 
	 	
	i = paste(i,"\n geneListData.out <- geneListData() ")
	i = paste(i,"\n volcanoPlot() ")
  	i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	 
	i = paste(i,"\n  scatterPlot() ")
  	i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	i = paste(i,"\n  MAplot() ")
  	i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------		
	i = paste(i,"\n  geneListGOTable.out <- geneListGOTable() ")
	i = paste(i, "\n # Read pathway data again \n GeneSets.out <-readGeneSets( geneSetFile,
    convertedData.out, input_selectGO2,input_selectOrg,
    c(input_minSetSize, input_maxSetSize)  ) ")	
 	i = op(i, input$removeRedudantSets, "Remove highly redundant gene sets?") 

	i = paste(i, "\n results <- geneListGO()  #Enrichment analysis
 results$adj.Pval <- format( results$adj.Pval,digits=3 )
 kable( results, row.names=FALSE)")	
  	i = paste(i,"\n```\n\nSTRING-db API access. 
	We need to find the taxonomy id of your species, this used by STRING.
  First we try to guess the ID based on iDEP's database. Users can also skip this step and assign NCBI taxonomy id directly by
  findTaxonomyID.out = 10090 # mouse 10090, human 9606 etc.
	\n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	
#enrichmentPlot(geneListGOTable.out, rightMargin=25  )
#enrichmentNetwork(geneListGOTable.out )	
	i = paste(i,"\n\n " )
	i = paste(i,"\n STRING10_species = read.csv(STRING10_speciesFile) ")
	i = paste(i,"\n ix = grep('Mus musculus', STRING10_species$official_name )
 findTaxonomyID.out <- STRING10_species[ix,1] # find taxonomyID
 findTaxonomyID.out  
")
   	i = paste(i,"\n``` \nEnrichment analysis using STRING     \n```{r, message=FALSE  } ") 
	i = paste(i,"\n  STRINGdb_geneList.out <- STRINGdb_geneList() #convert gene lists")	
 	i = op(i, input$STRINGdbGO, "'Process', 'Component', 'Function', 'KEGG', 'Pfam', 'InterPro'", TRUE) 
	i = paste(i, "\n results <- stringDB_GO_enrichmentData()  # enrichment using STRING	
 results$adj.Pval <- format( results$adj.Pval,digits=3 )
 kable( results, row.names=FALSE)")	
	
   	i = paste(i,"\n``` \nPPI network retrieval and analysis    \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------			
 	i = op(i, input$nGenesPPI, "Number of top genes for PPI retrieval and analysis") 
	i = paste(i,"\n stringDB_network1(1) #Show PPI network")
   	i = paste(i,"\n``` \nGenerating interactive PPI   \n```{r, message=FALSE  } ") # R Markdown 	
	i = paste(i,"\n write(stringDB_network_link(), 'PPI_results.html') # write results to html file")
	i = paste(i,"\n browseURL('PPI_results.html') # open in browser" )
   	i = paste(i,"\n```\n\n## 8. Pathway analysis   \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------			

#	i = paste(i, "\n\n##########################\n# 8. Pathway analysis \n##########################" ) 
	if( is.null( input$selectContrast1 ))
		i = paste(i,"\n input_selectContrast1 = limma.out$comparisons[1]") else
		i = op(i, input$selectContrast1, "select Comparison", TRUE) 
	i = paste(i,"\n #input_selectContrast1 = limma.out$comparisons[3] # manually set")

	if( is.null( input$selectGO ))
		i = paste(i,"\n input_selectGO = 'GOBP'  # gene set category") else
		i = op(i, input$selectGO, "Gene set category",TRUE) 	
	i = paste(i,"\n #input_selectGO='custom' # if custom gmt file" ) 
  	i = op(i, input$minSetSize, "Min size for gene set")  
  	i = op(i, input$maxSetSize, "Max size for gene set")  
	i = paste(i, "\n # Read pathway data again \n GeneSets.out <-readGeneSets( geneSetFile,
    convertedData.out, input_selectGO,input_selectOrg,
    c(input_minSetSize, input_maxSetSize)  ) ")	
  	i = op(i, input$pathwayPvalCutoff, "FDR cutoff") 
  	i = op(i, input$nPathwayShow, "Top pathways to show") 	
  	i = op(i, input$absoluteFold, "Use absolute values of fold-change?") 	
  	i = op(i, input$GenePvalCutoff, "FDR to remove genes") 

	i = paste(i,"\n\n input_pathwayMethod = 1  # 1  GAGE
 gagePathwayData.out <- gagePathwayData()  # pathway analysis using GAGE  
  ")
 	i = paste(i, "\n results <- gagePathwayData.out  #Enrichment analysis for k-Means clusters	
 results$adj.Pval <- format( results$adj.Pval,digits=3 )
 kable( results, row.names=FALSE)")	
 
    i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------		
	i = paste(i,"\n pathwayListData.out = pathwayListData() 
 enrichmentPlot(pathwayListData.out, 25  )")
    i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	 
 	i = paste(i,"\n  enrichmentNetwork(pathwayListData.out ) ")
    i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
 	i = paste(i,"\n  enrichmentNetworkPlotly(pathwayListData.out)")
    i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	
	i = paste(i,"\n\n input_pathwayMethod = 3  # 1  fgsea 
 fgseaPathwayData.out <- fgseaPathwayData() #Pathway analysis using fgsea")
  	i = paste(i, "\n results <- fgseaPathwayData.out  #Enrichment analysis for k-Means clusters	
 results$adj.Pval <- format( results$adj.Pval,digits=3 )
 kable( results, row.names=FALSE)")
    i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------		
	i = paste(i,"\n  pathwayListData.out = pathwayListData() 
 enrichmentPlot(pathwayListData.out, 25  )")
    i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	 
 	i = paste(i,"\n  enrichmentNetwork(pathwayListData.out ) ")
    i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
 	i = paste(i,"\n  enrichmentNetworkPlotly(pathwayListData.out)")
    i = paste(i,"\n```     \n```{r, message=FALSE,fig.width=9, fig.height=8  } ") # R Markdown block-------------------------------------	
	

 	i = paste(i,"\n   PGSEAplot() # pathway analysis using PGSEA" )
	
# PGSEA for selected contrast not showing up
    i = paste(i,"\n```\n\n## 9. Chromosome   \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
#	i = paste(i, "\n\n##########################\n# 9. Chromosome \n##########################" )  
	if( is.null( input$selectContrast2 ))
		i = paste(i,"\n input_selectContrast2 = limma.out$comparisons[1]") else
		i = op(i, input$selectContrast2, "select Comparison", TRUE) 

	i = paste(i,"\n #input_selectContrast2 = limma.out$comparisons[3] # manually set")
  	i = op(i, input$limmaPvalViz, "FDR to filter genes") 
  	i = op(i, input$limmaFCViz, "FDR to filter genes") 
	i = paste(i, "\n genomePlotly() # shows fold-changes on the genome")
	
    i = paste(i,"\n```\n\n## 10. Biclustering   \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------		
#	i = paste(i, "\n\n##########################\n# 10. Bicluster \n##########################" )  
  	i = op(i, input$nGenesBiclust, "Top genes for biclustering") 
  	i = op(i, input$biclustMethod, "Method: 'BCCC', 'QUBIC', 'runibic' ...", TRUE) 
	i = paste(i, "\n biclustering.out = biclustering()  # run analysis\n")
  	i = op(i, input$selectBicluster, "select a cluster") 	
	i = paste(i, "\n biclustHeatmap()   # heatmap for selected cluster" )
    i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	
	if( is.null( input$selectGO4 ))
		i = paste(i,"\n input_selectGO4 = 'GOBP'  # gene set category") else
		i = op(i, input$selectGO4, "Gene set category",TRUE) 	
	i = paste(i, "\n # Read pathway data again \n GeneSets.out <-readGeneSets( geneSetFile,
    convertedData.out, input_selectGO4,input_selectOrg,
    c(input_minSetSize, input_maxSetSize)  ) ")		

	i = paste(i, "\n results <- geneListBclustGO()  #Enrichment analysis for k-Means clusters	
 results$adj.Pval <- format( results$adj.Pval,digits=3 )
 kable( results, row.names=FALSE)")	
	
	
	
    i = paste(i,"\n```\n\n## 11. Co-expression network    \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	
#	i = paste(i, "\n\n##########################\n# 11. Co-expression network \n##########################" )  
  	i = op(i, input$mySoftPower, "SoftPower to cutoff") 
  	i = op(i, input$nGenesNetwork, "Number of top genes") 
  	i = op(i, input$minModuleSize, "Module size minimum") 
	i = paste(i, "\n wgcna.out = wgcna()   # run WGCNA ")
	    i = paste(i,"\n```     \n```{r, message=FALSE  } ")
	i = paste(i,"\n softPower()  # soft power curve")
    i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	i = paste(i, "\n  modulePlot()  # plot modules ")
	i = paste(i, "\n  listWGCNA.Modules.out = listWGCNA.Modules() #modules\n")
    i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------	
	if( is.null( input$selectGO5 ))
		i = paste(i,"\n input_selectGO5 = 'GOBP'  # gene set category") else
		i = op(i, input$selectGO5, "Gene set category",TRUE) 
	i = paste(i, "\n # Read pathway data again \n GeneSets.out <-readGeneSets( geneSetFile,
    convertedData.out, input_selectGO5,input_selectOrg,
    c(input_minSetSize, input_maxSetSize)  ) ")	
  	i = op(i, input$selectWGCNA.Module, "Select a module", TRUE) 
  	i = op(i, input$topGenesNetwork, "SoftPower to cutoff") 
  	i = op(i, input$edgeThreshold, "Number of top genes") 
	i = paste(i, "\n moduleNetwork()	# show network of top genes in selected module\n")
    i = paste(i,"\n```     \n```{r, message=FALSE  } ") # R Markdown block-------------------------------------		

  	i = op(i, input$removeRedudantSets, "Remove redundant gene sets")

	i = paste(i, "\n results <- networkModuleGO()  #Enrichment analysis of selected module
 results$adj.Pval <- format( results$adj.Pval,digits=3 )
 kable( results, row.names=FALSE)")		
	
	
    i = paste(i,"\n```   ") # R Markdown block-------------------------------------	
	return(i)
 })
 
 
})  # shiny Server
