####################################################
# Author: Steven Ge Xijin.Ge@sdstate.edu
# Co-author: Eric Tulowetzke, eric.tulowetzke@jacks.sdstate.edu
# Lab: Ge Lab
# R version 4.0.5
# Project: iDEP v93
# File: utility_functions
# Purpose of file:hold all of utility_functions for iDEP project 
# Start data: 06-06-2022 (mm-dd-yyyy)
# Data last modified: 06-06-2021, 12:46 PM CST (mm-dd-yyyy,TIME) 
# to help with github merge 
#######################################################

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