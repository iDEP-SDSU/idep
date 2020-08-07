library('R6')
source('server.config')

UtilFuns <- R6Class("UtilFuns")

UtilFuns$set("public", "DistanceFuns", NULL)
UtilFuns$set("public", "HierarchicalClusteringFuns", NULL)


UtilFuns$set("public", "initialize", 
	function(){
		self$DistanceFuns <- list(
			Correlation=self$dist2, 
			Euclidean=dist, 
			AbsolutePCC=self$dist3)

		self$HierarchicalClusteringFuns <- list(
			average  = self$hclust2, 
            complete = hclust, 
            single   = self$hclust.single,
			median   = self$hclust.median, 
            centroid = self$hclust.centroid, 
            mcquitty = self$hclust.mcquitty
		)
	}
)








##############################################################################
##########					Distance Functions						##########
##############################################################################

# distance function = 1-PCC (Pearson's correlation coefficient)
UtilFuns$set("public", "dist2",
	function(x,...){
		return(as.dist(1-cor(t(x), method="pearson")))
	}
)

# distance function = 1-abs(PCC) (Pearson's correlation coefficient)
UtilFuns$set("public", "dist3",
	function(x,...){
		return(as.dist(1-abs(cor(t(x), method="pearson"))))
	}
)






##############################################################################
##########				Hierarchical Clustering Functions			##########
##############################################################################

# average linkage
UtilFuns$set("public", "hclust2",
	function(x, method="average", ...){
		return(hclust(x, method=method, ...))
	}
)

# ward.D linkage
UtilFuns$set("public", "hclust.ward.D",
	function(x, method="ward.D", ...){
		return(hclust(x, method=method, ...))
	}
)

# ward.D2 linkage
UtilFuns$set("public", "hclust.ward.D2",
	function(x, method="ward.D2", ...){
		return(hclust(x, method=method, ...))
	}
)

# single linkage
UtilFuns$set("public", "hclust.single",
	function(x, method="single", ...){
		return(hclust(x, method=method, ...))
	}
)

# mcquitty linkage
UtilFuns$set("public", "hclust.mcquitty",
	function(x, method="mcquitty", ...){
		return(hclust(x, method=method, ...))
	}
)

# median linkage
UtilFuns$set("public", "hclust.median",
	function(x, method="median", ...){
		return(hclust(x, method=method, ...))
	}
)

# centroid linkage
UtilFuns$set("public", "hclust.centroid",
	function(x, method="centroid", ...){
		return(hclust(x, method=method, ...))
	}
)


# Find enriched TF binding motifs in promoters
UtilFuns$set("public", "promoter",
    function(converted,selectOrg, radio){
        speciesChoice <- LogicManager$DB$GetSpeciesChoice()

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
)