library('R6')
library(gplots,verbose=FALSE)		# for hierarchical clustering
library(ggplot2,verbose=FALSE)	# graphics

source('server.config')
source('constant.config')

Kmeans.Logic <- R6Class("Kmeans.Logic")

Kmeans.Logic$set("public", "PickupAndNormalizeGene",
    function(Reactive_ConvertedTransformedData, GeneCount, NormalizationMethod){
        # Check GeneCount value
        if( GeneCount > CONST_KNN_MAX_GENE_CLUSTERING ){
            GeneCount = CONST_KNN_MAX_GENE_CLUSTERING # max
        }

        if( GeneCount < CONST_KNN_MIN_GENE_CLUSTERING ){
            GeneCount = CONST_KNN_MIN_GENE_CLUSTERING
        }

        DataRow <- dim(Reactive_ConvertedTransformedData)[1]
        if( GeneCount > DataRow){
            GeneCount = DataRow
        }
        
        # Pick up top GeneCount genes
        dat <- Reactive_ConvertedTransformedData[1:GeneCount,]

        # normalization
        if( NormalizationMethod == 'L1Norm')
            dat = 100* dat / apply(dat,1,function(y) sum(abs(y))) 
        else if( NormalizationMethod == 'geneMean')
            dat = dat - apply(dat,1,mean)  
        else if( NormalizationMethod == 'geneStandardization')	
            dat = (dat - apply(dat,1,mean) ) / apply(dat,1,sd)

        return(dat)
    }
)

Kmeans.Logic$set("public", "CalcKmeansCluster",
    function( Reactive_ConvertedTransformedData, 
        GeneCount, 
        NormalizationMethod, 
        RerunSeed,
        NumberOfCluster )
    {

        dat <- self$PickupAndNormalizeGene(Reactive_ConvertedTransformedData, GeneCount, NormalizationMethod)
        
        set.seed(RerunSeed)
        
        cl = kmeans(dat, NumberOfCluster, iter.max = CONST_KNN_ITERATION_MAX)
    
        hc <- LogicManager$UtilFuns$hclust2(LogicManager$UtilFuns$dist2(cl$centers-apply(cl$centers,1,mean) )  )# perform cluster for the reordering of samples
        tem = match(cl$cluster,hc$order) #  new order 
        dat = dat[order(tem),] 
        bar = sort(tem)
        
        return(list( x = dat , bar = bar) )
    }
)

Kmeans.Logic$set("public", "CalcKmeansNCluster",
    function(Reactive_ConvertedTransformedData, 
        GeneCount, 
        NormalizationMethod)
    {
        dat <- self$PickupAndNormalizeGene(Reactive_ConvertedTransformedData, GeneCount, NormalizationMethod)
        
        set.seed(2) ## need double check with Dr GE.

        k = CONST_KNN_MAX_CLUSTER_NUMBER
        wss <- (nrow(dat)-1)*sum(apply(dat,2,var))
        for (i in 2:k){
            wss[i] <- sum(kmeans(dat, centers = i, iter.max = CONST_KNN_ITERATION_MAX)$withinss)
        } 
        return(wss)
    }
)

Kmeans.Logic$set("public", "PlotWithinGroupSquareSum",
    function(wss){
        par(mar=c(4,5,4,4))
	    plot(1:CONST_KNN_MAX_CLUSTER_NUMBER, wss, type="b", xlab="Number of Clusters (k)",
		    ylab="Within groups sum of squares",
		    cex=2,cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2	,xaxt="n"	 )
	    axis(1, at = seq(1, 30, by = 2),cex.axis=1.5,cex=1.5)
    }
)


Kmeans.Logic$set("public", "MergeGenInfoWithClusterResult",
	function(x, bar, geneInfo, selectedOrg){
		Cluster <- toupper(letters)[bar]
		x <- cbind(Cluster, x)

		# add gene symbol
		if( selectedOrg != "NEW") { 
			ix <- match( rownames(x), geneInfo[,1])
			x <- cbind(as.character( geneInfo$symbol)[ix], x) 
		}

		return(x)
	}
)


Kmeans.Logic$set("public", "GetKmeansClusterHeatmapWithGeneBar",
	function(x, bar, color){
		centered_x <- x - apply(x,1,mean)
		return(LogicManager$Display$GetHeatmapWithGeneGroups(centered_x, bar, 1000, mycolor = color))
	}
)

Kmeans.Logic$set("public", "CalculateTSNEAndGeneratePlot",
	function( train, Cluster, seed, colorGenes ){
		library(Rtsne,verbose=FALSE)
		set.seed(seed)
		tsne <- Rtsne(train[,-1], dims = 2, perplexity=30, verbose=FALSE, max_iter = 400)
		nClusters = length(unique(Cluster) )

		if(colorGenes) {			
			plot(tsne$Y[,1], tsne$Y[,2], pch = (0:(nClusters-1))[Cluster], cex = 1.,col = mycolors[Cluster], xlab="X",ylab="Y")
			legend("topright",toupper(letters)[1:nClusters], pch = 0:(nClusters-1), col=mycolors, title="Cluster"  )
		} else {
			plot(tsne$Y[,1], tsne$Y[,2],  cex = 1., xlab="X",ylab="Y")
		}
	}
)


Kmeans.Logic$set("public", "CalculateGeneDistributionAndPlot",
	function(pickedGeneCount, data){
		SDs <- apply(data,1,sd)
		maxSD <- mean(SDs)+ 4*sd(SDs)
		SDs[ SDs > maxSD] = maxSD

		if(pickedGeneCount > length(SDs)){
			pickedGeneCount = length(SDs)
		}

		Cutoff=sort(SDs,decreasing=TRUE)[pickedGeneCount] 

		SDs = as.data.frame(SDs)

		p <- ggplot(SDs, aes(x=SDs)) + 
		  	geom_density(color="darkblue", fill="lightblue") +
		  	labs(x = "Standard deviations of all genes", y="Density")+
		  	geom_vline(
				aes(xintercept=Cutoff),
				color="red", 
				linetype="dashed", 
				size=1
			) +
			annotate("text", 
				x = Cutoff + 0.4*sd(SDs[,1]), 
				y = 1, 
				colour = "red", 
				label = paste0("Top ", pickedGeneCount)
			) +				
			theme(axis.text=element_text(size=14),
				axis.title=element_text(size=16,face="bold")
			)
		return(p)
	}
)


Kmeans.Logic$set("public", "GetKmeansGoData", 
    function(minFDR, selectedOrg, gmtFile, nCluster, GO, is_Kmeans_RemoveRedudantSets,
        Reactive_Kmeans, Reactive_ConvertedIDResult, Reactive_AllGeneInfo,
        Reactive_GeneSet ){
        pp = 0
        for( i in 1:nCluster) {
            
			query = rownames(Reactive_Kmeans$x)[which(Reactive_Kmeans$bar == i)]
			if(selectedOrg == "NEW" && !is.null( gmtFile) ){ 
				result <- LogicManager$DB$FindOverlapGMT( query, Reactive_GeneSet,1) 
                
			} else {
				convertedID <- Reactive_ConvertedIDResult
				convertedID$IDs <- query
				if(is_Kmeans_RemoveRedudantSets) 
                    reduced = CONST_redudantGeneSetsRatio 
                else 
                    reduced = FALSE
				result = LogicManager$DB$FindOverlap(convertedID, Reactive_AllGeneInfo, GO, selectedOrg,1,reduced) 
			}
			if( dim(result)[2] ==1) next;   # result could be NULL
			result$direction = toupper(letters)[i] 
			if (pp==0) 
            { 
                results <- result; 
                pp <- 1;
			} else {
				results <- rbind(results,result)
			}
		}

        if(pp == 0) {
            return(as.data.frame("No enrichment found."))
        }
        
		results= results[,c(6,1,2,4,5)]
		colnames(results)= c("Cluster","FDR","nGenes","Pathways","Genes")
		if(min(results$FDR) > minFDR ){
            results = as.data.frame("No signficant enrichment found.") 
        } else {
            results = results[which(results$FDR < minFDR),]
        }
		
        if( is.null(results) )  return ( as.matrix("No significant enrichment.") )	
        if( class(results) != "data.frame")  return ( as.matrix("No significant enrichment.") )
        if( dim(results)[2] ==1)  return ( as.matrix("No significant enrichment.") )
        colnames(results)[2] = "adj.Pval"

        return(results)
    }
)

Kmeans.Logic$set("public", "GetGeneSetByGOOption",
	function( ConvertedIDResult, ConvertedTransformedData, selectOrg, gmtFile, GO, maxSetSize){
		if(is.null(ConvertedTransformedData) | is.null(ConvertedIDResult) ) {
			return(NULL)
		}

		if( selectOrg == "NEW" | !is.null(gmtFile) ){ # new species 
			inFile <- gmtFile$datapath
			return( DB$readGMTRobust(inFile) )
		}

		return(DB$QueryGeneSetsFromPathway(ConvertedIDResult, ConvertedTransformedData, selectOrg, gmtFile, GO, maxSetSize))
	}
)

Kmeans.Logic$set("public", "CalculateKmeansGoTableData",
    function(KmeansGOData, is_Kmeans_RemoveRedudantSets){
        if(is.null(KmeansGOData)) {
            return(NULL)
        }
            
        results1 = KmeansGOData
        if(dim(results1)[2] == 1) {
            return(results1) 
        }else{
            results1$adj.Pval <- sprintf("%-2.1e",as.numeric(results1$adj.Pval) )
            results1[,1] <- as.character(results1[,1])
            results1[ duplicated (results1[,1] ),1 ] <- ""  
            
            return( results1[,-5])
	    }
    }
)

Kmeans.Logic$set("public", "GenerateEnrichmentPlot",
    function(selectedGO, selectedOrg, gmtFile, Reactive_Kmeans, Reactive_KmeansGOData ){
        if(is.null(Reactive_KmeansGOData)){
            return(NULL) 
        }

        if(is.null(selectedGO)){
            return (NULL)
        }

        if(selectedGO == "ID not recognized!"){
            return ( as.matrix("Gene ID not recognized.") )#No matching species
        } 

   	    if( is.null(Reactive_Kmeans) ) {
            return(NULL)
        }

	    if( selectedOrg == "NEW" && is.null(gmtFile) ) {
            return(NULL) # new but without gmtFile
        }

       	tem1 = KmeansGOdata()
        colnames(tem1)[1]="Direction"
        return( LogicManager$Display$GetEnrichmentPlot(tem1, 46) )
    }
)


Kmeans.Logic$set("public", "PromoterAnalysis",
    function(nClusters, selectOrg, selectGO2, promoterBP, Reactive_Kmeans){
        results1 <- NULL
        result <- NULL 
	    pp<- 0
	    for( i in 1:nClusters ) {
	        
            #query = rownames(x)[which(bar == i)]
            query = rownames(Reactive_Kmeans$x)[which(Reactive_Kmeans$bar == i)]	
            convertedID = LogicManager$PreProcessing$GetConvertID(query, selectOrg);#"gmax_eg_gene"
            result <- LogicManager$UtilFuns$promoter( convertedID, selectOrg, promoterBP )
            
            if( is.null(result) ){
                next   # result could be NULL
            }

            if( dim(result)[2] ==1) {
                next
            }

            result$List = toupper(letters)[i]    
            if (pp==0 ) { 
                results1 <- result; 
                pp <- 1 
            } else { 
                results1 = rbind(results1,result) 
            }
        }

        if( is.null(results1) ) {
            return( as.data.frame("No significant motif enrichment found.") )
        } else {
            results1[ duplicated (results1[,4] ),4 ] <- ""
            return( results1[,c(4,1:3,5)] )
        }
    }
)