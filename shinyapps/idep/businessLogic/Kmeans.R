library('R6')
library(gplots,verbose=FALSE)		# for hierarchical clustering
library(ggplot2,verbose=FALSE)	# graphics

source('server.config')

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
        
        cl = kmeans(x, NumberOfCluster, iter.max = CONST_KNN_ITERATION_MAX)
    
        hc <- hclust2(dist2(cl$centers-apply(cl$centers,1,mean) )  )# perform cluster for the reordering of samples
        tem = match(cl$cluster,hc$order) #  new order 
        x = x[order(tem),] 
        bar = sort(tem)
        
        return(list( x = x , bar = bar) )
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
		return(LogicManager$Display$GetHeatmapWithGeneGroups(x, bar, 1000, mycolor = color))
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

