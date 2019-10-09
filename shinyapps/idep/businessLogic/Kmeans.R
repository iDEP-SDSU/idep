library('R6')
library(gplots,verbose=FALSE)		# for hierarchical clustering
library(ggplot2,verbose=FALSE)	# graphics

source('server.config')

Kmeans.Manager <- R6Class("Kmeans.Manager")

Kmeans.Manager$set("public", "CalcKmeansCluster",
    function( Reactive_ConvertedTransformedData, 
        GeneCount, 
        NormalizationMethod, 
        RerunSeed,
        NumberOfCluster )
    {
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
            dat = 100* dat / apply(dat,1,function(y) sum(abs(y))) else # L1 norm
        if( NormalizationMethod == 'geneMean')
            dat = dat - apply(dat,1,mean)  else # this is causing problem??????
        if( NormalizationMethod == 'geneStandardization')	
            dat = (dat - apply(dat,1,mean) ) / apply(dat,1,sd)
        
        set.seed(RerunSeed)
        
        cl = kmeans(x, NumberOfCluster, iter.max = CONST_KNN_ITERATION_MAX)
    
        hc <- hclust2(dist2(cl$centers-apply(cl$centers,1,mean) )  )# perform cluster for the reordering of samples
        tem = match(cl$cluster,hc$order) #  new order 
        x = x[order(tem),] 
        bar = sort(tem)
        
        return(list( SortedData = x , SortedIndex = bar) )
    }
)

Kmeans.Manager$set("public", "GetKmeansClusterHeatmapWithGeneBar",
	function(x, bar, color){
		centered_x <- x - apply(x,1,mean)
		return(LogicManager$Display$GetHeatmapWithGeneGroups(x, bar, 1000, mycolor = color))
	}
)

