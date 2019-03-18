library('R6')
source('server.config')

HierarchicalClusteringFunctions.Manager <- R6Class("HierarchicalClusteringFunctions.Manager")

HierarchicalClusteringFunctions.Manager$set("public", "hclustFuns", NULL)

# average linkage
HierarchicalClusteringFunctions.Manager$set("public", "hclust2",
	function(x, method="average", ...){
		return(hclust(x, method=method, ...))
	}
)

# ward.D linkage
HierarchicalClusteringFunctions.Manager$set("public", "hclust.ward.D",
	function(x, method="ward.D", ...){
		return(hclust(x, method=method, ...))
	}
)

# ward.D2 linkage
HierarchicalClusteringFunctions.Manager$set("public", "hclust.ward.D2",
	function(x, method="ward.D2", ...){
		return(hclust(x, method=method, ...))
	}
)

# single linkage
HierarchicalClusteringFunctions.Manager$set("public", "hclust.single",
	function(x, method="single", ...){
		return(hclust(x, method=method, ...))
	}
)

# mcquitty linkage
HierarchicalClusteringFunctions.Manager$set("public", "hclust.mcquitty",
	function(x, method="mcquitty", ...){
		return(hclust(x, method=method, ...))
	}
)

# median linkage
HierarchicalClusteringFunctions.Manager$set("public", "hclust.median",
	function(x, method="median", ...){
		return(hclust(x, method=method, ...))
	}
)

# centroid linkage
HierarchicalClusteringFunctions.Manager$set("public", "hclust.centroid",
	function(x, method="centroid", ...){
		return(hclust(x, method=method, ...))
	}
)

HierarchicalClusteringFunctions.Manager$set("public", "initialize",
 	function(){
		self$hclustFuns <- list(
			average   = self$hclust2, 
            complete = hclust, 
            single   = self$hclust.single,
			median   = self$hclust.median, 
            centroid = self$hclust.centroid, 
            mcquitty = self$hclust.mcquitty
		)
	}
)
