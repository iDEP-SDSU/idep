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


