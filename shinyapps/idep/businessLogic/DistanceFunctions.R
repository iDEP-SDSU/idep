library('R6')
source('server.config')

DistanceFunctions.Manager <- R6Class("DistanceFunctions.Manager")

# distance function = 1-PCC (Pearson's correlation coefficient)
DistanceFunctions.Manager$set("public", "dist2",
	function(x,...){
		return(as.dist(1-cor(t(x), method="pearson")))
	}
)

# distance function = 1-abs(PCC) (Pearson's correlation coefficient)
DistanceFunctions.Manager$set("public", "dist3",
	function(x,...){
		return(as.dist(1-abs(cor(t(x), method="pearson"))))
	}
)


