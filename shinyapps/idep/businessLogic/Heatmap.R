library('R6')
library(reshape2,verbose=FALSE) 

source('server.config')

Heatmap.Logic <- R6Class("Heatmap.Logic")

Heatmap.Logic$set("public", "CutData",
	function(dat, geneCount, numHeatmapCutoff, isGeneCentering, 
		isGeneNormalize, isSampleCentering, isSampleNormalize)
	{
		if(geneCount>dim(dat)[1]){
			geneCount = dim(dat)[1]	# max	as data
		}

		if(geneCount < 10){
			geneCount = 10				# min
		} 

		# this will cutoff very large values, which could skew the color 
		if(isGeneCentering){
			dat <- as.matrix(dat[1:geneCount,]) - apply(dat[1:geneCount,], 1, mean)
		}

		# standardize by gene
		if(isGeneNormalize){
			dat <- dat/apply(dat,1,sd)
		}

		# row centering and normalize
		dat <- scale(dat, center = isSampleCentering, scale = isSampleNormalize)

		cutoff = median(unlist(dat)) + numHeatmapCutoff * sd(unlist(dat))
		dat[dat>cutoff] <- cutoff
		cutoff = median(unlist(dat)) - numHeatmapCutoff * sd(unlist(dat)) 
		dat[dat<cutoff] <- cutoff

		return(dat)
	}
)



Heatmap.Logic$set("public", "GenerateHeatmap",
	function(dat, sampleInfo, geneCount, isHaveSelectFactorHeatmap, 
		isSampleClustering, selectFactorsHeatmap, 
		selectedDistFunction, selectedhclustFunction,
		selectedHeatColor)
	{
		groups = LogicManager$PreProcessing$DetectGroups( colnames(dat) )
		mylist <- list(
			selectedDistFunction = selectedDistFunction, 
			selectedhclustFunction = selectedhclustFunction,
			selectedHeatColor = selectedHeatColor
		)

		saveRDS(file='mylist', mylist)

		if(!is.null(sampleInfo) && isHaveSelectFactorHeatmap ) {
			if(selectFactorsHeatmap == "Sample_Name" ){
				groups = LogicManager$PreProcessing$DetectGroups( colnames(dat) ) 
			}else{ 	
				ix = match( selectFactorsHeatmap, colnames(dat) ) 
				groups = dat[,ix]
			}
		}

		return(
			LogicManager$Display$GetHeatmap2Plot(dat, geneCount, groups, isSampleClustering,
				selectedDistFunction, selectedhclustFunction, selectedHeatColor)
		)
	}
)
