library('R6')
library(reshape2,verbose=FALSE) 

source('server.config')

Heatmap.Logic <- R6Class("Heatmap.Logic")

Heatmap.Logic$set("public", "distFuns", NULL)
Heatmap.Logic$set("public", "hclustFuns", NULL)

Heatmap.Logic$set("public", "initialize", 
	function(DistanceFunctions, HierarchicalClusteringFunctions){
		self$distFuns <- list(
			Correlation=DistanceFunctions$dist2, 
			Euclidean=dist, 
			AbsolutePCC=DistanceFunctions$dist3)
		self$hclustFuns <- HierarchicalClusteringFunctions$hclustFuns
	}
)



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

# adding sample legends to heatmap; this is for the main heatmap
# https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
Heatmap.Logic$set("public", "add_legend",
	function(...) {
		opar <- par(
			fig=c(0, 1, 0, 1), 
			oma=c(0, 0, 0, 0), 
		 	mar=c(0, 0, 0, 6), 
			new=TRUE
		)
		on.exit(par(opar))
		plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
		legend(...)
	}
)


Heatmap.Logic$set("public", "GenerateHeatmap",
	function(dat, sampleInfo, geneCount, isHaveSelectFactorHeatmap, 
		isSampleClustering, selectFactorsHeatmap, 
		selectedDistFunction, selectedhclustFunction,
		selectedHeatColor)
	{
		groups = LogicManager$PreProcessing$DetectGroups(colnames(dat))

		if(!is.null(sampleInfo) && isHaveSelectFactorHeatmap ) {
			if(selectFactorsHeatmap == "Sample_Name" ){
				groups = LogicManager$PreProcessing$DetectGroups(colnames(dat) ) 
			}else{ 	
				ix = match(selectFactorsHeatmap, colnames(dat) ) 
				groups = dat[,ix]
			}
		}

		groups.colors = rainbow(length(unique(groups)))

		lmat = rbind(c(0,4),c(0,1),c(3,2),c(5,0))
		lwid = c(2,6)
		lhei = c(1.5,.2,8,1.1)
	
		if(ncol(dat) < 20){
			cexFactor = 2 
		}else if(ncol(dat) < 31){
			cexFactor = 1.5
		}else{
			cexFactor = 1
		}

		par(mar = c(5, 4, 1.4, 0.2))

		if( geneCount>110 ){
			heatmap.2(dat, 
				distfun = self$distFuns[[as.integer(selectedDistFunction)]],
				hclustfun = self$hclustFuns[[selectedhclustFunction]],
				Colv = isSampleClustering,
				col = LogicManager$Display$heatColors[as.integer(selectedHeatColor),],
				density.info="none", 
				trace="none", 
				scale="none", 
				keysize=.5,
				key=TRUE, 
				symkey=FALSE,
				ColSideColors=groups.colors[ as.factor(groups)],
				labRow="",
				margins=c(10,0),
				srtCol=45,
				cexCol=cexFactor,  # size of font for sample names
				lmat = lmat, 
				lwid = lwid, 
				lhei = lhei
			)
		}else{
			heatmap.2(dat, 
				distfun = self$distFuns[[as.integer(selectedDistFunction)]],
				hclustfun=hclustFuns[[selectedhclustFunction]],
				Colv=isSampleClustering,
				col= LogicManager$Display$heatColors[as.integer(selectedHeatColor),],
				density.info="none", 
				trace="none", 
				scale="none", 
				keysize=.5,
				key=TRUE, 
				symkey=FALSE,
				ColSideColors=groups.colors[ as.factor(groups)],
				margins=c(18,12),
				cexRow=1,
				srtCol=45,
				cexCol=cexFactor,  # size of font for sample names
				lmat = lmat, 
				lwid = lwid, 
				lhei = lhei
			)
		}

		if(length(unique(groups) ) <= 30 ) {  # only add legend when there is less categories
			par(lend = 1)           # square line ends for the color legend
			self$add_legend(
				"topleft",
				legend = unique(groups), # category labels
				col = groups.colors[unique(as.factor(groups))],  # color key
				lty= 1,             # line style
				lwd = 10            # line width
			)
		}
	}
)
