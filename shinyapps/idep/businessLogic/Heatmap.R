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


Heatmap.Logic$set("public", "ClusterGeneAndSample",
	function(dat, isSampleClustering){

		# clustering genes------
		clust <- dat %>% 
		 	LogicManager$UtilFuns$dist2() %>% 
		  	LogicManager$UtilFuns$hclust2()
		# Get order
		ord_row <- clust$order

		# clustering samples --------
		if( !isSampleClustering ){
			ord_column = 1:ncol(dat)
		} else { 
			clust <- t(dat) %>% 
			  	LogicManager$UtilFuns$dist2() %>% 
			  	LogicManager$UtilFuns$hclust2()
			# Get order
			ord_column <- clust$order
		}

		return(list(ord_row=ord_row, ord_column=ord_column))
	}
)

Heatmap.Logic$set("public", "GenerateHeatmapPlotly",
	function(dat, clusteredOrder, selectedHeatColor, allGeneInfo){
		ord_row=clusteredOrder$ord_row
		ord_column=clusteredOrder$ord_column

		# rename data use symbol
		rownames(dat) <- LogicManager$PreProcessing$GetGenesSymbolByEnsemblIDs(allGeneInfo, rownames(dat), TRUE)

		# Re-arrange based on order
		df <- t( dat[ord_row,ord_column] )%>%
	   		melt()
		colnames(df)[1:2] <- c("X","Y")

		# Call Display manager
		return(LogicManager$Display$GetHeatmapPlotly(df, selectedHeatColor))
	}
)

Heatmap.Logic$set("public", "CutData_SD",
	function(geneCount, ConvertedTransformedData){
		SDs=apply(ConvertedTransformedData,1,sd)
		maxSD = mean(SDs)+ 4*sd(SDs)
		SDs[ SDs > maxSD] = maxSD
		if(geneCount > length(SDs)){
				geneCount = length(SDs)
		}
		Cutoff=sort(SDs,decreasing=TRUE)[geneCount]
		SDs = as.data.frame(SDs)
		return(list(SDs=SDs, Cutoff=Cutoff))
	}
)

Heatmap.Logic$set("public", "GenerateSDHeatmapPlot",
	function(SDs, Cutoff){
		return(LogicManager$Display$GetSDHeatmapPlot(SDs, Cutoff))
	}
)


# Data cut function for CorrelationMatrixPlot
# remove bottom 25% lowly expressed genes, which inflate the PPC
# round the result by 2 digit
# melt and return the data
Heatmap.Logic$set("public", "CutDataForCorrelationMatrixPlot",
	function(transformedData){
		maxGene <- apply(transformedData,1,max)
		dat <- transformedData[which(maxGene > quantile(maxGene)[1] ) ,] # remove bottom 25% lowly expressed genes, which inflate the PPC
		melted_cormat <- melt(round(cor(dat),2), na.rm = TRUE)
		return(melted_cormat)
	}
)

Heatmap.Logic$set("public", "GenerateCorrelationPlot",
	function(dat, isLabelWithPCC){
		return(LogicManager$Display$GetHeatmapOfCorrelationMatrix(dat, isLabelWithPCC))
	}
)


# Data cut function for Sample Tree
# 1. remove bottom 25% lowly expressed genes, which inflate the PPC
# 2. center and standardize by gene
# 3. row centering and normalize
Heatmap.Logic$set("public", "CutDataForSampleTreePlot",
	function(transformedData, isGeneCentering, isGeneNormalize, isSampleCentering, isSampleNormalize){
		maxGene <- apply(transformedData,1,max)
		dat <- transformedData[which(maxGene > quantile(maxGene)[1] ) ,] # remove bottom 25% lowly expressed genes, which inflate the PPC

		# center by gene
		if(isGeneCentering){
			dat <- as.matrix(dat)-apply(dat,1,mean)	
		}
		
		# standardize by gene
		if(isGeneNormalize){
			dat <- dat / apply(dat,1,sd)
		}

		# row centering and normalize
		dat <- scale(dat, center = isSampleCentering, scale = isSampleNormalize) 
		
		return(dat)
	}
)

Heatmap.Logic$set("public", "GenerateSampLeTreePlot",
	function(dat, hclustfunName, distfunName){
		return(LogicManager$Display$GetSampLeTreePlot(dat, hclustfunName, distfunName))
	}
)



