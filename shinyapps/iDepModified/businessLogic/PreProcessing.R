library('R6')

library(e1071,verbose=FALSE) 		# computing kurtosis

source('server.config')

PreProcessing.Logic <- R6Class("PreProcessing.Logic")


PreProcessing.Logic$set("public","ReadDataPreprocess",
	# Data file preprocessing function
	# Major steps:
	# 	1. Remove non-numerical columns
	#		Note: The first column won't be remove. We by default consider the first column as Gene ID
	#	2. Remove all missing rows
	#		If a row (a Gene) does not contain any information, we don't need it shows in the dataset
	#	3. Clean Gene ID and set it as row name
	#		Gene Ids sometime contain invalid symbol, which should be clean.
	#		For most bioinfo R package, Gene ID is the rowname of the dataset. 
	#		So we need set our dataset rowname as gene id
	#	4. Clean Sample name.
	#		This related to some later operations.
	#	5. Impute missing data
	#	6. Calculate Kurtosis
	function(rawData, imputateMethod){
		tmp.data <- self$RemoveNonNumericalColumns(rawData)
		tmp.data <- self$RemoveAllMissingRows(tmp.data)
		dataSizeOriginal = dim(tmp.data); 
		dataSizeOriginal[2] = dataSizeOriginal[2] - 1
		
		tmp.data <- self$CleanGeneIDs(tmp.data)
		tmp.data <- self$SetGeneIDtoRowname(tmp.data)
		tmp.data <- self$CleanSampleNames(tmp.data)

		# sort by SD
		tmp.data <- tmp.data[order(- apply(tmp.data[,2:dim(tmp.data)[2]],1,sd) ),]
				
		if(sum(is.na(tmp.data))>0){
			tmp.data <- self$ImputateMissingValue(tmp.data,imputateMethod)
		}

		result <- CalcKurtosis(dat, dataType, minCounts, NminSamples, CountsTransform, 
				LogStart, LowFilter, isTransform, isNoFDR)

		dataSizeAfter <- dim(result$dat)
		dataSizeAfter[2] = dataSizeAfter[2] - 1

		return(list(data = result$dat, mean.kurtosis = result$mean.kurtosis,
					rawCount = result$rawCount, dataTypeWarning = result$warning,
					dataSizeOriginal = dataSizeOriginal, dataSizeAfter = dataSizeAfter,
					pvals=result$pvals))
	}
)

PreProcessing.Logic$set("public", "RemoveNonNumericalColumns",
	function(rawData){
	# Remove non numeric column
		isNumeric <- apply(rawData,2,is.numeric)

		if( sum(isNumeric) <= 1 ){
			return(NULL)
		}

		idxKeep <- isNumeric
		idxKeep[1] <- TRUE

		return(rawData[,idxKeep])
	}
)


PreProcessing.Logic$set("public", "RemoveAllMissingRows",
	## Remove rows with only missing value
	function(dat){
		idxKeep = which( apply(dat[,-1],1, function(y) sum( is.na(y) ) ) != dim(dat)[2]-1 )
		
		return(dat[idxKeep,])
	}
)


PreProcessing.Logic$set("public", "CleanGeneIDs",
	function(dat){
		# remove spaces in gene ids
		# remove " in gene ids, mess up SQL query				
		# remove ' in gene ids		
		# remove one or two digits after "." at the end.
    	# A35244.1 -> A35244  or A35244.23 -> A35244, but not more than two.  GLYMA.18G52160 stays the same.
		dat[,1] <- toupper(dat[,1])
		dat[,1] <- gsub(" |\"|\'|\\.[0-9]{1,2}$", "", dat[,1])

		# remove duplicated genes
		dat <- dat[!duplicated(dat[,1]) ,]

		return(dat)
	}
)

PreProcessing.Logic$set("public", "SetGeneIDtoRowname",
	function(dat){	
		rownames(dat) <- dat[,1]
		return(as.matrix(dat[,c(-1)]))
	}
)


PreProcessing.Logic$set("public", "CleanSampleNames",
	function(dat){	
		# remove "-" or "." from sample names
		colnames(dat) = gsub("-","",colnames(dat))
		colnames(dat) = gsub("\\.","",colnames(dat))
		return(dat)
	}
)

PreProcessing.Logic$set("public", "ImputateMissingValue",
	function(dat, method){
		avaiableMethod = c("geneMedian", "treatAsZero", "geneMedianInGroup")
		
		methodIdx <- which(avaiableMethod, method)

		if(length(methodIdx) != 1 ){
			return(dat)
		}

		dat <- switch(	methodIdx, 
			# case "geneMedian":
			self$FillViaGeneMedian(dat),	
			
			# case "treatAsZero":
			self$TreatAsZero(dat),
			
			# case "geneMedianInGroup":
			self$FillViaGroupMedian(dat)
		)
		
		return(dat)
	}
)


PreProcessing.Logic$set("public", "FillViaGeneMedian",
	function(dat){
		rowMedians <- apply(dat,1, median,na.rm=T)
		for( i in 1:dim(dat)[2] ) {
			ix = which(is.na(dat[,i]) )
			dat[ix,i] <- rowMedians[ix]						
		}
		return(dat)
	}
)

PreProcessing.Logic$set("public", "TreatAsZero",
	function(dat){
		dat[is.na(dat)] <- 0	
		return(dat)
	}
)

PreProcessing.Logic$set("public", "FillViaGroupMedian",
	function(dat){
		sampleGroups = detectGroups( colnames(dat))
		for (group in unique(sampleGroups) ){		
			samples = which( sampleGroups == group )
			rowMedians <- apply(dat[,samples, drop=F],1, median, na.rm=T)
			for( i in samples ) { 
				idx = which(is.na(dat[ ,i] ) )	
				if(length(idx) >0 )
					dat[idx, i]  <- rowMedians[idx]
			}										
		}
						
		# missing for entire sample group, use median for all samples
		if(sum(is.na(dat))>0 ) { 
			rowMedians <- apply(dat,1, median,na.rm=T)
			for( i in 1:dim(dat)[2] ) {
				idx = which(is.na(dat[,i]) )
				dat[idx,i] <- rowMedians[idx]			
			}						
		}
		return(dat)
	}
)


PreProcessing.Logic$set("public", "DetectGroups",
	function(x){
		# x are col names
		# Define sample groups based on column names
		# Args:
		#   x are vector of characters, column names in data file
		# Returns: 
		#   a character vector, representing sample groups.
		tem <- gsub("[0-9]*$","",x) # Remove all numbers from end
		#tem = gsub("_Rep|_rep|_REP","",tem)
		tem <- gsub("_$","",tem); # remove "_" from end
		tem <- gsub("_Rep$","",tem); # remove "_Rep" from end
		tem <- gsub("_rep$","",tem); # remove "_rep" from end
		tem <- gsub("_REP$","",tem)  # remove "_REP" from end
		return( tem )
	}
)

PreProcessing.Logic$set("public", "CalcKurtosis",
	function(dat, dataType, minCounts=NULL, NminSamples=NULL, CountsTransform=NULL, 
				LogStart=NULL, LowFilter=NULL, isTransform=FALSE, isNoFDR=TRUE ){
		mean.kurtosis = mean(apply(x,2, kurtosis),na.rm=T)
		rawCounts = NULL
		pvals= NULL		
		result <- switch(dataType,
			self$CalcKurtosisForReadCount(dat, mean.kurtosis, minCounts, NminSamples, CountsTransform, LogStart),
			self$CalcKurtosisForFPKM(dat, mean.kurtosis, LowFilter, NminSamples, isTransform, LogStart),
			self$CalcKurtosisForOtherDatatype(dat, isNoFDR)
		)
		return(result)
	}

)

PreProcessing.Logic$set("public", "CalcKurtosisForReadCount",
	function(dat, mean.kurtosis, minCounts, NminSamples, CountsTransform, CountsLogStart){
		if(!is.integer(dat) & (mean.kurtosis < CONST_KRUTOSIS_LOG ) ) {
			dataTypeWarning = -1
		}

		dat <- round(dat,0)
		dat <- dat[ which( apply( cpm(DGEList(counts = dat)), 1,  
							function(x) sum(x >=minCounts)) >= NminSamples ) , ]

		# construct DESeqExpression Object
		# colData = cbind(colnames(x), as.character(detectGroups( colnames(x) )) )
		tem = rep("A",dim(x)[2]); tem[1] <- "B"   # making a fake design matrix to allow process, even when there is no replicates
		colData = cbind(colnames(x), tem )
		colnames(colData)  = c("sample", "groups")
		dds <- DESeqDataSetFromMatrix(countData = x, colData = colData, design = ~ groups)
		dds <- estimateSizeFactors(dds) # estimate size factor for use in normalization later for started log method

		dat <- switch(CountsTransform,
			log2( counts(dds, normalized=TRUE) + CountsLogStart ) ,
			assay( rlog(dds, blind=TRUE) ),
			assay( vst(dds, blind=TRUE) )
		)
		
		return(list(dat=dat, warning=dataTypeWarning))
	}
)


PreProcessing.Logic$set("public", "CalcKurtosisForFPKM",
	function(dat, mean.kurtosis, LowFilter, NminSamples, isTransform, LogStart){
		if ( is.integer(dat) ) {
			dataTypeWarning = 1  # Data appears to be read counts
		}
		
		#-------------filtering
		#tem <- apply(x,1,max)
		#x <- x[which(tem > input$lowFilter),]  # max by row is at least 
		dat <- dat[ which( apply( dat, 1,  function(dat) sum(x >= LowFilter)) >= NminSamples ) , ] 
		
		dat <- dat[which( apply(dat,1, function(x) max(x))- min(x) ) > 0  ),]  # remove rows with all the same levels
		#--------------Log transform
		# Takes log if log is selected OR kurtosis is big than 100
		if ( (isTransform) | ( mean.kurtosis > CONST_KRUTOSIS_LOG ) ) {
			rawCount <- dat
			dat = log( dat + abs( LogStart ),2)
		}
			
		tem <- apply(dat,1,sd) 
		dat <- dat[order(-tem),]  # sort by SD

		return(list(dat=dat, warning=dataTypeWarning, rawCount=rawCount))
	}
)

PreProcessing.Logic$set("public", "CalcKurtosisForOtherDatatype",
	function(dat, isNoFDR){
		n2 = ( dim(dat)[2] %/% 2) # 5 --> 2
		if(!isNoFDR) { 	 # if we have FDR
			pvals = dat[,2*(1:n2 ),drop=FALSE ]  # 2, 4, 6
			dat = dat[, 2*(1:n2 )-1,drop=FALSE]   # 1, 3, 5				
		}
		return(list(dat = dat, pvals=pvals))
	}
)



