library('R6')
source('server.config')

PreProcessing.Logic <- R6Class("PreProcessing.Logic")



PreProcessing.Logic$set("public","ReadCountPreprocess",
	function(rawData, imputateMethod){
		tmp.data <- self$RemoveNonNumericalColumns(rawData)
		tmp.data <- self$RemoveAllMissingRows(tmp.data)
		tmp.data <- self$CleanGeneIDs(tmp.data)
		tmp.data <- self$SetGeneIDtoRowname(tmp.data)
		tmp.data <- self$CleanSampleNames(tmp.data)

		# sort by SD
		tmp.data <- tmp.data[order(- apply(tmp.data[,2:dim(tmp.data)[2]],1,sd) ),]
				
		if(sum(is.na(tmp.data))>0){
			tmp.data <- self$ImputateMissingValue(tmp.data,imputateMethod)
		}

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



