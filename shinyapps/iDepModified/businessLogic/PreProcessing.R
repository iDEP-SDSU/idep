library('R6')
source('server.config')

PreProcessing.Logic <- R6Class("PreProcessing.Logic")



PreProcessing.Logic$set("public","ReadCountPreprocess",
	function(rawData){
		tmp.data <- self$RemoveNonNumericalColumns(rawData)
		tmp.data <- self$RemoveAllMissingRows(tmp.data)
		tmp.data <- self$CleanGeneIDs(tmp.data)
		tmp.data <- self$SetGeneIDtoRowname(tmp.data)
		tmp.data <- self$CleanSampleNames(tmp.data)

		# sort by SD
		tmp.data <- tmp.data[order(- apply(tmp.data[,2:dim(tmp.data)[2]],1,sd) ),]
		
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


PreProcessing.Logic$set("public", "CleanSampleNames",
	function(dat){

	}
)
