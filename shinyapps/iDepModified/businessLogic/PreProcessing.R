library('R6')
source('server.config')

PreProcessing.Logic <- R6Class("PreProcessing.Logic")



PreProcessing.Logic$set("public","ReadCountPreprocess",
	function(rawData){
		tmp.data <- self$RemoveNonNumericalColumns(rawData)
		tmp.data <- self$RemoveAllMissingRows(tmp.data)
		tmp.data <- self$CleanGeneIDs(tmp.data)
		tmp.data <- self$CleanSampleNames(tmp.data)
		tmp.data <- 
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
		return(dat)
	}
)







