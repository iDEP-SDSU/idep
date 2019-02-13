library('R6')
source('server.config')

PreProcessing.Logic <- R6Class("PreProcessing.Logic")



PreProcessing.Logic$set("public","ReadCountPreprocess",
	function(rawData){
		tmp.dataExcludeNonNumeric <- self$RemoveNonNumericalColumns(rawData)
		tmp.dataExcludeAllMissingColumn <- 
		tmp.dataRefineGeneId <-
		tmp.dataRefineGeneId <-
		tmp.dataAfterImpute <- 
	}
)

PreProcessing.Logic$set("public", "RemoveNonNumericalColumns",
	function(rawData){

		isNumeric <- apply(rawData,2,is.numeric)
		
		if( sum(isNumeric) <= 1 ){
			return(NULL)
		}

		isKeep <- isNumeric
		isKeep[1] <- TRUE

		return(rawData[,isKeep])
	}
)


