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



