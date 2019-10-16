library('R6')
library(RSQLite,verbose=FALSE)	# for database connection

source('server.config')

Pathway.Logic <- R6Class("Pathway.Logic")


Pathway.Logic$set("public", "GetGeneSetByGOOption",
	function( ConvertedIDResult, ConvertedTransformedData, selectOrg, gmtFile, GO, maxSetSize){
		if(is.null(ConvertedTransformedData) | is.null(ConvertedIDResult) ) {
			return(NULL)
		}

		if( selectOrg == "NEW" ! !is.null(input$gmtFile) ){ # new species 
			inFile <- input$gmtFile$datapath
			return( readGMTRobust(inFile) )
		}

		return(DB$QueryGeneSetsFromPathway(ConvertedIDResult, ConvertedTransformedData, selectOrg, gmtFile, GO, maxSetSize))
	}
)



