library('R6')
library(RSQLite,verbose=FALSE)	# for database connection

source('server.config')

DB.Manager <- R6Class("DB.Manager")

DB.Manager$set("public","sqlite", NULL)
DB.Manager$set("public","dbConvert", NULL)
DB.Manager$set("public","OrgInfo", NULL)	# this field does not allow direct access, use 'QueryOrgInfo' to get org info
DB.Manager$set("public","Quotes", NULL)

## Init function
DB.Manager$set("public","initialize",
	function(){
		self$sqlite <- dbDriver("SQLite")
		
		self$dbConvert <- dbConnect( self$sqlite, paste0(CONFIG_DATA_DATAPATH, "convertIDs.db"), flags=SQLITE_RO) 
		
		dbQueryResult <- dbGetQuery(self$dbConvert, paste("select distinct * from orgInfo " ))
		self$OrgInfo <- dbQueryResult[order(dbQueryResult$name),]

		TmpQuotes <- dbGetQuery(self$dbConvert, " select * from quotes")
		self$Quotes <- paste0("\"",TmpQuotes$quotes,"\"", " -- ",TmpQuotes$author,".       ")

	}
)


# get id, ens, species from convertID mapping table
DB.Manager$set("public", "QuerySpeciesInfoFromConvertDBMapping",
	function(querySet){
		querySTMT <- paste( "select distinct id,ens,species from mapping where id IN ('", paste(querySet,collapse="', '"),"')",sep="")
		result <- dbGetQuery(self$dbConvert, querySTMT)
		return(result)
	}
)





