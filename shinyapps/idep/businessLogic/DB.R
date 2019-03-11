
sqlite  <- dbDriver("SQLite")
convert <- dbConnect( sqlite, paste0(CONFIG_DATA_DATAPATH, "convertIDs.db"), flags=SQLITE_RO) 


QuerySpeciesInfoFromConvertDBMapping(querySet){
	querySTMT <- paste( "select distinct id,ens,species from mapping where id IN ('", paste(querySet,collapse="', '"),"')",sep="")
	result <- dbGetQuery(convert, querySTMT)
	return(result)
}

# Prepare species list
QuerySpeciesListFromConvertDBOrgInfo(){

	# Create a list for Select Input options
	orgInfo <- dbGetQuery(convert, paste("select distinct * from orgInfo " ))
	orgInfo <- orgInfo[order(orgInfo$name),]
	speciesChoice <- setNames(as.list( orgInfo$id ), orgInfo$name2 )
	# add a defult element to list    # new element name       value
	speciesChoice <- append( setNames( "NEW","**NEW SPECIES**"), speciesChoice  )
	speciesChoice <- append( setNames( "BestMatch","Best matching species"), speciesChoice  )

}
