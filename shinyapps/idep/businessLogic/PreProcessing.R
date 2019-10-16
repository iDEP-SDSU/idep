library('R6')
library('edgeR')
library('DESeq2')
library('dplyr')
library(e1071,verbose=FALSE) 		# computing kurtosis

source('server.config')

PreProcessing.Logic <- R6Class("PreProcessing.Logic")

PreProcessing.Logic$set("public","RawDataPreprocess",
	# Data file preprocessing function
	# Convert from Server.R$readData
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
	function(rawData, imputateMethod, dataType, minCounts, minSamples, 
			CountsTransform, LogStart, isTransform, isNoFDR){

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

		result <- self$CalcKurtosis(tmp.data, dataType, minCounts, minSamples, 
						CountsTransform, LogStart, isTransform, isNoFDR)

		
		dataSizeAfter <- dim(result$dat)
		dataSizeAfter[2] = dataSizeAfter[2] - 1

		return(list(dat = result$dat, mean.kurtosis = result$mean.kurtosis,
					rawCount = result$rawCount, dataTypeWarning = result$warning,
					dataSizeOriginal = dataSizeOriginal, dataSizeAfter = dataSizeAfter,
					pvals=result$pvals))
	}
)

#	support function for RawDataPreprocess	
PreProcessing.Logic$set("public", "RemoveNonNumericalColumns",
	function(rawData){
	# Remove non numeric column
		# isNumeric <- apply(rawData,2,is.numeric)
		# apply does not work :
		# if any of the columns in your dataframe is not numeric, apply will try 
		# to coerce all of them to the least common supertype, and you'll get 
		# FALSE for each column;  this is not the case with sapply. 

		rawData <- as.data.frame(rawData)
		isNumeric <- sapply(rawData, is.numeric)

		if( sum(isNumeric) <= 1 ){
			return(NULL)
		}

		idxKeep <- isNumeric
		idxKeep[1] <- TRUE

		return(rawData[,idxKeep])
	}
)

#	support function for RawDataPreprocess
PreProcessing.Logic$set("public", "RemoveAllMissingRows",
	## Remove rows with only missing value
	function(dat){

		idxKeep = which( apply(dat[,-1],1, function(y) sum( is.na(y) ) ) != dim(dat)[2]-1 )
		
		return(dat[idxKeep,])
	}
)

#	support function for RawDataPreprocess
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

#	support function for RawDataPreprocess
PreProcessing.Logic$set("public", "SetGeneIDtoRowname",
	function(dat){	
		rownames(dat) <- dat[,1]
		return(as.matrix(dat[,c(-1)]))
	}
)

#	support function for RawDataPreprocess
PreProcessing.Logic$set("public", "CleanSampleNames",
	function(dat){	
		# remove "-" or "." from sample names
		colnames(dat) = gsub("-","",colnames(dat))
		colnames(dat) = gsub("\\.","",colnames(dat))
		return(dat)
	}
)

#	support function for RawDataPreprocess
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

#	support function for ImputateMissingValue
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

#	support function for ImputateMissingValue
PreProcessing.Logic$set("public", "TreatAsZero",
	function(dat){
		dat[is.na(dat)] <- 0	
		return(dat)
	}
)

#	support function for ImputateMissingValue
PreProcessing.Logic$set("public", "FillViaGroupMedian",
	function(dat){
		sampleGroups = self$DetectGroups(colnames(dat))
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
	function(dat){
		# x are col names
		# Define sample groups based on column names
		# Args:
		#   x are vector of characters, column names in data file
		# Returns: 
		#   a character vector, representing sample groups.
		tem <- gsub("[0-9]*$","",dat) # Remove all numbers from end
		#tem = gsub("_Rep|_rep|_REP","",tem)
		tem <- gsub("_$","",tem); # remove "_" from end
		tem <- gsub("_Rep$","",tem); # remove "_Rep" from end
		tem <- gsub("_rep$","",tem); # remove "_rep" from end
		tem <- gsub("_REP$","",tem)  # remove "_REP" from end
		return( tem )
	}
)

#	support function for RawDataPreprocess	
PreProcessing.Logic$set("public", "CalcKurtosis",
	function(dat, dataType, minCounts=NULL, minSamples=NULL, CountsTransform=NULL, 
				LogStart=NULL, isTransformFPKM=FALSE, isNoFDR=TRUE ){
		mean.kurtosis = mean(apply(dat, 2, kurtosis),na.rm=T)
		rawCounts = NULL
		pvals= NULL
				
		result <- switch(as.numeric(dataType),
			self$CalcKurtosisForReadCount(dat, mean.kurtosis, minCounts, minSamples, CountsTransform, LogStart),
			self$CalcKurtosisForFPKM(dat, mean.kurtosis, minCounts, minSamples, isTransformFPKM, LogStart),
			self$CalcKurtosisForOtherDatatype(dat, isNoFDR)
		)
		return(result)
	}

)
#	support function for CalcKurtosis	
PreProcessing.Logic$set("public", "CalcKurtosisForReadCount",
	function(dat, mean.kurtosis, minCounts, NminSamples, CountsTransform, CountsLogStart){
		if(!is.integer(dat) & (mean.kurtosis < CONST_KRUTOSIS_LOG ) ) {
			dataTypeWarning = -1
		}else{
			dataTypeWarning = NULL
		}
		
		dat <- round(dat,0)
		dat <- dat[ which( apply( cpm(DGEList(counts = dat)), 1,  
							function(x) sum(x >=minCounts)) >= NminSamples ) , ]
		
		rawCount = dat

		# construct DESeqExpression Object
		# colData = cbind(colnames(x), as.character(detectGroups( colnames(x) )) )
		tem = rep("A",dim(dat)[2]); tem[1] <- "B"   # making a fake design matrix to allow process, even when there is no replicates
		colData = cbind(colnames(dat), tem )
		colnames(colData)  = c("sample", "groups")
		dds <- DESeqDataSetFromMatrix(countData = dat, colData = colData, design = ~ groups)
		dds <- estimateSizeFactors(dds) # estimate size factor for use in normalization later for started log method

		dat <- switch(CountsTransform,
			log2( counts(dds, normalized=TRUE) + CountsLogStart ) ,
			assay( vst(dds, blind=TRUE) ),
			assay( rlog(dds, blind=TRUE) )
		)
		
		return(list(dat=dat, warning=dataTypeWarning, 
					mean.kurtosis = mean.kurtosis, rawCount=rawCount,
					pvals = NULL))
	}
)

#	support function for CalcKurtosis	
PreProcessing.Logic$set("public", "CalcKurtosisForFPKM",
	function(dat, mean.kurtosis, LowFilter, NminSamples, isTransform, LogStart){
		if ( is.integer(dat) ) {
			dataTypeWarning = 1  # Data appears to be read counts
		}
		
		#-------------filtering
		#tem <- apply(x,1,max)
		#x <- x[which(tem > input$lowFilter),]  # max by row is at least 
		dat <- dat[ which( apply( dat, 1,  function(x) sum(x >= LowFilter)) >= NminSamples ), ] 
		
		dat <- dat[ which( apply( dat, 1, function(x) max(x)- min(x) ) > 0 ), ]  # remove rows with all the same levels
		#--------------Log transform
		# Takes log if log is selected OR kurtosis is big than 100
		if ( (isTransform) | ( mean.kurtosis > CONST_KRUTOSIS_LOG ) ) {
			rawCount <- dat
			dat = log( dat + abs( LogStart ),2)
		}
			
		tem <- apply(dat,1,sd) 
		dat <- dat[order(-tem),]  # sort by SD

		return(list(dat=dat, warning=dataTypeWarning, 
					mean.kurtosis = mean.kurtosis, rawCount=rawCount,
					pvals = NULL))
	}
)

#	support function for CalcKurtosis	
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


PreProcessing.Logic$set("public","RawSampleInfoPreprocess", 
	function(rawInfo, rawData){
		if(is.null(rawInfo)){
			return(NULL)
		}
		
		colnameOfData <- colnames(rawData)
		temp.info <- self$CleanSampleNames(rawInfo)

		matchedInfo <- match(toupper(colnameOfData), toupper(colnames(temp.info)))
		matchedInfo <- matchedInfo[which(!is.na(matchedInfo))]

		validate(
				  need(	
					length(unique(matchedInfo) ) == dim(rawData)[2] 
				       	& dim(rawInfo)[1]>=1 & dim(rawInfo)[1] <500, # at least one row, it can not be more than 500 rows
					"Error!!! Sample information file not recognized. Sample names must be exactly the same. Each row is a factor. Each column represent a sample.  Please see documentation on format.")
		)

		#-----------Double check factor levels, change if needed
		# remove "-" or "." from factor levels
		for( i in 1:dim(temp.info)[1]) {
		   temp.info[i,] = gsub("-","",temp.info[i,])
		   temp.info[i,] = gsub("\\.","",temp.info[i,])				
		}
		
		if( length(unique(matchedInfo) ) == dim(rawData)[2]) { # matches exactly
			temp.info = temp.info[,matchedInfo]
			# if the levels of different factors are the same, it may cause problems
			if( sum( apply(temp.info, 1, function(x) length(unique(x)))) > length(unique(unlist(temp.info) ) ) ) {
				tem2 =apply(temp.info,2, function(x) paste0( names(x),x)) # factor names are added to levels
				rownames(tem2) = rownames(temp.info)
				temp.info <- tem2				
			}
			return(t(temp.info))			
		}else{
			retrun(NULL)
		}
	}
)


# convert gene IDs to ensembl gene ids and find species
# steps:
#	1. clean gene name
#	2. query from convertID.db/mapping table
#	3. get matched species
#	4. clean result table by removing duplicates and changing colnames
#
PreProcessing.Logic$set("public", "GetConvertID",
	function(geneNames, selectOrg) {
		querySet <- self$cleanGeneSet( unlist( strsplit( toupper(geneNames),'\t| |\n|\\,')))
		# querySet is ensgene data for example, ENSG00000198888, ENSG00000198763, ENSG00000198804

		result <- LogicManager$DB$QuerySpeciesInfoFromConvertDBMapping(querySet)
		# this will return a table with id,ens,species


		if( dim(result)[1] == 0  ){
		# if we didn't have a table return null
			return(NULL)
		} 

		# get matched species
		if(selectOrg == 'BestMatch') {
			speciesMatched <- self$ConvertIDAutoMatch(result)
		} else { # if species is selected
			result <- result[which(result$species == selectOrg ) ,]
			if( dim(result)[1] == 0  ) return(NULL) #stop("ID not recognized!")
			speciesMatched <- as.data.frame(paste("Using selected species ", LogicManager$DB$findSpeciesNameById(selectOrg) )  )
		}


		result <- result[which(!duplicated(result[,2]) ),] # remove duplicates in ensembl_gene_id
		result <- result[which(!duplicated(result[,1]) ),] # remove duplicates in user ID
		colnames(speciesMatched) = c("Matched Species (genes)") 
		conversionTable <- result[,1:2]
		colnames(conversionTable) = c("User_input","ensembl_gene_id")
		conversionTable$Species = sapply(result[,3], LogicManager$DB$findSpeciesNameById )


		return(	
			list(
				originalIDs = querySet,
				ensemblIDs=unique( result[,2]),			
				species = LogicManager$DB$findSpeciesById(result$species[1]), 
				speciesMatched = speciesMatched,
				conversionTable = conversionTable
			)
		)
	}
)

#		sub function for convertID
PreProcessing.Logic$set("public", "ConvertIDAutoMatch",
	function(result){
		comb = paste(result$species,result$idType)
		sortedCounts = sort(table(comb),decreasing=T)
		recognized =names(sortedCounts[1])
		result <- result[which(comb == recognized),]
		speciesMatched=sortedCounts
		names(speciesMatched )= sapply(as.numeric(gsub(" .*","",names(sortedCounts) ) ), LogicManager$DB$findSpeciesNameById  ) 
		speciesMatched <- as.data.frame( speciesMatched )

		if(length(sortedCounts) == 1) { # if only  one species matched
			speciesMatched[1,1] <-paste( rownames(speciesMatched), "(",speciesMatched[1,1],")",sep="")
		} else {# if more than one species matched
			speciesMatched[,1] <- as.character(speciesMatched[,1])
			speciesMatched[,1] <- paste( speciesMatched[,1]," (",speciesMatched[,2], ")", sep="") 
			speciesMatched[1,1] <- paste( speciesMatched[1,1],"   ***Used in mapping***  To change, select from above and resubmit query.") 	
			speciesMatched <- as.data.frame(speciesMatched[,1])
		}

		return(speciesMatched)
	}	
)


# Apply convert ID result to given data
# this data can be raw read count or transformed read count data

#	When given data is raw read count:
	# refer to convertedCounts() in 0.81 code
	# Convert Raw read count data based on Convert ID result.
	# If no conversion applied on ID, then use raw read count directly

#	When given data is transformed data:
	# refer to convertedData() in 0.81 code
	# Convert Transformed data based on Convert ID result. 
	# If no conversion applied on ID, then use transformed data directly

#	When given data is Pval:
	# refer to convertedPvals() in 0.81 code

PreProcessing.Logic$set("public", "ApplyConvertIDToGivenData",
	function(GivenData, ConversionTable){
		
		rownames(GivenData) = toupper(rownames(GivenData))

		# any gene not recognized by the database is disregarded
		# the 3 lines keeps the unrecogized genes using original IDs
		dat = merge(ConversionTable[,1:2], GivenData,  by.y = 'row.names', by.x = 'User_input', all.y=TRUE)

		# original IDs used if ID is not matched in database
		ix = which(is.na(dat[,2]) )
		dat[ix,2] = dat[ix,1] 	

		#multiple matched IDs, use the one with highest SD
		tem = apply(dat[,3:(dim(dat)[2])],1,sd)
		dat = dat[order(dat[,2],-tem),]
		dat = dat[!duplicated(dat[,2]) ,]
		rownames(dat) = dat[,2]
		dat = as.matrix(dat[,c(-1,-2)])
		tem = apply(dat,1,sd)
		dat = dat[order(-tem),]  # sort again by SD

		return(dat)
	}
)


# convert sorted species:idType combs into a list for repopulate species choice
PreProcessing.Logic$set("public", "matchedSpeciesInfo",
	function(x) {
  		a<- c()
  		for( i in 1:length(x)) {
  		  	a = c(a,paste( gsub("genes.*","",LogicManager$DB$findSpeciesNameById( as.numeric(gsub(" .*","",names(x[i])) ))), " (",
  		  	               x[i]," mapped from ",findIDtypeById( gsub(".* ","",names(x[i]) ) ),")",sep="") 
  		  	) 
		}      
  		return(a)
	}
)


# Clean up gene sets. Remove spaces and other control characters from gene names  
PreProcessing.Logic$set("public", "cleanGeneSet",
	function(x){
  		# remove duplicate; upper case; remove special characters
  		x <- unique( toupper( gsub("\n| ","",x) ) )
  		x <- x[which( nchar(x)>1) ]  # genes should have at least two characters
  		return(x)
	}
)

# allGeneInfo()
PreProcessing.Logic$set("public", "GetGenesInfomationByEnsemblIDs",
	function(ensemblIDs, species, selectOrg){
		if(selectOrg != 'BestMatch'){
			idxGeneInfoFile = grep(LogicManager$DB$findSpeciesById(selectOrg)[1,1], geneInfoFiles )
		}else{
			idxGeneInfoFile = grep(species[1,1],CONFIG_DATA_LIST_GENEINFOFILES)
		}
		
		if(length(idxGeneInfoFile) == 0 ) {
			return(as.data.frame("No matching gene info file found") ) 
		}


		if(length(idxGeneInfoFile) == 1){ 
			# if only one file           
			# read in the chosen file 
			geneInfo = read.csv(as.character(CONFIG_DATA_LIST_GENEINFOFILES[idxGeneInfoFile]) )
			geneInfo[,1]= toupper(geneInfo[,1]) #WBGene0000001 some ensembl gene ids in lower case
		} else { 
			return(as.data.frame("Multiple geneInfo file found!") )   
		}

		Set = match(geneInfo$ensembl_gene_id, ensemblIDs)
		Set[which(is.na(Set))]="Genome"
		Set[which(Set!="Genome")] ="List"

		return( cbind(geneInfo,Set) )
	}
)

# GetGeneSymbolUsingEnsembl 
PreProcessing.Logic$set("public", "GetGenesSymbolByEnsemblIDs",
	function(allGeneInfo, EnsemblIDs, isUseEnsemblIfNoSymbol){
		if(length(EnsemblIDs) == 0){
			return(NULL)
		}

		tb.QueryData <- as.data.frame(EnsemblIDs) 
		colnames(tb.QueryData) <- c('ensembl_gene_id')

		tb.result <- tb.QueryData %>% 
			left_join(allGeneInfo, by='ensembl_gene_id') %>%
			select('ensembl_gene_id', 'symbol')	
		
		tb.result$symbol <- as.character(tb.result$symbol)
		tb.result$ensembl_gene_id <- as.character(tb.result$ensembl_gene_id)
		if(isUseEnsemblIfNoSymbol){
			# if we want use ensembl id for the symbol missed gene
			symbol <- ifelse(is.na(tb.result$symbol), tb.result$ensembl_gene_id, tb.result$symbol)
		}else{
			# if not, just pull the symbol out
			symbol <- tb.result$symbol
		}
		
		return(symbol)
	}
)

# GenerateDataForSingleGenePlot
PreProcessing.Logic$set("public", "GenerateDataForSingleGenePlot",
	function(ConvertedTransformedData, AllGeneInfo, SelectedOrg, GeneID){
		Symbols <- rownames(ConvertedTransformedData)

		if( SelectedOrg != "NEW" &&  ncol(AllGeneInfo) != 1 ) {
			ix = match( rownames(ConvertedTransformedData), AllGeneInfo[,1])
			if( sum( is.na(AllGeneInfo$symbol) ) != dim(AllGeneInfo)[1] ) {  
				# symbol really exists?
				Symbols = as.character( AllGeneInfo$symbol[ix] )
				Symbols[which( nchar(Symbols) <= 2 ) ] <- rownames(ConvertedTransformedData) [which( nchar(Symbols) <= 2 )]
			}
		}

		dat <- as.data.frame(ConvertedTransformedData)

		dat$Genes = Symbols

		searchWord = gsub("^ ","",GeneID )
		ix = which(regexpr(  paste("^" , toupper(searchWord),sep="")   ,toupper(dat$Genes)) > 0)

		if(grepl(" $", searchWord)){
			# if there is space character, do exact match
			ix = match(gsub(" ","", toupper(searchWord)), toupper(dat$Genes) )
		}  
		
		# too few or too many genes found
		if(length(ix) == 0 | length(ix) > 50 ){
			return(NULL)
		}

		mdf = melt(dat[ix,],id.vars="Genes", value.name="value", variable.name="samples")

		return(mdf)
	}
)

# GenerateDataForAllSamplesSingleGenePlot
PreProcessing.Logic$set("public", "GenerateDataForAllSamplesSingleGenePlot",
	function(mdf){

		# Barplot with error bars
		mdf$count = 1
		g = self$DetectGroups(mdf$samples)
		Means = aggregate(mdf$value,by=list( g, mdf$Genes ), FUN = mean, na.rm=TRUE  )
		SDs = aggregate(mdf$value,by=list( g, mdf$Genes ), FUN = sd, na.rm=TRUE  )
		Ns = aggregate(mdf$count, by= list(g, mdf$Genes) , FUN = sum  )
		summarized = cbind(Means,SDs[,3],Ns[,3])
		colnames(summarized)= c("Samples","Genes","Mean","SD","N")
		summarized$SE = summarized$SD / sqrt(summarized$N)			

		return(summarized)
	}
)

# Logic used to generate 'processedData' 
# Basically, this function will attach ensemble_id and gene_symbol onto the converted transformed data
PreProcessing.Logic$set("public", "FormatProcessedTransformedDataForDownload",
	function(AllGeneInfo, ConvertedTransformedData, ConvertedIDResult){
		dat <- merge(
			AllGeneInfo[,c('ensembl_gene_id','symbol')],
			round(ConvertedTransformedData, 4),
			by.x = "ensembl_gene_id",
			by.y = "row.names",
			all.y = TRUE
		)

		dat[, 2] =  paste(" ", dat[, 2]) # add space to gene symbol to avoid auto convertion of symbols to dates by Excel 
		dat <- merge(ConvertedIDResult$conversionTable, dat, 
			by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.y = TRUE)
		dat <- dat[,-3]						# remove species column
	
		ix <- which( dat[,3] == "  NA" )	# find NA index 
		dat[ix,3] <- ""   					# remove NA's 

		# Legacy code:
		# ix <- which(is.na(dat[,2]))			# record not in mapping table
		# uid <- dat[,2]
		# uid[ix] = dat[ix,1]
		# Above code equal to: 
		uids <- ifelse(is.na(dat[,2]), dat[,1], dat[,2]) 
		uids = paste0(" ",uids)	# prevents Excel auto conversion	
		dat[,2] <- ifelse(is.na(dat[,2]), "", dat[,1])
		dat[,1] <- uids

		colnames(dat)[1:2]= c("User_ID sorted by SD","Ensembl_gene_id")
		
		# sort by sd
		dat = dat[ order( -apply(dat[,-3:-1],1,sd )   )  ,]
		
		rownames(dat)=1:nrow(dat)
		
		return(dat)
	}
)

# Logic used to generate 'processedCountsData' 
# Basically, this function will attach ensemble_id and gene_symbol onto the converted read count data
PreProcessing.Logic$set("public", "FormatProcessedRawReadcountDataForDownload",
	function(AllGeneInfo, ConvertedRawReadcountData, ConvertedIDResult){
		dat <- merge(
			AllGeneInfo[,c('ensembl_gene_id','symbol')],
			round(ConvertedRawReadcountData, 4),
			by.x = "ensembl_gene_id",
			by.y = "row.names",
			all.y = TRUE
		)

		dat[, 2] =  paste(" ", dat[, 2]) # add space to gene symbol to avoid auto convertion of symbols to dates by Excel 
		dat <- merge(ConvertedIDResult$conversionTable, dat, 
			by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.y = TRUE)
		dat <- dat[,-3]						# remove species column
	
		ix <- which( dat[,3] == "  NA" )	# find NA index 
		dat[ix,3] <- ""   					# remove NA's 

		# Legacy code:
		# ix <- which(is.na(dat[,2]))			# record not in mapping table
		# uid <- dat[,2]
		# uid[ix] = dat[ix,1]
		# Above code equal to: 
		uids <- ifelse(is.na(dat[,2]), dat[,1], dat[,2]) 
		uids = paste0(" ",uids)	# prevents Excel auto conversion	
		dat[,2] <- ifelse(is.na(dat[,2]), "", dat[,1])
		dat[,1] <- uids

		colnames(dat)[1:2]= c("User_ID sorted by SD","Ensembl_gene_id")
		
		# sort by sd
		dat = dat[ order( -apply(log2(10+ dat[,-3:-1]),1,sd )   )  ,]

		rownames(dat)=1:nrow(dat)
		
		return(dat)
	}
)



