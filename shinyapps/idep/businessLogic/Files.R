library('R6')
source('server.config')

File.Manager <- R6Class("File.Manager")


## Public Methods:
File.Manager$set("public", "LoadDemoData", 
	function(){
		rawData <- read.csv(CONFIG_DATA_DEMO_READCOUNT)
		if(dim(rawData)[2] <= 2 ){
			# if less than 3 columns, try tab-deliminated
			rawData <- read.table(CONFIG_DATA_DEMO_READCOUNT, sep="\t",header=TRUE,quote = "",comment.char="")
		}   

		rawSampleInfo <- read.csv(CONFIG_DATA_DEMO_SAMPLEINFO, row.names=1,header=T,colClasses="character")
		if(dim(rawSampleInfo)[2] <= 2 ){
			# if less than 3 columns, try tab-deliminated
			rawSampleInfo <- read.table(CONFIG_DATA_DEMO_SAMPLEINFO, row.names=1,sep="\t",header=TRUE,colClasses="character")
		}

		result <- list(dat=rawData, sampleInfo=rawSampleInfo)
		return(result)
	}
)

File.Manager$set("public", "LoadUploadedData", 
	function(dataUrl, sampleInfoUrl){
		rawData <- read.csv(dataUrl)
		if(dim(rawData)[2] <= 2 ){
			# if less than 3 columns, try tab-deliminated
			rawData <- read.table(dataUrl, sep="\t",header=TRUE,quote = "",comment.char="")
		}   

		if(is.null(sampleInfoUrl)){
			rawSampleInfo <- NULL
		}else {
		   	rawSampleInfo <- read.csv(sampleInfoUrl, row.names=1,header=T,colClasses="character")
			if(dim(rawData)[2] <= 2 ){
				# if less than 3 columns, try tab-deliminated
				rawData <- read.table(sampleInfoUrl, row.names=1,sep="\t",header=TRUE,colClasses="character")
			}  
		}
		
		result <- list(dat=rawData, sampleInfo=rawSampleInfo)
		return(result)
	}
)















## Private Methods:







