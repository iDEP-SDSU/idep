library('R6')
source('server.config')

File.Manager <- R6Class("File.Manager")


## Public Methods:
File.Manager$set("public", "LoadDemoData", 
	function(){
		rawData <- read.csv(CONFIG_DATA_DEMO_READCOUNT)
		rawTestDesign <- read.csv(CONFIG_DATA_DEMO_SAMPLEINFO)
		result <- list(dat=rawData, design=rawTestDesign)
		return(result)
	}
)

File.Manager$set("public", "LoadUploadedData", 
	function(dataUrl, designUrl){
		rawData <- read.csv(dataUrl)

		if(is.null(designUrl)){
			rawTestDesign <- NULL
		}else {
		   	rawTestDesign <- read.csv(designUrl)
		}
		
		result <- list(dat=rawData, design=rawTestDesign)
		return(result)
	}
)















## Private Methods:







