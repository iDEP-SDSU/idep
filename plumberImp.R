
source('classes/librarySetup.R', echo=T)
source('classes/DataClass.R', echo=T)

#* @get /mean
normalMean <- function(samples=10){
  data <- rnorm(samples)
  mean(data)
}


#* Get a graph of the values
#* @get /plots
plots <- function(){
  # load data
  dataObj <- DataClass$new(filePath="GSE37704_sailfish_genecounts.csv")
  localData = dataObj$data
  #colsToSum <- names(localData[,-1])
  dt2 = localData[, lapply(.SD, sum, na.rm=TRUE), .SDcols=colsToSum]
  #barplot( as.matrix(dt2/1e6), col="green",las=3, main="Total read counts (millions)")
  #plot(1:10, 1:10, main="test")
  
  list(result =as.matrix(dt2/1e6))
}
