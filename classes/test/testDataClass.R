source('classes/librarySetup.R', echo=F)
source('classes/DataClass.R', echo=F)

pdf(NULL) # this prevents error Cannot open file 'Rplots.pdf'
Min_overlap <- 2
minSetSize = 3;
mappingCoverage = 0.60 # 60% percent genes has to be mapped for confident mapping
mappingEdge = 0.5  # Top species has 50% more genes mapped
PvalGeneInfo = 0.05; minGenes = 10 # min number of genes for ploting
kurtosis.log = 50  # log transform is enforced when kurtosis is big
kurtosis.warning = 10 # log transformation recommnded
minGenesEnrichment = 2 # perform GO or promoter analysis only if more than this many genes
PREDA_Permutations =1000

#sqlite  <- dbDriver("SQLite")
#convert <- dbConnect(sqlite,"../go/convertIDs.db")
#set.seed(2)
mycolors = sort(rainbow(20))[c(1,20,10,11,2,19,3,12,4,13,5,14,6,15,7,16,8,17,9,18)] # 20 colors for kNN clusters
#keggSpeciesID = read.csv("KEGG_Species_ID.csv")
#dim(keggSpeciesID)
#str(keggSpeciesID)

# load data
dataObj <- DataClass$new(filePath="GSE37704_sailfish_genecounts.csv")
dataObj$data

# Pre-Process
colsToSum <- names(dataObj$data[,-1])
dt2 = dataObj$data[, lapply(.SD, sum, na.rm=TRUE), .SDcols=colsToSum]
barplot( as.matrix(dt2/1e6), col="green",las=3, main="Total read counts (millions)")

## EDA -  required data transform 
# Distribution of transformed Data 
## TODO: transform data 

colors = dim(dataObj$data)[2]
myColors = rainbow(colors)
startedLog = 2
x = log(dataObj$data[,-1]+abs(startedLog),2)
x = as.data.frame(x)
tem <- apply(x,1,sd)
x <- x[order(-tem),]  # sort by SD


par(mfrow=c(3,1))
plot(density(x[,1]),col = myColors[1],xlab="Expresson values", ylab="Density", main= "Distribution of transformed data")
for( i in 2:dim(x)[2] ){
  lines(density(as.matrix(x[,i])),col=myColors[i] )
}
png("r.png", width=10, height=5)
boxplot(x, las = 2, ylab="Transformed expression levels", main="Distribution of transformed data")
dev.off()
plot(x[,1:2],xlab=colnames(x)[1],ylab=colnames(x)[2], main="Scatter plot of first two samples")



