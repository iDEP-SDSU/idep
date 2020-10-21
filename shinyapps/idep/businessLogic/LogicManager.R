library('R6')
source('server.config')

source('businessLogic/UtilityFunctions.R')
source('businessLogic/Files.R') ## File Manager
source('businessLogic/Display.R') ## Plot and Image Manager
source('businessLogic/PreProcessing.R')		## Preprocessing logic
source('businessLogic/DB.R')		## database manager
source('businessLogic/Heatmap.R')
source('businessLogic/Kmeans.R')
source('businessLogic/PCA.R')

Logic.Manager <- R6Class("Logic.Manager")


Logic.Manager$set("public", "Files", "")
Logic.Manager$set("public", "Display", "")
Logic.Manager$set("public", "PreProcessing", "")
Logic.Manager$set("public", "DB", "")
Logic.Manager$set("public", "UtilFuns", "")
Logic.Manager$set("public", "Heatmap", "")
Logic.Manager$set("public", "Kmeans", "")
Logic.Manager$set("public", "PCA", "")




Logic.Manager$set("public", "initialize",
	function(){
		self$UtilFuns <- UtilFuns$new()
		self$Files <- File.Manager$new()
		self$Display <- Display.Manager$new()
		self$DB <- DB.Manager$new()
		self$PreProcessing <- PreProcessing.Logic$new()
		self$Heatmap <- Heatmap.Logic$new()
        self$Kmeans <- Kmeans.Logic$new()
        self$PCA <- PCA.Logic$new()
	}
)

