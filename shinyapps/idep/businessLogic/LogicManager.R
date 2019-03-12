library('R6')
source('server.config')

source('businessLogic/Files.R') ## File Manager
source('businessLogic/Display.R') ## Plot and Image Manager
source('businessLogic/PreProcessing.R')		## Preprocessing logic
source('businessLogic/DB.R')		## database manager

Logic.Manager <- R6Class("Logic.Manager")


Logic.Manager$set("public", "Files", "")
Logic.Manager$set("public", "Display", "")
Logic.Manager$set("public", "PreProcessing", "")
Logic.Manager$set("public", "DB", "")



Logic.Manager$set("public", "initialize",
	function(){
		self$Files <- File.Manager$new()
		self$PreProcessing <- PreProcessing.Logic$new()
		self$Display <- Display.Manager$new()
		self$DB <- DB.Manager$new()
	}
)

