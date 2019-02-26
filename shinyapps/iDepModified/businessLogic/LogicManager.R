library('R6')
source('server.config')

source('businessLogic/Files.R') ## File Manager
source('businessLogic/PreProcessing.R')		## Preprocessing logic

Logic.Manager <- R6Class("Logic.Manager")


Logic.Manager$set("public", "Files", "")
Logic.Manager$set("public", "PreProcessing", "")



Logic.Manager$set("public", "initialize",
	function(){
		self$Files <- File.Manager$new()
		self$PreProcessing <- PreProcessing.Logic$new()
	}
)

