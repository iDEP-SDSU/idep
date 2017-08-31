## DataClass
# library("data.table")

DataClass <- setRefClass("DataClass",
  fields = list(
    data = "data.table",
    train = "data.table",
    test = "data.table"
  ),
  methods = list(
    initialize = function(filePath){
      print("DataClass is initiated")
      demoDataFile = filePath
      .self$data = data.table(read.csv(demoDataFile))
    },
    finalize = function(){
      print("DataClass is terminated")
    }
  )
)
