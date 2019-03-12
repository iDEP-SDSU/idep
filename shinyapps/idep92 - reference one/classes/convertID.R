library(dplyr)
library(RPostgreSQL)
Postgres()

ConvertDB <- setRefClass("ConvertDB", 
  fields = list(dataSet = "ANY", cnn = "ANY", level2Terms = "ANY", idIndex="ANY", quotes="ANY", speciesChoice="ANY"),
  methods = list(
    ## Constructor 
    initialize = function(){
      sqlite  <- dbDriver("SQLite")
      .self$cnn = dbConnect(PostgreSQL(),
        host = "192.241.138.85", port = "5432",
        user = "idep", password = "sdsu57007", dbname="idep")
      .self$go_levels()
      .self$proc_idIndex()
      .self$proc_orgInfo()
      print("all init function done - ConvertDB")
    },
    go_levels = function(){
      query <- "select distinct id,level from GO WHERE GO = 'biological_process'"
      GO_levels = dbGetQuery(.self$cnn, query)
      .self$level2Terms = GO_levels %>% 
        filter(level %in% c(2,3)) %>% 
        select(id)
      return(GO_levels)
    },
    proc_idIndex = function(){
      queryidIndex = "select distinct * from idIndex"
      queryQuotes = "select * from quotes"
      
      .self$idIndex = dbGetQuery(.self$cnn, queryidIndex)
      quotesTemp <- dbGetQuery(.self$cnn, queryQuotes)
      
      .self$quotes = paste0("\"",quotesTemp$quotes,"\"", " -- ",quotesTemp$author,".       ")
    },
    proc_orgInfo = function() {
      # prepare species list
      # Create a list for Select Input options
      orgInfo = dbGetQuery(.self$cnn, paste("select distinct * from orgInfo " ))
      orgInfo = orgInfo[order(orgInfo$name),]
      speciesChoiceTemp = setNames(as.list( orgInfo$id ), orgInfo$name2 )
      # add a defult element to list    # new element name       value
      speciesChoiceTemp = append( setNames( "NEW","**NEW SPECIES**"), speciesChoiceTemp  )
      .self$speciesChoice = append( setNames( "BestMatch","Best matching species"), speciesChoiceTemp  )
      # move one element to the 2nd place
      # move2 <- function(i) c(speciesChoice[1:2],speciesChoice[i],speciesChoice[-c(1,i)])
      # i= grep("Glycine max" ,names(speciesChoice)); speciesChoice <- move2(i)
      # i= grep("Zea mays" ,names(speciesChoice)); speciesChoice <- move2(i)
      # i= grep("Arabidopsis thaliana",names(speciesChoice)); speciesChoice <- move2(i)
      # i= grep("Saccharomyces cerevisiae" ,names(speciesChoice)); speciesChoice <- move2(i)
      # i= grep("Caenorhabditis elegans",names(speciesChoice)); speciesChoice <- move2(i)
      # i= grep("Danio rerio" ,names(speciesChoice)); speciesChoice <- move2(i)
      # i= grep("Bos taurus" ,names(speciesChoice)); speciesChoice <- move2(i)
      # i= grep("Rattus norvegicus" ,names(speciesChoice)); speciesChoice <- move2(i)
      # i= grep("Mus musculus",names(speciesChoice)); speciesChoice <- move2(i)
      # i= grep("Homo sapiens",names(speciesChoice)); speciesChoice <- move2(i)
    }
  )
)

