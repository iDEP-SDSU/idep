####################################################
# Author: Eric Tulowetzke
# Lab: Ge Lab
# R version 4.0.5 (2020-10-10) -- "Bunny-Wunnies Freak Out"
# Project: ShinyGO Restructure
# Purpose of project: rewrite ShinyGO61 code to be easier to develop
# File: allsamplePreprocess.R
# Purpose of file: Preprocess data from SQLite database, into a feather file
#  where user can easily find example of gene IDs when species picked
# Start data: 02-02-2021 (mm-dd-yyyy)
# Data last modified: (mm-dd-yyyy)
#######################################################

Rlibs <- c('RSQLite', 'feather')
notInstalled <- setdiff(Rlibs, rownames(installed.packages()))
if(length(notInstalled)>0) {install.packages(notInstalled)}
lapply(Rlibs, require, character.only = TRUE)

path <- 'convertIDs.db'
filename <- 'example_of_id.feather'
convertIDsDatabase <- dbConnect(RSQLite::SQLite(), path)
specie <- dbGetQuery(convertIDsDatabase, paste('SELECT * FROM orginfo'))
specieList <- unique(c("Human", sort(specie$name2)))
RSQLite::dbDisconnect(convertIDsDatabase)


work <- function(sp, done, total) {
  file <- paste(sp, done, total, collapse = '-')
  file <- paste(file, '.txt')
  m <- matrix(c(2,2, 2,2), 2, 2)
  convertIDsDatabase <- dbConnect(RSQLite::SQLite(), path)
  query4Specie <- paste('SELECT * FROM orginfo WHERE name2 =', shQuote(sp))
  specie <- dbGetQuery(convertIDsDatabase, query4Specie)
  
  
  if (specie$totalGenes == 0) {
    matchIDtypeDf <- data.frame('index' = sp,
                                'id' = '0 Genes found',
                                'gene' = '0 Genes found')
  } else {
    query4IDmap <- paste('SELECT * FROM mapping WHERE species =', as.numeric(specie$id))
    foundGenesDf <- dbGetQuery(convertIDsDatabase, query4IDmap)
    idTypes <- unique(foundGenesDf$idType)
    idTypes4SQL <- paste(idTypes, collapse = ", ")
    query4idTypeMatch <- paste('SELECT * FROM idIndex WHERE id IN (', idTypes4SQL, ')')
    idIndexDf <- dbGetQuery(convertIDsDatabase, query4idTypeMatch)
    
    matchIDtypeDf <- data.frame('index' = vector(length = length(idIndexDf$idType)),
                                'id' = idIndexDf$idType)
    for (indexR in 1:length(idIndexDf$idType)) {
      matchIDtypeDf$index[indexR] <- sp
      tmpDf <- foundGenesDf[foundGenesDf$idType == idTypes[indexR],]
      tmpVec <- c(as.character(tmpDf$id[1]))
      for (index in 2:11) {
        tmpVec <- c(tmpVec, tmpDf$id[index])
      } # end of inner for
      tmpVec <- paste(tmpVec, collapse = ", ")
      matchIDtypeDf$gene[indexR] <- tmpVec
    } # end of outer for
  }
  
  RSQLite::dbDisconnect(convertIDsDatabase)
  write.table(m, file)
  return(matchIDtypeDf)
}

allSample <- data.frame('index' = vector(), 'id' = vector(),
                        'gene' = vector())
number <- 1
countS <- length(specieList)
for (s in specieList) {
  tmp <- work(sp = s, done = number, total = countS)
  allSample <- rbind(allSample, tmp)
  number <- number + 1
}
feather::write_feather(allSample, filename)





