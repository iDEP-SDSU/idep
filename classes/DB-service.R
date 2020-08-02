# Migrate SQLite to PostgreSQL
## 1. sqlite3 convertIDs.db .dump > d96.db.bak
## 2. Manually edit data types
## 3. Import .bak file on PostgreSQL
# (Alternative)
## 1. brew install --HEAD pgloader
## 2. createdb DB
## 3. pgloader sqliteDB.db postgresql://USER:PASSWORD@HOST:PORT/DB

# Installing drivers [on macOS]
## 1. brew install unixodbc
## 2. brew install psqlodbc
## 3. brew services start postgresql

# pgloader --root-dir /scratch/opt/tmp convertIDs.db postgresql://idep:iDEP666@localhost/idep96
# Create connection
con <- DBI::dbConnect(odbc::odbc(),
                      driver   = "/usr/local/lib/psqlodbcw.so",
                      database = "idep96",
                      UID      = "idep", #rstudioapi::askForPassword("Database user")
                      PWD      = rstudioapi::askForPassword("Database password"),
                      server     = "pt.jacks.local", 
                      Port     = 5432)

convert <- con
quotes <- DBI::dbGetQuery(convert, " select * from quotes")
quotes <- paste0("\"", quotes$quotes, "\"", " -- ", quotes$author,".       ")

# Create a list for Select Input options
orgInfo <- DBI::dbGetQuery(convert, paste("select distinct * from orgInfo " ))
orgInfo <- orgInfo[order(orgInfo$name), ]
annotatedSpeciesCounts <- sort(table(orgInfo$group)) # total species, Ensembl, Plants, Metazoa, STRINGv10
speciesChoice <- setNames(as.list(orgInfo$id), orgInfo$name2 )
# add a defult element to list    # new element name       value

speciesChoice <- append(setNames( "NEW","**NEW SPECIES**"), speciesChoice)
speciesChoice <- append(setNames( "BestMatch","Best matching species"), speciesChoice)

# move one element to the 2nd place
move2 <- function(i) c(speciesChoice[1:2], speciesChoice[i], speciesChoice[-c(1, 2, i)])
i <- which(names(speciesChoice) == "Vitis vinifera"); speciesChoice <- move2(i)
i <- which(names(speciesChoice) == "Oryza sativa Japonica Group"); speciesChoice <- move2(i)
i <- which(names(speciesChoice) == "Oryza sativa Indica Group"); speciesChoice <- move2(i)
i <- which(names(speciesChoice) == "Glycine max"); speciesChoice <- move2(i)
i <- which(names(speciesChoice) =="Zea mays"); speciesChoice <- move2(i)
i <- which(names(speciesChoice) =="Arabidopsis thaliana"); speciesChoice <- move2(i)
i <- which(names(speciesChoice) == "Saccharomyces cerevisiae"); speciesChoice <- move2(i)
i <- which(names(speciesChoice)  == "Caenorhabditis elegans"); speciesChoice <- move2(i)
i <- which(names(speciesChoice)  == "Drosophila melanogaster"); speciesChoice <- move2(i)
i <- which(names(speciesChoice) =="Dog"); speciesChoice <- move2(i)
i <- which(names(speciesChoice) =="Macaque"); speciesChoice <- move2(i)
i <- which(names(speciesChoice) =="Chicken"); speciesChoice <- move2(i)
i <- which(names(speciesChoice) =="Pig"); speciesChoice <- move2(i)
i <- which(names(speciesChoice) =="Zebrafish" ); speciesChoice <- move2(i)
i <- which(names(speciesChoice) == "Cow" ); speciesChoice <- move2(i)
i <- which(names(speciesChoice) == "Rat" ); speciesChoice <- move2(i)
i <- which(names(speciesChoice) == "Mouse"); speciesChoice <- move2(i)
i <- which(names(speciesChoice) == "Human"); speciesChoice <- move2(i)

GO_levels <- DBI::dbGetQuery(convert, "select distinct id,level from GO  
                                WHERE GO = 'biological_process'")
level2Terms <- GO_levels[which(GO_levels$level %in% c(2,3)), 1]  # level 2 and 3
idIndex <- DBI::dbGetQuery(convert, paste("select distinct * from idIndex "))
quotes <- DBI::dbGetQuery(convert, " select * from quotes")
quotes = paste0("\"", quotes$quotes, "\"", " -- ", quotes$author, ".       ")

# Clean up gene sets. Remove spaces and other control characters from gene names  
cleanGeneSet <- function (x) {
  # remove duplicate; upper case; remove special characters
  x <- unique(toupper(gsub("\n| ", "", x)))
  x <- x[which(nchar(x) > 1)]  # genes should have at least two characters
  return(x)
}

#find idType based on index 
findIDtypeById <- function(x){ # find 
  return( idIndex$idType[ as.numeric(x)] )
}
findSpeciesById <- function (speciesID){ # find species name use id
  return( orgInfo[which(orgInfo$id == speciesID),]  )
}
# just return name
findSpeciesByIdName <- function (speciesID){ # find species name use id
  return( orgInfo[which(orgInfo$id == speciesID),3]  )
}

# convert gene IDs to ensembl gene ids and find species
convertID <- function (query, selectOrg, selectGO) {
  querySet <- cleanGeneSet(unlist(strsplit(toupper(query), '\t| |\n|\\,')))
  # querySet is ensgene data for example, ENSG00000198888, ENSG00000198763, ENSG00000198804
  
  if(selectOrg == "BestMatch") { # query all species
    querySTMT <- paste0("select distinct id,ens,species from mapping where id IN ('", paste(querySet, collapse = "', '"), "')")
  } else {  # organism has been selected query specific one
    querySTMT <- paste0("select distinct id,ens,species from mapping where species = '", selectOrg,
                        "' AND id IN ('", paste(querySet, collapse="', '"), "')")
  }
  tictoc::tic("Big query [Convert Gene IDs]")
  result <- DBI::dbGetQuery(convert, querySTMT)
  tictoc::toc()
  
  if(dim(result)[1] == 0) return(NULL)
  if(selectOrg == speciesChoice[[1]]) {
    comb <- paste(result$species, result$idType)
    sortedCounts <- sort(table(comb), decreasing = TRUE)
    # Try to use Ensembl instead of STRING-db genome annotation
    if(sortedCounts[1] <= sortedCounts[2] * 1.1  # if the #1 species and #2 are close
        && as.numeric(names(sortedCounts[1])) > sum(annotatedSpeciesCounts[1:3])  # 1:3 are Ensembl species
        && as.numeric(names(sortedCounts[2])) < sum(annotatedSpeciesCounts[1:3])) { # and #2 come earlier (ensembl) than #1
      tem <- sortedCounts[2]
      sortedCounts[2] <- sortedCounts[1]
      names(sortedCounts)[2] <- names(sortedCounts)[1]
      sortedCounts[1] <- tem
      names(sortedCounts)[1] <- names(tem)    
    } 
    recognized <- names(sortedCounts[1])
    result <- result[which(comb == recognized),]
    speciesMatched <- sortedCounts
    names(speciesMatched) <- sapply(as.numeric(gsub(" .*", "", names(sortedCounts))), findSpeciesByIdName)
    speciesMatched <- as.data.frame(speciesMatched)
    if(length(sortedCounts) == 1) { # if only  one species matched
      speciesMatched[1,1] <- paste0(rownames(speciesMatched), "(", speciesMatched[1, 1], ")")
    } else {# if more than one species matched
      speciesMatched[,1] <- as.character(speciesMatched[,1])
      speciesMatched[,1] <- paste0(speciesMatched[, 1], " (", speciesMatched[, 2], ")") 
      speciesMatched[1,1] <- paste(speciesMatched[1, 1], "   ***Used in mapping***  To change, select from above and resubmit query.") 	
      speciesMatched <- as.data.frame(speciesMatched[, 1])
    }
  } else { # if species is selected
    result <- result[which(result$species == selectOrg), ]
    if(dim(result)[1] == 0) return(NULL) #stop("ID not recognized!")
    speciesMatched <- as.data.frame(paste("Using selected species ", findSpeciesByIdName(selectOrg)))
  }
  result <- result[which(!duplicated(result[,2]) ),] # remove duplicates in ensembl_gene_id
  result <- result[which(!duplicated(result[,1]) ),] # remove duplicates in user ID
  colnames(speciesMatched) = c("Matched Species (genes)" ) 
  conversionTable <- result[,1:2]; colnames(conversionTable) = c("User_input","ensembl_gene_id")
  conversionTable$Species = sapply(result[,3], findSpeciesByIdName )
  if(0){
    # generate a list of gene set categories
    ix = grep(findSpeciesById(result$species[1])[1,1],gmtFiles)
    if (length(ix) == 0 ) {categoryChoices = NULL}
    # If selected species is not the default "bestMatch", use that species directly
    if(selectOrg != speciesChoice[[1]]) {  
      ix = grep(findSpeciesById(selectOrg)[1,1], gmtFiles )
      if (length(ix) == 0 ) {categoryChoices = NULL}
      totalGenes <- orgInfo[which(orgInfo$id == as.numeric(selectOrg)),7]
    }
    pathway <- DBI::dbConnect(sqlite,gmtFiles[ix],flags=SQLITE_RO)
    # Generate a list of geneset categories such as "GOBP", "KEGG" from file
    geneSetCategory <-  DBI::dbGetQuery(pathway, "select distinct * from categories " ) 
    geneSetCategory  <- geneSetCategory[,1]
    categoryChoices <- setNames(as.list( geneSetCategory ), geneSetCategory )
    categoryChoices <- append( setNames( "All","All available gene sets"), categoryChoices  )
    #change GOBO to the full description for display
    names(categoryChoices)[ match("GOBP",categoryChoices)  ] <- "GO Biological Process"
    names(categoryChoices)[ match("GOCC",categoryChoices)  ] <- "GO Cellular Component"
    names(categoryChoices)[ match("GOMF",categoryChoices)  ] <- "GO Molecular Function"
    dbDisconnect(pathway)
  } #if (0)
  
  return(list(originalIDs = querySet,IDs=unique( result[,2]), 
              species = findSpeciesById(result$species[1]), 
              #idType = findIDtypeById(result$idType[1] ),
              speciesMatched = speciesMatched,
              conversionTable = conversionTable
  ))
}
# References
## [1] Drivers: https://db.rstudio.com/best-practices/drivers/