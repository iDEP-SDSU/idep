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

# References
## [1] Drivers: https://db.rstudio.com/best-practices/drivers/