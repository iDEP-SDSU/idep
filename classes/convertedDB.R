# Data Read 
# sqlite3 -header -csv -separator ',' 'select * from orgInfo;' > orgInfo.csv

convert <- dbConnect(sqlite,"./shinyapps/idep/data/convertIDs.db")
table <- dbGetQuery(convert, paste("SELECT name FROM sqlite_temp_master WHERE type='table'" ))
allTables <- dbListTables(convert)
GO <- dbGetQuery(convert, paste("SELECT * FROM ",allTables[1], sep="" ))
mapping <- dbGetQuery(convert, paste("SELECT * FROM ",allTables[2], sep="" ))
orgInfo <- dbGetQuery(convert, paste("SELECT * FROM ",allTables[3], sep="" ))
quotes <- dbGetQuery(convert, paste("SELECT * FROM ",allTables[4], sep="" ))


require("RPostgreSQL")

# loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")

pw <- {"example"}

install.packages("RPostgreSQL")
require("RPostgreSQL")
library(DBI)
sqlite  <- dbDriver("SQLite")
# creates a connection to the postgres database
# note that "con" will be used later in each connection to the database
con <- dbConnect(PostgreSQL(),
  host = "192.241.138.85", port = "5432",
  user = "idep", password = "sdsu57007", dbname="idep")
dtab = dbGetQuery(con, "select * from mapping limit 50")
dtab

