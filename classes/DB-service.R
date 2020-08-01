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


# References
## [1] Drivers: https://db.rstudio.com/best-practices/drivers/