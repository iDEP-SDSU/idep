# Migrate SQLite to PostgreSQL
## 1. sqlite3 convertIDs.db .dump > d96.db.bak
## 2. Manually edit data types
## 3. Import .bak file on PostgreSQL

# Installing drivers [on macOS]
## 1. brew install unixodbc
## 2. brew install psqlodbc
## 3. brew services start postgresql

# Create connection
con <- DBI::dbConnect(odbc::odbc(),
                      driver   = "/usr/local/lib/psqlodbcw.so",
                      database = "postgres",
                      UID      = "postgres", #rstudioapi::askForPassword("Database user")
                      PWD      = rstudioapi::askForPassword("Database password"),
                      host     = "localhost", 
                      Port     = 5432)


# References
## [1] Drivers: https://db.rstudio.com/best-practices/drivers/