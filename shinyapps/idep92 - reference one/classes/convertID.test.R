source('./classes/convertID.R')

cvObj = ConvertDB$new()
cvObj$speciesChoice

library(RPostgreSQL)
cnn = dbConnect(PostgreSQL(),
    host = "192.241.138.85", port = "5432",
    user = "idep", password = "sdsu57007", dbname="idep")

sqlInterpolate
sql <- "SELECT * FROM X WHERE name = ?name"
sqlInterpolate(ANSI(), sql, name = "Hadley")
sqlInterpolate(ANSI(), sql, name = "H'); DROP TABLE--;")
postgresqlExecStatement

res <- dbGetQuery(cnn, "SELECT * FROM idIndex WHERE 'idType' = 'ensembl_exon_id'")
str(res)
head(res)

src_postgres
install.packages("dbplyr")
library(dplyr)
require(dbplyr)
localdb <- src_postgres(host = "192.241.138.85", port = "5432",
  user = "idep", password = "sdsu57007", dbname="idep")
devtools::install_github("rstudio/pool")

library(pool)
library(DBI)
library(RPostgreSQL)
install.packages("RPostgres")
pool <- dbPool(
  drv = RPostgreSQL::PostgreSQL(max.con=20),
  dbname = "idep",
  host = "192.241.138.85",
  user = "idep", password = "sdsu57007"
)

dbGetInfo(pool)
result <- pool %>% tbl('KEGG_Species_ID') %>% filter(start_position < 2000)
sp <- result %>% select(start_position) %>% collect()
nrow(sp)

#https://shiny.rstudio.com/articles/sql-injections.html
sql <- "SELECT * FROM City WHERE ID = ?id1 OR ID = ?id2 OR ID = ?id3;"
query <- sqlInterpolate(conn, sql, id1 = "123", id2 = "bad123", id3 = "dfd")
dbGetQuery(conn, query)

query1 <- "SELECT kegg, array_to_string(array_agg(distinct KEGG_Species_ID.name) OVER (PARTITION BY kegg) FROM KEGG_Species_ID"
res <- dbGetQuery(cnn, query1)
res

query1 <- "SELECT distinct level, count(*) OVER(PARTITION BY level) as num_level FROM go order by level"
res <- dbGetQuery(cnn, query1)
barplot(t(as.matrix(res)), main ="Number of level records")

query2 <- "SELECT count(*) OVER (PARTITION BY idType) FROM mapping LIMIT 50"
result3
query3 <- "select t1.\"GO\", count(*) as cnt from go as t1 group by t1.\"GO\""
result3 <- dbGetQuery(cnn, query3)
barplot(t(as.matrix(result3$cnt)), main ="Number of level records for GO")

query4 <- "select t1.kegg, count(*) as cnt from \"KEGG_Species_ID\" as t1 group by t1.kegg"
result4 <- dbGetQuery(cnn, query4)

nonNull <- result4 %>% filter(kegg != "") %>% select(kegg)
nonNull


library(tidyr)
