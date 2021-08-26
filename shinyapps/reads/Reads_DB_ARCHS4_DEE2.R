################################################################################
# for processing read counts meta data from ARCHRS4 and DEE2
# This code creates file "GEO.db" with two tables, "GSEinfo" and "sampleInfo".
# The tables contain the following fields:
#     SampleInfo:
# [1] "GSMs"             "tissue"           "sample_title"
# [4] "sample_series_id" "species"          "SRR_accession"
# [7] "QC_summary"
#
#     GSEinfo
# [1] "GSEID"   "Species" "Samples" "Title"
################################################################################


# Required libraries
library("rhdf5")
library(GEOmetadb)
library("getDEE2")
library(dplyr)
setwd("~/idep-master/idep-master/shinyapps/reads")




#######################################################
# extract series and sample information for ARCHRS4
#######################################################


# Human and mouse matrices downloaded from https://maayanlab.cloud/archs4/download.html
destination_fileH <- "../../data/readCounts/human_matrix_v10.h5"
destination_fileM <- "../../data/readCounts/mouse_matrix_v10.h5"


# Check if ARCHS4 gene expression files was already downloaded, if not in current directory download file form repository
if (!file.exists(destination_fileH)) {
  print("Downloading compressed gene expression matrix.")
  url <- "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
  download.file(url, destination_file, quiet = FALSE)
} else {
  print("Local file already exists.")
}

# Check if mouse file exists
if (!file.exists(destination_fileM)) {
  print("Downloading compressed gene expression matrix.")
  url <- "https://s3.amazonaws.com/mssm-seq-matrix/mouse_matrix_v10.h5"
  download.file(url, destination_file, quiet = FALSE)
} else {
  print("Local file already exists.")
}


# Human meta data
destination_file <- destination_fileH

# View contents of .h5 file
# h5ls(destination_file)

# Retrieve information from compressed data
GSMs <- h5read(destination_file, "meta/samples/geo_accession")
tissue <- h5read(destination_file, "meta/samples/source_name_ch1")
sample_title <- h5read(destination_file, "meta/samples/title")
sample_series_id <- h5read(destination_file, "meta/samples/series_id")
humanNSample <- length(unique(sample_series_id))
GSEs <- unique(sample_series_id)
species <- rep("human", length(GSMs))

sample_infoH <- cbind(GSMs, tissue, sample_title, sample_series_id, species)
H5close()

# Mouse meta data
destination_file <- destination_fileM

# Retrieve information from compressed data
GSMs <- h5read(destination_file, "meta/samples/geo_accession")
tissue <- h5read(destination_file, "meta/samples/source_name_ch1")
sample_title <- h5read(destination_file, "meta/samples/title")
sample_series_id <- h5read(destination_file, "meta/samples/series_id")
mouseNSample <- length(unique(sample_series_id))
GSEs <- c(GSEs, unique(sample_series_id))
species <- rep("mouse", length(GSMs))

sample_infoM <- cbind(GSMs, tissue, sample_title, sample_series_id, species) # combine mouaw sample info
H5close()

# combine sample info for both human and mouse
sample_info_ARCHRS4 <- as.data.frame(rbind(sample_infoH, sample_infoM))
remove(sample_infoM, sample_infoH, sample_series_id, sample_title, species, tissue, GSEs, GSMs)

sample_info_ARCHRS4$SRR_accession <- "" # placeholder for DEE2
sample_info_ARCHRS4$QC_summary <- "" # placeholder for DEE2

# write.table(sample_info_ARCHRS4, "sampleInfo.txt", sep="\t",row.names=F)
# x = read.table("sampleInfo.txt", sep="\t",header=T )

library(Hmisc) # for capitalize function: human --> Human
sample_info_ARCHRS4$species <- capitalize(as.character(sample_info_ARCHRS4$species))
sample_info_ARCHRS4$species <- paste0("ARCHS4_", sample_info_ARCHRS4$species)

sample_info_ARCHRS4$sample_series_id <- as.character(sample_info_ARCHRS4$sample_series_id)
#sample_info_ARCHRS4$sample_series_id <- gsub("\t.*", "", sample_info_ARCHRS4$sample_series_id)
#sample_info_ARCHRS4$sample_series_id <- gsub("c\\(\\\"", "", sample_info_ARCHRS4$sample_series_id)
#sample_info_ARCHRS4$sample_series_id <- gsub("\\\".*", "", sample_info_ARCHRS4$sample_series_id)


#####################################################
## Add GSE ids from GEOmetadb.sqlite 
#####################################################


# unique gsms
archs4_gsm <- unique(sample_info_ARCHRS4$GSMs)
archs4_gsm <- subset(archs4_gsm, grepl("GSM", archs4_gsm))


archs4_gsm <- as.data.frame(archs4_gsm)


archs4_gsm <- as.vector(archs4_gsm)

s <- paste(archs4_gsm) # takes a long time ~30 sec
# s

s <- gsub("c\\(\\\"", "", s) # removes c\
s <- gsub("\"", "'", s) # replace \ with '

# download "GEOmetabd.sqlite" , may take ~20 minutes
if (!file.exists("GEOmetadb.sqlite")) {
  getSQLiteFile()
  print("File Downloaded")
} else {
  "File Already Exists"
}

# connect to file
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")

dbListTables(con)
dbListFields(con, "gse_gsm")
#dbGetQuery(con, "sql script")
query <- paste("select gse, gsm from gse_gsm where gsm IN ('",
               paste(s, collapse = "','"),
               sep = ""
)



archs4_gse_gsm <- dbGetQuery(con, query)


dbDisconnect(con)
remove(s)
#####################





#tem <- as.data.frame(table(sample_info_ARCHRS4[, 4])) # count samples
tem <- as.data.frame(table(archs4_gse_gsm[,1]))
sample_info_ARCHRS4 <- merge(archs4_gse_gsm, sample_info_ARCHRS4, by.x = "gsm", by.y = "GSMs")


GSEinfo1 <- as.data.frame(unique(sample_info_ARCHRS4[, c(2,6)]))
GSEinfo1 <- merge(GSEinfo1, tem, by.x = "gse", by.y = "Var1")

# clear environment
remove(tem, x, destination_file, destination_fileH, destination_fileM,s)

#####################################################
## Add DEE2 datasets
#####################################################

DEE2Species <- c(
  "athaliana",
  "celegans",
  "dmelanogaster",
  "drerio",
  "ecoli",
  "hsapiens",
  "mmusculus",
  "rnorvegicus",
  "scerevisiae"
)

DEE2Species.Display.Name <- c(
  "Arabidopsis",
  "Worm",
  "Fly",
  "Zebrafish",
  "E coli",
  "Human",
  "Mouse",
  "Rat",
  "Yeast"
)

# Reads DEE2 metadata for all species, may take several minutes to run
for (i in 1:length(DEE2Species)) {
  cat("\n", i, DEE2Species.Display.Name[i], "\n")
  mdat <- getDEE2Metadata(DEE2Species[i], outfile = paste0("metadata_", DEE2Species[i]))
  mdat$species <- DEE2Species.Display.Name[i]
  if (i == 1) {
    metaAll <- mdat
  } else {
    metaAll <- rbind(metaAll, mdat)
  }
}

## species --> DEE2_species
metaAll$species <- paste0("DEE2_", metaAll$species)
# find missing gses using gsms from metaAll and  gses in GEOmetabd


#################


# unique gsms
dee2_gsm <- unique(metaAll$Sample_name)
dee2_gsm <- subset(dee2_gsm, grepl("GSM", dee2_gsm))
dee2_gsm <- as.data.frame(dee2_gsm)


dee2_gsm <- as.vector(dee2_gsm)

s <- paste(dee2_gsm) # takes a long time ~30 sec
#s #takes a long time to run

s <- gsub("c\\(\\\"", "", s) # removes c\
s <- gsub("\"", "'", s) # replace \ with '

################


# download "GEOmetabd.sqlite" , may take ~20 minutes
if (!file.exists("GEOmetadb.sqlite")) {
  getSQLiteFile()
  print("File Downloaded")
} else {
  "File Already Exists"
}

# connect to file
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")

dbListTables(con)
dbListFields(con, "gse_gsm")
#dbGetQuery(con, message)
query <- paste("select gse, gsm from gse_gsm where gsm IN ('",
               paste(s, collapse = "','"),
               sep = ""
)

query_gse <- paste("select distinct gse from gse_gsm where gsm IN ('",
                   paste(s, collapse = "','"),
                   sep = ""
)


dee2_gse_gsm <- dbGetQuery(con, query)


dee2_gse <- dbGetQuery(con, query_gse)
dbDisconnect(con)
#####################

dee2 <- metaAll[metaAll$Sample_name %in% dee2_gsm$dee2_gsm, ]



dee2[dee2$Sample_name == "GSM1026903", ]


dee2 <- merge(dee2, dee2_gse_gsm, by.x = "Sample_name", by.y = "gsm")
table(dee2$species)
# tem = as.data.frame( table(metaAll$GEO_series))	# count samples\
tem <- as.data.frame(table(dee2$gse))
# remove missing GSEids and GSEids with 2nd number
# tem <- subset(tem, nchar(as.character(tem$Var1)) >0)
summary(nchar(as.character(tem$Var1)))
tem$Var1 <- gsub(";.*", "", tem$Var1)
GSEinfo2 <- as.data.frame(unique(dee2[, c("gse", "species")]))



GSEinfo2 <- merge(GSEinfo2, tem, by.x = "gse", by.y = "Var1")
colnames(GSEinfo2)[1] <- "sample_series_id"


##################
######################



# merge Dee2 frequencies with ARCHRS4 frequencies
colnames(GSEinfo1) <- colnames(GSEinfo2)
GSEinfo <- rbind(GSEinfo1, GSEinfo2)
remove(GSEinfo1, GSEinfo2)
head(GSEinfo)

#---------------------------
# sample level information table:
# GSMs        tissue    sample_title     sample_series_id species
# sample_info_DEE2 = metaAll[, c("GSM_accession", "experiment_title","GSE_accession","species","SRR_accession","QC_summary" )]
# sample_info_DEE2 = metaAll#[, c("GEO_series","species","SRR_accession","QC_summary", "Sample_name", "Library_name" )]
sample_info_DEE2 <- dee2
# clean up missing GSE info; this removes some samples from DEE2  588,140 --> 419,814
# find geo series from geometadb
# sample_info_DEE2 <- sample_info_DEE2[ which(nchar(sample_info_DEE2$GEO_series ) >0 ) , ]
# sample_info_DEE2$GEO_series<-gsub(";.*","",sample_info_DEE2$GEO_series)
head(sample_info_DEE2)
# remove failed samples??

remove(mdat, metaAll,s)

################################################################################################
# Use GEOmetadb library to retrieve experiment info from the NCBI Gene Expression Omnibus (GEO)

################################################################################################

# download "GEOmetabd.sqlite" , may take ~20 minutes
if (!file.exists("GEOmetadb.sqlite")) {
  getSQLiteFile()
  print("File Downloaded")
} else {
  "File Already Exists"
}

# connect to file
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")


# dbGetQuery(con, "select count(gse), count(gsm) from gse_gsm")

# dbGetQuery(con, "select * from gse_gsm where gsm in ('GSM1000400')")


# view tables and fields
dbListTables(con)
dbListFields(con, "gse")
dbListFields(con, "gsm")
dbListFields(con, "gse_gsm")

dbListFields(con, "gds")


# dbGetQuery(con,"select source_name_ch1, gsm from gsm where gsm IN ('GSM399100', 'GSM399101')")
colnames(GSEinfo)[1] <- "GSEID"


# #get gse and title fields from gse table for included GSE ids
rs <- dbGetQuery(con, paste("select gse,title, summary from gse where gse IN ('",
                            paste(GSEinfo$GSEID, collapse = "','"), "')",
                            sep = ""
))
# rs[,3] <- gsub("\t","",rs[,3] )

# add title column to GSEinfo
GSEinfo <- merge(GSEinfo, rs, by.x = "GSEID", by.y = "gse")



# sort GSEinfo by species and GSEID
GSEinfo <- GSEinfo[order(GSEinfo$species, GSEinfo$GSEID), ]


colnames(GSEinfo) <- c("GSEID", "Species", "Samples", "Title", "Summary")

GSEinfo$Summary_short <- strtrim(GSEinfo$Summary, 1100) # shorten summary length
GSEinfo$summary_short2 <- ifelse(nchar(GSEinfo$Summary_short) == 1100, gsub("\\s*\\w*$", "...", GSEinfo$Summary_short), GSEinfo$Summary_short)

GSEinfo <- GSEinfo[, c(1:4, 7)]
colnames(GSEinfo) <- c("GSEID", "Species", "Samples", "Title", "Summary")
write.table(GSEinfo, "GSEinfo.txt", sep = "\t", row.names = F)


#########################################################################
# Add tissue information for DEE2 samples from "GEOmetabd.sqlite"
#########################################################################

# Preview gsm table
### 	rs1 <- dbGetQuery(con, "select * from gsm LIMIT 3")
# rs1
# Get gsm and source_name_ch1 fields from gsm table
rs1 <- dbGetQuery(con, paste("select gsm, source_name_ch1, title from gsm where gsm IN ('",
                             paste(sample_info_DEE2$Sample_name, collapse = "','"), "')",
                             sep = ""
)) # long time to run  ~ 15 minutes

# disconect from .sqlite file
dbDisconnect(con)

# add source_name_ch1 to sampleInfo_dee2
sampleInfo_dee2 <- merge(sample_info_DEE2, rs1, by.x = "Sample_name", by.y = "gsm")
# sampleInfo_dee2 <- sampleInfo_dee2[ , c(1,7, 2:6)] #What? select entire table?
remove(sample_info_DEE2)


# colnames(sampleInfo_dee2)  <- colnames(sample_info)
colnames(sample_info_ARCHRS4)
colnames(sampleInfo_dee2)


# arrange columns


sample_info_ARCHRS4$title <- ""
# arrange columns of dee2 to match archs4
# a=sampleInfo_dee2 %>% select(Sample_name,source_name_ch1, title, GEO_series, species, SRR_accession, QC_summary, title)

# colnames(sampleInfo_dee2)  <- colnames(sample_info)
colnames(sample_info_ARCHRS4)

colnames(sampleInfo_dee2)

# arrange columns of dee2 to match archs4
a <- sampleInfo_dee2 %>% select(Sample_name, source_name_ch1, title, gse, species, SRR_accession, QC_summary)

sampleInfo_dee2 <- as.data.frame(a)


a <- sample_info_ARCHRS4 %>% select(gsm, tissue, sample_title, gse, species, SRR_accession, QC_summary)

sample_info_ARCHRS4 <- as.data.frame(a)
remove(a)
# change column names for each sampleINfo to match

colnames(sample_info_ARCHRS4)
colnames(sampleInfo_dee2)

colnames(sampleInfo_dee2) <- colnames(sample_info_ARCHRS4)
colnames(sampleInfo_dee2)


# combine sample infos
sampleInfoFinal <- rbind(sample_info_ARCHRS4, sampleInfo_dee2)
# colnames(sampleInfoFinal) <- c()

remove(sample_info_ARCHRS4, sampleInfo_dee2, rs, rs1, tem, dee2_gse, dee2_gse_gsm, dee2_gsm)

###############################################################################
# Create "GEO.db" database
###############################################################################


sqlite <- dbDriver("SQLite")
convert <- dbConnect(sqlite, "GEO.db")

# Write  GSEinfo.txt

dbWriteTable(convert, "GSEinfo", GSEinfo, overwrite = TRUE)

# add index to GSEinfo.txt
dbSendQuery(
  convert,
  "CREATE INDEX index2 ON GSEinfo (GSEID, Species);"
)

# Write sampleInfo.txt

dbWriteTable(convert, "sampleInfo", sampleInfoFinal, overwrite = TRUE)


# Create index on idType
dbSendQuery(
  convert,
  "CREATE INDEX index1 ON sampleInfo (species, sample_series_id);"
)

# View tables/fields
dbListTables(convert)
dbListFields(convert, "sampleInfo")
dbListFields(convert, "GSEinfo")

# disconnect
dbDisconnect(convert)



# debugging
dbGetQuery(convert, "Select * from sampleInfo limit 3")

dbGetQuery(convert, "Select Summary_short from GSEinfo limit 3")
dbGetQuery(convert, "Select count(GSMs) from sampleInfo where sample_series_id = 'GSE100058'")
dbGetQuery(convert, "Select Sample_series_id from sampleInfo where sample_series_id = 'GSE100035'")
dbGetQuery(convert, "Select * from sampleInfo where sample_series_id = 'GSE100035'")
h <- dbGetQuery(convert, "select * from GSEinfo where Samples < 1000 and Samples > 600")
h <- as.data.frame(h)
hist(h$Samples)
boxplot(dbGetQuery(convert, "select Samples from GSEinfo"))
