 # for processing read counts meta data from ARCHRS4 and DEE2

library("rhdf5")
setwd("I:/READS/shinyapps/reads")
library(GEOmetadb)  
library("getDEE2")
 
   destination_fileH = "../../data/readCounts/human_matrix.h5"
   destination_fileM = "../../data/readCounts/mouse_matrix.h5"

    infoFile = "GSM_sample_info.csv"
    
    # Check if gene expression file was already downloaded, if not in current directory download file form repository
    if(!file.exists(destination_fileH)){
      print("Downloading compressed gene expression matrix.")
      url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
      download.file(url, destination_file, quiet = FALSE)
    } else{
      print("Local file already exists.")
    }

    # Check if mouse file exists
    if(!file.exists(destination_fileM)){
      print("Downloading compressed gene expression matrix.")
      url = "https://s3.amazonaws.com/mssm-seq-matrix/mouse_matrix.h5"
      download.file(url, destination_file, quiet = FALSE)
    } else{
      print("Local file already exists.")
    }
    
    #######################################################
    # extract series and sample information for ARCHRS4
    #######################################################

   # if(!file.exists(infoFile)){ 
      destination_file = destination_fileH 
      # Retrieve information from compressed data
      GSMs = h5read(destination_file, "meta/Sample_geo_accession")
      tissue = h5read(destination_file, "meta/Sample_source_name_ch1")
      #genes = h5read(destination_file, "meta/genes")
      sample_title = h5read(destination_file, "meta/Sample_title")
      sample_series_id = h5read(destination_file, "meta/Sample_series_id")
      humanNSample=length(unique(sample_series_id) )
	  GSEs = unique(sample_series_id)
      species = rep("human",length(GSMs))
      sample_info = cbind(GSMs, tissue, sample_title,sample_series_id,species)
      H5close()
	  
      destination_file = destination_fileM 
      # Retrieve information from compressed data
      GSMs = h5read(destination_file, "meta/Sample_geo_accession")
      tissue = h5read(destination_file, "meta/Sample_source_name_ch1")
      #genes = h5read(destination_file, "meta/genes")
      sample_title = h5read(destination_file, "meta/Sample_title")
      sample_series_id = h5read(destination_file, "meta/Sample_series_id")
      mouseNSample=length(unique(sample_series_id) )
	  GSEs = c(GSEs, unique(sample_series_id))
      species = rep("mouse",length(GSMs))
      sample_infoM = cbind(GSMs, tissue, sample_title,sample_series_id,species)
      H5close()
	  
	  #sample info for both human and mouse
      sample_info = as.data.frame( rbind(sample_info, sample_infoM) )
      sample_info$SRR_accession = "" # placeholder for DEE2
      sample_info$QC_summary = ""
	  #write.table(sample_info, "sampleInfo.txt", sep="\t",row.names=F)
	  #x = read.table("sampleInfo.txt", sep="\t",header=T )

      library(Hmisc) # for capitalize function: human --> Human
      sample_info$species <- capitalize( as.character(sample_info$species ) )    
      sample_info$species <- paste0( "ARCHS4_", sample_info$species)

	  sample_info$sample_series_id <- as.character( sample_info$sample_series_id )
  	  sample_info$sample_series_id <- gsub("\t.*", "", sample_info$sample_series_id )   
 #   some samples are connected with two c(\"GSE52279\", \"GSE52623\")
  	  sample_info$sample_series_id <- gsub("c\\(\\\"", "", sample_info$sample_series_id )   	
   	  sample_info$sample_series_id <- gsub("\\\".*", "", sample_info$sample_series_id )   	 
	 tem = as.data.frame( table(sample_info[,4]))	# count samples
	 GSEinfo = as.data.frame( unique(sample_info[,4:5]))
     GSEinfo =  merge( GSEinfo, tem, by.x = "sample_series_id", by.y= "Var1")


   
	#####################################################
	## Add DEE2 datasets
	#####################################################
    # https://github.com/markziemann/dee2/blob/master/AccessDEEfromR.md

	DEE2Species = c("athaliana"
					,"celegans"
					,"dmelanogaster"
					,"drerio"
					,"ecoli"
					,"hsapiens"
					,"mmusculus"
					,"rnorvegicus"
					,"scerevisiae")
	DEE2Species.Display.Name = c("Arabidopsis"
					,"Worm"
					,"Fly"
					,"Zebrafish"
					,"E coli"
					,"Human"
					,"Mouse"
					,"Rat"
					,"Yeast")

	for ( i in 1:length(DEE2Species) ) {
	  cat ("\n", i, DEE2Species.Display.Name[i],"\n")
	  mdat <- getDee2Metadata( DEE2Species[i], outfile = paste0("metadata_",DEE2Species[i] ) )
	  mdat$species <- DEE2Species.Display.Name[i]
	  if(i == 1) 
		metaAll <- mdat 
	  else metaAll = rbind(metaAll, mdat)
	  }

	# Human --> Human_DEE2
	metaAll$species = paste0( "DEE2_", metaAll$species )

	tem = as.data.frame( table(metaAll$GSE_accession))	# count samples
	GSEinfo2 = as.data.frame( unique( metaAll[ ,c("GSE_accession", "species")]))
	GSEinfo2 =  merge( GSEinfo2, tem, by.x = "GSE_accession", by.y= "Var1")
	colnames(GSEinfo2)[1] = "sample_series_id"

	# merge with ARCHRS4 data

	GSEinfo <- rbind (GSEinfo, GSEinfo2)     

    #---------------------------
    # sample level information table:          
    #GSMs        tissue    sample_title     sample_series_id species
    sample_info_DEE2 = metaAll[, c("GSM_accession", "experiment_title","GSE_accession","species","SRR_accession","QC_summary" )]
    # clean up missing GSE info; this removes some samples from DEE2  588,140 --> 419,814
    sample_info_DEE2 <- sample_info_DEE2[ which(!is.na(sample_info_DEE2$GSM_accession ) ) , ]
    sample_info_DEE2 <- sample_info_DEE2[ which(!is.na(sample_info_DEE2$GSE_accession ) ) , ]    
     
    # remove failed samples?????

	 ########################################################
     # Use GEOmetadb to retrieve experiment info	 
	 # file can be 
     if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
	 con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
	 dbListTables(con)
	 dbListFields(con,'gse')
	 
	 rs <- dbGetQuery(con,paste("select gse,title from gse where gse IN ('",
                           paste(GSEinfo$sample_series_id,collapse="','"),   "')",sep=""))
    #rs[,3] <- gsub("\t","",rs[,3] )
	GSEinfo =  merge(GSEinfo, rs, by.x = "sample_series_id", by.y= "gse")	
    	

	GSEinfo = GSEinfo[order(GSEinfo$species, GSEinfo$sample_series_id),]
	colnames(GSEinfo) <- c("GSEID", "Species","Samples","Title")
	write.table(GSEinfo, "GSEinfo.txt", sep="\t",row.names=F)

  
    #########################################################################
    # Add tissue information for DEE2 samples
   	 #rs1 <- dbGetQuery(con, "select * from gsm LIMIT 3")

   	 rs1 <- dbGetQuery(con,paste("select gsm,source_name_ch1 from gsm where gsm IN ('",
                           paste(sample_info_DEE2$GSM_accession,collapse="','"),   "')",sep=""))
     
     dbDisconnect(con)
     sampleInfoFinal <- merge(sample_info_DEE2, rs1, by.x ="GSM_accession", by.y = "gsm" )
     sampleInfoFinal <- sampleInfoFinal[ , c(1,7, 2:6)]
     colnames(sampleInfoFinal)  <- colnames(sample_info)
   
     sampleInfoFinal <- rbind(sample_info, sampleInfoFinal)


    ###############################################################################
    # Create database
	
	sqlite  <- dbDriver("SQLite")
	convert <- dbConnect(sqlite,"GEO.db")
	dbWriteTable(convert,"GSEinfo",GSEinfo)
	dbSendQuery( convert,
	"CREATE INDEX index2 ON GSEinfo (GSEID, Species);")

	dbWriteTable(convert,"sampleInfo",sampleInfoFinal)
	# create index on idType
	dbSendQuery( convert,
	"CREATE INDEX index1 ON sampleInfo (species, sample_series_id);")

	dbListTables(convert)
	dbListFields(convert, 'sampleInfo')
	dbListFields(convert, 'GSEinfo')
	dbDisconnect(convert)
	