 # for processing read counts meta data from ARCHRS4

library("rhdf5")
 #setwd("C:/Users/Xijin.Ge/Google Drive/research/expression data cell line")
 library(GEOmetadb)  

 
 
	destination_fileH = "/docker2/idep/countsData/human_matrix.h5"
    destination_fileM = "/docker2/idep/countsData/mouse_matrix.h5"
    infoFile = "GSM_sample_info.csv"
    
    # Check if gene expression file was already downloaded, if not in current directory download file form repository
    if(!file.exists(destination_fileH)){
      print("Downloading compressed gene expression matrix.")
      url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
      download.file(url, destination_file, quiet = FALSE)
    } else{
      print("Local file already exists.")
    }
    
    # Check if gene expression file was already downloaded, if not in current directory download file form repository
    if(!file.exists(destination_fileM)){
      print("Downloading compressed gene expression matrix.")
      url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
      download.file(url, destination_file, quiet = FALSE)
    } else{
      print("Local file already exists.")
    }
    
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
	  
	        # sample info for both human and mouse
      sample_info = as.data.frame( rbind(sample_info, sample_infoM) )
	  	write.table(sample_info, "sampleInfo.txt", sep="\t",row.names=F)
	x = read.table("sampleInfo.txt", sep="\t",header=T )
	  
	  
	 tem = as.data.frame( table(sample_info[,4]))	# count samples
	 GSEinfo = as.data.frame( unique(sample_info[,4:5]))
     GSEinfo =  merge( GSEinfo, tem, by.x = "sample_series_id", by.y= "Var1")
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
    	
	dbDisconnect(con) 
	GSEinfo = GSEinfo[order(GSEinfo$species, GSEinfo$sample_series_id),]
	colnames(GSEinfo) <- c("GEO ID", "Species","Samples","Title")
	write.table(GSEinfo, "GSEinfo.txt", sep="\t",row.names=F)
	x = read.table("GSEinfo.txt", sep="\t",header=T )
	  
	  
	  
	  
	  
	  
	  
	  
	  