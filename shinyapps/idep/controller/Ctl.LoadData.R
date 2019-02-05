library('R6')
library(shiny)
library(shinyBS)
library(DT)
library("rhdf5")

source('server.config')

Ctl.LoadData <- R6Class("Ctl.LoadData")

## Fields 
Ctl.LoadData$set("public",	"destination_fileH", "")
Ctl.LoadData$set("public",	"destination_fileH", "")
Ctl.LoadData$set("public",	"sampleInfoFile", "")
Ctl.LoadData$set("public",	"GSEInfoFile", "")
Ctl.LoadData$set("public",	"sample_info", NULL)







## Init
Ctl.LoadData$set(  
	"public",
	"initialize",
	function(){
		self$destination_fileH = paste0(CONFIG_DATA_READCOUNT_PATH, "human_matrix.h5")
		self$destination_fileH = paste0(CONFIG_DATA_READCOUNT_PATH, "mouse_matrix.h5")
		self$destination_fileH = paste0(CONFIG_DATA_READCOUNT_PATH, "sampleInfo.txt")
		self$destination_fileH = paste0(CONFIG_DATA_READCOUNT_PATH, "GSEinfo.txt")

		
		if(file.exists(sampleInfoFile)) {
			self$sample_info =read.table(sampleInfoFile, sep="\t",header=T )
		} else {   # create sample info
			# Check if gene expression file was already downloaded, if not in current directory download file form repository
			if(!file.exists(destination_fileH)) {
				print("Downloading compressed gene expression matrix.")
				url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
				download.file(url, destination_file, quiet = FALSE)
			} else{
				print("Local file already exists.")
			}

			# Check if gene expression file was already downloaded, if not in current directory download file form repository
			if(!file.exists(destination_fileM)){
				print("Downloading compressed gene expression matrix.")
				url = "https://s3.amazonaws.com/mssm-seq-matrix/mouse_matrix.h5"
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
			
			species = rep("mouse",length(GSMs))
			sample_infoM = cbind(GSMs, tissue, sample_title, sample_series_id,species)
			H5close()
			# sample info for both human and mouse
			sample_info = as.data.frame( rbind(sample_info, sample_infoM) )
			#write.table(sample_info, sampleInfoFile, sep="\t",row.names=F)
		}
		
		self$humanNDataset <<- length(unique(sample_info$sample_series_id[which(sample_info$species == "human")]) )
		self$mouseNDataset <<- length(unique(sample_info$sample_series_id[which(sample_info$species == "mouse")]) )
	}
)








## Methods

Ctl.LoadData$set(  
	"public",
	"samples",
	function(){
		renderTable({
        	if (is.null(input$SearchData_rows_selected))   
			{	return(NULL)	}
        	if (is.null(Search()))   
			{	return(as.matrix("No dataset found!")) }
        	Search()$info
      	},bordered = TRUE)
	}
)

Ctl.LoadData$set(  
	"public",
	"DatasetInfo",
	function(){
		reactive({
			dataset.info <- read.table(self$GSEInfoFile, sep="\t",header=T )
			dataset.info$GEO.ID = as.character(dataset.info$GEO.ID)
			return(dataset.info)
      })
	}
)

Ctl.LoadData$set(
	"public",
	"Search",
	function(input){
		reactive({
			if (is.null(input$SearchData_rows_selected)) { 
				return(NULL)
			}

			withProgress(message = "Searching ...", {

				iy = which( self$DatasetInfo()$Species == input$selected.species.archs4 )
				ix = iy[input$SearchData_rows_selected]

				keyword = self$DatasetInfo()$GEO.ID[ix]
          		keyword = gsub(" ","",keyword)
          		ix = which(self$sample_info[,4]== keyword)

          		if(length(ix) == 0){
					return(NULL)
				}else{
					samp = sample_info[ix,1]
					if( names(sort(table(sample_info[ix,5]),decreasing=T))[1] == "human" ){	
						destination_file = self$destination_fileH
					}

            		if( names(sort(table(sample_info[ix,5]),decreasing=T))[1] == "mouse" ){
						destination_file = self$destination_fileM
					}

					samples = h5read(destination_file, "meta/Sample_geo_accession")
            		sample_locations = which(samples %in% samp)

					genes = h5read(destination_file, "meta/genes")
					expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
					tissue = h5read(destination_file, "meta/Sample_source_name_ch1")
					sample_title = h5read(destination_file, "meta/Sample_title")
					H5close()
					incProgress(1/2)
					
					rownames(expression) <-paste(" ",genes)
					colnames(expression) <- paste( samples[sample_locations], sample_title[sample_locations], sep=" ")
					expression <- expression[,order(colnames(expression))]
					tem = sample_info[ix,c(5,1:3)]
					tem = tem[order(tem[,4]),]
					colnames(tem) <- c("Species", "Sample ID","Tissue","Sample Title")
					incProgress(1)

					if(dim(tem)[1]>50) {
						tem = tem[1:50,]
					}
            		
					return( list(info=tem, counts = expression ) )
				}		
			}
		})
	}
)

Ctl.LoadData$set(
	"public",
	"RenderSampleTable",
	function(input){
		renderTable(
			{
				if (is.null(input$SearchData_rows_selected)){
					return(NULL)
				}   

				if (is.null(self$Search())){
					return(as.matrix("No dataset found!"))
				}

				self$Search()$info
      		},
			bordered = TRUE
		)
	}
)



Ctl.LoadData$set(
	"public",
	"",
	function(){

	}
)


