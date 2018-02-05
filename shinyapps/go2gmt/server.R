# This app converts gene-GO id mappings to a GMT file. 
library("GO.db")
library("GSEABase")

MinList = 2; MaxList = 5000
speciesCode = "Species"
##########################################################
# Function to generate gene set for a GO term
      gset <- function( x ){
        a=""
        if(length(geneIds(x)) >= MinList &&length(geneIds(x)) <= MaxList ) {
          GOID <- setName(x)
          setID <-paste("GO",Ontology(GOID),sep="")
          setID <- paste(setID, speciesCode, GOID,gsub(" ", "_",Term(GOID)), sep="_")
          setID <- substr( setID,1,255)  # limit label/ID to 255 characters
          
          
          detailsUrl = paste("http://amigo.geneontology.org/cgi-bin/amigo/term_details?term=",GOID,sep="")
          a = paste(setID,"\t",detailsUrl)
          a=paste(a,"\t",paste( cleanGeneSet(geneIds(x)),collapse="\t"),"\n",sep="")
        }
        a <- gsub("'"," ", a)   # Some GO terms have ' character in them 5', this causes problem
        return(a)
      }
      cleanGeneSet <- function (x){
        # remove duplicate; upper case; remove special characters
        x <- unique( toupper( gsub("\n| ","",x) ) )
        x <- x[ which( nchar(x)>1) ]  # genes should have at least two characters
        return(x)
      }
  
# Define server logic required to draw a histogram
server <- function(input, output) {
  
  convertData <- reactive({
    if (is.null(input$files))   return(NULL)
    inFile <- input$files
    inFile <- inFile$datapath
    
    x = read.table(inFile, sep="\t")
    
    withProgress(message="Converting ", {
      
	x[,1] = gsub("\\..*","",x[,1])
	  x[,1] <- as.character(x[,1])
      x[,2] <- as.character(x[,2])
      genes=c(); goids = c()
      for(i in 1:dim(x)[1]){
        tem = unlist(strsplit(x[i,2],"; ") )
        tem = gsub(" .*","",tem)
        
        for(tem2 in tem){
          genes=c(genes,x[i,1]); goids = c(goids, tem2)
        }
      }
      
    
      mapping = as.data.frame( cbind(goids, rep("IEA",length(genes)),  genes))
      
      colnames(mapping)=c("go","evidence","gene")

      incProgress(1/5)
	# convert to GMT
      ensembl_dataset="atauschii_eg_gene"
      idType ="ensembl_gene_id"
      
      
      goFrame=GOFrame(mapping);
      goAllFrame=GOAllFrame(goFrame); rm(goFrame)
      #build gene set collection
	  incProgress(1/3)
      gsc <-  GeneSetCollection(goAllFrame, setType = GOCollection())
      # toGmt(gsc,paste(species,"_GO_",idType,".gmt",sep=""));  # save gene set to file  toGmt() is a GSEABase method
      # paste lots of string inside a loop is very slow, using lapply
      geneSet <- paste(lapply(gsc,gset),collapse="",sep="")
      
	  incProgress(1)
      
    })
    
    return(list( geneSet=geneSet, n=length(gsc)) )
  })
  
  output$stat <- renderText({
    if (is.null(input$files))   return(NULL)
    if (is.null(convertData()))   return(NULL)   
    
    paste(convertData()$n , " gene sets converted.")
    
    
  })
  
  output$download.File <- downloadHandler(
    filename = function() {"converted.GMT"},
    content = function(file) {
      write(convertData()$geneSet, file)
    }
  )
}
