library('R6')

Ctl.PCA <- R6Class("Ctl.PCA")

###############################################################################
###################			Load/Reload UI Functions		###################
###############################################################################

Ctl.PCA$set("public", "InitColorSelection",
    function(PreprocessSampleInfoResult){
        renderUI({
            if (is.null(PreprocessSampleInfoResult) ){ 
                return(HTML("Upload a sample info file to customize this plot.") ) 
            }else{
                selectInput(
                    "select_PCA_ColorSelection", 
                    label="Color:",
                    choices=c( colnames(PreprocessSampleInfoResult, "Sample_Name"),
                    selected = "Sample_Name")  
            } 
        })
    }
)

Ctl.PCA$set("public", "InitShapeSelection",
    function(PreprocessSampleInfoResult){
        renderUI({
            if (is.null(PreprocessSampleInfoResult) ){ 
                return(NULL) 
            }else{ 
	            tem <- c( colnames(PreprocessSampleInfoResult), "Sample_Name")
	            if(length(tem)>1) { 
                    tem2 = tem[1]
                    tem[1] <- tem[2]
                    tem[1] = tem2
                } # swap 2nd factor with first

	            selectInput(
                    "select_PCA_ShapeSelection", 
                    label="Shape:",
                    choices=tem, 
                    selected = "Sample_Name"
                )
	        } 
        })
    }
)


###############################################################################
###################			Result Ouput Functions			###################
###############################################################################

Ctl.PCA$set("public", "GetMainPlot",  ## not done yet
    function( input, ConvertedTransformedData, PreprocessSampleInfoResult, GeneSetsPCA ){

        withProgress(message=sample(quotes,1), detail ="Running ", {

            if(input$select_PCA_Methods == 1){
                incProgress(1/3,detail="PCA")
                p <- LogicManager$PCA$GetPCAPlot(ConvertedTransformedData, PreprocessSampleInfoResult)
                incProgress(1,detail="Done")
                return(p)
            }else if (input$select_PCA_Methods == 2) {
                incProgress(1/8, detail="PGSEA")
                p <- LogicManager$PCA$plotPathway(ConvertedTransformedData, GeneSetsPCA)
                incProgress(1,detail="Done")
                return(p)
            }else if (input$select_PCA_Methods == 3) {
                incProgress(1/3,detail = " MDS")
                p <- LogicManager$PCA$plotMDS(ConvertedTransformedData, PreprocessSampleInfoResult)
                incProgress(1,detail="Done")
                return(p)
            }else{
                incProgress(1/4, detail=" t-SNE")
                p <- LogicManager$PCA$plotTSNE(ConvertedTransformedData, PreprocessSampleInfoResult, input$btn_PCA_RecalcTSNE)
                incProgress(1,detail="Done")
                return(p)
            }
        }
    }
)


Ctl.PCA$set("public", "ShowCorrelationBetweenPCs",
    function(){
        renderUI({
            if (is.null(input$file1)&& input$goButton == 0)   return(NULL)
            if(is.null(readSampleInfo())) return(NULL)
            npc = 5 # number of Principal components
            x <- readData()$data
            y <- readSampleInfo()

            pca.object <- prcomp(t(x))
            pcaData = as.data.frame(pca.object$x[,1:npc]); 
            pvals = matrix(1,nrow=npc,ncol=ncol(y))
            for (i in 1:npc ){
                for (j in 1:ncol(y) )
                    pvals[i,j] =summary( aov(pcaData[,i] ~ as.factor(y[,j])))[[1]][["Pr(>F)"]][1]
            }
            
            pvals = pvals * npc* ncol(y)   # correcting for multiple testing
            pvals[pvals>1] = 1

            colnames(pvals) = colnames(y)
            rownames(pvals) = paste0("PC",1:npc)
            a="<h4>Correlation between Principal Components (PCs) with factors </h4>"
            nchar0 = nchar(a)
            for ( i in 1:npc){
                j = which.min(pvals[i,])
                if(pvals[i,j]< 0.05) a=paste0(a,rownames(pvals)[i], 
                            " is correlated with ", colnames(pvals)[j],
                            " (p=",   sprintf("%-3.2e",pvals[i,j]),").<br>")
            }
            if(nchar(a) == nchar0 ) return(NULL) else 
            return( HTML(a) )
                
        }) 
    }
)


Ctl.PCA$set("public","")





