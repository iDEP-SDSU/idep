library('R6')
library(RSQLite,verbose=FALSE)	# for database connection

source('server.config')

PCA.Logic <- R6Class("PCA.Logic")

PCA.Logic$set("public", "GetPCAPlot",
    function(ConvertedTransformedData,PreprocessSampleInfoResult, colorSelection, shapeSelection){
        ####################################
        # for showing shapes on ggplot2. The first 6 are default. Default mapping can only show 6 types.
        shapes = c(16,17,15,3,7,8,   1,2,4:6,9:15,18:25  )
        
        pca.object <- prcomp(t(ConvertedTransformedData))
        pcaData = as.data.frame(pca.object$x[,1:2]); 
        pcaData = cbind(pcaData,LogicManager$PreProcessing$DetectGroups(colnames(ConvertedTransformedData)) )
		colnames(pcaData) = c("PC1", "PC2", "Sample_Name")
		percentVar=round(100*summary(pca.object)$importance[2,1:2],0)
		if(is.null(PreprocessSampleInfoResult)) { 
			p=ggplot(pcaData, aes(PC1, PC2, color=Sample_Name, shape = Sample_Name)) 
			} else {
			pcaData = cbind(pcaData,PreprocessSampleInfoResult )
			p=ggplot(pcaData, aes_string("PC1", "PC2", color=colorSelection,shape=shapeSelection))  

			}
		 if(ncol(ConvertedTransformedData)<20) # change size depending of # samples
			p <- p + geom_point(size=5)  else if(ncol(ConvertedTransformedData)<50)
			 p <- p + geom_point(size=3)  else 
			 p <- p + geom_point(size=2)
			 

		p <- p+	 scale_shape_manual(values= shapes)	 
		
		p=p + xlab(paste0("PC1: ",percentVar[1],"% variance")) 
		p=p + ylab(paste0("PC2: ",percentVar[2],"% variance")) 
		p=p + ggtitle("Principal component analysis (PCA)") +
            coord_fixed(ratio=1.0)+ 
		    theme(plot.title = element_text(size = 16,hjust = 0.5)) + 
            theme(aspect.ratio=1) +
		    theme(
                axis.text.x = element_text( size = 16),
			    axis.text.y = element_text( size = 16),
			    axis.title.x = element_text( size = 16),
			    axis.title.y = element_text( size = 16) 
            ) +
		    theme(legend.text=element_text(size=16))
	    return(p)
    }
)

PCA.Logic$set("public", "GetPGSEAPlot",
    function(ConvertedTransformedData, GeneSetPCA){
        library(PGSEA,verbose=FALSE)
		pca.object <- prcomp(t(ConvertedTransformedData))
		pca = 100*pca.object$rotation 
		Npca = 5
		Nterms = 6
		if (Npca > dim(pca)[2]) { Npca = dim(pca)[2] } else pca <-  pca[,1:Npca]
		#pca = pca[,1:5]
		if(is.null(GeneSetPCA ) ) return(NULL)  # no species recognized
		if(length(GeneSetPCA ) <= 1 ) return(NULL)
		#cat("\n\nGene Sets:",length( GeneSets()))
		pg = self$myPGSEA(pca,cl=GeneSetPCA,range=c(15,2000),p.value=TRUE, weighted=FALSE,nPermutation=1)
		incProgress(2/8)
		# correcting for multiple testing
		p.matrix = pg$p.result
		tem = p.adjust(as.numeric(p.matrix),"fdr")
		p.matrix = matrix(tem, nrow=dim(p.matrix)[1], ncol = dim(p.matrix)[2] )
		rownames(p.matrix) = rownames(pg$p.result); colnames(p.matrix) = colnames(pg$p.result)


		selected =c()
		for( i in 1:dim(p.matrix)[2]) {
		  tem = which( rank(p.matrix[,i],ties.method='first') <= Nterms)  # rank by P value
		 #tem = which( rank(pg$result[,i],ties.method='first') >= dim(p.matrix)[1]-3.1) # rank by mean
		 names(tem) = paste("PC",i," ", rownames(p.matrix)[tem], sep="" )
		 selected = c(selected, tem)
		}
		rowids = gsub(" .*","",names(selected))
		rowids = as.numeric( gsub("PC","",rowids) )
		pvals = p.matrix[ cbind(selected,rowids) ]
		a=sprintf("%-1.0e",pvals)
		tem = pg$result[selected,]
		rownames(tem) = paste(a,names(selected)); #colnames(tem)= paste("PC",colnames(tem),sep="")
		
		tem = tem[!duplicated(selected),] 
		#tem = t(tem); tem = t( (tem - apply(tem,1,mean)) ) #/apply(tem,1,sd) )

		p <- smcPlot(
            tem,
            scale =  c(-max(tem), max(tem)), 
            show.grid = T, 
            margins = c(1,1, 3, 20),
			col = .rwb,
            cex.lab=0.8
        )
        
        return(p)
    }
)


PCA.Logic$set("public", "GetMDSPlot",
    function(ConvertedTransformedData, PreprocessSampleInfoResult, colorSelection, shapeSelection ){
        ####################################
        # for showing shapes on ggplot2. The first 6 are default. Default mapping can only show 6 types.
        shapes = c(16,17,15,3,7,8,   1,2,4:6,9:15,18:25  )
        fit = cmdscale( LogicManager$UtilFuns$dist2(t(ConvertedTransformedData) ), eig=T, k=2)

		pcaData = as.data.frame(fit$points[,1:2]); pcaData = cbind(pcaData,LogicManager$PreProcessing$DetectGroups(colnames(ConvertedTransformedData)) )
		colnames(pcaData) = c("x1", "x2", "Sample_Name")
		

		if(is.null(PreprocessSampleInfoResult)){ 
		    p=ggplot(pcaData, aes(x1, x2, color=Sample_Name, shape = Sample_Name))  
		} else {
			pcaData = cbind(pcaData,PreprocessSampleInfoResult )
			p=ggplot(pcaData, aes_string("x1", "x2", color=colorSelection,shape=shapeSelection))  
		}
			
		if(ncol(ConvertedTransformedData)<20) # change size depending of # samples
			p <- p + geom_point(size=5)  else if(ncol(ConvertedTransformedData)<50)
			 p <- p + geom_point(size=3)  else 
			 p <- p + geom_point(size=2)
		p <- p+	 scale_shape_manual(values= shapes)	 
		
		p=p+xlab("Dimension 1") 
		p=p+ylab("Dimension 2") 
		p=p+ggtitle("Multidimensional scaling (MDS)")+ coord_fixed(ratio=1.)+ 
		 theme(plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio=1) +
			 theme(axis.text.x = element_text( size = 16),
			   axis.text.y = element_text( size = 16),
			   axis.title.x = element_text( size = 16),
			   axis.title.y = element_text( size = 16) ) +
		theme(legend.text=element_text(size=16))
		   
        return(p)
    }
)

PCA.Logic$set("public", "GetTSNEPlot",
    function(ConvertedTransformedData, PreprocessSampleInfoResult, tsneseed, colorSelection, shapeSelection){
        ####################################
        # for showing shapes on ggplot2. The first 6 are default. Default mapping can only show 6 types.
        shapes = c(16,17,15,3,7,8,   1,2,4:6,9:15,18:25  )
        library(Rtsne,verbose=FALSE)
        set.seed(tsneseed)
        tsne <- Rtsne(t(ConvertedTransformedData), dims = 2, perplexity=1, verbose=FALSE, max_iter = 400)
    
        pcaData = as.data.frame(tsne$Y); 
        pcaData = cbind(pcaData,LogicManager$PreProcessing$DetectGroups(colnames(ConvertedTransformedData)) )
		colnames(pcaData) = c("x1", "x2", "Sample_Name")
		
		#pcaData$Sample_Name = as.factor( pcaData$Sample_Name)

		if(is.null(PreprocessSampleInfoResult)) { 
			p=ggplot(pcaData, aes(x1, x2, color=Sample_Name, shape = Sample_Name)) 
		} else {
			pcaData = cbind(pcaData,PreprocessSampleInfoResult)
			p=ggplot(pcaData, aes_string("x1", "x2", color=colorSelection,shape=shapeSelection)) 
			}
			
		if(ncol(ConvertedTransformedData)<20) # change size depending of # samples
			p <- p + geom_point(size=5)  else if(ncol(ConvertedTransformedData)<50)
			 p <- p + geom_point(size=3)  else 
			 p <- p + geom_point(size=2)
		p <- p+	 scale_shape_manual(values= shapes)	  
		p=p+xlab("Dimension 1") 
		p=p+ylab("Dimension 2") 
		p=p+ggtitle("t-SNE plot")+ coord_fixed(ratio=1.)+ 
		 theme(plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio=1) +
			 theme(axis.text.x = element_text( size = 16),
			   axis.text.y = element_text( size = 16),
			   axis.title.x = element_text( size = 16),
			   axis.title.y = element_text( size = 16) ) +
		theme(legend.text=element_text(size=16))
		
        return(p)
    }
)


# Generate the downloaded dataset
PCA.Logic$set("public", "GeneratePCAData",
    function(PreprocessResultData, TSNESeed){
        x <- PreprocessResultData
        
        result = prcomp(t(x))$x[,1:2]
        fit = cmdscale( LogicManager$UtilFuns$dist2(t(x) ), eig=T, k=2)
        result = cbind( result, fit$points[,1:2] )
        library(Rtsne,verbose=FALSE)
        set.seed(TSNESeed)
        tsne <- Rtsne(t(x), dims = 2, perplexity=1, verbose=FALSE, max_iter = 400)
        
        result = cbind( result, tsne$Y)
        
        colnames(result) = c("PCA.x","PCA.y","MDS.x", "MDS.y", "tSNE.x", "tSNE.y")
        return(result)	
    }
)

# Runs pathway analysis using PGSEA; this is copied and revised from PGSEA package
PCA.Logic$set("public", "myPGSEA",
    function(exprs, cl, range = c(25, 500), ref = NULL, center = TRUE, 
        p.value = 0.005, weighted = TRUE, nPermutation=100, enforceRange = TRUE, ...) {
        if (is(exprs, "ExpressionSet")) 
            exprs <- exprs(exprs)
        if (!is.list(cl)) 
            stop("cl need to be a list")
        if (!is.null(ref)) {
            if (!is.numeric(ref)) 
                stop("column index's required")
        }
        if (!is.null(ref)) {
            if (options()$verbose) 
                cat("Creating ratios...", "\n")
            ref_mean <- apply(exprs[, ref], 1, mean, na.rm = TRUE)
            exprs <- sweep(exprs, 1, ref_mean, "-")
        }
        if (center) 
            exprs <- scale(exprs, scale = FALSE)         # column centering is done
        results <- matrix(NA, length(cl), ncol(exprs))
        rownames(results) <- names(cl)
        colnames(results) <- colnames(exprs)
        mode(results) <- "numeric"
        Setsize = c(rep(0,length(cl)))     # gene set size vector
        mean2 = c(rep(0,length(cl)))     # mean of the range of means 
        meanSD = c(rep(0,length(cl)))     # SD of the range of means	
        if (is.logical(p.value)) 
            { p.results <- results; mean.results <- results;}
        for (i in 1:length(cl)) {              # for each gene list
            #cat("\nProcessing gene set",i);
            if (class(cl[[i]]) == "smc") {
                clids <- cl[[i]]@ids
            }
            else if (class(cl[[i]]) %in% c("GeneColorSet", "GeneSet")) {
                clids <- cl[[i]]@geneIds
            }
            else {
                clids <- cl[[i]]
            }
            if (options()$verbose) 
                cat("Testing region ", i, "\n")
            ix <- match(clids, rownames(exprs))
            ix <- unique(ix[!is.na(ix)])
            present <- sum(!is.na(ix))
            Setsize[i] <- present 
            if (present < range[1]) {
                if (options()$verbose) 
                    cat("Skipping region ", i, " because too small-", 
                    present, ",\n")
                next
            }
            if (present > range[2]) {
                if (options()$verbose) 
                    cat("Skipping region ", i, " because too large-", 
                    present, "\n")
                next
            }
            texprs <- exprs[ix, ]           # expression matrix for genes in gene set
            if (any(is.na(texprs))) 
                cat("Warning - 'NA' values within expression data, enrichment scores are estimates only.\n")
            if (!is.matrix(texprs)) 
                texprs <- as.matrix(texprs)
                                
            stat <- try(apply(texprs, 2, t.test, ...))
            means <- try(apply(texprs, 2, mean,trim=0.1))   # trim mean
            ps <- unlist(lapply(stat, function(x) x$p.value))
            stat <- unlist(lapply(stat, function(x) x$statistic))
            p.results[i, ] <- ps
            mean.results[i,] <- means
            results[i, ] <- as.numeric(stat)
            
            # permutation of gene sets of the same size
            if(nPermutation > 2 )  { # no permutation if <=2
                meansRanges = c(0, rep(nPermutation))
                for( k in 1:nPermutation ) {
                    ix <- sample.int( dim(exprs)[1], length(ix) )
                    texprs <- exprs[ix, ] 
                    means <- try(apply(texprs, 2, mean,trim=0.1))   # trim mean
                    meansRanges[k] = dynamicRange(means)
                }
                mean2[i] = mean(meansRanges)
                meanSD[i]= sd(meansRanges,na.rm=TRUE)   # NA are removed before calculating standard deviation
            }
        }
        return(list(results = results, p.results = p.results, means = mean.results, size=Setsize, mean2=mean2, meanSD=meanSD))
    }
)

