library('R6')
library(gplots,verbose=FALSE)		# for hierarchical clustering
library(ggplot2,verbose=FALSE)	# graphics

source('server.config')

Display.Manager <- R6Class("Display.Manager")
Display.Manager$set("public", "HeatColors" , NULL)

Display.Manager$set("public", "initialize",
	function(){

		hmcols <- colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF",
			"#E0F3F8", "#91BFDB", "#4575B4")))(75)

		heatColors = list(  
			"Green-Black-Red" = greenred(75),     
			"Blue-White-Red" = bluered(75),     
	        "Green-Black-Magenta" = colorpanel(75,"green", "black","magenta"),
	        "Blue-Yellow-Red" = colorpanel(75,"blue", "yellow","red"), 
			"Blue-white-brown" = hmcols 
		)

		self$HeatColors <- heatColors
	}
)

Display.Manager$set("public", "GetReadcountBarPlot",
	function(readCount, groups){
		memo = ""
		dat <- readCount

		if( ncol(dat) > CONST_EAD_PLOT_MAX_SAMPLE_COUNT ){
			dat <- dat[,1:CONST_EAD_PLOT_MAX_SAMPLE_COUNT]
			groups <- groups[1:CONST_EAD_PLOT_MAX_SAMPLE_COUNT]
			memo = paste(" (only showing", CONST_EAD_PLOT_MAX_SAMPLE_COUNT, "samples)")
		}

		if(nlevels(groups)<=1 | nlevels(groups) >20){
			columnColor = "green"
		}else{
			columnColor = rainbow(nlevels(groups))[ groups ]	
		}
	   	
		tbl <- data.frame( x=colnames(dat), y= colSums(dat)/1e6 )
		tbl$x <- factor(tbl$x, levels = c(as.character(tbl$x))) 

		p <- plot_ly(data = tbl, x = ~x, y = ~y, type = 'bar',
					marker = list(color = columnColor) ) %>%
				layout(	title = paste("Total read counts (millions) of each sample", memo),
						xaxis = list(title="Sample"),
						yaxis = list(title="Total read couns (millions)") )

		return(p)
	}
)

Display.Manager$set("public", "GetTransformedDataBoxPlot",
	function(transformData, groups){
		# 1. Init/load vars
		memo = ""
		dat <- transformData

		# 2. Check sample count. If count is to much, then show first CONST_EAD_PLOT_MAX_SAMPLE_COUNT samples.
		if( ncol(dat) > CONST_EAD_PLOT_MAX_SAMPLE_COUNT ){
			dat <- dat[, 1:CONST_EAD_PLOT_MAX_SAMPLE_COUNT]
			groups <- groups[1:CONST_EAD_PLOT_MAX_SAMPLE_COUNT]
			memo = paste("(only showing", CONST_EAD_PLOT_MAX_SAMPLE_COUNT, "samples)")
		}

		# 3. Define colors
		if(nlevels(groups)<=1 | nlevels(groups) >20){
			columnColor = "green"
		}else{
			columnColor = rainbow(nlevels(groups))[ groups ]	
		}

		# 4. Build plot
		p <- plot_ly( type="box" )
		legGroupsList <- c()			# this stores the groups that already showing in the legend

		for(i in 1:ncol(dat)){
			
			groupName = as.character(groups[i])

			if( groups[i] %in% legGroupsList ){
				# if this group name already in list, then don't need show the leg
				showLeg = FALSE
			}else{
				# if this group name is not in list, then show the leg as 'group leg'
				showLeg = TRUE
				legGroupsList <- c(legGroupsList, groupName)
			}


			p <- p %>% 
					add_boxplot(y = dat[,i], line = list(color = columnColor[i]), 
						marker = list(color = columnColor[i]), legendgroup = groupName, showlegend = showLeg, 
						hoverinfo = colnames(dat)[i], name = groupName, x = colnames(dat)[i]
					)  
					# legendgroup = groupName defines which group it belongs to  
					# name = groupName defines the name showing on the legend
		}

		p <- p %>% layout(title = paste("Distribution of transformed data", memo) )

		return(p)
	}
)

Display.Manager$set("public", "GetTransformedDataDensityPlot",
	function(transformData, groups){
		# 1. Init/load vars
		memo = ""
		dat <- transformData

		# 2. Check sample count. If count is to much, then show first CONST_EAD_PLOT_MAX_SAMPLE_COUNT samples.
		if( ncol(dat) > CONST_EAD_PLOT_MAX_SAMPLE_COUNT ){
			dat <- dat[, 1:CONST_EAD_PLOT_MAX_SAMPLE_COUNT]
			groups <- groups[1:CONST_EAD_PLOT_MAX_SAMPLE_COUNT]
			memo = paste("(only showing", CONST_EAD_PLOT_MAX_SAMPLE_COUNT, "samples)")
		}

		# 3. Define colors
		if(nlevels(groups)<=1 | nlevels(groups) >20){
			columnColor = "green"
		}else{
			columnColor = rainbow(nlevels(groups))[ groups ]	
		}

		# 4. Build plot
		p <- plot_ly(type = "scatter")
		legGroupsList <- c()			# this stores the groups that already showing in the legend

		for( i in 1:ncol(dat) ){
			# Step 1. check this column belong to a new group or not
			groupName = as.character(groups[i])

			if( groups[i] %in% legGroupsList ){
				# if this group name already in list, then don't need show the leg
				showLeg = FALSE
			}else{
				# if this group name is not in list, then show the leg as 'group leg'
				showLeg = TRUE
				legGroupsList <- c(legGroupsList, groupName)
			}

			# Step 2. add this column to plot
			dens <- density(dat[,i])
			p <- p %>% 
				add_trace(x = dens$x, y = dens$y,  legendgroup = groupName, showlegend = showLeg,
				hoverinfo = colnames(dat)[i], name = groupName, line = list(color = columnColor[i]), mode = 'lines')
		}

		p <- p %>% 
			layout(title = paste("Density plot of transformed data", memo),
				xaxis = list(title = "Expression values")			
			)

		return(p)
	}
)

Display.Manager$set("public", "GetTransformedDataScatterPlot",
	function(transformData){
		# 1. Init/load vars
		memo = ""
		dat <- transformData

		# 2. Select first two columns
		dat <- dat[,1:2]

		# 3. Build Plot 
		p <- plot_ly(x = dat[,1], y = dat[,2], mode = 'markers') %>%
				layout(
					title = "Scatter plot of first two samples",
					xaxis = list(title = colnames(dat)[1] ),
					yaxis = list(title = colnames(dat)[2] )
				)

		return(p)
	}
)


Display.Manager$set("public", "GetHeatmap2Plot",
	function(dat, geneCount, groups, isSampleClustering,
		distfunName, hclustfunName, heatColorName){
		#http://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
		groups.colors = rainbow(length(unique(groups)))

		lmat = rbind(c(0, 4), c(0, 1), c(3, 2), c(5, 0))
		lwid = c(2, 6)
		lhei = c(1.5, .2, 8, 1.1)
	
		if(ncol(dat) < 20){
			cexFactor = 2 
		}else if(ncol(dat) < 31){
			cexFactor = 1.5
		}else{
			cexFactor = 1
		}

		par(mar = c(5, 4, 1.4, 0.2))

		if( geneCount>110 ){
			heatmap.2(dat, 
				distfun = LogicManager$UtilFuns$DistanceFuns[[distfunName]],
				hclustfun = LogicManager$UtilFuns$HierarchicalClusteringFuns[[hclustfunName]],
				Colv = isSampleClustering,
				col = self$HeatColors[[heatColorName]],
				density.info="none", 
				trace="none", 
				scale="none", 
				keysize=.5,
				key=TRUE, 
				symkey=FALSE,
				ColSideColors=groups.colors[ as.factor(groups)],
				labRow="",
				margins=c(10,0),
				srtCol=45,
				cexCol=cexFactor,  # size of font for sample names
				lmat = lmat, 
				lwid = lwid, 
				lhei = lhei
			)
		}else{
			heatmap.2(dat, 
				distfun = LogicManager$UtilFuns$DistanceFuns[[distfunName]],
				hclustfun = LogicManager$UtilFuns$HierarchicalClusteringFuns[[hclustfunName]],
				Colv = isSampleClustering,
				col = self$HeatColors[[heatColorName]],
				density.info="none", 
				trace="none", 
				scale="none", 
				keysize=.5,
				key=TRUE, 
				symkey=FALSE,
				ColSideColors=groups.colors[ as.factor(groups)],
				margins=c(18,12),
				cexRow=1,
				srtCol=45,
				cexCol=cexFactor,  # size of font for sample names
				lmat = lmat, 
				lwid = lwid, 
				lhei = lhei
			)
		}

		if(length(unique(groups) ) <= 30 ) {  # only add legend when there is less categories
			par(lend = 1)           # square line ends for the color legend
			self$add_legend(
				"topleft",
				legend = unique(groups), # category labels
				col = groups.colors[unique(as.factor(groups))],  # color key
				lty= 1,             # line style
				lwd = 10            # line width
			)
		}
	}
)

# This is a support function for GetHeatmap2Plot
# adding sample legends to heatmap; this is for the main heatmap
# https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
Display.Manager$set("public", "add_legend",
	function(...) {
		opar <- par(
			fig=c(0, 1, 0, 1), 
			oma=c(0, 0, 0, 0), 
		 	mar=c(0, 0, 0, 6), 
			new=TRUE
		)
		on.exit(par(opar))
		plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
		legend(...)
	}
)

Display.Manager$set("public", "GetHeatmapPlotly",
	function(dat, selectedHeatColor){

		# Generate Color names
		colorNames = unlist(strsplit(tolower(selectedHeatColor),"-" ) )

		# Build ggplot
		p <- dat %>%
	  		ggplot(aes(X, Y, fill = value)) + 
			  	geom_tile() + 
			  	scale_fill_gradient2(low = colorNames[1], mid = colorNames[2],high = colorNames[3]) + 
			  	theme(
					axis.title.y=element_blank(),   # remove y labels
				   	#axis.text.y=element_blank(),  # keep gene names for zooming
					axis.ticks.y=element_blank(),
					axis.title.x=element_blank()
				) + 
				theme(axis.text.x = element_text(size=10,angle = 45, hjust = 1))
		
		p <- ggplotly(p) %>% 
			layout(margin = list(b = 150,l=200))

		return(p)
	}
)

Display.Manager$set("public", "GetSDHeatmapPlot",
	function(dat, Cutoff){
		p <- ggplot(dat, aes(x=SDs)) + 
			geom_density(color="darkblue", fill="lightblue") +
		  	labs(x = "Standard deviations of all genes", y="Density")+
		  	geom_vline(
				aes(xintercept=Cutoff),
				color="red", linetype="dashed", size=1
			) + 
			annotate(
				"text", 
				x = Cutoff + 0.4*sd(dat[,1]), 
				y = 1,
				color = "red", 
				label = paste0("Top ", Cutoff)
			) +
			theme(
				axis.text=element_text(size=14),
				axis.title=element_text(size=16,face="bold")
			)
		return(p)
	}
)

# heatmap of correlation matrix
Display.Manager$set("public", "GetHeatmapOfCorrelationMatrix",
	function(dat,isLabelWithPCC){
		# Create a ggheatmap
		ggheatmap <- ggplot(dat, aes(Var2, Var1, fill = value))+
			geom_tile(color = "white")+
			scale_fill_gradient2(low = "green", high = "red",  mid = "white", 
			space = "Lab",  limit = c(min(dat[,3]) ,max(dat[,3])), midpoint = median(dat[,3]),
			name="Pearson's \nCorrelation") +
			theme_minimal()+ # minimal theme
			theme(axis.text.x = element_text(angle = 45, vjust = 1, size=14,hjust = 1))+
			theme(axis.text.y = element_text( size = 14 ))+
			coord_fixed()

		if(isLabelWithPCC && ncol(dat)<20){
			ggheatmap <- ggheatmap +  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4)
		}
		
		ggheatmap <- ggheatmap + 
		  	theme(
				axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				panel.grid.major = element_blank(),
				panel.border = element_blank(),
				panel.background = element_blank(),
				axis.ticks = element_blank(),
				legend.justification = c(1, 0),
				legend.position = c(0.6, 0.7),
				legend.direction = "horizontal"
			) +
			guides(fill = FALSE)
		
		return(ggheatmap)
	}
)

# Sample Tree plot for heatmap tab
Display.Manager$set("public", "GetSampLeTreePlot",
	function(dat, hclustfunName, distfunName){
		hclustfun <- LogicManager$UtilFuns$HierarchicalClusteringFuns[[hclustfunName]]
		distfun <- LogicManager$UtilFuns$DistanceFuns[[distfunName]]

		p <- plot(
			as.dendrogram(
				hclustfun(
					distfun( t(dat) )
				)
			),
			xlab = "",
			ylab = paste(distfunName, "(", hclustfunName, "linkage", ")" ),
			type = "rectangle"
		)

		return(p)
	}
)

# bar plot of Single Gene for individual samples
Display.Manager$set("public", "GetBarPlotSingleGeneOfIndividualSamples",
	function(dat, ymax){
		p <- ggplot(data=dat, aes(x=samples, y=value, group = Genes, shape=Genes, colour = Genes)) +
			geom_line() +
			geom_point( size=5,  fill="white") + 
			labs(y="Transformed expression level") +
			coord_cartesian(ylim = c(0, ymax))

		p <- p + theme(plot.title = element_text(size = 16,hjust = 0.5)) +
	 		theme(axis.text.x = element_text(angle=45, size = 16, hjust=1),
	    		axis.text.y = element_text(size = 16),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 16) 
			) +
			theme(legend.text=element_text(size=12))

		return(p)
	}
)

# bar plot of Single gene for all samples
Display.Manager$set("public", "GetBarPlotSingleGeneOfAllSamples",
	function(dat, isUseSD){
		#http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization
		if(isUseSD){
			p <- ggplot(dat, aes(x=Genes, y=Mean,fill=Samples) ) + # data & aesthetic mapping
				geom_bar(stat="identity", position=position_dodge()) + # bars represent average
				geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.2,position=position_dodge(.9)) +
				labs(y="Expression Level")
		}else{
			p <- ggplot(dat, aes(x=Genes, y=Mean,fill=Samples) ) + # data & aesthetic mapping
				geom_bar(stat="identity", position=position_dodge()) + # bars represent average
				geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.2,position=position_dodge(.9)) +
				labs(y="Expression Level") 
		}

		p <- p + theme(plot.title = element_text(size = 16,hjust = 0.5)) + # theme(aspect.ratio=1) +
	 		theme(axis.text.x = element_text(angle=45, size = 16, hjust=1),
	       		axis.text.y = element_text( size = 16),
		   		axis.title.x = element_blank(),
		   		axis.title.y = element_text( size = 16) 
			) +
			theme(legend.text=element_text(size=16))

		return(p)
	}
)




# heatmap with color bar define gene groups
Display.Manager$set("public", "GetHeatmapWithGeneGroups",
	function(x,bar=NULL,n=-1,mycolor=1,clusterNames=NULL,sideColors=NULL ) {
		# number of genes to show
		ngenes = as.character( table(bar))
		if(length(bar) >n && n != -1) {ix = sort( sample(1:length(bar),n) ); bar = bar[ix]; x = x[ix,]  }
		if(! is.null(bar) )
			if(is.null(sideColors) ) 
				sideColors = mycolors

		# this will cutoff very large values, which could skew the color 
		x=as.matrix(x)-apply(x,1,mean)
		cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
		x[x>cutoff] <- cutoff
		cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
		x[x< cutoff] <- cutoff
		#colnames(x)= detectGroups(colnames(x))
		if(is.null(bar)) # no side colors
			heatmap.2(x,  Rowv =F,Colv=F, dendrogram ="none",
				col=self$HeatColors[as.integer(mycolor)], density.info="none", trace="none", scale="none", keysize=.3
				,key=F, labRow = F,
				#,RowSideColors = mycolors[bar]
				,margins = c(8, 24)
				,srtCol=45
			) else
			heatmap.2(x,  Rowv =F,Colv=F, dendrogram ="none",
				col=self$HeatColors[as.integer(mycolor)], density.info="none", trace="none", scale="none", keysize=.3
				,key=F, labRow = F,
				,RowSideColors = sideColors[bar]
				,margins = c(8, 24)
				,srtCol=45
			)
			
		if(!is.null(bar)) { 
			legend.text = paste("Cluster ", toupper(letters)[unique(bar)], " (N=", ngenes,")", sep="") 
			if( !is.null( clusterNames ) && length(clusterNames)>= length( unique(bar) ) )  
				legend.text = paste(clusterNames[ 1:length( unique(bar) )  ], " (N=", ngenes,")", sep="") 
			
			par(lend = 1)           # square line ends for the color legend
			legend("topright",      # location of the legend on the heatmap plot
			legend = legend.text, # category labels
			col = sideColors,  # color key
			lty= 1,             # line style
			lwd = 10 )           # line width
		}
	}
)


# a program for ploting enrichment results by highlighting the similarities among terms
# must have columns: Direction, adj.Pval   Pathways Genes
#  Direction	adj.Pval	nGenes	Pathways		Genes
#Down regulated	3.58E-59	131	Ribonucleoprotein complex biogenesis	36	Nsun5 Nhp2 Rrp15 
#Down regulated	2.55E-57	135	NcRNA metabolic process	23	Nsun5 Nhp2 Rrp15 Emg1 Ddx56 Rsl1d1 enrichmentPlot <- function( enrichedTerms){
# Up or down regulation is color-coded
# gene set size if represented by the size of marker
Display.Manager$set("public", "GetEnrichmentPlot",
    function( enrichedTerms, rightMargin=33) {
        if(class(enrichedTerms) != "data.frame") return(NULL)
        if(nrow(enrichedTerms) <=1 ) return(NULL)  # only one term or less
        library(dendextend) # customizing tree
        
        geneLists = lapply(enrichedTerms$Genes, function(x) unlist( strsplit(as.character(x)," " )   ) )
        names(geneLists)= enrichedTerms$Pathways

        # compute overlaps percentage--------------------

        n = length(geneLists)
        w <- matrix(NA, nrow = n, ncol = n)
        # compute overlaps among all gene lists
            for (i in 1:n) {
                for (j in i:n) {
                    u <- unlist(geneLists[i])
                    v <- unlist(geneLists[j])
                    w[i, j] = length(intersect(u, v))/length(unique(c(u,v)))
                }
            }
        # the lower half of the matrix filled in based on symmetry
            for (i in 1:n) 
                for (j in 1:(i-1)) 
                    w[i, j] = w[j,i] 
        

        # compute overlaps P value---------------------
        if(0) {
        total_elements = 30000
        n = length(geneLists)
        w <- matrix(rep(0,n*n), nrow = n, ncol = n)
        # compute overlaps among all gene lists
            for (i in 1:n) {
                for (j in (i+1):n) {
                    u <- unlist(geneLists[i])
                    v <- unlist(geneLists[j])
                    xx= length( intersect(u, v) )
                    if(xx == 0)
                        next;
                    mm = length(u)
                    nn <- total_elements - mm	
                    kk = length(v)
                    w[i,j] = -sqrt( -phyper(xx-1,mm,nn,kk, lower.tail=FALSE,log.p = TRUE ));
                    
                }
            }
            

        # the lower half of the matrix filled in based on symmetry
            for (i in 1:n) 
                for (j in 1:(i-1)) 
                    w[i, j] = w[j,i] 
                    
            # w =  w-min(w) 			
            # for( i in 1:n) 		w[i,i] = 0;
        
        }

        Terms = paste( sprintf("%-1.0e",as.numeric(enrichedTerms$adj.Pval)), 
                        names(geneLists))
        rownames(w) = Terms
        colnames(w) = Terms
        par(mar=c(0,0,1,rightMargin)) # a large margin for showing 

        dend <- as.dist(1-w) %>%
            hclust (method="average") 
        ix = dend$order # permutated order of leaves

        leafType= as.factor( gsub(" .*","", enrichedTerms$Direction[ix] ) )
        if(length(unique(enrichedTerms$Direction)  ) ==2 )
            leafColors = c("green","red")  else  # mycolors
            leafColors = mycolors
            
        #leafSize = unlist( lapply(geneLists,length) ) # leaf size represent number of genes
        #leafSize = sqrt( leafSize[ix] )  
        leafSize = -log10(as.numeric( enrichedTerms$adj.Pval[ix] ) ) # leaf size represent P values
        leafSize = 1.5*leafSize/max( leafSize ) + .2
        
            dend %>% 
            as.dendrogram(hang=-1) %>%
            set("leaves_pch", 19) %>%   # type of marker
            set("leaves_cex", leafSize) %>% #Size
            set("leaves_col", leafColors[leafType]) %>% # up or down genes
            plot(horiz=TRUE)
            
        #legend("top",pch=19, col=leafColors[1:2],legend=levels(leafType),bty = "n",horiz =T  )
        # add legend using a second layer
            par(lend = 1)           # square line ends for the color legend
            add_legend("top",pch=19, col=leafColors,legend=levels(leafType),bty = "n",horiz =T 

            )
    }
)
