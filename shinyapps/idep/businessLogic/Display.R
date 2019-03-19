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




