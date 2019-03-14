library('R6')
library(gplots,verbose=FALSE)		# for hierarchical clustering
library(ggplot2,verbose=FALSE)	# graphics

source('server.config')

Display.Manager <- R6Class("Display.Manager")
Display.Manager$set("public", "heatColors" , NULL)
Display.Manager$set("public", "hmcols" , NULL)

Display.Manager$set("public", "initialize",
	function(){

		self$hmcols <- colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF",
			"#E0F3F8", "#91BFDB", "#4575B4")))(75)

		heatColors = rbind(  greenred(75),     bluered(75),     
	                     colorpanel(75,"green", "black","magenta"),
	                     colorpanel(75,"blue", "yellow","red"), self$hmcols )
		rownames(heatColors) = c("Green-Black-Red", "Blue-White-Red", "Green-Black-Magenta",
	                         "Blue-Yellow-Red", "Blue-white-brown")

		self$heatColors <- heatColors
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
	function(transedData, groups){
		# 1. Init/load vars
		memo = ""
		dat <- transedData

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
	function(transedData, groups){
		# 1. Init/load vars
		memo = ""
		dat <- transedData

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
	function(transedData){
		# 1. Init/load vars
		memo = ""
		dat <- transedData

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

