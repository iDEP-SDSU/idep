library('R6')

source('server.config')

Display.Manager <- R6Class("Display.Manager")


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
	   	
		p <- plot_ly(x = colnames(dat), y = colSums(dat)/1e6, type = 'bar',
					marker = list(color = columnColor) ) %>%
				layout(title = paste("Total read counts (millions)", memo) )

		return(p)
	}
)

Display.Manager$set("public", "GetTransedDataBoxPlot",
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

		for(i in 1:ncol(dat)){
			p <- p %>% add_boxplot(y = dat[,i], line = list(color = columnColor[i]), 
						marker = list(color = columnColor[i]), name = colnames(dat)[i])
		}

		p <- p %>% layout(title = paste("Distribution of transformed data (millions)", memo) )

		return(p)
	}
)

Display.Manager$set("public", "GetTransedDataDensityPlot",
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

		for( i in 1:ncol(dat) ){
			dens <- density(dat[,i])
			p <- p %>% 
				add_trace(x = dens$x, y = dens$y, name = colnames(dat)[i], color = columnColor[i], mode = 'lines')
		}

		p <- p %>% 
			layout(title = paste("Distribution of transformed data (millions)", memo),
				xaxis = list(title = "Expression values")			
			)

		return(p)
	}
)


