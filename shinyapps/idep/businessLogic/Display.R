library('R6')

source('server.config')

Display.Manager <- R6Class("Display.Manager")


Display.Manager$set("public", "GetReadcountBarPlot",
	function(readCount, groups){
		memo = ""
		dat <- readCount

		if( ncol(dat) > CONST_EAD_PLOT_MAX_SAMPLE_COUNT ){
			part = 1:CONST_EAD_PLOT_MAX_SAMPLE_COUNT
			dat <- dat[,part]
			memo = paste(" (only showing", CONST_EAD_PLOT_MAX_SAMPLE_COUNT, "samples)")
		}

		if(nlevels(groups)<=1 | nlevels(groups) >20){
			columnColor = "green"
		}else{
			columnColor = rainbow(nlevels(groups))[ groups ]	
		}
	   	
		p <- plot_ly(y = colSums(dat)/1e6, x = colnames(dat), type = 'bar',
					marker = list(color = columnColor) ) %>%
				layout(title = paste("Total read counts (millions)", memo) )

		return(p)
	}
)