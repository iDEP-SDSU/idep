library('R6')

Ctl.LoadData <- R6Class("Ctl.LoadData")

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
	"",
	function(){

	}
)


