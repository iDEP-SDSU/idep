library('R6')
library(shiny)
library(shinyBS)
library(plotly)


Ctl.Report <- R6Class("Ctl.Report")



Ctl.Report$set("public", "SaveReportInTempFile",
	function(
		file,
		input,
		ReactVars,
		PreProcessResult,
		PreprocessSampleInfoResult,
		ConvertedIDResult,
		AllGeneInfo,
		ConvertedPvals,
		ConvertedRawReadcountData,
		ConvertedTransformedData
	){
		params <- list(
			htmlInput = input,
			ReactVars = ReactVars,
			PreProcessResult = PreProcessResult,
			PreprocessSampleInfoResult = PreprocessSampleInfoResult,
			ConvertedIDResult = ConvertedIDResult,
			AllGeneInfo = AllGeneInfo,
			ConvertedPvals = ConvertedPvals,
			ConvertedRawReadcountData = ConvertedRawReadcountData,
			ConvertedTransformedData = ConvertedTransformedData
		)

		rmarkdown::render(
			"reports/main.Rmd", 
			output_file = file,
			output_format = input$selectReportOutputFormat,
			params = params
		)
	}
)

Ctl.Report$set("public", "GetReportFileName",
	function(input){
		fileType = gsub( "_document", "", input$selectReportOutputFormat )
		
		if(fileType == 'html'){
			return("iDepReport.html")
		}

		if(fileType == 'pdf'){
			return("iDepReport.pdf")
		}

		return("iDepReport.docx")
	}
)








