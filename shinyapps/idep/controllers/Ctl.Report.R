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
			htmlinput = input,
			ReactVars = ReactVars,
			PreProcessResult = PreProcessResult,
			PreprocessSampleInfoResult = PreprocessSampleInfoResult,
			ConvertedIDResult = ConvertedIDResult,
			AllGeneInfo = AllGeneInfo,
			ConvertedPvals = ConvertedPvals,
			ConvertedRawReadcountData = ConvertedRawReadcountData,
			ConvertedTransformedData = ConvertedTransformedData
		)
		self$CleanImages()
		self$SaveRequiredPlotsInReportFolder(params)
		rmarkdown::render(
			"reports/main.Rmd", 
			output_file = file,
			output_format = input$selectReportOutputFormat,
			params = params
		)
		self$CleanImages()
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

Ctl.Report$set("public", "CleanImages",
	function(){
		files <- dir( "reports/img/", pattern = ".png$", full.names = TRUE )
		if(length(files) != 0)
		file.remove(files)
	}
)

Ctl.Report$set("public", "SaveRequiredPlotsInReportFolder",
	function(params){
		ReactVars <- params$ReactVars
		PreProcessResult <- params$PreProcessResult
		PreprocessSampleInfoResult <- params$PreprocessSampleInfoResult
		ConvertedIDResult <- params$ConvertedIDResult
		AllGeneInfo <- params$AllGeneInfo
		ConvertedPvals <- params$ConvertedPvals
		ConvertedRawReadcountData <- params$ConvertedRawReadcountData
		ConvertedTransformedData <- params$ConvertedTransformedData
		input <- params$htmlinput

		png("reports/img/mainheatmap.png", width = 800, height = 900)
		HeatmapCtrl$GetMainHeatmap(input, ReactVars, PreProcessResult, PreprocessSampleInfoResult)
		dev.off()



	}
)






