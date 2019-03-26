library('R6')

View.Report <- R6Class("View.Report")




###############################################################################
########################		Main Structure		###########################
###############################################################################

# An empty side panel 
View.Report$set("public", "sidebarPanel",
	function(){
		sidebarPanel(
			"Nothing here"
		)
	}
)

# A main panel with
#		selections of sections
#		selections of output file format
View.Report$set("public", "mainPanel",
	function(){
		mainPanel(
			selectInput("selectReportOutputFormat", "Select Report Format",
				choices = list("html_document", "pdf_document", "word_document")
			),
			checkboxInput("isReportIncludeDataDescription", "Include Data Description", value=TRUE),
			checkboxInput("isReportIncludePreProcess", "Include Data Preprocess", value=TRUE),
			downloadButton('download_Report_MainReport', "Download Report")
		)
	}
)
