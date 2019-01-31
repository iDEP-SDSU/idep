library('R6')

View.LoadData <- R6Class("View.LoadData")

# side bar support function: sideLoadDemoDataSection
View.LoadData$set(  
	"public",
	"sideLoadDemoDataSection",
	function(){
		wellPanel(				   
			actionButton("goButton", TxtMsg$`Click here to load demo data`),
			tags$head(tags$style("#goButton{color: red;
								font-size: 16px;
								font-style: italic;
								}")),					
			h5(TxtMsg$` and just click the tabs for some magic!`, style = "color:red")
		)
	}
)

# side bar support function: sideChooseDataTypeSection
View.LoadData$set(  
	"public",
	"sideChooseDataTypeSection",
	function(){
		wellPanel( 
			radioButtons(   
				"dataFileFormat", 
				label = "1. Choose data type", 
				choices = list(
						"Read counts data (recommended)" = 1, 
						"Normalized expression values (RNA-seq FPKM, microarray, etc.)" = 2,
						"Fold-changes and corrected P values from CuffDiff or any other program" = 3
					),
				selected = 1
			),

			conditionalPanel(
				"input.dataFileFormat == 3",
				checkboxInput("noFDR", "Fold-changes only, no corrected P values", value = FALSE)
			)
		)
	}
)

# side bar support function: sideUploadFileSection
View.LoadData$set(  
	"public",
	"sideUploadFileSection",
	function(){
		wellPanel(
			fileInput(
				'file1', 
				'2. Upload expression data (CSV or text)',
				accept = c(
					'text/csv',
					'text/comma-separated-values',
					'text/tab-separated-values',
					'text/plain',
					'.csv',
					'.tsv'		  
				) 
			)
		)	
	}
)

# side bar support function: sideUploadDataFromPublicSource
View.LoadData$set(  
	"public",
	"sideUploadDataFromPublicSource",
	function(){
		wellPanel(
			a("New! Analyze public RNA-seq data", href="http://bioinformatics.sdstate.edu/reads/"),
			fileInput('file2', h5('Optional: Upload an experiment design file(CSV or text)'),
				accept = c(
					'text/csv',
					'text/comma-separated-values',
					'text/tab-separated-values',
					'text/plain',
					'.csv',
					'.tsv'		  
				)
			)
		)
	}
)

# side bar support function: sideGuessSpeciesSection
View.LoadData$set(  
	"public",
	"sideGuessSpeciesSection",
	function(){
		wellPanel(
			strong("3. Verify guessed species. Change if neccessary."),
			selectInput("selectOrg", label = NULL,"Best matching species",width='100%')
		)
	}
)

# side bar support function: side.CondPanel.GuessSpeciesSection
View.LoadData$set(  
	"public",
	"side.CondPanel.GuessSpeciesSection",
	function(){
		conditionalPanel("input.selectOrg == 'NEW'",
			fileInput('gmtFile', 'Upload a geneset .GMT file for enrichment analysis (optional)',
				accept = c(
						'text/csv',
						'text/comma-separated-values',
						'text/tab-separated-values',
						'text/plain',
						'.csv',
						'.tsv',
						'.gmt'
				)
			)
		)
	}
)


# Main side bar
View.LoadData$set(  
	"public", 
	"sidebarLayout", 
	function(){
		sidebarPanel(

			self$sideLoadDemoDataSection(),

			p(HTML("<div align=\"right\"> <A HREF=\"javascript:history.go(0)\">Reset</A></div>" )),

			self$sideChooseDataTypeSection(),

			self$sideUploadFileSection(),

			self$sideUploadDataFromPublicSource(),

			self$sideGuessSpeciesSection(),

			self$side.CondPanel.GuessSpeciesSection(),

			tableOutput('species'),
			
			a( h5("?",align = "right"), href="https://idepsite.wordpress.com/data-format/",target="_blank")
		)
	}
)




