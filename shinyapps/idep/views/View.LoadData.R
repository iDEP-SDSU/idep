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
				label = TxtMsg$`1. Choose data type`,
				choiceNames = list(
					TxtMsg$`Read counts data (recommended)`,
					TxtMsg$`Normalized expression values (RNA-seq FPKM, microarray, etc.)`,
					TxtMsg$`Fold-changes and corrected P values from CuffDiff or any other program`), 
				choiceValues = list(
						1, 
						2,
						3
					),
				selected = 1
			),

			conditionalPanel(
				"input.dataFileFormat == 3",
				checkboxInput("noFDR", TxtMsg$`Fold-changes only, no corrected P values`, value = FALSE)
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
				TxtMsg$`2. Upload expression data (CSV or text)`,
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
			a(TxtMsg$`New! Analyze public RNA-seq data`, href="http://bioinformatics.sdstate.edu/reads/"),
			fileInput('file2', h5(TxtMsg$`Optional: Upload an experiment design file(CSV or text)`),
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
			strong(TxtMsg$`3. Verify guessed species. Change if neccessary.`),
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
			fileInput('gmtFile', TxtMsg$`Upload a geneset .GMT file for enrichment analysis (optional)`,
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

			p(HTML(paste0("<div align=\"right\"> <A HREF=\"javascript:history.go(0)\">", TxtMsg$Reset, "</A></div>" ))),

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




