library('R6')

View.LoadData <- R6Class("View.LoadData")


###############################################################################
########################		Main Structure		###########################
###############################################################################

# Side Bar
View.LoadData$set( "public", "sidebarPanel", 
	function(){
		sidebarPanel(

			actionButton("btn_LoadData_DemoData", TxtLibrary$btn_label_LoadDemoData),

			div(actionButton("btn_Reset_Data", TxtLibrary$btn_label_ResetData, onclick="javascript:history.go(0)"),
				style="float:right"),

						
			br(),
			br(),
			
			# shows when data source is not 'demo'
			self$side_ConPanel_NonDemoDataBuild(),

			# shows when data source is 'demo'
			self$side_ConPanel_ExplainOfDemoData(),

#			self$sideGuessSpeciesSection(),
#
#			self$side.CondPanel.GuessSpeciesSection(),
#
#			tableOutput('species'),
			
			a( h5("?",align = "right"), href="https://idepsite.wordpress.com/data-format/",target="_blank")

			
		)
	}
)

# Main Panel
View.LoadData$set("public", "mainPanel",
	function(){
		mainPanel(
    		tableOutput('tbl_TestDesign'),
			tableOutput('tbl_RawDataTop20'),
    		
			br(),
			
			img(src='flowchart.png', align = "center",width="562", height="383"),
    		self$main_UpdateNote(),
			conditionalPanel( condition = "!output.txt_PackageLoadedMessage",
    			h3("Loading R packages ... ...")
			),
			
    		htmlOutput('txt_PackageLoadedMessage'),

			self$Pop_DownloadPublicData()
    	)
	}
)

###############################################################################
###################				Components				#######################
###############################################################################



###################				NonDemoDataBuild			###################
View.LoadData$set("public", "side_ConPanel_NonDemoDataBuild",
	#	Shows when output.DataSource != Demo
	#	This panel helps user to build up their own data:
	#		1. User could upload data from their PC, or search and use public
	#			dataset
	#		2. User could upload their own experiment design
	function(){
		conditionalPanel( condition = "output.DataSource != 'Demo'",
			tags$b("Build your own data set"),
			
			br(),

			self$side_ChooseDataTypeSection(),

			self$side_ChooseNonDemoData(),

			self$side_UploadExperimentDesign()
		)
	}
)

###################				ExplainOfDemoData			###################
View.LoadData$set("public", "side_ConPanel_ExplainOfDemoData",
	#	Shows when output.DataSource == Demo
	#	A brief description about the demo data.
	function(){
		conditionalPanel( condition = "output.DataSource == 'Demo'",
			h5('Some description of our demo data')
		)
	}
)

# side bar: side_ChooseDataTypeSection
View.LoadData$set(  
	"public",
	"side_ChooseDataTypeSection",
	function(){
		wellPanel(
			radioButtons( 
				"dataFileFormat", 
				label = '1. Choose Data Type',
				choiceNames = list(
					TxtLibrary$`Read counts data (recommended)`,
					TxtLibrary$`Normalized expression values (RNA-seq FPKM, microarray, etc.)`,
					TxtLibrary$`Fold-changes and corrected P values from CuffDiff or any other program`
				), 
				choiceValues = list(
					1, 
					2,
					3
				),
				selected = 1
			),
			conditionalPanel(
				"input.dataFileFormat == 3",
				checkboxInput("noFDR", TxtLibrary$`Fold-changes only, no corrected P values`, value = FALSE)
			)
		)
	}
)


# side bar: side_ChooseNonDemoData
View.LoadData$set(  
	"public",
	"side_ChooseNonDemoData",
	function(){
		wellPanel(
			fileInput(
				'fileUploadedData', 
				label = '2. Upload expression data (CSV or text)',
				accept = c(
					'text/csv',
					'text/comma-separated-values',
					'text/tab-separated-values',
					'text/plain',
					'.csv',
					'.tsv'		  
				) 
			),
			tags$b('Or, Search public data'),
			actionButton("btn_LoadData_SearchDataFromPublic", TxtLibrary$btn_label_SearchPublicData)
		)	
	}
)

# side bar: test design
View.LoadData$set(  
	"public",
	"side_UploadExperimentDesign",
	function(){
		wellPanel(
			fileInput('fileUploadedTestDesign', 
				label = h5(TxtLibrary$`Optional: Upload an experiment design file(CSV or text)`),
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

# side bar: sideGuessSpeciesSection
View.LoadData$set(  
	"public",
	"sideGuessSpeciesSection",
	function(){
		wellPanel(
			strong(TxtLibrary$`3. Verify guessed species. Change if neccessary.`),
			selectInput("selectOrg", label = NULL,"Best matching species",width='100%')
		)
	}
)

# side bar: side_CondPanel_GuessSpeciesSection

######## This part may need move to pre process
View.LoadData$set(  
	"public",
	"side_CondPanel_GuessSpeciesSection",
	function(){
		conditionalPanel("input.selectOrg == 'NEW'",
			fileInput('gmtFile', TxtLibrary$`Upload a geneset .GMT file for enrichment analysis (optional)`,
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

# popout: Pop_DownloadPublicData
View.LoadData$set("public", "Pop_DownloadPublicData",
	function(){
		bsModal(
			"modalSearchPublicDataTab", 
			"Public RNA-Seq and ChIP-Seq data", 
			"btn_LoadData_SearchDataFromPublic", 
			size="large",
			fluidPage(			
				fluidRow( 
					column(	9, 
							h5("Search and click on a dataset to see more information before download.") 
					),
        	       	column(	3, 
							selectInput(
								"selected.species.archs4", 
								"", 
        	        			choices = list("Human"= "human", "Mouse" = "mouse"), 
								selected = "human"
							)
					)
				),
				DT::dataTableOutput('SearchData'),
				HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />'),
				br(),
				fluidRow( 	
					column(4, 	textOutput('selectedDataset')),
	               	column(3, 	actionButton('btn_LoadData_UseSelectedPublicData', 'Use This Sample') )	
	     		),
				br(),
				tableOutput('ViewData_SelectFromPublic_Samples'),
				h4("Loading data and R packages ... ..."),
				htmlOutput('DoneLoading')
			)
		)
	}
)

# main panel: Update Notes
View.LoadData$set(
	"public",
	"main_UpdateNote",
	function(){
		wellPanel(
			h5("v0.81 Enabled downloading of publication-ready vector graphics files"),
    		h5("New v0.80  Updated annotation database. Comprehensive pathway 
    		    database for human. TF binding motifs for 200+ speceis. Old version made available."
    		),
    		h5("New v0.70  iDEP generates R and R Markdown codes for users to run in stand-alone!"),
    		a("R Markdown example.",align = "left", href="http://rpubs.com/ge600/R",target="_blank"),
    		h5("New v0.68! Try the STRING-db API access on the DEG2 page that offer 
    		    protein interaction networks and GO enrichment for thousands species, including bacteria."),
    		h5(	"Integrated Differential Expression and Pathway analysis (iDEP) of transcriptomic data. See ",
    		  	a(" documentation", href="https://idepsite.wordpress.com/", target="_blank"), 
    		  	"and",
    		  	a(" manuscript.", href="https://www.biorxiv.org/content/early/2018/04/20/148411",target="_blank"),
    		  	"Based on annotation of",
    		  	a( "167 animal and 53 plant genomes ",href="https://idepsite.wordpress.com/species/",target="_blank"), 
    		  	"in Ensembl BioMart as of 12/15/2017.",
    		  	a("STRING-db ", href="https://string-db.org/",target="_blank"),
    		  	"offer API access to protein interaction networks and annotations 
    		  	for 115 archaeal, 1678 bacterial, and 238 eukaryotic species.",
    			" Additional  data from",
    			a("KEGG, ", href="www.genome.jp/kegg/",target="_blank"),
    			a("Reactome, ", href="http://www.reactome.org/",target="_blank"),
    			a("MSigDB (human),", href="https://doi.org/10.1093/bioinformatics/btr260",target="_blank"),
    			a("GSKB (mouse)", href="http://biorxiv.org/content/early/2016/10/24/082511",target="_blank"), 
    			"and",
    			a("  araPath (arabidopsis).", href="https://doi.org/10.1093/bioinformatics/bts421", target="_blank"),
    			" For feedbacks or data contributions (genes and GO mapping of any species), please",
    			a("contact us, ",href="mailto:xijin.ge@sdstate.edu?Subject=iDEP" ),
    			"or visit our",a(" homepage.", href="http://ge-lab.org/",target="_blank"),
    			"Send us suggestions or any error message to help improve iDEP.",
    			a("Email",href="mailto:Xijin.Ge@SDSTATE.EDU?Subject=iDEP suggestions")
    		)
		)
	}
)


