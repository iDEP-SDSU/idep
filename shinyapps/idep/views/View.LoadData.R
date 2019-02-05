library('R6')

View.LoadData <- R6Class("View.LoadData")

# side bar: sideLoadDemoDataSection
View.LoadData$set(  
	"public",
	"sideLoadDemoDataSection",
	function(){
		wellPanel(				   
			actionButton("goButton", TxtMsg$`Click here to load demo data`),
			actionButton("tabBut2", "View Search Page"),
			tags$head(tags$style("#goButton{color: red;
								font-size: 16px;
								font-style: italic;
								}")),					
			h5(TxtMsg$` and just click the tabs for some magic!`, style = "color:red")
		)
	}
)

# side bar: sideChooseDataTypeSection
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

# side bar: sideUploadFileSection
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

# side bar: sideUploadDataFromPublicSource
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

# side bar: sideGuessSpeciesSection
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

# side bar: side.CondPanel.GuessSpeciesSection
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

View.LoadData$set(
	"public",
	"side.CondPanel.PublicSampleId",
	function(){
		wellPanel(textOutput('selectedSampleIds'))
	}
)

# popout: popSearchSeqData
View.LoadData$set(
	"public",
	"popSearchSeqData",
	function(){
		bsModal(
			"modalSearchTab", 
			"Public RNA-Seq and ChIP-Seq data", 
			"tabBut2", 
			size="large",
			fluidPage(
				h5("Download gene-level read counts data for 7,793 human and mouse datasets from",
					a("ARCHS4 ", href="http://amp.pharm.mssm.edu/archs4/help.html"),
					"on Nov. 5, 2018.",
					"HiSeq 2000 and HiSeq 2500 raw data from NCBI's SRA database are 
					aligned with Kallisto against human GRCh38 or mouse GRCm38 cDNA reference."
				),
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
					column(4, 	textOutput('selectedDataset') ),
	               	column(3, 	actionButton('downloadSearchedData', 'Download') ),
	               	column(5, 	h5(	"Edit and uplodad file to ", 
									a("iDEP", href="http://bioinformatics.sdstate.edu/idep/"), 
									"for analysis"
								)
					)	
	     		),
				br(),
				tableOutput('samples'),
				h4("Loading data and R packages ... ..."),
				htmlOutput('DoneLoading') 
			)	
		)
	}
)



# Main side bar
View.LoadData$set(  
	"public", 
	"sidebarPanel", 
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

			self$side.CondPanel.PublicSampleId(),
			
			a( h5("?",align = "right"), href="https://idepsite.wordpress.com/data-format/",target="_blank")

			
		)
	}
)

# Main Panel
View.LoadData$set(
	"public",
	"mainPanel",
	function(){
		mainPanel(
    		tableOutput('sampleInfoTable'),
			tableOutput('contents'),
    		#,h3("Service will not be available starting 6:30 am (US central time) on June 21 (Friday) 
    		#due to scheduled maintenance. It should take less than 45 minutes. ",  style = "color:red")
    		#,h3("We will re-sbumit our grant proposal to NIH. If you didn't send us a support letter last time, 
    		#    please consider sending us a brief email/letter before Nov. 15th, with your 
    		#    broad area of research and how iDEP helps your work. Thanks!"
    		#,a("Email",href="mailto:Xijin.Ge@SDSTATE.EDU?Subject=iDEP suggestions"), style = "color:red"),
    		
			br(),
			
			img(src='flowchart.png', align = "center",width="562", height="383"),
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
    		), #h5
    		h3("Loading R packages ... ..."),
    		htmlOutput('fileFormat'),

			self$popSearchSeqData()
    	)
	}
)


