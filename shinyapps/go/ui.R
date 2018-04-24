library(shiny,verbose=FALSE)
library(shinyBS,verbose=FALSE) # for popup figures
shinyUI(
  fluidPage(
    # Application title
    titlePanel("ShinyGO v0.4: Gene Ontology Enrichment Analysis & More "),
     sidebarLayout(
      sidebarPanel(
	  
	  				 fluidRow(
					column(6, actionButton("goButton", "Submit"))					
					,column(6, a(h5("Reset all",align = "right"), href="http://ge-lab.org:3838/go/") )					
				),
	  
	  
#	  a(h5("Reset all",align = "right"), href="http://ge-lab.org:3838/go/"), 
      tags$style(type="text/css", "textarea {width:100%}"),
      tags$textarea(id = 'input_text', placeholder = 'Just paste gene lists and click Submit. Most types of gene IDs accepted. Double check the guessed species,  adjust options and resubmit if needed. Any feedback is appreciated. ', rows = 8, ""),
 #      actionButton("goButton", "Submit"),
      h6(" "),
      htmlOutput("selectGO1"),
	  #htmlOutput("selectGO1_test"),
      selectInput("selectOrg", label = NULL,"Best matching species",width='100%'),  
      numericInput("minFDR", label = h5("P-value cutoff (FDR)"), value = 0.05),
     tags$style(type='text/css', "#minFDR { width:100%;   margin-top:-15px}"),  
     radioButtons("radio", label = NULL, choices = list("Upstream 300bp as promoter" = 300, "Upstream 600bp as promoter" = 600),selected = 300),
     tableOutput('species' ),
      h5( textOutput("text1") )
      ), # sidebarPanel
     mainPanel(
       tabsetPanel(
        tabPanel("Enrichment", tableOutput('EnrichmentTable'), 
		    conditionalPanel("input.goButton != 0", 
		                     downloadButton('downloadEnrichment', 'Download') )
            ,br(),br()	
			,conditionalPanel("input.selectGO == 'KEGG'", 
				htmlOutput('listSigPathways')	
				,imageOutput("KeggImage", width = "100%", height = "100%") )
				
			  
		    )
        #,tabPanel("Details", tableOutput('tableDetail')  )
		,tabPanel("Tree", downloadButton('GOTermsTree4Download','Figure' ),plotOutput('GOTermsTree')   )
		,tabPanel("Network", plotOutput('enrichmentNetworkPlot')   )

        ,tabPanel("Genes", tableOutput("conversionTable"), downloadButton('downloadGeneInfo', 'Download')  )
        ,tabPanel("Groups", tableOutput("grouping"), downloadButton('downloadGrouping', 'Download')   )          
        ,tabPanel("Plots", plotOutput("genePlot")  )
        , tabPanel("Genome", plotOutput("genomePlot", width = "100%")  )
        ,tabPanel("Promoter", tableOutput("promoter"), downloadButton('downloadPromoter', 'Download')   )  
		,tabPanel("STRING API", 
						h5("ShinyGO tries to match your species with the 115 archaeal, 1678 bacterial, 
						and 238 eukaryotic species in the",
							a(" STRING server", href="https://string-db.org/",target="_blank"),
								" and send the genes. If it is running, please wait until it finishes. This can take 5 minutes, especially for the first time when iDEP downloads large annotation files.")
						,htmlOutput("STRINGDB_species_stat") 
						,tags$head(tags$style("#STRINGDB_species_stat{color: blue;font-size: 15px;}"))						
						, selectizeInput('speciesName', label=NULL,choices = " ",
							   multiple = TRUE, options = list(maxItems = 1,							 
								 placeholder = 'Species name (e.g. Homo sapiens)',
								 onInitialize = I('function() { this.setValue(""); }')
							   )
							 )
						,textOutput('STRINGDB_mapping_stat')
						,tags$head(tags$style("#STRINGDB_mapping_stat{color: blue;font-size: 15px;}"))
						,br()
				,actionButton("ModalPPI","PPI network of DEGs"),br(),br()
						,selectInput("STRINGdbGO", label="Functional Enrichment", choices = list("GO Biological Process"= "Process"
																			  ,"GO Cellular Component"= "Component"
																			  ,"GO Molecular Function"= "Function"
																			  ,"KEGG" = "KEGG"
																			  ,"Pfam" = "Pfam"
																			  ,"InterPro" = "InterPro")
																, selected = "Process")							
						#,actionButton("submit2STRINGdb", "Submit")
						,downloadButton("STRING_enrichmentDownload")
						,tableOutput("stringDB_GO_enrichment")		
		
		)	
		,tabPanel("?", 	      h5("Based on annotation of 163 animal and 45 plant genomes in Ensembl BioMart as of 12/15/2017."
            ," Additional  data from",a("MSigDB (human),", href="https://doi.org/10.1093/bioinformatics/btr260",target="_blank") 
         ,a("GSKB (mouse)", href="http://biorxiv.org/content/early/2016/10/24/082511",target="_blank") 
         ,"and",a("  araPath (arabidopsis).", href="https://doi.org/10.1093/bioinformatics/bts421",target="_blank") 
         ," For feedbacks or data contributions, please"
         ,a("contact us, ",href="mailto:xijin.ge@sdstate.edu?Subject=ShinyGO" )
         , "or visit our",a(" homepage.", href="http://ge-lab.org/",target="_blank")
		 ,br(),br()
		 , a("Full list of supported organisms", href="https://idepsite.wordpress.com/species/",target="_blank")
		 ,br(),br()
		 , "ShinyGO Shares many functionalities and databases with ", a("iDEP.", href="http://ge-lab.org/idep/",target="_blank")
 		 ,br(),br()
		 ,h4("Changes:")
		 ,"4/24/2018: V0.4 Add STRING API, KEGG diagram, tree display and network."
	) ) 	
       ) #tabsetPanel
	   ,bsModal("ModalExamplePPI", "Protein-protein interaction(PPIs) networks ", "ModalPPI", size = "large"
		,h5("By sending your genes to the STRING website, 
			shinyGO is retrieving a sub-network, calculating PPI enrichment, 
		  and generating custom URLs to the STRING website containing your genes. This can take 5 minutes. Patience will pay off! ")
		,sliderInput("nGenesPPI", label = h5("Genes to include:"), min = 0, max = 400, value = 20,step=10) 
		,htmlOutput("stringDB_network_link")
		,tags$head(tags$style("#stringDB_network_link{color: blue; font-size: 15px;}"))

		,plotOutput("stringDB_network1")		 	   
   
     )# bsModal
     ) # mainPanel
    ) #sidebarLayout
    ,tags$head(includeScript("google_analytics.js")) # tracking usage
  ) #fluidPage
)



