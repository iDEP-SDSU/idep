library(shiny,verbose=FALSE)
library(shinyBS,verbose=FALSE) # for popup figures
shinyUI(
  fluidPage(
    # Application title
     sidebarLayout(
	 
      sidebarPanel(
	    titlePanel("ShinyGO v0.41: Gene Ontology Enrichment Analysis + more"),  
	  	p(HTML("<div align=\"right\"> <A HREF=\"javascript:history.go(0)\">Reset</A></div>" )),					
      tags$style(type="text/css", "textarea {width:100%}"),
      tags$textarea(id = 'input_text', placeholder = 'Just paste gene lists and click Submit. Most types of gene IDs accepted. Double check the guessed species, and adjust options if needed. Any feedback is appreciated. ', rows = 8, ""),
       actionButton("goButton", "Submit"),
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
        tabPanel("Enrichment" 
			,conditionalPanel("input.goButton == 0 "  # welcome screen
				,br(),br(),h4("Welcome to ShinyGO! Just paste your gene list to get enriched GO terms and othe pathways for over 200 plant and animal species. In addition, it also produces
				KEGG pathway diagrams with your genes highlighted, hierarchical clustering trees and networks summarizing 
				overlapping terms/pathways, protein-protein interaction networks, gene characterristics plots, and enriched promoter motifs. See example outputs below:")			
				,br(),img(src='enrich.png', align = "center",width="660", height="339")
				,br(),br(),img(src='KEGG2.png', align = "center",width="541", height="360")
				,br(),br(),img(src='GOtree3.png', align = "center",width="500", height="258")
				,br(),br(),img(src='GOnetwork2.png', align = "center",width="500", height="248")
				,br(),br(),img(src='PPInetwork2.png', align = "center",width="500", height="391")	
				,br(),br(),img(src='chr.png', align = "center",width="444", height="338")			
				,br(),br(),img(src='promoter.png', align = "center",width="717", height="288")					
			)

			,tableOutput('EnrichmentTable')	
		    ,conditionalPanel("input.goButton != 0", 
		        downloadButton('downloadEnrichment', 'Download')			
				,br(),br(),"Select KEGG pathways in the left to display pathway diagrams."
			)			
			,conditionalPanel("input.selectGO == 'KEGG'", 
				htmlOutput('listSigPathways')
				,br(),br(),imageOutput("KeggImage", width = "100%", height = "100%")				
				,br(),br()	 
			)
		) # enrichment tab
        #,tabPanel("Details", tableOutput('tableDetail')  )
		,tabPanel("Tree", downloadButton('GOTermsTree4Download','Figure' ),plotOutput('GOTermsTree')   )
		,tabPanel("Network", actionButton("layoutButton", "Change layout"),downloadButton('enrichmentNetworkPlotDownload')
		,plotOutput('enrichmentNetworkPlot')   )

        ,tabPanel("Genes", tableOutput("conversionTable"), downloadButton('downloadGeneInfo', 'Download')  )
        ,tabPanel("Groups", tableOutput("grouping"), downloadButton('downloadGrouping', 'Download')   )          
        ,tabPanel("Plots", plotOutput("genePlot", inline = TRUE,width='auto',height='auto')
			, plotOutput("genePlot2", inline = TRUE,width='auto',height='auto')  )
        #,tabPanel("Plots2", plotOutput("genePlot2")  )
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
		,tabPanel("?", 
		h5("Based on annotation of 163 animal and 45 plant genomes in Ensembl BioMart as of 12/15/2017."
            ," Additional  data from",a("MSigDB (human),", href="https://doi.org/10.1093/bioinformatics/btr260",target="_blank") 
         ,a("GSKB (mouse)", href="http://biorxiv.org/content/early/2016/10/24/082511",target="_blank") 
         ,"and",a("  araPath (arabidopsis).", href="https://doi.org/10.1093/bioinformatics/bts421",target="_blank") 
         ," For feedbacks or data contributions, please"
         ,a("contact us, ",href="mailto:xijin.ge@sdstate.edu?Subject=ShinyGO" )
         , "or visit our",a(" homepage.", href="http://ge-lab.org/",target="_blank")
		 ,br(),br()
		 , a("Full list of supported organisms", href="https://idepsite.wordpress.com/species/",target="_blank")
		 ,br(),br()
		 , "ShinyGO shares many functionalities and databases with ", a("iDEP.", href="http://ge-lab.org/idep/",target="_blank")
 		 ,br(),br()

		 #,includeHTML("help.htm")
		 ,h4("Input:")
		 ,"A list of gene ids, separated by tab, space, comma or new line."
		 ,h4("Output:")
		 ,"Enriched GO terms and pathways:"
		 ,br()
		 ,img(src='enrich.png', align = "center",width="660", height="339")		 
		 ,br(),br()		 
		 ,"In addition to the enrichment table, a set of plots are produced. If KEGG database is choosen, then enriched pathway diagrams are shown, with user's genes highlighted. Like this one below:"
		 ,br()
		 ,img(src='KEGG.png', align = "center",width="696", height="494")
		 ,br(),br()
		 ,"Many GO terms are related. Some are even redundant, like \"cell cycle\" and \"cell cycle process\". 
		 To visualize such relatedness in enrichment results, we use a hierarchical clustering tree and network. 
		 In this tree below, related GO terms are grouped together based on how many genes they share. The size of the solid circle corresponds to the enrichment FDR."
		 ,br()
		 ,img(src='GOtree.png', align = "center",width="700", height="361")
		 ,br(),br()
		 ,"In this network below, each node represent a enriched GO term. Related GO terms are connected by a line, whose thickness reflect percent of overlapping genes. Size of the node corresponds to number of genes." 
		 ,br(),img(src='GOnetwork.png', align = "center",width="500", height="248")
		 ,br(),br()		 
		 ,"Through API access to STRING-db, we also retrieve protein-protein interaction (PPI) network. In addition to a static network image, users can also get access to an interactive graphics at the www.string-db.org web server."
		 ,br(),img(src='PPInetwork.png', align = "center",width="700", height="547")
		 ,br(),br()	
		 ,"ShinyGO also detects transcription factor (TF) binding motifs enriched in the promoters of user's genes."
		,br(),br(),img(src='promoter.png', align = "center",width="717", height="288")			 
		 
		 
		 ,br(),h4("Changes:")
		 ,"4/24/2018: V0.4 Add STRING API, KEGG diagram, tree display and network."	
		 ,"4/27/2018: V0.41 Change to ggplot2, add grid and gridExtra packages"
		 ,"4/30/2018: V0.42 changed figure configurations for tree."
		 
		 ) ) 	
       ) #tabsetPanel
	   ,bsModal("ModalExamplePPI", "Protein-protein interaction(PPIs) networks ", "ModalPPI", size = "large"
		,h5("By sending your genes to the STRING website, 
			shinyGO is retrieving a sub-network, calculating PPI enrichment, 
		  and generating custom URLs to the STRING website containing your genes. This can take 5 minutes. Patience will pay off! ")
		,sliderInput("nGenesPPI", label = h5("Genes to include:"), min = 0, max = 400, value = 50,step=10) 
		,htmlOutput("stringDB_network_link")
		,tags$head(tags$style("#stringDB_network_link{color: blue; font-size: 15px;}"))

		,plotOutput("stringDB_network1")		 	   
   
     )# bsModal
     ) # mainPanel
    ) #sidebarLayout
    ,tags$head(includeScript("google_analytics.js")) # tracking usage
  ) #fluidPage
)



