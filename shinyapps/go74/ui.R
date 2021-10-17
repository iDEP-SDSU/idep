###################################################
# Author: Steven Ge Xijin.Ge@sdstate.edu
# co-author: Eric Tulowetzke, eric.tulowetzke@jacks.sdstate.edu
# Lab: Ge Lab
# R version 4.0.5
# Project: ShinyGO v65
# File: ui.R
# Purpose of file:ui logic of app
# Start data: NA (mm-dd-yyyy)
# Data last modified: 06-22-2021
#######################################################
library(shiny,verbose=FALSE)
library(shinyBS,verbose=FALSE) # for popup figures
library(plotly) # interactive network plot
library(visNetwork)
library('shinyjs', verbose = FALSE)
library('reactable', verbose = FALSE)
ui <- fluidPage(
  # Application title
  sidebarLayout(
    
    sidebarPanel(
      titlePanel("ShinyGO v0.74: Gene Ontology Enrichment Analysis + more"),  
      
      fluidRow(
        column(8,   actionButton("useDemo1", "Demo genes"),	  	  ),
        #column(4,   actionButton("useDemo2", "Demo 2"),	  	  ),
        column(4,  p(HTML("<div align=\"right\"> <A HREF=\"javascript:history.go(0)\">Reset</A></div>" ))	  	  )
      ),
      tags$style(type="text/css", "textarea {width:100%}"),	  	
                
      tags$textarea(id = 'input_text', placeholder = 'Just paste a list of genes and click Submit. More adjustments below. Most types of gene IDs accepted. Double check the guessed species, and adjust if needed. ', rows = 8, ""),
      
      
      fluidRow(
        column(8,  actionButton("backgroundGenes", "Background (recommended)")	),
        column(4,   actionButton("goButton", strong("Submit"))	  	  )
      ),
      br(),
      htmlOutput("selectGO1"),
      
      fluidRow(
        column(6,
               numericInput(inputId = "minFDR",
                            label = h5("P-value cutoff (FDR)"),
                            value = 0.05, step = 0.01)
        ),
        column(6,
               selectInput("maxTerms", h5("# of top pathways to show"),
                           choices = list("10" = 10,
                                          "20" = 20,
                                          "30" = 30,
                                          "40" = 40,
                                          "50" = 50,
                                          "100" = 100,
                                          "500" = 500),
                           selected = "30")
        )
      ),
      #tags$style(type='text/css', "#minFDR { width:100%;   margin-top:-15px}"),  
      # selectInput("selectOrg", label = NULL,"Best matching species",width='100%'),  
      
      h5("Select or search your species. Or use our best guess."),
      selectizeInput('selectOrg', 
                     label    = NULL,
                     choices  = " ",
                     multiple = TRUE,
                     options  = list( maxItems     = 1,               
                                      placeholder  = 'Best matching species',
                                      onInitialize = I('function() { this.setValue(""); }'))  
                     #,selected = "Best matching species"                                                  
      ), 
      tableOutput('species' ),
      #h5("Check this out if you want example of our gene ids,
      #     or download gene mapping."),
      actionButton(inputId = "geneIdButton",
                   label =  "Gene ID Examples")
    ), # sidebarPanel
    mainPanel(
      tabsetPanel(
        tabPanel("Enrichment" 
                 ,conditionalPanel("input.goButton == 0 "  # welcome screen

                                   ,h4("Oct. 15, 2021: Version 0.74. Database updated to Ensembl Release 104 and STRING v11. We now recommends the use of background genes in enrichment analysis. V.0.74 is much faster with even large set of background genes.")
                                   ,h4("We recently hired Jenny for database updates and user support.",
                                       a("Email Jenny ",href="mailto:gelabinfo@gmail.com?Subject=ShinyGO"),
                                       "for questions, suggestions or data contributions.") 

                                   
                                   #,h4("If your gene IDs are not recognized, please let us know. We might be able to add customized gene mappings to Ensembl gene IDs.")
                                   
                                   ,h4("2/3/2020: Now published by", a("Bioinformatics.", href="https://doi.org/10.1093/bioinformatics/btz931",target="_blank"))
                                   ,h4("Just paste your gene list to get enriched GO terms and othe pathways for over 400 plant and animal species, 
				    based on annotation from Ensembl , Ensembl plants and Ensembl Metazoa. An additional 5000 genomes 
				    (including bacteria and fungi) are   annotated based on STRING-db (v.11). In addition, it also produces
				    KEGG pathway diagrams with your genes highlighted, hierarchical clustering trees and networks summarizing 
				    overlapping terms/pathways, protein-protein interaction networks, gene characterristics plots, and enriched promoter motifs. 
                 See example outputs below:")			
                                   ,br(),img(src='enrich.png', align = "center",width="660", height="339")
                                   ,br(),br(),img(src='KEGG2.png', align = "center",width="541", height="360")
                                   ,br(),br(),img(src='GOtree3.png', align = "center",width="500", height="258")
                                   ,br(),br(),img(src='GOnetwork2.png', align = "center",width="500", height="248")
                                   ,br(),br(),img(src='PPInetwork2.png', align = "center",width="500", height="391")	
                                   ,br(),br(),img(src='chr.png', align = "center",width="444", height="338")			
                                   ,br(),br(),img(src='promoter.png', align = "center",width="717", height="288")					
                 )
                 ,div(style="display:inline-block", 
                          selectInput(inputId = "SortPathways",
                              label = NULL,
                              choices = c("Sort by FDR", "Sort by Fold Enriched", "Sort by Genes", "Sort by Category Name" ),
                              selected = "Sort by Fold Enriched" ),
                     style="algn:right")
                 ,tableOutput('EnrichmentTable')	
                 ,conditionalPanel("input.goButton != 0", 
                                   downloadButton('downloadEnrichment', 'Download table with gene IDs')			
                                   ,br(),br()
                 )			
                 ,conditionalPanel("input.selectGO == 'KEGG'", 
                                   htmlOutput('listSigPathways')
                                   ,br(),br(),imageOutput("KeggImage", width = "100%", height = "100%")				
                                   ,br(),br()	 
                 )
                 ,h5("Enrichment analysis based on hypergeometric distribution followed by FDR correction.")
                 ,h5("Select KEGG pathways in the left to display your genes in pathway diagrams.")
        ) # enrichment tab
        #,tabPanel("Details", tableOutput('tableDetail')  )
        ,tabPanel("Tree" 
                  ,h5("A hierarchical clustering tree summarizing the correlation among significant pathways 
       listed in the Enrichment tab. Pathways with many shared genes are clustered together. 
       Bigger dots indicate more significant P-values.")
                  ,downloadButton('GOTermsTree4Download','Figure' )
                  ,plotOutput('GOTermsTree')
        )
        
        ,tabPanel("Network" 
                  ,fluidRow(
                    column(2, actionButton("layoutButton", "Change layout") ),
                    column(2, actionButton("GONetwork", "Static plot") ),
                    column(2, h5("Edge cutoff:"), align="left" ) ,
                    column(2, numericInput("edgeCutoff", label = NULL, value = 0.20, min = 0, max = 1, step = .1), align="right"  ), 
                    column(2, checkboxInput("wrapTextNetwork", "Wrap text", value = FALSE)) 
                  ),      
                  visNetworkOutput("enrichmentNetworkPlotInteractive",height = "800px", width = "800px"),
                  downloadButton("enrichmentNetworkPlotInteractiveDownload","Download HTML"),
                  downloadButton("downloadEdges", "Edges"),
                  downloadButton("downloadNodes", "Nodes"),
                  h5("Similar to the Tree tab, this interactive plot also shows the relationship between enriched pathways. 
       Two pathways (nodes) are connected if they share 20% (default) or more genes. 
       You can move the nodes by dragging them, zoom in and out by scrolling, 
       and shift the entire network by click on an empty point and drag. 
       Darker nodes are more significantly enriched gene sets. 
       Bigger nodes represent larger gene sets.  
       Thicker edges represent more overlapped genes.")               		   
        )
        
        ,tabPanel("Genes"
                  , tableOutput("conversionTable")
                  , downloadButton('downloadGeneInfo', 'Download')
                  , h5("This table shows how your genes are converted to Ensembl gene IDs.")          
        )
        
        ,tabPanel("Groups"
                  , tableOutput("grouping")
                  , downloadButton('downloadGrouping', 'Download')
                  , h5("If Gene Ontology is chosen, your genes are grouped by functional categories defined by high-level GO terms. ")          
        ) 
        
        ,tabPanel("Plots"
                  , h5("The characteristics of your genes are compared with the rest in the genome. Chi-squared and Student's 
              t-tests are run to see if your genes have special characteristics when compared with all the other genes or, if uploaded, a customized background.")
                  , plotOutput("genePlot", inline = TRUE,width='auto',height='auto')
                  , plotOutput("genePlot2", inline = TRUE,width='auto',height='auto')  
        )
        
        ,tabPanel("Genome"
                  , h5("Your genes are marked in each of the chromosomes. 
                Note that the scale for each chromosomes are different.")
                  , plotOutput("genomePlot", width = "100%")  
                  
        )
        
        ,tabPanel("Promoter", 
                  radioButtons("radio", label = NULL, choices = list("Upstream 300bp as promoter" = 300, 
                                                                     "Upstream 600bp as promoter" = 600),selected = 300),
                  tableOutput("promoter"), 
                  downloadButton('downloadPromoter', 'Download'),
                  h5("The promoter sequences of your genes are compared with those of the 
              other genes in the genome in terms of transcription factor (TF) binding motifs. 
              \"*Query gene\" indicates a transcription factor coded by a gene included in 
              your list.")
                  
        )  
        ,tabPanel("STRING", 
                  h5("Your genes are sent to STRING-db website for enrichment analysis 
            and retrieval of a protein-protein network. We tries 
            to match your species with the archaeal, bacterial, 
						and eukaryotic species in the",
                     a(" STRING server", href="https://string-db.org/",target="_blank"),
                     " and send the genes. If it is running, please wait until it finishes. This can take 5 minutes, especially for the first time when shinyGO downloads large annotation files.")
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
                  ,selectInput("STRINGdbGO",
                               label="Functional Enrichment",
                               choices = list("GO Biological Process"= "Process",
                                            "GO Cellular Component"= "Component",
                                            "GO Molecular Function"= "Function",
                                            "KEGG" = "KEGG",
                                            "Pfam" = "Pfam",
                                            "InterPro" = "InterPro")
                               , selected = "Process")
                  ,downloadButton("STRING_enrichmentDownload")
                  ,tableOutput("stringDB_GO_enrichment")		
                  
        )	
        ,tabPanel("?"
                  ," For feedbacks, please"
                  ,a("contact us, ",href="mailto:gelabinfo@gmail?Subject=ShinyGO" )
                  , "or visit our",a(" homepage.", href="http://ge-lab.org/",target="_blank")
                  ,"For details, please see our", a("paper", href="https://doi.org/10.1093/bioinformatics/btz931", target="_blank")
                  ,"and a detailed", a("demo.", href="https://www.biorxiv.org/content/biorxiv/suppl/2018/05/04/315150.DC1/315150-1.pdf",target="_blank") 
                  , "ShinyGO shares many functionalities and databases with ", a("iDEP.", href="http://ge-lab.org/idep/",target="_blank")
                  ," Source code at", a(" GitHub. ", href="https://github.com/iDEP-SDSU/idep/tree/master/shinyapps/go",target="_blank")
                  ,br(),br()
                  ,strong("Citation:")
                  ,br()
                  ,"Ge SX, Jung D & Yao R,", a(" Bioinformatics 2020", href="https://doi.org/10.1093/bioinformatics/btz931", target="_blank")
                  ,br(),br()
                  ,strong("Previous versions (still functional):")
                  ,br()
                  ,a("ShinyGO V0.65, "
                     , href="http://bioinformatics.sdstate.edu/go65/")
                  ,"based on database derived from Ensembl Release 103, archived on Oct. 15, 2021"
                  ,br()
                  ,a("ShinyGO V0.61, "
                     , href="http://bioinformatics.sdstate.edu/go61/")
                  ,"based on database derived from Ensembl Release 96, archived on May 23, 2020"
                  ,br()		 
                  ,a("ShinyGO V0.60, "
                     , href="http://bioinformatics.sdstate.edu/go60/")
                  ,"based on database derived from Ensembl BioMart version 96, archived on Nov 6, 2019"
                  ,br()
                  ,a("ShinyGO V0.51, "
                     , href="http://bioinformatics.sdstate.edu/go51/")
                  ,"based on database derived from Ensembl BioMart version 95, archived on May 20, 2019"
                  ,br()
                  ,a("ShinyGO V0.50, "
                     , href="http://bioinformatics.sdstate.edu/go50/")
                  ,"based on database derived from Ensembl BioMart version 92, archived on March 29, 2019"
                  ,br()		
                  ,a("ShinyGO V0.41, "
                     , href="http://bioinformatics.sdstate.edu/go41/")
                  ,"based on database derived from Ensembl BioMart version 91, archived on July 11, 2018"
                  ,br()
                  
                  ,h5( "Based on gene onotlogy (GO) annotation and gene ID mapping of ",
                       a( "315 animal and  plant genomes ",href="https://idepsite.wordpress.com/species/",target="_blank"), 
                       "in Ensembl BioMart release 96 as of 5/20/2019."	 
                       , "In addition, 115 archaeal, 1678 bacterial, and 238 eukaryotic genomes are annotated based on STRING-db v10. 
            Additional pathway data are being collected for some model species from difference sources."
                       ,h5( "Genomes based on STRING-db is marked as STRING-db. If the same genome is included in both Ensembl and STRING-db, users should 
            use Ensembl annotation, as it is more updated and is supported in more functional modules. ")
                       ,includeHTML("human_mouse_source.html")
                       
                       ,br()
                       ,br(),br()
                       
                       #,includeHTML("help.htm")
                       ,h4("Input:")
                       ,"A list of gene ids, separated by tab, space, comma or the newline characters.
		   Ensembl gene IDs are used internally to identify genes. Other types of IDs will be mapped to Ensembl 
		   gene IDs using ID mapping information available in Ensembl BioMart. "
                       ,h4("Output:")
                       ,"Enriched GO terms and pathways:"
                       ,br()
                       ,img(src='enrich.png', align = "center",width="660", height="339")		 
                       ,br(),br()		 
                       ,"In addition to the enrichment table, a set of plots are produced. If KEGG database is choosen, then enriched pathway diagrams are shown, with user's genes highlighted, like this one below:"
                       ,br()
                       ,img(src='KEGG.png', align = "center",width="696", height="494")
                       ,br(),br()
                       ,"Many GO terms are related. Some are even redundant, like \"cell cycle\" and \"cell cycle process\". 
		 To visualize such relatedness in enrichment results, we use a hierarchical clustering tree and network. 
		 In this hierarchical clustering tree, related GO terms are grouped together based on how many genes they share. The size of the solid circle corresponds to the enrichment FDR."
                       ,br()
                       ,img(src='GOtree.png', align = "center",width="700", height="361")
                       ,br(),br()
                       ,"In this network below, each node represents an enriched GO term. Related GO terms are connected by a line, whose thickness reflects percent of overlapping genes. The size of the node corresponds to number of genes." 
                       ,br(),img(src='GOnetwork.png', align = "center",width="500", height="248")
                       ,br(),br()		 
                       ,"Through API access to STRING-db, we also retrieve a protein-protein interaction (PPI) network. In addition to a static network image, users can also get access to interactive graphics at the www.string-db.org web server."
                       ,br(),img(src='PPInetwork.png', align = "center",width="700", height="547")
                       ,br(),br()	
                       ,"ShinyGO also detects transcription factor (TF) binding motifs enriched in the promoters of user's genes."
                       ,br(),br(),img(src='promoter.png', align = "center",width="717", height="288")			 
                       
                       
                       ,br(),h4("Changes:")
                       ,h5("6/6/2021: V0.66 Adjusted interface. ")
                       ,h5("6/2/2021: V0.66 add customized background genes.")
                       ,h5("5/23/2021: V0.65 Database update to Ensembl 103 and STRING-db v11.")
                       ,h5("11/3/2019: V 0.61 Improved visualization based on suggestions from reviewers. Interactive networks.")
                       ,h5("5/20/2019: V0.60 Upgraded to Ensembl Biomart 96. Add annotation from STRING-db v10")
                       ,h5("3/29/2019: V0.51 Update annotation to Ensembl release 95. Interface change. Demo gene lists. Error messages.")
                       ,h5("9/10/2018: V0.5 Upgraded to Ensembl Biomart 92")
                       ,h5("4/30/2018: V0.42 changed figure configurations for tree.")
                       ,h5("4/27/2018: V0.41 Change to ggplot2, add grid and gridExtra packages")		
                       ,h5("4/24/2018: V0.4 Add STRING API, KEGG diagram, tree display and network.")
                       
                  ) ) 	
      ) #tabsetPanel
      ,bsModal("ModalExamplePPI", "Protein-protein interaction(PPIs) networks ", "ModalPPI", size = "large"
               ,h5("By sending your genes to the STRING website, 
			shinyGO is retrieving a sub-network, calculating PPI enrichment, 
		  and generating custom URLs to the STRING website containing your genes. This can take 5 minutes. Patience will pay off! ")
               ,sliderInput("nGenesPPI", label = h5("Genes to include:"), min = 0, max = 400, value = 50,step=10) 
               #,htmlOutput("stringDB_network_link")
               #,tags$head(tags$style("#stringDB_network_link{color: blue; font-size: 15px;}"))
               
               ,plotOutput("stringDB_network1")		 	   
               
      )# bsModal 1
      
      ,bsModal("InteractiveNetwork", "Interactive enrichment networks ", "GONetwork", size = "large"
               ,fluidRow(
                 column(2, actionButton("layoutButtonStatic", "Change layout") ),
                 column(2, downloadButton('enrichmentNetworkPlotDownload', "Download")),
                 column(2, checkboxInput("wrapTextNetworkStatic", "Wrap text", value = FALSE) )  
               )       
               ,plotOutput('enrichmentNetworkPlot')              
               
      )# bsModal 2
      
      ,bsModal("BackgroundGenes", "Customized background genes (Optional but highly recommended)", "backgroundGenes", size = "large"
               ,tags$textarea(id = 'input_text_b', 
                              placeholder = 'Paste all genes from which the gene list is derived. These are all genes whose expression or other activity that you measured. This could be all the genes on a special DNA microarray or all the genes detected by a proteomics experiment. 
                    ', rows = 15, "")  
               ,h4("By default, we compare your gene list with a background of all protein-coding genes in the genome.
     		         When your genes are not selected from genome-wide data, customized background genes might yield more accurate 
		             results for enrichment analysis. For gene lists derived from a typical RNA-seq dataset, 
                     many  use the subset of genes with detectable expression, typically the genes passed a minimum filter. 
                      We can also customize background genes to overcome bias in selection. Currently only less than 30,000 genes are accepted.
		              ")
               
      ),# bsModal 3	
      
      bsModal(id = "geneIDBs", title = "Gene ID Examples",
              trigger = "geneIdButton", size = "large",
              fluidPage(shinyjs::useShinyjs(),
                        sidebarLayout(fluid = TRUE,
                                      sidebarPanel(#Side panel
                                        #see server updateSelectizeInput 
                                        selectizeInput(inputId = "userSpecie",
                                                       label = "What's your specie name?", choices = NULL),
                                        shiny::tags$h5("Can erase and type in box"),
                                        
                                        #see server updateSelectizeInput
                                        selectizeInput(inputId = "userIDtype",
                                                       label = "What's your ID type? (Optional)", choices = NULL),
                                        shiny::tags$h5("Can erase and type in box"),
                                        actionButton(inputId = "submitIDPage", label = "submit"),
                                        actionButton(inputId = "resetIDPage", label = "reset"),
                                        downloadButton(outputId = "downloadIDPage", label = "Download mapping.csv")
                                      ),##End of side panel
                                      mainPanel(reactable::reactableOutput(outputId = "tableResult"),
                                                ##Instructions for the user
                                                shiny::tags$div(
                                                  shiny::tags$h1("Instructions for Usage"),
                                                  shiny::tags$h4("This page's purpose is to give the user some interactive tools to look at our database IDs.
                               There are two different uses for this page, see explanation below:"),
                                                  shiny::tags$ul(
                                                    shiny::tags$li(#Bullet point 1
                                                      shiny::tags$h4("If you only pick a species,
                                       you are receiving table with all the different IDs
                                       related to that species. (Showen below)")
                                                    ),#end of bullet point 1
                                                    shiny::tags$li(#Bullet point 2
                                                      shiny::tags$h4("If you pick a species and an ID type,
                                       a table with all the IDs of the ID type you pick and how they map to ensembl IDs(our preferred ID database)
                                       , and you can download a csv file of the mapping ID.")
                                                    )#end of bullet point 2
                                                  )#end of ul
                                                ), ##End of instructions
                                                reactable::reactableOutput(outputId = "tableDefault")
                                      )##End of main panel
                        )#END of sidebarLayout
              )) #end of gene id ui
      
    ) # mainPanel
  ) #sidebarLayout
  ,tags$head(includeScript("google_analytics.js")) # tracking usage
) #fluidPage




