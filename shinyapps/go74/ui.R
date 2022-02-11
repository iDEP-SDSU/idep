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
columnSelection = list("-log10(FDR)" = "EnrichmentFDR", 
                        "Fold Enrichment" = "FoldEnrichment", 
                        "Genes" =  "nGenes", 
                        "Category Name" = "Pathway")
ui <- fluidPage(
  # Application title
  sidebarLayout(
    
    sidebarPanel(
      titlePanel("ShinyGO v0.741: Gene Ontology Enrichment Analysis + more"),  
       h5("Select or search your species. Or use our best guess."),
       fluidRow( 
         column(9, selectizeInput('selectOrg', 
                     label    = NULL,
                     choices  = " ",
                     multiple = TRUE,
                     options  = list( maxItems     = 1,               
                                      placeholder  = 'Best matching species',
                                      onInitialize = I('function() { this.setValue(""); }')) 
                  )), 
        column(3, actionButton("MorgInfo", "Info"))  
      ),  
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
      

      tableOutput('species' ),
      actionButton("MGeneIDexamples", "Gene IDs examples"),
      h5("Try ", a(" iDEP", href="https://bioinformatics.sdstate.edu/idep/",target="_blank"), "for RNA-Seq data analysis")


    ), # sidebarPanel
    mainPanel(
      tabsetPanel(
        tabPanel("Enrichment" 
                 ,conditionalPanel("input.goButton == 0 "  # welcome screen
                                   ,br()
                                   ,p("Please use the new version  ", a("ShinyGO 0.75.", href="http://bioinformatics.sdstate.edu/go/"))
                                   ,p("Feb. 11, 2022: Like ShinyGO but your genome is not covered?", 
                                   a("Customized ShinyGO", href="http://bioinformatics.sdstate.edu/goc/"), " is now available. 
                                   Its database includes several custom genomes requested by users. To request your genome, fill in this ", 
                                   a("Form.", href="https://forms.gle/zLtLnqxkW187AgT76"), style = "color:red")
                                   ,p("Nov. 15, 2021: Database update. ShinyGO v0.75 available in testing mode. It includes Ensembl database update, new species from Ensembl Fungi and Ensembl Protists, and STRINGdb (5090 species) update to 11.5.")
                                   ,p("Oct. 25, 2021: Interactive genome plot. Identificantion of genomic regions signficantly enriched with user genes.")
                                   ,p("Oct. 23, 2021: Version 0.741 A fully customizable enrichment chart! Switch between bar, dot or lollipop plots.  Detailed gene informations with links on the Genes tab.")
                                   ,p("Oct. 15, 2021: Version 0.74. Database updated to Ensembl Release 104 and STRING v11. We now recommends the use of background genes in enrichment analysis. V.0.74 is much faster with even large set of background genes.")
                                   ,p("We recently hired Jenny for database updates and user support.",
                                       a("Email Jenny ",href="mailto:gelabinfo@gmail.com?Subject=ShinyGO"),
                                       "for questions, suggestions or data contributions.") 
                                   ,p("ShinyGO has not been tested thoroughly! Please let us know if you find any problems.", style = "color:red")
                                   
                                   
                                   #,h4("If your gene IDs are not recognized, please let us know. We might be able to add customized gene mappings to Ensembl gene IDs.")
                                   
                                   ,p("2/3/2020: Now published by", a("Bioinformatics.", href="https://doi.org/10.1093/bioinformatics/btz931",target="_blank"))
                                   ,p("Just paste your gene list to get enriched GO terms and othe pathways for over 520 species, 
				    based on annotation from Ensembl (Plants, Metazoa, Protists, Fungi). An additional 5090 genomes 
				    (including bacteria and fungi) are  annotated based on STRING-db (v.11.5). In addition, it also produces
				    KEGG pathway diagrams with your genes highlighted, hierarchical clustering trees and networks summarizing 
				    overlapping terms/pathways, protein-protein interaction networks, gene characterristics plots, and enriched promoter motifs. 
                 See example outputs below:")			
                                   ,br(),img(src='enrich.png', align = "center",width="660", height="339")
                                   ,br(),img(src='enrichmentChart.png', align = "center",width="700", height="400")
                                   ,br(),br(),img(src='KEGG2.png', align = "center",width="541", height="360")
                                   ,br(),br(),img(src='GOtree3.png', align = "center",width="500", height="258")
                                   ,br(),br(),img(src='GOnetwork2.png', align = "center",width="500", height="248")
                                   ,br(),br(),img(src='PPInetwork2.png', align = "center",width="500", height="391")	
                                   ,br(),br(),img(src='chr.png', align = "center",width="444", height="338")
                                   ,br(),br(),img(src='downSyndrome.png', align = "center",width="371", height="276")				
                                   ,br(),br(),img(src='promoter.png', align = "center",width="717", height="288")					
                 )
                 ,br()
                 ,p("FDR is adjusted from the hypergeometric test. Fold Enrichment is defined as the percentage of genes in your list belonging to a pathway, divided by the corresponding percentage in the background. FDR tells us how likely the enrichment is by chance; Fold Enrichment indicates how drastically genes of a certain pathway is overrepresented. We think the latter deserves at least some attention.")
                 ,div(style="display:inline-block", 
                          selectInput(inputId = "SortPathways",
                              label = NULL,
                              choices = c("Sort by FDR" = "Sort by FDR", 
                                          "Sort by Fold Enrichment" = "Sort by Fold Enrichment", 
                                          "Sort by Genes" =  "Sort by Genes", 
                                          "Sort by Category Name" = "Sort by Category Name"),
                              selected = "Sort by Fold Enrichment" ),
                     style="algn:right")

                 ,tableOutput('EnrichmentTable')	
                 ,conditionalPanel("input.goButton != 0", 
                                   downloadButton('downloadEnrichment', 'Download table with gene IDs')			
                                   ,br(),br()
                 )			

        ) # enrichment tab
        #,tabPanel("Details", tableOutput('tableDetail')  )

 #---Enrichment Chart-----------------------------------------------------------
        ,tabPanel("Chart" 
                  ,plotOutput('enrichChart')
                  ,br(), h4("Adjusting the width of your browser window to resize figure.")
                  ,fluidRow(
                    column(3, selectInput(inputId = "SortPathwaysPlot",
                                          label = h5("Sort Pathway by"),
                                          choices = columnSelection,
                                          selected = columnSelection[2] ) )
                    ,column(3, selectInput(inputId = "SortPathwaysPlotX",
                                           label = h5("x-axis"),
                                           choices = columnSelection[1:3],
                                           selected = columnSelection[2] )  )
                     ,column(3, selectInput(inputId = "SortPathwaysPlotColor",
                                           label = h5("Color"),
                                           choices = columnSelection[1:3],
                                           selected = columnSelection[1] )  )
                     ,column(3, selectInput(inputId = "SortPathwaysPlotSize",
                                           label = h5("Size"),
                                           choices = columnSelection[1:3],
                                           selected = columnSelection[3] )  )
                  ) # first row

                  ,fluidRow(
                    column(3, numericInput(inputId = "SortPathwaysPlotFontSize",
                                           label = h5("Font Size"),
                                           value = 12,
                                           min = 3,
                                           max = 18,
                                           step = 1 ) )
                    ,column(3, numericInput(inputId = "SortPathwaysPlotMarkerSize",
                                           label = h5("Circle Size"),
                                           value = 4,
                                           min = 0,
                                           max = 10,
                                           step = 1 ))
                    ,column(3, selectInput(inputId = "SortPathwaysPlotHighColor",
                                           label = h5("Color:High"),
                                           choices = c("red", "orange", "yellow", "green", "blue", "purple"),
                                           selected = "red"
                                           ))
                    ,column(3, selectInput(inputId = "SortPathwaysPlotLowColor",
                                           label = h5("Color:Low"),
                                           choices = c("red", "orange", "yellow", "green", "blue", "purple"),
                                           selected = "blue"
                                           ))
                  ) # 2nd row
                  ,selectInput(inputId = "enrichChartType",
                                           label = h5("Chart type"),
                                           choices = c("lollipop", "dotplot", "barplot"),
                                           selected = "lollipop"
                                           )
        )

 #---Tree-----------------------------------------------------------
        ,tabPanel("Tree" 
                  ,h5("A hierarchical clustering tree summarizing the correlation among significant pathways 
       listed in the Enrichment tab. Pathways with many shared genes are clustered together. 
       Bigger dots indicate more significant P-values.")
                  ,downloadButton('GOTermsTree4Download','Figure' )
                  ,plotOutput('GOTermsTree')
        )

 #---Enrichment Chart-----------------------------------------------------------        
        ,tabPanel("Network" 
                  ,fluidRow(
                    column(2, actionButton("layoutButton", "Change layout") ),
                    column(2, actionButton("GONetwork", "Static plot") ),
                    column(2, h5("Edge cutoff:"), align="left" ) ,
                    column(2, numericInput("edgeCutoff", label = NULL, value = 0.30, min = 0, max = 1, step = .1), align="right"  ), 
                    column(2, checkboxInput("wrapTextNetwork", "Wrap text", value = TRUE)) 
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
 #---Genes-----------------------------------------------------------        
        ,tabPanel("KEGG"
                 ,h5("Select KEGG pathways in the left to display your genes in pathway diagrams. Selected species only. Takes 1-3 minutes.")
                 ,conditionalPanel("input.selectGO == 'KEGG'", 
                                   htmlOutput('listSigPathways')
                                   ,br(),br(),imageOutput("KeggImage", width = "100%", height = "100%")				
                                   ,br(),br()	 
                 )

      
        )
 #---Genes-----------------------------------------------------------        
        ,tabPanel("Genes"
                  , downloadButton('downloadGeneInfo', 'Download')
                  , tableOutput("conversionTable")
      
        )
 #---Groups-----------------------------------------------------------        
        ,tabPanel("Groups"
                  , downloadButton('downloadGrouping', 'Download')
                  , h5("If Gene Ontology is chosen, your genes are grouped by functional categories defined by high-level GO terms. ")  
                  , tableOutput("grouping")
        
        ) 
 #---Plots-----------------------------------------------------------        
        ,tabPanel("Plots"
                  , h5("The characteristics of your genes are compared with the rest in the genome. Chi-squared and Student's 
              t-tests are run to see if your genes have special characteristics when compared with all the other genes or, if uploaded, a customized background.")
                  , plotOutput("genePlot", inline = TRUE,width='auto',height='auto')
                  , plotOutput("genePlot2", inline = TRUE,width='auto',height='auto')  
        )

 #---Genome-----------------------------------------------------------        
        ,tabPanel("Genome"
                  , h5("The genes are represented by red dots. The purple lines indicate regions where these genes are statistically enriched, compared to the density of genes in the background. We scanned the genome with a sliding window. Each window is further divided into several equal-sized steps for sliding. Within each window we used the hypergeometric test to determine if the presence of your genes are significant. Essentially, the genes in each window define a gene set/pathway, and we carried out enrichment analysis. The chromosomes may be only partly shown as we use the last gene's location to draw the line. Mouse over to see gene symbols. Zoom in regions of interest.")
                  ,plotlyOutput("genomePlotly",height = "900px")
                  ,fluidRow(
                    column(3, selectInput(inputId = "MAwindowSize",
                                           label = h5("Window Size(Mb)"),
                                           selected = 6,
                                           choices = c(1, 2, 4, 6, 8, 10, 15, 20) ))
                    ,column(3, selectInput(inputId = "MAwindowSteps",
                                           label = h5("Steps in a window"),
                                           selected = 2,
                                           choices = c(1, 2, 3, 4)))
                    ,column(3, selectInput(inputId = "chRegionPval", 
                                           label = h5("FDR cutoff for windows"),
                                           selected = 0.00001,
                                           choices = c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))))
                  ,fluidRow(  
                     column(4, checkboxInput("labelGeneSymbol", "Label genes", value = FALSE) )
                     ,column(4, checkboxInput("ignoreNonCoding", "Coding genes only", value = TRUE) )  
                    ,column(4, actionButton("gPlotstatic", "Static plot") ) )                
        )

 #---Genome-----------------------------------------------------------                
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
 #---STRING-----------------------------------------------------------
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

 #---?-----------------------------------------------------------
        ,tabPanel("?"
                  ," For feedbacks, please"
                  ,a("contact us, ",href="mailto:gelabinfo@gmail.com?Subject=ShinyGO" )
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
      
      ,bsModal("BackgroundGenes", "Customized background genes (Highly recommended)", "backgroundGenes", size = "large"
               ,tags$textarea(id = 'input_text_b', 
                              placeholder = 'Paste all genes from which the gene list is derived. These are all genes whose expression or other activity that you measured. This could be all the genes on a special DNA microarray or all the genes detected by a proteomics experiment. 
                    ', rows = 15, "")  
               ,h4("By default, we compare your gene list with a background of all protein-coding genes in the genome.
     		         When your genes are not selected from genome-wide data, customized background genes might yield more accurate 
		             results for enrichment analysis. For gene lists derived from a typical RNA-seq dataset, 
                     many  use the subset of genes with detectable expression, typically the genes passed a minimum filter. 
                      We can also customize background genes to overcome bias in selection. Currently only less than 30,000 genes are accepted.
		              ")
               
      )# bsModal 3	
      
      ,bsModal("geneIDexamples", "What the gene IDs in our database look like?", "MGeneIDexamples", size = "large"
               ,selectizeInput(inputId = "userSpecieIDexample",
                               label = "Select or search for species", choices = NULL)
               ,tableOutput("showGeneIDs4Species")

       )# bsModal 4	

      ,bsModal("orgInfoButton", "Supported species (Search by common and scientific names, or NCBI taxonomy IDs)", "MorgInfo", size = "large"
               ,DT::dataTableOutput('orgInfoTable')

       )# bsModal 5	
      ,bsModal("genomePlotStaticButton", "Static genome plot", "gPlotstatic", size = "large"
                  , h5("Your genes are marked in each of the chromosomes. 
                Note that the scale for each chromosomes are different.")
                  , plotOutput("genomePlot", width = "100%")  

       )# bsModal 6	      
    ) # mainPanel
  ) #sidebarLayout
  ,tags$head(includeScript("google_analytics.js")) # tracking usage
) #fluidPage




