library(shiny)
shinyUI(
  fluidPage(
    # Application title
    titlePanel("ShinyGO: Gene Ontology Enrichment Analysis & More v0.1"),
      h5("Based on annotation of 69 metazoa and 42 plant genomes in Ensembl BioMart as of 11/15/2016."
            ," Additional  data from",a("MSigDB (human),", href="https://doi.org/10.1093/bioinformatics/btr260") 
         ,a("GSKB (mouse)", href="http://biorxiv.org/content/early/2016/10/24/082511") 
         ,"and",a("  araPath (arabidopsis).", href="https://doi.org/10.1093/bioinformatics/bts421") 
         , "Built with R and", a("Shiny!",href="http://shiny.rstudio.com/")
         ," For feedbacks or data contributions (genes and GO mapping of any species), please"
         ,a("contact us, ",href="mailto:xijin.ge@sdstate.edu?Subject=ShinyGO" )
         , "or visit our",a(" homepage.", href="http://ge-lab.org/")
         ),
     sidebarLayout(
      sidebarPanel(
	  a(h5("Reset all",align = "right"), href="http://ge-lab.org:3838/go/"), 
      tags$style(type="text/css", "textarea {width:100%}"),
      tags$textarea(id = 'input_text', placeholder = 'Just paste gene lists and click Submit. Most types of gene IDs accepted. Double check the guessed species,  adjust options and resubmit if needed. Any feedback is appreciated. ', rows = 8, ""),
      actionButton("goButton", "Submit"),
      h6(" "),
      htmlOutput("selectGO1"),
      selectInput("selectOrg", label = NULL,"Best matching species",width='100%'),  
      numericInput("minFDR", label = h5("P-value cutoff (FDR)"), value = 0.05),
     tags$style(type='text/css', "#minFDR { width:100%;   margin-top:-15px}"),  
     radioButtons("radio", label = NULL, choices = list("Upstream 300bp as promoter" = 300, "Upstream 600bp as promoter" = 600),selected = 300),
     tableOutput('species' ),
      h5( textOutput("text1") )
      ), # sidebarPanel
     mainPanel(
       tabsetPanel(
        tabPanel("Enrichment", tableOutput('table'), 
		    conditionalPanel("input.goButton != 0", 
		                     downloadButton('downloadEnrichment', 'Download') )  
		    )
        ,tabPanel("Details", tableOutput('tableDetail')  )
        ,tabPanel("Genes", tableOutput("conversionTable"), downloadButton('downloadGeneInfo', 'Download')  )
        ,tabPanel("Groups", tableOutput("grouping"), downloadButton('downloadGrouping', 'Download')   )          
        ,tabPanel("Plots", plotOutput("genePlot")  )
        , tabPanel("Genome", plotOutput("genomePlot", width = "100%")  )
        ,tabPanel("Promoter", tableOutput("promoter"), downloadButton('downloadPromoter', 'Download')   )  
       ) #tabsetPanel
     ) # mainPanel
    ), #sidebarLayout
    tags$head(includeScript("google_analytics.js")) # tracking usage
  ) #fluidPage
)



