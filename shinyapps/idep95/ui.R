# iDEP user interface, By Steven Ge, Xijin.Ge@sdstate.edu
# Integrated Differential Gene Expression and Pathway analysis
# hosted at http://ge-lab.org/idep/

library(shiny,verbose=FALSE)
library("shinyAce",verbose=FALSE) # for showing text files, code
library(shinyBS,verbose=FALSE) # for popup figures
library(plotly,verbose=FALSE)
library('shinyjs', verbose = FALSE)
library('shinyjs', verbose = FALSE)
library('reactable', verbose = FALSE)
library(visNetwork) # interative network graphs
iDEPversion = "iDEP.951"

shinyUI(
navbarPage(
iDEPversion,
  id='navBar', 
  
#================================================================================================== 
#   Load Data
#================================================================================================== 
  
  tabPanel("Load Data", value = 1,
 # titlePanel(h5("Upload Files")),
  sidebarLayout(
    # sidebar---------------------------------
    sidebarPanel(
  
      actionButton("goButton", "Click here to load demo data"),
      tags$head(tags$style("#goButton{color: red;
                                 font-size: 16px;
                                 font-style: italic;
                                 }"),
        # Fix excessive padding around the body 
        tags$style("
          body {
            padding: 0 !important;
          }"
                ))                    
      ,h5(" and just click the tabs for some magic!", style = "color:red")
      ,p(HTML("<div align=\"right\"> <A HREF=\"javascript:history.go(0)\">Reset</A></div>" ))
      ,strong("1. Recommended: Select or search for your species. ")
      ,fluidRow( 
         column(9, selectizeInput('selectOrg', 
                     label    = NULL,
                     choices  = " ",
                     multiple = TRUE,
                     options  = list( maxItems     = 1,               
                                      placeholder  = 'Best matching species',
                                      onInitialize = I('function() { this.setValue(""); }')) 
                  )), 
        column(3, actionButton("MorgInfo", "Info"))  
      )  
      ,conditionalPanel("input.selectOrg == 'NEW'",
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
      ) # conditionalPanel


      ,radioButtons("dataFileFormat", 
                     label = "2. Choose data type", 
                     choices = list("Read counts data (recommended)"                                          = 1, 
                                     "Normalized expression values (RNA-seq FPKM, microarray, etc.)"          = 2,
                                     "Fold-changes and corrected P values from CuffDiff or any other program" = 3),
                     selected = 1
      )      
      ,conditionalPanel("input.dataFileFormat == 3",
        checkboxInput("noFDR", "Fold-changes only, no corrected P values", value = FALSE)
      )
      
      ,fileInput('file1', '3. Upload expression data (CSV or text)',
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    'text/tab-separated-values',
                    'text/plain',
                    '.csv',
                    '.tsv'          
                  ) 
      )
      ,a(h4("Analyze public RNA-seq datasets for 9 species"), href="http://bioinformatics.sdstate.edu/reads/")
      ,fileInput('file2', h5('Optional: Upload an experiment design file(CSV or text)'),
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    'text/tab-separated-values',
                    'text/plain',
                    '.csv',
                    '.tsv'          
                  )
      )
      ,tableOutput('species' )

      ,actionButton("MGeneIDexamples", "Example gene IDs")
      ,bsModal("geneIDexamples", "What the gene IDs in our database look like?", "MGeneIDexamples", size = "large"
               ,selectizeInput(inputId = "userSpecieIDexample",
                               label = "Select or search for species", choices = NULL)
               ,tableOutput("showGeneIDs4Species")

       )# bsModal 4	

      ,bsModal("orgInfoButton", "Search annotated species by common and scientific names, or NCBI taxonomy id. For other species, you can still use iDEP, except enrichment analysis.", "MorgInfo", size = "large"
               ,DT::dataTableOutput('orgInfoTable')

       )# bsModal 4	
      ,h5("Try ", a(" ShinyGO", href="https://bioinformatics.sdstate.edu/go/",target="_blank"), "for GO enrichment analysis")
      ,a( h5("?",align = "right"), href="https://idepsite.wordpress.com/data-format/",target="_blank")
                                                                                       # new window
    ), #sidebarPanel
  
  # Main Panel -------------------------------------
    mainPanel(
      shinyjs::useShinyjs(),
      tableOutput('sampleInfoTable')
      ,tableOutput('contents')

      ,div(id='loadMessage',
           h4('Loading R packages, please wait ... ... ...'))
      ,htmlOutput('fileFormat')
                                   ,br(),h4("Trusted charities providing Aid in Ukraine selected by ", a("CharityWatch,",
                                   href="https://www.charitywatch.org/charity-donating-articles/top-rated-charities-providing-aid-in-ukraine"),
                                   "charities and individuals verified by",
                                   a(" GoFundMe.", href=("https://www.gofundme.com/en-ie/c/act/donate-to-ukraine-relief?utm_source=email&utm_medium=marketing&utm_content=annoucement&utm_campaign=022522_helpukraine_send14_dedicatedpage"))
                                   ), br()
      ,p("April 25, 2022: Gene ID conversion is much faster now, even when species has to be guessed. So is the DEG2 tab.")
      ,p("April 24, 2022: Add a tab for visualizing the fold-change of all genes in all KEGG diagrams across all comparisons!")

      ,p("Feb. 11, 2022: Like iDEP but your genome is not covered?", 
      a("Customized iDEP", href="http://bioinformatics.sdstate.edu/idepc/"), " is now available. 
      Its database includes several custom genomes requested by users. To request to add new species/genome, fill in this ", 
      a("Form.", href="https://forms.gle/zLtLnqxkW187AgT76"), style = "color:red")
  ,p( a("Email Jenny for questions.",href="mailto:gelabinfo@gmail.com?Subject=iDEP"), "Dr. Ge is notorisly slow in responding to emails.") 
      ,p("If it is slow, restart from a new browser window (not a new tab). You will be assigned to a new worker computer.")
      ,p("iDEP has not been thoroughly tested. Please let us know if you find any issue/bug.", style = "color:red")
      ,h3("iDEP: Integrated Differential Expression and Pathway analysis")
      ,br(),img(src='flowchart.png', align = "center",width="562", height="383")  
      ,br(),img(src='figs.gif', align = "center",width="640", height="480")

    
       #,img(src='flowchart.png', align = "center",width="562", height="383")
     # ) # conditionalPanel

    ) # main panel
  ) #sidebarLayout
) #tabPanel
       
     
#================================================================================================== 
#   Pre-Process
#================================================================================================== 
  ,tabPanel("Pre-Process", value = 2,
    sidebarLayout(
  # sidebar of pre-process -----------------------------------
      sidebarPanel(
    
      # If FPKM data
        conditionalPanel("input.dataFileFormat == 2",   
          strong("Only keep genes above this level in at least n samples:" )
          ,fluidRow(
            column(6, numericInput("lowFilter", label = h5(" Min. level"), value = -1000))
            ,column(6, numericInput("NminSamples2", label = h5("n samples"), value = 1) )          
          )
          ,tags$style(type='text/css', "#lowFilter { width:100%;   margin-top:-12px}")
          ,tags$style(type='text/css', "#NminSamples2 { width:100%;   margin-top:-12px}")      
          ,radioButtons("transform", "Log Transformation",c("No"=FALSE,"Yes"=TRUE) )
          ,numericInput("logStart", label = h5("Constant c for started log: log(x+c)"), value = 1)
          ,tags$style(type='text/css', "#logStart { width:100%;   margin-top:-12px}")
          ,textOutput("textTransform") 
          ,tags$head( tags$style("#textTransform{color: blue;
                      font-size: 16px;
                      font-style: italic;}"
                     )
          )

    ) # conditionalPanel

   
    # If read count data
        ,conditionalPanel("input.dataFileFormat == 1",            
          strong("Keep genes with minimal counts per million (CPM) in at least n libraries:")
          ,fluidRow(
            column(6, numericInput("minCounts", label = h5("Min. CPM"), value = 0.5) )
            ,column(6, numericInput("NminSamples", label = h5("n libraries"), value = 1) )
          ) # fluidRow
          ,tags$style(type='text/css', "#minCounts { width:100%;   margin-top:-12px}")
          ,tags$style(type='text/css', "#NminSamples { width:100%;   margin-top:-12px}")
              
          ,radioButtons("CountsTransform", "Transform counts data for clustering & PCA.",  
                        c("VST: variance stabilizing transform" = 2, 
                          "rlog: regularized log (slow) "       = 3,
              "EdgeR: log2(CPM+c)"                  = 1
              ),                          
                        selected = 1 
          )
      
          ,conditionalPanel("input.CountsTransform == 1",
            fluidRow(
              column(5, h5("Pseudo count c:")  )
              ,column(7, numericInput("countsLogStart", label = NULL, value = 4) )
            )        
          )        
        ) # conditionalPanel
    
        ,selectInput("missingValue", 
                     label   = "Missing values imputation:",
                     choices = list("Gene median"                  = "geneMedian",
                                     "Treat as zero"               = "treatAsZero", 
                                     "Median within sample groups" = "geneMedianInGroup"),
                     selected = "geneMedian")
        ,actionButton("genePlot1", "Plot one or more genes")
        ,br(),br()
        ,actionButton("examineDataB", "Search processed data")
        ,br(),br()
        ,checkboxInput("noIDConversion", "Do not convert gene IDs to Ensembl.", value = FALSE)
        ,downloadButton('downloadProcessedData', 'Processed data') 
        ,conditionalPanel("input.dataFileFormat == 1", 
           downloadButton('downloadConvertedCounts', 'Converted counts data') )
        ,downloadButton('downloadEDAplot', 'High-resolution figure')  
        ,br(),br()
        ,textOutput('nGenesFilter')
        ,tags$head(tags$style("#nGenesFilter{color: blue;
                       font-size: 16px;
                       font-style: italic;
                       }"))                 
        ,textOutput("readCountsBias") 
        ,tags$head(tags$style("#readCountsBias{color: red;
                 font-size: 16px;
                 font-style: italic;
                 }" ) )
        ,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/pre-process/",target="_blank")
    ), # sidebarPanel

    #main panel of pre-process tab ------------------------------------------------------------------------
    
      mainPanel(
        h5("Aspect ratios of figures can be adjusted by changing the width of browser window.")
        ,conditionalPanel("input.dataFileFormat == 1", plotOutput("totalCounts") )
        ,fluidRow( 
          column(4, selectInput("scatterX", "Select a sample for x-axis", choices = 1:5, selected = 1))  
          ,column(4, selectInput("scatterY", "Select a sample for y-axis", choices = 1:5, selected = 2) )
          ) 
        ,plotOutput("EDA")
        ,bsModal("modalExample10", "Converted data (Most variable genes on top)", 
                 "examineDataB", size = "large", DT::dataTableOutput('examineData'))
         
        ,bsModal("modalExample1021", "Search for genes", "genePlot1", size = "large", 
          textInput("geneSearch", "Enter full or partial gene ID, or list of genes separated by semicolon:", "SNCA;Robo3;GAPDH"),
          checkboxInput("genePlotBox", label = "Show individual samples", value = FALSE),
          plotOutput("genePlot"),
          conditionalPanel("input.genePlotBox == 0", 
          checkboxInput("useSD", label = "Use standard deviation instead of standard error", value = FALSE))
          ,downloadButton('downloadGenePlot', 'Figure')  
        )
    
      ) # mainPanel
    )  #sidebarLayout     
  ) #tabPanel
  
  
#================================================================================================== 
#   Heat map
#================================================================================================== 
  ,tabPanel("Heatmap", value = 3,
    sidebarLayout(

      # sidebar of heatmap -----------------------------------
      sidebarPanel(
        sliderInput("nGenes", 
                    label = h4("Most variable genes to include:"), 
                    min   = 0, 
                    max   = 12000, 
                    value = 1000,
                    step  = 100) 
        ,actionButton("showGeneSD_heatmap", "Gene SD distribution")    
        ,actionButton("showStaticHeatmap", "Interactive heatmap")
        ,br()
        ,actionButton("showCorrelation", "Correlation matrix")  
        ,actionButton("showSampleTree", "Sample Tree")      
        ,HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />') # a solid line
        ,strong("Customize hierarchical clustering (Default values work well):")
        ,fluidRow(
          column(3, h5("Color")  )
          ,column(9, selectInput("heatColors1", label = NULL,"green-black-red",width='100%') )
        )
        ,fluidRow(
          column(4, h5("Distance")  )
          ,column(8, selectInput("distFunctions", label = NULL,"Correlation",width='100%') )
        )          
        ,fluidRow(
           column(4, h5("Linkage")  )
           ,column(8, selectInput("hclustFunctions", label = NULL,"average",width='100%') )
        )
        ,fluidRow(
           column(8, h5("Cut-off Z score")  )
           ,column(4, numericInput("heatmapCutoff", label = NULL, value = 4,min=2,step=1) )
        )
        ,checkboxInput("geneCentering", "Center genes (substract mean)", value = TRUE)          
        ,checkboxInput("geneNormalize", "Normalize genes (divide by SD)", value = FALSE)
        ,checkboxInput("sampleCentering", "Center samples (substract mean)", value = FALSE)  
        ,checkboxInput("sampleNormalize", "Normalize samples(divide by SD)", value = FALSE)
        ,checkboxInput("noSampleClustering", "Do not re-order or cluster samples", value = FALSE)
        ,htmlOutput('listFactorsHeatmap')
        ,downloadButton('downloadData', 'Heatmap data')
        ,downloadButton('downloadHeatmap1', 'High-resolution figure')       
        ,br(),a(h5("?",align = "right"), href="https://idepsite.wordpress.com/heatmap/",target="_blank")
      ), # sidebarPanel
    
      # Main Panel of Heat map tab ---------------------------------------------------------------------
      mainPanel(          
        plotOutput("heatmap1")      
        ,bsModal("modalExample8", "Correlation matrix using top 75% genes", "showCorrelation", size = "large"
          ,downloadButton("downloadCorrelationMatrix","Data")
          ,downloadButton('downloadCorrelationMatrixPlot', 'Figure') 
          ,checkboxInput("labelPCC", "Label w/ Pearson's correlation coefficients", value = TRUE)
          ,plotOutput("correlationMatrix")

        )  
    
        ,bsModal("modalExample228", "Hierarchical clustering tree.", "showSampleTree", size = "large"
          ,h4("Using genes with maximum expression level at the top 75%. Data is transformed 
           and clustered as specified in the main page. ")
          ,plotOutput("sampleTree") 
          ,downloadButton('downloadSampleTree', 'Figure') 
        )
    
        ,bsModal("modalExample28", "Heatmap with hierarchical clustering tree", "showStaticHeatmap", 
           size = "large",
           sliderInput("nGenesPlotly", 
                       label = h4("Most variable genes to include:"),
                       min   = 0, 
                       max   = 12000, 
                       value = 50,
                       step  = 100),
           h5("Mouse over to see gene names. To zoom, click and drag up or downward and release."),
           plotlyOutput("heatmapPlotly", width = "100%", height = "800px")
        )
       
        ,bsModal("modalExample1233", "Distribution of variations among genes", "showGeneSD_heatmap", 
           size = "large"
        ,downloadButton("downloadDistributionSD1","Figure")
        ,plotOutput('distributionSD_heatmap'))
            
      ) # mainPanel
    )  #sidebarLayout     
  ) #tabPanel

  
  
#================================================================================================== 
#  k-Means
#================================================================================================== 
  ,tabPanel("k-Means",value = 4,
    sidebarLayout(
  
  # sidebar of k-Means -----------------------------------
      sidebarPanel(
        sliderInput("nGenesKNN", 
                    label = h4("Most variable genes to include "),
                    min   = 0, 
                    max   = 12000, 
                    value = 2000,
                    step  = 100) 
          
        ,sliderInput("nClusters", 
                     label = h4("Number of Clusters"), 
                     min   = 2, 
                     max   = 20, 
                     value = 4,
                     step  = 1) 
        ,actionButton("KmeansReRun", "Re-Run")
        ,actionButton("NClusters", "How many clusters?")
        ,actionButton("showGeneSD", "Gene SD distribution")
        ,actionButton("geneTSNE", "t-SNE map")
        ,selectInput("kmeansNormalization", h5("Normalize by gene:"), 
                     choices = list("Mean center"       = "geneMean",
                                     "Standardization"  = "geneStandardization",
                                      "L1 Norm"         = "L1Norm"),
                     selected = "Standardization")        
        ,tags$style(type='text/css', "#kmeansNormalization { width:100%;   margin-top:-9px}")
        ,actionButton("showMotifKmeans", "Enriched TF binding motifs")
        ,br(),downloadButton('downloadDataKmeans', 'K-means data')
        ,downloadButton('downloadKmeansHeatmap', 'High-resolution figure')
    
        # a solid line as divider
        ,HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />') 
        ,h5("Pathway database")
        ,htmlOutput("selectGO3")
        ,tags$style(type='text/css', "#selectGO3 { width:100%;   margin-top:-9px}")
        ,checkboxInput("removeRedudantSets", "Remove redudant genesets", value = TRUE)
        ,checkboxInput("useFilteredAsBackground", "Filtered genes as background for enrichment", value = TRUE)
        ,actionButton("ModalEnrichmentPlotKmeans", "Visualize enrichment")
        ,downloadButton('downloadKmeansGO',"Enrichment details")             
        ,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/k-means/",target="_blank")
      ), # sidebar ends
    
    # main panel of k-Means -------------------------------------------------------------------------
      mainPanel(
        plotOutput("KmeansHeatmap",inline=TRUE)
        ,br()
        ,h4("Enriched pathways for each cluster")
        ,tableOutput("KmeansGO")
    
        ,bsModal("modalExample2", "Enriched TF binding motifs in promoters of Kmeans clusters", 
             "showMotifKmeans", size = "large"
          ,radioButtons("radioPromoterKmeans", 
                        label    = NULL, 
                        choices  = list("Upstream 300bp as promoter" = 300, 
                                        "Upstream 600bp as promoter" = 600),
                        selected = 300)
          ,tableOutput("KmeansPromoter")
        )
       
        ,bsModal("modalExample9", "Determining the number of clusters (k)", "NClusters", size = "large",
          h5("Following the elbow method, one should choose k so that adding another cluster 
          does not substantially reduce the within groups sum of squares."
          ,a("Wikipedia", href="https://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set",
          target="_blank"))
          ,plotOutput("KmeansNclusters")
        )
      
        ,bsModal("modalExample229", "t-SNE plot of genes", "geneTSNE", size = "large",
          h5("We use the dimension reduction algorith "
            ,a("t-SNE", href="https://lvdmaaten.github.io/tsne/",target="_blank"), 
            "to map the top genes. Examine the distribution can help choose the nubmer of clusters in k-Means. ")          
          ,checkboxInput("colorGenes", "Color genes by the results of k-Means", value = TRUE)
          ,actionButton("seedTSNE", "Re-calculate using different random numbers")
          ,plotOutput("tSNEgenePlot")
        )
      
        ,bsModal("modalExample233", "Distribution of variations among genes", "showGeneSD",
                 size = "large"
                 ,downloadButton('downloadDistributionSD',"Figure") 
                 ,plotOutput('distributionSD'))
    
        ,bsModal("ModalEnrichmentPlotKmeans1", "Visualize enrichment", "ModalEnrichmentPlotKmeans",
                 size="large"
          ,h5("Gene sets closer on the tree share more genes. Sizes of dot correspond to adjuested Pvalues")
          ,downloadButton('enrichmentPlotKmeans4Download',"Figure")
      ,plotOutput('enrichmentPlotKmeans')
        )
    
      ) # mainPanel
    )  #sidebarLayout     
  ) #tabPanel

  
#================================================================================================== 
#  PCA
#================================================================================================== 
  ,tabPanel("PCA",value = 5,
    sidebarLayout(
  
  # sidebar of PCA ----------------------------------------------------------------------------------
      sidebarPanel(
        radioButtons("PCA_MDS", "Methods", 
                     c("Principal Component Analysis"      = 1,            
                        "Multidimensional Scaling"         = 3, 
                        "t-SNE"                            = 4,
                        "Pathway Analysis of PCA rotation" = 2             
            ))
            
        ,conditionalPanel("input.PCA_MDS == 2" # only show if PCA or MDS (not pathway)
          ,htmlOutput('selectGO6')
        )
      
        ,conditionalPanel("input.PCA_MDS != 2" # only show if PCA or MDS (not pathway)
          ,htmlOutput('listFactors')
          ,htmlOutput('listFactors2')
        )
      
        ,br()
        ,downloadButton('downloadPCAData', 'Coordinates')
        ,downloadButton('downloadPCA', 'High-resolution figure')
        ,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/pca/",target="_blank")        
      ), # end of sidebar
    
      # main Panel of PCA ------------------------------------------------------------------------------
      mainPanel(
        plotOutput("PCA", inline=TRUE)
        ,conditionalPanel("input.PCA_MDS == 1", # only show if t-SNE
          fluidRow( 
            column(4, selectInput("PCAx", "Principal component for x-axis", choices = 1:5, selected = 1))  
            ,column(4, selectInput("PCAy", "Principal component for y-axis", choices = 1:5, selected = 2) )
          ) )
#        ,tags$style(type='text/css', "#PCAx { width:100%;   margin-top:-12px}") 
#        ,tags$style(type='text/css', "#FirstPCA { width:100%;   margin-top:-12px}")    
        ,br()
        ,conditionalPanel("input.PCA_MDS == 4", # only show if t-SNE
          actionButton("tsneSeed2", "Re-calculate t-SNE"),br(),br() )
        
        ,conditionalPanel("input.PCA_MDS == 1 |input.PCA_MDS == 2 " # only show if PCA or MDS (not pathway)
          ,htmlOutput('PCA2factor')
          ,br(),br() )
         
      
      ) # mainPanel
    )  #sidebarLayout     
  ) #tabPanel

#================================================================================================== 
#  DEG1: Differentially expressed gene 1
#================================================================================================== 
  ,tabPanel("DEG1",value = 6,
    sidebarLayout(
  
    # sidebar of DEG1 --------------------------------------------------------------------------------
      sidebarPanel(
        h5("Identifying Differential Expressed Genes (DEGs). See next tab for details."),
        conditionalPanel("input.dataFileFormat == 1",
          selectInput("CountsDEGMethod", "Method:", 
                      choices = list("DESeq2"      = 3,
                                     "limma-voom"  = 2,
                                     "limma-trend" = 1),
                      selected = 3)        
          ,tags$style(type='text/css', "#CountsDEGMethod { width:100%;   margin-top:-12px}")
        )
    
        ,conditionalPanel( "input.dataFileFormat == 2", h5("Using the limma package")        )        
        ,fluidRow(
          column(5,numericInput("limmaPval",  # FDR cutoff
                 label = h5("FDR cutoff"), 
                 value = 0.1,
                 min   = 1e-5,
                 max   = 1,
                 step  =.05)  )
          ,column(7, numericInput("limmaFC",   # fold change cutoff
                  label = h5("Min fold change"),
                  value = 2,
                  min   = 1,
                  max   = 100,
                  step  = 0.5) )
        ) # fluidRow
        ,tags$style(type='text/css', "#limmaPval { width:100%;   margin-top:-12px}")
        ,tags$style(type='text/css', "#limmaFC { width:100%;   margin-top:-12px}")
        ,actionButton("modelAndComparisons", "Select factors & comparisons")
        ,tags$head(tags$style("#modelAndComparisons{color: blue;font-size: 15px;}"))        
        ,br(),br() 
        ,fluidRow(
           column(6,actionButton("showVenn", "Venn Diagram") )
           #, column(6,actionButton("showDEGstats", "Barplot"))
        )
        ,br(),fluidRow(          
        column(8, downloadButton('downloadGeneListsGMT', 'Gene lists') ) )
        ,br()
# the fold change data is wrong when using LIMMA from DEG1 tab
        ,fluidRow(
           column(8, downloadButton('download.DEG.data', 'FDR & fold-changes for all genes(Only use with DESeq2) ') )
        )
        ,downloadButton('downloadSigGeneStats', 'Figure')
        ,br(),h4( textOutput("textLimma") )
        ,tags$head(tags$style("#textLimma{color: blue;font-size: 15px;}"))  
        ,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/degs/",target="_blank")
    
        , width = 4
      ), # end of sidebar
  
            
      mainPanel(
        plotOutput('sigGeneStats')
        ,br(),br(),
        h4("Numbers of differentially expressed genes for all comparisons.
           \"B-A\" means B vs. A. Interaction terms start with \"I:\" ")
        ,tableOutput('sigGeneStatsTable')
    
        ,bsModal("modalExample", "Venn Diagram", "showVenn", size = "large",
           checkboxInput("UpDownRegulated", label = "Split gene lists by up- or down-regulation", value = FALSE)
           ,htmlOutput('listComparisonsVenn')
           ,downloadButton("DownloadVenn","Figure")
           ,plotOutput("vennPlot")           
           )
      
        ,bsModal("modalExample21", "Build model and/or select comparisons", "modelAndComparisons", size = "large",
          fluidRow(
            column(6, htmlOutput('listFactorsDE'))
            ,column(6, htmlOutput('listBlockFactorsDE') ) 
          )
      
          ,fluidRow(
           column(6, htmlOutput('selectReferenceLevels1'))
           ,column(6, htmlOutput('selectReferenceLevels2') ) 
          )
          ,fluidRow(
            column(6, htmlOutput('selectReferenceLevels3'))
            ,column(6, htmlOutput('selectReferenceLevels4') ) 
          )
          ,fluidRow(
            column(6, htmlOutput('selectReferenceLevels5'))
            ,column(6, htmlOutput('selectReferenceLevels6') ) 
          )            
          ,htmlOutput('listInteractionTerms')
          ,textOutput('experimentDesign')
          ,tags$head(tags$style("#experimentDesign{color: red;font-size: 16px;}"))    
          ,htmlOutput('listModelComparisons')
          ,actionButton("submitModelButton", "Submit & re-calculate",style="float:center")
          ,tags$head(tags$style("#submitModelButton{color: blue;font-size: 20px;}"))
          ,h5("Close this window to see results.")
          ,a(h5("More info on DESeq2 experiment design",align = "right"),
            href="http://rpubs.com/ge600/deseq2",target="_blank")
        ) # end of model construction window
      
      ) # mainPanel
    )  #sidebarLayout     
  ) #tabPanel

  
#================================================================================================== 
#  DEG2: Differentially expressed genes 2
#================================================================================================== 
  ,tabPanel("DEG2",value = 7,
    sidebarLayout(
  
      # sidebar of DEG2 --------------------------------------------------------------------------------
      sidebarPanel(
        h5("Examine the results of DEGs for each comparison")
        ,htmlOutput("listComparisons")
        ,br()
        ,actionButton("showVolcano", "Volcano Plot")
        ,actionButton("showMAplot", "MA Plot")  
        ,actionButton("showScatter", "Scatter Plot") 
        ,actionButton("showMotif", "TF binding motifs in promoters")
        # ,tags$style(type='text/css', "#showMotif { width:100%;   margin-top:-12px}")
        ,downloadButton('download.selectedHeatmap.data', "Gene list & data" )
        ,tags$style(type='text/css', "#download.selectedHeatmap.data { width:100%;   margin-top:-12px}")         
        ,downloadButton('downloadSelectedHeatmap',"High-resolution figure")

        ,HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />') 
        ,h5("Enrichment analysis for DEGs:")
        ,htmlOutput("selectGO2")
        ,checkboxInput("UseFilteredGenesEnrich", 
                         label = "Filtered genes as background for enrichment", 
                         value = TRUE)
        ,tags$style(type='text/css', "#selectGO2 { width:100%;   margin-top:-9px}")
        ,actionButton("ModalEnrichmentPlot", "Enrichment tree")
     
        ,actionButton("ModalVisNetworkDEG", "Network (New!)" )  
        ,tags$head(tags$style("#ModalVisNetworkDEG{color: red}"))   
        ,downloadButton('downloadGOTerms', "Enrichment details" )
        ,actionButton("STRINGdb_GO", "Enrichment using STRING API")
        ,h5("Also try",  a("ShinyGO", href="http://ge-lab.org/go/", target="_blank") )  
        ,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/degs/", target="_blank")        
        
        ,width=4
      ), # end of sidebar
      
      # main Panel of KEG2 -----------------------------------------------------------------------------    
      mainPanel(
        #h4("Expression pattern of DEGs for selected comparison:")
        plotOutput("selectedHeatmap", inline=TRUE)
        ,br()
        ,h4("Enriched pathways in DEGs for the selected comparison:")            
        ,tableOutput("geneListGO")
        ,h4("Top Genes for selected comparison:"),tableOutput('geneList')
    
      # Volcano Plot
        ,bsModal("modalExample4", "Volcano plot", "showVolcano", size = "large",
          checkboxInput("volcanoPlotBox", 
                         label = "Show interactive version w/ gene symbols", 
                         value = FALSE)
          ,conditionalPanel("input.volcanoPlotBox == 0",  
            downloadButton('downloadVolcanoPlot',"Figure" )
            ,plotOutput("volcanoPlot"))          
          ,conditionalPanel("input.volcanoPlotBox == 1",
            plotlyOutput("volcanoPlotly", width = "550px", height = "550px"))           
        ) #bsModal
    
      # Scatter plot
        ,bsModal("modalExample5", "Scatter plot", "showScatter", size = "large",
          checkboxInput("scatterPlotBox", 
                        label = "Show interactive version w/ gene symbols", 
                        value = FALSE)
          ,conditionalPanel("input.scatterPlotBox == 0",  
            downloadButton('downloadScatterPlot',"Figure" ),  
            plotOutput("scatterPlot") )
          ,conditionalPanel("input.scatterPlotBox == 1",
            plotlyOutput("scatterPlotly", width = "550px", height = "550px") )
        ) #bsModal
    
      #Enrichment plot: Tree
        ,bsModal("ModalEnrichmentPlot1", "Visualize enrichment", "ModalEnrichmentPlot", size="large"
          ,h5("Gene sets closer on the tree share more genes. Sizes of dot correspond to adjuested Pvalues")
          ,downloadButton('enrichmentPlotDEG24Download',"Figure")
          ,plotOutput('enrichmentPlotDEG2') 
        ) #bsModal
      
      # visNetwork ------------------------------
        ,bsModal("ModalVisNetworkDEG1", "Network of pathways", "ModalVisNetworkDEG", size="large"
          ,h5("Connected gene sets share more genes. Color of node correspond to adjuested Pvalues.")            
          ,fluidRow(
            column(2, actionButton("layoutVisDEG", "Change layout") ),
            column(1, h5("Cutoff:"), align="right" ) ,
            column(2, numericInput("edgeCutoffDEG", label = NULL, value = 0.30, min = 0, max = 1, step = .1), align="left"  ), 
            column(2, checkboxInput("wrapTextNetworkDEG", "Wrap text", value = TRUE)), 
            column(1, downloadButton("visNetworkDEGDownload","Network") ),
            column(1, downloadButton("downloadEdgesDEG", "Edges")) ,
            column(1, downloadButton("downloadNodesDEG", "Nodes"))
           ) 
 
		   ,selectInput("upORdownRegDEG", NULL,
			 c("Both Up & Down" = "Both",
			   "Up regulated" = "Up",
			   "Down regulated" = "Down")) 

           ,h6("Two pathways (nodes) are connected if they share 30% (default, adjustable) or more genes. 
			   Green and red represents down- and up-regulated pathways. You can move the nodes by dragging them, zoom in and out by scrolling, 
			   and shift the entire network by click on an empty point and drag. 
			   Darker nodes are more significantly enriched gene sets. 
			   Bigger nodes represent larger gene sets.  
			   Thicker edges represent more overlapped genes. ")    
           ,visNetworkOutput("visNetworkDEG",height = "800px", width = "800px")
        ) #bsModal
       # end VisNetwork --------------------------------    
           
        # M-A plot
        ,bsModal("modalExample5555", "M-A plot", "showMAplot", size = "large",
           checkboxInput("MAPlotBox", label = "Show interactive version w/ gene symbols", value = FALSE)
           ,conditionalPanel("input.MAPlotBox == 0",  
             downloadButton('downloadMAPlot',"Figure" ),
             plotOutput("MAplot") )
           ,conditionalPanel("input.MAPlotBox == 1",
             plotlyOutput("MAplotly",width = "550px", height = "550px") )
        ) #bsModal
       
      # STRING API
        ,bsModal("modalExampleSTRING", "Enrichment and network visualization using STRING API", 
                 "STRINGdb_GO", size = "large"
                 ,h5("iDEP tries to match your species with the archaeal, 
                     bacterial, and eukaryotic species in the",
                     a(" STRING server", 
                       href="https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html",target="_blank"),
                     " and send the DEGs. If it is running, please wait until it finishes. This can take 5 minutes", 
                     " especially for the first time when iDEP downloads large annotation files.")
                  ,htmlOutput("STRINGDB_species_stat") 
                  ,tags$head(tags$style("#STRINGDB_species_stat{color: blue;font-size: 15px;}"))            
                  ,selectizeInput('speciesName', 
                                  label    = NULL,
                                  choices  = " ",
                                  multiple = TRUE,
                                  options  = list( maxItems     = 1,               
                                                   placeholder  = 'Species name (e.g. Homo sapiens)',
                                                   onInitialize = I('function() { this.setValue(""); }'))                 
                  )
                  ,textOutput('STRINGDB_mapping_stat')
                  ,tags$head(tags$style("#STRINGDB_mapping_stat{color: blue;font-size: 15px;}"))
                  ,br()
                  ,actionButton("ModalPPI","PPI network of DEGs"),br(),br()
                  ,selectInput("STRINGdbGO", 
                               label   = "Functional Enrichment", 
                               choices = list("GO Biological Process"  = "Process"
                                              ,"GO Cellular Component" = "Component"
                                              ,"GO Molecular Function" = "Function"
                                              ,"KEGG"                  = "KEGG"
                                              ,"Pfam"                  = "Pfam"
                                              ,"InterPro"              = "InterPro")
                                ,selected = "Process"),
                 numericInput(inputId = "STRINGFDR",
                              label = "FDR cutoff",
                              value = 0.01,
                              step = 0.01),
                 actionButton("submit2STRINGdb", "Submit"),
                 downloadButton("STRING_enrichmentDownload")
                  ,tableOutput("stringDB_GO_enrichment")
        ) #bsModal
    
        ,bsModal("ModalExamplePPI", "Analyze top DEGs on protein interaction networks (PPIs)", 
                 "ModalPPI", size = "large"
                 ,h5("By sending top DEGs (ranked by fold-change) to the STRING website, 
                      iDEP is retrieving a sub-network, calculating PPI enrichment, 
                      and generating custom URLs to the STRING website containing your genes. 
                      This can take 5 minutes. Patience will pay off! ")
                 ,sliderInput("nGenesPPI", 
                              label = h5("Top up- or down-regulated genes to include:"), 
                              min   = 0, 
                              max   = 400, 
                              value = 100,
                              step  = 10) 
                 #,htmlOutput("stringDB_network_link")  #no longer working due to update in STRING-db
                 #,tags$head(tags$style("#stringDB_network_link{color: blue; font-size: 15px;}"))
                 ,br(),br()
                 ,h5("Interactions among proteins encoded by top up-regulated proteins", align = "center")
                 ,plotOutput("stringDB_network1")    
        ) #bsModal

        ,bsModal("modalExample1", "Enriched TF binding motifs in promoters of DEGs", "showMotif", size = "large"
                  ,radioButtons("radio.promoter", 
                                label = NULL, 
                                choices = list("Upstream 300bp as promoter" = 300, 
                                               "Upstream 600bp as promoter" = 600),
                                selected = 300)
                 ,tableOutput("DEG.Promoter")
        ) #bsModal
           
      ) # mainPanel
    )  #sidebarLayout     
  ) #tabPanel  

#================================================================================================== 
#  Pathway Analysis
#================================================================================================== 
  ,tabPanel("Pathway",value = 8,
    sidebarLayout(
  
      # sidebar of Pathway --------------------------------------------------------------------------------
      sidebarPanel(
        htmlOutput("listComparisonsPathway")
        ,tags$style(type='text/css', "#listComparisonsPathway { width:100%;   margin-top:-12px}")
        ,selectInput("pathwayMethod", 
                     label = "Select method:", 
                     choices = list("GAGE"                   = 1, 
                                    "GSEA (preranked fgsea)" = 3,
                                    "PGSEA"                  = 2, 
                                    "PGSEA w/ all samples"   = 4, 
                                    "ReactomePA"             = 5),
                     selected = 1) #
        ,tags$style(type='text/css', "#pathwayMethod { width:100%;   margin-top:-12px}")
        #,h5("Select genesets (use KEGG to show pathway diagrams w/ fold-change):")

        ,htmlOutput("selectGO1")
        ,tags$style(type='text/css', "#selectGO1 { width:100%;   margin-top:-12px}")
        #,h5("Filter genesets by size:")
        ,fluidRow( 
          column(6, numericInput( "minSetSize", 
                                  label = h5("Geneset size: Min."), 
                                  min   = 1, 
                                  max   = 30, 
                                  value = 5,
                                  step  = 1) ),
          column(6, numericInput( "maxSetSize", 
                                  label = h5("Max."), 
                                  min   = 1000, 
                                  max   = 2000, 
                                  value = 2000,
                                  step  = 100) )
        ) # fluidRow
        ,tags$style(type='text/css', "#minSetSize { width:100%;   margin-top:-12px}")
        ,tags$style(type='text/css', "#maxSetSize { width:100%;   margin-top:-12px}")

        ,numericInput("pathwayPvalCutoff", 
                      label = h5("Pathway signifiance cutoff (FDR)"),
                      value = 0.2,
                      min   = 1e-20,
                      max   = 1,
                      step  = .05)
        ,tags$style(type='text/css', "#pathwayPvalCutoff { width:100%;   margin-top:-12px}")

        ,numericInput("nPathwayShow", 
                      label = h5("Number of top pathways to show"),
                      value = 30, 
                      min   = 5,
                      max   = 100,
                      step  = 5)
        ,tags$style(type='text/css', "#nPathwayShow { width:100%;   margin-top:-12px}")

        ,checkboxInput("absoluteFold", label = "Use absolute values of fold changes for GSEA and GAGE", value = FALSE)
        ,numericInput("GenePvalCutoff", 
                      label = h5("Remove genes with big FDR before pathway analysis:"),
                      value = 1,
                      min   = 1e-20,
                      max   = 1,
                      step  = .05)
        ,tags$style(type='text/css', "#GenePvalCutoff { width:100%;   margin-top:-12px}")

        # if pathway analysis using methods other than ReactomePA
        ,conditionalPanel("input.pathwayMethod == 1 | input.pathwayMethod == 2 |
                           input.pathwayMethod == 3| input.pathwayMethod == 4" 
          ,actionButton("ModalEnrichmentPlotPathway", "Pathway tree") 
          ,actionButton("ModalVisNetworkPA", "Network(New!)" )
        ,tags$head(tags$style("#ModalVisNetworkPA{color: red}"))   
          #,actionButton("ModalExaminePathways", "Gene expression by pathway")
          ,downloadButton('downloadPathwayListData', "Pathway list w/ genes")          
        )

        ,conditionalPanel("input.pathwayMethod == 4", 
                          downloadButton("PGSEAplotAllSamples.Download","High-resolution figure") )
        ,conditionalPanel("input.pathwayMethod == 2", 
                          downloadButton("PGSEAplot.Download", "High-resolution figure") )

        #,actionButton("examinePathway", "Examine individual pathways")
        #,conditionalPanel("input.pathwayMethod == 2",downloadButton('download.PGSEAplot.data', 'Download PGSEA pathway data'))
        ,h5("* Warning! The many combinations can lead to false positives in pathway analyses.")        
        ,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/pathways/",target="_blank")    
      ), # end of sidebar



      # main panel of pathway analysis -----------------------------------------------------------------------------------
      mainPanel(        
        conditionalPanel("input.pathwayMethod == 2",
                         h5("Red and blue indicates activated and suppressed pathways, respectively. ")          
                         ,plotOutput("PGSEAplot", inline=TRUE))
          
        ,conditionalPanel("input.pathwayMethod == 1", tableOutput("gagePathway") )
        ,conditionalPanel("input.pathwayMethod == 3", tableOutput("fgseaPathway") )
        ,conditionalPanel("input.pathwayMethod == 5", tableOutput("ReactomePAPathway") )
        ,conditionalPanel("input.pathwayMethod == 4",
                          h5("Red and blue indicates activated and suppressed pathways, respectively.  ") 
                          ,plotOutput("PGSEAplotAllSamples", inline=TRUE)
        )
          
        #,bsModal("ModalExaminePathway", "Gene expression on pathways", "ModalExaminePathways", size="large"          
        ,conditionalPanel("input.pathwayMethod == 1 | input.pathwayMethod == 2
                           | input.pathwayMethod == 3|input.pathwayMethod == 4" 
                          ,htmlOutput("listSigPathways")
                          ,downloadButton('downloadSelectedPathwayData', 'Expression data for genes in selected pathway')  )

        ,conditionalPanel(" (input.pathwayMethod == 1 | input.pathwayMethod == 2| 
                           input.pathwayMethod == 3 |input.pathwayMethod == 4) & input.selectGO == 'KEGG'"
                          ,h5("Red and green represent up- and down-regulated genes, respectively.")
                          ,imageOutput("KeggImage", width = "100%", height = "100%"))

        ,conditionalPanel(" (input.pathwayMethod == 1 | input.pathwayMethod == 2 
                           |input.pathwayMethod == 3 |input.pathwayMethod == 4) & input.selectGO != 'KEGG'",
                          plotOutput("selectedPathwayHeatmap"))        

        ,bsModal("ModalEnrichmentPlotPathway1", "Significant pathways", "ModalEnrichmentPlotPathway", size="large"
                 ,h5("Gene sets closer on the tree share more genes. Sizes of dot correspond to adjuested Pvalues")
                 ,downloadButton('enrichmentPlotPathway4Download',"Figure")
                 ,plotOutput('enrichmentPlotPathway') )

        ,bsModal("ModalEnrichmentPlotPahtway2", "Significant pathways", "ModalEnrichmentNetworkPathway", size="large"
                 ,h5("Connected gene sets share more genes. Color of node correspond to adjuested Pvalues.")            
                 ,actionButton("layoutButton3", "Change layout")
                 ,plotOutput('enrichmentNetworkPlotPathway')
      
        )#bsModal    


      # visNetwork ------------------------------
        ,bsModal("ModalVisNetworkPA1", "Related pathways", "ModalVisNetworkPA", size="large"
          ,h5("Connected gene sets share more genes. Color of node correspond to adjuested Pvalues.")            
          ,fluidRow(
            column(2, actionButton("layoutVisPA", "Change layout") ),
            column(1, h5("Cutoff:"), align="right" ) ,
            column(2, numericInput("edgeCutoffPA", label = NULL, value = 0.30, min = 0, max = 1, step = .1), align="left"  ), 
            column(2, checkboxInput("wrapTextNetworkPA", "Wrap text", value = TRUE)), 
            column(1, downloadButton("visNetworkPADownload","Network") ),
            column(1, downloadButton("downloadEdgesPA", "Edges")) ,
            column(1, downloadButton("downloadNodesPA", "Nodes"))
           )    
		   ,selectInput("upORdownRegPA", NULL,
						 c("Both Up & Down" = "Both",
						   "Up regulated" = "Up",
						   "Down regulated" = "Down")) 

           ,h6("This interactive plot also shows the relationship between enriched pathways. 
			   Two pathways (nodes) are connected if they share 30% (default, adjustable) or more genes. 
			   Green and red represents down- and up-regulated pathways. You can move the nodes by dragging them, zoom in and out by scrolling, 
			   and shift the entire network by click on an empty point and drag. 
			   Darker nodes are more significantly enriched gene sets. 
			   Bigger nodes represent larger gene sets.  
			   Thicker edges represent more overlapped genes. ")    
           ,visNetworkOutput("visNetworkPA",height = "800px", width = "800px")
        ) #bsModal
       # end VisNetwork --------------------------------    


       
      ) # mainPanel
    )  #sidebarLayout     
  ) #tabPanel  

#==================================================================================================
#  KEGG    4/24/2022
#==================================================================================================
  ,tabPanel("KEGG", value = 9,
    sidebarLayout(
  
      # sidebar of Bicluster --------------------------------------------------------------------------------
      sidebarPanel(
        h5("Visualize your fold-changes of all genes on all KEGG pathways"),
        htmlOutput('keggPathwaysAll'),        
        htmlOutput("listComparisonsKEGG"),        
        selectInput("maxFCKEGG", "Fold-change (log2) cutoff in color code", 
            choices = c(0.5, 1, 1.5, 2, 3, 4),
            selected = 2),  
        "Please cite the papers for ", 
                  a("pathview, ", href="https://doi.org/10.1093/bioinformatics/btt285"), "and ",
                  a("KEGG.", href="https://doi.org/10.1093/nar/gkaa970") 
      ), # end of sidebar


      # main panel of Biocluster ---------------------------------------------------------------------------
      mainPanel(
        h5("Bright red indicates most upregulated; bright green, most downregulated."),
       imageOutput("KeggImage_temp", width = "100%", height = "100%")

        
      ) # mainPanel
    )  #sidebarLayout     
  ) #tabPanel  


#================================================================================================== 
#  Genome view
#================================================================================================== 
  ,tabPanel("Genome",value = 10,
    sidebarLayout(  
      # sidebar of Genome --------------------------------------------------------------------------------
      sidebarPanel( 
        htmlOutput("listComparisonsGenome")
        ,tags$style(type='text/css', "#listComparisonsPathway { width:100%;   margin-top:-12px}")  
        ,fluidRow(
          column(6, numericInput("limmaPvalViz", 
                                 label = h5("Genes: FDR "), 
                                 value = 0.1,
                                 min   = 1e-5, 
                                 max   = 1,
                                 step  = .05) )
          ,column(6, numericInput( "limmaFCViz", 
                                   label = h5("Fold change"), 
                                   value = 2,
                                   min   = 1,
                                   max   = 100,
                                   step  = 0.5 ) )
        ) # fluidRow  
        ,tags$style(type='text/css', "#limmaPvalViz { width:100%;   margin-top:-12px}")
        ,tags$style(type='text/css', "#limmaFCViz { width:100%;   margin-top:-12px}")    
        
     

        ,fluidRow( 
                  column(6, checkboxInput("labelGeneSymbol", "Label Genes", value = FALSE) ) 
                  ,column(6, checkboxInput("ignoreNonCoding", "Coding genes only", value = TRUE) )  
                  )
        ,HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />') 

                  ,fluidRow(
                    column(6, selectInput(inputId = "MAwindowSize",
                                           label = h5("Window Size (Mb)"),
                                           selected = 6,
                                           choices = c(1, 2, 4, 6, 8, 10, 15, 20) ))
                    ,column(6, selectInput(inputId = "MAwindowSteps",
                                           label = h5("Steps"),
                                           selected = 2,
                                           choices = c(1, 2, 3, 4)))
                   )
                  ,selectInput(inputId = "chRegionPval", 
                                           label = h5("FDR cutoff for window"),
                                           selected = 0.0001,
                                           choices = c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001)) 
        ,HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />')  

        ,actionButton("runPREDA", "Run PREDA (5 mins)")
        ,br()
        ,div(style = "display:inline-block; float:right", actionButton("genomeTabButton", "Info"))
      ), # end of sidebar
    
      
      # main panel of Genome tab -----------------------------------------------------------------------------------
      mainPanel(                   
                  bsModal("genomeTabButtonID", "View your genes on chromosomes", "genomeTabButton", size = "large"
                          ,h3("Where are your differentially expressed genes (DEGs) located on the genome?" )
                          ,p("Red and blue dots represent significantly up- and down-regulated genes, respectively, according to the criteria on the side panel. These criteria could differ from the one in DEG1 tab. The distance of the dots from the closest chromosome is proportional to the log2 fold-change (FC).")

                         ,h3("Are there regions of the genome where genes are coherently up- or down-regulated?")
                         ,p("To answer this question, we scan the genome with sliding windows. Within each window we take several steps to slide forward. For example if you choose a window size = 6Mbps and steps = 3, the first window is from 0 to 6 Mbps, the 2nd  from 2 to 8Mbps, and the third from 4 to 10 Mbps, and so on.")
                         ,p("For all genes in a window/region, we test whether the mean of FC of these genes is zero using a t-test. All genes analyzed by DESeq2 or limma, significant or otherwise, are included in this analysis. Hence this result is indepdent of our DEG cutoffs. P values from the test of the mean are adjusted to FDR. Essentially, we considered genes located in a genomic region as a gene set or pathway, and we performed simple pathway analysis by asking whether these genes are behaving non-randomly.")

                         ,p("Based on an FDR cutoff for the windows, red and blue segments indicate genomic regions with genes coherently up- or down-regulated, respectively.  Below you can adjust the window size, and steps in a window, and FDR cutoff for windows.  Mouse over to see gene symbols or IDs. Zoom in regions of interest. The chromosomes may be only partly shown as we use the last gene's location to draw the line.")

                         ,p("As an alternative approach, you can use "
                             ,a("PREDA.", href="https://doi.org/10.1093/bioinformatics/btr404",target="_blank"),
                             "Very slow (5 mins), but may be useful in studying cancer or
                             other diseases that might involve chromosomal gain or loss."  ) 
                   )
           

                  ,plotlyOutput("genomePlotly",height = "900px")

  
        ,bsModal("modalExample111", "Differentially expressed genomic loci", "runPREDA", size="large"
                 ,fluidRow( 
                           column(3, numericInput("RegionsPvalCutoff", 
                                                  label = h5("Min. FDR"), 
                                                  value = 0.01,
                                                  min   = 1e-20,
                                                  max   = 1,
                                                  step  = .05) ),
                           column(3, numericInput("StatisticCutoff",
                                                  label = h5("Min. Statistic"),
                                                  min   = .2, 
                                                  max   = 1.5, 
                                                  value = .5,
                                                  step  =.1) ),

                           column(3,actionButton("showRegions", "Significant Loci")),
                           column(3,actionButton("showGenesRegions", "Genes"))
                 ) #fluidRow
                 ,tags$style(type='text/css', "#showRegions { width:100%; margin-top: 40px;}")
                 ,tags$style(type='text/css', "#showGenesRegions { width:100%; margin-top: 40px;}")
                 ,plotOutput("genomePlot")
        ) #bsModal

        ,bsModal("modalExample15", "Diff. expressed regions", "showRegions", size = "large"
                 ,downloadButton('downloadRegions', 'Download')
                 ,dataTableOutput("chrRegions"))

        ,bsModal("modalExample16", "Diff. expressed regions", "showGenesRegions", size = "large"
                 ,downloadButton('downloadGenesInRegions', 'Download')
                 ,dataTableOutput("genesInChrRegions"))

      ) # mainPanel
    )  #sidebarLayout     
  ) #tabPanel 

 
#==================================================================================================
#  Bicluster
#==================================================================================================
  ,tabPanel("Bicluster",value = 11,
    sidebarLayout(
  
      # sidebar of Bicluster --------------------------------------------------------------------------------
      sidebarPanel(
        h5("Biclustering can discover genes correlated on subset of samples. 
            Only useful when  sample size is large(>10). 
            Uses methods implemented in the biclust R package. ")

        ,numericInput("nGenesBiclust", 
                      label = h5("Most variable genes to include "), 
                      min   = 10, 
                      max   = 2000, 
                      value = 1000) 

        ,selectInput("biclustMethod", "Method:", 
                      choices = list( 
                                      "BCCC"        = "BCCC()"
                                      ,"QUBIC"      = "BCQU()"
                                      ,"runibic"    = "BCUnibic()"                                    
                                      ,"BCXmotifs"  ="BCXmotifs()"
                                      , "BCPlaid"   = "BCPlaid()"
                                      ,"BCSpectral" = "BCSpectral()"
                                      ,"BCBimax"    = "BCBimax()"
                                      ,"BCQuest"    = "BCQuest()" ),
                     selected = "BCCC()")        
  
        ,htmlOutput('listBiclusters')
        ,h5("Enrichment database")
        ,htmlOutput("selectGO4"),tags$style(type='text/css', "#selectGO4 { width:100%;   margin-top:-9px}")              
        ,downloadButton('download.biclust.data', 'Download all biclusters')
        ,br(),br(),textOutput('biclusterInfo')
        ,tags$head(tags$style("#biclusterInfo{color: blue;
                       font-size: 16px;
                       font-style: italic;
                       }"))          
        ,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/biclustering/",target="_blank")
      ), # end of sidebar


      # main panel of Biocluster ---------------------------------------------------------------------------
      mainPanel(  
        plotOutput('biclustHeatmap')
        ,h3("Enriched gene sets in selected bicluster")
        ,tableOutput('geneListBclustGO')
        ,br(),br()
        ,h3("Genes in this cluster")
        ,tableOutput('geneListBicluster')
        
      ) # mainPanel
    )  #sidebarLayout     
  ) #tabPanel  

#==================================================================================================
#  Network
#==================================================================================================
  ,tabPanel("Network",value = 12,
    sidebarLayout(
  
      # sidebar of Network --------------------------------------------------------------------------
      sidebarPanel(
        h5("Identify co-expression networks and sub-modules using",
           a( "WGCNA.", href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559",
             target="_blank" )
           ,"Only useful when  sample size is large(>15). ")
        ,numericInput("nGenesNetwork", 
                      label = h5("Most variable genes to include (<3001) "), 
                      min   = 10, 
                      max   = 3000, 
                      value = 1000) 
        ,fluidRow(
          column(6, numericInput("mySoftPower", 
                                 label = h5("Soft Threshold"), 
                                 min = 1, 
                                 max = 20, 
                                 value = 5))
          ,column(6, numericInput("minModuleSize",
                                  label = h5("Min. Module Size"),
                                  min = 10, 
                                  max = 100, 
                                  value = 20)  )
        ) # fluidRow  
        ,tags$style(type='text/css', "#mySoftPower { width:100%;   margin-top:-12px}")
        ,tags$style(type='text/css', "#minModuleSize { width:100%;   margin-top:-12px}")
        
        ,actionButton("chooseSoftThreshold", "Choose soft threshold")  
        ,actionButton("showModuleHeatmap", "Heatmap")          

        ,textOutput('moduleStatistics')
        ,tags$head(tags$style("#moduleStatistics{color: blue;
                       font-size: 14px;
                       font-style: italic;
                       }" ))           
                          
        ,HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />') # a solid line
        ,htmlOutput('listWGCNA.Modules')
        ,fluidRow(
          column(6, numericInput("edgeThreshold", 
                                 label = h5("Edge Threshold"), 
                                 min   = 0, 
                                 max   = 1, 
                                 value = .4, 
                                 step  = .1))
          ,column(6, numericInput("topGenesNetwork", 
                                  label = h5("Top genes"), 
                                  min   = 10, 
                                  max   = 2000, 
                                  value = 10, 
                                  step  = 10)  )
        ) # fluidRow  
        ,tags$style(type='text/css', "#edgeThreshold { width:100%;   margin-top:-12px}")
        ,tags$style(type='text/css', "#topGenesNetwork { width:100%;   margin-top:-12px}")

        ,actionButton("networkLayout", "Change network layout", style="float:center")        
        ,h5("Enrichment database:")
        ,htmlOutput('selectGO5')
        ,downloadButton('download.WGCNA.Module.data',"Download all modules")
        ,downloadButton('downloadSelectedModule',"Download network for selected module")  
        ,h5("The network file can be imported to", 
            a(" VisANT", href="http://visant.bu.edu/",target="_blank"),
            " or ", 
            a("Cytoscape.",href="http://www.cytoscape.org/",target="_blank" )  )  

        ,a(h5("?",align = "right"), href="https://idepsite.wordpress.com/network/",target="_blank")          
      ), # end of sidebar


      mainPanel(            
        plotOutput('modulePlot')
        ,br(),br()
        ,plotOutput('moduleNetwork')
        ,h3("Enriched pathways among all genes in selected module")          
        ,tableOutput('networkModuleGO')

        ,bsModal("modalExample112", "Choose soft threshold", "chooseSoftThreshold", 
                 size="large"
                 ,plotOutput('softPower') )

        ,bsModal("modalExample116", "Heatmap of identified modules", "showModuleHeatmap",
                 size="large"
                 ,plotOutput('networkHeatmap') )

       ) # mainPanel
    )  #sidebarLayout     
  ) #tabPanel  

#==================================================================================================
#  Reproducibility
#================================================================================================== 
  ,tabPanel("R",value = 13,
    fluidRow(    
      column(12,
        h4( a("Email us", href= "mailto:Xijin.Ge@SDSTATE.EDU?Subject=iDEP", target="_top") , 
           " for questions, suggestions, or data contributions. Stay connected via ", 
           a("user group", href="https://groups.google.com/d/forum/idep",target="_blank"),
           " or ",a("Twitter.", href="https://twitter.com/StevenXGe", target="_blank"),
           " Visit our ", a("GitHub", href="https://github.com/iDEP-SDSU/idep", target="_blank"), 
           " page to see source code, install a local version, or report bugs and request features.", 
           "iDEP is being developed by a very small team: Dr. Xijin Ge and a graduate student ",
           a("(Homepage).", href="http://ge-lab.org/")
           ),
        h2("R as in Reproducibility"),
        
        h5(a("Documentation site.",href="https://idepsite.wordpress.com/", target="_blank"),
        " Source code on ", a("GitHub,", href="https://github.com/iDEP-SDSU/idep", target="_blank"),
           "where users can also report bugs or requet features."
           ),
        

        h5("To improve reproducibility, iDEP generates custom  R code  
            based on your data and choices of parameters. Users with some
            R coding experience should be able to re-run most analyses 
            by downloading all of the files below. If Ensembl IDs is not used in users' original file, 
            we should use the converted data file. Click through all the tabs and then download all
            these file to a folder. Run the Customized R code or the Markdown file. "
            ,a("R Markdown example.",align = "left", href="https://gex.netlify.com/post/reproducing-idep-analyses-with-auto-generated-r-markdown/"
             ,target="_blank"))
    
        ,downloadButton('downloadRcode',"Customized R code")   
        ,downloadButton('downloadRcodeMarkdown',"Customized R code(Markdown)")  
        ,downloadButton('downloadRfunctions',"iDEP core functions")
        ,downloadButton('downloadGeneInfo',"Gene Info file")
        ,downloadButton('downloadPathwayFile',"Pathway file (large)")
        ,conditionalPanel("input.dataFileFormat == 1",     
                           downloadButton('downloadConvertedCountsRtab'
                          ,"Converted counts"))
        ,conditionalPanel("input.dataFileFormat == 2",     
                          downloadButton('downloadProcessedDataRtab'
                          ,"Converted data"))
        ,conditionalPanel("input.file2 > 0 | input.goButton > 0",        
                          downloadButton('downloadSampleInfoFile'
                          ,"Experiment design file"))
        ,br()
        ,h4("Previous versions of iDEP")
        ,a("iDEP 0.94 with Ensembl Release 104, archived on Feb. 8, 2022 "
        , href="http://bioinformatics.sdstate.edu/idep94/")

        ,a("iDEP 0.93 with Ensembl Release 103, archived on Oct. 15, 2021 "
           , href="http://bioinformatics.sdstate.edu/idep93/")  
        ,br()
        ,a("iDEP 0.92 with Ensembl Release 100, archived on May 20, 2021 "
           , href="http://bioinformatics.sdstate.edu/idep92/")  
        ,br()
        ,a("iDEP 0.90 with Ensembl Release 96, archived on May 20, 2021 "
           , href="http://bioinformatics.sdstate.edu/idep90/")  
        ,br()
        ,a("iDEP 0.85 with Ensembl Release 95, archived on May 19, 2019 "
           , href="http://bioinformatics.sdstate.edu/idep85/")  
        ,br()
        ,a("iDEP 0.82 with Ensembl  Release 92, archived on March 29, 2019 "
           , href="http://bioinformatics.sdstate.edu/idep82/")  
        ,br()
        ,a("iDEP 0.73 with Ensembl  Release 91, archived on July 11, 2018 "
            , href="http://bioinformatics.sdstate.edu/idep73/")  
        ,br()
        ,h4("Citation")  
        ,h5("Please cite: Ge SX, Son EW, Yao R: iDEP: an integrated web application for differential 
            expression and pathway analysis of RNA-Seq data. BMC Bioinformatics 2018, 19(1):534. PMID:30567491 ",
            a(" Full text",href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2486-6", target="_blank"))
        ,h4("Usage Statistics")
        ,h5("As of May 11 2021, iDEP website has been visited 167832 times by 43,220 users from 100+ countries. 
            iDEP has been used by researchers to analyze transcriptomic data 
            from sun flowers to primates to human, 
            as demonstrated by the 191 papers citing iDEP.")
       ,br()  
       ,htmlOutput('RsessionInfo')

       ,br(), br(), h4("Changes")
       ,h5("v 0.38 Gene ID conversion, remove redudancy;  rlog option set to blind=TRUE.")
       ,h5("v 0.39 Reorganized code. Updated to Bioconductor 3.5; solved problems with PREDA 9/8/17.")
       ,h5("v 0.40 Moved libraries from the beginning to different places to save loading time.")
       ,h5("v 0.50 Enable changing colors for heatmap. Fixed error with voom; also changed method to TMM, see", a( "here.",href=
         "https://www.bioconductor.org/help/workflows/RNAseq123/#transformations-from-the-raw-scale"))
 
       ,h5("2/5/2018 V 0.62 :  Add  tree and networks to visualize 
            overlaps between enriched gene sets, in K-means, DEG2, and pathway tags.
            Downloads of pathway analysis results and high-resolution figures.")
       ,h5("2/6/2018: Fixed errors caused by gene symbol matching for unknown species.
             More user control of hierarchical clustering tree")
       ,h5("2/9/2018: V 0.65 Added API access to STRINGdb website on the DEG2 tab. 
            Supports thousands of bacterial species")
       ,h5("2/10/2018: V 0.66 Improved API access to STRINGdb, by adding automatic species matching.")
       ,h5("2/11/2018: V 0.67 Tested with larger dataset of 259 samples. Changed figure configurations.")
       ,h5("2/14/2018: V 0.68 Fixed Pathview loading code. Connected Pathview to PGSEA.")
       ,h5("2/25/2018: V 0.68 Fixed Fold-change, FDR data upload and parsing.
            Figure resolution using the res=150 option help improve readability of labels.")
       ,h5("2/28/2018: V0.69 Change interactive heat maps to re-order columns")
       ,h5("3/12/2018: V0.70 Generating R code and downloading annotation files used in analysis.")
       ,h5("3/14/2018: V0.71 Improve R markdown file; add color to EDA plots; detect bias in sequencing depth ")
       ,h5("3/18/2018: v0.711 Fixed error caused by gene names containing characters such as \' or \" ")
       ,h5("3/28/2018: v0.712 Fine tuned EDA plots")
       ,h5("4/3/2018: V0.713 add permutations for fgsea. Expand quotes.")
       ,h5("4/13/2018: V0.72 fixed bug caused error in R Markdown file when users choose species other than the default.")
       ,h5("4/25/2018: V0.73 Fixed bug for GO terms tree in k-Means. Added layout button for network viz of GO terms")
       ,h5("5/28/2018: V0.73 Enabled download of all FDR and fold-changes at the DEG1 tab. 
           Problem with the FDR and fold-change data with only one columns. ")
       ,h5("7/28/2018: v0.80 New v0.80  Updated annotation database. Comprehensive pathway 
          database for human. TF binding motifs for 200+ speceis. Old version made available.")
       ,h5("7/29/2018: v0.80 Orgized UI.R according to", a("Google R style guide.", 
            href="https://google.github.io/styleguide/Rguide.xml") )   
       ,h5("7/30/2018: V0.80 Fixed error in downloaded up- and down-regulated gene lists. We thank Juan Xie for pointing this out. ")    
       ,h5("12/2/2018: v0.81 High resolution figure download with eps format, which can be eidted with Adobe Illustrator.")
       ,h5("3/5/2019:  v0.82 Fix a bug regarding limma for identification of D.E.Gs. Up- and down-regulation are opposite in some cases.")
       ,h5("3/29/2019: v0.85 Annotation database upgrade. Ensembl v 95. Ensembl plants v.42, and Ensembl Metazoa v.42.")
       ,h5("5/19/2019: v0.90 Annotation database upgrade. Ensembl v 96. Ensembl plants v.43, and Ensembl Metazoa v.43. STRING-db v10")
       ,h5("2/3/2020: v0.90 customizable PCA plot and scatter plot")
       ,h5("5/10/2021: V0.93 updated to Ensembl Release 103 and String-DB v11.")
       ,h5("10/15/21: For GO enrichment analysis, we now recommend using background genes, instead of all genes on the genome. In the KNN and DEG2 tabs, it is now the default that all genes passed the filter in Pre-Process tab are used as a customized background.")
      ,h4("10/18/20: Interactive network enabling users to easily visualize the relatedness of significant pathways similar to EnrichmentMap. ")
     ,p("April 23, 2022: KEGG diagram on the Pathway tab is now working! Please email Jenny if you see any issues!")
      ,p("April 23, 2022: DEG2 tab is now much faster! I hope people complained about it earlier!")
      ,p("April 20, 2022: Plot one or more genes is now working. Bug fixed.")
      ,p("Mar. 7, 2022: Fixed an R library issue affected KEGG diagrams for some organisms.")
      ,p("Feb. 19, 2022: R upgraded from 4.05 to 4.1.2. This solved the STRING API issues. Some Bioconductor packages are also upgraded.", style = "color:red")
      ,p("Feb. 8, 2022: iDDEP v0.95 becomes default. Old versions are still available. See the last tab.")
      ,p("Nov. 15, 2021: iDEP v0.95 available in testing mode. It includes Ensembl database update, new species from Ensembl Fungi and Ensembl Protists, and STRINGdb (5090 species) update from v11 to 11.5.")
      ,p("10/26/2021: The Genome view is now much improved! Automatically detects chromosomal regions enriched with genes having abnormaly high and low fold-changes.")

      ,p("10/15/21: For GO enrichment analysis, we now recommend using background genes, instead of all genes on the genome. In the KNN and DEG2 tabs, it is now the default that all genes passed the filter in Pre-Process tab are used as a customized background.")

      ,p("We updated instruction for local installation", a("here.", href="https://github.com/iDEP-SDSU/idep#readme"), 
          "The most recent database file is now publically available, free of charge for non-profit organizations.")

      ,p("Check out the 50,000+ datasets of uniformly processed public RNA-seq data ", a("here!", href="http://bioinformatics.sdstate.edu/reads/" ))
     ,p("10/18/20: Interactive network enables users to easily visualize the relatedness 
           of pathways, similar to EnrichmentMap. Using the Network buttons on DEG2 and Pathway tabs,
           you can generate and export interactive networks like this one below. You can move the nodes by dragging them, zoom in by scrolling, 
			   and shift the entire network by click on an empty point and drag. ")    
       ,br(),br()
       ,h5("In loving memory of my parents. X.G.")

       ) #column
     ) # fluidRow
   ) # tabPanel

  ,tags$head(includeScript("ga.js")) # tracking usage  
  )# Navibar

) #shinyUI
