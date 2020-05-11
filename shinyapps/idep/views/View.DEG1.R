library('R6')
library(shiny)
library(shinyBS)
library(DT,verbose=FALSE) 		# for renderDataTable


View.DEG1 <- R6Class("View.DGE1")


###############################################################################
###################				Main Structure				###################
###############################################################################

View.DGE1$set("public", "sidebarPanel",
    function(){
        sidebarPanel(
            h5("Identifying Differential Expressed Genes (DEGs). See next tab for details."),
            self$Cond_DataIsReadCount(),
            self$Cond_DataIsNormalizedExpressionValues(),
            self$OtherSettings(),
            width = 4
        )
    }
)

View.DEG1$set("public", "mainPanel",
	function(){
		mainPanel(
            plotOutput('sigGeneStats'),
            br(),br(),
            h4("Numbers of differentially expressed genes for all comparisons.
            \"B-A\" means B vs. A. Interaction terms start with \"I:\" "),
            tableOutput('sigGeneStatsTable'),
            self$PopShowVennDiagram(),
            self$PopShowModelAndComparisons()
        )
    }
)



###############################################################################
###################			Component Functions				###################
###############################################################################

View.DEG1$set("public", "PopShowVennDiagram",
    function(){
        bsModal(
            "modalExample", 
            "Venn Diagram", 
            "showVenn", 
            size = "large",
            checkboxInput("UpDownRegulated", label = "Split gene lists by up- or down-regulation", value = FALSE),
            htmlOutput('listComparisonsVenn'),
            downloadButton("DownloadVenn", "Figure"),
            plotOutput("vennPlot"),
        )
    }
)

View.DEG1$set("public", "PopShowModelAndComparisons",
    function(){
        bsModal(
            "modalExample21", 
            "Build model and/or select comparisons", 
            "modelAndComparisons", 
            size = "large",
            fluidRow(
                column(6, htmlOutput('listFactorsDE')),
                column(6, htmlOutput('listBlockFactorsDE') ) 
            ),

            fluidRow(
                column(6, htmlOutput('selectReferenceLevels1')),
                column(6, htmlOutput('selectReferenceLevels2') ) 
            ),
            fluidRow(
                column(6, htmlOutput('selectReferenceLevels3')),
                column(6, htmlOutput('selectReferenceLevels4') ) 
            ),
            fluidRow(
                column(6, htmlOutput('selectReferenceLevels5')),
                column(6, htmlOutput('selectReferenceLevels6') ) 
            ),            
            htmlOutput('listInteractionTerms'),
            textOutput('experimentDesign')
            tags$head(tags$style("#experimentDesign{color: red;font-size: 16px;}")),
            htmlOutput('listModelComparisons'),
            actionButton("submitModelButton", "Submit & re-calculate",style="float:center"),
            tags$head(tags$style("#submitModelButton{color: blue;font-size: 20px;}")),
            h5("Close this window to see results."),
            a(  h5("More info on DESeq2 experiment design",
                align = "right"),
                href="http://rpubs.com/ge600/deseq2",
                target="_blank"
            )
        )
    }
)

View.DEG1$set("public", "Cond_DataIsReadCount",
    function(){
        conditionalPanel(
            "input.dataFileFormat == 1",
            selectInput("CountsDEGMethod", "Method:", 
                choices = list(
                    "DESeq2"      = 3,
                    "limma-voom"  = 2,
                    "limma-trend" = 1
                ),
                selected = 3
            ),       
            tags$style(type='text/css', "#CountsDEGMethod { width:100%;   margin-top:-12px}")
        )
    }
)

View.DEG1$set("public", "Cond_DataIsNormalizedExpressionValues",
    function(){
        conditionalPanel( 
            "input.dataFileFormat == 2", 
            h5("Using the limma package")        
        )
    }
)

View.DEG1$set("public", "OtherSettings",
    function(){

        fluidRow(
            column( 5,
                numericInput("limmaPval",  # FDR cutoff
                    label = h5("FDR cutoff"), 
                    value = 0.1,
                    min   = 1e-5,
                    max   = 1,
                    step  =.05
                ) 
            ),
            column(7,
                numericInput("limmaFC",   # fold change cutoff
                    label = h5("Min fold change"),
                    value = 2,
                    min   = 1,
                    max   = 100,
                    step  = 0.5
                ) 
            )
        ),
        # style of limmaPval and limmaFC
        tags$style(type='text/css', "#limmaPval { width:100%;   margin-top:-12px}"),
        tags$style(type='text/css', "#limmaFC { width:100%;   margin-top:-12px}"),

        actionButton("modelAndComparisons", "Select factors & comparisons"),
        tags$head(tags$style("#modelAndComparisons{color: blue;font-size: 15px;}")),
        br(),br(),

        fluidRow(
           column(6,actionButton("showVenn", "Venn Diagram") )
        ),
        br(),
        fluidRow(          
            column(8, downloadButton('downloadGeneListsGMT', 'Gene lists') ) 
        ),
        br(),
        fluidRow(
           column(8, downloadButton('download.DEG.data', 'FDR & fold-changes for all genes') )
        ),
        downloadButton('downloadSigGeneStats', 'Figure'),
        br(),
        h4( textOutput("textLimma") ),
        tags$head(tags$style("#textLimma{color: blue;font-size: 15px;}"))  ,
        a(h5("?",align = "right"), href="https://idepsite.wordpress.com/degs/",target="_blank")
    }
)



