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
)

