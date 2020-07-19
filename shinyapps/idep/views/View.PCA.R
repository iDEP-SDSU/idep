library('R6')
library(shiny)
library(shinyBS)
library(DT,verbose=FALSE) 		# for renderDataTable


View.PCA <- R6Class("View.PCA")

###############################################################################
###################				Main Panel					###################
###############################################################################
View.PCA$set("public", "mainPanel",
	function(){
		mainPanel(
            plotOutput("PCA_mainplot", inline=TRUE),
            self$ConPanel_TSNEMainPanel(),
            br(),
            br(),
            self$ConPanel_MDSorPCAMainPanel()
        )
    }
)





###############################################################################
###################				sidebar Panel				###################
###############################################################################

View.PCA$set("public", "sidebarPanel",
    function(){
        sidebarPanel(
            radioButtons("select_PCA_Methods", "Methods", 
                c(  "Principal Component Analysis"     = 1,            
                    "Multidimensional Scaling"         = 3, 
                    "t-SNE"                            = 4,
                    "Pathway Analysis of PCA rotation" = 2  
                )
            ),
            self$ConPanel_PathwaySettings(),
            self$ConPanel_NonPathwaySettings(),
            br(),
            downloadButton('download_PCA_Data', 'Coordinates'),
            downloadButton('download_PCA_Mainplot', 'High-resolution figure'),
            a(h5("?",align = "right"), href="https://idepsite.wordpress.com/pca/",target="_blank")
        )
    }
)


###############################################################################
###################			Component Functions				###################
###############################################################################

# side panel setting, conditional
View.PCA$set("public", "ConPanel_PathwaySettings",
    function(){
        conditionalPanel("input.select_PCA_Methods == 2",
            htmlOutput('select_PCA_GO'),  ## selectGO6 in previous version
            fluidRow( 
                column(6, numericInput( "num_PCA_minGeneSetSize", 
                    label = h5("Geneset size: Min."), 
                    min   = 5, 
                    max   = 30, 
                    value = 15,
                    step  = 1)
                ),
                column(6, numericInput( "num_PCA_maxGeneSetSize", 
                    label = h5("Max."), 
                    min   = 1000, 
                    max   = 2000, 
                    value = 2000,
                    step  = 100) 
                )
            ), # fluidRow
            tags$style(type='text/css', "#num_PCA_minGeneSetSize { width:100%;   margin-top:-12px}"),
            tags$style(type='text/css', "#num_PCA_maxGeneSetSize { width:100%;   margin-top:-12px}")
        )
    }
)

# side panel setting, conditional
View.PCA$set("public", "ConPanel_NonPathwaySettings",
    function(){
        conditionalPanel("input.select_PCA_Methods != 2",
            htmlOutput('div_PCA_MainPlotColorSelection'), 
            htmlOutput('div_PCA_MainPlotShapeSelection') 
        )
    }
)

# main panel conditional panel, for TSNE only
View.PCA$set("public", "ConPanel_TSNEMainPanel",
    function(){
        conditionalPanel("input.select_PCA_Methods == 4",
            actionButton("btn_PCA_RecalcTSNE", "Re-calculate t-SNE"),
            br(),
            br()
        )
    }
)

View.PCA$set("public", "ConPanel_MDSorPCAMainPanel",
    function(){
        conditionalPanel("input.select_PCA_Methods == 1 | input.select_PCA_Methods == 2",
            htmlOutput('div_PCA_CorrelationBetweenPCs')
        )
    }
)

