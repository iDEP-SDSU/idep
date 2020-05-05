library('R6')
library(shiny)
library(shinyBS)
library(DT,verbose=FALSE) 		# for renderDataTable


View.PCA <- R6Class("View.PCA")

###############################################################################
###################				Main Panel					###################
###############################################################################
View.PCA$set("public", "sidebarPanel",
	function(){
		mainPanel(
            plotOutput("PCA_mainplot", inline=TRUE),
            self$ConPanel_tSNEMainPanel(),
            br(),
            br(),
            self$ConPanel_Need a name()
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
            downloadButton('download_PCA_Mainplot', 'Coordinates'),
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
            htmlOutput('select_PCA_GO')  ## selectGO6 in previous version
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
View.PCA$set("public", "ConPanel_tSNEMainPanel",
    function(){
        conditionalPanel("input.select_PCA_Methods == 4",
            actionButton("btn_PCA_RecalcTSNE", "Re-calculate t-SNE"),
            br(),
            br()
        )
    }
)

View.PCA$set("public", "ConPanel_Need a name",
    function(){
        conditionalPanel("CONDITIONNEEDCONFIRM",
            htmlOutput('PCA2factor') ## need rename
        )
    }
)

