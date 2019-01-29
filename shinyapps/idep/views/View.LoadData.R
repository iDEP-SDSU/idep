library('R6')

View.LoadData <- R6Class("View.LoadData")

View.LoadData$set(  "public", 
                    "sidebarLayout", 
                    function(){
                            # sidebar---------------------------------
                            sidebarPanel(
                                actionButton("goButton", "Click here to load demo data"),
                                tags$head(tags$style("#goButton{color: red;
                                                        font-size: 16px;
                                                        font-style: italic;
                                                        }"))                    
                                ,h5(" and just click the tabs for some magic!", style = "color:red")
                                ,p(HTML("<div align=\"right\"> <A HREF=\"javascript:history.go(0)\">Reset</A></div>" ))
                                ,radioButtons("dataFileFormat", 
                                            label = "1. Choose data type", 
                                            choices = list("Read counts data (recommended)"                                          = 1, 
                                                            "Normalized expression values (RNA-seq FPKM, microarray, etc.)"          = 2,
                                                            "Fold-changes and corrected P values from CuffDiff or any other program" = 3),
                                            selected = 1
                                )      
                                ,conditionalPanel("input.dataFileFormat == 3",
                                    checkboxInput("noFDR", "Fold-changes only, no corrected P values", value = FALSE)
                                )
                            
                                ,fileInput('file1', '2. Upload expression data (CSV or text)',
                                        accept = c(
                                            'text/csv',
                                            'text/comma-separated-values',
                                            'text/tab-separated-values',
                                            'text/plain',
                                            '.csv',
                                            '.tsv'          
                                        ) 
                                )
                                ,a("New! Analyze public RNA-seq data", href="http://bioinformatics.sdstate.edu/reads/")
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

                                ,strong("3. Verify guessed species. Change if neccessary.")
                                ,selectInput("selectOrg", label = NULL,"Best matching species",width='100%') 
                                    
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

                                ,tableOutput('species' )
                                ,a( h5("?",align = "right"), href="https://idepsite.wordpress.com/data-format/",target="_blank")
                                                                                                                # new window
                            )
                    }
                )