#' mod_07_genome UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_07_genome_ui <- function(id){
  ns <- NS(id)
  tabPanel(
    "Genome",
    sidebarLayout(  
      sidebarPanel( 
        htmlOutput(outputId = ns("list_comparisons_genome")),
        tags$style(
          type = "text/css",
          "#genome-list_comparisons_genome{ width:100%;   margin-top:-12px}"
        ),
        fluidRow(
          column(
            width = 6,
            numericInput(
              inputId = ns("limma_p_val_viz"), 
              label = h5("Genes: FDR "), 
              value = 0.1,
              min = 1e-5, 
              max = 1,
              step = .05
            )
          ),
          column(
            width = 6,
            numericInput(
              inputId = ns("limma_fc_viz"), 
              label = h5("Fold change"), 
              value = 2,
              min = 1,
              max = 100,
              step = 0.5
            )
          )
        ),
        tags$style(
          type = "text/css",
          "#genome-limma_p_val_viz{ width:100%;   margin-top:-12px}"
        ),
        tags$style(
          type = "text/css",
          "#genome-limma_fc_viz{ width:100%;   margin-top:-12px}"
        ),
        fluidRow( 
          column(
            width = 6,
            checkboxInput(
              inputId = ns("label_gene_symbol"),
              label = "Label Genes",
              value = FALSE
            )
          ),
          column(
            width = 6,
            checkboxInput(
              inputId = ns("ignore_non_coding"),
              label = "Coding genes only",
              value = TRUE
            )
          )  
        ),
        HTML(
          "<hr style='height:1px;border:none;
          color:#333;background-color:#333;' />"
        ),
        fluidRow(
          column(
            width = 6,
            selectInput(
              inputId = ns("ma_window_size"),
              label = h5("Window Size (Mb)"),
              selected = 6,
              choices = c(1, 2, 4, 6, 8, 10, 15, 20)
            )
          ),
          column(
            width = 6,
            selectInput(
              inputId = ns("ma_window_steps"),
              label = h5("Steps"),
              selected = 2,
              choices = c(1, 2, 3, 4)
            )
          )
        ),
        selectInput(
          inputId = ns("ch_region_p_val"), 
          label = h5("FDR cutoff for window"),
          selected = 0.0001,
          choices = c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001)
        )
      ),
      mainPanel(    
        tabsetPanel(
          tabPanel(
            "Info",
            h3(
              "Where are your differentially expressed genes (DEGs)
              located on the genome?"
            ),
            p(
              "Red and blue dots represent significantly up- and
              down-regulated genes, respectively, according to the
              criteria on the side panel. These criteria could
              differ from the one in DEG1 tab. The distance of the
              dots from the closest chromosome is proportional to
              the log2 fold-change (FC)."
            ),
            h3(
              "Are there regions of the genome where genes are
              coherently up- or down-regulated?"
            ),
            p(
              "To answer this question, we scan the genome with
              sliding windows. Within each window we take several
              steps to slide forward. For example if you choose a
              window size = 6Mbps and steps = 3, the first window is
              from 0 to 6 Mbps, the 2nd  from 2 to 8Mbps, and the
              third from 4 to 10 Mbps, and so on."
            ),
            p(
              "For all genes in a window/region, we test whether the
              mean of FC of these genes is zero using a t-test. All
              genes analyzed by DESeq2 or limma, significant or
              otherwise, are included in this analysis. Hence this
              result is indepdent of our DEG cutoffs. P values from
              the test of the mean are adjusted to FDR. Essentially,
              we considered genes located in a genomic region as a
              gene set or pathway, and we performed simple pathway
              analysis by asking whether these genes are behaving
              non-randomly."
            ),
            p(
              "Based on an FDR cutoff for the windows, red and blue
              segments indicate genomic regions with genes coherently
              up- or down-regulated, respectively. Below you can
              adjust the window size, and steps in a window, and FDR
              cutoff for windows.  Mouse over to see gene symbols or
              IDs. Zoom in regions of interest. The chromosomes may
              be only partly shown as we use the last gene's location
              to draw the line."
            ),
            p(
              "As an alternative approach, you can use ",
              a(
                "PREDA.",
                href = "https://doi.org/10.1093/bioinformatics/btr404",
                target = "_blank"
              ),
              "Very slow (5 mins), but may be useful in studying
              cancer or other diseases that might involve chromosomal
              gain or loss."
            ) 
          ),
          tabPanel(
            "Chromosome Plot",
            plotly::plotlyOutput(
              outputId = ns("genome_plotly"),
              height = "900px"
            )
          ),
          tabPanel(
            "PREDA (5 Mins)",
            fluidRow( 
              column(
                width = 3,
                numericInput(
                  inputId = ns("regions_p_val_cutoff"), 
                  label = h5("Min. FDR"), 
                  value = 0.01,
                  min = 1e-20,
                  max = 1,
                  step = .05
                )
              ),
              column(
                width = 3,
                numericInput(
                  inputId = ns("statistic_cutoff"),
                  label = h5("Min. Statistic"),
                  min = .2, 
                  max = 1.5, 
                  value = .5,
                  step = .1
                )
              )
            ),
            plotOutput(
              outputId = ns("genome_plot"),
              height = "700px",
              width = "100%"
            )
          ),
          tabPanel(
            "(PREDA) Significant Loci",
            DT::dataTableOutput(outputId = ns("chr_regions"))
          ),
          tabPanel(
            "(PREDA) Genes",
            DT::dataTableOutput(outputId = ns("genes_chr_regions"))
          )
        )
      )
    )   
  )
}
    
#' mod_07_genome Server Functions
#'
#' @noRd 
mod_07_genome_server <- function(id, pre_process, deg, idep_data){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    output$list_comparisons_genome <- renderUI({
      if(is.null(deg$limma()$comparisons)) {
        selectInput(
          inputId = ns("select_contrast"),
          label = NULL,
          choices = list("All" = "All"),
          selected = "All"
        )  
			}	else {
        selectInput(
          inputId = ns("select_contrast"),
          label = 
            "Select a comparison to examine. \"A-B\" means A vs. B (See heatmap).
            Interaction terms start with \"I:\"",
          choices = deg$limma()$comparisons
	     )
      } 
	  })
    
    # visualizing fold change on chrs. 
    output$genome_plotly <- plotly::renderPlotly({
      req(!is.null(deg$limma()))

      chromosome_plotly(
        limma = deg$limma(),
        select_contrast = input$select_contrast,
        all_gene_info = pre_process$all_gene_info(),
        ignore_non_coding = input$ignore_non_coding,
        limma_p_val_viz = input$limma_p_val_viz,
        limma_fc_viz = input$limma_fc_viz,
        label_gene_symbol = input$label_gene_symbol,
        ma_window_size = input$ma_window_size,
        ma_window_steps = input$ma_window_steps,
        ch_region_p_val = input$ch_region_p_val
      )
    })

    # Pre-calculating PREDA, so that changing FDR cutoffs does not trigger entire calculation
    genome_plot_data_pre <- reactive({
      req(!is.null(deg$limma()))

      get_genome_plot_data_pre(
        select_contrast = input$select_contrast,
        limma = deg$limma(),
        all_gene_info = pre_process$all_gene_info()
      ) 
    })

    # Results from PREDA
    genome_plot_data <- reactive({
      req(!is.null(genome_plot_data_pre()))

      get_genome_plot_data(
        genome_plot_data_pre = genome_plot_data_pre(),
        all_gene_info = pre_process$all_gene_info(),
        select_contrast = input$select_contrast,
        limma = deg$limma(),
        regions_p_val_cutoff = input$regions_p_val_cutoff,
        statistic_cutoff = input$statistic_cutoff
      )
    })

    # Using PREDA to identify significant genomic regions 
    output$genome_plot <- renderPlot({
      req(!is.null(genome_plot_data()))

	    get_genome_plot(
        genome_plot_data = genome_plot_data(),
        regions_p_val_cutoff = input$regions_p_val_cutoff,
        statistic_cutoff = input$statistic_cutoff
      )
    })

    output$chr_regions <- DT::renderDataTable({
	    req(!is.null(genome_plot_data()))

      region_data <- genome_plot_data()$Regions[, c(8, 1:6, 9)]
	    colnames(region_data)[1] <- "RegionID"

      DT::datatable(
        region_data,
        options = list(
          pageLength = 20,
          scrollX = "400px"
        ),
        rownames = FALSE
      )
    })

    output$genes_chr_regions <- DT::renderDataTable({
	    req(!is.null(genome_plot_data()))

      genes <- genome_plot_data()$Genes[, -c(5, 10, 12)]
	    genes$Fold <- round(genes$Fold, 3)
	    colnames(genes)[2] <- "Dir"

      DT::datatable(
        genes,
        options = list(
          pageLength = 20,
          scrollX = "400px"
        ),
        rownames = FALSE
      )
    })
  })
}
    
## To be copied in the UI
# mod_07_genome_ui("mod_07_genome_ui_1")
    
## To be copied in the server
# mod_07_genome_server("mod_07_genome_ui_1")
