#' 08_bicluster UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_08_bicluster_ui <- function(id){
  ns <- NS(id)
  tabPanel(
    "Bicluster",
    sidebarLayout(
      sidebarPanel(
        h5(
          "Biclustering can discover genes correlated on subset of samples. 
           Only useful when sample size is large(N>15)
           and more than 2 sample groups. 
           Based on the", 
           a("biclust", 
             href="https://cran.r-project.org/web/packages/biclust/index.html"
           ),
           " and ", 
           a("QUBIC", 
             href="https://www.bioconductor.org/packages/release/bioc/html/QUBIC.html"
           ),
           " R packages."
        ),
        numericInput(
          inputId = ns("n_genes"), 
          label = h5("Most variable genes to include: "), 
          min = 10, 
          max = 2000, 
          value = 1000
        ),
        selectInput(
          inputId = ns("biclust_method"),
          label = "Method:", 
          choices = list( 
            "BCCC" = "biclust::BCCC()",
            "QUBIC" = "QUBIC::BCQU()",
            "runibic" = "runibic::BCUnibic()",
            "BCXmotifs" = "biclust::BCXmotifs()",
            "BCPlaid" = "biclust::BCPlaid()",
            "BCSpectral" = "biclust::BCSpectral()",
            "BCBimax" = "biclust::BCBimax()",
            "BCQuest" = "biclust::BCQuest()"
          ),
          selected = "BCCC()"
        ),
        selectInput(
          inputId = ns("heatmap_color_select"),
          label = "Select Heatmap Color: ",
          choices = "green-black-red",
          width = "100%"
        ),
        htmlOutput(outputId = ns("list_biclusters")),
        textOutput(ns("bicluster_info")),
        a(
          h5(
            "Questions?",
            align = "right"
          ),
          href = "https://idepsite.wordpress.com/biclustering/",
          target = "_blank"
        )
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Heatmap",
            fluidRow(
              column(
                width = 3,
                plotOutput(
                  outputId = ns("biclust_main_heatmap"),
                  height = "450px",
                  width = "100%",
                  brush = ns("ht_brush")
                ),
                br(),
                h5("Selected Cell (Submap):"),
                uiOutput(
                  outputId = ns("ht_click_content")
                )
              ),
              column(
                width = 9,
                plotOutput(
                  outputId = ns("biclust_sub_heatmap"),
                  height = "650px",
                  width = "100%",
                  click = ns("ht_click")
                )
              )
            )
          ),
          tabPanel(
            "Enrichment",
            fluidRow(
              column(
                width = 4,
                style = "margin-top: 15px;",
                htmlOutput(outputId = ns("select_go_selector")),
              ),
              column(
                width = 4,
                checkboxInput(
                  inputId = ns("filtered_background"), 
                  label = "Use filtered data as background in enrichment (slow)", 
                  value = TRUE
                )
              ),
              column(
                width = 4,
                checkboxInput(
                  inputId = ns("remove_redudant"),
                  label = "Remove Redudant Gene Sets",
                  value = FALSE
                )
              )
            ),
            tableOutput(ns("pathway_data_biclust"))

          ),
          tabPanel(
            "Genes",
            DT::dataTableOutput(
              outputId = ns("gene_list_bicluster"),
            ),
            br(),
            uiOutput(ns("download_biclust_button"))
          )
        )
      )
    )
  )
}
    




#------------------------------------------------------------------------



#' 08_bicluster Server Functions
#'
#' @noRd 
mod_08_bicluster_server <- function(id, pre_process, idep_data, tab){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    # Interactive heatmap environment
    biclust_env <- new.env()

    # GMT choices for enrichment ----------
    output$select_go_selector <- renderUI({
	    req(!is.null(pre_process$gmt_choices()))

	    selectInput(
        inputId = ns("select_go"),
        label = NULL,
        choices = pre_process$gmt_choices(),
        selected = "GOBP"
      )
    })

    # all clusters
    biclustering <- reactive({
      req(!is.null(pre_process$data()))

      get_biclustering(
        data = pre_process$data(),
        n_genes = input$n_genes,
        biclust_method = input$biclust_method
      )
    })

    output$list_biclusters <- renderUI({
      req(tab() == "Bicluster")
		  req(!is.null(biclustering()))
      req(biclustering()$res@Number != 0)
		
			selectInput(
        inputId = ns("select_bicluster"), 
				label = "Select a cluster",
        selected = 1,
			  choices = 1:biclustering()$res@Number  
			)		
	  })

    # Heatmap Colors ----------
    heatmap_colors <- list(
      "Green-Black-Red" = c("green", "black", "red"),
      "Blue-White-Red" = c("blue", "white", "red"),
      "Green-Black-Magenta" = c("green", "black", "magenta"),
      "Blue-Yellow-Red" = c("blue", "yellow", "red"),
      "Blue-White-Brown" = c("blue", "white", "brown")
    )
    heatmap_choices <- c(
      "Green-Black-Red",
      "Blue-White-Red",
      "Green-Black-Magenta",
      "Blue-Yellow-Red",
      "Blue-White-Brown"
    )
    observe({
      updateSelectInput(
        session = session,
        inputId = "heatmap_color_select",
        choices = heatmap_choices
      )
    })

    # information for a specific cluster
    biclust_data <- reactive({
      req(!is.null(biclustering()) && !is.null(input$select_bicluster))
      req(biclustering()$res@Number != 0)

      return(
        biclust::bicluster(
          biclustering()$data,
          biclustering()$res,
          as.numeric(input$select_bicluster)
        )[[1]]
      )
    })

    output$biclust_main_heatmap <- renderPlot({
      req(!is.null(biclust_data()))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Creating Heatmap",
        color = "#000000"
      )

      # Assign heatmap to be used in multiple components
      biclust_env$ht <- basic_heatmap(
        data = biclust_data(),
        heatmap_color_select = heatmap_colors[[input$heatmap_color_select]]
      )

      # Use heatmap position in multiple components
      biclust_env$ht_pos_main <- InteractiveComplexHeatmap::htPositionsOnDevice(biclust_env$ht)

      shinybusy::remove_modal_spinner()

      return(biclust_env$ht)
    })

    output$biclust_sub_heatmap <- renderPlot({
      if (is.null(input$ht_brush)) {
        grid::grid.newpage()
        grid::grid.text("No region is selected.", 0.5, 0.5)
      } else {
        biclust_heat_return <- basic_heat_sub(
          ht_brush = input$ht_brush,
          ht = biclust_env$ht,
          ht_pos_main = biclust_env$ht_pos_main,
          heatmap_data = biclust_data()
        )

        biclust_env$ht_select <- biclust_heat_return$ht_select
        biclust_env$submap_data <- biclust_heat_return$submap_data
        biclust_env$group_colors <- biclust_heat_return$group_colors
        biclust_env$column_groups <- biclust_heat_return$column_groups
        
        biclust_env$ht_sub <- ComplexHeatmap::draw(
          biclust_env$ht_select,
          annotation_legend_side = "top",
          heatmap_legend_side = "top"
        )

        biclust_env$ht_pos_sub <- InteractiveComplexHeatmap::htPositionsOnDevice(
          biclust_env$ht_sub
        )

        return(biclust_env$ht_sub)
      }
    })

    # Sub Heatmap Click Value ---------
    output$ht_click_content <- renderUI({
      if (is.null(input$ht_click)) { 
        "Click for Info."
      } else {
        heat_click_info(
          click = input$ht_click,
          ht_sub = biclust_env$ht_sub,
          ht_sub_obj = biclust_env$ht_select,
          ht_pos_sub = biclust_env$ht_pos_sub,
          sub_groups = biclust_env$column_groups,
          group_colors = biclust_env$group_colors,
          data = biclust_env$submap_data
        )
      }
    })

    # Biclustering summary message -----------
    output$bicluster_info <- renderText({		
      req(!is.null(biclustering()) && !is.null(input$select_bicluster))

      bicluster_summary_message(
        biclustering = biclustering(),
        select_bicluster = input$select_bicluster
      )	
	  }) 

    pathway_table_biclust <- reactive({
      req(!is.null(biclust_data()))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Running Analysis",
        color = "#000000"
      )

      gene_names <- merge_data(
        all_gene_names = pre_process$all_gene_names(),
        data = biclust_data(),
        merge_ID = "ensembl_ID"
      )
      # Only keep the gene names and scrap the data
      gene_names_query <- dplyr::select_if(gene_names, is.character)

      req(!is.null(input$select_go))

      gene_sets <- read_pathway_sets(
        all_gene_names_query = gene_names_query,
        converted = pre_process$converted(),
        go = input$select_go,
        select_org = pre_process$select_org(),
        gmt_file = pre_process$gmt_file(),
        idep_data = idep_data,
        gene_info = pre_process$all_gene_info()
      )

      pathway_info <- find_overlap(
        pathway_table = gene_sets$pathway_table,
        query_set = gene_sets$query_set,
        total_genes = gene_sets$total_genes,
        processed_data = pre_process$data(),
        gene_info = pre_process$all_gene_info(),
        go = input$select_go,
        idep_data = idep_data,
        select_org = pre_process$select_org(),
        sub_pathway_files = gene_sets$pathway_files,
        use_filtered_background = input$filtered_background,
        reduced = input$remove_redudant
      )

      shinybusy::remove_modal_spinner()

      return(pathway_info)
    })

    # Enrichment Data Table ----------
    output$pathway_data_biclust <- renderTable({
      req(!is.null(pathway_table_biclust()))

      if(ncol(pathway_table_biclust()) > 1) {
        pathway_table <- pathway_table_biclust()[, 1:4]
      } else {
        pathway_table <- pathway_table_biclust()
      }
      pathway_table
    }, digits = -1, 
     spacing="s", 
     striped=TRUE,
     bordered = TRUE, 
     width = "auto",
     hover=TRUE, 
     sanitize.text.function = function(x) x )

    # list of genes and symbols in the currently selected cluster
    genes_in_selected_cluster <- reactive({
      req(!is.null(biclust_data()))
      req(!is.null(biclustering()) && !is.null(input$select_bicluster) && !is.null(input$select_go))

      get_biclust_table_data(
        res = biclustering()$res,
        biclust_data = biclust_data(),
        select_go = input$select_go,
        select_org = pre_process$select_org(),
        all_gene_info = pre_process$all_gene_info()
      )

    })

    output$gene_list_bicluster <- DT::renderDataTable({
      req(!is.null(biclust_data()) && !is.null(genes_in_selected_cluster()))
      req(!is.null(biclustering()) && !is.null(input$select_bicluster) && !is.null(input$select_go))

      DT::datatable(
        genes_in_selected_cluster(),
        options = list(
          pageLength = 20,
          scrollX = "400px"
          #dom = 'ftg' # hide "Show 20 entries"
        ),
        class = 'cell-border stripe',
        rownames = FALSE
      )
    })

    output$download_biclust <- downloadHandler(
      filename = "bicluster.csv", 
      content = function(file) {
        write.csv(genes_in_selected_cluster(), file, row.names = FALSE)
      }
    )
    
    output$download_biclust_button <- renderUI({
      req(!is.null(biclust_data()) && !is.null(genes_in_selected_cluster()))
      req(!is.null(biclustering()) && !is.null(input$select_bicluster) && !is.null(input$select_go))
      downloadButton(
        outputId = ns("download_biclust"), 
        "All genes"
      )
    })



  })
}
