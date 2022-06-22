#' 09_network UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_09_network_ui <- function(id){
  ns <- NS(id)
  tabPanel(
    "Network",
    sidebarLayout(
      sidebarPanel(
        h5(
          "Identify co-expression networks and sub-modules using",
          a(
            "WGCNA.",
            href = "https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559",
            target = "_blank"
          ),
          "Only useful when  sample size is large(> 15)."
        ),
        numericInput(
          inputId = ns("n_genes_network"), 
          label = h5("Most variable genes to include (< 3001)"), 
          min = 10, 
          max = 3000, 
          value = 1000
        ),
        fluidRow(
          column(
            width = 6,
            numericInput(
              inputId = ns("soft_power"), 
              label = h5("Soft Threshold"), 
              min = 1, 
              max = 20, 
              value = 5
            )
          ),
          column(
            width = 6,
            numericInput(
              inputId = ns("min_module_size"),
              label = h5("Min. Module Size"),
              min = 10, 
              max = 100, 
              value = 20
            )
          )
        ),
        tags$style(
          type = "text/css",
          "#network-soft_power{ width:100%;   margin-top:-12px}"
        ),
        tags$style(
          type = "text/css",
          "#network-min_module_size{ width:100%;   margin-top:-12px}"
        ),
        uiOutput(outputId = ns("heatmap_color_ui")),
        HTML(
          "<hr style='height:1px;border:none;color:#333;background-color:#333;' />"
        ),
        htmlOutput(outputId = ns("list_wgcna_modules")),
        fluidRow(
          column(
            width = 6,
            numericInput(
              inputId = ns("edge_threshold"), 
              label = h5("Edge Threshold"), 
              min = 0, 
              max = 1, 
              value = .4, 
              step = .1
            )
          ),
          column(
            width = 6,
            numericInput(
              inputId = ns("top_genes_network"), 
              label = h5("Top genes"), 
              min = 10, 
              max = 2000, 
              value = 10, 
              step = 10
            )
          )
        ),
        tags$style(
          type = "text/css",
          "#network-edge_threshold{ width:100%;   margin-top:-12px}"
        ),
        tags$style(
          type = "text/css",
          "#network-top_genes_network{ width:100%;   margin-top:-12px}"
        ),
        br(),
        htmlOutput(
          outputId = ns("select_go_selector")
        ),
        h5("The network file can be imported to", 
          a("VisANT", href = "http://visant.bu.edu/", target = "_blank"),
          " or ", 
          a("Cytoscape.", href = "http://www.cytoscape.org/", target = "_blank")
        ),
        # Show Network message
        actionButton(
          inputId = ns("show_messages"),
          label = "Show Network Summary"
        ),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/network/",
          target = "_blank"
        )
      ),
      mainPanel(
        tabsetPanel(
          id = ns("network_tabs"),
          tabPanel(
            "Module Plot",
            plotOutput(
              outputId = ns("module_plot"),
              width = "100%",
              height = "500px"
            )
          ),
          tabPanel(
            "Network Plot",
            br(),
            actionButton(
              inputId = ns("network_layout"),
              label = "Change network layout",
              style = "float:center"
            ),
            plotOutput(outputId = ns("module_network"))
          ),
          tabPanel(
            "Scale Independence",
            br(),
            plotOutput(
              outputId = ns("scale_independence_plot"),
              width = "100%",
              height = "500px"
            )
          ),
          tabPanel(
            "Mean Connectivity",
            br(),
            plotOutput(
              outputId = ns("mean_connectivity_plot"),
              width = "100%",
              height = "500px"
            )
          ),
          tabPanel(
            "Enrichment Table",
            h3("Enriched gene sets in selected module"),
            fluidRow(
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
            DT::dataTableOutput(
              outputId = ns("pathway_data_network")
            )
          ),
          tabPanel(
            "Heatmap",
            fluidRow(
              column(
                width = 3,
                plotOutput(
                  outputId = ns("network_main_heatmap"),
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
                  outputId = ns("network_sub_heatmap"),
                  height = "650px",
                  width = "100%",
                  click = ns("ht_click")
                )
              )
            )
          )
        )
      )
    )
  )
}
    
#' 09_network Server Functions
#'
#' @noRd 
mod_09_network_server <- function(id, pre_process, idep_data, tab){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    # Interactive heatmap environment
    network_env <- new.env()

    output$list_wgcna_modules <- renderUI({
      req(!is.null(wgcna()))
      module_list <- get_wgcna_modules(wgcna = wgcna())
      req(!is.null(module_list))
      selectInput(
        inputId = ns("select_wgcna_module"), 
				label = "Select a module",
				choices = module_list
			)
	  })

    # GMT choices for enrichment ----------
    output$select_go_selector <- renderUI({
	    req(!is.null(pre_process$gmt_choices()))

	    selectInput(
        inputId = ns("select_go"),
        label = "Select Geneset:",
        choices = pre_process$gmt_choices(),
        selected = "GOBP"
      )
    })

    wgcna <- reactive({
      req(!is.null(pre_process$data()))

      get_wgcna(
        data = pre_process$data(),
        n_genes = input$n_genes_network,
        soft_power = input$soft_power,
        min_module_size = input$min_module_size
      )
	  })

    output$module_plot <- renderPlot({
      req(!is.null(wgcna()))

      get_module_plot(wgcna())
		})

    network <- reactiveValues(network_plot = NULL)

    observe({
      req(!is.null(input$select_wgcna_module))
      req(!is.null(wgcna()))

      network$network_plot <- get_network_plot(
        select_wgcna_module = input$select_wgcna_module,
        wgcna = wgcna(),
        top_genes_network = input$top_genes_network,
        select_go = input$select_go,
        select_org = pre_process$select_org(),
        all_gene_info = pre_process$all_gene_info(),
        edge_threshold = input$edge_threshold
      )
    })

    observeEvent(
      input$network_layout, {
        req(!is.null(input$select_wgcna_module))
        req(!is.null(wgcna()))

        network$network_plot <- get_network_plot(
          select_wgcna_module = input$select_wgcna_module,
          wgcna = wgcna(),
          top_genes_network = input$top_genes_network,
          select_go = input$select_go,
          select_org = pre_process$select_org(),
          all_gene_info = pre_process$all_gene_info(),
          edge_threshold = input$edge_threshold
        )
	  })

    output$module_network <- renderPlot({
      network$network_plot()
	  })

    network_query <- reactive({
      req(!is.null(input$select_wgcna_module))
      req(!is.null(wgcna()))

      network_query <- network_enrich_data(
        select_wgcna_module = input$select_wgcna_module,
        wgcna = wgcna()
      )
    })

    pathway_table_network <- reactive({
      req(!is.null(network_query()))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Running Analysis",
        color = "#000000"
      )

      gene_names_query <- dplyr::filter(
        pre_process$all_gene_names(),
        ensembl_ID %in% network_query()
      )

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

    # Pathway Data Table ----------
    output$pathway_data_network <- DT::renderDataTable({
      req(!is.null(pathway_table_network()))

      if(ncol(pathway_table_network()) > 1) {
        pathway_table <- pathway_table_network()[, 1:4]
      } else {
        pathway_table <- pathway_table_network()
      }

      DT::datatable(
        pathway_table,
        options = list(
          pageLength = 20,
          scrollX = "400px"
        ),
        rownames = FALSE
      )
    })

    output$scale_independence_plot <- renderPlot({
      req(!is.null(wgcna()))

      plot_scale_independence(
        wgcna = wgcna()
      )
	  })

    output$mean_connectivity_plot <- renderPlot({
      req(!is.null(wgcna()))
      
      plot_mean_connectivity(
        wgcna = wgcna()
      )
	  })

    module_statistic <- reactive({
      req(!is.null(wgcna()))
      paste(
        "A network of",
        wgcna()$n_genes,
        "genes was divided into ",
        wgcna()$n_modules,
        "modules."
      )
    })

    # Show messages when on the Network tab or button is clicked
    observe({
      req(input$show_messages || tab() == "Network")
      req(!is.null(module_statistic()))

      showNotification(
        ui = module_statistic(),
        id = "network_summary",
        duration = NULL,
        type = "default"
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
    output$heatmap_color_ui <- renderUI({
      req(!is.null(input$network_tabs == "Heatmap"))

      selectInput(
        inputId = ns("heatmap_color_select"),
        label = "Select Heatmap Color: ",
        choices = heatmap_choices,
        width = "100%"
      )
    })

    network_data <- reactive({
      req(!is.null(network_query()))

      data <- pre_process$data()[rownames(pre_process$data()) %in% network_query(), ]

      if(ncol(pre_process$all_gene_names()) > 2) {
        data <- rowname_id_swap(
          data_matrix = data,
          all_gene_names = pre_process$all_gene_names(),
          select_gene_id = "symbol"
        )
      }
    })

    output$network_main_heatmap <- renderPlot({
      req(!is.null(network_data()))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Creating Heatmap",
        color = "#000000"
      )

      # Assign heatmap to be used in multiple components
      network_env$ht <- basic_heatmap(
        data = network_data(),
        heatmap_color_select = heatmap_colors[[input$heatmap_color_select]]
      )

      # Use heatmap position in multiple components
      network_env$ht_pos_main <- InteractiveComplexHeatmap::htPositionsOnDevice(network_env$ht)

      shinybusy::remove_modal_spinner()

      return(network_env$ht)
    })

    output$network_sub_heatmap <- renderPlot({
      if (is.null(input$ht_brush)) {
        grid::grid.newpage()
        grid::grid.text("No region is selected.", 0.5, 0.5)
      } else {
        network_heat_return <- basic_heat_sub(
          ht_brush = input$ht_brush,
          ht = network_env$ht,
          ht_pos_main = network_env$ht_pos_main,
          heatmap_data = network_data()
        )

        network_env$ht_select <- network_heat_return$ht_select
        network_env$submap_data <- network_heat_return$submap_data
        network_env$group_colors <- network_heat_return$group_colors
        network_env$column_groups <- network_heat_return$column_groups
        
        network_env$ht_sub <- ComplexHeatmap::draw(
          network_env$ht_select,
          annotation_legend_side = "top",
          heatmap_legend_side = "top"
        )

        network_env$ht_pos_sub <- InteractiveComplexHeatmap::htPositionsOnDevice(
          network_env$ht_sub
        )

        return(network_env$ht_sub)
      }
    })

    # Sub Heatmap Click Value ---------
    output$ht_click_content <- renderUI({
      if (is.null(input$ht_click)) { 
        "Click for Info."
      } else {
        heat_click_info(
          click = input$ht_click,
          ht_sub = network_env$ht_sub,
          ht_sub_obj = network_env$ht_select,
          ht_pos_sub = network_env$ht_pos_sub,
          sub_groups = network_env$column_groups,
          group_colors = network_env$group_colors,
          data = network_env$submap_data
        )
      }
    })

  })
}
    
## To be copied in the UI
# mod_09_network_ui("09_network_ui_1")
    
## To be copied in the server
# mod_09_network_server("09_network_ui_1")
