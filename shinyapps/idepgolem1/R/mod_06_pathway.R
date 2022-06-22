#' 06_pathway UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_06_pathway_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Pathway",
    sidebarLayout(
      sidebarPanel(
        htmlOutput(
          outputId = ns("list_comparisons_pathway")
        ),
        tags$style(
          type = "text/css",
          "#pathway-list_comparisons_pathway { width:100%;   margin-top:-12px}"
        ),

        selectInput(
          inputId = ns("pathway_method"), 
          label = "Select method:", 
          choices = list(
            "GAGE" = 1, 
            "GSEA (preranked fgsea)" = 3,
            "PGSEA" = 2, 
            "PGSEA w/ all samples" = 4, 
            "ReactomePA" = 5
          ),
          selected = 1
        ),
        tags$style(
          type = "text/css",
          "#pathway-pathway_method { width:100%;   margin-top:-12px}"
        ),
        htmlOutput(
          outputId = ns("select_go_selector")
        ),
        tags$style(
          type = "text/css",
          "#pathway-select_go { width:100%;   margin-top:-12px}"
        ),
        fluidRow( 
          column(
            width = 6,
            numericInput(
              inputId = ns("min_set_size"), 
              label = h5("Geneset size: Min."), 
              min = 5, 
              max = 30, 
              value = 15,
              step = 1
            )
          ),
          column(
            width = 6,
            numericInput(
              inputId = ns("max_set_size"), 
              label = h5("Max."), 
              min = 1000, 
              max = 2000, 
              value = 2000,
              step = 100
            )
          )
        ),
        tags$style(
          type = "text/css",
          "#pathway-min_set_size { width:100%;   margin-top:-12px}"
        ),
        tags$style(
          type = "text/css",
          "#pathway-max_set_size { width:100%;   margin-top:-12px}"
        ),
        numericInput(
          inputId = ns("pathway_p_val_cutoff"), 
          label = h5("Pathway signifiance cutoff (FDR)"),
          value = 0.2,
          min = 1e-20,
          max = 1,
          step = .05
        ),
        tags$style(
          type = "text/css",
          "#pathway-pathway_p_val_cutoff { width:100%;   margin-top:-12px}"
        ),
        numericInput(
          inputId = ns("n_pathway_show"), 
          label = h5("Number of top pathways to show"),
          value = 30, 
          min = 5,
          max = 100,
          step = 5
        ),
        tags$style(
          type = "text/css",
          "#pathway-n_pathway_show { width:100%;   margin-top:-12px}"
        ),
        checkboxInput(
          inputId = ns("absolute_fold"),
          label = "Use absolute values of fold changes for GSEA and GAGE",
          value = FALSE
        ),
        numericInput(
          inputId = ns("gene_p_val_cutoff"), 
          label = h5("Remove genes with big FDR before pathway analysis:"),
          value = 1,
          min = 1e-20,
          max = 1,
          step = .05
        ),
        tags$style(
          type = "text/css",
          "#pathway-gene_p_val_cutoff { width:100%;   margin-top:-12px}"
        ),
        h5("* Warning! The many combinations can lead to false positives in pathway analyses."),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/pathways/",
          target = "_blank"
        )    
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Results",
            br(),
            conditionalPanel(
              condition = "input.pathway_method == 1",
              DT::dataTableOutput(outputId = ns("gage_pathway_table")),
              ns = ns
            ),
            conditionalPanel(
              condition = "input.pathway_method == 2",
              h5("Red and blue indicates activated and suppressed pathways, respectively."),
              plotOutput(
                outputId = ns("pgsea_plot"),
                inline = TRUE
              ),
              ns = ns
            ),
            conditionalPanel(
              condition = "input.pathway_method == 3",
              DT::dataTableOutput(outputId = ns("fgsea_pathway")),
              ns = ns
            ),
            conditionalPanel(
              condition = "input.pathway_method == 4",
              h5("Red and blue indicates activated and suppressed pathways, respectively."),
              plotOutput(
                outputId = ns("pgsea_plot_all_samples"),
                inline = TRUE
              ),
              ns = ns
            ),
            conditionalPanel(
              condition = "input.pathway_method == 5",
              DT::dataTableOutput(outputId = ns("reactome_pa_pathway")),
              ns = ns
            )
          ),
          tabPanel(
            "Heatmap",
            conditionalPanel(
              condition = "input.pathway_method == 1 | input.pathway_method == 2 |
                           input.pathway_method == 3 | input.pathway_method == 4",
              fluidRow(
                br(),
                column(
                  width = 5,
                  htmlOutput(outputId = ns("list_sig_pathways")),
                ),
                column(
                  width = 4, 
                  conditionalPanel(
                    condition =  "input.select_go != 'KEGG'", 
                    selectInput(
                      inputId = ns("heatmap_color_select"),
                      label = "Select Heatmap Color: ",
                      choices = "green-black-red",
                      width = "100%"
                    ), 
                    ns = ns
                  ), 
                  conditionalPanel(
                    condition = "input.select_go == 'KEGG'", 
                    selectInput(
                      inputId = ns("kegg_color_select"), 
                      label = "Select colors (low-high)", 
                      choices = "green-red", 
                      width = "100%"
                    ), 
                    ns = ns
                  )
                )
              ),
              ns = ns
            ),
            conditionalPanel(
              condition = "(input.pathway_method == 1 | input.pathway_method == 2 |
                            input.pathway_method == 3 | input.pathway_method == 4) &
                            input.select_go != 'KEGG'",
              h5("Brush for sub-heatmap, click for value. (Shown Below)"),
              br(),
              fluidRow(
                column(
                  width = 3,
                  plotOutput(
                    outputId = ns("path_main_heatmap"),
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
                    outputId = ns("path_sub_heatmap"),
                    height = "650px",
                    width = "100%",
                    click = ns("ht_click")
                  )
                )
              ),
              ns = ns
            ),
            conditionalPanel(
              condition = "(input.pathway_method == 1 | input.pathway_method == 2 | 
                            input.pathway_method == 3 | input.pathway_method == 4) &
                            input.select_go == 'KEGG'",
              h5("Red and green represent up- and down-regulated genes, respectively."),
              imageOutput(
                outputId = ns("kegg_image"),
                width = "100%",
                height = "100%"
              ),
              ns = ns
            )
          ),
          tabPanel(
            "Tree",
            plotOutput(
              outputId = ns("enrichment_tree"),
              width = "100%"
            )
          ),
          tabPanel(
            "Network",
            h5("Connected gene sets share more genes. Color of node correspond to adjuested Pvalues."),
            fluidRow(
              column(
                width = 2,
                actionButton(
                  inputId = ns("layout_vis_deg"),
                  label = "Change layout"
                )
              ),
              column(
                width = 1,
                h5("Cutoff:"),
                align="right"
              ),
              column(
                width = 2,
                numericInput(
                  inputId = ns("edge_cutoff_deg"),
                  label = NULL,
                  value = 0.30,
                  min = 0,
                  max = 1,
                  step = .1
                ),
                align="left"
              ),
              column(
                width = 2,
                checkboxInput(
                  inputId = ns("wrap_text_network_deg"),
                  label = "Wrap text",
                  value = TRUE
                )
              )
            ),
            selectInput(
              inputId = ns("up_down_reg_deg"),
              NULL,
              choices = c(
                "Both Up & Down" = "Both",
                "Up regulated" = "Up",
                "Down regulated" = "Down"
              )
            ),
            h6(
              "Two pathways (nodes) are connected if they share 30% (default, adjustable) or more genes.
              Green and red represents down- and up-regulated pathways. You can move the nodes by 
              dragging them, zoom in and out by scrolling, and shift the entire network by click on an 
              empty point and drag. Darker nodes are more significantly enriched gene sets. Bigger nodes
              represent larger gene sets. Thicker edges represent more overlapped genes."
            ),
            visNetwork::visNetworkOutput(
              outputId = ns("vis_network_path"),
              height = "800px",
              width = "100%"
            )
          )
        )
      )
    )
  )
}

#' mod_06_pathway Server Functions
#'
#' @noRd
mod_06_pathway_server <- function(id, pre_process, deg, idep_data, tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Interactive heatmap environment
    path_env <- new.env()

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

    output$list_comparisons_pathway <- renderUI({
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

    output$list_sig_pathways <- renderUI({
	    if(tab() != "Pathway") {
        selectInput(
          inputId = ns("sig_pathways"),
          label = NULL, 
          choices = list("All" = "All"),
          selected = "All"
        )
      }	else {
        req(!is.null(input$pathway_method))
        # Default, sometimes these methods returns "No significant pathway found"
		    choices <- "All"  
		    if(input$pathway_method == 1) { 
			    if(!is.null(gage_pathway_data())) {
            if(dim(gage_pathway_data())[2] > 1) {
              choices <- gage_pathway_data()[, 2]
            } 
          }
		    } else if(input$pathway_method == 2) {
          if(!is.null(pgsea_plot_data())) {
            if(dim(pgsea_plot_data())[2] > 1) {
					 	  pathways <- as.data.frame(pgsea_plot_data())
						  choices <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
					  }
          }
				} else if(input$pathway_method == 3) {
          if(!is.null(fgsea_pathway_data())) {
            if(dim(fgsea_pathway_data())[2] > 1) {
              choices <- fgsea_pathway_data()[, 2]
            }
          }
        } else if(input$pathway_method == 4) {
          if(!is.null(pgsea_plot_all_samples_data())) {
            if(dim(pgsea_plot_all_samples_data())[2] > 1) {
					    pathways <- as.data.frame(pgsea_plot_all_samples_data())
					    choices <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
				    }
          }
        } else if(input$pathway_method == 5) {
			    if(!is.null(reactome_pa_pathway_data())) {
            if(dim(reactome_pa_pathway_data())[2] > 1) {
              choices <- reactome_pa_pathway_data()[, 2]
            }  
          }
        }
        
        selectInput(
          inputId = ns("sig_pathways"),
          label = 
            "Select a pathway to show expression pattern of related genes on a heatmap or a
             KEGG pathway diagram:",
          choices = choices
        )
	    } 
	  })

    gene_sets <- reactive({
      req(tab() == "Pathway")
      req(!is.null(input$select_go))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Finding Gene Sets",
        color = "#000000"
      )

      if(pre_process$select_org() == "NEW" && !is.null(pre_process$gmt_file())) {
        in_file <- pre_process$gmt_file()
        in_file <- in_file$datapath
				gene_sets <- read_gmt_robust(in_file)
      } else {
        gene_sets <- read_gene_sets(
          converted = pre_process$converted(),
          all_gene_names = pre_process$all_gene_names(),
          go = input$select_go,
          select_org = pre_process$select_org(),
          idep_data = idep_data,
          my_range = c(input$min_set_size, input$max_set_size)
        )
      }

      shinybusy::remove_modal_spinner()

      return(gene_sets)
    })

    gage_pathway_data <- reactive({  
      req(input$pathway_method == 1)
      req(!is.null(deg$limma()))
      req(!is.null(gene_sets()))

      gage_data(
        select_go = input$select_go,
        select_contrast = input$select_contrast,
        min_set_size = input$min_set_size,
        max_set_size = input$max_set_size,
        limma = deg$limma(),
        gene_p_val_cutoff = input$gene_p_val_cutoff,
        gene_sets = gene_sets(),
        absolute_fold = input$absolute_fold,
        pathway_p_val_cutoff = input$pathway_p_val_cutoff,
        n_pathway_show = input$n_pathway_show
      )
    })

    output$gage_pathway_table <- DT::renderDataTable({
      req(!is.null(gage_pathway_data()))
      
      DT::datatable(
        gage_pathway_data(),
        options = list(
          pageLength = 15,
          scrollX = "500px"
        ),
        rownames = FALSE
      )
    })

    contrast_samples <- reactive({
      req(!is.null(input$select_contrast))
      req(!is.null(pre_process$data()))

      find_contrast_samples(
        select_contrast = input$select_contrast, 
		    all_sample_names = colnames(pre_process$data()),
		    sample_info = pre_process$sample_info(),
		    select_factors_model = deg$select_factors_model(),
		    select_model_comprions = deg$select_model_comprions(), 
		    reference_levels = deg$reference_levels(),
		    counts_deg_method = deg$counts_deg_method(),
		    data_file_format = pre_process$data_file_format()
	    )
    })

    output$pgsea_plot <- renderPlot({
      req(input$pathway_method == 2)

      plot_pgsea(
        my_range = c(input$min_set_size, input$max_set_size),
        processed_data = pre_process$data(),
        contrast_samples = contrast_samples(),
        gene_sets = gene_sets(),
        pathway_p_val_cutoff = input$pathway_p_val_cutoff,
        n_pathway_show = input$n_pathway_show
      )
    }, 
      height = 800,
      width = 800
    )

    pgsea_plot_data <- reactive({
      req(!is.null(gene_sets()))

      get_pgsea_plot_data(
        my_range = c(input$min_set_size, input$max_set_size),
        data = pre_process$data(),
        select_contrast = input$select_contrast,
        gene_sets = gene_sets(),
        sample_info = pre_process$sample_info(),
        select_factors_model = deg$select_factors_model(),
        select_model_comprions = deg$select_model_comprions(),
        pathway_p_val_cutoff = input$pathway_p_val_cutoff,
        n_pathway_show = input$n_pathway_show
      )
    })

    fgsea_pathway_data <- reactive({
      req(input$pathway_method == 3)
      req(!is.null(deg$limma()))
      req(!is.null(gene_sets()))

      fgsea_data(
        select_contrast = input$select_contrast,
        my_range = c(input$min_set_size, input$max_set_size),
        limma = deg$limma(),
        gene_p_val_cutoff = input$gene_p_val_cutoff,
        gene_sets = gene_sets(),
        absolute_fold = input$absolute_fold,
        pathway_p_val_cutoff = input$pathway_p_val_cutoff,
        n_pathway_show = input$n_pathway_show
      )
    })

    output$fgsea_pathway <- DT::renderDataTable({
      req(!is.null(fgsea_pathway_data()))
      
      DT::datatable(
        fgsea_pathway_data(),
        options = list(
          pageLength = 15,
          scrollX = "500px"
        ),
        rownames = FALSE
      )
	  })

    reactome_pa_pathway_data <- reactive({
      req(!is.null(deg$limma()))

      reactome_data(
        select_contrast = input$select_contrast,
        my_range = c(input$min_set_size, input$max_set_size),
        limma = deg$limma(),
        gene_p_val_cutoff = input$gene_p_val_cutoff,
        converted = pre_process$converted(),
        idep_data = idep_data,
        pathway_p_val_cutoff = input$pathway_p_val_cutoff,
        n_pathway_show = input$n_pathway_show,
        absolute_fold = input$absolute_fold
      )
    })

    output$pgsea_plot_all_samples <- renderPlot({
      req(input$pathway_method == 4)

      pgsea_plot_all(
        go = input$select_go,
        my_range = c(input$min_set_size, input$max_set_size),
        data = pre_process$data(),
        select_contrast = input$select_contrast,
        gene_sets = gene_sets(),
        pathway_p_val_cutoff = input$pathway_p_val_cutoff,
        n_pathway_show = input$n_pathway_show
      )
    }, height = 800, width = 800)

    pgsea_plot_all_samples_data <- reactive({
      req(!is.null(gene_sets()))

      get_pgsea_plot_all_samples_data(
        data = pre_process$data(),
        select_contrast = input$select_contrast,
        gene_sets = gene_sets(),
        my_range = c(input$min_set_size, input$max_set_size),
        pathway_p_val_cutoff = input$pathway_p_val_cutoff,
        n_pathway_show = input$n_pathway_show
      )
    })

    output$reactome_pa_pathway <- DT::renderDataTable({
      req(!is.null(reactome_pa_pathway_data()))
	    
      DT::datatable(
        reactome_pa_pathway_data(),
        options = list(
          pageLength = 15,
          scrollX = "500px"
        ),
        rownames = FALSE
      )
	  })

    selected_pathway_data <- reactive({
      req(!is.null(input$sig_pathways))
      req(!is.null(gene_sets()))
      
      pathway_select_data(
        sig_pathways = input$sig_pathways,
        gene_sets = gene_sets(),
        contrast_samples = contrast_samples(),
        data = pre_process$data(),
        select_org = pre_process$select_org(),
        all_gene_names = pre_process$all_gene_names()
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
    
    # Kegg Colors --------
    kegg_colors <- list(
      "Green-Red" = c("green", "red"), 
      "Blue-Orange" = c("blue", "orange")
    )
    kegg_choices <- c(
      "Green-Red", 
      "Blue-Orange"
    )
    observe({
      updateSelectInput(
        session = session, 
        inputId = "kegg_color_select", 
        choices = kegg_choices
      )
    })

    output$path_main_heatmap <- renderPlot({
      req(!is.null(selected_pathway_data()))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Creating Heatmap",
        color = "#000000"
      )

      # Assign heatmap to be used in multiple components
      path_env$ht <- basic_heatmap(
        data = selected_pathway_data(),
        heatmap_color_select = heatmap_colors[[input$heatmap_color_select]]
      )

      # Use heatmap position in multiple components
      path_env$ht_pos_main <- InteractiveComplexHeatmap::htPositionsOnDevice(path_env$ht)

      shinybusy::remove_modal_spinner()

      return(path_env$ht)
    })

    output$path_sub_heatmap <- renderPlot({
      if (is.null(input$ht_brush)) {
        grid::grid.newpage()
        grid::grid.text("No region is selected.", 0.5, 0.5)
      } else {
        path_heat_return <- basic_heat_sub(
          ht_brush = input$ht_brush,
          ht = path_env$ht,
          ht_pos_main = path_env$ht_pos_main,
          heatmap_data = selected_pathway_data()
        )

        path_env$ht_select <- path_heat_return$ht_select
        path_env$submap_data <- path_heat_return$submap_data
        path_env$group_colors <- path_heat_return$group_colors
        path_env$column_groups <- path_heat_return$column_groups
        
        path_env$ht_sub <- ComplexHeatmap::draw(
          path_env$ht_select,
          annotation_legend_side = "top",
          heatmap_legend_side = "top"
        )

        path_env$ht_pos_sub <- InteractiveComplexHeatmap::htPositionsOnDevice(path_env$ht_sub)

        return(path_env$ht_sub)
      }
    })

    # Sub Heatmap Click Value ---------
    output$ht_click_content <- renderUI({
      if (is.null(input$ht_click)) { 
        "Click for Info."
      } else {
        heat_click_info(
          click = input$ht_click,
          ht_sub = path_env$ht_sub,
          ht_sub_obj = path_env$ht_select,
          ht_pos_sub = path_env$ht_pos_sub,
          sub_groups = path_env$column_groups,
          group_colors = path_env$group_colors,
          data = path_env$submap_data
        )
      }
    })

    output$kegg_image <- renderImage({

      kegg_pathway(
        go = input$select_go,
        gage_pathway_data = gage_pathway_data(),
        sig_pathways = input$sig_pathways,
        select_contrast = input$select_contrast,
        limma = deg$limma(),
        converted = pre_process$converted(),
        idep_data = idep_data,
        select_org = pre_process$select_org(), 
        low_color = kegg_colors[[input$kegg_color_select]][1], 
        high_color = kegg_colors[[input$kegg_color_select]][2]
      )
    }, deleteFile = TRUE)
    
    # List of pathways with details
    pathway_list_data <- reactive({
      get_pathway_list_data(
        pathway_method = input$pathway_method,
        gage_pathway_data = gage_pathway_data(),
        fgsea_pathway_data = fgsea_pathway_data(),
        pgsea_plot_data = pgsea_plot_data(),
        pgsea_plot_all_samples_data = pgsea_plot_all_samples_data(),
        go = input$select_go,
        select_org = pre_process$select_org(),
        gene_info = pre_process$all_gene_info(),
        gene_sets = gene_sets()
      )
    })

    # Enrichment Tree -----------
    output$enrichment_tree <- renderPlot({
      req(!is.null(pathway_list_data()))

      enrichment_plot(
        go_table = pathway_list_data(),
        45
      )
    })

    # Define a Network
    network_data_path <- reactive({
      req(!is.null(pathway_list_data()))

      network_data(
        network = pathway_list_data(),
        up_down_reg_deg = input$up_down_reg_deg,
        wrap_text_network_deg= input$wrap_text_network_deg,
        layout_vis_deg = input$layout_vis_deg,
        edge_cutoff_deg = input$edge_cutoff_deg
      )
    })

    # Interactive vis network plot
    output$vis_network_path <- visNetwork::renderVisNetwork({
      req(!is.null(network_data_path()))
      
      vis_network_plot(
        network_data = network_data_path()
      )
    })
  })
}

## To be copied in the UI
# mod_08_pathway_ui("08_pathway_ui_1")

## To be copied in the server
# mod_08_pathway_server("08_pathway_ui_1")
