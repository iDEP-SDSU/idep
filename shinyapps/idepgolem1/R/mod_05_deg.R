#' 05_deg1 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_05_deg_1_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "DEG1",
    sidebarLayout(
      sidebarPanel(
        # DEG analysis methods for read counts data
        conditionalPanel(
          condition = "output.data_file_format == 1",
          selectInput(
            inputId = ns("counts_deg_method"),
            label = "Method:", 
            choices = list(
              "DESeq2" = 3,
              "limma-voom" = 2,
              "limma-trend" = 1
            ),
            selected = 3 
            
          ),
          tags$style(
            type = 'text/css',
            "#deg-counts_deg_method {width:100%;   margin-top:-12px}"
          ), 
          ns = ns
        ),
        # Label when the limma method is selected
        conditionalPanel(
          condition = "output.data_file_format == 2",
          h5("Using the limma package"), 
          ns = ns
        ),
        fluidRow(
          column(
            width = 5,
            # Adjusted significant p-value to use
            numericInput(
              inputId = ns("limma_p_val"),
              label = h5("FDR cutoff"), 
              value = 0.1,
              min = 1e-5,
              max = 1,
              step = .05
            )
          ),
          column(
            width = 7,
            # Min fold change to use
            numericInput(
              inputId = ns("limma_fc"),
              label = h5("Min fold change"),
              value = 2,
              min = 1,
              max = 100,
              step = 0.5
            )
          ),
          # Style both numeric inputs
          tags$style(
            type = "text/css",
            "#deg-limma_p_val { width:100%;   margin-top:-12px}"
          ),
          tags$style(
            type = "text/css",
            "#deg-limma_fc { width:100%;   margin-top:-12px}"
          )
        ),
        # Button to run DEG analysis for the specified model
        actionButton(
          inputId = ns("submit_model_button"),
          label = "Submit & Calculate",
          style = "float:center"
        ),
        tags$head(tags$style(
          "#deg-submit_model_button{font-size: 20px;}"
        )),
        tags$br(),
        tags$br(),
        uiOutput(ns("download_lfc_button")),
        uiOutput(ns("note")),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/degs/",
          target = "_blank"
        )
      ),
      mainPanel(
        tabsetPanel(
          id = ns("step_1"),
          tabPanel(
            title = "Experiment Design",
            fluidRow(
              column(
                width = 6,
                htmlOutput(outputId = ns("list_factors_deg"))
              ),
              column(
                width = 6,
                htmlOutput(outputId = ns("list_block_factors_deg"))
              ) 
            ),
            fluidRow(
              htmlOutput(outputId = ns("select_reference_levels"))
            ),
            htmlOutput(outputId = ns("list_interaction_terms")),
            textOutput(outputId = ns("experiment_design")),
            tags$head(tags$style(
              "#deg-experiment_design{color: red;font-size: 16px;}"
            )),
            htmlOutput(outputId = ns("list_model_comparisons")),
            h3("Use the submit button in the sidebar once the desired 
               design is selected!"
            ),
            a(
              h5("More info on DESeq2 experiment design", align = "right"),
              href = "http://rpubs.com/ge600/deseq2",
              target = "_blank"
            )
          ),
          tabPanel(
            title = "Results",
            value = ("results_tab"),
            plotOutput(
              outputId = ns("sig_gene_stats")
            ),
            br(),
            br(),
            h4(
              "Numbers of differentially expressed genes for all comparisons.
              \"B-A\" means B vs. A. Interaction terms start with \"I:\" "
            ),
            tableOutput(
              outputId = ns("sig_gene_stats_table")
            )
          ),
          tabPanel(
            title = "Venn Diagram",
            checkboxInput(
              inputId = ns("up_down_regulated"),
              label = "Split gene lists by up- or down-regulation",
              value = TRUE
            ),
            htmlOutput(outputId = ns("list_comparisons_venn")),
            plotOutput(outputId = ns("venn_plot")), 
            ottoPlots::mod_download_figure_ui(
              id = ns("dl_venn"), 
              label = "Dowload venn diagram"
            )
          )
        )
      )
    )
  )
}

mod_05_deg_2_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "DEG2",
    sidebarLayout(
      sidebarPanel(
        h5("Examine the results of DEGs for each comparison"),
        htmlOutput(outputId = ns("list_comparisons")),
        # Heatmap customizing features ----------
        conditionalPanel(
          condition = "input.step_2 == 'Heatmap'",
          fluidRow(
            column(width = 3, h5("Color:")),
            column(
              width = 9,
              selectInput(
                inputId = ns("heatmap_color_select"),
                label = NULL,
                choices = "green-black-red",
                width = "100%"
              )
            )
          ),
          ns = ns
        ),
        conditionalPanel(
          condition = "input.step_2 == 'Volcano Plot' | 
            input.step_2 == 'MA Plot'", 
          fluidRow(
            column(width = 3, h5("Plot colors:")), 
            column(
              width = 9, 
              selectInput(
                inputId = ns("plot_color_select"), 
                label = NULL, 
                choices = "Red-Green"
              )
            )
          ), 
          ns = ns
        ),
        HTML(
          "<hr style='height:1px;border:none;color:
          #333;background-color:#333;' />"
        ),
        h5("Enrichment analysis for DEGs:"),
        htmlOutput(outputId = ns("select_go_selector")),
        tags$style(
          type = "text/css",
          "#deg-select_go { width:100%; margin-top:-9px}"
        ),
        checkboxInput(
          inputId = ns("filtered_background"), 
          label = "Use filtered data as background in enrichment (slow)", 
          value = TRUE
        ),
        checkboxInput(
          inputId = ns("remove_redudant"),
          label = "Remove Redudant Gene Sets",
          value = FALSE
        ), 
        conditionalPanel(
          condition = "input.step_2 == 'Enrich Table'", 
          downloadButton(
            outputId = ns("dl_enrich_up"), 
            label = "Up Enrichment"
          ),
          downloadButton(
            outputId = ns("dl_enrich_down"), 
            label = "Down Enrichment"
          ),
          ns = ns
        )
      ),
      mainPanel(
        tabsetPanel(
          id = ns("step_2"),
          tabPanel(
            title = "Heatmap",
            h5("Brush for sub-heatmap, click for value. (Shown Below)"),
            br(),
            fluidRow(
              column(
                width = 3,
                plotOutput(
                  outputId = ns("deg_main_heatmap"),
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
                  outputId = ns("deg_sub_heatmap"),
                  height = "650px",
                  width = "100%",
                  click = ns("ht_click")
                )
              )
            )
          ),
          tabPanel(
            title = "Volcano Plot",
            br(),
            plotOutput(
              outputId = ns("volcano_plot"),
              height = "500px",
              width = "100%"
            ),
            ottoPlots::mod_download_figure_ui(ns("download_volcano"))
          ),
          tabPanel(
            title = "MA Plot",
            br(),
            plotOutput(
              outputId = ns("ma_plot"),
              height = "500px",
              width = "100%"
            ), 
            ottoPlots::mod_download_figure_ui(ns("download_ma"))
          ),
          
          tabPanel(
            title = "Scatter Plot",
            br(),
            plotOutput(
              outputId = ns("scatter_plot"),
              height = "500px",
              width = "100%"
            )
          ),
          tabPanel(
            title = "Enrich Table",
            br(),
            tags$p("To see the list of genes that were differently expressed in
                    each category, download the data in the left-hand panel."),
            br(),
            strong(h3("Up Regulated Genes")),
            br(),
            DT::dataTableOutput(
              outputId = ns("pathway_data_up")
            ),
            br(),
            strong(h3("Down Regulated Genes")),
            br(),
            DT::dataTableOutput(
              outputId = ns("pathway_data_down")
            )
          ),
          tabPanel(
            title = "Enrich Tree",
            plotOutput(
              outputId = ns("enrichment_tree"),
              width = "100%"
            )
          ),
          tabPanel(
            title = "Pathway Network",
            h5("Connected gene sets share more genes. 
               Color of node correspond to adjuested Pvalues."),
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
              outputId = ns("vis_network_deg"),
              height = "800px",
              width = "100%"
            )
          )
        )
      )        
    )
  )
}

#' 05_deg1 Server Functions
#'
#' @noRd
mod_05_deg_server <- function(id, pre_process, idep_data, load_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Interactive heatmap environment
    deg_env <- new.env()
    
    # DEG STEP 1 ----------
    output$data_file_format <- reactive({
      pre_process$data_file_format()
    })
    outputOptions(
      x = output,
      name = "data_file_format",
      suspendWhenHidden = FALSE
    )

    # Experiment Design UI Elements ------------
    output$list_factors_deg <- renderUI({
      list_factors <- list_factors_ui(
        sample_info = pre_process$sample_info(),
        data_file_format = pre_process$data_file_format(),
        counts_deg_method = input$counts_deg_method
      )
      if(class(list_factors)[1] == "list") {
        return(
          checkboxGroupInput(
            inputId = ns("select_factors_model"), 
            h5(list_factors$title), 
            choices = list_factors$choices,
            selected = NULL
          )
        )
      } else {
        return(list_factors)
      }
	  })

    output$list_block_factors_deg <- renderUI({ 
      choices <- list_block_factors_ui(
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model
      )
      req(!is.null(choices))
      return(
        checkboxGroupInput(
          inputId = ns("select_block_factors_model"), 
          h5("Select a factor for batch effect or paired samples, if needed."), 
          choices = choices,
          selected = NULL
        )
      )
	  })

    output$list_model_comparisons <- renderUI({
      req(pre_process$data())
		  model_comparisons <- list_model_comparisons_ui(
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model,
        processed_data = pre_process$data()
      )
      req(!is.null(model_comparisons))
      checkboxGroupInput(
        inputId = ns("select_model_comprions"), 
			  label = h5(model_comparisons$title),
        choices = model_comparisons$choices,
        selected = model_comparisons$choices[[1]]
      )
	  })

    output$list_interaction_terms <- renderUI({
      interactions <- list_interaction_terms_ui(
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model
      )
      req(!is.null(interactions))
      checkboxGroupInput(
        inputId = ns("select_interactions"), 
				label = h5(
          "Interaction terms between factors(e.g. genotypes repond differently
          to treatment?):"
        ),
				choices = interactions,
        selected = NULL
      )
	  }) 

	
	  # Set limits for selections of factors
    observe({
    	if(length(input$select_factors_model) > 6) {
        updateCheckboxGroupInput(
          session,
          ns("select_factors_model"),
          selected = tail(input$select_factors_model, 6)
        )
      }
      if(input$counts_deg_method !=3 ) {
        if(length(input$select_factors_model) > 2) {
          updateCheckboxGroupInput(
            session,
            ns("select_factors_model"),
            selected = tail(input$select_factors_model, 2)
          )
        }
        if(length(input$select_block_factors_model) > 1) {
          updateCheckboxGroupInput(
            session,
            ns("select_block_factors_model"),
            selected = tail(input$select_block_factors_model, 1)
          )
        }
      }
      
      if(length(input$select_comparisons_venn) > 5) {
        updateCheckboxGroupInput(
          session,
          ns("select_comparisons_venn"),
          selected = tail(input$select_comparisons_venn, 5)
        )
      }
		})

    output$experiment_design <- renderText({
      experiment_design_txt(
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model,
        select_block_factors_model = input$select_block_factors_model,
        select_interactions = input$select_interactions
      )
    })

    output$select_reference_levels <- renderUI({
      select_choices <- select_reference_levels_ui(
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model,
        data_file_format = pre_process$data_file_format(),
        counts_deg_method = input$counts_deg_method
      )
      req(!is.null(select_choices))
      lapply(names(select_choices), function(x) {
        tagList(
          column(
            width = 4,
            selectInput(
              inputId = ns(
                paste0(
                  "reference_level_factor_",
                  which(names(select_choices) == x)
                )
              ), 
							label = h5(paste0("Reference/baseline level for ", x)),
							choices = setNames(
                as.list(paste0(x, ":", select_choices[[x]])),
                select_choices[[x]]
              )
            )
          )
        )
      })
    })

    factor_reference_levels <- reactive(
      return(
        c(
          input$reference_level_factor_1,
				  input$reference_level_factor_2,
				  input$reference_level_factor_3,
				  input$reference_level_factor_4,
				  input$reference_level_factor_5,
				  input$reference_level_factor_6
        )
      )
    )

      # Observe submit button ------ 
    deg <- reactiveValues(limma = NULL)
    observeEvent(
      input$submit_model_button, {
        req(!is.null(pre_process$raw_counts()))
      
        shinybusy::show_modal_spinner(
          spin = "orbit",
          text = "Running Analysis",
          color = "#000000"
        )

        deg$limma <- limma_value(
          data_file_format = pre_process$data_file_format(),
          counts_deg_method = input$counts_deg_method,
          raw_counts = pre_process$raw_counts(),
          limma_p_val = input$limma_p_val,
          limma_fc = input$limma_fc,
          select_model_comprions = input$select_model_comprions,
          sample_info = pre_process$sample_info(),
          select_factors_model = input$select_factors_model,
          select_interactions = input$select_interactions,
          select_block_factors_model = input$select_block_factors_model,
          factor_reference_levels = factor_reference_levels(),
          processed_data = pre_process$data(),
          counts_log_start = pre_process$counts_log_start(),
          p_vals = pre_process$p_vals()
        )
        
        updateTabsetPanel(
          session = session, 
          inputId = "step_1", 
          selected = "results_tab"
        )
        shinybusy::remove_modal_spinner()
      }  
    )
    
    deg_info <- reactive({
      req(!is.null(deg$limma$results))
      
      deg_information(
        limma_value = deg$limma, 
        gene_names = pre_process$all_gene_names(),
        processed_data = pre_process$data(), 
        no_id_conversion = load_data$no_id_conversion()
      )[[1]]
    })
    
    deg_method <- c(
      "limma_trend", 
      "limma_voom", 
      "DESeq2"
    )
    
    name <- reactive({
      paste0(
        "deg_values_", 
        deg_method[as.numeric(input$counts_deg_method)], 
        ".csv"
      )
    })
    
    output$download_lfc <- downloadHandler(
      filename = function(){
        name()
      }, 
      content = function(file) {
        write.csv(deg_info(), file, row.names = FALSE)
      }
    )
    
    output$download_lfc_button <- renderUI({
      req(!is.null(deg_info()))
      downloadButton(
        outputId = ns("download_lfc"), 
        "Download DEG Data"
      )
    })
    
    output$note <- renderUI({
      req(!is.null(deg_info()))
      tippy::tippy_this(
        elementId = ns("download_lfc"),
        tooltip = "This data includes log fold change, adjusted p-value and 
            processed data from Pre-Process tab.",
        theme = "light-border"
      )
    })

    output$sig_gene_stats <- renderPlot({
      req(!is.null(deg$limma$results))
      sig_genes_plot(
        results = deg$limma$results
      )
    })

    output$sig_gene_stats_table <- renderTable({
      req(!is.null(deg$limma))
      
		  genes_stat_table(limma = deg$limma)
    },
      digits = 0,
      spacing = "s",
      include.rownames = FALSE,
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = T
    )

    output$list_comparisons_venn <- renderUI({
      req(!is.null(deg$limma))

      venn_comp <- list_comp_venn(
        limma = deg$limma,
        up_down_regulated = input$up_down_regulated
      )
      if(is.null(venn_comp$choices)) {
        selectInput(
          inputId = ns("select_comparisons_venn"),
          label = NULL,
          choices = list("All" = "All"),
          selected = "All"
        )
      } else {
        checkboxGroupInput(
          inputId = ns("select_comparisons_venn"), 
			    label = h4("Select up to 5 comparisons"), 
			    choices = venn_comp$choices,
			    selected = venn_comp$choices_first_three
        )
      }

	  })

    # venn diagram ----- 
    venn <- reactive({
      req(!is.null(deg$limma))
      req(!is.null(input$select_comparisons_venn))
      
      venn <- plot_venn(
        limma = deg$limma,
        up_down_regulated = input$up_down_regulated,
        select_comparisons_venn = input$select_comparisons_venn
      )
      p <- recordPlot()
      return(p)
    })
    output$venn_plot <- renderPlot({
      print(venn())
    })
    dl_venn <- ottoPlots::mod_download_figure_server(
      id = "dl_venn", 
      filename = "venn_diagram", 
      figure = reactive({ venn() })
    )
    

    # DEG STEP 2 --------
    output$list_comparisons <- renderUI({
      
      if(is.null(deg$limma$comparisons)) {
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
          choices = deg$limma$comparisons
	     )
      } 
	  })

    contrast_samples <- reactive({
      req(!is.null(input$select_contrast))

      find_contrast_samples(
        select_contrast = input$select_contrast, 
		    all_sample_names = colnames(pre_process$data()),
		    sample_info = pre_process$sample_info(),
		    select_factors_model = input$select_factors_model,
		    select_model_comprions = input$select_model_comprions, 
		    reference_levels = factor_reference_levels(),
		    counts_deg_method = input$counts_deg_method,
		    data_file_format = pre_process$data_file_format()
	    )
    })

    heat_data <- reactive({
      req(!is.null(deg$limma))
      req(!is.null(input$select_contrast))

      deg_heat_data(
        limma = deg$limma,
        select_contrast = input$select_contrast,
        processed_data = pre_process$data(),
        contrast_samples = contrast_samples()
      )
    })
    
    # Plot colors ------- 
    plot_colors <- list(
      "Green-Red" = c("green", "grey45", "red"), 
      "Red-Green" = c("red", "grey45", "green"), 
      "Blue-Red" = c("blue", "grey45", "red"), 
      "Green-Magenta" = c("green", "grey45", "magenta"),
      "Orange-Blue" = c("orange", "grey45", "blue")
    )
    
    plot_choices <- c(
      "Green-Red", 
      "Red-Green", 
      "Blue-Red", 
      "Green-Magenta", 
      "Orange-Blue"
    )
    
    observe({
      updateSelectInput(
        session = session, 
        inputId = "plot_color_select", 
        choices = plot_choices
      )
    })

    # Heatmap Colors ----------
    heatmap_colors <- list(
      "Green-Black-Red" = c("green", "black", "red"),
      "Red-Black-Green" = c("red", "black", "red"), 
      "Blue-White-Red" = c("blue", "white", "red"),
      "Green-Black-Magenta" = c("green", "black", "magenta"),
      "Blue-Yellow-Red" = c("blue", "yellow", "red"),
      "Blue-White-Brown" = c("blue", "white", "brown"), 
      "Orange-White-Blue" = c("orange", "white", "blue")
    )
    heatmap_choices <- c(
      "Green-Black-Red",
      "Red-Black-Green", 
      "Blue-White-Red",
      "Green-Black-Magenta",
      "Blue-Yellow-Red",
      "Blue-White-Brown", 
      "Orange-White-Blue"
    )
    observe({
      updateSelectInput(
        session = session,
        inputId = "heatmap_color_select",
        choices = heatmap_choices
      )
    })

    output$deg_main_heatmap <- renderPlot({
      req(!is.null(heat_data()$genes))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Creating Heatmap",
        color = "#000000"
      )

      # Assign heatmap to be used in multiple components
      deg_env$ht <- deg_heatmap(
        data = heat_data()$genes,
        bar = heat_data()$bar,
        heatmap_color_select = heatmap_colors[[input$heatmap_color_select]]
      )

      # Use heatmap position in multiple components
      deg_env$ht_pos_main <- InteractiveComplexHeatmap::htPositionsOnDevice(deg_env$ht)

      shinybusy::remove_modal_spinner()

      return(deg_env$ht)
    })

    output$deg_sub_heatmap <- renderPlot({
      if (is.null(input$ht_brush)) {
        grid::grid.newpage()
        grid::grid.text("No region is selected.", 0.5, 0.5)
      } else {
        deg_heat_return <- deg_heat_sub(
          ht_brush = input$ht_brush,
          ht = deg_env$ht,
          ht_pos_main = deg_env$ht_pos_main,
          heatmap_data = heat_data(),
          all_gene_names = pre_process$all_gene_names()
        )

        deg_env$ht_select <- deg_heat_return$ht_select
        deg_env$submap_data <- deg_heat_return$submap_data
        deg_env$group_colors <- deg_heat_return$group_colors
        deg_env$column_groups <- deg_heat_return$column_groups
        deg_env$bar <- deg_heat_return$bar
        
        deg_env$ht_sub <- ComplexHeatmap::draw(
          deg_env$ht_select,
          annotation_legend_side = "top",
          heatmap_legend_side = "top"
        )

        deg_env$ht_pos_sub <- InteractiveComplexHeatmap::htPositionsOnDevice(deg_env$ht_sub)

        return(deg_env$ht_sub)
      }
    })

    # Sub Heatmap Click Value ---------
    output$ht_click_content <- renderUI({
      if (is.null(input$ht_click)) { 
        "Click for Info."
      } else {
        deg_click_info(
          click = input$ht_click,
          ht_sub = deg_env$ht_sub,
          ht_sub_obj = deg_env$ht_select,
          ht_pos_sub = deg_env$ht_pos_sub,
          sub_groups = deg_env$column_groups,
          group_colors = deg_env$group_colors,
          bar = deg_env$bar,
          data = deg_env$submap_data
        )
      }
    })
    
    # volcano plot -----
    vol_plot <- reactive({
      req(!is.null(deg$limma$top_genes))
      
      vol <- plot_volcano(
        select_contrast = input$select_contrast,
        comparisons = deg$limma$comparisons,
        top_genes = deg$limma$top_genes,
        limma_p_val = input$limma_p_val,
        limma_fc = input$limma_fc,
        plot_colors = plot_colors[[input$plot_color_select]]
      )
    })
    
    

    output$volcano_plot <- renderPlot({
      print(vol_plot())
    })
    
    download_volcano <- ottoPlots::mod_download_figure_server(
      "download_volcano", 
      filename = "volcano_plot", 
      figure = reactive({ vol_plot() }) # stays as a reactive variable
    )
    
    # ma plot----------------
    ma_plot <- reactive({
      req(!is.null(deg$limma$top_genes))
      
      plot_ma(
        select_contrast = input$select_contrast,
        comparisons = deg$limma$comparisons,
        top_genes = deg$limma$top_genes,
        limma_p_val = input$limma_p_val,
        limma_fc = input$limma_fc,
        contrast_samples = contrast_samples(),
        processed_data = pre_process$data(), 
        plot_colors = plot_colors[[input$plot_color_select]]
      )
    })
    
    output$ma_plot <- renderPlot({
	    print(ma_plot())
    })
    
    download_ma <- ottoPlots::mod_download_figure_server(
      "download_ma", 
      filename = "ma_plot", 
      figure = reactive({ ma_plot() })
    )

    output$scatter_plot <- renderPlot({
      req(!is.null(deg$limma$top_genes))

      plot_deg_scatter(
        select_contrast = input$select_contrast,
        comparisons = deg$limma$comparisons,
        top_genes = deg$limma$top_genes,
        limma_p_val = input$limma_p_val,
        limma_fc = input$limma_fc,
        contrast_samples = contrast_samples(),
        processed_data = pre_process$data(),
        sample_info = pre_process$sample_info()
      )
    })

    # Split up and down genes into two data bases ---------
    up_reg_data <- reactive({
      req(!is.null(heat_data()))

      return(
        heat_data()$genes[heat_data()$bar == 1, ]
      )
    })
    down_reg_data <- reactive({
      req(!is.null(heat_data()))

      return(
        heat_data()$genes[heat_data()$bar == -1, ]
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

    # Enrichment Analysis Up Data -----------
    pathway_table_up <- reactive({
      req(!is.null(up_reg_data()))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Running Analysis",
        color = "#000000"
      )

      gene_names <- merge_data(
        all_gene_names = pre_process$all_gene_names(),
        data = up_reg_data(),
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

    # Subheatmap Data Table ----------
    output$pathway_data_up <- DT::renderDataTable({
      req(!is.null(pathway_table_up()))
      
      DT::datatable(
        if (ncol(pathway_table_up()) < 5){
          data = pathway_table_up()
        } else {
          data = pathway_table_up()[ ,1:4]
        },
        options = list(
          pageLength = 20,
          scrollX = "400px"
        ),
        rownames = TRUE
      )
    })
    output$dl_enrich_up <- downloadHandler(
      filename = function() {
        "deg_up_enrichment.csv"
      }, 
      content = function(file) {
        write.csv(data_frame_with_list(pathway_table_up()), file)
      }
    )

    # Enrichment Analysis Down Data -----------
    pathway_table_down <- reactive({
      req(!is.null(down_reg_data()))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Running Analysis",
        color = "#000000"
      )

      gene_names <- merge_data(
        all_gene_names = pre_process$all_gene_names(),
        data = down_reg_data(),
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

    # Subheatmap Data Table ----------
    output$pathway_data_down <- DT::renderDataTable({
      req(!is.null(pathway_table_down()))

      DT::datatable(
        if (ncol(pathway_table_down()) < 5){
          data = pathway_table_down()
        } else {
          data = pathway_table_down()[ ,1:4]
        },
        options = list(
          pageLength = 20,
          scrollX = "400px"
        ),
        rownames = TRUE
      )
    })
    output$dl_enrich_down <- downloadHandler(
      filename = function() {
        "deg_down_enrichment.csv"
      }, 
      content = function(file) {
        write.csv(data_frame_with_list(pathway_table_down()), file)
      }
    )

    go_table <- reactive({
      req(!is.null(pathway_table_up()) || !is.null(pathway_data_down()))

      go_table_data(
        up_enrich_data = pathway_table_up(),
        down_enrich_data = pathway_table_down()
      )
    })

    # Enrichment Tree -----------
    output$enrichment_tree <- renderPlot({
      req(!is.null(go_table()))
	    
      enrichment_plot(
        go_table = go_table(),
        45
      )
    })

    # Define a Network
    network_data_deg <- reactive({
      req(!is.null(go_table()))

      network_data(
        network = go_table(),
        up_down_reg_deg = input$up_down_reg_deg,
        wrap_text_network_deg= input$wrap_text_network_deg,
        layout_vis_deg = input$layout_vis_deg,
        edge_cutoff_deg = input$edge_cutoff_deg
      )
    })

    # Interactive vis network plot
    output$vis_network_deg <- visNetwork::renderVisNetwork({
      req(!is.null(network_data_deg()))
      
      vis_network_plot(
        network_data = network_data_deg()
      )
    })

    list(
      limma = reactive(deg$limma),
      select_factors_model = reactive(input$select_factors_model),
      select_model_comprions = reactive(input$select_model_comprions),
      reference_levels = reactive(factor_reference_levels()),
      counts_deg_method = reactive(input$counts_deg_method)
    )
  })
}

## To be copied in the UI
# mod_06_deg_ui("06_deg_ui_1")

## To be copied in the server
# mod_06_deg_server("06_deg_ui_1")
