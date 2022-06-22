#' 01_pre_process UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_02_pre_process_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Pre-Process",
    sidebarLayout(

      # Pre-Process Panel Sidebar ----------
      sidebarPanel(

        # Conditional panel for read count data -----------
        conditionalPanel(
          condition = "output.data_file_format == 1",
          strong("Keep genes with minimal counts per million (CPM) in at
                  least n libraries:"),
          fluidRow(
            column(
              width = 6,

              # Min counts per million (works with min samples)
              numericInput(
                inputId = ns("min_counts"),
                label = h5("Min. CPM"),
                value = 0.5
              )
            ),
            column(
              width = 6,

              # Min samples per row to have min CPM
              numericInput(
                inputId = ns("n_min_samples_count"),
                label = h5("n libraries"),
                value = 1
              )
            )
          ),
          tags$style(
            type = "text/css",
            "#pre_process-min_counts { width:100%;   margin-top:-12px}"
          ),
          tags$style(
            type = "text/css",
            "#pre_process-n_min_samples_count { width:100%;   margin-top:-12px}"
          ),

          # Type of transformation to perform on the counts data
          radioButtons(
            inputId = ns("counts_transform"),
            label = "Transform counts data for clustering & PCA.",
            choices = c(
              "VST: variance stabilizing transform" = 2,
              "rlog: regularized log (slow) " = 3,
              "EdgeR: log2(CPM+c)" = 1
            ),
            selected = 1
          ),

          # Conditional panel for EdgeR transformation -----------
          conditionalPanel(
            condition = "input.counts_transform == 1",
            fluidRow(
              column(
                width = 5,
                h5("Pseudo count c:")
              ),
              column(
                width = 7,

                # Constant to add for a log transform
                numericInput(
                  inputId = ns("counts_log_start"),
                  label = NULL,
                  value = 4
                )
              ) 
            ),
            ns = ns
          ),
          ns = ns
        ),

        # Conditional panel for FPKM data (2)----------
        conditionalPanel(
          condition = "output.data_file_format == 2",
          strong("Only keep genes above this level in at least n samples:"),
          fluidRow(
            column(
              width = 6,

              # Fold counts min (works with min samples)
              numericInput(
                inputId = ns("low_filter_fpkm"),
                label = h5("Min. level"),
                value = -1000
              )
            ),
            column(
              width = 6,

              # Min samples per row to have the low filter
              numericInput(
                inputId = ns("n_min_samples_fpkm"),
                label = h5("n samples"),
                value = 1
              )
            )
          ),
          tags$style(
            type = "text/css",
            "#pre_process-low_filter_fpkm { width:100%;margin-top:-12px}"
          ),
          tags$style(
            type = "text/css",
            "#pre_process-n_min_samples_fpkm { width:100%;margin-top:-12px}"
          ),

          # Perform a log transform or not
          radioButtons(
            inputId = ns("log_transform_fpkm"),
            label = "Log Transformation",
            choices = c("No" = FALSE, "Yes" = TRUE)
          ),

          # Constant to add if yes to a log transform
          numericInput(
            inputId = ns("log_start_fpkm"),
            label = h5("Constant c for started log: log(x+c)"),
            value = 1
          ),
          tags$style(
            type = "text/css",
            "#pre_process-log_start { width:100%;   margin-top:-12px}"
          ),
          ns = ns
        ),

        # Select input for missing value ------------
        selectInput(
          inputId = ns("missing_value"),
          label = "Missing values imputation:",
          choices = list(
            "Gene median" = "geneMedian",
            "Treat as zero" = "treatAsZero",
            "Median within sample groups" = "geneMedianInGroup"
          ),
          selected = "geneMedian"
        ),
        br(),

        strong("Download Processed Data"),
        br(),
        # Download button for processed data -----------
        downloadButton(
          outputId = ns("download_processed_data"),
          label = "Processed data"
        ),

        # Conditional panel for read count data ------------
        conditionalPanel(
          condition = "output.data_file_format == 1",

          # Download the counts data with converted IDs
          downloadButton(
            outputId = ns("download_converted_counts"),
            label = "Converted counts data"
          ),
          ns = ns
        ),
        br(),
        
        br(),
        strong("Download R-Markdown Report"),
        br(),
        
        downloadButton(
          outputId = ns("report"),
          label = "Generate report"
        ),  
        br(),
        br(),
        # Show transform messages
        actionButton(
          inputId = ns("show_messages"),
          label = "Show Conversion Messages"
        ),

        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/pre-process/",
          target = "_blank"
        ),
      ),
      

      # Pre-Process Panel Main -----------
      mainPanel(
        h5(
          "Aspect ratios of figures can be adjusted by changing
           the width of browser window. (To save a plot, right-click)"
        ),
       
        tabsetPanel(
          id = ns("eda_tabs"),

          # Barplot for read counts data ----------
          tabPanel(
            title = "Barplot",
            br(),
            plotOutput(
              outputId = ns("total_counts_gg"),
              width = "100%",
              height = "500px"
            ), 
            ottoPlots::mod_download_figure_ui(
              id = ns("dl_total_counts"), 
              label = "Download barplot"
            )
          ),

          # Scatterplot with interactive axes ----------
          tabPanel(
            title = "Scatterplot",
          # Axis selectors -----------
            br(),
            fluidRow(
              column(
                width = 4,
                selectInput(
                  inputId = ns("scatter_x"),
                  label = "Select a sample for x-axis",
                  choices = 1:5,
                  selected = 1
                )
              ),
              column(
                width = 4,
                selectInput(
                  inputId = ns("scatter_y"),
                  label = "Select a sample for y-axis",
                  choices = 1:5,
                  selected = 2
                )
              )
            ),
            br(),
            plotOutput(
              outputId = ns("eda_scatter"),
              width = "100%",
              height = "500px"
            ), 
            ottoPlots::mod_download_figure_ui(
              id = ns("dl_eda_scatter"), 
              label = "Download scatterplot"
            )
          ),

          # Boxplot of transformed data ----------
          tabPanel(
            title = "Boxplot",
            br(),
            plotOutput(
              outputId = ns("eda_boxplot"),
              width = "100%",
              height = "500px"
            ), 
            ottoPlots::mod_download_figure_ui(
              id = ns("dl_eda_boxplot"),
              label = "Download boxplot"
            )
          ),

          # Density plot of transformed data ---------
          tabPanel(
            title = "Density Plot",
            br(),
            plotOutput(
              outputId = ns("eda_density"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(
              id = ns("dl_eda_density"), 
              label = "Download density plot"
            )
          ),

          # Density plot of transformed data ---------
          tabPanel(
            title = "SD vs. Mean Plot",
            br(),
            fluidRow(
              column(
                width = 4,
                selectInput(
                  inputId = ns("heat_color_select"),
                  label = "Select Heat Colors",
                  choices = NULL
                )
              ),
              column(
                width = 4,
                checkboxInput(
                  inputId = ns("rank"),
                  label = "Use rank of mean values"
                )
              ),
            ),
            plotOutput(
              outputId = ns("dev_transfrom"),
              width = "100%",
              height = "500px"
            ), 
            ottoPlots::mod_download_figure_ui(
              id = ns("dl_dev_transform"), 
              label = "Download transformed plot"
            )
          ),

          # Searchable table of transformed converted data ---------
          tabPanel(
            title = "Converted Data",
            br(),
            DT::dataTableOutput(outputId = ns("examine_data"))
          ),

          # Plot panel for individual genes ---------
          tabPanel(
            title = "Individual Genes",
            br(),
            fluidRow(
              column(
                4, 
                # Gene ID Selection -----------
                selectInput(
                  inputId = ns("select_gene_id"),
                  label = "Select Gene ID Label",
                  choices = NULL,
                  selected = NULL
                ),
                selectizeInput(
                  inputId = ns("selected_gene"),
                  label = "Select/Search for Genes",
                  choices = "",
                  selected = NULL,
                  multiple = TRUE
                )
              ),
              column(
                4,
                checkboxInput(
                  inputId = ns("gene_plot_box"),
                  label = "Show individual samples",
                  value = FALSE
                ),
                uiOutput(ns("sd_checkbox"))
              ),
              column(
                4,
                radioButtons(
                  inputId = ns("angle_ind_axis_lab"),
                  label = "Angle Axis Labels",
                  choices = c(0, 45, 90),
                  selected = 0
                )
              )
            ),
            plotOutput(
              outputId = ns("gene_plot"),
              width = "100%",
              height = "500px"
            ), 
            ottoPlots::mod_download_figure_ui(
              id = ns("dl_gene_plot"), 
              label = "Download gene plot"
            )
          )
        )
      )
    )
    )
}

#' 01_pre_process Server Functions
#'
#' @noRd
mod_02_pre_process_server <- function(id, load_data, tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Data file format for conditional panels ----------
    # outputOptions required otherwise the value can only be used
    # if it is rendered somewhere else in the UI
    output$data_file_format <- reactive({
      load_data$data_file_format()
    })
    outputOptions(output, "data_file_format", suspendWhenHidden = FALSE)

    # Update Variable Selection for the Scatter Plots ----------
    observe({
      req(!is.null(load_data$converted_data()))

      updateSelectInput(
        session,
        inputId = "scatter_x",
        choices = colnames(processed_data()$data),
        selected = colnames(processed_data()$data)[1]
      #load_data$converted_data())[1]
      )
      updateSelectInput(
        session,
        inputId = "scatter_y",
        choices = colnames(processed_data()$data),
        selected = colnames(processed_data()$data)[2]
      )
    })

    # Dynamic Barplot Tab ----------
    observe({
      if (load_data$data_file_format() != 1) {
        hideTab(inputId = "eda_tabs", target = "Barplot")
        updateTabsetPanel(session, "eda_tabs", selected = "Scatterplot")
      } else if (load_data$data_file_format() == 1) {
        showTab(inputId = "eda_tabs", target = "Barplot")
        updateTabsetPanel(session, "eda_tabs", selected = "Barplot")
      }
    })

    # Process the data with user defined criteria ----------
    processed_data <- reactive({
      req(!is.null(load_data$converted_data()))
      req(input$n_min_samples_count)
      req(input$min_counts)
      req(input$low_filter_fpkm)
      req(input$n_min_samples_fpkm)
      req(input$counts_log_start)
      req(input$log_start_fpkm)

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Pre-Processing Data",
        color = "#000000"
      )

      processed_data <- pre_process(
        data = load_data$converted_data(),
        missing_value = input$missing_value,
        data_file_format = load_data$data_file_format(),
        low_filter_fpkm = input$low_filter_fpkm,
        n_min_samples_fpkm = input$n_min_samples_fpkm,
        log_transform_fpkm = input$log_transform_fpkm,
        log_start_fpkm = input$log_start_fpkm,
        min_counts = input$min_counts,
        n_min_samples_count = input$n_min_samples_count,
        counts_transform = input$counts_transform,
        counts_log_start = input$counts_log_start,
        no_fdr = load_data$no_fdr()
      )
      shinybusy::remove_modal_spinner()
      
      return(processed_data)
    })

    # Counts barplot ------------
    total_counts <- reactive({
      req(!is.null(processed_data()$data))
      
      total_counts_ggplot(
        counts_data = processed_data()$raw_counts,
        sample_info = load_data$sample_info()
      )
    })
    output$total_counts_gg <- renderPlot({
      print(total_counts())
    })
    dl_total_counts <- ottoPlots::mod_download_figure_server(
      id = "dl_total_counts", 
      filename = "total_counts_barplot", 
      figure = reactive({ total_counts() })
    )

    # Scatter eda plot ----------
    scatter <- reactive({
      req(!is.null(processed_data()$data))
      
      eda_scatter(
        processed_data = processed_data()$data,
        plot_xaxis = input$scatter_x,
        plot_yaxis = input$scatter_y
      )
    })
    output$eda_scatter <- renderPlot({
      print(scatter())
    })
    dl_eda_scatter <- ottoPlots::mod_download_figure_server(
      id = "dl_eda_scatter", 
      filename = "scatter_plot", 
      figure = reactive({ scatter() })
    )

    # Box eda plot ----------
    eda_box <- reactive({
      req(!is.null(processed_data()$data))
      
      eda_boxplot(
        processed_data = processed_data()$data,
        sample_info = load_data$sample_info()
      )
    })
    output$eda_boxplot <- renderPlot({
      print(eda_box())
    })
    dl_eda_boxplot <- ottoPlots::mod_download_figure_server(
      id = "dl_eda_boxplot", 
      filename = "transformed_boxplot", 
      figure = reactive({ eda_box() })
    )

    # Density eda plot ----------
    density <- reactive({
      req(!is.null(processed_data()$data))
      
      eda_density(
        processed_data = processed_data()$data,
        sample_info = load_data$sample_info()
      )
    })
    output$eda_density <- renderPlot({
      print(density())
    })
    dl_eda_density <- ottoPlots::mod_download_figure_server(
      id = "dl_eda_density", 
      filename = "density_plot", 
      figure = reactive({ density() })
    )

    # Standard deviation vs mean plot ----------
    # Heatmap Colors ----------
    heat_colors <- list(
      "Green" = c("green"),
      "Red" = c("red"),
      "Magenta" = c("magenta"),
      "Blue" = c("blue"),
      "Brown" = c("brown")
    )
    heat_choices <- c(
      "Green",
      "Red",
      "Magenta",
      "Blue",
      "Brown"
    )
    observe({
      updateSelectInput(
        session = session,
        inputId = "heat_color_select",
        choices = heat_choices
      )
    })
    
    # Mean vs SD plot --------
    dev <- reactive({
      req(!is.null(processed_data()$data))
      
      mean_sd_plot(
        processed_data = processed_data()$data,
        heat_cols = heat_colors[[input$heat_color_select]],
        rank = input$rank
      )
    })
    output$dev_transfrom <- renderPlot({
      print(dev())
    })
    dl_dev_transform <- ottoPlots::mod_download_figure_server(
      id = "dl_dev_transform", 
      filename = "transform_plot", 
      figure = reactive({ dev() })
    )

    # Merge Data Sets with Gene names ----------
    merged_processed_data <- reactive({
      req(!is.null(processed_data()$data))

      merged_data <- merge_data(
        load_data$all_gene_names(),
        processed_data()$data,
        merge_ID = "ensembl_ID"
      )
    })
    merged_raw_counts_data <- reactive({
      req(!is.null(processed_data()$data))

      merged_data <- merge_data(
        load_data$all_gene_names(),
        processed_data()$raw_counts,
        merge_ID = "ensembl_ID"
      )
    })

    # Pre-Process Data Table ----------
    output$examine_data <- DT::renderDataTable({
      req(!is.null(merged_processed_data()))

      DT::datatable(
        merged_processed_data(),
        options = list(
          pageLength = 10,
          scrollX = "400px"
        ),
        rownames = FALSE
      )
    })

    # Gene ID Name Choices ----------
    observe({
      req(!is.null(load_data$all_gene_names()))

      updateSelectInput(
        session = session,
        inputId = "select_gene_id",
        choices = colnames(load_data$all_gene_names())
      )
    })

    # Individual plot data ------------
    individual_data <- reactive({
      req(!is.null(processed_data()$data))

      rowname_id_swap(
        data_matrix = processed_data()$data,
        all_gene_names = load_data$all_gene_names(),
        select_gene_id = input$select_gene_id
      )
    })

    # Individual genes selection ----------
    observe({
      req(tab() == "Pre-Process")

      updateSelectizeInput(
        session,
        inputId = "selected_gene",
        choices = rownames(individual_data()),
        selected = NULL,
        server = TRUE
      )
    })

    # Dynamic individual gene checkbox ----------
    output$sd_checkbox <- renderUI({
      req(input$gene_plot_box == FALSE)

      checkboxInput(
        inputId = ns("use_sd"),
        label = "Use standard deviation instead of standard error",
        value = FALSE
      )
    })

    # Individual gene plot ---------
    gene_plot <- reactive({
      req(!is.null(individual_data()))
      req(!is.null(input$selected_gene))
      
      individual_plots(
        individual_data = individual_data(),
        sample_info = load_data$sample_info(),
        selected_gene = input$selected_gene,
        gene_plot_box = input$gene_plot_box,
        use_sd = input$use_sd,
        lab_rotate = input$angle_ind_axis_lab
      )
    })
    output$gene_plot <- renderPlot({
      print(gene_plot())
    })
    dl_gene_plot <- ottoPlots::mod_download_figure_server(
      id = "dl_gene_plot", 
      filename = "gene_plot", 
      figure = reactive({ gene_plot() })
    )
    
    

    # Download buttons ----------
    output$download_processed_data <- downloadHandler(
      filename = function() {
        "processed_data.csv"
      },
      content = function(file) {
        write.csv(merged_processed_data(), file)
      }
    )
    output$download_converted_counts <- downloadHandler(
      filename = function() {
        "converted_counts_data.csv"
      },
      content = function(file) {
        write.csv(merged_raw_counts_data(), file)
      }
    )

    # Markdown report
    output$report <- downloadHandler(
      
      # For PDF output, change this to "report.pdf"
      filename ="pre_process_report.pdf",
      content = function(file) {
        #Show Loading popup
        shinybusy::show_modal_spinner(
          spin = "orbit",
          text = "Generating Report",
          color = "#000000"
        )
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "test_workflow2.Rmd")
        tempReport
        tempReport<-gsub("\\", "/",tempReport,fixed = TRUE)

        #This should retrieve the project location on your device:
        #"C:/Users/bdere/Documents/GitHub/idepGolem"
        wd <- getwd()

        markdown_location <-paste0(wd, "/vignettes/Reports/test_workflow2.Rmd")
        file.copy(from=markdown_location,to = tempReport, overwrite = TRUE)

        # Set up parameters to pass to Rmd document
        params <- list(
          loaded_data = load_data$converted_data(),
          sample_info = load_data$sample_info(),
          data_file_format = load_data$data_file_format(),
          no_id_conversion = input$no_id_conversion,
          min_counts = input$min_counts,
          n_min_samples_count = input$n_min_samples_count,
          counts_transform = input$counts_transform,
          counts_log_start = input$counts_log_start,
          log_transform_fpkm = input$log_transform_fpkm,
          log_start_fpkm = input$log_start_fpkm,
          low_filter_fpkm = input$low_filter_fpkm,
          missing_value = input$missing_value,
          scatter_x = input$scatter_x,
          scatter_y = input$scatter_y,
          sd_color = heat_colors[[input$heat_color_select]],
          rank = input$rank,
          no_fdr = load_data$no_fdr()

        )

        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(
          input = tempReport,#markdown_location, 
          output_file = file,
          params = params,
          envir = new.env(parent = globalenv())
        )
        shinybusy::remove_modal_spinner()
        
      }
      
    )

    # Number of converted IDs ---------
    n_matched <- reactive({
      req(!is.null(processed_data))

      match_process_ids <-
        load_data$matched_ids() %in% rownames(processed_data()$data)

      return(sum(match_process_ids))
    })

    # Bias detected message -------
    read_counts_bias <- reactive({
      req(!is.null(processed_data()$raw_counts))

      counts_bias_message(
        raw_counts = processed_data()$raw_counts,
        sample_info = load_data$sample_info()
      )
    })

    # Text Output Information -----------
    converted_message <- reactive({
      req(processed_data()$data_size)

      conversion_counts_message(
        data_size = processed_data()$data_size,
        all_gene_names = load_data$all_gene_names(),
        n_matched = n_matched()
      )
    })

    # Show messages when on the Pre-Process tab or button is clicked
    observe({
      req(input$show_messages || tab() == "Pre-Process")

      showNotification(
        ui = converted_message(),
        id = "conversion_counts",
        duration = NULL,
        type = "default"
      )

      req(!is.null(read_counts_bias()))
      showNotification(
        ui = read_counts_bias(),
        id = "read_counts_message",
        duration = NULL,
        type = "error"
      )
    })
    # Data type warning -------
    observe({
      req(input$show_messages || tab() == "Pre-Process")
      req(processed_data()$data_type_warning != 0)
      
      message <- switch(processed_data()$data_type_warning,
             "-1" = "Integers detected. Did you mean to select 'read counts'?",
             "1" = "Non count values detected. Did you mean select 'Normalized expression values'?"
  )
  
      
      showNotification(
        ui = message,
        id = "data_type_warning",
        duration = NULL,
        type = "error"
      )
    })

    
    # Remove messages if the tab changes --------
    observe({
      req(tab() != "Pre-Process")

      removeNotification("conversion_counts")
      removeNotification("read_counts_message")
      removeNotification("data_type_warning")
    })

    all_gene_info <- reactive({
      req(!is.null(load_data$converted()))

      return(
        get_gene_info(
          load_data$converted(),
          load_data$select_org(),
          all_gene_info = idep_data$gene_info_files
        )
      ) 
		})


    # Return Values -----------
    list(
      raw_counts = reactive(processed_data()$raw_counts),
      data = reactive(processed_data()$data),
      p_vals = reactive(processed_data()$p_vals),
      sample_info = reactive(load_data$sample_info()),
      all_gene_names = reactive(load_data$all_gene_names()),
      gmt_choices = reactive(load_data$gmt_choices()),
      converted = reactive(load_data$converted()),
      select_org = reactive(load_data$select_org()),
      gmt_file = reactive(load_data$gmt_file()),
      all_gene_info = reactive(load_data$all_gene_info()),
      data_file_format = reactive(load_data$data_file_format()),
      counts_log_start = reactive(input$counts_log_start),
      all_gene_info = reactive(all_gene_info())
    )
  })
}

## To be copied in the UI
# mod_02_pre_process_ui("pre_process") #nolint

## To be copied in the server
# mod_02_pre_process_server("pre_process_ui") #nolint
