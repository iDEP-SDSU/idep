#' 03_heatmap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_03_clustering_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Clustering",
    sidebarLayout(

      # Heatmap Panel Sidebar ----------
      sidebarPanel(
        conditionalPanel(
          condition = "input.cluster_panels == 'Heatmap/Enrichment' | 
          input.cluster_panels == 'Gene SD Distribution' ",
          
          numericInput(
            inputId = ns("n_genes"), 
            label = h4("Top n most variable genes to include:"), 
            min = 10, 
            max = 12000, 
            value = 100, 
            step = 10
          ), 
          ns = ns
        ),
        
        conditionalPanel(
          condition = "(input.cluster_panels == 'Heatmap/Enrichment' | 
            input.cluster_panels == 'Sample Tree') &&  input.cluster_meth == 2",
          
          # k- means slidebar -----------
            
            sliderInput(
              inputId = ns("k_clusters"),
              label = "Number of Clusters:",
              min   = 2,
              max   = 20,
              value = 4,
              step  = 1
            ),
            
            # Re-run k-means with a different seed
            actionButton(
              inputId = ns("k_means_re_run"),
              label = "Re-Run"
            ),
            
            # Elbow plot pop-up 
            actionButton(
              inputId = ns("elbow_pop_up"),
              label = "How many clusters?"
            ),
            # Line break ---------
            HTML(
              '<hr style="height:1px;border:none;
           color:#333;background-color:#333;" />'
            ),
          
          ns = ns 
        ),

      

        # Select Clustering Method ----------
        conditionalPanel(
          condition = "input.cluster_panels == 'Heatmap/Enrichment' | 
            input.cluster_panels == 'Sample Tree'",
          
          selectInput(
            inputId = ns("cluster_meth"),
            label = "Select Clustering Method:",
            choices = list(
              "Hierarchical" = 1,
              "k-Means" = 2
            ),
            selected = 1
          ),

          ns = ns
        ),

        # Heatmap customizing features ----------
        conditionalPanel(
          condition = "input.cluster_panels == 'Heatmap/Enrichment' ",
          
          # Gene ID Selection -----------
          selectInput(
            inputId = ns("select_gene_id"),
            label = "Select Gene ID Label (<= 50 genes):",
            choices = NULL,
            selected = NULL
          ),
          
          # Sample coloring bar -----------
          htmlOutput(ns("list_factors_heatmap")),

          strong("Customize heatmap (Default values work well):"),
          fluidRow(
            br(),
            column(width = 3, h5("Color")),
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

        # Clustering methods for hierarchical ----------
        conditionalPanel(
          condition = "input.cluster_meth == 1 && 
            (input.cluster_panels == 'Heatmap/Enrichment' | 
            input.cluster_panels == 'Sample Tree')",
          fluidRow(
            column(width = 4, h5("Distance")),
            column(
              width = 8,
              selectInput(
                inputId = ns("dist_function"),
                label = NULL,
                choices = NULL,
                width = "100%"
              )
            )
          ),
          fluidRow(
            column(width = 4, h5("Linkage")),
            column(
              width = 8,
              selectInput(
                inputId = ns("hclust_function"),
                label = NULL,
                choices = c(
                  "average", "complete", "single",
                  "median", "centroid", "mcquitty"
                ),
                width = "100%"
              )
            )
          ),
          fluidRow(
            column(width = 8, h5("Cut-off Z score")),
            column(
              width = 4,
              numericInput(
                inputId = ns("heatmap_cutoff"),
                label = NULL,
                value = 4,
                min = 2,
                step = 1
              )
            )
          ),
          ns = ns
        ),

        # Checkbox features ------------
        conditionalPanel(
          condition = "input.cluster_panels == 'Heatmap/Enrichment' | input.cluster_panels == 'Sample Tree' ",
          
          checkboxInput(
            inputId = ns("gene_centering"),
            label = "Center genes (substract mean)",
            value = TRUE
          ),
          checkboxInput(
            inputId = ns("gene_normalize"),
            label = "Normalize genes (divide by SD)",
            value = FALSE
          ),
          ns = ns
        ),
        
        conditionalPanel(
          condition = "input.cluster_panels == 'Heatmap/Enrichment' ",
          
          checkboxInput(
            inputId = ns("no_sample_clustering"),
            label = "Do not re-order or cluster samples",
            value = FALSE
          ),
          checkboxInput(
            inputId = ns("show_row_dend"),
            label = "Show Row Dendogram",
            value = TRUE
          ),
          br(),
          downloadButton(
            outputId = ns("download_heatmap_data"),
            label = "Heatmap data"
          ),
          
          ns = ns
        ),

        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/heatmap/",
          target = "_blank"
        )
      ),
      mainPanel(
        tabsetPanel(
          id = ns("cluster_panels"),

          # Heatmap panel ----------
          tabPanel(
            title = "Heatmap/Enrichment",
            h3("Heatmap"), 
            h5("Brush for sub-heatmap, click for value. (Shown Below)"),
            br(),
            
            fluidRow(
              column(
                width = 3,
                plotOutput(
                  outputId = ns("heatmap_main"),
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
                  outputId = ns("sub_heatmap"),
                  height = "650px",
                  width = "100%",
                  click = ns("ht_click")
                )
              )
            ),
            h3("Enrichment"), 
            h5("Enrichment analysis is  based on selected genes from heatmap."),
            h6("List of genes included in each pathway can be found 
               in downloaded data."),
            fluidRow(
              column(
                width = 4,
                htmlOutput(outputId = ns("select_go_selector"))
              ),
              column(
                width = 8,
                checkboxInput(
                  inputId = ns("filtered_background"), 
                  label = "Use filtered data as background in enrichment (slow)", 
                  value = TRUE
                ),
                checkboxInput(
                  inputId = ns("remove_redudant"),
                  label = "Remove Redudant Gene Sets",
                  value = FALSE
                )
              ),
              tags$style(
                type='text/css',
                "#clustering-min_set_size {width:100%; margin-top:-12px}"
              ),
              tags$style(
                type='text/css',
                "#clustering-max_set_size {width:100%; margin-top:-12px}"
              )
            ),
            verbatimTextOutput(ns("test")),
            uiOutput(outputId = ns("pathway_data"))
          ),
          
          # Gene Standard Deviation Distribution ----------
          tabPanel(
            title = "Gene SD Distribution",
            br(),
            plotOutput(
              outputId = ns("sd_density_plot"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(ns("dl_gene_dist"))
          ),

          # Sample Tree Plot ---------
          tabPanel(
            title = "Sample Tree",
            h5(
              "Using genes with maximum expression level at the top 75%.
               Data is transformed and clustered as specified in the sidebar."
            ),
            br(),
            plotOutput(
              outputId = ns("sample_tree"),
              width = "100%",
              height = "400px"
            ),
            ottoPlots::mod_download_figure_ui(ns("dl_sample_tree"))
          )
        )
      )
    )
  )
}


#' 03_heatmap Server Functions
#'
#' @noRd
mod_03_clustering_server <- function(id, pre_process, idep_data, tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Interactive heatmap environment
    shiny_env <- new.env()

    # Update Slider Input ---------
    observe({
      req(tab() == "Clustering")
      req(!is.null(pre_process$data()))
      if (nrow(pre_process$data()) > 12000) {
        max_genes <- 12000
      } else {
        max_genes <- round(nrow(pre_process$data()), -2)
      }
      updateNumericInput(
        inputId = "n_genes", 
        value = 100, 
        max = max_genes
      )
    })

    # Heatmap Colors ----------
    heatmap_colors <- list(
      "Green-Black-Red" = c("green", "black", "red"),
      "Red-Black-Green" = c("red", "black", "green"), 
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

    # Distance functions -----------
    dist_funs <- dist_functions()
    dist_choices <- setNames(
      1:length(dist_funs),
      names(dist_funs)
    )
    observe({
      updateSelectInput(
        session = session,
        inputId = "dist_function",
        choices = dist_choices
      )
    })

    # Hclust Functions ----------
    hclust_funs <- hcluster_functions()

    # Gene ID Name Choices ----------
    observe({
      req(!is.null(pre_process$all_gene_names()))

      updateSelectInput(
        session = session,
        inputId = "select_gene_id",
        choices = colnames(pre_process$all_gene_names())
      )
    })

    # Sample color bar selector ----------
    output$list_factors_heatmap <- renderUI({
      selectInput(
        inputId = ns("select_factors_heatmap"),
        label = "Sample Color Bar:",
        choices = c("Sample_Name", colnames(pre_process$sample_info()))
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

    # Standard Deviation Density Plot ----------
    sd_density_plot <- reactive({
      req(!is.null(pre_process$data()))
      
      sd_density(
        data = pre_process$data(),
        n_genes_max = input$n_genes
      )
    })
    
    output$sd_density_plot <- renderPlot({
      print(sd_density_plot())
    })
    
    dl_gene_dist <- ottoPlots::mod_download_figure_server(
      id = "dl_gene_dist", 
      filename = "sd_density_plot", 
      figure = reactive({ sd_density_plot() })
    )
    
    

    # Heatmap Data -----------
    heatmap_data <- reactive({
      req(!is.null(pre_process$data()))

      process_heatmap_data(
        data = pre_process$data(),
        n_genes_max = input$n_genes,
        gene_centering = input$gene_centering,
        gene_normalize = input$gene_normalize,
        sample_centering = FALSE,
        sample_normalize = FALSE,
        all_gene_names = pre_process$all_gene_names(),
        select_gene_id = input$select_gene_id
      )
    })

    # HEATMAP -----------
    # Information on interactivity
    # https://jokergoo.github.io/2020/05/15/interactive-complexheatmap/
    output$heatmap_main <- renderPlot({
      req(!is.null(heatmap_data()))
      req(!is.null(input$select_factors_heatmap))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Creating Heatmap",
        color = "#000000"
      )

      # Assign heatmap to be used in multiple components
      shiny_env$ht <- heatmap_main(
        data = heatmap_data(),
        cluster_meth = input$cluster_meth,
        heatmap_cutoff = input$heatmap_cutoff,
        sample_info = pre_process$sample_info(),
        select_factors_heatmap = input$select_factors_heatmap,
        dist_funs = dist_funs,
        dist_function = input$dist_function,
        hclust_function = input$hclust_function,
        no_sample_clustering = input$no_sample_clustering,
        heatmap_color_select = heatmap_colors[[input$heatmap_color_select]],
        row_dend = input$show_row_dend,
        k_clusters = input$k_clusters,
        re_run = input$k_means_re_run
      )

      # Use heatmap position in multiple components
      shiny_env$ht_pos_main <- InteractiveComplexHeatmap::htPositionsOnDevice(shiny_env$ht)

      shinybusy::remove_modal_spinner()

      return(shiny_env$ht)
    })

    # Heatmap Click Value ---------
    output$ht_click_content <- renderUI({
      if (is.null(input$ht_click)) { 
        "Click for Info."
      } else {
        cluster_heat_click_info(
          click = input$ht_click,
          ht_sub = shiny_env$ht_sub,
          ht_sub_obj = shiny_env$ht_sub_obj,
          ht_pos_sub = shiny_env$ht_pos_sub,
          sub_groups = shiny_env$sub_groups,
          group_colors = shiny_env$group_colors,
          cluster_meth = input$cluster_meth,
          click_data = shiny_env$click_data
        )
      }
    })

    # Subheatmap creation ---------
    output$sub_heatmap <- renderPlot({
      if (is.null(input$ht_brush)) {
        grid::grid.newpage()
        grid::grid.text("No region is selected.", 0.5, 0.5)
      } else {
        submap_return <- heat_sub(
          ht_brush = input$ht_brush,
          ht = shiny_env$ht,
          ht_pos_main = shiny_env$ht_pos_main,
          heatmap_data = heatmap_data(),
          sample_info = pre_process$sample_info(),
          select_factors_heatmap = input$select_factors_heatmap,
          cluster_meth = input$cluster_meth
        )

        # Objects used in other components ----------
        shiny_env$ht_sub_obj <- submap_return$ht_select
        shiny_env$submap_data <- submap_return$submap_data
        shiny_env$sub_groups <- submap_return$sub_groups
        shiny_env$group_colors <- submap_return$group_colors
        shiny_env$click_data <- submap_return$click_data
        
        shiny_env$ht_sub <- ComplexHeatmap::draw(
          shiny_env$ht_sub_obj,
          annotation_legend_list = submap_return$lgd,
          annotation_legend_side = "top"
        )

        shiny_env$ht_pos_sub <- InteractiveComplexHeatmap::htPositionsOnDevice(shiny_env$ht_sub)

        return(shiny_env$ht_sub)
      }
    })

    # Enrichment Analysis ----------
    # Gene sets reactive
    pathway_table <- reactive({
      req(!is.null(input$select_gene_id))
      req(!is.null(input$ht_brush))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Running Analysis",
        color = "#000000"
      )
      
      pathway_info <- list()
      
      if (input$cluster_meth == 1) {
        gene_names <- merge_data(
          all_gene_names = pre_process$all_gene_names(),
          data = shiny_env$submap_data,
          merge_ID = input$select_gene_id
        )
      
        # Only keep the gene names and scrap the data
        gene_names_query <- dplyr::select_if(gene_names, is.character)

        req(!is.null(pre_process$all_gene_names()))
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

        pathway_info[["Hierarchical_Selection"]] <- find_overlap(
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
      } else if (input$cluster_meth == 2) {
        # Get the cluster number and Gene IDs
        row_ord <- ComplexHeatmap::row_order(shiny_env$ht)
        for (i in 1:length(row_ord)) {
          if (i == 1) {
          clusts <- data.frame(
            "cluster" = rep(names(row_ord[i]), length(row_ord[[i]])),
            "row_order" = row_ord[[i]]
          )
          } else {
          tem <- data.frame(
            "cluster" = rep(names(row_ord[i]), length(row_ord[[i]])),
            "row_order" = row_ord[[i]]
          )
          clusts <- rbind(clusts, tem)
          }
        }
        clusts$id <- rownames(heatmap_data()[clusts$row_order, ]) 

        for (i in 1:length(shiny_env$click_data)) {
          cluster_data <- shiny_env$click_data[[i]]

          gene_names <- merge_data(
            all_gene_names = pre_process$all_gene_names(),
            data = cluster_data,
            merge_ID = input$select_gene_id
          )
      
          # Only keep the gene names and scrap the data
          gene_names_query <- dplyr::select_if(gene_names, is.character)

          req(!is.null(pre_process$all_gene_names()))
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

          pathway_sub_info <- find_overlap(
            pathway_table = gene_sets$pathway_table,
            query_set = gene_sets$query_set,
            total_genes = gene_sets$total_genes,
            processed_data = pre_process$data(),
            gene_info = pre_process$all_gene_info(),
            go = input$select_go,
            idep_data = idep_data,
            select_org = pre_process$select_org(),
            sub_pathway_files = gene_sets$pathway_files,
            use_filtered_background = TRUE,
            reduced = FALSE
          )

          # Get cluster by matching gene ID from query to cluster number
          clust_num <- clusts$cluster[clusts$id == gene_names_query[1, 2]]

          pathway_info[[paste0("Cluster_", clust_num)]] <- pathway_sub_info
        }
      }

      shinybusy::remove_modal_spinner()
      


      return(pathway_info)
    })
    
    

    # Pathway Data Table ----------
    output$pathway_data <- renderUI({
      req(!is.null(pathway_table()))
      
      
      #exclude gene list column from displayed table, but keep for download
      lapply(names(pathway_table()), function(x) {
        output[[x]] = DT::renderDataTable({
          DT::datatable(
            if (ncol(pathway_table()[[x]]) < 5){
              data = pathway_table()[[x]]
            } else {
              data = pathway_table()[[x]][,1:4]
            },
            options = list(
              pageLength = 20,
              scrollX = "400px",
              dom = 'ft'
            ),
            rownames = FALSE,
            selection = 'single'
          )
        })

        down_data <- data_frame_with_list(pathway_table()[[x]])

        output[[paste0("table_", x)]] = downloadHandler(
          filename = function() {
            paste0(x, ".csv")
          },
          content = function(file) {
            write.csv(down_data, file)
          }
        )
      })
  
      return(lapply(names(pathway_table()), function(x) {
        tagList(
          br(),
          downloadButton(
            outputId = ns(paste0("table_", x)),
            label = paste0("Enrichment: ", gsub("_", " ", x))
          ),
          br(),
          strong(h3(gsub("_", " ", x))),
          DT::dataTableOutput(ns(x))
        )
      })
      )
    })

  
    # Sample Tree ----------
    sample_tree <- reactive({
      req(!is.null(pre_process$data()))
      
      draw_sample_tree(
        tree_data = pre_process$data(),
        gene_centering = input$gene_centering,
        gene_normalize = input$gene_normalize,
        sample_centering = FALSE,
        sample_normalize = FALSE,
        hclust_funs = hclust_funs,
        hclust_function = input$hclust_function,
        dist_funs = dist_funs,
        dist_function = input$dist_function
      )  
      p <- recordPlot() 
      return(p)
    })
    
    output$sample_tree <- renderPlot({
      print(sample_tree())
    })
    
    dl_sample_tree <- ottoPlots::mod_download_figure_server(
      id = "dl_sample_tree", 
      filename = "sample_tree", 
      figure = reactive({ sample_tree() })
    )
    
    # k-Cluster elbow plot ----------
    output$k_clusters <- renderPlot({
      req(!is.null(heatmap_data()))

      k_means_elbow(
        heatmap_data = heatmap_data()
      )
    })
    # pop-up modal 
    observeEvent(input$elbow_pop_up, {
      showModal(modalDialog(
        plotOutput(ns("k_clusters")), 
        footer = NULL, 
        easyClose = TRUE, 
        title = tags$h5(
          "Following the elbow method, one should choose k so that adding 
          another cluster does not substantially reduce the within groups sum of squares.",
          tags$a(
            "Wikipedia",
            href = "https://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set",
            target = "_blank"
          )
        ),  
      ))
    })

     # Heatmap Download Data -----------
    heatmap_data_download <- reactive({
      req(!is.null(pre_process$all_gene_names()))
      req(!is.null(heatmap_data()))

      merged_data <- merge_data(
        pre_process$all_gene_names(),
        heatmap_data(),
        merge_ID = input$select_gene_id
      )
    })

    output$download_heatmap_data <- downloadHandler(
      filename = function() {
        "heatmap_data.csv"
      },
      content = function(file) {
        write.csv(heatmap_data_download(), file)
      }
    )

  })
}

## To be copied in the UI
# mod_03_heatmap_ui("03_heatmap_ui_1")

## To be copied in the server
# mod_03_heatmap_server("03_heatmap_ui_1")
