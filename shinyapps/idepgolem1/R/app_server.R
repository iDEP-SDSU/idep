#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # file size is 5MB by default. This changes it to 30MB
  #options(shiny.maxRequestSize=30*1024^2)

  # load static data files such as list of species, gmt files, etc
  # This could be moved to run_app as global variable, as in global.R
  # see https://github.com/ThinkR-open/golem/issues/6
  idep_data <- get_idep_data() 

  # Tab Variable to control reactivity
  tab <- reactive(input$navbar)

  load_data <- mod_01_load_data_server(
    id = "load_data",
    idep_data = idep_data,
    tab = tab
  )
  pre_process <- mod_02_pre_process_server(
    id = "pre_process",
    load_data = load_data,
    tab = tab
  )
  mod_03_clustering_server(
    id = "clustering",
    pre_process = pre_process,
    idep_data = idep_data,
    tab = tab
  )
  mod_04_pca_server(
    id = "pca",
    pre_process = pre_process,
    idep_data = idep_data
  )
  deg <- mod_05_deg_server(
    id = "deg",
    pre_process = pre_process, 
    idep_data = idep_data, 
    load_data = load_data
  )
  mod_06_pathway_server(
    id = "pathway",
    pre_process = pre_process,
    deg = deg,
    idep_data = idep_data,
    tab = tab
  )
  mod_07_genome_server(
    id = "genome",
    pre_process = pre_process,
    deg = deg,
    idep_data = idep_data
  )
  mod_08_bicluster_server(
    id = "bicluster",
    pre_process = pre_process,
    idep_data = idep_data,
    tab = tab
  )
  mod_09_network_server(
    id = "network",
    pre_process = pre_process,
    idep_data = idep_data,
    tab = tab
  )
  mod_10_doc_server(
    id = "doc",
    pre_process = pre_process,
    idep_data = idep_data,
    tab = tab
  )
}
