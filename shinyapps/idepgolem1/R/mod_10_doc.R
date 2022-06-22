#' 10_doc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_10_doc_ui <- function(id){
  ns <- NS(id)
  tabPanel(
    "About",
    fluidPage(
      h3("Citation:"),
      p("Ge, Son & Yao, iDEP: an integrated web application for differential 
        expression and pathway analysis of RNA-Seq data, ", 
        a("BMC Bioinformatics 19:1-24, 2018.", 
           href="https://doi.org/10.1186/s12859-018-2486-6",
           target="_blank"
        )
      ),
      p("Consider citing other tools that form the foundation of iDEP, such as ", 
        a("ENSEMBL, ", 
           href="https:/doi.org/10.1093/nar/gkab1049",
           target="_blank"
        ),
        a("STRING-db,", 
           href="https://doi.org/10.1093/nar/gky1131",
           target="_blank"
        ),
        a("DESeq2", 
           href="https://doi.org/10.1186/s13059-014-0550-8",
           target="_blank"
        ),
        " and many others.",
        " If you use the KEGG diagram, please also cite ", 
        a("pathview, ", 
           href="https://doi.org/10.1093/bioinformatics/btt285",
           target="_blank"
        ),
        "and ",
        a("KEGG.", 
           href="https://doi.org/10.1093/nar/gkaa970",
           target="_blank"
        )
      )
    )
  )

}
    
#' 10_doc Server Functions
#'
#' @noRd 
mod_10_doc_server <- function(id, pre_process, idep_data, tab){
  moduleServer( id, function(input, output, session){
    ns <- session$ns 
  })
}
    
## To be copied in the UI
# mod_09_doc_ui("10_doc_ui_1")
    
## To be copied in the server
# mod_09_doc_server("10_doc_ui_1")
