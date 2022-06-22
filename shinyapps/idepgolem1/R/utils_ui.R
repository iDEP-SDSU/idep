#' utils_ui.R These are helper functions for ui code,
#' shiny code can go in this file
#'
#'
#' @section utils_ui.R functions:
#' \code{make_pull_down_menu} makes it easy to make pull down menu
#'
#'
#' @name utils_ui.R
NULL


#' make_pull_down_menu makes it easy to make pull down menu
#'
#'
#' @description
#'
#'
#' @param funs
#'
#'
#' @return
make_pull_down_menu <- function(funs) {
  return(setNames(1:length(funs), names(funs)))
}
