#' 10_doc 
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#'
#' @noRd
NULL


#' Density plot for the processed data
#'
#' This function takes in the processed data and sample info
#' and creates a density plot for the distribution of sequences
#' that are mapped to each sample.
#'
#' @param processed_data Data that has gone through the pre-processing
#' @param sample_info Sample_info from the experiment file
#'
#' @return Returns a formatted gg density plot
doc_hist <- function(
  processed_data,
  sampleID
) {
  hist(processed_data[, as.integer(sampleID)])
}
