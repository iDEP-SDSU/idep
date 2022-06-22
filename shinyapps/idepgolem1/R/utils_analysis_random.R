#' utils_analysis_random.R Miscellaneous helper functions here,
#' we find a place for them later.These are functions don't manipulate data,
#' but there were used in data analysis functions for other purposes.
#'
#'
#' @section utils_analysis_random.R functions:
#' add later....
#'
#'
#' @name utils_analysis_random.R
NULL

#' hcluster_functions
#'
#'
#' @description
#'
#'
#' @return
#'
#'
hcluster_functions <- function() {
  # Functions for hierarchical clustering
  hclust_average <- function(x, method = "average", ...) { # average linkage
    hclust(x, method = method, ...)
  }
  hclust_ward_D <- function(x, method = "ward.D", ...) { # ward.D linkage # nolint
    hclust(x, method = method, ...)
  }
  hclust_ward_D2 <- function(x, method = "ward.D2", ...) { # ward.D2 linkage # nolint
    hclust(x, method = method, ...)
  }
  hclust_single <- function(x, method = "single", ...) { # single linkage
    hclust(x, method = method, ...)
  }
  hclust_mcquitty <- function(x, method = "mcquitty", ...) { # mcquitty linkage # nolint
    hclust(x, method = method, ...)
  }
  hclust_median <- function(x, method = "median", ...) { # median linkage
    hclust(x, method = method, ...)
  }
  hclust_centroid <- function(x, method = "centroid", ...) { # centroid linkage # nolint
    hclust(x, method = method, ...)
  }

  return(list(
    complete = hclust,
    average = hclust_average,
    ward_D = hclust_ward_D,
    ward_D2 = hclust_ward_D2,
    single = hclust_single,
    mcquitty = hclust_mcquitty,
    median = hclust_median,
    centroid = hclust_centroid
  ))
}


#' dist_functions
#'
#'
#' @description
#'
#'
#' @return
#'
#'
dist_functions <- function() {
  dist_pcc <- function(x, ...) {
    # distance function = 1-PCC (Pearson's correlation coefficient)
    as.dist(1 - cor(t(x), method = "pearson"))
  }

  dist_abs_pcc <- function(x, ...) {
    # distance function = 1-abs(PCC) (Pearson's correlation coefficient)
    as.dist(1 - abs(cor(t(x), method = "pearson")))
  }

  return(list(
    euclidean = dist,
    pearson_correlation = dist_pcc,
    absolute_pcc = dist_abs_pcc
  ))
}


#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param num_set DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
dynamic_range <- function(num_set) {
  # Given a set of numbers,
  # find the difference between 2nd largest and 2nd smallest
  sorted_set <- sort(num_set)
  if (length(num_set) >= 4) {
    k <- 2
  } else {
    k <- 1
  }
  return(sorted_set[length(num_set) - k + 1] - sorted_set[k])
}


#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param sample_names DESCRIPTION.
#' @param sample_info DESCRIPTION.
#'
#' @export
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
detect_groups <- function(sample_names, sample_info = NULL) {
  # sample_names are col names parsing samples by either the name
  # or using a data frame of sample infos.
  # Note that each row of the sample_info data frame represents a sample.
  sample_group <- NULL
  if (is.null(sample_info)) {
    # Remove all numbers from end
    # remove "_" from end
    # remove "_Rep" from end
    # remove "_rep" from end
    # remove "_REP" from end
    sample_group <- gsub(
      "[0-9]*$", "",
      sample_names
    )
    sample_group <- gsub("_$", "", sample_group)
    sample_group <- gsub("_Rep$", "", sample_group)
    sample_group <- gsub("_rep$", "", sample_group)
    sample_group <- gsub("_REP$", "", sample_group)
  } else {
    # the orders of samples might not be the same.
    # The total number of samples might also differ
    match_sample <- match(sample_names, row.names(sample_info))
    sample_info2 <- sample_info[match_sample, , drop = FALSE]
    if (ncol(sample_info2) == 1) {
      # if there's only one factor
      sample_group <- sample_info2[, 1]
    } else {
      # multiple columns/factors
      foo <- function(y) paste(y, collapse = "_")
      sample_group <- unlist(apply(
        X = sample_info2,
        MARGIN = 1,
        FUN = foo
      ))
      names(sample_group) <- row.names(sample_info2)
      if (min(table(sample_group)) == 1) { # no replicates?
        sample_group <- sample_info2[, 1]
      }
    }
  }
  return(as.character(sample_group))
}


# Clean up gene sets. Remove spaces and other control characters from gene names
#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param gene_set DESCRIPTION.
#'
#' @export
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
clean_gene_set <- function(gene_set) {
  # remove duplicate; upper case; remove special characters
  gene_set <- unique(toupper(gsub("\n| ", "", gene_set)))
  # genes should have at least two characters
  gene_set <- gene_set[which(nchar(gene_set) > 1)]
  return(gene_set)
}


#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param query_input DESCRIPTION.
#'
#' @export
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
clean_query <- function(query_input) {
  return(clean_gene_set(unlist(strsplit(
    x = toupper(query_input),
    split = "\t| |\n|\\,"
  ))))
}


# Read gene sets GMT file
# This functions cleans and converts to upper case
#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param file_path DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
read_gmt <- function(file_path) { # size restriction
  # Read in the first file
  gmt_data <- scan(file = file_path, what = "", sep = "\n")
  gmt_data <- gsub(
    pattern = " ",
    replacement = "",
    x = gmt_data
  )
  gmt_data <- toupper(gmt_data)

  #----Process the first file
  # Separate elements by one or more whitespace
  gmt_data <- strsplit(x = gmt_data, split = "\t")
  # Extract the first vector element and set it as the list element name
  extract <- function(x) x[[1]]
  names(gmt_data) <- sapply(X = gmt_data, FUN = extract)
  # Remove the first vector element from each list element
  extract2 <- function(x) x[-1]
  gmt_data <- lapply(X = gmt_data, FUN = extract2)
  # remove duplicated elements
  ## need to try to get this to lappy later
  for (i in 1:length(gmt_data)) {
    gmt_data[[i]] <- clean_gene_set(gmt_data[[i]])
  }
  # check the distribution of the size of gene lists sapply(y, length) hold a vector of sizes
  gene_size <- max(sapply(X = gmt_data, FUN = length))
  if (gene_size < 5) {
    cat("Warning! Gene sets have very small number of genes!\n
     Please double check format.")
  }
  # gene sets smaller than 1 is ignored!!!
  gmt_data <- gmt_data[which(sapply(X = gmt_data, FUN = length) > 1)]
  return(gmt_data)
}

### edit later
# This function convert gene set names
# x="GOBP_mmu_mgi_GO:0000183_chromatin_silencing_at_rDNA"
# chromatin silencing at rDNA
#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param x DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
proper <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))


#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param word_list DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
extract_word <- function(word_list) {
  words <- unlist(strsplit(word_list, "_"))
  if (length(words) <= 4) {
    return(gsub("_", " ", word_list))
  } else {
    words <- words[-c(1:4)]
    return(proper(paste(words, collapse = " ")))
  }
}


#' check_object_state an Utility function to simplify checking object states.
#'
#'
#' This function is a booling utility function,
#' which helps evaluates the state of objects,
#' and helps with sending messages from server to UI logic of shiny app.
#'
#'
#' @param check_exp An expression that should be evaluated
#'
#' @param true_message Message displayed on console and
#'  returned if \code{check_exp} is true
#'
#' @param false_message Optional message returned if \code{check_exp} is false,
#'  default to blank string
#'
#'
#' @return A list is returned, with the following elements:
#'  \code{bool} is either true or false depending on \code{check_exp} evaluation
#'
#'  \code{content} either message and depended on \code{check_exp} evaluation
#' @examples
#' check <- check_object_state(
#'   check_exp = (is.null(NULL)),
#'   true_message = as.data.frame("This is NULL")
#' )
#' # will eval to true and display message to console
#' if (check$bool) {
#'   return(check)
#' }
#'
#' check <- check_object_state(
#'   check_exp = (length(c(1, 2)) == 0),
#'   true_message = "this has 0 elements",
#'   false_message = "this doesn't have 0 elements"
#' )
#' # This will not return check and check$content will be false_message
#' if (check$bool) {
#'   return(check)
#' }
check_object_state <- function(
  check_exp,
  true_message,
  false_message = ""
) {
  if (check_exp) {
    message(true_message)
    return(list(
      bool = TRUE,
      content = true_message
    ))
  } else {
    return(list(
      bool = FALSE,
      content = false_message
    ))
  }
}

#' ggplot colors function
#'
#' This function will return the colors that
#' the ggplot2 package uses for plots.
#'
#' @param n Number of colors to return
#' 
#' @return Vector of hex color codes for a plot.
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Swap rowname IDs for data matrix
#'
#' This function uses the all_gene_names dataframe
#' to swap the current rownames of a data matrix with
#' the desired gene ID. For instance, if the rownames
#' were currently ensembl, this function is able to
#' switch them back to the original form.
#'
#' @param data_matrix Data matrix with ensembl or user gene ID rownames
#' @param all_gene_names Data frame of gene names
#' @param select_gene_id Desired ID type for rownames
#'   (User_ID, ensembl_ID, symbol)
#'
#' @export
#' @return Data matrix with changed rownames
rowname_id_swap <- function(
  data_matrix,
  all_gene_names,
  select_gene_id
) {
  if (select_gene_id == "User_ID" && ncol(all_gene_names) == 1) {
    return(data_matrix)
  } else if (select_gene_id == "User_ID") {
    data_matrix <- as.data.frame(data_matrix)
    data_matrix$order <- seq(1, nrow(data_matrix), 1)
    new_data <- merge(
      data_matrix,
      all_gene_names,
      by.x = "row.names",
      by.y = "ensembl_ID",
      all.x = T
    )
    rownames(new_data) <- new_data$User_ID
    nums <- unlist(lapply(new_data, is.numeric))
    new_data <- new_data[, nums]
    new_data <- new_data[order(new_data$order), ]
    new_data <- dplyr::select(new_data, -order)
    new_data <- as.matrix(new_data)
    return(new_data)
  } else if(select_gene_id == "ensembl_ID") {
    data_matrix <- as.data.frame(data_matrix)
    data_matrix$order <- seq(1, nrow(data_matrix), 1)
    new_data <- merge(
      data_matrix,
      all_gene_names,
      by.x = "row.names",
      by.y = "ensembl_ID",
      all.x = T
    )
    rownames(new_data) <- new_data$ensembl_ID
    nums <- unlist(lapply(new_data, is.numeric))
    new_data <- new_data[, nums]
    new_data <- new_data[order(new_data$order), ]
    new_data <- dplyr::select(new_data, -order)
    new_data <- as.matrix(new_data)
    return(new_data)
  } else if (select_gene_id == "symbol") {
    data_matrix <- as.data.frame(data_matrix)
    data_matrix$order <- seq(1, nrow(data_matrix), 1)
    new_data <- merge(
      data_matrix,
      all_gene_names,
      by.x = "row.names",
      by.y = "ensembl_ID",
      all.x = T
    )
    rownames(new_data) <- new_data$symbol
    nums <- unlist(lapply(new_data, is.numeric))
    new_data <- new_data[, nums]
    new_data <- new_data[order(new_data$order), ]
    new_data <- dplyr::select(new_data, -order)
    new_data <- as.matrix(new_data)
    return(new_data)
  }
}

#' Merge data from gene info and processed data
#'
#' This function takes in the gene info data and merges it
#' with the data that has gone through the processing fcn.
#' The returned data contains the gene names as well as the
#' ensembl id in the first two columns.
#'
#' @param all_gene_names All matched gene names from idep data
#' @param data Data matrix with rownames to merge with gene names
#' 
#' @export
#' @return Inputted data with all gene name information.
merge_data <- function(
  all_gene_names,
  data,
  merge_ID
) {
  isolate({
    if (dim(all_gene_names)[2] == 1) {
      new_data <- round(data, 2)
      new_data <- as.data.frame(new_data)
      new_data$User_id <- rownames(new_data)
      new_data <- dplyr::select(
        new_data,
        User_id,
        tidyselect::everything()
      )
      rownames(new_data) <- seq(1, nrow(new_data), 1)
      tmp <- apply(new_data[, 2:dim(new_data)[2]], 1, sd)
      new_data <- new_data[order(-tmp), ]

      return(new_data)
    } else if (dim(all_gene_names)[2] == 2) {
      new_data <- merge(
        all_gene_names,
        round(data, 2),
        by.x = merge_ID,
        by.y = "row.names",
        all.y = T
      )
      new_data <- dplyr::select(
        new_data,
        User_ID,
        ensembl_ID,
        tidyselect::everything()
      )
      rownames(new_data) <- seq(1, nrow(new_data), 1)
      tmp <- apply(new_data[, 3:dim(new_data)[2]], 1, sd)
      new_data <- new_data[order(-tmp), ]

      return(new_data)
    } else {
      new_data <- merge(
        all_gene_names,
        round(data, 2),
        by.x = merge_ID,
        by.y = "row.names",
        all.y = T
      )
      new_data <- dplyr::select(
        new_data,
        User_ID,
        ensembl_ID,
        symbol,
        tidyselect::everything()
      )
      rownames(new_data) <- seq(1, nrow(new_data), 1)
      tmp <- apply(new_data[, 4:dim(new_data)[2]], 1, sd)
      new_data <- new_data[order(-tmp), ]

      return(new_data)
    }
  })
}

# adding sample legends to heatmap; this is for the main heatmap
# https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
add_legend <- function(...) {
  opar <- par(
    fig = c(0, 1, 0, 1),
    oma = c(0, 0, 0, 0),
    mar = c(0, 0, 0, 6),
    new = TRUE
  )
  on.exit(par(opar))
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend(...)
}

#' Wrapping long text by adding \n 
#' "Mitotic DNA damage checkpoint"  --> "Mitotic DNA damage\ncheckpoint"
#' https://stackoverflow.com/questions/7367138/text-wrap-for-plot-titles
wrap_strings <- function(
  vector_of_strings,
  width = 30
  ) { 
  as.character(
    sapply(vector_of_strings, FUN = function(x) {
      paste(strwrap(x, width = width), collapse = "\n")
    })
  )
}

#' ENRICHMENT NETWORK FUNCTION
enrichment_network <- function(
  go_table,
  layout_button = 0,
  edge_cutoff = 5
){
	gene_lists <- lapply(go_table$Genes, function(x) unlist(strsplit(as.character(x), " ")))
	names(gene_lists) <- go_table$Pathways
	go_table$Direction <- gsub(" .*", "", go_table$Direction)

	g <- enrich_net(
    data = go_table,
    gene_set = gene_lists,
    node_id = "Pathways",
    num_char = 100, 
	  p_value = "adj.Pval",
    p_value_cutoff = 1,
    degree_cutoff = 0,
	  n = 200,
    group = go_table$Direction,
    vertex.label.cex = 1, 
    vertex.label.color = "black",
    show_legend = FALSE, 
    layout_button = layout_button,
    edge_cutoff = edge_cutoff
  )
}

#' numChar=100 maximum number of characters
#' n=200  maximum number of nodes
#' degree.cutoff = 0    Remove node if less connected
#' from PPInfer
enrich_net <-  function(
  data,
  gene_set,
  node_id,
  node_name = node_id,
  p_value, 
  n = 50,
  num_char = NULL,
  p_value_cutoff = 0.05,
  edge_cutoff = 0.05, 
  degree_cutoff = 0,
  edge_width = function(x) {5 * x^2},
  node_size = function(x) {2.5 * log10(x)},
  group = FALSE,
  group_color = c("green", "red"),
  group_shape = c("circle", "square"),
  legend_parameter = list("topright"),
  show_legend = TRUE,
  plotting = TRUE, 
  layout_button = 0,
  ...
) {
	set.seed(layout_button)
  data <- data.frame(data, group)
  colnames(data)[length(colnames(data))] <- "Group"
  data <- data[as.numeric(data[, "adj_p_val"]) < p_value_cutoff, ]
  data <- data[order(data[, "adj_p_val"]), ]
  n <- min(nrow(data), n)
  if (n == 0) {
    stop("no enriched term found...")
  }
  data <- data[1:n, ]
  index <- match(data[, node_id], names(gene_set))
  gene_sets_list <- list()
  for (i in 1:n) {
    gene_sets_list[[i]] <- gene_set[[index[i]]]
  }
  names(gene_sets_list) <- data[, node_name]
    
  if(is.null(num_char)) {
    num_char <- max(nchar(as.character(data[, node_name])))
  } else {
    if(length(unique(substr(data[, node_name], 1, num_char))) < nrow(data)) {
      num_char <- max(nchar(as.character(data[, node_name])))
      message("Note : numChar is too small.", "\n")
    }
  }
  data[, node_name] <- paste(
    substr(data[, node_name], 1, num_char), 
    ifelse(nchar(as.character(data[, node_name])) > num_char, "...", ""),
    sep = ""
  )
  w <- matrix(NA, nrow = n, ncol = n)

  for (i in 1:n) {
    for (j in i:n) {
      u <- unlist(gene_sets_list[i])
      v <- unlist(gene_sets_list[j])
      w[i, j] <- length(intersect(u, v)) / length(unique(c(u, v)))
    }
  }
  list_edges <- stack(data.frame(w))
  list_edges <- cbind(
    list_edges[, 1],
    rep(data[, node_name], n),
    rep(data[, node_name], each = n)
  )
  list_edges <- list_edges[list_edges[, 2] != list_edges[,3], ]
  list_edges <- list_edges[!is.na(list_edges[, 1]), ]
  g <- igraph::graph.data.frame(list_edges[, -1], directed = FALSE)
  igraph::E(g)$width <- edge_width(as.numeric(list_edges[, 1]))
  igraph::V(g)$size <- node_size(lengths(gene_sets_list))
  g <- igraph::delete.edges(g, igraph::E(g)[as.numeric(list_edges[, 1]) < edge_cutoff])
  index_deg <- igraph::degree(g) >= degree_cutoff
  g <- igraph::delete.vertices(g, igraph::V(g)[!index_deg])
  data <- data[index_deg, ]
  index <- index[index_deg]
  if(length(igraph::V(g)) == 0) {
    stop("no categories greater than degree_cutoff...")
  }
  n <- min(nrow(data), n)
  data <- data[1:n, ]
  group_level <- sort(unique(group))
  p_values <- log10(as.numeric(data[, "adj_p_val"]))
    
  for(i in 1:length(group_level)) {
    index <- data[, "Group"] == group_level[i]
    igraph::V(g)$shape[index] <- group_shape[i]
    group_p_values <- p_values[index]
    
    if(length(group_p_values) > 0) {
      if(max(group_p_values) == min(group_p_values)) {
        igraph::V(g)$color[index] <- grDevices::adjustcolor(
          group_color[i], 
          alpha.f = 0.5
        )
      } else {
        igraph::V(g)$color[index] <- sapply(
          1 - .9 * (group_p_values - min(group_p_values)) / 
          (max(group_p_values) - min(group_p_values)), 
          function(x) {
            grDevices::adjustcolor(group_color[i], alpha.f =  .1 + x)
          }
        )
      }
    }
  }
	if(plotting) {
    plot(g, , vertex.label.dist = 1.2, ...)
    if(show_legend) {
			legend.parameter$legend <- group.level
			legend.parameter$text.col <- group.color
			legend.parameter$bty <- "n"	
			do.call(legend, legend.parameter)
		}
  }
  
  return(g)
}

#' VIS NETWORK FUNCTION
vis_network_plot <- function(
  network_data
) {
  visNetwork::visNetwork(
    nodes = network_data$nodes,
    edges = network_data$edges,
    height = "700px",
    width = "700px"
  ) |> 
    visNetwork::visIgraphLayout(layout = "layout_with_fr") |>
    visNetwork::visNodes(
      color = list(
        border = "#000000",
        highlight = "#FF8000"
      ),
      font = list(
        color = "#000000",
        size = 20
      ),
      borderWidth = 1,
      shadow = list(
        enabled = TRUE,
        size = 10)
    ) |>
    visNetwork::visEdges(
      shadow = FALSE,
      color = list(
        color = "#A9A9A9",
        highlight = "#FFD700"
      )
    ) |> 
    visNetwork::visExport(
      type = "jpeg", 
      name = "export-network", 
      float = "left", 
      label = "Export as an image (only what's visible on the screen!)", 
      background = "white", 
      style = ""
    ) 
}	

#' EXTRAT UNDERSCORE
extract_under <- function(x) {
  words <- unlist(strsplit(x, "_"))
  if(length(words) <= 4) {
    return(gsub("_", " ", x))
  } else {
    words <- words[-c(1:4)]
    return(loose.rock::proper(paste(words, collapse = " ")))
  }
}

#' PATHVIEW FUNCTION
#' https://rdrr.io/bioc/pathview/src/R/colorpanel2.R
colorpanel2 <- function (
  n,
  low,
  mid,
  high
) {
  if (missing(mid) || missing(high)) {
    low <- grDevices::col2rgb(low)
    if (missing(high))
      high <- grDevices::col2rgb(mid)
    else high <- grDevices::col2rgb(high)
    red <- seq(low[1, 1], high[1, 1], length = n)/255
    green <- seq(low[3, 1], high[3, 1], length = n)/255
    blue <- seq(low[2, 1], high[2, 1], length = n)/255
  }
  else {
    isodd <- n%%2 == 1
    if (isodd) {
      n <- n + 1
    }
    low <- grDevices::col2rgb(low)
    mid <- grDevices::col2rgb(mid)
    high <- grDevices::col2rgb(high)
    lower <- floor(n/2)
    upper <- n - lower
    red <- c(seq(low[1, 1], mid[1, 1], length = lower),
             seq(mid[1, 1], high[1, 1], length = upper))/255
    green <- c(seq(low[3, 1], mid[3, 1], length = lower),
               seq(mid[3, 1], high[3, 1], length = upper))/255
    blue <- c(seq(low[2, 1], mid[2, 1], length = lower),
              seq(mid[2, 1], high[2, 1], length = upper))/255
    if (isodd) {
      red <- red[-(lower + 1)]
      green <- green[-(lower + 1)]
      blue <- blue[-(lower + 1)]
    }
  }
  grDevices::rgb(red, blue, green)
}

#' PATHVIEW SOURCE FUNCTION
#' https://rdrr.io/bioc/pathview/src/R/render.kegg.node.R
render.kegg.node <- function(
  plot.data,
  cols.ts,
  img,
  same.layer = TRUE,
  type = c("gene","compound")[1],
  text.col = "black",
  cex = 0.25
){
  width=ncol(img)
  height=nrow(img)
  nn=nrow(plot.data)
  pwids=plot.data$width
  if(!all(pwids==max(pwids))){
    message("Info: ", "some node width is different from others, and hence adjusted!")
    wc=table(pwids)
    pwids=plot.data$width=as.numeric(names(wc)[which.max(wc)])
  }

  if(type=="gene"){
  if(same.layer!=T){
    rect.out=sliced.shapes(plot.data$x+0.5, height-plot.data$y, plot.data$width/2-0.5, plot.data$height/2-0.25,  cols=cols.ts, draw.border=F, shape="rectangle")
    text(plot.data$x+0.5, height-plot.data$y, labels = as.character(plot.data$labels),
         cex = cex, col = text.col)
    return(invisible(1))
  } else{
    img2=img
    pidx=cbind(ceiling(plot.data$x-plot.data$width/2)+1,
      floor(plot.data$x+plot.data$width/2)+1,
      ceiling(plot.data$y-plot.data$height/2)+1,
      floor(plot.data$y+plot.data$height/2)+1)
    cols.ts=cbind(cols.ts)
    ns=ncol(cols.ts)
    brk.x= sapply(plot.data$width/2, function(wi) seq(-wi, wi, length = ns+1))
    for(k in 1:ns){
      col.rgb=col2rgb(cols.ts[,k])/255
      pxr=t(apply(pidx[,1:2], 1, function(x) x[1]:x[2]))-plot.data$x-1
      sel=pxr>=ceiling(brk.x[k,]) & pxr<=floor(brk.x[k+1,])
      for(i in 1:nn){
      sel.px=(pidx[i,1]:pidx[i,2])[sel[i,]]
      node.rgb=img[pidx[i,3]:pidx[i,4],sel.px, 1:3]
      node.rgb.sum=apply(node.rgb,c(1,2), sum)
      blk.ind=which(node.rgb.sum==0|node.rgb.sum==1,arr.ind=T)
      node.rgb=array(col.rgb[,i],dim(node.rgb)[3:1])
      node.rgb=aperm(node.rgb, 3:1)
      for(j in 1:3) node.rgb[cbind(blk.ind,j)]=0
      img2[pidx[i,3]:pidx[i,4],sel.px, 1:3]=node.rgb
    }
  }
    return(img2)
  }
} else if(type=="compound"){
  if(same.layer!=T){
    nc.cols=ncol(cbind(cols.ts))
    if(nc.cols>2){#block the background circle
      na.cols=rep("#FFFFFF", nrow(plot.data))
      cir.out=sliced.shapes(plot.data$x, height-plot.data$y, plot.data$width[1], plot.data$width[1], cols=na.cols, draw.border=F, shape="ellipse", lwd=0.2)
    }
    cir.out=sliced.shapes(plot.data$x, height-plot.data$y, plot.data$width[1], plot.data$width[1], cols=cols.ts, shape="ellipse", blwd=0.2)
    return(invisible(1))
  } else{
#    col.rgb=col2rgb(cols.ts)/255
    blk=c(0,0,0)
    img2=img
    w=ncol(img) #repeat
    h=nrow(img) #repeat
    cidx=rep(1:w, each=h)
    ridx=rep(1:h, w)
    pidx=lapply(1:nn, function(i){
      ii=which((cidx-plot.data$x[i])^2+(ridx-plot.data$y[i])^2<(plot.data$width[i])^2)
      imat=cbind(cbind(ridx, cidx)[rep(ii,each=3),],1:3)
      imat[,1:2]=imat[,1:2]+1
      ib=which(abs((cidx-plot.data$x[i])^2+(ridx-plot.data$y[i])^2-(plot.data$width[i])^2)<=8)
      ibmat=cbind(cbind(ridx, cidx)[rep(ib,each=3),],1:3)
      ibmat[,1:2]=ibmat[,1:2]+1
      return(list(fill=imat,border=ibmat))
    })

    cols.ts=cbind(cols.ts)
    ns=ncol(cols.ts)
    brk.x= sapply(plot.data$width, function(wi) seq(-wi, wi, length = ns+1))
    for(i in 1:nn){
      pxr=pidx[[i]]$fill[,2]-1-plot.data$x[i]
      col.rgb=col2rgb(cols.ts[i,])/255
      for(k in 1:ns){
        sel=pxr>=brk.x[k,i] & pxr<=brk.x[k+1,i]
        img2[pidx[[i]]$fill[sel,]]=col.rgb[,k]
      }
      img2[pidx[[i]]$border]=blk
    }
    return(img2)
  }
} else stop("unrecognized node type!")
}

#' PATHVIEW SOURCE FUNCTION
#' 
pathview.stamp <- function(
  x = NULL,
  y = NULL,
  position = "bottomright",
  graph.sizes,
  on.kegg = TRUE,
  cex = 1
){
  if(on.kegg)    labels ="Data on KEGG graph\nRendered by Pathview"
  else labels="-Data with KEGG pathway-\n-Rendered  by  Pathview-"
  if(is.null(x)| is.null(y)){
    x=graph.sizes[1]*.80
    y=graph.sizes[2]/40
    if(length(grep('left',position))==1)  x=graph.sizes[1]/40
    if(length(grep('top', position))==1)  y=graph.sizes[2]-y
  }
  text(x=x, y=y, labels=labels, adj=0, cex = cex, font=2)
}

#' Heatmap of the data
#' 
#' Create a ComplexHeatmap from a data matrix.
#' 
#' @param data A basic data matrix
#' @param heatmap_color_select Vector of colors to use for the fill
#'  in the heatmap
#' @export
#' @return A drawn ComplexHeatmap.
basic_heatmap <- function(
  data,
  heatmap_color_select
) {
  # Number of genes to show
	n_genes <- nrow(data)

  data <- as.matrix(data) - apply(data, 1, mean)
  cutoff <- median(unlist(data)) + 3 * sd(unlist(data)) 
	data[data > cutoff] <- cutoff
	cutoff <- median(unlist(data)) - 3 * sd(unlist(data)) 
	data[data < cutoff] <- cutoff
	
	data <- data[which(apply(data, 1, sd) > 0), ]

  # Color scale
  if(min(data) < 0) {
    col_fun <- circlize::colorRamp2(
      c(min(data), 0, max(data)),
      heatmap_color_select
    )
  } else {
    col_fun <- circlize::colorRamp2(
      c(min(data), median(data), max(data)),
      heatmap_color_select
    )
  }
  
	groups <- detect_groups(colnames(data))
	group_count <- length(unique(groups))
  groups_colors <- gg_color_hue(group_count)
  
  top_ann <- ComplexHeatmap::HeatmapAnnotation(
    Group = groups,
    col = list(
      Group = setNames(groups_colors, unique(groups))
    ),
    annotation_legend_param = list(
      Group = list(nrow = 1, title = NULL)
    ),
    show_annotation_name = list(Group = FALSE),
    show_legend = FALSE
  )

  heat <- ComplexHeatmap::Heatmap(
    data,
    name = "Expression",
    col = col_fun,
    cluster_rows = TRUE,
    clustering_method_rows = "average",
    clustering_distance_rows = function(x) as.dist(
      1 - cor(t(x), method = "pearson")
    ),
    cluster_columns = TRUE,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    top_annotation = top_ann,
    show_row_names = FALSE,
    show_column_names = FALSE,
    heatmap_legend_param = list(
      direction = "horizontal",
      legend_width = grid::unit(6, "cm"),
      title = "Color Key",
      title_position = "topcenter"
    )
  )
  
  return(
    heatmap = ComplexHeatmap::draw(
      heat,
      heatmap_legend_side = "bottom"
    )
  )
}

#' Heatmap of User brush selection
#' 
#' Create a heatmap from the brush selection of the main heatmap.
#' Used in iDEP to create an interactive heatmap and enable the
#' User to zoom in on areas they find interesting.
#' 
#' @param ht_brush Brush information from the User on the main
#'  heatmap
#' @param ht Main heatmap to create the sub-heatmap from
#' @param ht_pos_main Position information from the main heatmap
#'  to use for the sub-heatmap
#' @param heatmap_data Data matrix that is being plotted in the
#'  main heatmap
#' 
#' @export
#' @return A ComplexHeatmap object that will be inputted into the
#'  draw function in the server, the sub-heatmap data matrix, the
#'  group color mapping for the annotation, and the groups that
#'  the columns fall into.
basic_heat_sub <- function(
  ht_brush,
  ht,
  ht_pos_main,
  heatmap_data
) {
  lt <- InteractiveComplexHeatmap::getPositionFromBrush(ht_brush)
  pos1 <- lt[[1]]
  pos2 <- lt[[2]]

  pos <- InteractiveComplexHeatmap::selectArea(
    ht,
    mark = FALSE,
    pos1 = pos1,
    pos2 = pos2,
    verbose = FALSE,
    ht_pos = ht_pos_main
  )

  # Annotations ----------
    column_groups <- detect_groups(colnames(heatmap_data))
    groups_colors <- gg_color_hue(length(unique(column_groups)))
  
    top_ann <- ComplexHeatmap::HeatmapAnnotation(
      Group = column_groups,
      col = list(
        Group = setNames(
          groups_colors,
          unique(column_groups)
        )
      ),
      annotation_legend_param = list(
        Group = list(nrow = 1, title = NULL)
      ),
      show_annotation_name = list(Group = FALSE),
      show_legend = TRUE
    )
    
    group_col_return <- setNames(
      groups_colors,
      c(unique(column_groups))
    )
  # End annotation ---------

  column_index <- unlist(pos[1, "column_index"])
  row_index <- unlist(pos[1, "row_index"])
  top_ann <- top_ann[column_index]
  column_groups <- column_groups[column_index]
  m <- ht@ht_list[[1]]@matrix

  if (length(row_index) > 50) {
    show_rows <- FALSE
  } else {
    show_rows <- TRUE
  }

  submap_data <- m[row_index, column_index, drop = FALSE]

  ht_select <- ComplexHeatmap::Heatmap(
    submap_data,
    col = ht@ht_list[[1]]@matrix_color_mapping@col_fun,
    show_heatmap_legend = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = show_rows,
    top_annotation = top_ann,
    name = "heat_1"
  )

  return(list(
    ht_select = ht_select,
    submap_data = submap_data,
    group_colors = group_col_return,
    column_groups = column_groups
  ))
}


#' HTML code for sub-heatmap selected cell
#' 
#' Create HTML code for a cell of information on the cell of the
#' sub-heatmap that the User clicks on. The cell contains the
#' expression value, the sample, the gene, and the group.
#' 
#' @param click Information from the User clicking on a cell of
#'  the sub-heatmap
#' @param ht_sub The drawn sub-heatmap
#' @param ht_sub_obj The sub-heatmap ComplexHeatmap object
#' @param ht_pos_sub Position information for the sub-heatmap
#' @param sub_groups Vector of the groups that the samples
#'  belong to
#' @param group_colors color of the top annotation that
#'  is used for each group
#' @param data Sub data matrix that is plotted in the sub-heatmap
#' 
#' @return HTML code that will be used in the shiny UI to tell
#'  the user the information of the cell they selected.
heat_click_info <- function(
  click,
  ht_sub,
  ht_sub_obj,
  ht_pos_sub,
  sub_groups,
  group_colors,
  data
) {
  pos1 <- InteractiveComplexHeatmap::getPositionFromClick(click)
    
  pos <- InteractiveComplexHeatmap::selectPosition(
    ht_sub,
    mark = FALSE,
    pos = pos1,
    verbose = FALSE,
    ht_pos = ht_pos_sub
  )
  
  row_index <- pos[1, "row_index"]
  column_index <- pos[1, "column_index"]

  if (is.null(row_index)) {
    return("Select a cell in the heatmap.")
  }

  value <- data[row_index, column_index]
  col <- ComplexHeatmap::map_to_colors(ht_sub_obj@matrix_color_mapping, value)
  sample <- colnames(data)[column_index]
  gene <- rownames(data)[row_index]
  group_name <- sub_groups[column_index]
  group_col <- group_colors[[group_name]]

  # HTML for info table
  # Pulled from https://github.com/jokergoo/InteractiveComplexHeatmap/blob/master/R/shiny-server.R
  # Lines 1669:1678
  html <- GetoptLong::qq("
<div>
<pre>
Value: @{round(value, 2)} <span style='background-color:@{col};width=50px;'>    </span>
Sample: @{sample}
Gene: @{gene} 
Group: @{group_name} <span style='background-color:@{group_col};width=50px;'>    </span>   
</pre></div>"
)

 return(HTML(html))
}