#' fct_05_pca.R This file holds all of the main data analysis functions
#' associated with fifth tab of the iDEP website.
#'
#'
#' @section fct_05_pca.R functions:
#' \code{my_pgsea}
#'
#'
#' @name fct_05_pca.R
NULL



#' Principal Component Analysis
#'
#' Draw a PCA plot where user selects which PCs on axes
#'
#' @param data Data that has been through pre-processing
#' @param sample_info Matrix array with experiment info
#' @param PCAx PC on x axis
#' @param PCAy PC on y axis
#'
#' @export
#' @return Formatted PCA plot 
#' 
PCA_plot <- function(
  data,
  sample_info,
  PCAx = 1,
  PCAy = 2,
  selected_color = "Sample_Name",
  selected_shape = "Sample_Name"
) {

  #no design file
  if(is.null(selected_color)){
    selected_color <- "Sample_Name"
  }
  if(is.null(selected_shape)){
    selected_shape <- "Sample_Name"
  }
  counts <- data
  memo <- ""
  
  if (ncol(counts) > 100) {
    part <- 1:100
    counts <- counts[, part]
    memo <- paste("(only showing 100 samples)")
  }
  
  if (ncol(counts) < 31) {
    x_axis_labels <- 16
  } else {
    x_axis_labels <- 12
  }
  
  x <- data
  y <- sample_info
  pca.object <- prcomp(t(x))
  
  #5 pc's or number of columns if <5
  npc <- min(5, ncol(data))
  pcaData <- as.data.frame(pca.object$x[, 1:npc])

  groups <- detect_groups(sample_names = colnames(data), sample_info = sample_info)
  #Missing design clause
  if( is.null(sample_info)){
    pcaData <- cbind(pcaData, detect_groups(colnames(data), sample_info))
  } else {
    pcaData <- cbind(pcaData, detect_groups(colnames(x), y), sample_info)
  }
  #dim(pcaData)[2]
  colnames(pcaData)[npc + 1] <- "Sample_Name"
  if (nlevels(groups) <= 1 | nlevels(groups) > 20) {
    group_fill <- NULL
    legend <- "none"
  } else {
    group_fill <- groups
    legend <- "right"
  }
  
  if (ncol(counts) < 31) {
    x_axis_labels <- 16
  } else {
    x_axis_labels <- 12
  }
  
  P1 <- paste("PC", PCAx, sep = "")
  P2 <- paste("PC", PCAy, sep = "")
  
  
  # Set point & text size based on number of sample
  point_size <- 6
  if (ncol(x) >= 40) {
    point_size <- 3
  }
  
  plot_PCA <- ggplot2::ggplot(
    data = pcaData, 
    ggplot2::aes_string(
      x = paste0("PC",PCAx),
      y = paste0("PC",PCAy),
      color = selected_color,
      shape = selected_shape
    )
  )   +
  ggplot2::geom_point(size = point_size)+
  ggplot2::theme_light() +
  ggplot2::theme(
      legend.position = "right", # TODO no legend for large data
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = x_axis_labels
      ),
      axis.text.y = ggplot2::element_text(
        size = 16
      ),
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      )
    ) +
    ggplot2::labs(
      title = paste("Principal Component Analysis (PCA) ", memo),
      y = "Dimension 2",
      x = "Dimension 1"
    ) +
    ggplot2::guides(color=ggplot2::guide_legend(override.aes=list(shape=15)))
  # selected principal components
  PCAxy <- c(as.integer(PCAx), as.integer(PCAy))
  percentVar <- round(100 * summary(pca.object)$importance[2, PCAxy], 0)
  plot_PCA <- plot_PCA + ggplot2::xlab(paste0("PC", PCAx, ": ", percentVar[1], "% Variance"))
  plot_PCA <- plot_PCA + ggplot2::ylab(paste0("PC", PCAy, ": ", percentVar[2], "% Variance"))
  return(plot_PCA)
}

#' TSNE FUNCTION 
#'
#' Draw a t-sne plot where user selects which PCs on axes
#'
#'
#' @param data Data that has been through pre-processing
#' @param sample_info Matrix array with experiment info
#'
#' @export
#' @return Formatted T-sne plot
#'
#'
t_SNE_plot <- function(
  data,
  sample_info,
  selected_color,
  selected_shape
) {
  
  #no design file
  if(is.null(selected_color)){
    selected_color <- "Sample_Name"
  }
  if(is.null(selected_shape)){
    selected_shape <- "Sample_Name"
  }
  
  counts <- data
  memo <- ""
  
  if (ncol(counts) > 100) {
    part <- 1:100
    counts <- counts[, part]
    memo <- paste("(only showing 100 samples)")
  }

  if (ncol(counts) < 31) {
    x_axis_labels <- 16
  } else {
    x_axis_labels <- 12
  }

  
  x <- data
  y <- sample_info
  tsne <- Rtsne::Rtsne(t(x), dims = 2, perplexity = 1, verbose = FALSE, max_iter = 400)
  pcaData <- as.data.frame(tsne$Y)
  
  #Missing design clause
  if( is.null(sample_info)){
    pcaData <- cbind(pcaData, detect_groups(colnames(x), y))
  } else {
    pcaData <- cbind(pcaData, detect_groups(colnames(x), y), sample_info)
  }

  colnames(pcaData)[1:3] <- c("x1", "x2", "Sample_Name")
  
  # Set point size based on number of sample
  point_size <- 6
  if (ncol(x) >= 40) {
    point_size <- 3
  }
  
  #Generate plot
  plot_t_SNE <- ggplot2::ggplot(
    data = pcaData,
    ggplot2::aes_string(
      x = "x1",
      y = "x2",
      color = selected_color,
      shape = selected_shape)
    ) +
  ggplot2::geom_point(size = point_size) +
  ggplot2::theme_light() +
  ggplot2::theme(
    legend.position = "right",
    axis.title.y = ggplot2::element_text(
      color = "black",
      size = 14
    ),
    axis.title.x = ggplot2::element_text(
      color = "black",
      size = 14
    ),
    axis.text.x = ggplot2::element_text(
      angle = 90,
      size = x_axis_labels
    ),
    axis.text.y = ggplot2::element_text(
      size = 16
    ),
    plot.title = ggplot2::element_text(
      color = "black",
      size = 16,
      face = "bold",
      hjust = .5
    )
  ) +
    ggplot2::labs(
      title = paste("T-SNE ", memo),
      y = "Dimension 2",
      x = "Dimension 1"
    )

  return(plot_t_SNE)
}

#' MDS FUNCTION
#'
#' Draw a MDS plot
#'
#'
#' @param data Data that has been through pre-processing
#' @param sample_info Matrix array with experiment info
#'
#' @export
#' @return Formatted PCA plot
#'
MDS_plot <- function(
  data,
  sample_info,
  selected_shape = "Sample_Name",
  selected_color = "Sample_Name"
) {

  #no design file
  if(is.null(selected_color)){
    selected_color <- "Sample_Name"
  }
  if(is.null(selected_shape)){
    selected_shape <- "Sample_Name"
  }
  
  counts <- data
  memo <- ""
  
  if (ncol(counts) > 100) {
    part <- 1:100
    counts <- counts[, part]
    memo <- paste("(only showing 100 samples)")
  }
  
  if (ncol(counts) < 31) {
    x_axis_labels <- 16
  } else {
    x_axis_labels <- 12
  }
  
  x <- data
  y <- sample_info

  fit <- cmdscale(
    dist_functions()$pearson_correlation(t(x)),
    eig = T,
    k = 2
  )
  pcaData <- as.data.frame(fit$points[, 1:2])
  
  #Missing design clause
  if( is.null(sample_info)){
    pcaData <- cbind(pcaData, detect_groups(colnames(x), y))
  } else {
    pcaData <- cbind(pcaData, detect_groups(colnames(x), y), sample_info)
  }
  colnames(pcaData)[1:3] <- c("x1", "x2", "Sample_Name")

  # Set point & text size based on number of sample
  point_size <- 6

  if (ncol(x) >= 40) {
    point_size <- 3
    #text_size <- 16
  }

  
  p <- ggplot2::ggplot(
    data = pcaData,
    ggplot2::aes_string(
      x = "x1",
      y = "x2",
      color = selected_color,
      shape = selected_shape 
    )
  )
  p <- p + ggplot2::geom_point(size = point_size)+ 
  ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right",
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = x_axis_labels
      ),
      axis.text.y = ggplot2::element_text(
        size = 16
      ),
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      ) 
    ) +
    ggplot2::labs(
      title = paste("Multi-Dimensional Scaling (MDS) ", memo),
      y = "Dimension 2",
      x = "Dimension 1"
    )
  return(p)
}


#' Correlations Between Principle Components and Factors
#'
#' Desc
#'
#'
#' @param data Data that has been through pre-processing
#' @param sample_info Matrix array with experiment info
#'
#' @export
#' @return text with correlation
#'
pc_factor_correlation <- function(
  data,
  sample_info
) {
  x <- data
  y <- sample_info
  if (is.null(y)){
    y <- as.matrix(detect_groups(colnames(data)))
  }
  
  if(dim(y)[2] == 1)
  {
    return("No design file uploaded")
  }
  pca.object <- prcomp(t(x))
  
  #5 pc's or number of columns if <5
  npc <- min(5, ncol(data))
  
  pcaData <- as.data.frame(pca.object$x[, 1:npc])
  pvals <- matrix(1, nrow = npc, ncol = ncol(y))
  for (i in 1:npc) {
    for (j in 1:ncol(y)) {
      pvals[i, j] <- summary(
        aov(
          pcaData[, i] ~ as.factor(y[, j])
          )
        )[[1]][["Pr(>F)"]][1]
    }
  }
  pvals <- pvals * npc * ncol(y) # correcting for multiple testing
  pvals[pvals > 1] <- 1
  colnames(pvals) <- colnames(y)
  rownames(pvals) <- paste0("PC", 1:npc)
  a <- "Correlation between Principal Components (PCs) with factors: "
  nchar0 <- nchar(a)
  for (i in 1:npc) {
    j <- which.min(pvals[i, ])
    if (pvals[i, j] < 0.05) {
      a <- paste0(
        a, rownames(pvals)[i],
        " is correlated with ", colnames(pvals)[j],
        " (p=", sprintf("%-3.2e", pvals[i, j]), ")."
      )
    }
  }
  return(a)
}



#' Principal Component Analysis with PCAtools package
#'
#' Draw a PCA plot using PCAtools package
#'
#' @param data Data that has been through pre-processing
#' @param sample_info Matrix array with experiment info
#'
#' @export
#' @return Formatted PCA plot using PCAtools package
#' 
PCA_biplot <- function(
  data,
  sample_info,
  select_gene_id = "symbol",
  all_gene_names, # = pre_process$all_gene_names(),
  selected_x = "PC1",
  selected_y = "PC2",
  encircle = TRUE,
  encircleFill = TRUE,
  showLoadings = TRUE,
  pointlabs = TRUE,
  point_size = 4.0,
  ui_color = NULL,
  ui_shape = NULL
) {
  #missing design
  if(is.null(sample_info)) {
    meta_data <- as.data.frame(colnames(data))
    rownames(meta_data) <- colnames(data)
  } else {
    meta_data <- sample_info
  }

  #Swap rownames
  data <- rowname_id_swap(
    data_matrix = data,
    all_gene_names = all_gene_names,
    select_gene_id = select_gene_id
  )
  
  pca_obj <- PCAtools::pca(data, metadata = meta_data, removeVar = 0.1)
  
  if(pointlabs == TRUE){
    show_point_labels <- rownames(pca_obj$metadata)
  } else{
    show_point_labels <- NULL
  }

   PCAtools::biplot(
    pcaobj = pca_obj,
    x = selected_x,
    y = selected_y,
    colby = ui_color,
    shape = ui_shape,
    #colLegendTitle = 'Color?',
    encircle = encircle,
    encircleFill = encircleFill,
    showLoadings = showLoadings,
    lab = show_point_labels,
    legendPosition = 'right',
    legendLabSize = 16,
    legendIconSize = 8.0,
    pointSize = point_size,
    title = "Principal Component Scores"
    )

}

#' Principal Component Analysis with PCAtools package
#'
#' Draw a Scree plot with Horn's and Elbow suggestion for cutoffs using PCAtools package
#'
#' @param data Data that has been through pre-processing
#'
#' @export
#' @return Formatted Scree plot using PCAtools package
PCA_Scree <- function(
  processed_data 
) {

  suppressWarnings(
  pca_obj <- PCAtools::pca(mat = processed_data, removeVar = 0.1)
  )
  suppressWarnings(
  horn <- PCAtools::parallelPCA(processed_data)
  )
  
  suppressWarnings(
  elbow <- PCAtools::findElbowPoint(pca_obj$variance)
  )

  suppressWarnings(
  p <- plot(PCAtools::screeplot(
    pca_obj,
    vline = c(horn$n, elbow)) +
    ggplot2::geom_label(
      ggplot2::aes(
        x = horn$n + .1, 
        y = 60,
        label = 'Horn\'s', 
        vjust = .5,
        hjust = .5,
        size = 8)
      ) +
    ggplot2::geom_label(
      ggplot2::aes(x = elbow + .1,
                   y = 70,
                   label = 'Elbow',
                   vjust = .5,
                   hjust = .5,
                   size = 8))
  )

  )
  return(p)
}

#' Principal Component Analysis with PCAtools package
#'
#' Generates a plot showing correlations between Principal Components and design factors
#'
#' @param data Data that has been through pre-processing
#' @param sample_info Design Matrix
#' @return Formatted plot generated with PCAtools package
PCAtools_eigencorplot <- function(
  processed_data,
  sample_info
)
{
  #missing design
  if(is.null(sample_info)) {
    return("Upload Design file to see EigenCor plot.")
    #meta_data <- as.data.frame(colnames(processed_data))
    #colnames(meta_data)[1] <- "Sample_Name"
  } else {
    meta_data <- sample_info
    
    #Design Factors must be converted to numeric
    meta_data <- as.data.frame(meta_data)
    meta_data <- sapply(meta_data, function(x) as.numeric(factor(x)))
    meta_data <- as.data.frame(meta_data)
    
    #maintain rownames
    rownames(meta_data) <- rownames(sample_info)
    
    #create PCA object
    pca_obj <- PCAtools::pca(processed_data, metadata = meta_data, removeVar = 0.1)
    
    #plot
    p <- PCAtools::eigencorplot(pca_obj, metavars = colnames(meta_data))
    return(p)
  }
  }
