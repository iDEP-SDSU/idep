#' fct_02_pre_process.R This file holds all of the main data analysis functions
#' associated with second tab of the iDEP website.
#'
#'
#' @section fct_02_pre_process.R functions:
#' \code{plot_genes}
#'
#'
#' @name fct_02_pre_process.R
NULL


#' @title Pre-Process the data
#'
#' @description This function takes in user defined values to
#' process the data for the EDA that occurs on the second page.
#' All the filtering and transformation that the app does will
#' occur on this page.
#'
#' @param data Data that has already gone through the convert_data fcn
#' @param missing_value Method to deal with missing data
#' @param data_file_format Type of data being examined
#' @param low_filter_fpkm Low count filter for the fpkm data
#' @param n_min_samples_fpkm Min samples for fpkm data
#' @param log_transform_fpkm Type of transformation for fpkm data
#' @param log_start_fpkm Value added to log transformation for fpkm
#' @param min_counts Low count filter for count data
#' @param n_min_samples_count Min sample for count data
#' @param counts_transform Type of transformation for counts data
#' @param counts_log_start Value added to log for counts data
#' @param no_fdr Fold changes only data with no p values
#' 
#' @export
#' @return A list containing the transformed data, the mean kurtosis,
#' the raw counts, a data type warning, the size of the original data,
#' and p-values.
pre_process <- function(
  data,
  missing_value,
  data_file_format,
  low_filter_fpkm,
  n_min_samples_fpkm,
  log_transform_fpkm,
  log_start_fpkm,
  min_counts,
  n_min_samples_count,
  counts_transform,
  counts_log_start,
  no_fdr
) {
  data_type_warning <- 0
  data_size_original <- dim(data)
  kurtosis_log <- 50

  # Sort by standard deviation -----------
  data <- data[order(-apply(
    data[, 1:dim(data)[2]],
    1,
    sd
  )), ]

  # Missng values in expression data ----------
  if (sum(is.na(data)) > 0) {
    if (missing_value == "geneMedian") {
      row_medians <- apply(data, 1, function(y) median(y, na.rm = T))
      for (i in 1:ncol(data)) {
        val_miss_row <- which(is.na(data[, i]))
        data[val_miss_row, i] <- row_medians[val_miss_row]
      }
    } else if (missing_value == "treatAsZero") {
      data[is.na(data)] <- 0
    } else if (missing_value == "geneMedianInGroup") {
      sample_groups <- detect_groups(colnames(data))
      for (group in unique(sample_groups)) {
        samples <- which(sample_groups == group)
        row_medians <- apply(
          data[, samples, drop = F],
          1,
          function(y) median(y, na.rm = T)
        )
        for (i in samples) {
          missing <- which(is.na(data[, i]))
          if (length(missing) > 0) {
            data[missing, i] <- row_medians[misssing]
          }
        }
      }
      if (sum(is.na(data)) > 0) {
        row_medians <- apply(
          data,
          1,
          function(y) median(y, na.rm = T)
        )
        for (i in 1:ncol(data)) {
          missing <- which(is.na(data[, i]))
          data[missing, i] <- row_medians[missing]
        }
      }
    }
  }
  # Compute kurtosis ---------
  mean_kurtosis <- mean(apply(data, 2, e1071::kurtosis), na.rm = TRUE)
  raw_counts <- NULL
  pvals <- NULL

  # Pre-processing for each file format ----------
  if (data_file_format == 2) {
    if (is.integer(data)) {
      data_type_warning <- 1
    }

    # Filters ----------
    # Not enough counts
    data <- data[which(apply(
      data,
      1,
      function(y) sum(y >= low_filter_fpkm)
    ) >= n_min_samples_fpkm), ]

    # Same levels in every entry
    data <- data[which(apply(
      data,
      1,
      function(y) max(y) - min(y)
    ) > 0), ]

    # Takes log if log is selected OR kurtosis is bigger than 50
    if (
      (log_transform_fpkm == TRUE) ||
        (mean_kurtosis > kurtosis_log)
    ) {
      data <- log(data + abs(log_start_fpkm), 2)
    }

    std_dev <- apply(data, 1, sd)
    data <- data[order(-std_dev), ]
  } else if (data_file_format == 1) {
    if (!is.integer(data) && mean_kurtosis < kurtosis_log) {
      data_type_warning <- -1
    }

    data <- round(data, 0)

    data <- data[which(apply(
      edgeR::cpm(edgeR::DGEList(counts = data)),
      1,
      function(y) sum(y >= min_counts)
    ) >= n_min_samples_count), ]

    raw_counts <- data

    # Construct DESeqExpression Object ----------
    tem <- rep("A", dim(data)[2])
    tem[1] <- "B"
    col_data <- cbind(colnames(data), tem)
    colnames(col_data) <- c("sample", "groups")
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = data,
      colData = col_data,
      design = ~groups
    )
    dds <- DESeq2::estimateSizeFactors(dds)

    # Counts Transformation ------------
    if (counts_transform == 3) {
      data <- DESeq2::rlog(dds, blind = TRUE)
      data <- SummarizedExperiment::assay(data)
    } else if (counts_transform == 2) {
      data <- DESeq2::vst(dds, blind = TRUE)
      data <- SummarizedExperiment::assay(data)
    } else {
      data <- log2(BiocGenerics::counts(
        dds,
        normalized = TRUE
      ) + counts_log_start)
    }
  } else if (data_file_format == 3) { # LFC and P-values
    n2 <- (ncol(data) %/% 2)
    raw_counts <- data
    if (!no_fdr) {
      pvals <- data[, 2 * (1:n2), drop = FALSE]
      data <- data[, 2 * (1:n2) - 1, drop = FALSE]
      if (ncol(data) == 1) {
        placeholder <- rep(1, dim(data)[1])
        pvals <- cbind(pvals, placeholder)
        zero_placeholder <- rep(0, dim(data)[1])
        data <- cbind(data, zero_placeholder)
      }
    }
  }
  data_size <- dim(data)

  validate(
    need(
      nrow(data) > 5 && ncol(data) >= 1,
      "Data file not recognized. Please double check."
    )
  )

  data <- data[order(-apply(
    data[, 1:dim(data)[2]],
    1,
    sd
  )), ]

  results <- list(
    data = as.matrix(data),
    mean_kurtosis = mean_kurtosis,
    raw_counts = raw_counts,
    data_type_warning = data_type_warning,
    data_size = c(data_size_original, data_size),
    p_vals = pvals
  )

  return(results)
}

#' Creates a barplot of the count data
#'
#' This function takes in the rount count data
#' and creates a formatted gg barplot that shows the number
#' of genes mapped to each sample in millions.
#'
#' @param counts_data Raw counts from gene expression data
#' @param sample_info Experiment file information for grouping
#'  samples
#'
#' @export
#' @return formatted ggbarplot
total_counts_ggplot <- function(
  counts_data,
  sample_info
) {
  counts <- counts_data
  memo <- ""

  if (ncol(counts) > 100) {
    part <- 1:100
    counts <- counts[, part]
    memo <- paste("(only showing 100 samples)")
  }
  groups <- as.factor(
    detect_groups(colnames(counts_data), sample_info)
  )

  if (nlevels(groups) <= 1 || nlevels(groups) > 20) {
    grouping <- NULL
  } else {
    grouping <- groups
  }

  if (ncol(counts) < 31) {
    x_axis_labels <- 16
  } else {
    x_axis_labels <- 12
  }
  plot_data <- data.frame(
    sample = as.factor(colnames(counts)),
    counts = colSums(counts) / 1e6,
    group = groups,
    grouping = grouping
  )

  plot <- ggplot2::ggplot(
    data = plot_data,
    ggplot2::aes(x = sample, y = counts, fill = grouping)
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right",
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(
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
      title = paste("Total Raw Read Counts (Millions)", memo),
      y = "Raw Counts (Millions)"
    )

  return(plot)
}

#' Scatterplot for EDA on processed data
#'
#' This function takes the data after it has been pre-processed and
#' creates a scatterplot of the counts for two samples that are
#' selected by the user.
#'
#' @param processed_data Data that has gone through the pre-processing
#' @param plot_xaxis Sample to plot on the x-axis
#' @param plot_yaxis Sample to plot on the y axis
#'
#' @export 
#' @return Returns a formatted gg scatterplot
eda_scatter <- function(
  processed_data,
  plot_xaxis,
  plot_yaxis
) {
  plot_data <- as.data.frame(processed_data)
  scatter <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = get(plot_xaxis),
      y = get(plot_yaxis)
    )
  ) +
    ggplot2::geom_point(size = 1) +
    ggplot2::labs(
      title = paste0(
        "Scatter for ", plot_xaxis, " and ", plot_yaxis,
        " transfromed expression"
      ),
      x = paste0(plot_xaxis),
      y = paste0(plot_yaxis)
    ) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "none",
      axis.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 16),
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      )
    )

  return(scatter)
}


#' Boxplot for processed data
#'
#' This function takes the processed data and creates
#' a boxplot of number of sequences mapped to each
#' tissue sample.
#'
#' @param processed_data Data that has gone through the pre-processing
#' @param sample_info Sample_info from the experiment file
#'
#' @export 
#' @return Formatted gg boxplot of the distribution of counts for each
#'  sample
eda_boxplot <- function(
  processed_data,
  sample_info
) {
  counts <- as.data.frame(processed_data)
  memo <- ""

  if (ncol(counts) > 40) {
    part <- 1:40
    counts <- counts[, part]
    memo <- paste(" (only showing 40 samples)")
  }
  groups <- as.factor(
    detect_groups(colnames(processed_data), sample_info)
  )

  if (nlevels(groups) <= 1 | nlevels(groups) > 20) {
    grouping <- NULL
  } else {
    grouping <- groups
  }
  if (ncol(counts) < 31) {
    x_axis_labels <- 16
  } else {
    x_axis_labels <- 12
  }

  longer_data <- tidyr::pivot_longer(
    data = counts,
    colnames(counts),
    names_to = "sample",
    values_to = "expression"
  )
  longer_data$groups <- rep(groups, nrow(counts))
  longer_data$grouping <- rep(grouping, nrow(counts))

  plot <- ggplot2::ggplot(
    data = longer_data,
    ggplot2::aes(x = sample, y = expression, fill = grouping)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right",
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(
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
      title = paste("Distribution of Transformed Data", memo),
      y = "Transformed Expression"
    )

  return(plot)
}

#' Density plot for the processed data
#'
#' This function takes in the processed data and sample info
#' and creates a density plot for the distribution of sequences
#' that are mapped to each sample.
#'
#' @param processed_data Data that has gone through the pre-processing
#' @param sample_info Sample_info from the experiment file
#'
#' @export 
#' @return Returns a formatted gg density plot
eda_density <- function(
  processed_data,
  sample_info
) {
  counts <- as.data.frame(processed_data)
  memo <- ""

  if (ncol(counts) > 40) {
    part <- 1:40
    counts <- counts[, part]
    memo <- paste(" (only showing 40 samples)")
  }
  groups <- as.factor(
    detect_groups(colnames(processed_data), sample_info)
  )

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

  longer_data <- tidyr::pivot_longer(
    data = counts,
    colnames(counts),
    names_to = "sample",
    values_to = "expression"
  )
  longer_data$groups <- rep(groups, nrow(counts))
  longer_data$group_fill <- rep(group_fill, nrow(counts))

  plot <- ggplot2::ggplot(
    data = longer_data,
    ggplot2::aes(x = expression, color = group_fill, group = sample)
  ) +
    ggplot2::geom_density(size = 1) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = legend,
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.text.x = ggplot2::element_text(
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
      title = paste("Density Plot of Transformed Data", memo),
      x = "Transformed Expression",
      y = "Density",
      color = "Sample"
    )

  return(plot)
}



#' Individual plotting function for genes
#'
#' Takes in the merged data and other data to provide
#' plots on individual gene names. Depening on the selections
#' this function will either return a barplot from the
#' grouped data or a lineplot from the individual sample.
#'
#' @param merged_data Data that has been merged with the gene info
#' @param sample_info Experiment file information for grouping samples
#' @param select_dene Gene(s) to be plotted
#' @param gene_plot_box T/F for individual sample plot or grouped data plot
#' @param use_sd T/F for standard error or standard deviation bars on bar plot
#' @param select_org Species the expression data is from
#'
#' @export 
#' @return A formatted ggplot. For gene_plot_box = TRUE the return will be a
#'  lineplot for the expression of each individual sample for the selected gene.
#'  If gene_plot_box = FALSE the return will be a barplot for the groups provided
#'  in the sample information.
individual_plots <- function(
  individual_data,
  sample_info,
  selected_gene,
  gene_plot_box,
  use_sd,
  lab_rotate
) {
  individual_data <- as.data.frame(individual_data)
  individual_data$symbol <- rownames(individual_data)

  plot_data <- individual_data |>
    dplyr::filter(symbol %in% selected_gene) |>
    tidyr::pivot_longer(!symbol, names_to = "sample", values_to = "value")
  if (ncol(plot_data) < 31) {
    x_axis_labels <- 14
  } else {
    x_axis_labels <- 10
  }

  if (gene_plot_box == TRUE) {
    ind_line <- ggplot2::ggplot(
      data = plot_data,
      ggplot2::aes(x = sample, y = value, group = symbol, color = symbol)
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 5, fill = "white") +
      ggplot2::labs(
        title = "Transformed Expression Level",
        y = "Transformed Expression"
      ) +
      ggplot2::coord_cartesian(ylim = c(0, max(plot_data$value))) +
      ggplot2::theme_light() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          color = "black",
          size = 16,
          face = "bold",
          hjust = .5
        ),
        axis.text.x = ggplot2::element_text(
          angle = as.numeric(lab_rotate),
          size = x_axis_labels,
          vjust = .5
        ),
        axis.text.y = ggplot2::element_text(size = 16),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(
          color = "black",
          size = 14
        ),
        legend.text = ggplot2::element_text(size = 12)
      )

    return(ind_line)
  } else if (gene_plot_box == FALSE) {
    plot_data$groups <- detect_groups(plot_data$sample, sample_info)

    summarized <- plot_data |>
      dplyr::group_by(groups, symbol) |>
      dplyr::summarise(Mean = mean(value), SD = sd(value), N = dplyr::n())

    summarized$SE <- summarized$SD / sqrt(summarized$N)

    gene_bar <- ggplot2::ggplot(
      summarized,
      ggplot2::aes(x = symbol, y = Mean, fill = groups)
    ) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
      ggplot2::labs(
        title = "Expression Level",
        y = "Grouped Transformed Expression"
      ) +
      ggplot2::theme_light() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          color = "black",
          size = 16,
          face = "bold",
          hjust = .5
        ),
        axis.text.x = ggplot2::element_text(
          angle = as.numeric(lab_rotate),
          size = x_axis_labels,
          vjust = .5
        ),
        axis.text.y = ggplot2::element_text(size = 16),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(
          color = "black",
          size = 14
        ),
        legend.text = ggplot2::element_text(size = 12)
      )

    if (use_sd == TRUE) {
      gene_bar <- gene_bar + ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = Mean - SD,
          ymax = Mean + SD
        ),
        width = 0.2,
        position = ggplot2::position_dodge(.9)
      )
    } else {
      gene_bar <- gene_bar + ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = Mean - SE,
          ymax = Mean + SE
        ),
        width = 0.2,
        position = ggplot2::position_dodge(.9)
      )
    }

    return(gene_bar)
  }
}

#' Data processing message
#'
#' Creates a message about the size of the counts
#' data and the amount of IDs that were converted.
#'
#' @param data_size Data size matrix from processing function
#' @param all_gene_names Data frame with all gene names
#' @param n_matched Count of matched IDs after processing
#'
#' @export 
#' @return Message about processed data
conversion_counts_message <- function(
  data_size,
  all_gene_names,
  n_matched
) {
  if (ncol(all_gene_names) == 1) {
    return(paste(
      data_size[1], "genes in", data_size[4], "samples.",
      data_size[3], " genes passed filter. Original gene IDs used."
    ))
  } else {
    return(paste(
      data_size[1], "genes in", data_size[4], "samples.",
      data_size[3], " genes passed filter, ", n_matched,
      " were converted to Ensembl gene IDs in our database.
      The remaining ", data_size[3] - n_matched, " genes were
      kept in the data using original IDs."
    ))
  }
}

#' Create message for sequencing depth bias
#'
#' This function creates a warning message for the
#' UI to present to the user regarding the sequencing
#' depth bias.
#'
#' @param raw_counts Raw counts data from the processing function
#' @param sample_info Experiment file information about each sample
#'
#' @export 
#' @return Message for the UI
counts_bias_message <- function(
  raw_counts,
  sample_info
) {
  total_counts <- colSums(raw_counts)
  groups <- as.factor(
    detect_groups(
      colnames(raw_counts),
      sample_info
    )
  )
  message <- NULL
  # ANOVA of total read counts vs sample groups parsed by sample name
  pval <- summary(aov(total_counts ~ groups))[[1]][["Pr(>F)"]][1]
  means <- aggregate(total_counts, by = list(groups), mean)
  max_min_ratio <- max(means[, 2]) / min(means[, 2])

  if (is.null(pval)) {
    message <- NULL
  } else if (pval < 0.05) {
    message <- paste(
      "Warning! Sequencing depth bias detected. Total read counts are
       significantly different among sample groups
       (p=", sprintf("%-3.2e", pval), ") based on ANOVA.
       Total read counts max/min =", round(max_min_ratio, 2)
    )
  }
  # ANOVA of total read counts vs factors in experiment design
  if (!is.null(sample_info)) {
    y <- sample_info
    for (j in 1:ncol(y)) {
      pval <- summary(aov(
        total_counts ~ as.factor(y[, j])
      ))[[1]][["Pr(>F)"]][1]

      if (pval < 0.01) {
        message <- paste(
          message, " Total read counts seem to be correlated with factor",
          colnames(y)[j], "(p=", sprintf("%-3.2e", pval), ").  "
        )
      }
    }
  }
  return(message)
}

#' Mean vs. Standard Deviation plot
#' 
#' Create a plot that shows the standard deviation as the
#' Y-axis across the mean of the counts data on the X-axis.
#' Option to make the X-axis the rank of the mean which
#' does a better job showing the spread of the data.
#' 
#' @param processed_data Data that has gone through the pre-processing
#' @param rank TRUE/FALSE whether to use the rank of the mean or not
#' @param heat_cols Heat color to use with black in the plot
#' 
#' @export 
#' @return A formatted ggplot hexplot of the mean and standard
#'  deviation of the processed data
mean_sd_plot <- function(
  processed_data,
  rank,
  heat_cols
) {
  table_data <- data.frame(
    "x_axis" = apply(
      processed_data,
      1,
      mean
    ),
    "y_axis" = apply(
      processed_data,
      1,
      sd
    )
  )

  if (rank) {
    table_data$x_axis <- rank(table_data$x_axis) 
  }
  low_col <- "black"
  high_col <- heat_cols[[1]]

  hex_plot <- ggplot2::ggplot(
      table_data,
      ggplot2::aes(x = x_axis, y = y_axis)
  ) +
  ggplot2::geom_hex() +
  ggplot2::geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs")
  ) +
  ggplot2::scale_fill_gradient2(
    mid = low_col,
    high = high_col
  ) +
  ggplot2::theme_light() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      color = "black",
      size = 16,
      face = "bold",
      hjust = .5
    ),
    axis.text.x = ggplot2::element_text(size = 14),
    axis.text.y = ggplot2::element_text(size = 14),
    axis.title.x = ggplot2::element_text(
      color = "black",
      size = 14
    ),
    axis.title.y = ggplot2::element_text(
      color = "black",
      size = 14
    ),
    legend.text = ggplot2::element_text(size = 12)
  ) +
  ggplot2::labs(
    title = "Mean vs. Standard Deviation",
    y = "Standard Deviation"
  )

  if (rank) {
    hex_plot <- hex_plot +
      ggplot2::labs(x = "Rank of Mean")
  } else {
    hex_plot <- hex_plot +
      ggplot2::labs(x = "Mean")
  }

  return(hex_plot)
}