#' fct_03_heatmap.R This file holds all of the main data analysis functions
#' associated with third tab of the iDEP website.
#'
#'
#' @section fct_03_heatmap.R functions:
#'
#'
#'
#' @name fct_03_heatmap.R
NULL


#' Density plot of data standard deviation
#'
#' Draw a denisty plot of the standard deviation in the
#' data. Add vertical red lines for a range of genes.
#'
#' @param data Data that has been through pre-processing
#' @param n_genes_max Upper limit of gene range
#'
#' @export
#' @return Formatted density plot of the standard deviation
#' distribution.
sd_density <- function(
  data,
  n_genes_max 
) {
  n_genes_min <- 5
  sds <- apply(data[, 1:dim(data)[2]], 1, sd)
  max_sd <- mean(sds) + 4 * sd(sds)
  sds[sds > max_sd] <- max_sd

  if (n_genes_max - n_genes_min == 0) {
    n_genes_max <- n_genes_min + 10
  }

  if (n_genes_max > length(sds)) {
    n_genes_max <- length(sds)
  }

  if (n_genes_min > length(sds)) {
    n_genes_min <- length(sds) - 10
  }

  sds <- as.data.frame(sds)

  plot <- ggplot2::ggplot(sds, ggplot2::aes(x = sds)) +
    ggplot2::geom_density(color = "darkblue", fill = "lightblue") +
    ggplot2::labs(
      title = "Standard Deviations of All Genes",
      y = "Density",
      x = "Standard Deviation"
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
        size = 14
      ),
      axis.text.y = ggplot2::element_text(size = 16),
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      legend.text = ggplot2::element_text(size = 12)
    )

  if (n_genes_max == 10 && n_genes_min == 0) {
    return(plot)
  } else if (n_genes_min == 0) {
    cutoff <- sort(sds$sds, decreasing = TRUE)[n_genes_max]

    plot <- plot +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = cutoff),
        color = "red",
        linetype = "dashed",
        size = 1
      ) +
      ggplot2::annotate(
        "text",
        x = cutoff + 0.4 * sd(sds[, 1]),
        y = 1,
        colour = "red",
        label = paste0("Upper: ", n_genes_max)
      )

    return(plot)
  } else {
    cutoff_max <- sort(sds$sds, decreasing = TRUE)[n_genes_max]
    cutoff_min <- sort(sds$sds, decreasing = TRUE)[n_genes_min]

    plot <- plot +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = cutoff_max),
        color = "red",
        linetype = "dashed",
        size = 1
      ) +
      ggplot2::annotate(
        "text",
        x = cutoff_max - 0.4 * sd(sds[, 1]),
        y = 1,
        colour = "red",
        label = paste0("Upper: ", n_genes_max)
      ) +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = cutoff_min),
        color = "red",
        linetype = "dashed",
        size = 1
      ) +
      ggplot2::annotate(
        "text",
        x = cutoff_min + 0.4 * sd(sds[, 1]),
        y = 1,
        colour = "red",
        label = paste0("Lower: ", n_genes_min)
      )

    return(plot)
  }
}


#' Heatmap data process
#'
#' This function prepares the data from pre-processing
#' to be displayed in a heatmap. It takes in limits for
#' what genes to subset, what centering and standardizing
#' to perform, and what gene ID label to use.
#'
#' @param data Processed data matrix
#' @param n_genes_max Row number upper limit to display in heatmap
#' @param gene_centering TRUE/FALSE subtract mean from gene rows
#' @param gene_normalize TRUE/FALSE divide by SD in gene rows
#' @param sample_centering TRUE/FALSE subtract mean from sample columns
#' @param sample_normalize TRUE/FALSE divide by SD in sample columns
#' @param all_gene_names Data frame of gene names
#' @param select_gene_id Desired ID type for heatmap labels
#'  (User_ID, ensembl_ID, symbol)
#'
#' @export
#' @return Subsetted data matrix ([n_genes_min:n_genes_max, ]) with
#'   gene IDs as the select_gene_id
process_heatmap_data <- function(
  data, 
  n_genes_max,
  gene_centering,
  gene_normalize,
  sample_centering,
  sample_normalize,
  all_gene_names,
  select_gene_id
) {
  
  data <- rowname_id_swap(
    data_matrix = data,
    all_gene_names = all_gene_names,
    select_gene_id = select_gene_id
  )

  data <- data[order(-apply(
    data[, 1:dim(data)[2]],
    1,
    sd
  )), ]

  if (n_genes_max > nrow(data)) {
    n_genes_max <- nrow(data)
  } else if (n_genes_max < 10) {
    n_genes_max <- 10
  }

  # Center by genes, cutoff really large values ---------
  if (gene_centering) {
    data <-
      data[1:n_genes_max, ] -
      apply(data[1:n_genes_max, ], 1, mean)
  }

  # Standardize by gene ---------
  if (gene_normalize) {
    data <- data / apply(data, 1, sd)
  }

  # Row centering and normalize ----------
  data <- scale(
    data,
    center = sample_centering,
    scale = sample_normalize
  )

  if (gene_centering) {
    return(round(data, 3))
  } else {
    data <- data[1:n_genes_max, ]
  }

  return(round(data, 3))
}

#' Draw a heatmap of processed data
#'
#' Uses the package ComplexHeatmaps to draw a heatmap of the
#' processed data that has been prepared for the heatmap. The
#' returned heatmap visualizes the processed expression of the
#' gene range provided in process_heatmap_data.
#'
#' @param data Data returned from process_heatmap_data
#' @param cluster_meth Type of clustering to use (Hierarchical/k-Means)
#' @param heatmap_cutoff Z score max to filter data
#' @param sample_info Experiment design information from load data
#' @param select_factors_heatmap Factor group for annotation legend
#' @param dist_funs List of distance functions to use in heatmap
#' @param dist_function The selected distance function to use
#' @param hclust_function Type of clustering to perform
#' @param no_sample_clustering TRUE/FALSE Specify whehter to cluster columns
#' @param heatmap_color_select Vector of colors for heatmap scale
#' @param row_dend TRUE/FALSE Hide row dendogram,
#' @param k_clusters Number of clusters to use for k-means
#' @param re_run Re-run k-means with a different seed
#'
#' @export
#' @return Heatmap of the processed data.
heatmap_main <- function(
  data,
  cluster_meth,
  heatmap_cutoff,
  sample_info,
  select_factors_heatmap,
  dist_funs,
  dist_function,
  hclust_function,
  no_sample_clustering,
  heatmap_color_select,
  row_dend,
  k_clusters,
  re_run
) {
  # Filter with max z-score
  cutoff <- median(unlist(data)) + heatmap_cutoff * sd(unlist(data))
  data[data > cutoff] <- cutoff
  cutoff <- median(unlist(data)) - heatmap_cutoff * sd(unlist(data))
  data[data < cutoff] <- cutoff

  # Color scale
  if (min(data) < 0) {
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

  # Annotation for groups
  groups <- detect_groups(colnames(data))

  if (!is.null(sample_info) && !is.null(select_factors_heatmap)) {
    if (select_factors_heatmap == "Sample_Name") {
      groups <- detect_groups(colnames(data))
    } else {
      ix <- match(select_factors_heatmap, colnames(sample_info))
      groups <- sample_info[, ix]
    }
  }

  groups_colors <- gg_color_hue(length(unique(groups)))

  if (length(groups) < 30) {
    show_group_leg <- TRUE
  } else {
    show_group_leg <- FALSE
  }
  
  heat_ann <- ComplexHeatmap::HeatmapAnnotation(
    Group = groups,
    col = list(Group = setNames(groups_colors, unique(groups))),
    annotation_legend_param = list(
      Group = list(nrow = 1, title = NULL)
    ),
    show_annotation_name = list(Group = FALSE),
    show_legend = FALSE
  )

  # Different heatmaps for hierarchical and k-means
  if (cluster_meth == 1) {
    heat <- ComplexHeatmap::Heatmap(
      data,
      name = "Expression",
      col = col_fun,
      clustering_method_rows = hclust_function,
      clustering_method_columns = hclust_function,
      clustering_distance_rows = function(x) {
        dist_funs[[as.numeric(dist_function)]](x)
      },
      clustering_distance_columns = function(x) {
        dist_funs[[as.numeric(dist_function)]](x)
      },
      cluster_rows = TRUE,
      cluster_columns = !(no_sample_clustering),
      show_column_dend = TRUE,
      show_row_dend = row_dend,
      row_dend_side = "left",
      row_dend_width = grid::unit(1, "cm"),
      top_annotation = heat_ann,
      show_row_names = FALSE,
      show_column_names = FALSE,
      heatmap_legend_param = list(
        direction = "horizontal",
        legend_width = grid::unit(6, "cm"),
        title = "Color Key",
        title_position = "topcenter"
      )
    )
  } else if (cluster_meth == 2) {
    set.seed(re_run)
    if(k_clusters > 10){
      row_title = 8
    } else {
      row_title = 10
    }

    heat <- ComplexHeatmap::Heatmap(
      data,
      name = "Expression",
      col = col_fun,
      row_km = k_clusters,
      cluster_rows = TRUE,
      cluster_columns = !(no_sample_clustering),
      show_column_dend = TRUE,
      show_row_dend = row_dend,
      row_dend_side = "left",
      row_dend_width = grid::unit(1, "cm"),
      top_annotation = heat_ann,
      show_row_names = FALSE,
      show_column_names = FALSE,
      heatmap_legend_param = list(
        direction = "horizontal",
        legend_width = grid::unit(6, "cm"),
        title = "Color Key",
        title_position = "topcenter"
      ),
      row_title_gp = grid::gpar(fontsize = row_title)
    )
  }

  return(
    ComplexHeatmap::draw(
      heat,
      heatmap_legend_side = "bottom"
    )
  )  
}

#' Draw a dendogram of data samples
#' 
#' Create a clustered tree of the samples in the dataset.
#' 
#' @param tree_data Data that has been through pre-processing
#' @param gene_centering TRUE/FALSE subtract mean from gene rows
#' @param gene_normalize TRUE/FALSE divide by SD in gene rows
#' @param sample_centering TRUE/FALSE subtract mean from sample columns
#' @param sample_normalize TRUE/FALSE divide by SD in sample columns
#' @param hclust_funs Clustering functions defined in idep
#' @param hclust_function String of chosen clustering method
#' @param dist_funs Distance functions defined in idep
#' @param dist_function Selected distance function
#' 
#' @export
#' @return Dendogram plot of dataset samples
draw_sample_tree <- function(
  tree_data,
  gene_centering,
  gene_normalize,
  sample_centering,
  sample_normalize,
  hclust_funs,
  hclust_function,
  dist_funs,
  dist_function
) {
  max_gene <- apply(tree_data, 1, max)
  # Remove bottom 25% lowly expressed genes, which inflate the PPC
  tree_data <- tree_data[which(max_gene > quantile(max_gene)[1]), ]
  # Center by gene
  if (gene_centering) {
    tree_data <- tree_data - apply(tree_data, 1, mean)
  }
  # Normalize by gene
  if (gene_normalize) {
    tree_data <- tree_data / apply(tree_data, 1, sd)
  }
  # Center and normalize by sample
  tree_data <- scale(
    tree_data,
    center = sample_centering,
    scale = sample_normalize
  )
  
  par(mar = c(5.1,4.1,4.1,20))
  
  plot(
    stats::as.dendrogram(
      hclust_funs[[hclust_function]]
      (dist_funs[[as.numeric(dist_function)]](t(tree_data)))
    ),
    xlab = "",
    ylab = paste(
      names(dist_funs)[as.numeric(dist_function)], "(",
      hclust_function, "linkage", ")"
    ),
    type = "rectangle", 
    #leaflab = "textlike"
    horiz = TRUE
  )

}

#' Draw an elbow plot for k-cluster selection
#' 
#' This function takes in the processed heatmap data and
#' creates an elbow plot to guide the selection of the
#' number of clusters to create
#' 
#' @param heatmap_data Processed heatmap data
#' 
#' @export
#' @return Formatted elbow plot for the data
k_means_elbow <- function(
  heatmap_data
) {
  k.max = 20
  
  validate(
    need(nrow(heatmap_data) > k.max, 
    message = paste(paste("To create the elbow plot, please select at least", k.max  + 1), "genes."))
  )
  

  factoextra::fviz_nbclust(
    heatmap_data,
    kmeans,
    method = "wss",
    k.max = k.max
  ) +
  ggplot2::theme_light() +
  ggplot2::theme(
    legend.position = "none",
    axis.title.x = ggplot2::element_text(
      color = "black",
      size = 14
    ),
    axis.title.y = ggplot2::element_text(
      color = "black",
      size = 14
    ),
    axis.text.x = ggplot2::element_text(
      size = 12
    ),
    axis.text.y = ggplot2::element_text(
      size = 12
    ),
    plot.title = ggplot2::element_text(
      color = "black",
      size = 16,
      face = "bold",
      hjust = .5
    )
  ) +
  ggplot2::labs(
    title = "k-Means Elbow Plot",
    y = "Within Sum of Square",
    x = "Clusters"
  )
}


#' Create annotation for shiny subheatmap
#' 
#' Use the heatmap data to make an annotation for the
#' submap that will also show the legend
#' 
#' @param data Heatmap data
#' @param sample_info Experiment design file from load data
#' @param select_factors_heatmap Factor to group by in the samples
#' 
#' @export
#' @return A list containing a ComplexHeatmap annotation object,
#'  a ComplexHeatmap legend, list of groups, and list of group colors.
sub_heat_ann <- function(
  data,
  sample_info,
  select_factors_heatmap
) {
  groups <- detect_groups(colnames(data))

  if (!is.null(sample_info) && !is.null(select_factors_heatmap)) {
    if (select_factors_heatmap == "Sample_Name") {
      groups <- detect_groups(colnames(data))
    } else {
      ix <- match(select_factors_heatmap, colnames(sample_info))
      groups <- sample_info[, ix]
    }
  }

  groups_colors <- gg_color_hue(length(unique(groups)))
  group_colors <- setNames(groups_colors, unique(groups))

  if (length(unique(groups)) < 10) {
    lgd <- ComplexHeatmap::Legend(
      at = unique(groups),
      legend_gp = grid::gpar(fill = groups_colors),
      nrow = 1
    )
  } else {
    lgd <- NULL
  }
  
  heat_sub_ann <- ComplexHeatmap::HeatmapAnnotation(
    Group = groups,
    col = list(Group = setNames(groups_colors, unique(groups))),
    show_annotation_name = list(Group = FALSE),
    show_legend = FALSE
  )

  return(list(
    heat_sub_ann = heat_sub_ann,
    lgd = lgd,
    groups = groups,
    group_colors = group_colors
  ))
}

#' Interactive click text for subheatmap
#' 
#' Create a text output to tell the user the cell
#' information for their click.
#' 
#' @param click Click input from subheatmap
#' @param ht_sub Drawn subheatmap
#' @param ht_sub_obj Heatmap object with mapping info
#' @param ht_pos_sub Position information from submap
#' @param sub_groups Vector of group labels from submap
#' @param group_colors Colors for the group annotation
#' @param cluster_meth Type of clustering being performed
#' @param click_data Data matrix to get the data value from
#' 
#' @export
#' @return HTML code to produce a table with information
#'  about the selected cell.
cluster_heat_click_info <- function(
  click,
  ht_sub,
  ht_sub_obj,
  ht_pos_sub,
  sub_groups,
  group_colors,
  cluster_meth,
  click_data
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

  if (cluster_meth == 1) {

    value <- click_data[row_index, column_index]
    col <- ComplexHeatmap::map_to_colors(ht_sub_obj@matrix_color_mapping, value)
    sample <- colnames(click_data)[column_index]
    gene <- rownames(click_data)[row_index]

  } else if (cluster_meth == 2) {

    clust <- pos[1, "heatmap"]
    sub_click_data <- click_data[[clust]]
    value <- sub_click_data[row_index, column_index]
    col <- ComplexHeatmap::map_to_colors(
      ht_sub_obj@ht_list$heat_1@matrix_color_mapping,
      value
    )
    sample <- colnames(sub_click_data)[column_index]
    gene <- rownames(sub_click_data)[row_index]

  }
  
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
  HTML(html)
}

#' Draw sub heatmap from brush input
#' 
#' Use the brush input from the main heatmap to
#' create a larger subheatmap. 
#' 
#' @param ht_brush Brush input from the main heatmap
#' @param ht Main heatmap object
#' @param ht_pos_main Position of brush on main heatmap
#' @param heatmap_data Data for the heatmap
#' @param sample_info Experiment design file
#' @param select_factors_heatmap Group design to label by
#' @param cluster_meth Type of clustering being performed
#' 
#' @export
#' @return A list containing a Heatmap from the brush selection
#'  of the main heatmap, the submap data matrix, the groups for
#'  the submap, the submap legend, and data for the click info.
heat_sub <- function(
  ht_brush,
  ht,
  ht_pos_main,
  heatmap_data,
  sample_info,
  select_factors_heatmap,
  cluster_meth
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

  column_index <- unlist(pos[1, "column_index"])

  # Annotation, groups, and legend
  sub_heat <- sub_heat_ann(
    data = heatmap_data,
    sample_info = sample_info,
    select_factors_heatmap = select_factors_heatmap
  )
  sub_ann <- sub_heat$heat_sub_ann[column_index]
  sub_groups <- sub_heat$groups[column_index]
  lgd <- sub_heat$lgd
  group_colors <- sub_heat$group_colors
  
  if (cluster_meth == 1) {
    row_index <- unlist(pos[1, "row_index"])
    m <- ht@ht_list[[1]]@matrix
    if (length(row_index) > 50) {
      show_rows <- FALSE
    } else {
      show_rows <- TRUE
    }
    submap_data <- m[row_index, column_index, drop = FALSE]
    click_data <- submap_data

    ht_select <- ComplexHeatmap::Heatmap(
      m[row_index, column_index, drop = FALSE],
      col = ht@ht_list[[1]]@matrix_color_mapping@col_fun,
      show_heatmap_legend = FALSE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = show_rows,
      top_annotation = sub_ann,
      name = "heat_1"
    )
  } else if (cluster_meth == 2) {
    sub_heats = c()
    all_rows = c()
    click_data = c()

    for (i in 1:nrow(pos)) {
      all_rows <- c(all_rows, unlist(pos[i, "row_index"]))
    }

    if (length(all_rows) > 50) {
        show_rows <- FALSE
      } else {
        show_rows <- TRUE
    }
    m <- ht@ht_list[[1]]@matrix
    submap_data <- m[all_rows, column_index, drop = FALSE]
    
    for (i in 1:nrow(pos)) {
      row_index <- unlist(pos[i, "row_index"])
      
      click_data[[paste0("heat_", i)]] <- m[row_index, column_index, drop = FALSE]

      sub_heats[[i]] <- ComplexHeatmap::Heatmap(
        m[row_index, column_index, drop = FALSE],
        col = ht@ht_list[[1]]@matrix_color_mapping@col_fun,
        show_heatmap_legend = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = show_rows,
        name = paste0("heat_", i)
      )
      if (i == 1) {
        sub_heats[[i]] <- ComplexHeatmap::add_heatmap(
          sub_ann,
          sub_heats[[i]],
          direction = "vertical"
        )
      } else if (i >= 2) {
        sub_heats[[i]] <- ComplexHeatmap::add_heatmap(
          sub_heats[[i-1]],
          sub_heats[[i]],
          direction = "vertical"
        )
      }
      ht_select <- sub_heats[[i]]
    }
  }

  return(list(
    ht_select = ht_select,
    submap_data = submap_data,
    sub_groups = sub_groups,
    lgd = lgd,
    group_colors = group_colors,
    click_data = click_data
  ))
}

#' Create a correlation plot from heatmap data
#' 
#' Creates a correlation matrix heatmap from the
#' heatmap data to demonstrate the correlation
#' between samples.
#' 
#' @param data Heatmap data
#' @param label_pcc Label with correlation coefficient when TRUE
#' @param heat_cols Heat colors to use in the correlation matrix
#' @param text_col Color to make the text labels in the plot
#' 
#' @export
#' @return ggplot heatmap of correlation matrix
cor_plot <- function(
  data,
  label_pcc,
  heat_cols,
  text_col
) {
  # remove bottom 25% lowly expressed genes, which inflate the PPC 
	max_gene <- apply(data, 1, max)
	data <- data[which(max_gene > quantile(max_gene)[1] ), ]
  low_col <- heat_cols[[1]]
  mid_col <- heat_cols[[2]]
  high_col <- heat_cols[[3]]
		
	melted_cormat <- reshape2::melt(round(cor(data), 2), na.rm = TRUE)

	ggheatmap <- ggplot2::ggplot(
    melted_cormat,
    ggplot2::aes(Var2, Var1, fill = value)
  ) +
	ggplot2::geom_tile(color = text_col) +
	ggplot2::scale_fill_gradient2(
    low = low_col,
    high = high_col, 
    mid = mid_col, 
		space = "Lab",
    limit = c(
      min(melted_cormat[, 3]),
      max(melted_cormat[, 3])
    ),
    midpoint = median(melted_cormat[, 3]),
		name = "Pearson's \nCorrelation"
  ) +
	ggplot2::theme_minimal()+ # minimal theme
	ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      vjust = 1,
      size = 14,
      hjust = 1
    ),
    axis.text.y = ggplot2::element_text(size = 14 )
  ) +
	ggplot2::coord_fixed()

	if(label_pcc && ncol(data)<20) {
    ggheatmap <- ggheatmap +
    ggplot2::geom_text(
      ggplot2::aes(Var2, Var1, label = value),
      color = text_col,
      size = 4
    )
  }	
		
  ggheatmap + ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
		axis.title.y = ggplot2::element_blank(),
		panel.grid.major = ggplot2::element_blank(),
		panel.border = ggplot2::element_blank(),
		panel.background = ggplot2::element_blank(),
		axis.ticks = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(
      color = "black",
      size = 14
    ),
    legend.text = ggplot2::element_text(
      color = "black",
      size = 9,
      angle = 0,
      hjust = .5,
      vjust = .5
    ),
    legend.title.align = 0.5,
    legend.position = "right"
  )
}

#' GET RID OF LISTS IN A DATA FRAME
data_frame_with_list <- function(data_object) {
  set_lists_to_chars <- function(x) { 
    if(class(x) == 'list') {
      y <- paste(unlist(x[1]), sep='', collapse=', ')
    } else {
      y <- x
    } 
    return(y)
  }
  new_frame <- data.frame(
    lapply(data_object, set_lists_to_chars),
    stringsAsFactors = F
  )
  return(new_frame)
}