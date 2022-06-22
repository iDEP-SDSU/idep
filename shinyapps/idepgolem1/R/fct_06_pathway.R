#' fct_06_pathway.R This file holds all of the main data analysis functions
#' associated with eighth tab of the iDEP website.
#'
#'
#' @section fct_06_pathway.R functions:
#'
#'
#' @name fct_06_pathway.R
NULL

#' Pathway analysis with gage package
#' 
#' Run pathway analysis with the gage package using the results
#' from the limma_value function.
#' 
#' @param select_go The portion of the database to use for
#'  pathway analysis
#' @param select_contrast Comparison from DEG analysis to
#'  use the top genes from in the pathway analysis
#' @param min_set_size Minimum gene set size for a pathway
#' @param max_set_size Maximum gene set size for a pathway
#' @param limma Results list from the limma_value function
#' @param gene_p_val_cutoff Significant p-value to filter
#'  the top genes fold change by
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database (from read_gene_sets function)
#' @param absolute_fold Use the absolute value of the fold
#'  change
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of pathways to return in final
#'  result
#' 
#' @export
#' @return A data frame with the results of the pathway analysis.
#'  The data frame has five columns for the direction of the
#'  regulation, the pathway description, the stat value, the
#'  number of overlapping genes, and the p-value.
gage_data <- function(
  select_go,
  select_contrast,
  min_set_size,
  max_set_size,
  limma,
  gene_p_val_cutoff,
  gene_sets,
  absolute_fold,
  pathway_p_val_cutoff,
  n_pathway_show
) {
  if(select_go == "ID not recognized!") {
    return(as.data.frame("Gene ID not recognized."))
  }
  if(is.null(select_contrast)) {
    return(NULL)
  }
  my_range <- c(min_set_size, max_set_size)
  no_sig <- as.data.frame("No significant pathway found.")
  if(length(limma$top_genes) == 0 ) {
    return(no_sig)
  }
  
  if(length(limma$comparisons) == 1) {
    top_1 <- limma$top_genes[[1]]  
  } else {
    top <- limma$top_genes
    ix <- match(select_contrast, names(top))
    if(is.na(ix)) {
      return(no_sig)
    }
    top_1 <- top[[ix]] 
  }
  if(dim(top_1)[1] == 0) {
    return(no_sig)
  }
  colnames(top_1) <- c("Fold", "FDR")

  top_1 <- top_1[which(top_1$FDR < gene_p_val_cutoff), ]
 
  gmt <- gene_sets 
  if(length(gene_sets) == 0) {
    return(as.data.frame("No gene set found!"))
  }

  fold <- top_1[, 1]
  names(fold) <- rownames(top_1)
  if(absolute_fold) {
    fold <- abs(fold)
  }
  paths <- gage::gage(fold, gsets = gmt, ref = NULL, samp = NULL)

  paths <-  rbind(paths$greater, paths$less)
 
  if(dim(paths)[1] < 1 | dim(paths)[2] < 6 ) {
    return(no_sig)
  }
  top_1 <- paths[, c("stat.mean", "set.size", "q.val")]
  colnames(top_1) <- c("statistic", "Genes", "adj.Pval")
  top_1 <- top_1[order(top_1[, 3]), ]  
  if(length(which(top_1[, 3] <= pathway_p_val_cutoff)) == 0) {
    return(no_sig)
  }
  top_1 <- top_1[which(top_1[, 3] <= pathway_p_val_cutoff), , drop = FALSE]
  if(dim(top_1)[1] > n_pathway_show) {
    top_1 <- top_1[1:n_pathway_show, ,drop=FALSE]
  }		 
  top_1 <- as.data.frame(top_1)
  top_1 <- cbind(rep(select_contrast, dim(top_1)[1]), row.names(top_1), top_1) 
  top_1$statistic <- as.character(round(as.numeric(top_1$statistic), 4)) 
  top_1$adj.Pval <- sprintf("%-2.1e", as.numeric(top_1$adj.Pval))
  top_1[, 2] <- as.character(top_1[, 2])
  top_1[, 1] <- as.character(top_1[, 1])
  colnames(top_1)[1] <- "Direction"
  colnames(top_1)[2] <- paste("GAGE analysis:", gsub("-", " vs ", select_contrast))
  top_1[which(top_1[, 3] > 0), 1] <- "Up"
  top_1[which(top_1[, 3] < 0), 1] <- "Down"
  top_1 <- top_1[order(top_1[, 1], -abs(as.numeric(top_1[, 3]))), ]
  top_1[duplicated(top_1[, 1]), 1] <- ""
  rownames(top_1) <- seq(1, nrow(top_1), 1)

  return(top_1)
}

#' Pathway analysis with the PGSEA package
#' 
#' Run pathway analysis with the PGSEA package using the results
#' from the limma_value function.
#' 
#' @param processed_data Data that has been through the pre-processing
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database (from read_gene_sets function)
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of significant pathways to show
#' 
#' @export
#' @return A list with a data frame and a numeric value that is used
#'  in the plot_pgsea function to create a heatmap.
pgsea_data <- function(
  processed_data,
  gene_sets,
  my_range,
  pathway_p_val_cutoff,
  n_pathway_show
){
  subtype = detect_groups(colnames(processed_data))
	
  # Cut off to report in PGSEA. Otherwise NA
	p_value <- 0.01  
	if(length(gene_sets) == 0) {
    return(list(pg3 = NULL, best = 1))
  }
	
	pg <- PGSEA::PGSEA(
    processed_data - rowMeans(processed_data),
    cl = gene_sets,
    range = my_range,
    p.value = TRUE,
    weighted = FALSE
  )
	
	pg_results <- pg$results
  # Remove se/wrts with all missing(non-signficant)
	pg_results <- pg_results[rowSums(is.na(pg_results)) < ncol(pg_results), ]
	if(dim(pg_results)[1] < 2) {
    return()
  }
	best <- max(abs(pg_results))
	
	if(length(subtype) < 4 || length(unique(subtype)) <2 ||
     length(unique(subtype)) == dim(processed_data)[2]) { 
	  pg_results <- pg_results[order(-apply(pg_results, 1, sd)), ]
	  return(list(pg_data = pg_results[1:top, ], best <- best ))
	} 
    
	cat("\nComputing P values using ANOVA\n");
	path_p_value <- function (
    k,
    pg_results,
    subtype
  ){
	  return(summary(aov(pg_results[k, ] ~ subtype))[[1]][["Pr(>F)"]][1])
	}
	p_values <- sapply(1:dim(pg_results)[1], function(x) {
    path_p_value(
      k = x,
      pg_results = pg_results,
      subtype = subtype)
  })
	p_values <- stats::p.adjust(p_values, "fdr")
	

  if(sort(p_values)[2] > pathway_p_val_cutoff) {
    return(list(pg_data = NULL, best = best)) 
  } else {  
    n_sig_t <- rowSums(pg$p.results < p_value)
	
	  result <- cbind(as.matrix(p_values), n_sig_t, pg_results) 
	  result <- result[order(result[, 1]), ]
    result <- result[which(result[, 1] < pathway_p_val_cutoff), , drop = F]
	
	  pg_results <- result[, -2]

	  # When there is only 1 left in the matrix pg_results becomes a vector
	  if(sum(p_values < pathway_p_val_cutoff) == 1) {
      pg_data <- t(as.matrix(pg_results))
      pg_data <- rbind(pg_data, pg_data)
    } else {
      if(dim(pg_results)[1] > n_pathway_show) {
        pg_data <- pg_results[1:n_pathway_show, ]
      } else {
        pg_data <- pg_results
      }
    }

	  rownames(pg_data) <- sapply(rownames(pg_data) , extract_under)
	  a <- sprintf("%-3.2e", pg_data[, 1])
	  rownames(pg_data) <- paste(a, rownames(pg_data), sep = " ")
	  pg_data <- pg_data[, -1]
	
    # Sort by SD
	  pg_data <- pg_data[order(-apply(pg_data, 1, sd)), ]
  
    return(list(
      pg_data = pg_data,
      best = best
    ))
  }
 }

#' Heatmap of PGSEA pathway analysis
#' 
#' Create a heatmap from the pathway analysis using the PGSEA
#' package. The heatmap shows the expression in each group for
#' each significantly enriched pathway.
#' 
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param processed_data Data that has been through the pre-processing
#' @param contrast_samples Sample columns that correspond to the
#'  selected comparison
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database (from read_gene_sets function)
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of significant pathways to show
#' 
#' @export
#' @return A heatmap plot with the rows as the significant
#'  pathways and the columns corresponding to the samples.
plot_pgsea <- function(
  my_range,
  processed_data,
  contrast_samples,
  gene_sets,
  pathway_p_val_cutoff,
  n_pathway_show
) {
  genes <- processed_data[, contrast_samples]	
	if(length(gene_sets)  == 0)  {
    return(
      NULL
    )
  } else {
	  subtype <- detect_groups(colnames(genes))
	  result <- pgsea_data(
      processed_data = genes,
      gene_sets = gene_sets,
      my_range = my_range,
      pathway_p_val_cutoff = pathway_p_val_cutoff,
      n_pathway_show = n_pathway_show
    )
					 
	  if(is.null(result$pg_data)) {
      return(
        NULL
      )
    } else {
      PGSEA::smcPlot(
        result$pg_data,
        factor(subtype),
        scale = c(-max(result$pg_data), max(result$pg_data)),
        show.grid = T,
        margins = c(3, 1, 13, 38),
        col = PGSEA::.rwb,
        cex.lab = 0.5
      )
    }
  }
}

#' Pathway analysis with the FGSEA package
#' 
#' Run pathway analysis with the FGSEA package using the results
#' from the limma_value function.
#' 
#' @param select_contrast Comparison from DEG analysis to
#'  use the top genes from in the pathway analysis
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param limma Results list from the limma_value function
#' @param gene_p_val_cutoff Significant p-value to filter
#'  the top genes fold change by
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database (from read_gene_sets function)
#' @param absolute_fold Use the absolute value of the fold
#'  change
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of pathways to return in final
#'  result
#' 
#' @export
#' @return A data frame with the results of the pathway analysis.
#'  The data frame has five columns for the direction of the
#'  regulation, the pathway description, the stat value, the
#'  number of overlapping genes, and the p-value.
fgsea_data <- function(
  select_contrast,
  my_range,
  limma,
  gene_p_val_cutoff,
  gene_sets,
  absolute_fold,
  pathway_p_val_cutoff,
  n_pathway_show
) {
	no_sig <- as.data.frame("No significant pathway found.")
	if(length(limma$top_genes) == 0) {
    return(no_sig)
  }
	if(length(limma$comparisons) == 1) {
    top_1 <- limma$top_genes[[1]]  
	} else {
	  top <- limma$top_genes
	  ix <- match(select_contrast, names(top))
	  if(is.na(ix)) {
      return(no_sig)
    }
	  top_1 <- top[[ix]] 
	}
	if(dim(top_1)[1] == 0) {
    return(no_sig)
  }
	colnames(top_1) <- c("Fold","FDR")
	  
	# Remove some genes
	top_1 <- top_1[which(top_1$FDR < gene_p_val_cutoff), ]
	  
  if(length(gene_sets) == 0) {
    return(as.data.frame("No gene set found!"))
  }


	fold <- top_1[, 1]
  names(fold) <- rownames(top_1)
	
  # Use absolute value of fold change, disregard direction
  if(absolute_fold) {
    fold <- abs(fold)
  }
	 
  paths <- fgsea::fgsea(
    pathways = gene_sets, 
    stats = fold,
    minSize = my_range[1],
    maxSize = my_range[2],
    nPerm = 100000                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ,
		nproc = 6
  )
	
  if(dim(paths)[1] < 1) {
    return(no_sig)
  }
	paths <- as.data.frame(paths)
  # Sort by NES
  paths <- paths[order(-abs(paths[, 5])), ]
	top_1 <- paths[, c(1, 5, 7, 3)]
	colnames(top_1) <- c("Pathway", "NES", "Genes", "adj.Pval")
	  
	if(length(which(top_1[, 4] <= pathway_p_val_cutoff)) == 0) {
    return(no_sig)
  }
	top_1 <- top_1[which(top_1[, 4] <= pathway_p_val_cutoff), , drop = FALSE]
	
  if(dim(top_1)[1] > n_pathway_show) {
    top_1 <- top_1[1:n_pathway_show, , drop = FALSE]
  }
	
	top_1 <- as.data.frame(top_1)
	top_1 <- cbind(rep(select_contrast, dim(top_1)[1]), top_1) 
	top_1[, 4] <- as.character(round(as.numeric(top_1[, 4]), 4)) 
	top_1$adj.Pval <- sprintf("%-2.1e", as.numeric(top_1$adj.Pval))
	top_1[, 1] <- as.character(top_1[, 1])
	colnames(top_1)[1] <- "Direction"
	colnames(top_1)[2] <- paste("GSEA analysis:", gsub("-"," vs ", select_contrast))
	top_1[which(as.numeric(top_1[, 3]) > 0), 1] <- "Up"
  top_1[which(as.numeric(top_1[, 3]) < 0), 1] <- "Down"
	top_1 <- top_1[order(top_1[, 1], -abs(as.numeric(top_1[, 3]))), ]
	top_1[duplicated(top_1[, 1]), 1] <- ""	 
	top_1[, 3] <- as.character(round(as.numeric(top_1[, 3]), 4))

	return(top_1)
}

#' Pathway analysis with reactome package
#' 
#' Run pathway analysis with the reactome package using the results
#' from the limma_value function.
#' 
#' @param select_contrast Comparison from DEG analysis to
#'  use the top genes from in the pathway analysis
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param limma Results list from the \code{limma_value}
#'  function
#' @param gene_p_val_cutoff Significant p-value to filter
#'  the top genes fold change by
#' @param converted Return value from the convert_id function, contains
#'  information about the gene IDs for the matched species
#' @param idep_data Read data files from the database
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of pathways to return in final
#'  result
#' @param absolute_fold Use the absolute value of the fold
#'  change
#' 
#' @export
#' @return A data frame with the results of the pathway analysis.
#'  The data frame has five columns for the direction of the
#'  regulation, the pathway description, the stat value, the
#'  number of overlapping genes, and the p-value.
reactome_data <- function(
  select_contrast,
  my_range,
  limma,
  gene_p_val_cutoff,
  converted,
  idep_data,
  pathway_p_val_cutoff,
  n_pathway_show,
  absolute_fold
) {
  ensembl_species <- c(
    "hsapiens_gene_ensembl","rnorvegicus_gene_ensembl", "mmusculus_gene_ensembl",
	  "celegans_gene_ensembl","scerevisiae_gene_ensembl", "drerio_gene_ensembl",
    "dmelanogaster_gene_ensembl"
  )
  reactome_pa_species <- c("human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly" )
	no_sig <- as.data.frame("No significant pathway found.")
	if(length(limma$top_genes) == 0) {
    return(no_sig)
  }
	if(length(limma$comparisons) == 1) {
    top_1 <- limma$top_genes[[1]]  
	} else {
	  top <- limma$top_genes
	  ix <- match(select_contrast, names(top))
	  if(is.na(ix)) {
      return(no_sig)
    }
	  top_1 <- top[[ix]] 
	}
	if(dim(top_1)[1] == 0) {
    return(no_sig)
  }
	colnames(top_1) <- c("Fold", "FDR")

	# Remove some genes
	top_1 <- top_1[which(top_1$FDR < gene_p_val_cutoff), ]
	 
	fold <- top_1[, 1]
  names(fold) <- rownames(top_1)
	if(absolute_fold) {
    # Use absolute value of fold change, disregard direction
    fold <- abs(fold) 
  }
  
  species <- converted$species[1, 1]
  ix <- match(species, ensembl_species)	
  
  if(is.na(ix)) {
    return(as.data.frame("Species not coverted by ReactomePA package!"))
  }
	  
	fold <- convert_ensembl_to_entrez(
    query = fold,
    species = species,
    org_info = idep_data$org_info
  )  
	
	
  fold <- sort(fold, decreasing = T)
	paths <- ReactomePA::gsePathway(
    fold,
    nPerm = 5000,
    organism = reactome_pa_species[ix],
    minGSSize = my_range[1], 
		maxGSSize = my_range[2],
		pvalueCutoff = 0.5,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  
  paths <- as.data.frame(paths)
	  
	if(is.null(paths)) {
    return(no_sig)
  }
	if(dim(paths)[1] == 0) {
    return(no_sig)
  }

  if(dim(paths)[1] < 1) {
    return(no_sig)
  }
  paths <- as.data.frame(paths)
  paths <- paths[order(-abs(paths[, 5])), ]

	top_1 <- paths[, c(2, 5, 3, 7)]

	colnames(top_1) <- c("Pathway", "NES", "Genes", "adj.Pval")
	  
	if(length(which(top_1[, 4] <= pathway_p_val_cutoff)) == 0) {
    return(no_sig)
  }  
  top_1 <- top_1[which(top_1[, 4] <= pathway_p_val_cutoff), , drop = FALSE]
	
  if(dim(top_1)[1] > n_pathway_show) {
    top_1 <- top_1[1:n_pathway_show, , drop = FALSE]
  }
 	
	top_1 <- as.data.frame(top_1)
	top_1 <- cbind(rep(select_contrast, dim(top_1)[1]), top_1) 
	top_1[, 4] <- as.character(round(as.numeric(top_1[, 4]), 4)) 
	top_1$adj.Pval <- sprintf("%-2.1e", as.numeric(top_1$adj.Pval))
	top_1[, 1] <- as.character(top_1[, 1])
	colnames(top_1)[1] <- "Direction"
	colnames(top_1)[2] <- paste(
    "ReactomePA analysis:",
    gsub("-"," vs ", select_contrast)
  )
	top_1[which(as.numeric(top_1[, 3]) > 0), 1] <- "Up"
	top_1[which(as.numeric(top_1[, 3]) < 0), 1] <- "Down"
	top_1 <- top_1[order(top_1[, 1], -abs(as.numeric(top_1[, 3]))), ]
	top_1[duplicated(top_1[, 1]), 1] <- ""	 
  top_1[, 3] <- as.character(round(as.numeric(top_1[, 3]), 4))
	
  return(top_1)
}

#' Pathway analysis with the PGSEA package on all samples
#' 
#' Run pathway analysis with the PGSEA package using the results
#' from the limma_value function on all samples.
#' 
#' @param go Portion of the database to use for the pathway analysis
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param data Data that has been through the pre-processing
#' @param select_contrast Comparison from DEG analysis to
#'  use the top genes from in the pathway analysis
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database (from read_gene_sets function)
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of pathways to return in final
#'  result
#' 
#' @export
#' @return A data frame with the results of the pathway analysis.
#'  The data frame has five columns for the direction of the
#'  regulation, the pathway description, the stat value, the
#'  number of overlapping genes, and the p-value.
pgsea_plot_all <- function(
  go,
  my_range,
  data,
  select_contrast,
  gene_sets,
  pathway_p_val_cutoff,
  n_pathway_show
) {
  if(length(gene_sets)  == 0) {
    plot.new()
    text(0, 1, "No gene sets!")
  } else {
    subtype <- detect_groups(colnames(data)) 
	  result <- pgsea_data(
      processed_data = data,
      gene_sets = gene_sets,
      my_range = my_range,
      pathway_p_val_cutoff = pathway_p_val_cutoff,
      n_pathway_show = n_pathway_show
    )
    if(is.null(result$pg_data)) {
      plot.new()
      text(0.5, 1, "No significant pathway found!")
    } else {
      PGSEA::smcPlot(
        result$pg_data,
        factor(subtype),
        scale = c(-max(result$pg_data), max(result$pg_data)),
        show.grid = T,
        margins = c(3, 1, 13, 38),
        col = PGSEA::.rwb,
        cex.lab = 0.5
      )
    }
  }
}

#' Data from PGSEA plot
#' 
#' Get the data matrix that is plotted in the heatmap created by
#' the plot_pgsea function.
#' 
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param data Data that has been through the pre-processing
#' @param select_contrast Comparison from DEG analysis to
#'  use the top genes from in the pathway analysis
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database (from read_gene_sets function)
#' @param sample_info Experiment file information for grouping
#' @param select_factors_model The selected factors for the model
#'  expression
#' @param select_model_comprions Selected comparisons to analyze
#'  in the DEG analysis
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of pathways to return in final
#'  result
#' 
#' @export
#' @return Data matrix with the rownames the descriptions of pathways
#'  and the matrix the returned expression calculation from the PGSEA
#'  package.
get_pgsea_plot_data <- function(
  my_range,
  data,
  select_contrast,
  gene_sets,
  sample_info,
  select_factors_model,
  select_model_comprions,
  pathway_p_val_cutoff,
  n_pathway_show
) {	
	# Find sample related to the comparison
	iz <- match(detect_groups(colnames(data)), unlist(strsplit(select_contrast, "-")))
	iz <- which(!is.na(iz))
	
  if(!is.null(sample_info) & !is.null(select_factors_model) & length(select_model_comprions) > 0 ) {
		# Strings like: "groups: mutant vs. control"
    comparisons <- gsub(".*: ", "", select_model_comprions)
		comparisons <- gsub(" vs\\. ", "-", comparisons)
    # Corresponding factors
		factors_vector <- gsub(":.*", "", select_model_comprions)
    # Selected contrast lookes like: "mutant-control"
		ik <- match(select_contrast, comparisons)
		if(is.na(ik)) {
      iz <- 1:(dim(data)[2]) 
    } else {
      # Interaction term, use all samples	
      # Corresponding factors
			selected_factor <- factors_vector[ik]
			iz <- match(sample_info[, selected_factor], unlist(strsplit(select_contrast, "-")))
			iz <- which(!is.na(iz))
		}
	}

	if(grepl("I:", select_contrast)) {
    # If it is factor design use all samples
    iz <- 1:(dim(data)[2])
  } 
	if(is.na(iz)[1] | length(iz) <= 1) {
    iz <- 1:(dim(data)[2]) 
  }
	
	genes <- data
	genes <- genes[, iz]

	subtype <- detect_groups(colnames(genes)) 
  
  if(length(gene_sets)  == 0) {
    return(as.data.frame("No significant pathway!"))
  } else {
    result <- pgsea_data(
      processed_data = genes,
      gene_sets = gene_sets,
      my_range = my_range,
      pathway_p_val_cutoff = pathway_p_val_cutoff,
      n_pathway_show = n_pathway_show
    )
					 
	  if(is.null(result$pg_data)) {
      return(as.data.frame("No significant pathway!"))
    } else {
      return( as.data.frame(result$pg_data) )
    } 
  }
}

#' Data from PGSEA plot all samples
#' 
#' Get the data matrix that is plotted in the heatmap created by
#' the pgsea_plot_all function.
#' 
#' @param data Data that has been through the pre-processing
#' @param select_contrast Comparison from DEG analysis to
#'  use the top genes from in the pathway analysis
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database (from read_gene_sets function)
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of pathways to return in final
#'  result
#' 
#' @export
#' @return Data matrix with the rownames the descriptions of pathways
#'  and the matrix the returned expression calculation from the PGSEA
#'  package.
get_pgsea_plot_all_samples_data <- function(
  data,
  select_contrast,
  gene_sets,
  my_range,
  pathway_p_val_cutoff,
  n_pathway_show
) {
  genes <- data
	subtype <- detect_groups(colnames(genes)) 
  
  if(length(gene_sets)  == 0) {
    plot.new()
    text(0,1, "No gene sets!")
  } else {
	  result <- pgsea_data(
      processed_data = genes,
      gene_sets = gene_sets,
      my_range = my_range,
      pathway_p_val_cutoff = pathway_p_val_cutoff,
      n_pathway_show = n_pathway_show
    )
					 
	  if(is.null(result$pg_data)) {
      return(as.data.frame("No significant pathway!"))
    } else {
      return(as.data.frame(result$pg_data))
    }
  }	
}

#' Get data from genes in selected pathway
#' 
#' Return a data matrix that is a subset of the processed data and
#' only contains genes that are in the gene set of the desired
#' pathway.
#' 
#' @param sig_pathways Description of the pathway for which to
#'  obtain the gene expression data
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database (from read_gene_sets function)
#' @param contrast_samples Sample columns that correspond to the
#'  selected comparison
#' @param data Data that has been through the pre-processing
#' @param select_org The organism that the gene data is for
#' @param all_gene_names Matrix of all the matched and converted
#'  gene IDs
#' 
#' @export
#' @return Sub-data matrix from the processed data. Only contains
#'  genes from the selected pathway and samples that correspond to
#'  the comparison being analyzed.
pathway_select_data <- function(
  sig_pathways,
  gene_sets,
  contrast_samples,
  data,
  select_org,
  all_gene_names
) {
  if(sig_pathways == "All") {
    return (NULL) 
  }
  # Find the gene set
	ix <- which(names(gene_sets) == sig_pathways)
	if(length(ix) == 0) {
    return(NULL)
  }
  # Retrieve genes
  genes <- gene_sets[[ix]]

	# Find related samples	
	iz <-contrast_samples
	x <- data[which(rownames(data) %in% genes), iz]
	if(ncol(all_gene_names) == 3) {
    x <- rowname_id_swap(
      data_matrix = x,
      all_gene_names = all_gene_names,
      select_gene_id = "symbol"
    )
  }
	
	return(x)
}

#' Create pathway table with gene sets
#' 
#' Create a data frame of significant pathways and their analysis
#' values. Also add a column that contains the gene sets for the
#' pathway. 
#' 
#' @param pathway_method What method is being used for the pathway
#'  analysis
#' @param gage_pathway_data Return matrix from \code{gage_data}
#'  function
#' @param fgsea_pathway_data Return matrix from \code{fgsea_data}
#'  function
#' @param pgsea_plot_data Return matrix from \code{get_pgsea_plot_data}
#'  function
#' @param pgsea_plot_all_samples_data Return matrix from
#'  \code{get_pgsea_plot_all_samples_data} function
#' @param go Portion of the database to use for the pathway analysis
#' @param select_org Organism that the gene data is for
#' @param gene_info Return from \code{gene_info} function, all gene
#'  info from the database query with the User gene IDs
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database (from read_gene_sets function)
#' 
#' @export
#' @return A data frame with the pathway analysis statistics and 
#'  the gene sets for each significantly enriched pathway.
get_pathway_list_data <- function(
  pathway_method,
  gage_pathway_data,
  fgsea_pathway_data,
  pgsea_plot_data,
  pgsea_plot_all_samples_data,
  go,
  select_org,
  gene_info,
  gene_sets
) {
	pathways <- NULL
	if(pathway_method == 1) {
    if(!is.null(gage_pathway_data)) {
      if(dim(gage_pathway_data)[2] > 1) { 
				pathways <- gage_pathway_data
				colnames(pathways)[2] <- "Pathways" 	
				colnames(pathways)[4] <- "nGenes"
			}
    }
  }
	if(pathway_method == 3) {
    if(!is.null(fgsea_pathway_data)) {
      if(dim(fgsea_pathway_data)[2] > 1) {
				pathways <- fgsea_pathway_data
				colnames(pathways)[2] <- "Pathways" 	
				colnames(pathways)[4] <- "nGenes" 
			}
    }
  }
	if(pathway_method == 2) {
    if(!is.null(pgsea_plot_data)) {
      if(dim(pgsea_plot_data)[2] > 1) {
				pathways <- as.data.frame(pgsea_plot_data)
				pathways$Pathways <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
				pathways$adj.Pval <- gsub(" .*", "", rownames(pathways))
				pathways$Direction <- "Diff"
			}
    }
  }
	if(pathway_method == 4) {
    if(!is.null(pgsea_plot_all_samples_data)) {
      if(dim(pgsea_plot_all_samples_data)[2] >1 ) {
				pathways <- as.data.frame(pgsea_plot_all_samples_data)
				pathways$Pathways <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
				pathways$adj.Pval <- gsub(" .*", "", rownames(pathways))
				pathways$Direction <- "Diff"	
			}
    }
  }	
	if(is.null(pathways)) {
    return(NULL)
  }	
	# if no gene set data, return pathway list
	if(is.null(gene_sets)) {
    return(pathways)
	}
	
	pathways$adj_p_val <- as.numeric(pathways$adj.Pval)
  pathways <- subset(pathways, select = -c(adj.Pval))
  pathways$adj_p_val <- as.character(pathways$adj_p_val)

  # Sometimes only one pathway is in the table
	if(nrow(pathways) > 1) {
    for(i in 2:nrow(pathways)) {
      if(nchar(pathways$Direction[i]) <= 1) {
			  pathways$Direction[i] = pathways$Direction[i-1]
      }
    }
	}	

	# Gene symbol matching symbols 
	probe_to_gene <- NULL
	if(go != "ID not recognized!" & select_org != "NEW") {
    # If more than 50% genes has symbol
    if(sum(is.na(gene_info$symbol)) / dim(gene_info)[1] < .5) { 
		  probe_to_gene <- gene_info[, c("ensembl_gene_id", "symbol")]
		  probe_to_gene$symbol <- gsub(" ", "", probe_to_gene$symbol)

		  ix <- which(
        is.na(probe_to_gene$symbol) |
				nchar(probe_to_gene$symbol) < 2 | 
				toupper(probe_to_gene$symbol) == "NA" |  
				toupper(probe_to_gene$symbol) == "0"
      )
      # Use gene ID
		  probe_to_gene[ix, 2] <- probe_to_gene[ix, 1]  
	  }
  }

	pathways$Genes <- vector(mode = "list", length = nrow(pathways))
	# looking up genes for each pathway
	for(i in 1:nrow(pathways)) {
    # Find the gene set
		ix <- which(names(gene_sets) == pathways$Pathways[i])
		if(length(ix) != 0) {
      # Retrieve genes
			genes <- gene_sets[[ix]]
			
			if(!is.null(probe_to_gene)) { 
				iy <- match(genes, probe_to_gene[, 1])
				genes <- probe_to_gene[iy, 2]
			}
			pathways$Genes[[i]] <- c(genes)
		}
	}
  
	return(pathways)
}

#' Use KEGG to create a pathway diagram
#' 
#' In the databse, use the KEGG information to create an image
#' that is a diagram of the pathway that is being enriched.
#' 
#' @param go Portion of the database to use for the pathway analysis
#' @param gage_pathway_data Return matrix from \code{gage_data}
#'  function
#' @param sig_pathways Description of the pathway for which to
#'  obtain the gene expression data
#' @param select_contrast Comparison from DEG analysis to
#'  use the top genes from in the pathway analysis
#' @param limma Results list from the \code{limma_value}
#'  function
#' @param converted Return value from the convert_id function, contains
#'  information about the gene IDs for the matched species
#' @param idep_data Read data files from the database
#' @param select_org The organism that the gene data is for
#' @param low_color Color value for the low-ly expressed genes 
#' @param high_color Color vlaue for the high-ly expressed genes 
#' 
#' @export
#' @return Make an image and return the path to the image to be
#'  rendered in the server.
kegg_pathway <- function(
  go,
  gage_pathway_data,
  sig_pathways,
  select_contrast,
  limma,
  converted,
  idep_data,
  select_org, 
  low_color = "green", 
  high_color = "red"
) {
  # First generate a blank image. Otherwise return(NULL) gives us errors.
  out_file <- tempfile(fileext = '.png')
  png(out_file, width = 400, height = 300)

  frame()
	dev.off()
  blank <- list(
    src = out_file,
    contentType = 'image/png',
    width = 400,
    height = 300,
    alt = " "
  )	
  

	if(is.null(go) || go != "KEGG") {
    return(blank)
  }
	if(is.null(gage_pathway_data)) {
    return(blank)
  }
	if(is.null(sig_pathways)) {
    return(blank)
  }

  # These two functions are from the pathview package, modified to write to a designated folder: temp.
  mypathview <- function(
    gene.data = NULL,
    cpd.data = NULL,
    pathway.id,
    species = "hsa", 
    kegg.dir = ".",
    cpd.idtype = "kegg",
    gene.idtype = "entrez", 
    gene.annotpkg = NULL,
    min.nnodes = 3,
    kegg.native = TRUE, 
    map.null = TRUE,
    expand.node = FALSE,
    split.group = FALSE, 
    map.symbol = TRUE,
    map.cpdname = TRUE,
    node.sum = "sum", 
    discrete = list(gene = FALSE, cpd = FALSE),
    limit = list(gene = 1, cpd = 1),
    bins = list(gene = 10, cpd = 10),
    both.dirs = list(gene = T, cpd = T),
    trans.fun = list(gene = NULL, cpd = NULL), 
    low = list(gene = low_color, cpd = "blue"),
    mid = list(gene = "gray", cpd = "gray"),
    high = list(gene = high_color, cpd = "yellow"), 
    na.col = "transparent",
    ...
  ) {
    dtypes = !is.null(gene.data) + (!is.null(cpd.data))
    cond0 = dtypes == 1 & is.numeric(limit) & length(limit) > 1
    if (cond0) {
      if (limit[1] != limit[2] & is.null(names(limit))) 
        limit = list(gene = limit[1:2], cpd = limit[1:2])
    }
    if (is.null(trans.fun)) 
      trans.fun = list(gene = NULL, cpd = NULL)
    arg.len2 = c(
      "discrete", "limit", "bins", "both.dirs", "trans.fun", 
      "low", "mid", "high"
    )
    for (arg in arg.len2) {
      obj1 = eval(as.name(arg))
      if (length(obj1) == 1) 
        obj1 = rep(obj1, 2)
      if (length(obj1) > 2) 
        obj1 = obj1[1:2]
      obj1 = as.list(obj1)
      ns = names(obj1)
      if (length(ns) == 0 | !all(c("gene", "cpd") %in% ns)) 
        names(obj1) = c("gene", "cpd")
      assign(arg, obj1)
    }
    if(is.character(gene.data)) {
      gd.names = gene.data
      gene.data = rep(1, length(gene.data))
      names(gene.data) = gd.names
      both.dirs$gene = FALSE
      ng = length(gene.data)
      nsamp.g = 1
    } else if(!is.null(gene.data)) {
      if(length(dim(gene.data)) == 2) {
        gd.names = rownames(gene.data)
        ng = nrow(gene.data)
        nsamp.g = 2
      } else if(is.numeric(gene.data) & is.null(dim(gene.data))) {
        gd.names = names(gene.data)
        ng = length(gene.data)
        nsamp.g = 1
      } else stop("wrong gene.data format!")
    } else if(is.null(cpd.data)) {
      stop("gene.data and cpd.data are both NULL!")
    }
    gene.idtype = toupper(gene.idtype)
    bods <- pathview::bods
    if(species != "ko") {
      species.data = pathview::kegg.species.code(
        species,
        na.rm = T, 
        code.only = FALSE
      )
    } else {
      species.data = c(
        kegg.code = "ko",
        entrez.gnodes = "0", 
        kegg.geneid = "K01488",
        ncbi.geneid = NA,
        ncbi.proteinid = NA, 
        uniprot = NA
      )
      gene.idtype = "KEGG"
      msg.fmt = "Only KEGG ortholog gene ID is supported, make sure it looks like \"%s\"!"
      msg = sprintf(msg.fmt, species.data["kegg.geneid"])
      message("Note: ", msg)
    }
    if (length(dim(species.data)) == 2) {
      message("Note: ", "More than two valide species!")
      species.data = species.data[1, ]
    }
    species = species.data["kegg.code"]
    entrez.gnodes = species.data["entrez.gnodes"] == 1
    if(is.na(species.data["ncbi.geneid"])) {
      if(!is.na(species.data["kegg.geneid"])) {
        msg.fmt = "Mapping via KEGG gene ID (not Entrez) is supported for this species,\nit looks like \"%s\"!"
        msg = sprintf(msg.fmt, species.data["kegg.geneid"])
        message("Note: ", msg)
      } else {
        stop("This species is not annotated in KEGG!")
      }
    }
    if (is.null(gene.annotpkg)) 
      gene.annotpkg = bods[match(species, bods[, 3]), 1]
    if(
      length(grep("ENTREZ|KEGG|NCBIPROT|UNIPROT", gene.idtype)) < 1 & !is.null(gene.data)
    ) {
      if(is.na(gene.annotpkg)) 
        stop("No proper gene annotation package available!")
      if(!gene.idtype %in% gene.idtype.bods[[species]]) 
        stop("Wrong input gene ID type!")
      gene.idmap = pathview::id2eg(
        gd.names,
        category = gene.idtype, 
        pkg.name = gene.annotpkg,
        unique.map = F
      )
      gene.data = pathview::mol.sum(gene.data, gene.idmap)
      gene.idtype = "ENTREZ"
    }
    if(gene.idtype != "KEGG" & !entrez.gnodes & !is.null(gene.data)) {
      id.type = gene.idtype
      if(id.type == "ENTREZ") 
        id.type = "ENTREZID"
      kid.map = names(species.data)[-c(1:2)]
      kid.types = names(kid.map) = c(
        "KEGG", "ENTREZID", "NCBIPROT", "UNIPROT"
      )
      kid.map2 = gsub("[.]", "-", kid.map)
      kid.map2["UNIPROT"] = "up"
      if (is.na(kid.map[id.type])) 
        stop("Wrong input gene ID type for the species!")
      message("Info: Getting gene ID data from KEGG...")
      gene.idmap = KEGGREST::keggConv(kid.map2[id.type], species)
      message("Info: Done with data retrieval!")
      kegg.ids = gsub(paste(species, ":", sep = ""), "", names(gene.idmap))
      in.ids = gsub(paste0(kid.map2[id.type], ":"), "", gene.idmap)
      gene.idmap = cbind(in.ids, kegg.ids)
      gene.data = pathview::mol.sum(gene.data, gene.idmap)
      gene.idtype = "KEGG"
    }
    if(is.character(cpd.data)) {
      cpdd.names = cpd.data
      cpd.data = rep(1, length(cpd.data))
      names(cpd.data) = cpdd.names
      both.dirs$cpd = FALSE
      ncpd = length(cpd.data)
    } else if(!is.null(cpd.data)) {
      if (length(dim(cpd.data)) == 2) {
        cpdd.names = rownames(cpd.data)
        ncpd = nrow(cpd.data)
      } else if(is.numeric(cpd.data) & is.null(dim(cpd.data))) {
        cpdd.names = names(cpd.data)
        ncpd = length(cpd.data)
      }
      else stop("wrong cpd.data format!")
    }
    if(length(grep("kegg", cpd.idtype)) < 1 & !is.null(cpd.data)) {
      data(rn.list)
      cpd.types = c(names(rn.list), "name")
      cpd.types = tolower(cpd.types)
      cpd.types = cpd.types[-grep("kegg", cpd.types)]
      if (!tolower(cpd.idtype) %in% cpd.types) 
        stop("Wrong input cpd ID type!")
      cpd.idmap = pathview::cpd2kegg(cpdd.names, in.type = cpd.idtype)
      cpd.data = pathview::mol.sum(cpd.data, cpd.idmap)
    }
    warn.fmt = "Parsing %s file failed, please check the file!"
    if(length(grep(species, pathway.id)) > 0) {
      pathway.name = pathway.id
      pathway.id = gsub(species, "", pathway.id)
    } else pathway.name = paste(species, pathway.id, sep = "")
    kfiles = list.files(path = kegg.dir, pattern = "[.]xml|[.]png")
    npath = length(pathway.id)
    out.list = list()
    tfiles.xml = paste(pathway.name, "xml", sep = ".")
    tfiles.png = paste(pathway.name, "png", sep = ".")
    if(kegg.native) 
      ttype = c("xml", "png")
    else ttype = "xml"
    xml.file <- paste(kegg.dir, "/", tfiles.xml, sep = "")
    for(i in 1:npath) {
      if (kegg.native) 
        tfiles = c(tfiles.xml[i], tfiles.png[i])
      else tfiles = tfiles.xml[i]
      if (!all(tfiles %in% kfiles)) {
        dstatus = pathview::download.kegg(
          pathway.id = pathway.id[i], 
          species = species,
          kegg.dir = kegg.dir,
          file.type = ttype
        )
        if (dstatus == "failed") {
          warn.fmt = "Failed to download KEGG xml/png files, %s skipped!"
          warn.msg = sprintf(warn.fmt, pathway.name[i])
          message("Warning: ", warn.msg)
          return(invisible(0))
        }
      }
      if(kegg.native) {
        node.data = try(pathview::node.info(xml.file[i]), silent = T)
        if (class(node.data) == "try-error") {
          warn.msg = sprintf(warn.fmt, xml.file[i])
          message("Warning: ", warn.msg)
          return(invisible(0))
        }
        node.type = c("gene", "enzyme", "compound", "ortholog")
        sel.idx = node.data$type %in% node.type
        nna.idx = !is.na(
          node.data$x + node.data$y + node.data$width + node.data$height
        )
        sel.idx = sel.idx & nna.idx
        if(sum(sel.idx) < min.nnodes) {
          warn.fmt = "Number of mappable nodes is below %d, %s skipped!"
          warn.msg = sprintf(warn.fmt, min.nnodes, pathway.name[i])
          message("Warning: ", warn.msg)
          return(invisible(0))
        }
        node.data = lapply(node.data, "[", sel.idx)
      } else {
        gR1 = try(
          pathview::parseKGML2Graph2(
            xml.file[i],
            genes = F, 
            expand = expand.node,
            split.group = split.group
          ), 
          silent = T
        )
        node.data = try(
          pathview::node.info(gR1),
          silent = T
        )
        if (class(node.data) == "try-error") {
          warn.msg = sprintf(warn.fmt, xml.file[i])
          message("Warning: ", warn.msg)
          return(invisible(0))
        }
      }
      if(species == "ko") 
        gene.node.type = "ortholog"
      else gene.node.type = "gene"
      if ((
        !is.null(gene.data) | map.null) & sum(node.data$type == gene.node.type) > 1
      ) {
        plot.data.gene = pathview::node.map(
          gene.data,
          node.data,
          node.types = gene.node.type, 
          node.sum = node.sum,
          entrez.gnodes = entrez.gnodes
        )
        kng = plot.data.gene$kegg.names
        kng.char = gsub("[0-9]", "", unlist(kng))
        if (any(kng.char > "")) 
          entrez.gnodes = FALSE
        if (map.symbol & species != "ko" & entrez.gnodes) {
          if (is.na(gene.annotpkg)) {
            warn.fmt = "No annotation package for the species %s, gene symbols not mapped!"
            warn.msg = sprintf(warn.fmt, species)
            message("Warning: ", warn.msg)
          } else {
            # Try to fix this error: Error in $<-.data.frame: replacement has 97 rows, data has 103
			      plot.data.gene$labels = NA 
            plot.data.gene$labels = pathview::eg2id(
              as.character(plot.data.gene$kegg.names), 
              category = "SYMBOL",
              pkg.name = gene.annotpkg
            )[, 2]
            mapped.gnodes = rownames(plot.data.gene)
            node.data$labels[mapped.gnodes] = plot.data.gene$labels
          }
        }
        cols.ts.gene = pathview::node.color(
          plot.data.gene,
          limit$gene, 
          bins$gene,
          both.dirs = both.dirs$gene,
          trans.fun = trans.fun$gene, 
          discrete = discrete$gene,
          low = low$gene,
          mid = mid$gene, 
          high = high$gene,
          na.col = na.col
        )
      } else plot.data.gene = cols.ts.gene = NULL
      if ((
        !is.null(cpd.data) | map.null) & sum(node.data$type == "compound") > 1
      ) {
        plot.data.cpd = pathview::node.map(
          cpd.data,
          node.data,
          node.types = "compound", 
          node.sum = node.sum
        )
        if(map.cpdname & !kegg.native) {
          plot.data.cpd$labels = pathview::cpdkegg2name(plot.data.cpd$labels)[, 2]
          mapped.cnodes = rownames(plot.data.cpd)
          node.data$labels[mapped.cnodes] = plot.data.cpd$labels
        }
        cols.ts.cpd = pathview::node.color(
          plot.data.cpd,
          limit$cpd, 
          bins$cpd,
          both.dirs = both.dirs$cpd,
          trans.fun = trans.fun$cpd, 
          discrete = discrete$cpd,
          low = low$cpd,
          mid = mid$cpd, 
          high = high$cpd,
          na.col = na.col
        )
      } else plot.data.cpd = cols.ts.cpd = NULL
      if(kegg.native) {
        pv.pars = my.keggview.native(
          plot.data.gene = plot.data.gene, 
          cols.ts.gene = cols.ts.gene, 
          plot.data.cpd = plot.data.cpd, 
          cols.ts.cpd = cols.ts.cpd,
          node.data = node.data, 
          pathway.name = pathway.name[i],
          kegg.dir = kegg.dir, 
          limit = limit,
          bins = bins,
          both.dirs = both.dirs, 
          discrete = discrete,
          low = low,
          mid = mid,
          high = high, 
          na.col = na.col,
          ...
        )
      } else {
        pv.pars = pathview::keggview.graph(
          plot.data.gene = plot.data.gene, 
          cols.ts.gene = cols.ts.gene,
          plot.data.cpd = plot.data.cpd, 
          cols.ts.cpd = cols.ts.cpd,
          node.data = node.data, 
          path.graph = gR1,
          pathway.name = pathway.name[i], 
          map.cpdname = map.cpdname,
          split.group = split.group, 
          limit = limit,
          bins = bins,
          both.dirs = both.dirs, 
          discrete = discrete,
          low = low,
          mid = mid,
          high = high, 
          na.col = na.col,
          ...
        )
      }
      plot.data.gene = cbind(plot.data.gene, cols.ts.gene)
      if(!is.null(plot.data.gene)) {
        cnames = colnames(plot.data.gene)[-(1:8)]
        nsamp = length(cnames)/2
        if(nsamp > 1) {
          cnames[(nsamp + 1):(2 * nsamp)] = paste(
            cnames[(nsamp + 1):(2 * nsamp)], "col", sep = "."
          )
        } else cnames[2] = "mol.col"
        colnames(plot.data.gene)[-(1:8)] = cnames
      }
      plot.data.cpd = cbind(plot.data.cpd, cols.ts.cpd)
      if(!is.null(plot.data.cpd)) {
        cnames = colnames(plot.data.cpd)[-(1:8)]
        nsamp = length(cnames)/2
        if(nsamp > 1) {
          cnames[(nsamp + 1):(2 * nsamp)] = paste(
            cnames[(nsamp + 1):(2 * nsamp)], "col", sep = "."
          )
        } else cnames[2] = "mol.col"
        colnames(plot.data.cpd)[-(1:8)] = cnames
      }
      out.list[[i]] = list(
        plot.data.gene = plot.data.gene, 
        plot.data.cpd = plot.data.cpd
      )
    }
    if (npath == 1) 
      out.list = out.list[[1]]
    else names(out.list) = pathway.name
    return(invisible(out.list))
  }# <environment: namespace:pathview>
  my.keggview.native <- function (
    plot.data.gene = NULL,
    plot.data.cpd = NULL,
    cols.ts.gene = NULL, 
    cols.ts.cpd = NULL,
    node.data,
    pathway.name,
    out.suffix = "pathview", 
    kegg.dir = ".",
    multi.state = TRUE,
    match.data = TRUE,
    same.layer = TRUE, 
    res = 400,
    cex = 0.25,
    discrete = list(gene = FALSE, cpd = FALSE), 
    limit = list(gene = 1, cpd = 1),
    bins = list(gene = 10, cpd = 10), 
    both.dirs = list(gene = T, cpd = T),
    low = list(gene = "green", cpd = "blue"),
    mid = list(gene = "gray", cpd = "gray"), 
    high = list(gene = "red", cpd = "yellow"),
    na.col = "transparent", 
    new.signature = TRUE,
    plot.col.key = TRUE,
    key.align = "x", 
    key.pos = "topright", 
    ...
  ) {
    img <- png::readPNG(
      paste(kegg.dir, "/", pathway.name, ".png", sep = "")
    )
    width <- ncol(img)
    height <- nrow(img)
    cols.ts.gene = cbind(cols.ts.gene)
    cols.ts.cpd = cbind(cols.ts.cpd)
    nc.gene = max(ncol(cols.ts.gene), 0)
    nc.cpd = max(ncol(cols.ts.cpd), 0)
    nplots = max(nc.gene, nc.cpd)
    pn.suffix = colnames(cols.ts.gene)
    if (length(pn.suffix) < nc.cpd) 
      pn.suffix = colnames(cols.ts.cpd)
    if (length(pn.suffix) < nplots) 
      pn.suffix = 1:nplots
    if (length(pn.suffix) == 1) {
      pn.suffix = out.suffix
    }
    else pn.suffix = paste(out.suffix, pn.suffix, sep = ".")
    na.col = colorpanel2(1, low = na.col, high = na.col)
    if ((match.data | !multi.state) & nc.gene != nc.cpd) {
      if (nc.gene > nc.cpd & !is.null(cols.ts.cpd)) {
        na.mat = matrix(na.col, ncol = nplots - nc.cpd, nrow = nrow(cols.ts.cpd))
        cols.ts.cpd = cbind(cols.ts.cpd, na.mat)
      }
      if (nc.gene < nc.cpd & !is.null(cols.ts.gene)) {
        na.mat = matrix(
          na.col,
          ncol = nplots - nc.gene, 
          nrow = nrow(cols.ts.gene)
        )
        cols.ts.gene = cbind(cols.ts.gene, na.mat)
      }
      nc.gene = nc.cpd = nplots
    }
    out.fmt = "Working in directory %s"
    wdir = getwd()
    out.msg = sprintf(out.fmt, wdir)
    message("Info: ", out.msg)
    out.fmt = "Writing image file %s"
    multi.state = multi.state & nplots > 1
    if (multi.state) {
      nplots = 1
      pn.suffix = paste(out.suffix, "multi", sep = ".")
      if (nc.gene > 0) 
        cols.gene.plot = cols.ts.gene
      if (nc.cpd > 0) 
        cols.cpd.plot = cols.ts.cpd
    }
    for (np in 1:nplots) {
		  img.file = paste(
        kegg.dir,
        "/",
        pathway.name,
        ".",
        pn.suffix[np],
        ".png", 
			  sep = ""
      )
      out.msg = sprintf(out.fmt, img.file)
      message("Info: ", out.msg)
      png(img.file, width = width, height = height, res = res)
      op = par(mar = c(0, 0, 0, 0))
      plot(
        c(0, width),
        c(0, height),
        type = "n",
        xlab = "", 
        ylab = "",
        xaxs = "i",
        yaxs = "i"
      )
      if (new.signature) 
        img[height - 4:25, 17:137, 1:3] = 1
      if (same.layer != T) 
        rasterImage(img, 0, 0, width, height, interpolate = F)
      if (!is.null(cols.ts.gene) & nc.gene >= np) {
        if (!multi.state) 
          cols.gene.plot = cols.ts.gene[, np]
        if (same.layer != T) {
          render.kegg.node(
            plot.data.gene,
            cols.gene.plot, 
            img,
            same.layer = same.layer,
            type = "gene", 
            cex = cex
          )
        } else {
          img = render.kegg.node(
            plot.data.gene,
            cols.gene.plot, 
            img,
            same.layer = same.layer,
            type = "gene"
          )
        }
      }
      if (!is.null(cols.ts.cpd) & nc.cpd >= np) {
        if (!multi.state) 
          cols.cpd.plot = cols.ts.cpd[, np]
        if (same.layer != T) {
          render.kegg.node(
            plot.data.cpd,
            cols.cpd.plot, 
            img,
            same.layer = same.layer,
            type = "compound", 
            cex = cex
          )
        } else {
          img = render.kegg.node(
            plot.data.cpd,
            cols.cpd.plot, 
            img,
            same.layer = same.layer,
            type = "compound"
          )
        }
      }
      if (same.layer == T) 
        graphics::rasterImage(img, 0, 0, width, height, interpolate = F)
      pv.pars = list()
      pv.pars$gsizes = c(width = width, height = height)
      pv.pars$nsizes = c(46, 17)
      pv.pars$op = op
      pv.pars$key.cex = 2 * 72/res
      pv.pars$key.lwd = 1.2 * 72/res
      pv.pars$sign.cex = cex
      off.sets = c(x = 0, y = 0)
      align = "n"
      ucol.gene = unique(as.vector(cols.ts.gene))
      na.col.gene = ucol.gene %in% c(na.col, NA)
      if (plot.col.key & !is.null(cols.ts.gene) & !all(na.col.gene)) {
        off.sets = pathview::col.key(
          limit = limit$gene,
          bins = bins$gene, 
          both.dirs = both.dirs$gene,
          discrete = discrete$gene, 
          graph.size = pv.pars$gsizes,
          node.size = pv.pars$nsizes, 
          key.pos = key.pos,
          cex = pv.pars$key.cex,
          lwd = pv.pars$key.lwd, 
          low = low$gene,
          mid = mid$gene,
          high = high$gene, 
          align = "n"
        )
        align = key.align
      }
      ucol.cpd = unique(as.vector(cols.ts.cpd))
      na.col.cpd = ucol.cpd %in% c(na.col, NA)
      if (plot.col.key & !is.null(cols.ts.cpd) & !all(na.col.cpd)) {
        off.sets = pathview::col.key(
          limit = limit$cpd,
          bins = bins$cpd, 
          both.dirs = both.dirs$cpd,
          discrete = discrete$cpd, 
          graph.size = pv.pars$gsizes,
          node.size = pv.pars$nsizes, 
          key.pos = key.pos,
          off.sets = off.sets,
          cex = pv.pars$key.cex, 
          lwd = pv.pars$key.lwd,
          low = low$cpd,
          mid = mid$cpd, 
          high = high$cpd,
          align = align
        )
      }
      if (new.signature) 
        pathview.stamp(x = 17, y = 20, on.kegg = T, cex = pv.pars$sign.cex)
      par(pv.pars$op)
      dev.off()
    }
    return(invisible(pv.pars))
  }

  # Modify function in a package, change namespace
  # http://stackoverflow.com/questions/23279904/modifying-an-r-package-function-for-current-r-session-assigninnamespace-not-beh
  tmpfun <- get("keggview.native", envir = asNamespace("pathview"))
  environment(my.keggview.native) <- environment(tmpfun)
  # Don't know if this is really needed	
  attributes(my.keggview.native) <- attributes(tmpfun)  


	if (is.null(select_contrast)) {
    return(blank)
  }
	
	if(sig_pathways == "All") {
    return(blank)
  }
	
	if(length(limma$top_genes) == 0) {
    return(blank)
  }

	# Get fold change
	if(length(limma$comparisons)  == 1) {
    top_1 <- limma$top_genes[[1]]  
	} else {
	  top <- limma$top_genes
	  ix <- match(select_contrast, names(top))
	  if(is.na(ix)) {
      return (blank)
    }
	  top_1 <- top[[ix]] 
	}
	if(dim(top_1)[1] == 0 ) {
    return(blank)
	}
  
	 
	colnames(top_1) <- c("Fold","FDR")
  Species <- converted$species[1,1]
	  
	fold <- top_1[, 1]
  names(fold) <- rownames(top_1)
	fold <- convert_ensembl_to_entrez(
	  fold,
	  Species,
	  idep_data$org_info
	)
  kegg_species_id <- idep_data$kegg_species_id
	 
  kegg_species <- as.character(
    kegg_species_id[which(kegg_species_id[, 1] == Species), 3]
  )
	 
  if(nchar(kegg_species) <= 2) {
    return(blank)
  }

	 
	path_id <- kegg_pathway_id(
    sig_pathways,
    Species,
    "KEGG",
    select_org,
    idep_data$gmt_files,
    idep_data$org_info, 
    idep_data
  )
  
  # Kegg pathway id not found.
  if(is.null(path_id)) {
    return(blank)
  }
	if(nchar(path_id) < 3) {
    return(blank)
  }
  random_string <- gsub(".*file", "", tempfile()) 
	temp_folder <- tempdir() 
	out_file <- paste(
    temp_folder,
    "/",
    path_id,
    ".",
    random_string,
    ".png",
    sep = ""
  )
	
	pv.out <- mypathview(
    gene.data = fold,
    pathway.id = path_id,
    kegg.dir = temp_folder,
    out.suffix = random_string,
    species = kegg_species,
    kegg.native = TRUE
  )
	
  # Return a list containing the filename
  list(
    src = out_file,
    contentType = 'image/png',
    width = "100%",
    height = "100%",
    alt = "KEGG pathway image."
  )
}