#' fct_06_deg1.R This file holds all of the main data analysis functions
#' associated with sixth tab of the iDEP website.
#'
#'
#' @section fct_06_deg1.R functions:
#' \code{change_names}
#'
#'
#' @name fct_06_deg1.R
NULL


#' Change comparison names 
#'
#' Change comparison names in limma from "KO_ko-WT_ko" to    "KO-WT_for_ko"
#'
#' @param comparison Comparison to change the name for
#'
#' @export
#' @return The altered comparison string, see the description for example.
change_names <- function(comparison) {
  # check to see if work needs to be done
  # if no work needs to be done return input
  if (grepl(".*_.*-.*_.*", comparison)) {
    comparison_levels <- unlist(strsplit(comparison, "-"))
    comparison_levels <- unlist(strsplit(comparison_levels, "_"))
    if (length(comparison_levels) != 4) {
      comparison <- comparison_levels
    } else {
      comp_dup_index <- which(duplicated(comparison_levels))
      comparison <- paste0(
        comparison_levels[-c(comp_dup_index, comp_dup_index - 2)],
        collapse = "-"
      )
      comparison <- paste0(
        comparison,
        "_for_",
        comparison_levels[comp_dup_index]
      )
    }
  }
  return(comparison)
}

#' Create choices and title for model factors
#' 
#' Create a list of options for comparisons from the 
#' sample_info. Also create a title to be used in a 
#' UI checkbox based off the date file format and the
#' DEG method. The list of options will become
#' checkboxes for the user to create their experiment
#' design from.
#' 
#' @param sample_info Experiment file information for grouping
#' @param data_file_format Type of data being examined
#' @param counts_deg_method Method of DEG being performed (See
#'  DEG UI for options)
#' 
#' @export
#' @return A list containing a string title and a vector of
#'  comparisons to choose from.
list_factors_ui <- function(
  sample_info,
  data_file_format,
  counts_deg_method
) {
  if(is.null(sample_info)) {
    return(
      HTML(
        "<font size = \"2\">A <a href=\"https://idepsite.wordpress.com/data-format/\">
        sample information file</a> can be uploaded to build a linear model according
        to experiment design. </font>"
      )
    ) 
  } else {
		factors <- colnames(sample_info)
    choices <- setNames(factors, factors)
		title <- "1. Select 1 or 2 main factors. Or leave it blank and just choose pairs
              of sample groups below."	
    if (data_file_format == 1 & counts_deg_method==3) {
      title <- "1. Select 6 or less main factors. Or skip this step and just choose 
                pairs of sample groups below."
    }
    return(list(
      title = title,
      choices = choices
    ))	   	  
	}
}

#' Get the block factor choices
#' 
#' This function uses the sample_info file and the selected
#' factors for the model to create a selection for the batch
#' effect. Returns a vector that turns into a checkbox for
#' the User.
#' 
#' @param sample_info Experiment file information for grouping
#' @param select_factors_model The selected factors for the model
#'  expression
#' 
#' @export
#' @return This function returns a vector of choices for a batch
#'  effect or paired samples.
list_block_factors_ui <- function(
  sample_info,
  select_factors_model
) {
  if (is.null(sample_info)) {
		return(NULL)		   
	} else { 
		factors <- colnames(sample_info)
		factors <- setdiff(factors, select_factors_model)

		if (length(factors) == 0 ) {
      return(NULL)
    }

		choices = setNames(factors, factors)	

    return(choices)
	}
}

#' Create model comparisons choices
#' 
#' This function uses the sample_info file and the selected
#' factors to create a list of options for model comparisons.
#' Changes with the input of select_factors_model. If there
#' is no selected factor then it defaults to comparisons that
#' can be created from the processed data.
#' 
#' @param sample_info Experiment file information for grouping
#' @param select_factors_model The selected factors for the model
#'  expression
#' @param processed_data Data that has been through the pre-processing
#'  function
#' 
#' @export
#' @return Returns a list containing a vector of choices and a
#'  title for the UI element.
list_model_comparisons_ui <- function(
  sample_info,
  select_factors_model,
  processed_data
) {
	if (is.null(sample_info) | is.null(select_factors_model)) { 
    factors <- as.character (
      detect_groups(
        colnames(processed_data)
      )
    )
		factors <- unique(factors)
		comparisons <- apply(
      t(combn(factors, 2)),
      1,
      function(x) paste(x, collapse = " vs. ")
    )
	  comparisons <- c(
      comparisons,
      apply(
        t(combn(rev(factors),2)),
        1,
        function(x) paste(x, collapse = " vs. ")
      )
    )	
		comparisons <- sort(comparisons)
		choices <- stats::setNames(gsub(" vs\\. ","-",comparisons), comparisons)
    title <- "Select comparisons among sample groups:"

    return(list(
      choices = choices,
      title = title
    ))
  } else { 
		choices = list()
				
    for(selected_factors in select_factors_model) { 
			ix = match(selected_factors, colnames(sample_info))
				
      if(is.na(ix)) {
        next
      }

			factors <- unique(sample_info[, ix])
			comparisons <- apply(
        t(combn(factors, 2)),
        1,
        function(x) paste(x, collapse = " vs. ")
      )
				comparisons <- c(
          comparisons,
          apply(
            t(combn(rev(factors), 2)),
            1,
            function(x) paste(x, collapse = " vs. ")
          )
        )	
				comparisons <- sort(comparisons)
				comparisons <- paste0(selected_factors, ": ", comparisons)
				choices <- append(choices, stats::setNames(comparisons, comparisons))
        title <- "2. Select one or more comparisons:"
		}
			
    if(length(choices) == 0) {
      return(NULL)
    } else {
      return(list(
        choices = choices,
        title = title
      ))
    }
	} 
}

#' List model interaction terms
#' 
#' This functions uses the sample info file and the selected
#' model factors to create interaction terms to be used in the 
#' DEG process. 
#' 
#' @param sample_info Experiment file information for grouping
#' @param select_factors_model The selected factors for the model
#'  expression
#' 
#' @export
#' @return Returns a character string of an interaction term
#'  between the selected factors. Used in a checkbox for the
#'  User to create a model expression
list_interaction_terms_ui <- function(
  sample_info,
  select_factors_model
) {
	if (is.null(sample_info) | is.null(select_factors_model)) {
    return(NULL)
	}	 else { 
		selected_factors = select_factors_model[!grepl(":", select_factors_model)]
    
    if(length(selected_factors) <= 1) {
      return(NULL) 
    }
		interactions <- apply(
      t(combn(selected_factors, 2)),
      1,
      function(x) paste(x, collapse = ":")
    )
    return(interactions)
  }
}

#' Create a string of the model design
#' 
#' Use the model design selections to create a string of the
#' model design being used for the DEG analysis.
#' 
#' @param sample_info Experiment file information for grouping
#' @param select_factors_model The selected factors for the model
#'  expression
#' @param select_block_factors_model The selected factors for
#'  batch effect
#' @param select_interactions The interaction terms being used in
#'  the model design
#' 
#' @export
#' @return Returns a string of the model design being used for
#'  the DEG analysis
experiment_design_txt <- function(
  sample_info,
  select_factors_model,
  select_block_factors_model,
  select_interactions
) {
  if (is.null(sample_info) | is.null(select_factors_model)) {
    return(NULL)
	} else {
    model <- paste(
      "Model: expression ~ ",
      paste(select_factors_model, collapse = " + ")
    )
		if(!is.null(select_block_factors_model)) {
      model <- paste0(
        model,
        " + ",
        paste(select_block_factors_model, collapse = " + " )
      )
    }
    if(!is.null(select_interactions)) {
      model <- paste0(
        model,
        " + ",
        paste(select_interactions, collapse = " + ")
      )
    }
    
    return(model)									
	}
}

#' Refernce levels for selected factor
#' 
#' This function uses a vector of selected factors to create
#' choices for the reference level to use for the factor in
#' the DEG analysis. 
#' 
#' @param sample_info Experiment file information for grouping
#' @param select_factors_model The selected factors for the model
#'  expression
#' @param data_file_format Type of gene data being examined
#' @param counts_deg_method The method or package being used for
#'  the DEG analysis
#' 
#' @export
#' @return A list the same length as the vector of selected factors.
#'  Each entry in the list corresponds to the choice of group to
#'  use for the reference level.
select_reference_levels_ui <- function(
  sample_info,
  select_factors_model,
  data_file_format,
  counts_deg_method
) {
  if (is.null(sample_info) | is.null(select_factors_model)) {
    return(NULL)
  }	else {
    selected_factors <- select_factors_model[!grepl(":", select_factors_model)]
    
    if(length(selected_factors) == 0) {
      return(NULL)
    }
    if (!(data_file_format == 1 & counts_deg_method == 3)) {
      return(NULL)
    }
    select_choices <- c()
    for(i in selected_factors){
      if(is.na(match(i, colnames(sample_info)))) {
        select_choices[[i]] <- NULL
      } else {
        select_choices[[i]] <- unique(
          sample_info[, i]
        )      
      }
    }
    return(select_choices)	
	}
}

#' DEG analysis function
#' 
#' Use the limma or DESeq2 package to perform DEG analysis
#' with the specified model design. Core function for the
#' DEG panel of iDEP.
#' 
#' @param data_file_format Type of gene data being examined
#' @param counts_deg_method The method or package being used for
#'  the DEG analysis
#' @param raw_counts The matrix of counts before processing for
#'  gene expression data
#' @param limma_p_val Significant p-value to use for expressed
#'  genes
#' @param limma_fc Minimum fold-change cutoff for the DEG
#'  analysis
#' @param select_model_comprions Selected comparisons to analyze
#'  in the DEG analysis
#' @param sample_info Experiment file information for grouping
#' @param select_factors_model The selected factors for the model
#'  expression
#' @param select_interactions The interaction terms being used in
#'  the model design
#' @param select_block_factors_model The selected factors for
#'  batch effect
#' @param factor_reference_levels Vector of reference levels to
#'  use for the selected factors
#' @param processed_data Data that has been through the pre-processing
#' @param counts_log_start The constant added to the log transformation
#'  from pre-processing
#' @param p_vals The vector of p-vals calculated in pre-process for
#'  significant expression
#' 
#' @export
#' @return List with the results of the DEG analysis. When the function
#'  is successful there are four entries in the list. "results" is a
#'  matrix with the same dimensions as the processed data. The entries
#'  in "results" are c(-1, 0, 1) for (negative fold change, no significant
#'  change, positive fold change) respectively. The second entry is
#'  "comparisons" and is a character vector of the different comparisons
#'  that were analyzed in the function. Third is "exp_type" and details
#'  the model expression that was used for the DEG analysis. Lastly is
#'  "top_genes" which is itself a list. The "top_genes" list has an entry
#'  for each comparison. Each entry is a data frame with two columns. One
#'  column is the calculated fold change for the comparison and the other
#'  is the adjusted p-value for the fold change calculation.
limma_value <- function(
  data_file_format,
  counts_deg_method,
  raw_counts,
  limma_p_val,
  limma_fc,
  select_model_comprions,
  sample_info,
  select_factors_model,
  select_interactions,
  select_block_factors_model,
  factor_reference_levels,
  processed_data,
  counts_log_start,
  p_vals
) {
  if(data_file_format == 1) {
    if(counts_deg_method == 3) {
      return(
        deg_deseq2(
          raw_counts = raw_counts,
          max_p_limma = limma_p_val,
          min_fc_limma = limma_fc,
				  selected_comparisons = select_model_comprions,
          sample_info = sample_info,
					model_factors = c(select_factors_model, select_interactions), 
				  block_factor = select_block_factors_model,
          reference_levels = factor_reference_levels
        )
      )
		} else if(counts_deg_method < 3) {
      return(
        deg_limma(
          processed_data = processed_data,
          max_p_limma = limma_p_val,
          min_fc_limma = limma_fc,
					raw_counts = raw_counts,
          counts_deg_method = counts_deg_method,
					prior_counts = counts_log_start,
          data_file_format = data_file_format,
					selected_comparisons = select_model_comprions,
          sample_info = sample_info,
					model_factors = c(select_factors_model, select_interactions),
					block_factor = select_block_factors_model
        )
      )
    }
  } else if(data_file_format == 2) {
	  return(
      deg_limma(
        processed_data = processed_data,
        max_p_limma = limma_p_val,
        min_fc_limma = limma_fc,
				raw_counts = raw_counts,
        counts_deg_method = counts_deg_method,
				prior_counts = counts_log_start,
        data_file_format = data_file_format,
				selected_comparisons = select_model_comprions,
        sample_info = sample_info,
				model_factors = c(select_factors_model, select_interactions),
				block_factor = select_block_factors_model,
      )
    )
	} else {
		if(!is.null(p_vals)) {
		  ix <- match(rownames(processed_data), rownames(p_vals))
		  p_vals <- p_vals[ix, ]
		}

		# Looks like ratio data, take log2
		if(sum(
      round(apply(processed_data, 2, median) + .2) == 1
    ) == dim(x)[2] & min(x) > 0) {
      processed_data <- log2(processed_data)
    }
		
		exp_type <- "None standard data without replicates."
		all_calls <- processed_data
		for(i in 1:dim(all_calls)[2]) { 
			tem <- all_calls[, i]
			all_calls[which(
        tem <= log2(limma_fc) & tem >= -log2(limma_fc)
      ) , i] <- 0			
			all_calls[which(tem > log2(limma_fc)), i] <- 1
			all_calls[which(tem < -log2(limma_fc)), i] <- -1		
			if(!is.null(p_vals)) {
        all_calls[which(p_vals[, i] > limma_p_val), i] <- 0
      } 
		}
		comparisons <- colnames(all_calls)
		extract_column <- function(
      i,
      processed_data,
      p_vals,
      top_genes
    ) {
			top_genes <- as.data.frame(processed_data[, i, drop = FALSE])
			if(is.null(p_vals)) {
        top_genes$FDR <- 0
      } else {
        top_genes$FDR <- p_vals[, i]
      }
			colnames(top_genes) <- c("Fold","FDR")
			return(top_genes)	
		} 
		top_genes <- lapply(1:dim(processed_data)[2], function(x) {
      extract_column(
        i = x,
        processed_data = processed_data,
        p_vals = p_vals,
        top_genes = top_genes
      )
    })
		top_genes <- setNames(top_genes, colnames(processed_data))
		
		return(list(
      results = all_calls,
      comparisons = colnames(processed_data),
      exp_type = exp_type,
      top_genes = top_genes
    ))
  }
}

#' Differential expression using DESeq2 package
#' 
#' Used in the limma_value function to perform DEG analysis using the
#' DESeq2 package. It is not recommended to use this function on its own.
#' 
#' @param raw_counts The matrix of counts before processing for
#'  gene expression data
#' @param max_p_limma Significant p-value to use for the fold-change
#'  values
#' @param min_fc_limma Minimum fold-change to include in the results
#' @param selected_comparisons Comparisons being analyzed in the DEG
#'  analysis
#' @param sample_info Experiment file information for grouping
#' @param model_factors Vector of selected factors and interaction terms
#'  from the model design
#' @param block_factor The selected factors for batch effect
#' @param reference_levels Vector of reference levels to use for the
#'  selected factors
#' 
#' @export
#' @return The return value is the results of the DEG analysis. These
#'  results are filtered and formatted by the limma_value function.
deg_deseq2 <- function(
  raw_counts,
  max_p_limma = .05,
  min_fc_limma = 2,
  selected_comparisons = NULL,
  sample_info = NULL,
  model_factors = NULL,
  block_factor = NULL,
  reference_levels = NULL
){
  max_samples <- 100
  max_comparisons <- 20

	groups <- as.character(detect_groups(colnames(raw_counts), sample_info))
	unique_groups <- unique(groups)	
	
	# Check for replicates, removes samples without replicates
  # Number of replicates per biological sample
	reps <- as.matrix(table(groups))
  # Less than 2 samples with replicates
	if (sum(reps[, 1] >= 2) < 2) {
    return(list(
      results= NULL,
      comparisons = NULL,
      exp_type = 
        "Failed to parse sample names to define groups. Cannot perform DEGs
         and pathway analysis. Please double check column names! Use 
         WT_Rep1, WT_Rep2 etc. ",
      top_genes = NULL
    ))
  }
	# Remove samples without replicates
	unique_groups <- rownames(reps)[which(reps[, 1] > 1)]
	ix <- which(groups %in% unique_groups)  
	groups <- groups[ix]   
	raw_counts <- raw_counts[, ix] 
	exp_type <- paste(length(unique_groups)," sample groups detected.")	
	
	# Too many samples 
	if(ncol(raw_counts)  > max_samples) { 
		return(list(
      results = NULL,
      comparisons = NULL, 
			exp_type = paste(
        exp_type,
        "Too many samples for DESeq2. Please choose limma-voom or 
         limma-trend."
      ),
			top_genes = NULL
    ))
	}	
		
	# All pair-wise comparisons
	comparisons <- ""
	for(i in 1:(length(unique_groups) - 1)) {
    for (j in (i + 1):length(unique_groups)){
      comparisons <- c(
        comparisons,
        paste(unique_groups[j], "-", unique_groups[i], sep = "")
      )
    }
  }
	comparisons <- comparisons[-1]

   # Too many comparisons 
	if(length(comparisons)  > max_comparisons) { 
		exp_type = paste(
      exp_type,
      " Too many comparisons. Only first",
      max_comparisons,
      "of the ",
      length(comparisons), 
			"comparisons calculated. Please choose comparisons."
    )
		comparisons <- comparisons[1:max_comparisons]
	}	
	
	col_data <- cbind(colnames(raw_counts), groups)

	# No sample file, but user selected comparisons using column names
	if(is.null(model_factors) & length(selected_comparisons) > 0) {
		comparisons <- selected_comparisons
  }

	comparison_names <- comparisons
	# Set up the DESeqDataSet Object and run the DESeq pipeline
	dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = col_data,
    design = ~groups
  )								

	if(is.null(model_factors)) {
    dds = DESeq2::DESeq(dds)
  } else {
    # Using selected factors and comparisons ----------
		# Build model
    # Block factor is just added in
		model_factors <- c(model_factors, block_factor)  
    
    # Selected factors and interactions:
    # c( "strain", "treatment", "strain:treatment")
		factors <- model_factors  
    # Non-interaction terms
		factors <- factors[!grepl(":", factors)]
		# Interaction terms like strain:treatment
		interactions <- model_factors[grepl(":", model_factors)]
		
		col_data <- sample_info  
    # Factors are encoded as "A", "B", "C"; Avoids illegal letters
		factors_coded <- toupper(letters)[1: dim(col_data)[2]]
    # For look up; each column of sample_info
		names(factors_coded) <- colnames(col_data)
    # All columns named A B C D  
		colnames(col_data) <- factors_coded  

		col_data <- as.data.frame(col_data)
		
		# Set reference levels for factors
    # c("genotype:wt", "treatment:control")
		if(! is.null(reference_levels) ) {
			# First factor
			for(refs in reference_levels) {
        if(!is.null(refs)) {
					# Corresponding column id for factor
          ix <- match(
            gsub(":.*", "", refs),
            colnames(sample_info)
          )
					col_data[, ix] <- as.factor(col_data[, ix])
					col_data[, ix] <- relevel(
            col_data[, ix],
            gsub(".*:", "", refs)
          )
				}
      }
		}
				
		# Base model
    deseq2_object <- paste(
      "dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts,
      colData = col_data, design = ~ ", 
			paste(factors_coded[factors], collapse = "+")
    )	
		exp_type <- paste(
      "Model: ~",
      paste(model_factors, collapse = " + ")
    )

		
		# Create model
		if(length(interactions) > 0) {
			for(interaction_terms in interactions) {
				# Split strain:mutant as "strain" and "mutant"
        interacting_factors <- unlist(
          strsplit(interaction_terms, ":")
        )
        # Convert "strain:mutant" to "A:B"
				tem <- paste(
          factors_coded[interacting_factors],
          collapse = ":"
        )   
				deseq2_object <- paste(deseq2_object, " + ", tem)
			}			
		}

    # End the model
		deseq2_object <- paste(deseq2_object, ")") 

		eval(parse(text = deseq2_object))
	
		dds = DESeq2::DESeq(dds)  # main function		


		# Comparisons 
		# "group: control vs. mutant"
		comparisons <- gsub(".*: ", "", selected_comparisons)
		comparisons <- gsub(" vs\\. ", "-", comparisons)
    # Corresponding factors for each comparison
		factors_vector <- gsub(":.*", "", selected_comparisons)
		
		# comparison_names holds names for display with real factor names
		# comparisons is used in calculation it is A, B, C for factors
		comparison_names <- comparisons

		# Note that with interaction terms, not all meaningful comparisons is
    # listed for selection. This is complex. Only under reference level.
		
		# Comparisons due to interaction terms
		if(length(interactions) > 0) {
			interaction_comparisons <- DESeq2::resultsNames(dds)
			interaction_comparisons <- interaction_comparisons[grepl(
        "\\.", interaction_comparisons
      )]
	
			comparisons <- c(comparisons, interaction_comparisons)
		
			# Translate comparisons generated in interaction terms back to factor names
			interaction_comparison_names <- interaction_comparisons
			for(i in 1:length(interaction_comparison_names)) {
				tem <- unlist(strsplit(interaction_comparison_names[i], "\\."))
				tem_factors <- substr(tem, 1, 1) 
				
        # Get the first letter and translate into real factor names
				tem_factors[1] <- names(factors_coded)[factors_coded == tem_factors[1]]
        # Get the 2nd letters and translate into real factor names  
				tem_factors[2] <- names(factors_coded)[factors_coded == tem_factors[2]]  

				interaction_comparisons[i] <- paste0(
          "I:",
          tem_factors[1],
          "_",
          substr(tem[1], 2, nchar(tem[1])),
          ".",
          tem_factors[2],
          "_",
          substr(tem[2], 2, nchar(tem[2])) 
			  )				
			}
			comparison_names = c(comparison_names,interaction_comparisons)
		}
	}
	
	# Extract contrasts according to comparisons defined above
	result_first <- NULL
  all_calls <- NULL
	top_genes <- list()
  # Counter
  pk <- 1
  # First results
	pp <- 0 

	for(kk in 1:length(comparisons)) {
		tem = unlist(strsplit(comparisons[kk], "-"))
		
    # Group comparison using sample names
		if(is.null(model_factors)) {
      selected <- DESeq2::results(dds, contrast = c("groups", tem[1], tem[2]) )
    } else {
      # Not interaction term: they contain . interaction term
			if(!grepl("\\.", comparisons[kk])) {
        selected <- DESeq2::results(
          dds,
          contrast = c(factors_coded[factors_vector[kk]], tem[1], tem[2])
        )
      # either A, B, C ...
      } else {
        # Interaction term
        selected <- DESeq2::results(dds, name = comparisons[kk])
      }
		}

		selected$calls <- 0

		selected$calls[which(
      selected$log2FoldChange > log2(min_fc_limma) & selected$padj < max_p_limma
      )]  <-  1
		
    selected$calls[which(
      selected$log2FoldChange < -log2(min_fc_limma) & selected$padj < max_p_limma
    )] <- -1
		
    colnames(selected) <- paste(
      as.character(comparison_names[kk]),
      "___",
      colnames(selected),
      sep = ""
    )
		selected <- as.data.frame(selected)
    # First one with significant genes, collect gene list and Pval+ fold
		if(pp == 0) {
      result_first <- selected
      pp <- 1
		  top_genes[[1]] <- selected[, c(2, 6)] 
		  names(top_genes)[1] <- comparison_names[kk]
    } else {
      result_first <- merge(result_first, selected, by = "row.names") 
			rownames(result_first) <- result_first[, 1]
      result_first <- result_first[, -1]
      pk <- pk + 1
      top_genes[[pk]] <- selected[, c(2,6)]
      # Assign name to comparison
			names(top_genes)[pk] <- comparison_names[kk]
		}
	}

	interactions <- c()
	if(!is.null(model_factors)) {
    interactions <- model_factors[grepl(":", model_factors )]
  }
		
  # Add comprisons for non-reference levels. It adds to the results_first object.	
	if(length(interactions) > 0) {
    # Factor whose values are factors and names are factor and level combination
    # conditionTreated, genotypeWT
    factor_lookup <- c()
		level_lookup <- c()
		
		for(i in 1:dim(sample_info)[2]) {
      unique_sample_info <- unique(sample_info)
			tem <- rep(toupper(letters)[i], dim(unique_sample_info)[1])
			names(tem) <- paste0(toupper(letters)[i], unique_sample_info[, i])
			factor_lookup <- c(factor_lookup, tem)  
			
			tem <- as.character(unique_sample_info[, i])
			names(tem) <- paste0(toupper(letters)[i], unique_sample_info[, i])
			level_lookup <- c(level_lookup, tem)
		}
		
		# None interaction terms 
		none_inter_terms <- DESeq2::resultsNames(dds)[!grepl(
      "\\.", DESeq2::resultsNames(dds)
    )]
		none_inter_terms <- none_inter_terms[-1]
		all_interaction_terms <- DESeq2::resultsNames(dds)[grepl(
      "\\.", DESeq2::resultsNames(dds)
    )]

    # Each none interaction term
		for(kk in 1:length(none_inter_terms)) { 
			# Not just group comparison using sample names
      if(!is.null(model_factors)) {
				# Current factor
				c_factor <- gsub("_.*", "", none_inter_terms[kk])

				for(interaction_term in all_interaction_terms) {
          # 4 components
					splits <- split_interaction_terms(
            interaction_term,
            factor_lookup = factor_lookup,
            level_lookup = level_lookup
          )

					if (c_factor != splits[1] & c_factor != splits[3]) {
            next
          }

					selected <- DESeq2::results(
            dds,
            list(c(none_inter_terms[kk], interaction_term))
          ) 
					comparison_name <- paste0(
            none_inter_terms[kk],
            "__",
            gsub("\\.", "", interaction_term)
          )
						
					if(c_factor == splits[1]) {
            other_level <- splits[4]
          } else {
            other_level <- splits[2]
          }
							
					comparison_name = paste0(
            gsub(
              "_vs_",
              "-",
              substr(none_inter_terms[kk], 3, nchar(none_inter_terms[kk]))
            ), 
						"_for_",
            other_level
          )
					comparison_names <- c(comparison_names, comparison_name)
					selected$calls <- 0   
					selected$calls[which(
            selected$log2FoldChange > log2(min_fc_limma) & selected$padj < max_p_limma
          )] <- 1
					selected$calls[which(
            selected$log2FoldChange <  -log2(min_fc_limma) & selected$padj < max_p_limma
          )] <- -1
					colnames(selected) <- paste(comparison_name, "___", colnames(selected), sep = "")
					selected <- as.data.frame(selected)
          # First one with significant genes, collect gene list and Pval+ fold
					if(pp == 0) {
						result_first <- selected
            pp = 1 
						top_genes[[1]] <- selected[, c(2, 6)] 
						names(top_genes)[1] <- comparison_name
          } else {
            result_first = merge(result_first, selected, by = "row.names") 
						rownames(result_first) <- result_first[, 1] 
						result_first <- result_first[, -1]
						pk <- pk + 1 
						top_genes[[pk]] <- selected[, c(2, 6)]
            names(top_genes)[pk] <- comparison_name
					}
				}			
			} 
		}
	}

	if(!is.null(result_first)) { 
		# Note that when you only select 1 column from a data frame it automatically 
    # converts to a vector. drop = FALSE prevents that.
		all_calls <- as.matrix(
      result_first[, grep("calls", colnames(result_first)), drop = FALSE]
    )
		colnames(all_calls) <- gsub("___.*", "", colnames(all_calls))
    # Note that samples names should have no "."
		colnames(all_calls) <- gsub("\\.", "-", colnames(all_calls))
		colnames(all_calls) <- gsub("^I-", "I:", colnames(all_calls))
	}

	return(list(
    results= all_calls,
    comparisons = comparison_names,
    exp_type = exp_type,
    top_genes = top_genes
  )) 
}

#' SPLIT INTERACION TERMS
#' Split  genotypeI.conditionTrt --> c("genotype","I","conditoin","Trt")
split_interaction_terms <- function(
  term,
  factor_lookup,
  level_lookup
) {
	if(!grepl("\\.", term)) {
    return(NULL)
  }
	terms_split <- unlist(strsplit(term, "\\."))
	# factor1, level1, factor2, level2
	return(
    c(
      factor_lookup[terms_split[1]],
      level_lookup[terms_split[1]],
      factor_lookup[terms_split[2]],
      level_lookup[terms_split[2]]
    )
  )
}


#' Differential expression using limma package
#' 
#' Used in the limma_value function to perform DEG analysis using the 
#' limma package. It is not recommended to use this function on its own.
#' 
#' @param processed_data Data that has been through the pre-processing
#' @param max_p_limma Significant p-value to use for the fold-change
#'  values
#' @param min_fc_limma Minimum fold-change to include in the results
#' @param raw_counts The matrix of counts before processing for
#'  gene expression data
#' @param counts_deg_method The method or package being used for
#'  the DEG analysis
#' @param prior_counts The constant added to the log transformation
#'  from pre-processing
#' @param data_file_format Type of gene data being examined
#' @param selected_comparisons Selected comparisons to analyze
#'  in the DEG analysis
#' @param sample_info Experiment file information for grouping
#' @param model_factors Vector of selected factors and interaction terms
#'  from the model design
#' @param block_factor The selected factors for batch effect
#' 
#' @export
#' @return The return value is the results of the DEG analysis. These
#'  results are filtered and formatted by the limma_value function.
deg_limma <- function(
  processed_data,
  max_p_limma = .1,
  min_fc_limma = 2,
  raw_counts,
  counts_deg_method,
  prior_counts,
  data_file_format,
  selected_comparisons = NULL,
  sample_info = NULL,
  model_factors = NULL,
  block_factor = NULL
){
  # Many different situations:
  # 1. Just use sample names
  # 2. Just one factor
  # 3. Two factors no interaction
	# 4. Two factors with interaction
  # 5. Block factor 

	top_genes <- list()
  limma_trend <- FALSE
	if(data_file_format == 2) {
    eset <- methods::new("ExpressionSet", exprs = as.matrix(processed_data))
  } else {
    # Limma-trend method selected for counts data
		if(counts_deg_method == 1) {
      # Use transformed data for limma  
			eset <- methods::new("ExpressionSet", exprs = as.matrix(processed_data))
			limma_trend = TRUE
		}
	}

	groups <- colnames(processed_data)
	groups <-  detect_groups(groups, sample_info) 
	unique_groups <- unique(groups)  
	
	# Check for replicates, removes samples without replicates
  # Number of replicates per biological sample
	reps <- as.matrix(table(groups))
  # Less than 2 samples with replicates
	if(sum(reps[, 1] >= 2) < 2) {
    return(
      list(
        results = NULL,
        comparisons = NULL,
        exp_type = 
          "Failed to parse sample names to define groups. Cannot perform
          DEGs and pathway analysis. Please double check column names! Use
          WT_Rep1, WT_Rep2 etc. ",
        topGenes=NULL
      )
    )
  }
	 
	# Remove samples without replicates
	unique_groups <- rownames(reps)[which(reps[, 1] > 1)]
	ix <- which(groups %in% unique_groups)  
	groups <- groups[ix]   
	processed_data <- processed_data[, ix]
  raw_counts <- raw_counts[, ix] 
	
  # Just two groups
	if(length(unique_groups) == 2) {  
		unique_groups <- unique(groups)
    # "Mutant-WT"
		comparisons <-  paste(unique_groups[2], "-", unique_groups[1], sep = "")
		
		# No sample file, but user selected comparisons using column names
		if(is.null(model_factors) & length(selected_comparisons) > 0) {
      comparisons <- selected_comparisons
    }

    # Set reference level based on the order in which the levels appear
    # The first appearing level is set as reference; otherwise, we get
    # up and down-regulation reversed.
    groups <- factor(groups, levels = unique_groups) 

		design <- model.matrix(~0 + groups)
		colnames(design) <- unique_groups
		
    # Voom
		if(!is.null(raw_counts) && counts_deg_method == 2) {
			dge <- edgeR::DGEList(counts = raw_counts)
      # Normalization
			dge <- edgeR::calcNormFactors(dge, method = "TMM")
			voom_results <- limma::voom(dge, design)
      fit <- limma::lmFit(v, design)
    } else {
      # Regular limma
			fit <- limma::lmFit(eset, design)
    }

    cont_matrix <- limma::makeContrasts(contrasts = comparisons, levels = design)
		contrasts_fit <- limma::contrasts.fit(fit, cont_matrix)
		fit <- limma::eBayes(contrasts_fit, trend = limma_trend)

		# Calls differential gene expression 1 for up, -1 for down
		results <- limma::decideTests(
      fit,
      p.value = max_p_limma,
      lfc = log2(min_fc_limma)
    )

		top_genes_table <- limma::topTable(fit, number = 1e12, sort.by = "M")
		if(dim(top_genes_table)[1] != 0) {
      top_genes_table <- top_genes_table[, c('logFC', 'adj.P.Val')]
      top_genes[[1]] <- top_genes_table
    }

		# Log fold change is actually substract of means. So if the data is natral log
    # transformed, it should be natral log.
		exp_type = "2 sample groups."
  }
	
  # More than two sample groups
	if(length(unique_groups) > 2) {
	  # Set reference level based on the order in which the levels appear
    # The first appearing level is set as reference; otherwise, we get up and 
    # down-regulation reversed.
    groups <- factor(groups, levels = unique_groups)

		design <- model.matrix(~ 0 + groups)
		colnames(design) <- gsub("^groups", "", colnames(design))

		if(!is.null(raw_counts) && counts_deg_method == 2) {
			limma_voom <- limma::voom(raw_counts, design) 
			fit <- limma::lmFit(limma_voom, design) 
		} else {
      fit <- limma::lmFit(eset, design)
    }
		
		fit <- limma::eBayes(fit, trend = limma_trend)
		
		comparisons <- ""
		for(i in 1:(length(unique_groups) - 1)) {
      for (j in (i + 1):length(unique_groups)) {
        comparisons <- c(
          comparisons,
          paste(unique_groups[j], "-", unique_groups[i], sep = "")
        )
      }
    }
		comparisons <- comparisons[-1]

		# No sample file, but user selected comparisons using column names
		if(is.null(model_factors) & length(selected_comparisons) > 0) {
      comparisons <- selected_comparisons
    }
		
		make_contrast <- limma::makeContrasts(contrasts = comparisons[1], levels = design)
		if(length(comparisons) > 1) {
      for(kk in 2:length(comparisons) ) {
        make_contrast <- cbind(
          make_contrast,
          limma::makeContrasts(contrasts = comparisons[kk], levels = design)
        )
      }
    }
		exp_type = paste(length(unique_groups), " sample groups detected.")	 
		
		# Factorial design 2x2, 2x3, 3x5 etc.
		# All samples must be something like WT_control_rep1
		if(sum(sapply(strsplit(unique_groups, "_"), length) == 2) == length(unique_groups)) {
			comparisons <- ""
			for(i in 1:(length(unique_groups) - 1)) {
        for (j in (i + 1):length(unique_groups)) {
          # Only compare WT_control vs. WT_treatment
          if(strsplit(unique_groups[i], "_")[[1]][1] == strsplit(unique_groups[j], "_")[[1]][1] |
             strsplit(unique_groups[i], "_")[[1]][2] == strsplit(unique_groups[j], "_")[[1]][2]) {
            comparisons <- c(comparisons, paste(unique_groups[j], "-", unique_groups[i], sep = ""))
          }
        }
      }
			comparisons <- comparisons[-1]

			# Factors genotype treatment levels
			extract_treatment <- function(x) {
        paste(gsub(".*_", "", unlist(strsplit(x, "-"))), collapse = "-")
      }
			extract_genotype <- function (x) {
        gsub("_.*", "", unlist(strsplit(x, "-")))[1]
      }
			extract_treatment_counting <- unique(gsub(".*_", "", unlist(strsplit(unique_groups, "-"))))
			treatments <- sapply(comparisons, extract_treatment)
			genotypes <- sapply(comparisons, extract_genotype)
			exp_type <- paste(
        exp_type,
        "\nFactorial design:",
        length(unique(genotypes)),
        "X",
        length(extract_treatment_counting),
        sep = ""
      )

			# Pairwise contrasts
			make_contrast <- limma::makeContrasts(
        contrasts = comparisons[1],
        levels = design
      )
			for(kk in 2:length(comparisons)) {
        make_contrast <- cbind(
          make_contrast,
          limma::makeContrasts(contrasts = comparisons[kk], levels = design)
        )
      }
			contrast_names = colnames(make_contrast)

			# Interaction contrasts
			for (kk in 1:(length(comparisons) - 1)) {
        for(kp in (kk + 1):length(comparisons)) {
          if(treatments[kp] == treatments[kk]) {  
					  make_contrast <- cbind(
              make_contrast,
              make_contrast[, kp] - make_contrast[, kk]
            )
					  contrast_names <- c(
              contrast_names,
              paste(
                "I:", 
                genotypes[kp],
                "-",
                genotypes[kk],
                "(",
                gsub("-", ".vs.", treatments[kp]),
                ")",
                sep = ""
              )
            )
				  }
        }   
			}

			colnames(make_contrast) <- contrast_names
			comparisons <- contrast_names
		}
		

		# Sample information is uploaded and user selected factors and comparisons
		if(!is.null(model_factors) & length(selected_comparisons) > 0) {
			exp_type <- paste("Model: ~", paste(model_factors, collapse = " + "))
      # Default value to be re-write if needed
			interaction_term <- FALSE
			
			# Model factors that does not contain interaction terms
			# model_factors "genotype", "condition", "genotype:condition"
			key_model_factors <- model_factors[!grepl(":", model_factors)]
			
			# "genotype: control vs. mutant"
      # Corresponding factors for each comparison
			factors_vector <- gsub(":.*", "", selected_comparisons) 
			# Remove factors not used in comparison, these are batch effects/pairs/blocks, 	
			
			# A factor is selected both in block and main factors, then use it as block factor
			key_model_factors <- key_model_factors[
        is.na(match(key_model_factors, block_factor))
      ]	
		
      # Design matrix ----------
      # Remove factors not used.
			sample_info_filter <- sample_info[, key_model_factors, drop = F]
			groups <- apply(sample_info_filter, 1, function(x) paste(x, collapse = "_"))
			unique_groups <- unique(groups)  

      groups <- factor(groups, levels = unique_groups)

		  design <- stats::model.matrix(~ 0 + groups)  
		  colnames(design) <- gsub("^groups", "", colnames(design))
            	
			if(!is.null(raw_counts) && counts_deg_method == 2) {
				voom_results <- limma::voom(raw_counts, design) 
				fit <- limma::lmFit(voom_results, design) 
			} else {
        fit <- limma::lmFit(eset, design)
      }
		
			fit <- limma::eBayes(fit, trend = limma_trend)	
			
			# Making comaprisons-----------
      # Only one factor, or more than two then use all pairwise comparisons
			if(length(key_model_factors) != 2 | length(block_factor) > 1)  {
				comparisons <- gsub(".*: ", "", selected_comparisons)
				comparisons <- gsub(" vs\\. ", "-", comparisons)
      # Two key factors
			} else if(length(key_model_factors) == 2){ 
				if(sum(grepl(":", model_factors) > 0)) { 
					interaction_term <- TRUE
					# All pairwise comparisons
					comparisons <- ""
					for(i in 1:(length(unique_groups) - 1)) {
            for(j in (i + 1):length(unique_groups)) {
              # Only compare WT_control vs. WT_treatment
              if(strsplit(unique_groups[i], "_")[[1]][1] == strsplit(unique_groups[j], "_")[[1]][1] |
                 strsplit(unique_groups[i], "_")[[1]][2] == strsplit(unique_groups[j], "_")[[1]][2]) {
                comparisons <- c(comparisons, paste(unique_groups[j], "-", unique_groups[i], sep = ""))
              }
            }
          }
					comparisons <- comparisons[-1]
					
					# Pairwise contrasts
				  make_contrast <- limma::makeContrasts(
            contrasts = comparisons[1],
            levels = design
          )
					if(length(comparisons) > 1) {
            for(kk in 2:length(comparisons)) {
              make_contrast <- cbind(
                make_contrast,
                limma::makeContrasts(contrasts = comparisons[kk], levels = design)
              )
            }
          }
					
						
					contrast_names <- colnames(make_contrast)		
				
					# All possible interactions
					# Interaction contrasts
					
					contrast_compare <- NULL
					contrast_names <- ""
					for (kk in 1:(dim(make_contrast)[2] - 1)) {
					  for(kp in (kk+1):dim(make_contrast)[2]) {
              if(is.null(contrast_compare)) {
                contrast_compare <- make_contrast[, kp] - make_contrast[,kk]
              } else {
                contrast_compare <- cbind(
                  contrast_compare,
                  make_contrast[, kp] - make_contrast[, kk]
                )
              }
							contrast_names <- c(
                contrast_names,
                paste0(
                  "I:",
                  colnames(make_contrast)[kp],
                  ".vs.",
                  colnames(make_contrast)[kk]
                )
              )
						}   
					}
					colnames(contrast_compare) <- contrast_names[-1]
					
					# Remove nonsense contrasts from interactions
					contrast_compare <- contrast_compare[, which(
            apply(abs(contrast_compare), 2, max) == 1), drop = F
          ]
					contrast_compare <- contrast_compare[, which(
            apply(abs(contrast_compare), 2, sum) == 4), drop = F
          ]
          # Remove duplicate columns
					contrast_compare <- t(unique(t(contrast_compare)))		
					
					# Remove unwanted contrasts involving more than three levels in 
          # either factor
					keep <- c()
					for(i in 1:dim(contrast_compare)[2]) {
						tem <- rownames(contrast_compare)[contrast_compare[, i] != 0]
						tem1 <- unique(unlist(gsub("_.*", "", tem)))
						tem2 <- unique(unlist(gsub(".*_", "", tem)))
						if(length(tem1) == 2 & length(tem2) == 2) {
              keep <- c(keep, colnames(contrast_compare)[i])
            }
					}
					contrast_comapre <- contrast_compare[, keep, drop = F]
					comparison_names = colnames(contrast_compare) 				 
				}

				# "stage: MN vs. EN"  -->  c("MN_AB-EN_AB", "EN_Nodule-EN_AB") 
				#  comparisons in all levels of the other factor 
				transform_comparisons <- function(
          comparison,
          key_model_factors
        ) {
					tem <- gsub(".*: ", "", comparison)
          # control  mutant
					tem <- unlist(strsplit(tem, " vs\\. ")) 							
					factor <- gsub(":.*", "", comparison)
          
          # 1: first factor, 2: 2nd factor
					ix <- match(factor, key_model_factors)
          # 3-1 = 2; 3-1=1
					other_factor <- key_model_factors[3 - ix]
					other_factor_levels <- unique(sample_info_filter[, other_factor])				
					comparisons <- c()
					
					for(factor_levels in other_factor_levels) {
						if(ix == 1) {
							comparisons <- c(
                comparisons,
                paste(paste0(tem, "_", factor_levels), collapse = "-")
              )
						} else {
							comparisons <- c(
                comparisons,
                paste(paste0(factor_levels, "_", tem), collapse = "-")
              )
						}
					}
					return(comparisons)		
				}	
				
				comparisons <- unlist(sapply(selected_comparisons, transform_comparisons))
				comparisons <- as.vector(comparisons)
			}
			
			# make contrasts
			make_contrast <- limma::makeContrasts(contrasts = comparisons[1], levels = design)
			if(length(comparisons) > 1) {
        for(kk in 2:length(comparisons)) {
          make_contrast <- cbind(
            make_contrast,
            limma::makeContrasts(contrasts = comparisons[kk], levels = design)
          )
        }
      }
				
			if(interaction_term) {
				make_contrast <- cbind(make_contrast, contrast_compare)
				contrast_names <- c(colnames(make_contrast), colnames(contrast_compare))
				comparisons <- c(comparisons, comparison_names)
			}
			
      # Factor is selected as block
			if(length(block_factor) >= 1) { 
				if(length(block_factor) >= 1) {
          # If multiple use the first one
          block_factor <- block_factor[1] 
        }

				block <- sample_info[, block_factor]
			
				if(!is.null(raw_counts) && counts_deg_method == 2) {
					voom_results <- limma::voom(raw_counts, design)
					corfit <- limma::duplicateCorrelation(voom_results, design, block = block)			
					fit <- limma::lmFit(
            voom_results,
            design,
            block = block,
            correlation = corfit$consensus
          ) 
				} else {
					corfit <- limma::duplicateCorrelation(eset, design, block = block)
					fit <- limma::lmFit(
            eset,
            design,
            block = block,
            correlation = corfit$consensus
          )
				}
				fit <- limma::eBayes(fit, trend = limma_trend)	
			}
		}
		
		fit_contrast <- limma::contrasts.fit(fit, make_contrast)
		fit_contrast <- limma::eBayes(fit_contrast, trend = limma_trend)
		results <- limma::decideTests(
      fit_contrast,
      p.value = max_p_limma,
      lfc = log2(min_fc_limma )
    )
		# Extract fold change for each comparison
		# There is issues with direction of foldchange. Sometimes opposite
		top <- function(
      comp,
      fit_contrast,
      processed_data
    ) {
			tem <- limma::topTable(
        fit_contrast,
        number = 1e12,
        coef = comp,
        sort.by = "M"
      ) 
			if(dim(tem)[1] == 0) {
        return(1) 
			} else { 			
				# Compute fold change for the first gene (ranked by absolute value)
				tem2 <- as.numeric(
          processed_data[which(rownames(processed_data) == rownames(tem)[1]), ]
        )
				names(tem2) <- colnames(processed_data) 
					
				return(tem[, c(1,5)]) 
			}											
		}
		
		
		top_genes <- lapply(comparisons, function(x) {
      top(
        comp = x,
        fit_contrast = fit_contrast,
        processed_data = processed_data
      )
    })
		top_genes <- setNames(top_genes, comparisons)

		ix <- which(unlist(lapply(top_genes, class)) == "numeric")
		if(length(ix) > 0) {
      top_genes <- top_genes[-ix]
    } 
	}
    
	return(list(
    results = results,
    comparisons = comparisons,
    exp_type = exp_type,
    top_genes = top_genes
  )) 
}

#' Significant genes bar plot
#' 
#' Create a bar plot of the number of genes with a significant fold
#' change. The plot will break down all of the comparisons that
#' were analyzed and the count of significant genes for both up
#' and down changes.
#' 
#' @param results Results matrix from the limma_value function
#'  returned list
#' 
#' @export
#' @return Formatted gg barplot of the significantly expressed
#'  genes.
sig_genes_plot <- function(
  results
) {
	Up <-  apply(results, 2, function(x) sum(x == 1))
	Down <- apply(results, 2, function(x) sum(x == -1)) 
	stats <- rbind(Up, Down)
				 
	gg <- reshape2::melt(stats)

	colnames(gg) <- c("Regulation","Comparisons","Genes")
		 
	plot_bar <- ggplot2::ggplot(
    gg,
    ggplot2::aes(x = Comparisons, y = Genes, fill = Regulation)
  ) +
  ggplot2::geom_bar(position = "dodge", stat = "identity") +
  ggplot2::coord_flip() +
  ggplot2::theme(
    legend.position = "top",
    axis.title.y = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 14),
    legend.text = ggplot2::element_text(size = 12)
  ) +
  ggplot2::theme_light() +
  ggplot2::ylab("Number of differntially expressed genes") +
	ggplot2::geom_text(
    ggplot2::aes(label = Genes),
    position = ggplot2::position_dodge(width = 0.9),
    vjust = 0.5,
    hjust = 0
  )

	return(plot_bar)
}

#' Create a table from DEG results
#' 
#' Using the limma_value return list, create a table of the
#' number of significantly expressed genes for each analyzed
#' comparison.
#' 
#' @param limma Return list from the limma_value function
genes_stat_table <- function(
  limma
) {
  results <- limma$results
  
  # If only one comparison
  if(dim(results)[2] == 1) { 
    Up <- sum(results == 1)
    Down <- sum(results == -1)
    stats <- c(colnames(results), Up, Down)
    stats <- t(as.data.frame(stats))
    row.names(stats) <- colnames(results)
    colnames(stats) <- c("Comparison","Up", "Down")
  
  # More than one comparisons
  } else {  
		Up <-  apply(results, 2, function(x) sum(x == 1))
		Down <- apply(results, 2, function(x) sum(x == -1)) 
		stats <- rbind(Up, Down)
		stats <- t(stats)
		stats <- cbind(rownames(stats), stats)
		colnames(stats)[1] <- "Comparisons"
    # Reverse row order, to be the same with plot
		stats <- stats[dim(stats)[1]:1, ] 
  }
		 
  return(as.data.frame(stats))
}

#' List comparisons for venn diagram plot
#' 
#' Create a list of the comparisons that were detected and
#' analyzed by the limma_value function. These comparisons
#' can be used in the plot_venn function to find the number
#' of overlapping significantly expressed genes.
#' 
#' @param limma Returned list of results from the limma_value
#'  function
#' @param up_down_regulated Split the comparisons into either
#'  up or down regulated
#' 
#' @export
#' @return A character vector of the comparisons that were used
#'  in the DEG analysis and can be plotted with the venn_plot
#'  function.
list_comp_venn <- function(
  limma,
  up_down_regulated
) {
  if(is.null(limma$comparisons)) {
    return(NULL)  
	}	else {
    choices <- stats::setNames(limma$comparisons, limma$comparisons)
		
    if(up_down_regulated) {
      tem <- c(
        paste0("Up_", limma$comparisons),
        paste0("Down_", limma$comparisons)
      )
			choices <- stats::setNames(tem, tem)
		}
				
		choices_first_three <- choices
		
    if(length(choices_first_three) > 3) {
      # By default only 3 are selected
      choices_first_three <- choices[1:3]
    }
		return(list(
      choices = choices,
      choices_first_three = choices_first_three
    ))	
	} 
}

#' Create a venn diagram plot
#' 
#' Plot a venn diagram that illustrates the number of significantly
#' expressed genes that overlap for multiple comparisons. 
#' 
#' @param limma limma Returned list of results from the limma_value
#'  function
#' @param up_down_regulated Split the comparisons into either
#'  up or down regulated
#' @param select_comparisons_venn The comparisons to plot on the
#'  venn diagram
#' 
#' @export
#' @return A formatted venn diagram plot of the selected comparisons.
plot_venn <- function(
  limma,
  up_down_regulated,
  select_comparisons_venn
) {
  results <- limma$results

	# Split by up or down regulation
	if(up_down_regulated) {
    result_up <- results 
		result_up[result_up < 0] <- 0
		colnames(result_up) <- paste0("Up_", colnames(result_up))
		result_down <- results 
		result_down[result_down > 0] <- 0
		colnames(result_down) <- paste0("Down_", colnames(result_down))				
		results <- cbind(result_up, result_down)
	}			
	
  ixa <- c()
	for(comps in select_comparisons_venn) {
    # If not interaction term
		if(!grepl("^I:|^I-|^Up_I:|^Up_I-|^Down_I:|^Down_I-", comps) ) {  
			ix <- match(comps, colnames(results)) 
		} else {
			# Mismatch in comparison names for interaction terms for DESeq2
			# I:water_Wet.genetic_Hy 	in the selected Contrast
			# Diff-water_Wet-genetic_Hy  in column names
			tem <- gsub("^I-", "I:", colnames(results))
			tem <- gsub("-", "\\.", tem)
			ix <- match(comps, tem) 
      
      # This is for limma package
			if(is.na(ix)) {
        ix <- match(comps, colnames(results))
      } 						
		}
		ixa <- c(ixa,ix)
	}
  # Only use selected comparisons
	results <- results[, ixa, drop = FALSE]
	if(dim(results)[2] > 5) {
    results <- results[, 1:5]
  }
  colnames(results) <- gsub("^I-", "I:", colnames(results))	
	
  return(
    limma::vennDiagram(
      results,
      circle.col = rainbow(5),
      cex = c(1., 1, 0.7)
    )
  )
}

#' Find data for the heatmap
#' 
#' Filter the processed data into a submatrix that only contains
#' the genes that had a significant fold change. The samples will
#' also be subsetted to only contain samples that pertain to the
#' selected comparison or contrast. 
#' 
#' @param limma Return results list from the limma_value function
#' @param select_contrast Comparison from DEG analysis to filter
#'  for the significant genes
#' @param processed_data Data that has been through the pre-processing
#' @param contrast_samples Columns that are in the group of the
#'  selected comparison
#' 
#' @export
#' @return Submatrix of the processed data with only the significantly
#'  expressed genes and the columns that are in the selected contrast
#'  group.
deg_heat_data <- function(
  limma,
  select_contrast,
  processed_data,
  contrast_samples
) {
	genes <- limma$results

  if(is.null(genes)) {
    return(NULL)
  }
	
  # If not interaction term
  if(!grepl("I:", select_contrast)) {
		ix <- match(select_contrast, colnames(genes)) 
	} else {
		# Mismatch in comparison names for interaction terms for DESeq2
		# I:water_Wet.genetic_Hy in the selected Contrast
		# Diff-water_Wet-genetic_Hy in column names
	  tem <- gsub("I-", "I:", colnames(genes))
		tem <- gsub("-", "\\.", tem)
		ix <- match(select_contrast, tem) 
		
    # This is for limma package
		if(is.na(ix)) {
      ix <- match(select_contrast, colnames(genes)) 	
    }		
			
	}
	  
	if(is.null(ix) || is.na(ix)) {
    return(NULL)
  }
  # No significant genes for this comparison
	if(sum(abs(genes[, ix])) <= 1) {
    return(NULL) 
  }
  if(dim(genes)[2] < ix) {
    return(NULL)
  }	
  query <- rownames(genes)[which(genes[, ix] != 0)]
	if(length(query) == 0) {
    return(NULL)
  }
  iy <- match(query, rownames(processed_data))
		  

	iz <- contrast_samples
	
	# Color bar
	bar <- as.vector(genes[, ix])
	names(bar) <- row.names(genes)
	bar <- bar[bar != 0]

	# Retreive related data		 
	genes <- processed_data[iy, iz, drop = FALSE]

	genes <- genes[order(bar), , drop = FALSE]
	bar <- sort(bar)

	return(list(
    genes = genes,
    bar = bar
  ))
}

#' Heatmap for significant contrast genes
#' 
#' Create a ComplexHeatmap of the processed expression data for
#' the genes that were significantly expressed in the selected 
#' comparison. The data for this heatmap comes from the
#' deg_heat_data function.
#' 
#' @param data Submatrix of the processed data matrix from the
#'  deg_heta_data function
#' @param bar Vector to signify a positive (1) expression fold
#'  change or a negative (-1) change
#' @param heatmap_color_select Color vector to use for the
#'  heatmap expression scale
#' 
#' @export
#' @return A drawn heatmap from the filtered data.
deg_heatmap <- function(
  data,
  bar,
  heatmap_color_select
) {
  # Number of genes to show
	n_genes <- as.character(table(bar))

  data <- as.matrix(data) - apply(data, 1, mean)
  cutoff <- median(unlist(data)) + 3 * sd(unlist(data)) 
	data[data > cutoff] <- cutoff
	cutoff <- median(unlist(data)) - 3 * sd(unlist(data)) 
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
  
	groups <- detect_groups(colnames(data))
	group_count <- length(unique(groups))
  groups_colors <- gg_color_hue(2 + group_count)
  
  top_ann <- ComplexHeatmap::HeatmapAnnotation(
    Group = groups,
    col = list(
      Group = setNames(groups_colors[1:group_count], unique(groups))
    ),
    annotation_legend_param = list(
      Group = list(nrow = 1, title = NULL)
    ),
    show_annotation_name = list(Group = FALSE),
    show_legend = FALSE
  )
  
  bar[bar == -1] <- "Negative"
  bar[bar == 1]  <- "Positive"
  groups <- bar

  row_ann <- ComplexHeatmap::rowAnnotation(
    Change = groups,
    col = list(
      Change = setNames(
        groups_colors[(group_count + 1):length(groups_colors)], unique(groups)
      )
    ),
    annotation_legend_param = list(
      Change = list(nrow = 1, title = NULL)
    ),
    show_annotation_name = list(Change = FALSE),
    show_legend = FALSE
  )

  heat <- ComplexHeatmap::Heatmap(
    data,
    name = "Expression",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    left_annotation = row_ann,
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

#' Plot brush selection from main heatmap
#' 
#' Create a ComplexHeatmap object from the User brush selection
#' that is a sub plot of the main plot.
#' 
#' @param ht_brush Input from the user creating a brush selection
#'  on the main heatmap
#' @param ht Main heatmap from the deg_heatmap function
#' @param ht_pos_main Main heatmap position information to use
#'  for the sub heatmap
#' @param heatmap_data Original data matrix that was plotted in
#'  the main heatmap
#' @param all_gene_names Data matrix of all the mapped gene names
#' 
#' @export
#' @return A ComplexHeatmap object of the brushed selection from
#'  the main heatmap.
deg_heat_sub <- function(
  ht_brush,
  ht,
  ht_pos_main,
  heatmap_data,
  all_gene_names
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
    column_groups <- detect_groups(colnames(heatmap_data$genes))
	  group_count <- length(unique(column_groups))
    groups_colors <- gg_color_hue(2 + group_count)
  
    top_ann <- ComplexHeatmap::HeatmapAnnotation(
      Group = column_groups,
      col = list(
        Group = setNames(
          groups_colors[1:group_count],
          unique(column_groups)
        )
      ),
      annotation_legend_param = list(
        Group = list(nrow = 1, title = NULL)
      ),
      show_annotation_name = list(Group = FALSE),
      show_legend = TRUE
    )
    
    bar <- heatmap_data$bar
    bar[bar == -1] <- "Down"
    bar[bar == 1]  <- "Up"
    row_groups <- bar

    row_ann <- ComplexHeatmap::rowAnnotation(
      Change = row_groups,
      col = list(
        Change = setNames(
          groups_colors[(group_count + 1):length(groups_colors)],
          unique(row_groups)
        )
      ),
      annotation_legend_param = list(
        Change = list(nrow = 1, title = NULL)
      ),
      show_annotation_name = list(Change = FALSE),
      show_legend = TRUE
    )
    
    group_col_return <- setNames(
      groups_colors,
      c(unique(column_groups), unique(row_groups))
    )
  # End annotation ---------

  column_index <- unlist(pos[1, "column_index"])
  row_index <- unlist(pos[1, "row_index"])
  top_ann <- top_ann[column_index]
  row_ann <- row_ann[row_index]
  column_groups <- column_groups[column_index]
  m <- ht@ht_list[[1]]@matrix

  bar_return <- bar[row_index]

  if (length(row_index) > 50) {
    show_rows <- FALSE
  } else {
    show_rows <- TRUE
  }
  if(ncol(all_gene_names) == 3) {
    genes <- rowname_id_swap(
      data_matrix = m[row_index, column_index, drop = FALSE],
      all_gene_names = all_gene_names,
      select_gene_id = "symbol"
    )
  } else {
    genes <- m[row_index, column_index, drop = FALSE]
  }
  submap_data <- genes

  ht_select <- ComplexHeatmap::Heatmap(
    genes,
    col = ht@ht_list[[1]]@matrix_color_mapping@col_fun,
    show_heatmap_legend = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = show_rows,
    top_annotation = top_ann,
    left_annotation = row_ann,
    name = "heat_1"
  )

  return(list(
    ht_select = ht_select,
    submap_data = submap_data,
    group_colors = group_col_return,
    column_groups = column_groups,
    bar = bar_return
  ))
}

#' HTML code for sub-heatmap selected cell
#' 
#' Create HTML code for a cell of information on the cell of the
#' sub-heatmap that the User clicks on. The cell contains the
#' expression value, the sample, the gene, the group and the
#' direction of the fold change. 
#' 
#' @param click Information fro what cell is clicked in the
#'  sub-heatmap
#' @param ht_sub The drawn sub-heatmap
#' @param ht_sub_obj The sub-heatmap ComplexHeatmap object
#' @param ht_pos_sub Position information for the sub-heatmap
#' @param sub_groups Vector of the groups that the samples
#'  belong to
#' @param group_colors The color of the top annotation that
#'  is used for each group and the side annotation that denotes
#'  the direction of the expression regulation
#' @param bar Vector to signify a positive (1) expression fold
#'  change or a negative (-1) change
#' @param data Sub data matrix that is plotted in the sub-heatmap
#' 
#' @export
#' @return HTML code that will be used in the shiny UI to tell
#'  the user the information of the cell they selected.
deg_click_info <- function(
  click,
  ht_sub,
  ht_sub_obj,
  ht_pos_sub,
  sub_groups,
  group_colors,
  bar,
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
  up_down <- bar[row_index]
  up_down_col <- group_colors[[up_down]]
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
Regulation: @{up_down} <span style='background-color:@{up_down_col};width=50px;'>    </span>   
</pre></div>"
)

 return(HTML(html))
}

#' Volcano DEG plot
#' 
#' Use the results from limma-value to create a volcano plot to
#' illustrate the significantly expressed genes. Up and down
#' regulated genes are colored on the ggplot.
#' 
#' @param select_contrast Comparison from DEG analysis to filter
#'  for the significant genes
#' @param comparisons The comparisons vector from the results list
#'  of the limma_value function
#' @param top_genes top_genes list from results list of the
#'  limma_value function
#' @param limma_p_val Significant p-value to use to in determining
#'  the expressed genes
#' @param limma_fc Minimum fold change value to use in determining
#'  the expressed genes
#' @param plot_colors List containing three colors to differentiate between   
#'  the up-regulated, down-regulated, and other genes
#' 
#' @export
#' @return ggplot with the fold value as the X-axis and the log 10
#'  value of the adjusted p-value as the Y-axis.
plot_volcano <- function(
  select_contrast,
  comparisons,
  top_genes,
  limma_p_val,
  limma_fc,
  plot_colors
) {
  if(is.null(select_contrast) || is.null(comparisons) ||
     length(top_genes) == 0) {
    return(NULL)
  }
	if(length(comparisons)  == 1) { 
    top_1 = top_genes[[1]]  
	} else {
	  top <- top_genes
	  ix <- match(select_contrast, names(top))
	  if(is.na(ix)) {
      return (NULL)
    }
	  top_1 <- top[[ix]]
	}
  if(dim(top_1)[1] == 0 ) {
    return(NULL)
  } 
  colnames(top_1) <- c("Fold","FDR")
  # Convert to data frame
	top_1 <- as.data.frame(top_1)
  # Remove NA's
  top_1 <- top_1[which(!(is.na(top_1$Fold) | is.na(top_1$FDR))), ] 
	top_1$upOrDown <- "None"
	top_1$upOrDown[which(
    top_1$FDR <= limma_p_val & top_1$Fold  >= log2(limma_fc)
  )] <- "Up"
	top_1$upOrDown[which(
    top_1$FDR <= limma_p_val & top_1$Fold  <= -log2(limma_fc)
  )] <- "Down"

  return(
    ggplot2::ggplot(top_1, ggplot2::aes(x = Fold, y = -log10(FDR))) +
      ggplot2::geom_point(ggplot2::aes(color = upOrDown))	+
      ggplot2::scale_color_manual(values = plot_colors) + 
      ggplot2::theme_light() +
      ggplot2::theme(
        legend.position = "right",
        axis.title.x = ggplot2::element_text(
          color = "black",
          size = 14
        ),
        axis.title.y = ggplot2::element_text(
          color = "black",
          size = 14
        ),
        axis.text.x = ggplot2::element_text(
          size = 16
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
        title = "Fold Change vs. Adjusted p-Value",
        y = "-log10(Adjusted p-Val)",
        x = "Fold Change",
        color = "Regulated"
      )
  ) 
}

#' Plot mean expression and fold change
#' 
#' Draw a ggplot of the overal mean expression for each gene 
#' and the calculated fold change for the selected comparison.
#' 
#' @param select_contrast Comparison from DEG analysis to filter
#'  for the significant genes
#' @param comparisons The comparisons vector from the results list
#'  of the limma_value function
#' @param top_genes top_genes list from results list of the
#'  limma_value function
#' @param limma_p_val Significant p-value to use to in determining
#'  the expressed genes
#' @param limma_fc Minimum fold change value to use in determining
#'  the expressed genes
#' @param contrast_samples Samples that are included in the selected
#'  comparison
#' @param processed_data Data matrix that has gone through
#'  pre-processing
#' @param plot_colors List containing three colors to differentiate between   
#'  the up-regulated, down-regulated, and other genes
#' 
#' @export
#' @return A ggplot with the X-axis the mean expression value and
#'  the Y-axis the calculated fold-change from the DEG analysis.
plot_ma <- function(
  select_contrast,
  comparisons,
  top_genes,
  limma_p_val,
  limma_fc,
  contrast_samples,
  processed_data, 
  plot_colors
) {
  if(grepl("I:", select_contrast)) {
    return(NULL)
  }
  if(length(comparisons) == 1) {
    top_1 <- top_genes[[1]]  
	} else {
	  top <- top_genes
	  ix <- match(select_contrast, names(top))
	  if(is.na(ix)) {
      return(NULL)
    }
	  top_1 <- top[[ix]] 
	}
  if(dim(top_1)[1] == 0) {
    return(NULL)
  }
  colnames(top_1) <- c("Fold", "FDR")
  # Convert to data frame
  top_1 <- as.data.frame(top_1)
  # Remove NA's
  top_1 <- top_1[which(!(is.na(top_1$Fold) | is.na(top_1$FDR))), ]
  top_1$upOrDown <- "None"
	top_1$upOrDown[which(
    top_1$FDR <= limma_p_val & top_1$Fold  >= log2(limma_fc)
  )] <- "Up"
	top_1$upOrDown[which(
    top_1$FDR <= limma_p_val & top_1$Fold  <= -log2(limma_fc)
  )] <- "Down"

  iz <- contrast_samples

  average_data <- as.data.frame(apply(processed_data[, iz], 1, mean))
  colnames(average_data) <- "Average"
	rownames(average_data) <- rownames(processed_data)
		
	genes <-  merge(average_data, top_1, by = "row.names")


  return(
    ggplot2::ggplot(genes, ggplot2::aes(x = Average, y = Fold)) +
      ggplot2::geom_point(ggplot2::aes(color = upOrDown))	+
      ggplot2::scale_color_manual(values = plot_colors) +
      ggplot2::theme_light() +
      ggplot2::theme(
        legend.position = "right",
        axis.title.x = ggplot2::element_text(
          color = "black",
          size = 14
        ),
        axis.title.y = ggplot2::element_text(
          color = "black",
          size = 14
        ),
        axis.text.x = ggplot2::element_text(
          size = 16
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
        title = "Average Expression vs. Log2 Fold Change",
        y = "Log2 Fold Change",
        x = "Average Expression",
        color = "Regulated"
      )
  ) 
}

#' Scatter plot of the comparison groups
#' 
#' Create a scatter plot of the expression value for each gene 
#' in the two groups for the selected contrast. For the selected
#' contrast, the mean expression is calculated for a gene in both
#' group of samples in the contrast and plotted in a scatter plot.
#' 
#' @param select_contrast Comparison from DEG analysis to filter
#'  for the significant genes
#' @param comparisons The comparisons vector from the results list
#'  of the limma_value function
#' @param top_genes top_genes list from results list of the
#'  limma_value function
#' @param limma_p_val Significant p-value to use to in determining
#'  the expressed genes
#' @param limma_fc Minimum fold change value to use in determining
#'  the expressed genes
#' @param contrast_samples Samples that are included in the selected
#'  comparison
#' @param processed_data Data matrix that has gone through
#'  pre-processing
#' @param sample_info Experiment file information for grouping
#' 
#' @export
#' @return A formatted ggplot with the X-axis as the mean expression
#'  of one contrast group and the Y-axis as the mean expression of
#'  the other contrast group.
plot_deg_scatter <- function(
  select_contrast,
  comparisons,
  top_genes,
  limma_p_val,
  limma_fc,
  contrast_samples,
  processed_data,
  sample_info
) {
	if(grepl("I:", select_contrast)) {
    return(NULL)
  }

	if(length(comparisons)  == 1) {
    top_1 <- top_genes[[1]]  
	} else {
	  top <- top_genes
	  ix <- match(select_contrast, names(top))
	  if(is.na(ix)) {
      return(NULL)
    }
	  top_1 <- top[[ix]] 
	}
	  
  if(dim(top_1)[1] == 0) {
    return(NULL)
  }
	colnames(top_1) <- c("Fold", "FDR")
  # Convert to data frame
	top_1 <- as.data.frame(top_1)
  # Remove NA's
  top_1 <- top_1[which(!(is.na(top_1$Fold) | is.na(top_1$FDR))), ] 
	top_1$upOrDown <- "None"
	top_1$upOrDown[which(
    top_1$FDR <= limma_p_val & top_1$Fold  >= log2(limma_fc)
  )] <- "Up"
	top_1$upOrDown[which(
    top_1$FDR <= limma_p_val & top_1$Fold  <= -log2(limma_fc)
  )] <- "Down"

	iz <- contrast_samples
	
  genes <- processed_data[,iz]
	 
	g <- detect_groups(colnames(genes), sample_info)
	 
	if(length(unique(g)) > 2) {
    plot.new()
    text(0.5,0.5, "Not available for more than two groups.")
  } else{
		average_1 <- apply(genes[, which(g == unique(g)[1])], 1, mean)

		average_2 <- apply(genes[, which(g == unique(g)[2])], 1, mean)

		genes_1 <- cbind(average_1, average_2)
		rownames(genes_1) <- rownames(genes)
		genes_1 <-  merge(genes_1, top_1, by = "row.names")

    colors <- gg_color_hue(2)

    return(
      ggplot2::ggplot(genes_1, ggplot2::aes(x = average_1, y = average_2)) +
        ggplot2::geom_point(ggplot2::aes(color = upOrDown))	+
        ggplot2::scale_color_manual(values = c(colors[1], "grey45", colors[2])) +
        ggplot2::theme_light() +
        ggplot2::theme(
          legend.position = "right",
          axis.title.x = ggplot2::element_text(
            color = "black",
            size = 14
          ),
          axis.title.y = ggplot2::element_text(
            color = "black",
            size = 14
          ),
          axis.text.x = ggplot2::element_text(
            size = 16
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
          title = "Average Expression in Group",
          y = paste0("Average Expression: ", unique(g)[2]),
          x = paste0("Average Expression: ", unique(g)[1]),
          color = "Regulated"
        )
    )
	} 
}

#' Dendogram of enriched pathways
#' 
#' Create a dendogram plot of the enriched pathways to illustrate
#' which paths contain similar genes.
#' 
#' @param go_table Enrichment table from the pathway analysis
#'  functions
#' @param right_margin Control the size of the dendogram labels
#' 
#' @export
#' @return A dendogram plot that shows the users what pathways are
#'  that are enriched share genes.
enrichment_plot <- function(
  go_table,
  right_margin = 10
) {
  # a program for ploting enrichment results by highlighting the similarities among terms
  # must have columns: Direction, adj.Pval   Pathways Genes
  #  Direction	adj.Pval	nGenes	Pathways		Genes
  #Down regulated	3.58E-59	131	Ribonucleoprotein complex biogenesis	36	Nsun5 Nhp2 Rrp15 
  #Down regulated	2.55E-57	135	NcRNA metabolic process	23	Nsun5 Nhp2 Rrp15 Emg1 Ddx56 Rsl1d1
  # Up or down regulation is color-coded
  # gene set size if represented by the size of marker
  data <- go_table
  if(class(data) != "data.frame") {
    return(NULL)
  }
  # only one term or less
  if(nrow(data) <=1 || is.null(data)) {
    return(NULL)
  }
  
  gene_lists <- lapply(
    data$Genes,
    function(x) unlist(strsplit(as.character(x), " "))
  )
  names(gene_lists) <- data$Pathways

  # Compute overlaps percentage--------------------
  n <- length(gene_lists)
  w <- matrix(NA, nrow = n, ncol = n)
  
  # Compute overlaps among all gene lists
  for(i in 1:n) {
    for (j in i:n) {
      u <- unlist(gene_lists[i])
      v <- unlist(gene_lists[j])
      w[i, j] <- length(intersect(u, v)) / length(unique(c(u, v)))
    }
  }
  # The lower half of the matrix filled in based on symmetry
  for(i in 1:n) {
    for(j in 1:(i-1)) {
      w[i, j] <- w[j, i]
    } 
  }

  Terms <- paste(
    sprintf("%-1.0e",
    as.numeric(data$adj_p_val)), 
	  names(gene_lists)
  )
  rownames(w) <- Terms
  colnames(w) <- Terms

  # A large margin for showing 
  par(mar = c(0, 0, 1, right_margin))

  dend <- stats::as.dist(1 - w) |>
	  stats::hclust(method = "average")
  # Permutated order of leaves 
  ix <- dend$order 

  leaf_type <- as.factor(data$Direction[ix])
  leaf_colors <- gg_color_hue(length(unique(data$Direction)))
  # Leaf size represent P values
  leaf_size <- -log10(as.numeric(data$adj_p_val[ix]))
  leaf_size <- 1.5 * leaf_size / max(leaf_size) + .2
  
	dend |> 
	  stats::as.dendrogram(hang = -1) |>
    # Type of marker
	  dendextend::set("leaves_pch", 19) |>
    # Size
	  dendextend::set("leaves_cex", leaf_size) |>
    # up or down genes
	  dendextend::set("leaves_col", leaf_colors[leaf_type]) |>
	  plot(horiz = TRUE)
	
  # Add legend using a second layer
  par(lend = 1)
	add_legend(
    "top",
    pch = 19,
    col = leaf_colors,
    legend = levels(leaf_type),
    bty = "n",
    horiz = T 
  )
}

#' Create a single table from up and down enrichments
#' 
#' Use the enrichment table from the down genes analysis and
#' the up genes analysis to creeate a single table.
#' 
#' @param up_enrich_data Enrichment table from the up-regulated
#'  genes pathway analysis
#' @param down_enrich_data Enrichment table from the down-regulated
#'  genes pathway analysis.
#' 
#' @export
#' @return A combined enrichment analysis table.
go_table_data <- function(
 up_enrich_data,
 down_enrich_data 
) {
  if(nrow(up_enrich_data) >= 2 || is.null(up_enrich_data)) {
    up_data <- as.data.frame(up_enrich_data)
    up_data$direction <- rep("Up", nrow(up_enrich_data))
    up_data <- up_data[, c(6, 1, 2, 4, 5)]
    colnames(up_data) <- c("Direction", "adj_p_val", "n_genes", "Pathways", "Genes")
  } else {
    up_data <- NULL
  }
  if(nrow(down_enrich_data) >= 2 || is.null(down_enrich_data)) {
    down_data <- as.data.frame(down_enrich_data)
    down_data$direction <- rep("Down", nrow(down_enrich_data))
    down_data <- down_data[, c(6, 1, 2, 4, 5)]
    colnames(down_data) <- c("Direction", "adj_p_val", "n_genes", "Pathways", "Genes")
  } else {
    down_data <- NULL
  }
  data <- rbind(up_data, down_data)

  return(data)
}

#' VisNetwork data
#' 
#' Create VisNetwork data that can be inputted in the vis_network_plot
#' function to create an interactive network of enriched pathways.
#' 
#' @param network GO table from the pathway analysis
#' @param up_down_reg_deg Plot just up/down or both
#' @param wrap_text_network_deg Wrap the text from the pathway description
#' @param layout_vis_deg BUtton to reset the layout of the network
#' @param edge_cutoff_deg P-value to cutoff enriched pathways
#' 
#' @export
#' @return Data that can be inputted in the vis_network_plot function
#'  to create an interactive network.
network_data <- function(
  network,
  up_down_reg_deg,
  wrap_text_network_deg,
  layout_vis_deg,
  edge_cutoff_deg
) {
  if(up_down_reg_deg != "Both") {
    network <- network[network$Direction == up_down_reg_deg, ]
  }
  if(dim(network)[1] == 0) {
    return(NULL)
  }
  
  if(wrap_text_network_deg) {
    # Wrap long pathway names using default width of 30
    network$Pathways <- wrap_strings(network$Pathways)
  }

  g <- enrichment_network(
    network,
    layout_button = layout_vis_deg,
    edge_cutoff = edge_cutoff_deg
  )

  vis_net <- visNetwork::toVisNetworkData(g)
    
  # Color codes: https://www.rapidtables.com/web/color/RGB_Color.html
  vis_net$nodes$shape <- "dot"
    
  vis_net$nodes$size <- 5 + vis_net$nodes$size^2 
    
  return(vis_net)
}

deg_information <- function(
    limma_value, 
    gene_names, 
    processed_data, 
    no_id_conversion = FALSE 
){
  if (no_id_conversion){
    
    # get the first comparison level 
    degs_data <- limma_value$top_genes[[1]]
    degs_data$User_ID <- rownames(degs_data)
    
    # get the additional comparison levels if they exists 
    if (length(names(limma_value$top_genes)) > 1){
      for (i in 2:length(names(limma_value$top_genes))){
        temp <- limma_value$top_genes[[i]]
        temp$User_ID <- rownames(temp)
        degs_data <- dplyr::inner_join(degs_data, temp, by = "User_ID")
      }
    }
    
    # connect to gene symbols and original user id 
    processed_data <- as.data.frame(processed_data)
    processed_data$User_ID <- rownames(processed_data)
    
    degs_data <- dplyr::full_join(degs_data, gene_names, by = "User_ID")
    degs_data <- dplyr::full_join(degs_data, processed_data, by = "User_ID")
    
    
    degs_data <- degs_data %>%
      dplyr::relocate(User_ID)
    
  } else {
    
    # get the first comparison level 
    degs_data <- limma_value$top_genes[[1]]
    colnames(degs_data) <- c(
      (paste(limma_value$comparisons[[1]], "logFC", sep = "_")), 
      (paste(limma_value$comparisons[[1]], "adjPval", sep = "_"))
    )
    degs_data$ensembl_ID <- rownames(degs_data)
    
    # get the additional comparison levels if they exists 
    if (length(names(limma_value$top_genes)) > 1){
      for (i in 2:length(names(limma_value$top_genes))){
        temp <- limma_value$top_genes[[i]]
        colnames(temp) <- c(
          (paste(limma_value$comparisons[[i]], "logFC", sep = "_")),
          (paste(limma_value$comparisons[[i]], "adjPval", sep = "_"))
        )
        temp$ensembl_ID <- rownames(temp)
        degs_data <- dplyr::inner_join(degs_data, temp, by = "ensembl_ID")
      }
    }
    
    # connect to gene symbols and original user id 
    processed_data <- as.data.frame(processed_data)
    processed_data$ensembl_ID <- rownames(processed_data)
    
    degs_data <- dplyr::full_join(degs_data, gene_names, by = "ensembl_ID")
    degs_data <- dplyr::full_join(degs_data, processed_data, by = "ensembl_ID")
    
    
    degs_data <- degs_data %>%
      dplyr::relocate(User_ID) %>%
      dplyr::relocate(ensembl_ID) %>%
      dplyr::relocate(symbol)
  }
  return(list(degs_data, limma_value$Results))
}