#' fct_analysis_random.R Miscellaneous data analysis functions here,
#'  we find a place for them later
#'
#'
#' @section fct_analysis_random.R functions:
#' \code{gene_group_heatmap} heatmap with color bar define gene groups
#'
#' \code{find_overlap_gmt} Given a gene set,
#' finds significant overlaps with a gene set database.
#'
#'
#' @name fct_analysis_random.R
NULL

### work in process (Needs to be tested)
#' Calculate overlap for species with GMT file
#'
#' For a species not in the iDEP database, calculate the overlap
#'  for pathway analysis.
#'
#' @param query Vector of genes to calculate the overlap on
#' @param gene_set Gene sets from the custom GMT file to use in
#'  the overlap
#' @param min_fdr Significant p-value to determine significantly
#'  expressed genes
#' @param min_size Minimum size for a pathway gene set
#' @param max_size Maximum size for a pathway gene set
#' 
#' @export
find_overlap_gmt <- function(
  query, 
  gene_set,
  min_fdr = .2,
  min_size = 2,
  max_size = 10000
) {
  total_elements <- 30000 # why 3000?
  min_overlap <- 1 # nolint
  max_terms <- 10 # max number of enriched terms should be user input ????
  no_sig <- as.data.frame("No significant enrichment found!")
  query <- clean_gene_set(gene_set = query) # convert to upper case, unique()
  query_length <- length(query)

  if (query_length <= 2 || length(gene_set) < 1) {
    return(no_sig)
  }

  # gene sets smaller than 1 is ignored!!!
  gene_set <- gene_set[which(sapply(X = gene_set, FUN = length) > min_size)]
  # gene sets smaller than 1 is ignored!!!
  gene_set <- gene_set[which(sapply(X = gene_set, FUN = length) < max_size)]

  foo <- function(x, query) {
    length(intersect(query, x))
  }
  result <- unlist(lapply(X = gene_set, FUN = foo, query = query))
  result <- cbind(unlist(lapply(X = gene_set, FUN = length)), result)
  result <- result[which(result[, 2] > min_overlap), , drop = FALSE]
  if (nrow(result) == 0) {
    return(no_sig)
  }
  xx <- result[, 2]
  nn <- total_elements - query_length
  kk <- result[, 1]
  pval_enrich <- phyper(xx - 1, query_length, nn, kk, lower.tail = FALSE)
  fdr <- p.adjust(pval_enrich, method = "fdr", n = length(gene_set))
  result <- as.data.frame(cbind(fdr, result))
  result <- result[, c(1, 3, 2)]
  result$pathway <- rownames(result)
  result$Genes <- "" # place holder just
  colnames(result) <- c(
    "Corrected P value (FDR)",
    "Genes in list", "Total genes in category", "Functional Category", "Genes"
  )
  result <- result[which(result[, 1] < min_fdr), , drop = FALSE]
  if (nrow(result) == 0 || min(fdr) > min_fdr) {
    return(no_sig)
  }
  result <- result[order(result[, 1]), ]
  if (nrows(result) > max_terms) {
    result <- result[1:max_terms, ]
  }
  return(result)
}

#' Find overlap for pathway analysis
#' 
#' Use the pathway table from the read_pathway_sets function
#' to calculate adjusted p-values. Adjusted p-values determine
#' the enriched pathways from the selected qeury.
#' 
#' @param pathway_table Results from the read_pathway_sets
#'   query. If NULL or 0 rows there is no significant
#'   enrichment
#' @param query_set The vector of IDs that the enrichment
#'   analysis is being performed on
#' @param total_genes Length of the query set subtracted from
#'   the total number of genes in the database. Could change
#'   within the function if the background set changes to the
#'   filtered genes.
#' @param processed_data Data that has been filtered and
#'   transformed in the pre_process function
#' @param gene_info The gene info from the converted IDs and
#'   the function gene_info()
#' @param go Section of the database to query for pathway
#'   analysis
#' @param idep_data Data built in to idep
#' @param sub_pathway_files GMT files in the database that
#'   contain information for the matched species
#' @param use_filtered_background T/F Use the genes that
#'   passed the pre_process filter as the backgrounf
#' @param select_org Input for what organism the IDs are 
#'   pertaining to
#' @param reduced T/F Remove gene sets that are redudant
#'   from the final result
#'   
#' @export
#' @return If there is significant enrichment, the data frame
#'   that is returned has a pathway for each row with the
#'   total genes in the database mapping to it as well as the
#'   number of genes in the query that map to it. It also
#'   contains a column for the p-value and a list of the
#'   specific IDs included in the pathway from the query.
find_overlap <- function(
  pathway_table,
  query_set,
  total_genes,
  processed_data,
  gene_info,
  go,
  idep_data,
  sub_pathway_files,
  use_filtered_background,
  select_org,
  reduced = FALSE
) {
  max_pval_filter <- 0.3
  max_genes_background <- 30000
  min_genes_background <- 2000
  max_terms <- 15
  min_fdr <- .05
  if(reduced) {
    reduced <- .9
  }
  error_msg <- data.frame("Enrichment" = "No significant enrichment found!")

  if(select_org == "NEW" && is.null(pathway_table)) {
    return(data.frame("Enrichment" = "No GMT file provided!"))
  } else if (select_org == "NEW") {
    find_overlap_gmt(
      query = query_set,
      gene_set = pathway_table,
      min_fdr = .2,
      min_size = 2,
      max_size = 10000
    )
  }
  
  # pathway_table <- pathway_table[pathway_table$overlap > 1, ]

  if(dim(pathway_table)[1] == 0 || is.null(pathway_table)) {
    return(error_msg)
  }

  pathway_table <- pathway_table[which(
    pathway_table$overlap / length(query_set) /
    (as.numeric(pathway_table$n) / total_genes) > 1
  ), ]

  pathway_table$pval <- stats::phyper(
    pathway_table$overlap - 1,
		length(query_set),
		total_genes - length(query_set),   
		as.numeric(pathway_table$n), 
		lower.tail=FALSE
  )
  pathway_table <- subset(pathway_table, pathway_table$pval < max_pval_filter)

  # Background genes -----------
  if(!is.null(use_filtered_background)) {
    if(
      use_filtered_background && 
      length(row.names(processed_data)) > min_genes_background &&
      length(row.names(processed_data)) < max_genes_background + 1
    ) {
      pathway_table_bg <- background_pathway_sets(
        processed_data = processed_data,
        gene_info = gene_info,
        sub_query = query_set,
        go = go,
        pathway_table = pathway_table,
        idep_data = idep_data,
        sub_pathway_files = sub_pathway_files
      )

      pathway_table$pval <- phyper(
        pathway_table_bg$overlap - 1,
        length(query_set),
        length(rownames(processed_data)) - length(query_set),   
        as.numeric(pathway_table_bg$overlap_bg),
        lower.tail=FALSE
      ) 
    }
  }

  pathway_table$fdr <- stats::p.adjust(pathway_table$pval, method = "fdr")
  pathway_table <- pathway_table[order(pathway_table$fdr), ]

  if(dim(pathway_table)[1] > max_terms) {
    pathway_table <- pathway_table[1:max_terms, ]
  }
  
  if(min(pathway_table$fdr) > min_fdr) {
    pathway_table <- error_msg
  } else {
    pathway_table <- pathway_table[which(pathway_table$fdr < min_fdr), ]

    pathway_table <- subset(pathway_table, select = c(fdr, overlap, n, description, gene_sets))
    pathway_table$n <- as.numeric(pathway_table$n)
    pathway_table$fdr <-  formatC(pathway_table$fdr, format = "e", digits = 2)
    colnames(pathway_table) <- c(
      "Corrected P value (FDR)", "Genes in query", "Total genes in category",
      "Functional Category", "Genes"
    )

    # Remove redudant gene sets; only do it when there are more than 5. Error when there is only 1 or 2.
		if(reduced != FALSE && dim(pathway_table)[1] > 5){
			n <- nrow(pathway_table)
			tem <- rep(TRUE, n)
			gene_lists = pathway_table$Genes
			for(i in 2:n) {
        for(j in 1:(i-1)) { 
				  if(tem[j]) {
					  common_genes = length(intersect(gene_lists[[i]], gene_lists[[j]]))
					  if(common_genes/ length(gene_lists[[j]]) > reduced) {
              tem[j] = FALSE
            }	
				  }			
				}				
      }
								
			pathway_table <- pathway_table[which(tem), ]		
		}
  }

  return(pathway_table)
}

#' Determine samples in selected contrast
#' 
#' Find the samples that are in the group of the selected
#' contrast. Can be used to subset the data to only include
#' the samples that correspond to the chosen comparison.
#' 
#' @param select_contrast Comparison from DEG analysis to filter
#'  for the significant genes
#' @param all_sample_names The column names of the processed data
#' @param sample_info Experiment file information for grouping
#'  samples
#' @param select_factors_model The selected factors for the model
#'  expression
#' @param select_model_comprions Selected comparisons to analyze
#'  in the DEG analysis
#' @param reference_levels Vector of reference levels to use for the
#'  selected factors
#' @param counts_deg_method Method of DEG being performed (See
#'  DEG UI for options)
#' @param data_file_format Type of data being examined
#' 
#' @export
#' @return A numeric vector that can be used to index the processed
#'  data and subset to only include the columns from the selected
#'  contrast.
find_contrast_samples <- function(
  select_contrast,
  all_sample_names,
  sample_info = NULL,
  select_factors_model = NULL,
  select_model_comprions = NULL,
  reference_levels = NULL,
  counts_deg_method = NULL,
  data_file_format = NULL
){
	iz <- match(detect_groups(all_sample_names), unlist(strsplit(select_contrast, "-")))
	iz <- which(!is.na(iz))
	
	# Has design file, but didn't select factors
	if(!is.null(sample_info) & is.null(select_factors_model) &
     length(select_model_comprions) == 0) {
	  
    find_samples <- function(
      factor_level,
      sample_info
    ) { 
	    # Given a factor level such as "wt", return a vector indicating the samples 
      # with TRUE FALSE
	    #  p53_mock_1  p53_mock_2  p53_mock_3  p53_mock_4 ... p53_IR_4 null_mock_1 null_mock_2 
	    #  TRUE        TRUE        TRUE        TRUE       ... TRUE     FALSE       FALSE 
  	  tem <- apply(sample_info, 2, function(y) y == factorLevel )
  	  colSums(tem) > 0
  	  tem <- tem[, colSums(tem) > 0]
  	  return(tem)
	  }
	  
	  
	  sample_1 <- gsub("-.*", "", select_contrast)
	  level_1 <- gsub("_.*", "", sample_1)
	  level_2 <- gsub(".*_", "", sample_1)
	  iz <- which(find_samples(level_1, sample_info) & find_samples(level_2, sample_info))
	  
	  sample_2 <- gsub(".*-", "", select_contrast)
	  level_1 <- gsub("_.*", "", sample_2)
	  level_2 <- gsub(".*_", "", sample_2)
	  iz <- c(iz, which(find_samples(level_1, sample_info) & find_samples(level_2, sample_info)))	  
	}
	
	#Has design file and chose factors
	if(!is.null(sample_info) & !is.null(select_factors_model) &
     length(select_model_comprions) > 0) {
    
    # Strings like: "groups: mutant vs. control"
		comparisons <- gsub(".*: ", "", select_model_comprions)
    comparisons <- gsub(" vs\\. ", "-", comparisons)		
    # Corresponding factors
		factors_vector <- gsub(":.*", "", select_model_comprions)

	  # If read counts data and DESeq2
		if(data_file_format == 1 & counts_deg_method == 3) {
      # Could be "wt-mu" or "wt-mu_for_conditionB"
			contrast <- gsub("_for_.*", "", select_contrast)
      # Selected contrast lookes like: "mutant-control" 
			ik <- match(contrast, comparisons)

			other_factor_level <- gsub(".*_for_", "", select_contrast)
			# Find the corresponding factor for the other factor
			other_factor <- " "
			if(nchar(other_factor_level) > 0) {
				for(each_factor in colnames(sample_info)) {
          if (other_factor_level %in% sample_info[, each_factor]) {
            other_factor <- each_factor
          }
        }		
			}
			
			if(is.na(ik)) {
        iz <- 1:(length(all_sample_names))
      # Interaction term, use all samples
      } else {
        # Corresponding factors
				selected_factor <- factors_vector[ik] 

				iz <- which(sample_info[, selected_factor] %in% unlist(strsplit(contrast, "-")))

				# Filter by other factors: reference level
        # c("genotype:wt", "treatment:control")
				if(!is.null(reference_levels)) {   
					for(refs in reference_levels) {
            if(!is.null(refs) & gsub(":.*", "", refs) != selected_factor) {
							current_factor <- gsub(":.*", "", refs)
              # If not reference level
							if(nchar(other_factor_level) > 0 & current_factor == other_factor) {
								iz <- intersect(iz, which(sample_info[, current_factor] == other_factor_level))
							} else {
                iz <- intersect(iz, which(sample_info[, current_factor] == gsub(".*:", "", refs)))
              }	
						}
          }
				}			
				iz <- iz[which(!is.na(iz))]
			# Switching from limma to DESeq2 causes problem, as reference level is not defined.
			}
    # Not DESeq2
		} else {
		  # Given level find corresponding sample ids
			find_ids_from_level <- function(
        a_level,
        sample_info
      ){
				# Find factor
				current_factor = ""
				for(each_factor in colnames(sample_info)) {
          if (a_level %in% sample_info[, each_factor]) {
            current_factor <- each_factor
          }
        }			
				
        if(nchar(current_factor) > 0) {
          return(which(sample_info[, current_factor] %in% a_level)) 
        } else {
          return(NULL)
        }	
      }		
					
			if(!grepl(".*_.*-.*_.*", select_contrast)) {
        iz <- c()
      }
      # Double split!
			levels_4 <- unlist(strsplit(unlist(strsplit(select_contrast, "-")), "_")) 
			if(length(levels_4) != 4) { 
				iz <- c() 
			} else {
        # First sample
				iz <- intersect(
          find_ids_from_level(levels_4[1], sample_info),
          find_ids_from_level(levels_4[2], sample_info)
        )
        # 2nd sample
				iz <- c(
          iz,
          intersect(
            find_ids_from_level(levels_4[3], sample_info),
            find_ids_from_level(levels_4[4], sample_info)
          )
        ) 
			}
		}
	}	
	 
	if(grepl("I:", select_contrast)) {
    # If it is factor design use all samples
    iz <- 1:length(all_sample_names)
  }
	if(is.na(iz)[1] | length(iz) <= 1 ) {
    iz <- 1:length(all_sample_names)
  }

	return(iz)
}

#' Read gene sets GMT file
#' This functions cleans and converts to upper case
#' @export
read_gmt_robust <- function (in_file) {
	# Read in the first file 
	x <- scan(in_file, what = "", sep = "\n")
  # GMT files saved by Excel has a lot of empty cells "\t\t\t\t" "\t." means one or more tab
	# Remove white space
  x <- gsub(" ", "", x)  
	# Convert to upper case
  x <- toupper(x)

	#----Process the first file
	# Separate elements by one or more whitespace
	y <- strsplit(x, "\t")
	# Extract the first vector element and set it as the list element name
	names(y) <- sapply(y, "[[", 1)
	# Remove the first vector element from each list element
	y <- lapply(y, "[", -c(1, 2))
	# Remove duplicated elements
	for(i in 1:length(y)) {
    y[[i]] <- clean_gene_set(y[[i]])
  }
	# Check the distribution of the size of gene lists sapply(y, length) hold a vector of sizes
	if(max(sapply(y, length)) < 5) {
    cat("Warning! Gene sets have very small number of genes!\n Please double check format.")
  }
  # Gene sets smaller than 1 is ignored!!!
	y <- y[which(sapply(y, length) > 1)]

	return(y)
}

# retrieve detailed info on genes
get_gene_info <- function(
  converted,
  select_org,
  gene_info_files
) {
	if(is.null(converted)) {
    return(as.data.frame("ID not recognized!"))
  } 
	query_set <- converted$ids
	if(length(query_set) == 0) {
    return(as.data.frame("ID not recognized!"))
  }
	ix <- grep(converted$species[1, 1], gene_info_files)
	if(length(ix) == 0) {
    return(as.data.frame("No matching gene info file found") )
  } else {
	  # If selected species is not the default "bestMatch", use that species directly
	  if(select_org != "BestMatch") {  
		  ix <- grep(find_species_by_id(select_org)[1, 1], gene_info_files)
	  }
	  if(length(ix) == 1) {
      x <- read.csv(as.character(gene_info_files[ix])) 
      x[, 1] <- toupper(x[, 1]) 
      # If symbol is missing use Ensembl IDs
      x$symbol[is.na(x$symbol) ] <- x[, 1]
      # If duplicated symbol, paste Ensembl id to the end
      n_occur <- data.frame(table(x$symbol))
      # Rows with duplicated symbols
      ix_duplicated <- which(n_occur$Freq > 1) 
      x$symbol[ix_duplicated] <- paste(x$symbol[ix_duplicated], x[ix_duplicated, 1])
    } else {
      # Read in the chosen file 
      return(as.data.frame("Multiple geneInfo file found!"))
    }
	  Set <- match(x$ensembl_gene_id, query_set)
	  Set[which(is.na(Set))] = "Genome"
	  Set[which(Set!="Genome")] = "List"

	  return(cbind(x, Set))
  }
}
