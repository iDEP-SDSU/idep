#' fct_01_load_data.R This file holds all of the main data analysis functions
#' associated with first tab of the iDEP website.
#'
#'
#' @section fct_01_load_data.R functions:
#'
#'
#' @name fct_01_load_data.R
NULL

#' Retrieve detailed info on genes
#' 
#' This function retrieves detailed gene information from the
#' database for the matched species. 
#' 
#' @param converted Return value from the convert_id function. Contains
#'  information about the gene IDs for the matched species.
#' @param select_org The user selected organism for the expression data.
#'  Default is "BestMatch."
#' @param idep_data Read data files from the database. 
#' 
#' @export
#' @return A data frame containing information on all ensembl IDs for the
#'  matched species. Usually a very large data frame due to the amount
#'  of IDs that the data base contains.
gene_info <- function(
  converted,
  select_org,
  idep_data
) {
  check <- check_object_state(
    check_exp = (is.null(converted)),
    true_message = as.data.frame("ID not recognized!")
  )
  if (check$bool) {
    return(check)
  }

  query_set <- converted$ids
  check <- check_object_state(
    check_exp = (length(query_set) == 0),
    true_message = as.data.frame("ID not recognized!")
  )
  if (check$bool) {
    return(check)
  }

  ix <- grep(
    pattern = converted$species[1, 1],
    x = idep_data$gene_info_files
  )
  check <- check_object_state(
    check_exp = (length(ix) == 0),
    true_message = as.data.frame("No matching gene info file found")
  )
  if (check$bool) {
    return(check)
  }

  # If selected species is not the default "bestMatch",
  # use that species directly
  if (select_org != idep_data$species_choice[[1]]) {
    ix <- grep(
      pattern = find_species_by_id(
        species_id = select_org,
        org_info = idep_data$org_info
      )[1, 1],
      x = idep_data$gene_info_files
    )
  }

  check <- check_object_state(
    check_exp = (length(ix) != 1),
    true_message = as.data.frame("Multiple geneInfo file found!")
  )
  if (check$bool) {
    return(check)
  } else {
    gene_info_csv <- read.csv(as.character(idep_data$gene_info_files[ix]))
    gene_info_csv[, 1] <- toupper(gene_info_csv[, 1])
  }

  set <- match(gene_info_csv$ensembl_gene_id, query_set)
  set[which(is.na(set))] <- "Genome"
  set[which(set != "Genome")] <- "List"
  return(cbind(gene_info_csv, set))
}




#' Load basic data information
#'
#' This function does the immediate loading of the data and
#' sample info to present in the data preview table and the
#' sample info table. The data undergoes very basic filtering
#' and transformation before entering the table.
#'
#' @param expression_file The data path for the expression file, should be
#'  accessed with \code{expression_file$datapath}
#' @param experiment_file The data path for the experiment file, should be
#'  accessed with \code{experiment_file$datapath}
#' @param go_button TRUE/FALSE that tells the app to
#' load the demo data files
#' @param demo_data_file Expression demo data path (idep_data$demo_data_file)
#' @param demo_metadata_file Experiment demo data path 
#'  (idep_data$demo_metadata_file)
#'
#' @export
#' @return This returns a list that contains the expression data
#' and the sample information. If there is no experiment file it
#' only returns the expression data.
#'
input_data <- function(
  expression_file,
  experiment_file,
  go_button,
  demo_data_file,
  demo_metadata_file
) {
  in_file_data <- expression_file
  in_file_data <- in_file_data$datapath

  if (is.null(in_file_data) && go_button == 0) {
    return(NULL)
  } else if (go_button > 0) {    # use demo data
    in_file_data <- demo_data_file
  }

  isolate({
    # Read expression file -----------
    data <- read.csv(in_file_data, quote = "", comment.char = "")
    # Tab-delimented if not CSV
    if (ncol(data) <= 2) {
      data <- read.table(
        in_file_data,
        sep = "\t",
        header = TRUE,
        quote = "",
        comment.char = ""
      )
    }

    # Filter out non-numeric columns ---------
    num_col <- c(TRUE)
    for (i in 2:ncol(data)) {
      num_col <- c(num_col, is.numeric(data[, i]))
    }
    if (sum(num_col) <= 2) {
      return(NULL)
    }
    data <- data[, num_col]

    # Format gene ids --------
    data[, 1] <- toupper(data[, 1])
    data[, 1] <- gsub(" |\"|\'", "", data[, 1])

    # Remove duplicated genes ----------
    data <- data[!duplicated(data[, 1]), ]

    # Set gene ids as rownames and get rid of column ---------
    rownames(data) <- data[, 1]
    data <- as.matrix(data[, c(-1)])

    # Remove "-" or "." from sample names ----------
    colnames(data) <- gsub("-", "", colnames(data))
    colnames(data) <- gsub("\\.", "", colnames(data))

    # Order by SD ----------
    data <- data[order(-apply(
      data[, 1:ncol(data)],
      1,
      sd
    )), ]
  })

  # Read experiment file ----------
  in_file_expr <- experiment_file
  in_file_expr <- in_file_expr$datapath
  if (is.null(in_file_expr) && go_button == 0) {
    return(list(
      data = data,
      sample_info = NULL
    ))
  } else if (go_button > 0) {
    sample_info_demo <- t(read.csv(
      demo_metadata_file,
      row.names = 1,
      header = T,
      colClasses = "character"
    ))
    return(list(
      sample_info = sample_info_demo,
      data = data
    ))
  }

  isolate({
    # Read experiment file ----------
    expr <- read.csv(
      in_file_expr,
      row.names = 1,
      header = TRUE,
      colClasses = "character"
    )
    if (ncol(expr) <= 2) {
      expr <- read.table(
        in_file_expr,
        row.names = 1,
        sep = "\t",
        header = TRUE,
        colClasses = "character"
      )
    }
    # remove "-" or "." from sample names ----------
    colnames(expr) <- gsub("-", "", colnames(expr))
    colnames(expr) <- gsub("\\.", "", colnames(expr))

    # Matching with column names of expression file ----------
    matches <- match(
      toupper(colnames(data)), toupper(colnames(expr))
    )
    matches <- matches[which(!is.na(matches))] # remove NA
    validate(need(
      length(unique(matches)) == ncol(data) &&
        nrow(expr) >= 1 && nrow(expr) < 500,
      "Error!!! Sample information file not recognized. Sample names
       must be exactly the same. Each row is a factor. Each column
       represent a sample.  Please see documentation on format."
    ))

    # Check factor levels, change if needed ----------
    for (i in 1:nrow(expr)) {
      expr[i, ] <- gsub("-", "", expr[i, ])
      expr[i, ] <- gsub("\\.", "", expr[i, ])
      # Remove whitespace
      expr[i, ] <- gsub(" ", "", expr[i, ]) 
      # Convert to upper case to avoid mixing
      expr[i, ] <- toupper(expr[i, ])    
    }

    # Factor levels match ---------
    if (length(unique(matches)) == ncol(data)) {
      expr <- expr[, matches]
      if (
        sum(apply(expr, 1, function(y) length(unique(y)))) >
          length(unique(unlist(expr)))) {
        factor_names <- apply(
          expr,
          2,
          function(y) paste0(names(y), y)
        )
        rownames(factor_names) <- rownames(expr)
        expr <- factor_names
      }
      return(list(
        data = data,
        sample_info = t(expr)
      ))
    } else {
      return(list(
        data = data
      ))
    }
  })
}

#' Convert expression data rownames to ensembl
#'
#' This function works with the converted ids to take in the
#' expression data and swap the rownames to the converted ids.
#' It returns the data exactly the same except for the changed
#' ids.
#'
#' @param converted Data from convert_id function containing converted ids
#' @param no_id_conversion TRUE/FALSE for converting data ids or not
#' @param data Data from inputed expression file
#' @export
#' @return Returns original data with rownames converted to ensembl
convert_data <- function(
  converted,
  data,
  no_id_conversion
) {
  if (is.null(converted) || no_id_conversion) {
    return(list(
      data = data,
      mapped_ids = rownames(data)
    ))
  } else {
    mapping <- converted$conversion_table

    rownames(data) <- toupper(rownames(data))
    merged <- merge(
      mapping[, 1:2],
      data,
      by.y = "row.names",
      by.x = "User_input",
      all.y = TRUE
    )
    no_match <- which(is.na(merged[, 2]))
    merged[no_match, 2] <- merged[no_match, 1]

    mapped_ids <- merged[, 1:2]

    # Multiple matches use one with highest SD ----------
    tmp <- apply(merged[, 3:(ncol(merged))], 1, sd)
    merged <- merged[order(merged[, 2], -tmp), ]
    merged <- merged[!duplicated(merged[, 2]), ]
    rownames(merged) <- merged[, 2]
    merged <- as.matrix(merged[, c(-1, -2)])

    # Order by SD ----------
    merged <- merged[order(-apply(
      merged[, 1:ncol(merged)],
      1,
      sd
    )), ]

    return(list(
      data = merged,
      mapped_ids = mapped_ids
    ))
  }
}

#' Get all matched gene names
#'
#' This function will create a data frame with all the matched
#' gene names from the database.
#'
#' @param mapped_ids Matched IDs from the convert_data functions
#' @param all_gene_info Gene information matched species in idep database
#' @export
#' @return Data frame containing all the matched ID names from idep
#' database. Three columns denotes a recognized species for which
#' idep had gene names for. Two columns means the IDs were converted
#' to ensembl format, but no species was found for the gene names.
#' One means no conversion occurred.
get_all_gene_names <- function(
  mapped_ids,
  all_gene_info
) {
  if (is.null(dim(mapped_ids))) {
    return(data.frame("User_ID" = mapped_ids))
  } else if (!is.null(all_gene_info$bool)) {
    return(data.frame(
      "User_ID" = mapped_ids[, 1],
      "ensembl_ID" = mapped_ids[, 2]
    ))
  } else {
    mapped_ids <- data.frame(
      "User_ID" = mapped_ids[, 1],
      "ensembl_ID" = mapped_ids[, 2]
    )
    all_names <- merge(
      mapped_ids,
      all_gene_info[, c("ensembl_gene_id", "symbol")],
      by.x = "ensembl_ID",
      by.y = "ensembl_gene_id",
      all.x = T
    )
    all_names <- all_names[!duplicated(all_names$ensembl_ID), ]
    all_names$symbol[all_names$symbol == ""] <- NA
    all_names$symbol[is.na(all_names$symbol)] <- {
      all_names$ensembl_ID[is.na(all_names$symbol)]
    }
    duplicates <- all_names$symbol %in%
      all_names$symbol[duplicated(all_names$symbol)]
    all_names$symbol[duplicates] <- paste(
      all_names$symbol[duplicates],
      all_names$ensembl_ID[duplicates]
    )
    all_names <- dplyr::select(
      all_names,
      User_ID,
      ensembl_ID,
      symbol,
      tidyselect::everything()
    )

    return(all_names)
  }
}