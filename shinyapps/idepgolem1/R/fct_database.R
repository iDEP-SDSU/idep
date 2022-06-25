#' fct_database.R: A fie with all functions to deal with database
#'
#'
#' @section fct_database.R functions:
#' \code{connect_convert_db} connects to the converID.db
#' and returns the objects.
#'
#'
#'
#' @name fct_database.R
NULL


DATAPATH <- "../../data/data104/"
#DATAPATH <- "F:/data104b_final/"

#' connect_convert_db connects to the convertIDs.db and returns the
#' objects.
#'
#' Create a database connection with the DBI package.
#'
#' @param datapath Folder path to the data file
#'
#' @export
#' @return Database connection.
connect_convert_db <- function(datapath = DATAPATH) {
  return(DBI::dbConnect(
    drv = RSQLite::dbDriver("SQLite"),
    dbname = paste0(datapath, "convertIDs.db"),
    flags = RSQLite::SQLITE_RO
  ))
}


#' Load iDEP data
#'
#' Use this function call to load data that is
#' used in other functions.
#'
#' 
#' @param datapath Folder path to the iDEP data
#' @export
#' @return Large list of the iDEP data.
#' 1. kegg_species_id:  KEGG species list
#' 2. gmt_files: list of pathway files
#' 3. gene_info_files: list of geneInfo files
#' 4. demo_data_file: demo data file 
#' 5. demo_metadata_file: experimental design file for demo data
#' 6. quotes:  quotes
#' 7. string_species_go_data: List of STRING species
#' 8. org_info: orgInfo for Ensembl species
#' 9. annotated_species_count: total number of annotated species
#' 10. go_levels: GO levels
#' 11. go_level_2_terms: mapping of GO levels to terms
#' 12. id_index: idtype and index
#' 13. species_choice: list of species for populating selection.
get_idep_data <- function(datapath = DATAPATH) {

  # if prepared RData files exists, return the objects.
  # file is prepared with this command
  #  saveRDS(get_idep_data(), file="prepared_data.RData", compress=FALSE)
  # then this file needs to be copied to /data_go/
  # BE CAREFUL: the RDS file is not checked. It should be updated if the files change.
  if(file.exists(paste0(datapath, "data_go/prepared_data.RData"))) {
    return(readRDS(paste0(datapath, "data_go/prepared_data.RData")))
  } 
  # Code below will not be executed if RData file exists.

  kegg_species_id <- read.csv(paste0(datapath, "data_go/KEGG_Species_ID.csv"))

  gmt_files <- list.files(
    path = paste0(datapath, "pathwayDB"),
    pattern = ".*\\.db"
  )
  gmt_files <- paste(datapath,
    "pathwayDB/", gmt_files,
    sep = ""
  )

  gene_info_files <- list.files(
    path = paste0(datapath, "geneInfo"),
    pattern = ".*GeneInfo\\.csv"
  )
  gene_info_files <- paste(datapath,
    "geneInfo/", gene_info_files,
    sep = ""
  )

  demo_data_file <- paste0(datapath, "data_go/BcellGSE71176_p53.csv")
  demo_metadata_file <- paste0(
    datapath,
    "data_go/BcellGSE71176_p53_sampleInfo.csv"
  )

  conn_db <- connect_convert_db()

  quotes <- DBI::dbGetQuery(
    conn = conn_db,
    statement = "select * from quotes"
  )
  quotes <- paste0(
    "\"", quotes$quotes,
    "\"", " -- ", quotes$author, ".       "
  )

  string_species_go_data <- read.csv(paste0(
    datapath,
    "data_go/STRING11_species.csv"
  ))

  org_info <- DBI::dbGetQuery(
    con = conn_db,
    statement = "select distinct * from orgInfo"
  )
  org_info <- org_info[order(org_info$name), ]

  annotated_species_count <- sort(table(org_info$group))

  go_levels <- DBI::dbGetQuery(
    conn = conn_db,
    statement = "select distinct id, level from GO
         WHERE GO = 'biological_process'"
  )

  go_level_2_terms <- go_levels[which(go_levels$level %in% c(2, 3)), 1]

  id_index <- DBI::dbGetQuery(
    conn = conn_db,
    statement = "select distinct * from idIndex"
  )

  species_choice <- setNames(as.list(org_info$id), org_info$name2)
  species_choice <- append(
    setNames("NEW", "**NEW SPECIES**"),
    species_choice
  )
  species_choice <- append(
    setNames("BestMatch", "Best matching species"),
    species_choice
  )
  top_choices <- c(
    "Best matching species", "**NEW SPECIES**", "Human", "Mouse", "Rat", "Cow",
    "Zebrafish", "Pig", "Chicken", "Macaque", "Dog", "Drosophila melanogaster",
    "Caenorhabditis elegans", "Saccharomyces cerevisiae",
    "Arabidopsis thaliana", "Zea mays", "Glycine max",
    "Oryza sativa Indica Group", "Oryza sativa Japonica Group", "Vitis vinifera"
  )
  other_choices <- names(species_choice)[
    !(names(species_choice) %in% top_choices)
  ]
  species_choice <- species_choice[c(top_choices, other_choices)]
  

  DBI::dbDisconnect(conn = conn_db)

  return(list(
    kegg_species_id = kegg_species_id,
    gmt_files = gmt_files,
    gene_info_files = gene_info_files,
    demo_data_file = demo_data_file,
    demo_metadata_file = demo_metadata_file,
    quotes = quotes,
    string_species_go_data = string_species_go_data,
    org_info = org_info,
    annotated_species_count = annotated_species_count,
    go_levels = go_levels,
    go_level_2_terms = go_level_2_terms,
    id_index = id_index,
    species_choice = species_choice
  ))
}


#' Find the ID type
#'
#' Using an inputted ID, find the type of IDs
#'
#' @param id ID to find the type from
#' @param id_index Index of IDs to use
#' 
#' @export
find_id_type_by_id <- function(id, id_index) {
  # find idType based on index
  return(id_index$idType[as.numeric(id)])
}

#' Find a species by ID
#'
#' Find a species in the iDEP database with an
#' ID.
#'
#' @param species_id Species ID to search the database with
#' @param org_info iDEP data org_info file
#'
#' @export
#' @return Species information in \code{org_info} from the
#'  matched ID.
find_species_by_id <- function(species_id, org_info) {
  return(org_info[which(org_info$id == species_id), ])
}


#' Find a species by ID
#'
#' Find a species in the iDEP database with an
#' ID.
#' 
#' @param species_id Species ID to search the database with
#' @param org_info iDEP data org_info file
#'
#' @export
#' @return Only return the species name with this function.
find_species_by_id_name <- function(species_id, org_info) {
  # find species name use id
  return(org_info[which(org_info$id == species_id), 3])
}



#' convert_id This function takes gene IDs and converts them to ensembl data.
#'
#' FUNCTION_DESCRIPTION
#'
#' @param query A character vector of gene IDs
#' @param idep_data A instance of the output from get_idep_data
#'  (link to documentation)
#' @param select_org A character of the species that wants to be looked up,
#'  default to \code{"BestMatch"}
#'
#' @export
#' @return A large list of the conversion ID information that was gathered
#'  from querying the database with the original IDs.
convert_id <- function(
  query,
  idep_data,
  select_org = "BestMatch"
) {
  query <- gsub(pattern = "\"|\'", "", x = query)
  # remove " in gene ids, mess up SQL query
  # remove ' in gene ids
  # |\\.[0-9] remove anything after A35244.1 -> A35244
  #  some gene ids are like Glyma.01G002100
  query_set <- clean_query(query_input = query)
  query_string <- paste0("('", paste(query_set,collapse="', '"),"')")

	#use a small set of genes to guess species and idType; to improve speed
  # use a small set of gene ids, with the max #
  # when the query is small, use the quary
  max_sample_ids = 100 # default 
  # use the query by default
  n_sample_ids = length(query_set) # acutal number of samples, for calculating % later
  test_query_string <- query_string
	if(length(query_set) > max_sample_ids) {
    n_sample_ids <- max_sample_ids
	  test_query_set <- sample(query_set, max_sample_ids)
    test_query_string <- paste0("('", paste(test_query_set, collapse="', '"),"')")
  }

  conn_db <- connect_convert_db()

  # if best match species--------------------------------------------------
  if (select_org == idep_data$species_choice[[1]]) { 

	  #First send a query to determine the species
	  query_species <- paste0( "select species, idType, COUNT(species) 
      as freq from mapping where id IN ", 
	    test_query_string," GROUP by species,idType"
    )
	  species_ranked <- DBI::dbGetQuery(conn_db, query_species)
	  if(dim(species_ranked)[1] == 0) {
      return(NULL)
    }

    species_ranked <- species_ranked[order(-species_ranked$freq), ]
    #  species idType freq
    #     131      1   99
    #     131     90   87
    #     -10090    409   86    # negative is STRING species
    #     -10090    410   86

    # try to use Ensembl species if it ranked a close 2nd 
    if(nrow(species_ranked) > 1) { # if more than one matched
      #if 2nd is close to the first,       
      if (species_ranked$freq[1] <= species_ranked$freq[2] * 1.2 
        && species_ranked$species[1] < 0 # the first is STRING 
        && species_ranked$species[2] > 0 # the 2nd is Ensembl
      ) {
          #swap 2nd with 1st
          species_ranked[1:2, ] <- species_ranked[c(2, 1), ]
      }
    }
    
    # add species name "Mouse"
    species_ranked$name <- sapply(
      X = species_ranked$species,
      FUN = find_species_by_id_name,
      org_info = idep_data$org_info
    )

    #Bind species and % matched genes such as "Mouse(93)"
    species_ranked$name <- paste0(species_ranked$name,
      " (", round(species_ranked$freq / n_sample_ids * 100, 0), ")",
      sep = ""
    )

    species_matched <- as.data.frame(species_ranked$name)

    #Query for mapping
    query_statement <- paste0(
      "select distinct id,ens,species,idType from mapping where ",
      " species = '", species_ranked$species[1], "'",
      " AND idType = '", species_ranked$idType[1], "'",
      " AND id IN ", query_string
    )
		result <- DBI::dbGetQuery(conn_db, query_statement)
		if( dim(result)[1] == 0 ){ 
      return(NULL)
    }
  } else { 
    # if species is selected ---------------------------------------------------
    query_statement <- paste0(
      "select distinct id,ens,species,idType from mapping where species = '", 
      select_org,
	    "' AND id IN ", 
      query_string
    )
	  result <- DBI::dbGetQuery(conn_db, query_statement)

    if (nrow(result) == 0) {
      return(NULL)
    } # stop("ID not recognized!")
    species_matched <- as.data.frame(paste(
      find_species_by_id_name(
        species_id = select_org,
        org_info = idep_data$org_info
      )
    ))
  }

  # remove duplicates in query gene ids
  result <- result[which(!duplicated(result[, 1])), ]
  # remove duplicates in ensembl_gene_id
  result <- result[which(!duplicated(result[, 2])), ]
  colnames(species_matched) <- c("Matched Species (%genes)")
  conversion_table <- result[, 1:2]
  colnames(conversion_table) <- c("User_input", "ensembl_gene_id")
  conversion_table$Species <- sapply(
    X = result[, 3],
    FUN = find_species_by_id_name,
    org_info = idep_data$org_info
  )

  species <- find_species_by_id(
    species_id = result$species[1],
    org_info = idep_data$org_info
  )

  return(list(
    origninal_ids = query_set,
    ids = unique(result[, 2]),
    species = species,
    species_matched = species_matched,
    conversion_table = conversion_table
  ))
}

#' Read pathway sets for gene query
#' 
#' This function provides the gene set information
#' to perform enrichment analysis on the provided
#' query.
#' 
#' @param all_gene_names_query Subsetted rows of the 
#'   all_gene_names data frame to query the database with
#' @param converted Conversion information from the original
#'   IDs returned from the convert_id() function
#' @param go Section of the database to query for pathway
#'   analysis
#' @param select_org Input for what organism the IDs are 
#'   pertaining to
#' @param gmt_file For NEW species the gmt file to use
#'   for the pathway analysis
#' @param idep_data Data built in to idep
#' @param gene_info The gene info from the converted IDs and
#'   the function gene_info()
#' 
#' @export
#' @return This function returns a list with values that are
#'   used in the find_overlap function. The list contains
#'   pathway_table which is the overlap and total genes for
#'   each pathway that is enriched in the query. The list
#'   also contains the query_set of genes, the total_genes
#'   number which is used in the calculation of the p-values
#'   and the pathway files that contain gmt information on
#'   the mathced species.
read_pathway_sets <- function (
  all_gene_names_query,
  converted,
  go,
  select_org,
  gmt_file,
  idep_data,
  gene_info
) {
	id_not_recognized = as.data.frame("ID not recognized!")

  if(select_org == "NEW" && is.null(gmt_file)) {
    return(as.data.frame("No GMT file provided!"))
  } else if (select_org == "NEW" && !is.null(gmt_file)) {
    in_file <- gmt_file
    in_file <- in_file$datapath

    return(read_gmt(in_file))
  }

	if(ncol(all_gene_names_query) == 1) {
    return(id_not_recognized)
  }

  query_set <- all_gene_names_query[, 2]

  if(!is.null(gene_info)) {
    if(dim(gene_info)[1] > 1) {  
	    gene_info <- gene_info[which(
        gene_info$gene_biotype == "protein_coding"
      ), ]  
	    query_set <- intersect(query_set, gene_info[, 1])
	  }
  }

  if(length(query_set) == 0) {
    return(id_not_recognized)
  } 

	ix = grep(converted$species[1, 1], idep_data$gmt_files)
  total_genes <- converted$species[1, 7]

	# If selected species is not the default "bestMatch", use that species directly
	if(select_org != "BestMatch") {  
		ix = grep(
      find_species_by_id(select_org, idep_data$org_info)[1, 1],
      idep_data$gmt_files
    )
		total_genes <- idep_data$org_info[which(
      idep_data$org_info$id == as.numeric(select_org)
    ), 7]
	}

  if (length(ix) == 0) {
    return(id_not_recognized )
  }
  
  pathway_files <- idep_data$gmt_files[ix]
	pathway <- DBI::dbConnect(
    drv = RSQLite::dbDriver("SQLite"),
    dbname = pathway_files,
    flags= RSQLite::SQLITE_RO
  )
	
	if(is.null(go)) {
    go <- "GOBP"
  }

	sql_query <- paste(
    "select distinct gene, pathwayID from pathway where gene IN ('",
    paste(query_set, collapse = "', '"), "')",
    sep = ""
  )

	if(go != "All") {
    sql_query = paste0(sql_query, " AND category ='", go,"'")
  }
	result <- DBI::dbGetQuery(pathway, sql_query)

	if(dim(result)[1] == 0) {
    return(pathway_table <- NULL)
  }

	# List pathways and frequency of genes
	pathway_ids <- stats::aggregate(
    result$pathwayID,
    by = list(unique_values = result$pathwayID),
    FUN = length
  )
  colnames(pathway_ids) <- c("pathway_id", "overlap")

	if(dim(pathway_ids)[1] == 0) {
    return(pathway_table <- NULL)
  } else {
    # Convert pathways into lists
	  gene_sets <- lapply(
      pathway_ids[, 1],
      function(x)  result[which(result$pathwayID == x), 1]
    )
	  pathway_info <- DBI::dbGetQuery(
      pathway,
      paste(
        "select distinct id, n, Description from pathwayInfo where id IN ('", 
				paste(pathway_ids[, 1], collapse = "', '"), "') ", sep = ""
      ) 
    )
	  ix <- match(pathway_ids[, 1], pathway_info[, 1])
    pathway_merge <- merge(
      x = pathway_ids,
      y = pathway_info,
      by.x = "pathway_id",
      by.y = "id"
    )
	  names(gene_sets) <- pathway_info[ix, 1]

    pathway_merge$gene_sets <- gene_sets
  }
	
	DBI::dbDisconnect(pathway)

  # Gene sets and info for the enrichment analysis
	return(list(
    pathway_table = pathway_merge,
    query_set = query_set,
    total_genes = total_genes,
    pathway_files = pathway_files
  ))
}

#' Background gene pathway sets
#' 
#' This function reads the pathway sets for the filtered
#' gene IDs and performs pathway analysis for the query
#' with the filtered background.
#' 
#' @param processed_data The data matrix that has been through
#'   the pre-processing function
#' @param gene_info The gene info from the converted IDs and
#'   the function gene_info()
#' @param sub_query Query that the pathway analysis is being
#'   performed for
#' @param go Section of the database to query for pathway
#'   analysis
#' @param pathway_table Table from the initial querying for the
#'  sub-query
#' @param idep_data Data built in to idep
#' @param sub_pathway_files The subset of GMT files that contain
#'   information for the matched species
#' 
#' @export
#' @return Pathway gene set table for the background genes. Used
#'   in find_overlap to calculate pvals for the filtered background
background_pathway_sets <- function(
  processed_data,
  gene_info,
  sub_query,
  go,
  pathway_table,
  idep_data,
  sub_pathway_files
){
  query_set <- rownames(processed_data)
  
  if(!is.null(gene_info)) {
    if(dim(gene_info)[1] > 1) {
	    query_set <- intersect(
        query_set, 
        gene_info[which(gene_info$gene_biotype == "protein_coding"), 1]
      )
	  }  
  }

  pathway <- DBI::dbConnect(
    drv = RSQLite::dbDriver("SQLite"),
    dbname = sub_pathway_files,
    flags= RSQLite::SQLITE_RO
  )

  if(length(intersect(query_set, sub_query) ) == 0) {
    # If none of the selected genes are in background genes
    return(list(
      bg_result = as.data.frame(
        "None of the selected genes are in the background genes!"
      )
    ))
  }

  # Make sure the background set includes the query set
  query_set <- unique(c(query_set, sub_query))

  sql_query <- paste(
    "select distinct gene,pathwayID from pathway where gene IN ('", 
    paste(query_set, collapse = "', '"), "')", sep = ""
  )    
  # Restrict to pathways with genes matching sub_query (Improves speed)
  sql_query <- paste0(sql_query, " AND pathwayID IN ('",
    paste(pathway_table$pathway_id, collapse = "', '"),"')" 
  )

  if(go != "All") {
    sql_query <- paste0(sql_query, " AND category ='", go, "'")
  }

  results <- DBI::dbGetQuery(pathway, sql_query)

  if(dim(results)[1] == 0) {
    return(list(
      bg_result = as.data.frame("No matching species or gene ID file!")
    ))
  }
  bg_result <- table(results$pathwayID)
  bg_result <- as.data.frame(bg_result)
  colnames(bg_result) <- c("pathway_id", "overlap_bg")

  pathway_table_bg <- merge(
    pathway_table,
    bg_result,
    by = "pathway_id",
    all.x = TRUE
  )

  DBI::dbDisconnect(pathway)
  
  return(pathway_table_bg)
}

#' Database choices for the converted IDs
#' 
#' This function provides a list of the portions of the
#' database to query that contain genes from the converted
#' IDs. 
#' 
#' @param converted Conversion information from the original IDs returned
#'   from the convert_id() function
#' @param converted_data Data matrix with the converted ensembl_IDs
#' @param select_org Input for what organism the IDs are pertaining to
#' @param gmt_file Inputed gene set pathway file for NEW species
#' @param idep_data Data built in to idep
#' 
#' @export
#' @return A character vector of names of the section of data to perform
#'   pathway analysis on
gmt_category <- function(
  converted,
  converted_data,
  select_org,
  gmt_file,
  idep_data
) {
	if (select_org == "NEW" && !is.null(gmt_file)) {
    return(list(custom_gene_set ="Custom"))
  }
		
	id_not_recognized = as.data.frame("ID not recognized!")

	if(is.null(converted)) {
    return(id_not_recognized)
  }

	query_set <- rownames(converted_data)

	if(length(query_set) == 0){
    return(id_not_recognized )
  }

	ix = grep(converted$species[1,1], idep_data$gmt_files)

	if (length(ix) == 0) {
    return(id_not_recognized)
  }
	
	# If selected species is not the default "bestMatch", use that species directly
	if(select_org != idep_data$species_choice[[1]]) {  
		ix = grep(
      find_species_by_id(select_org, idep_data$org_info)[1,1],
      idep_data$gmt_files
    )
		if (length(ix) == 0) {
      return(id_not_recognized)
    }
	}

	pathway <- DBI::dbConnect(
    drv = RSQLite::dbDriver("SQLite"),
    dbname = idep_data$gmt_files[ix],
    flags= RSQLite::SQLITE_RO
  )

	# Generate a list of geneset categories such as "GOBP", "KEGG" from file
	gene_set_category <-  DBI::dbGetQuery(pathway, "select distinct * from categories") 
	gene_set_category  <- sort(gene_set_category[, 1])
	category_choices <- setNames(as.list(gene_set_category ), gene_set_category )
	
	# Set order of popular elements
  top_choices <- c("GOBP", "GOCC", "GOMF", "KEGG")
  other_choices <- names(category_choices)[
    !(names(category_choices) %in% top_choices)
  ]
  category_choices <- category_choices[c(top_choices, other_choices)]

	# Change names to the full description for display
	names(category_choices)[match("GOBP", category_choices)] <- "GO Biological Process"
	names(category_choices)[match("GOCC", category_choices)] <- "GO Cellular Component"
	names(category_choices)[match("GOMF", category_choices)] <- "GO Molecular Function"
	category_choices <- append(setNames("All", "All available gene sets"), category_choices)
	
	DBI::dbDisconnect(pathway)

	return(category_choices)
}

#' Read pathway gene sets
#' 
#' Use the IDs from the converted database to find the
#' gene sets for all the pathways in the database. Returns
#' a list with each entry a vector of IDs corresponding to
#' a description of a pathway in the database.
#' 
#' @param converted Return value from the \code{convert_id} function,
#'  contains information about the gene IDs for the matched species
#' @param all_gene_names Data frame with all gene names
#' @param go Portion of the database to use for the pathway analysis
#' @param select_org The user selected organism for the expression data,
#'  default is "BestMatch"
#' @param idep_data Read data files from the database
#' @param my_range Vector of the (min_set_size, max_set_size)
#' 
#' @export
#' @return A list with each entry a list of gene IDs that correspond to
#'  a pathway.
read_gene_sets <- function(
  converted,
  all_gene_names,
  go,
  select_org,
  idep_data,
  my_range
) {
	id_not_recognized <- as.data.frame("ID not recognized!")
	if(is.null(converted)) {
    return(id_not_recognized)
  }
  if(!is.null(all_gene_names[, 2])) {
    query_set <- all_gene_names[, 2]
  }
	if(is.null(query_set) || length(query_set) == 0) {
    return(id_not_recognized)
  }
	ix <- grep(converted$species[1, 1], idep_data$gmt_files)
	if(length(ix) == 0) {
    return(id_not_recognized)
  }
	
	# If selected species is not the default "bestMatch", use that species directly
	if(select_org != "BestMatch") {  
		ix <- grep(find_species_by_id(select_org)[1, 1], idep_data$gmt_files)
		if(length(ix) == 0) {
      return(id_not_recognized)
    }
	}
	pathway <- DBI::dbConnect(
    drv = RSQLite::dbDriver("SQLite"),
    dbname = idep_data$gmt_files[ix],
    flags= RSQLite::SQLITE_RO
  )
	
	if(is.null(go)) {
    go <- "GOBP"
  }  

	# Get Gene sets
	sql_query <- paste(
    "select distinct gene, pathwayID from pathway where gene IN ('",
    paste(query_set, collapse = "', '"),
    "')",
    sep = ""
  )
	
	if(go != "All") {
    sql_query <- paste0(sql_query, " AND category ='", go,"'")
  }
	result <- DBI::dbGetQuery(pathway, sql_query)
	if(dim(result)[1] == 0) {
    return(list(x = as.data.frame("No matching species or gene ID file!" )))
  }
	# List pathways and frequency of genes
	pathway_ids <- aggregate(
    result$pathwayID,
    by = list(unique.values = result$pathwayID),
    FUN = length
  )
	pathway_ids <- pathway_ids[which(pathway_ids[, 2] >= my_range[1]), ]
	pathway_ids <- pathway_ids[which(pathway_ids[, 2] <= my_range[2]), ]
	if(dim(pathway_ids)[1] == 0) {
    gene_sets = NULL
  }
	
	# Convert pathways into lists like those generated by readGMT
	gene_sets <- lapply(pathway_ids[, 1], function(x) result[which(result$pathwayID == x), 1])
	pathway_info <- DBI::dbGetQuery(
    pathway,
    paste(
      "select distinct id,Description from pathwayInfo where id IN ('", 
		   paste(pathway_ids[, 1], collapse = "', '"),
       "') ",
      sep = ""
    )
  )
	ix <- match(pathway_ids[, 1], pathway_info[, 1])
	names(gene_sets) <- pathway_info[ix, 2]
	DBI::dbDisconnect(pathway)
	return(gene_sets)
}

#' Convert IDs from ensembl to entrez
#' 
#' Convert an ID qeury for a species from the ensembl
#' ID type to entrez type.
#' 
#' @param query Vector of IDs to convert
#' @param species Species to map the IDs for
#' @param org_info org_info file from the \code{get_idep_data}
#'  function
#' 
#' @export
#' @return The qeuried genes with converted IDs.
convert_ensembl_to_entrez <- function(
  query,
  species,
  org_info
) { 
	query_set <- clean_gene_set(unlist(strsplit(toupper(names(query)), '\t| |\n|\\, ')))
  # Note uses species Identifying
	species_id <- org_info$id[which(org_info$ensembl_dataset == species)]  
	# idType 6 for entrez gene ID
  convert <- connect_convert_db()

  id_type_entrez <- DBI::dbGetQuery(
    convert,
    paste(
      "select distinct * from idIndex where idType = 'entrezgene_id'" 
    )
  )
  if(dim(id_type_entrez)[1] != 1) {
    cat("Warning! entrezgene ID not found!")
  }
  id_type_entrez <- as.numeric(id_type_entrez[1, 1])

	result <- DBI::dbGetQuery(
    convert,
		paste(
      "select  id,ens,species from mapping where ens IN ('",
      paste(query_set, collapse = "', '"),
			"') AND  idType ='",
      id_type_entrez,
      "'",
      sep = ""
    )
  )
	DBI::dbDisconnect(convert)
	if(dim(result)[1] == 0) {
    return(NULL)
  }
  result <- subset(result, species == species_id, select = -species)

	ix = match(result$ens, names(query))
  
	tem <- query[ix]
  names(tem) <- result$id
  return(tem)
}

#' Find pathway IDs for a KEGG description
#' 
#' From a pathway description find the KEGG ID.
#' 
#' @param pathway_description Description for a pathway in the
#'  iDEP database
#' @param Species Species that the pathway is for
#' @param GO Portion of the database to query ("KEGG")
#' @param select_org The organism that the gene data is for
#' @param gmt_files GMT files from the iDEP database
#' @param org_info Organism information files from the iDEP
#'  database
#' @param idep_data Background idep_data
#' 
#' @export
#' @return Return the KEGG ID for the pathway.
kegg_pathway_id <- function (
  pathway_description,
  Species,
  GO,
  select_org,
  gmt_files,
  org_info,
  idep_data
) {
	ix <- grep(Species, gmt_files)

	if(length(ix) == 0) {
    return(NULL)
  }
	
	# If selected species is not the default "bestMatch", use that species directly
	if(select_org != "BestMatch") {  
		ix <- grep(find_species_by_id(select_org)[1, 1], gmt_files)
		if(length(ix) == 0) {
      return(NULL)
    }
		total_genes <- org_info[which(org_info$id == as.numeric(select_org)), 7]
	}
	pathway <- DBI::dbConnect(
    drv = RSQLite::dbDriver("SQLite"),
    dbname = idep_data$gmt_files[ix],
    flags = RSQLite::SQLITE_RO
  )
	
	# change Parkinson's disease to Parkinson\'s disease    otherwise SQL 
	pathway_description <- gsub("\'", "\'\'", pathway_description)
							
	pathway_info <- DBI::dbGetQuery(
    pathway,
    paste(
      " select * from pathwayInfo where description =  '", 
			pathway_description,
      "' AND name LIKE '",
      GO,
      "%'",
      sep = ""
    )
  )
	DBI::dbDisconnect(pathway)
	
  if(dim(pathway_info)[1] != 1) {
    return(NULL)
  }
	tem <- gsub(".*:", "", pathway_info[1, 2])  
	
  return(gsub("_.*", "", tem))
}