#' fct_07_genome.R This file holds all of the main data analysis functions
#' associated with eighth tab of the iDEP website.
#'
#'
#' @section fct_07_genome.R functions:
#'
#'
#' @name fct_07_genome.R
NULL

#' Plotly of chromosome position
#' 
#' Calculate a plotly that shows the chromosome position of significant genes
#' and significantly enriched regions.
#' 
#' @param limma Return from \code{limma_value} function
#' @param select_contrast DEG contrast to examine
#' @param all_gene_info Gene information return from \code{get_gene_info}
#' @param ignore_non_coding When TRUE only use protein coding genes
#' @param limma_p_val_viz Adjusted p-value to use for significant genes
#' @param limma_fc_viz Minimum fold-change value to filter with
#' @param label_gene_symbol Paste the gene symbol label on the plot
#' @param ma_window_size Moving average window size for a chromosome 
#'   (1, 2, 4, 6, 8, 10, 15, 20)
#' @param ma_window_steps Number of moving average window steps (1, 2, 3, 4)
#' @param ch_region_p_val P-value to use for finding significant chromosome
#'  region enrichment
#'
#' @export
#' @return Plotly visualization of chromosomes and significantly enriched genes
chromosome_plotly <- function(
  limma,
  select_contrast,
  all_gene_info,
  ignore_non_coding,
  limma_p_val_viz,
  limma_fc_viz,
  label_gene_symbol,
  ma_window_size,
  ma_window_steps,
  ch_region_p_val
) {
  # Default plot
  fake <- data.frame(a = 1:3, b = 1:3)
  p <- ggplot2::ggplot(fake, ggplot2::aes(x = a, y = b)) +
    ggplot2::geom_blank() +
    ggplot2::ggtitle("No genes with position info.") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )

  if(length(limma$comparisons) == 1) {
    top_1 <- limma$top_genes[[1]]  
  } else {
	top <- limma$top_genes
	ix <- match(select_contrast, names(top))
	if(is.na(ix)) {
      return(plotly::ggplotly(p))
    }
	top_1 <- top[[ix]] 
  }
  if(dim(top_1)[1] == 0) {
    return(plotly::ggplotly(p))
  }

  # Species in STRING-db do not have chr. location data
  if(sum(!is.na(all_gene_info$start_position)) < 5) { 
   return(p) 
  }

  colnames(top_1) <- c("Fold", "FDR")
  
  x <- merge(
    top_1,
    all_gene_info,
    by.x = "row.names",
    by.y = "ensembl_gene_id"
  )
  
  colnames(x)[which(colnames(x) == "Row.names")] <- "ensembl_gene_id"

  # only coding genes? 
  if(ignore_non_coding) {
    x <- subset(x, gene_biotype == "protein_coding")
  }
  
  # If no chromosomes found. For example if user do not convert gene IDs.
  if(dim(x)[1] > 5) {
    x <- x[order(x$chromosome_name, x$start_position), ]
    x$ensembl_gene_id <- as.character(x$ensembl_gene_id)
    # If symbol is missing use Ensembl id
    x$symbol <- as.character(x$symbol)  
    ix <- which(is.na(x$symbol))
    ix2 <- which(nchar(as.character(x$symbol)) <= 2)
    ix3 <- which(duplicated(x$symbol))
    ix <- unique(c(ix, ix2, ix3))
    x$symbol[ix] <- x$ensembl_gene_id[ix] 
    
    x <- x[!is.na(x$chromosome_name), ]
    x <- x[!is.na(x$start_position), ]
            
    tem <- sort(table(x$chromosome_name), decreasing = T)
    # ch with less than 100 genes are excluded
    ch <- names(tem[tem >= 1])  
    if(length(ch) > 50) {
      # At most 50 ch
      ch <- ch[1:50]
    }
    # ch. name less than 10 characters
    ch <- ch[nchar(ch) <= 12]
    ch <- ch[order(as.numeric(ch) ) ]
    tem <- ch
    # The numbers are continous from 1 to length(ch)
    ch <- 1:(length(ch))
    # The names are real chr. names
    names(ch) <- tem
    x <- x[which(x$chromosome_name %in% names(ch)), ]
    x <- droplevels(x)
    # Numeric encoding
    x$chNum <- 1 
    x$chNum <- ch[x$chromosome_name]
      
    # Use max position as chr. length before filtering
    ch_length_table <- aggregate(start_position ~ chromosome_name, data = x, max)
    # Add chr. numer 
    ch_length_table$chNum <- ch[ch_length_table$chromosome_name]
    ch_length_table <- ch_length_table[!is.na(ch_length_table$chNum), ]
    ch_length_table <- ch_length_table[order(ch_length_table$chNum), c(3, 2)]
    ch_length_table <- ch_length_table[order(ch_length_table$chNum), ]
    ch_length_table$start_position <- ch_length_table$start_position / 1e6

    # Only keep significant genes
    ix <- which(
      (x$FDR< as.numeric(limma_p_val_viz)) &
      (abs(x$Fold) > log2(as.numeric(limma_fc_viz)))
    )

    if(length(ix) > 5) { 
      # Remove non-significant / not selected genes
      # Keep a copy
      x0 <- x
      x <- x[ix, ]
      # Prepare coordinates
      # Mbp
      x$start_position <- x$start_position / 1000000
      # Distance between chs.
      chD <- 30
      # Max log2 fold 
      fold_cutoff <- 4
        
      # log2 fold within -5 to 5
      x$Fold[which(x$Fold > fold_cutoff)] <- fold_cutoff
      x$Fold[which(x$Fold < -1 * fold_cutoff)] <- -1 * fold_cutoff 
      x$Fold <- 4 * x$Fold
    
      x$y <- x$chNum * chD + x$Fold
      ch_total <- dim(ch_length_table)[1] 
      x$R <- as.factor(sign(x$Fold))
        
      colnames(x)[which(colnames(x) == "start_position")] <- "x"
      x$R <- as.character(x$R)
      x$R[x$R == "-1"] <- "Down"
      x$R[x$R == "1"] <- "Up"
      x$R <- as.factor(x$R)
        
      # Don't define x and y, so that we could plot use two datasets
      p <- ggplot2::ggplot() +
        ggplot2::geom_point(
          data = x,
          ggplot2::aes(
            x = x,
            y = y,
            color = R,
            text = paste0(
              "Symbol: ",
              symbol,
              "\nRegulation: ",
              R,
              "\nChr Pos: ",
              x
            )
          ),
          shape = 20,
          size = 0.2
        )

      if(label_gene_symbol) {
        p <- p + ggplot2::geom_text(
          data = x,
          ggplot2::aes(x = x, y = y, label = symbol),
          check_overlap = FALSE,
          angle = 45,
          size = 2,
          vjust = 0,
          nudge_y = 4 
        )
      }
      # Label y with ch names
      p <- p + ggplot2::scale_y_continuous(
        labels = paste(
          "chr",
          names(ch[ch_length_table$chNum]), sep = ""
        ), 
        breaks = chD * (1:ch_total), 
        limits = c(0, chD * (ch_total + 1) + 5)
      )
        
      # Draw horizontal lines for each ch.
      for(i in 1:dim(ch_length_table)[1]) {
        p <- p + ggplot2::annotate(
          "segment",
          x = 0,
          xend = ch_length_table$start_position[i],
          y = ch_length_table$chNum[i] * chD,
          yend = ch_length_table$chNum[i] * chD
        )
      }
      # Change legend  http://ggplot2.tidyverse.org/reference/scale_manual.html
      # Customize legend text
      p <- p + ggplot2::scale_color_manual(
        name = "",   
        values = c("red", "blue"),
        breaks = c("Up", "Down"),
        labels = c("Up", "Dn")
      ) 
      p <- p + ggplot2::xlab("Position on chrs. (Mbp)") +
        ggplot2::theme(axis.title.y = ggplot2::element_blank())      
      p <- p + ggplot2::theme(legend.position = "none")

             
      x0 <- x0[x0$chromosome_name %in% unique(x$chromosome_name), ]
      # Numeric encoding
      x0$chNum <- 1
      x0$chNum <- ch[ x0$chromosome_name ]
      # Mbp
      x0$start_position <- x0$start_position / 1e6             

      window_size <- as.numeric(ma_window_size) 
      # Step size is then windowSize / steps         
      steps <- as.numeric(ma_window_steps)       

      # Centering 
      x0$Fold <-  x0$Fold - mean(x0$Fold)             

            
                             
      for(i in 0:(steps-1)) {
        # Step size is windowSize/steps   
        # If windowSize=10 and steps = 2; then step size is 5Mb
        # 1.3 becomes 5, 11.2 -> 15 for step 1
        # 1.3 -> -5
        x0$x <- (
          floor((x0$start_position - i * window_size / steps) / window_size)  
          + 0.5 + i / steps
        ) * window_size

        moving_average_start <- x0 |>
          dplyr::select(chNum, x, Fold) |>
          # Beginning bin can be negative for first bin in the 2nd step
          dplyr::filter( x >= 0) |>   
          dplyr::group_by(chNum, x) |>
          dplyr::summarize(
            ma = mean(Fold),
            n = dplyr::n(),
            pval = ifelse(
              dplyr::n() >= 3 && sd(Fold) > 0, stats::t.test(Fold)$p.value, 0
            )
          ) |>
          # na when only 1 data point?
          dplyr::filter(!is.na(pval))

        if(i == 0) {
          moving_average <- moving_average_start
        } else {
          moving_average <- rbind(moving_average, moving_average_start)  
        }        
      }      
      # Translate fold to y coordinates
      moving_average <- moving_average |>
        dplyr::filter(n >= 3) |>
        dplyr::mutate(pval = stats::p.adjust(pval, method = "fdr")) |>
        dplyr::filter(pval < as.numeric(ch_region_p_val)) |>
        dplyr::mutate(y = ifelse(ma > 0, 1, -1)) |>
        dplyr::mutate(y = chNum * chD + 3 * y) |>
        dplyr::mutate(ma = ifelse(ma > 0, 1, -1)) |>
        dplyr::mutate(ma = as.factor(ma))

      moving_average$ma <- as.character(moving_average$ma)
      moving_average$ma[moving_average$ma == "-1"] <- "Down"
      moving_average$ma[moving_average$ma == "1"] <- "Up"
      moving_average$ma <- as.factor(moving_average$ma)

      # Significant regions are marked as horizontal error bars 
      if(dim(moving_average)[1] > 0) {
        p <- p + ggplot2::geom_errorbarh(
          data = moving_average,
          ggplot2::aes(
            x = x, 
            y = y, 
            xmin = x - window_size / 2, 
            xmax = x + window_size / 2,
            colour = ma,
            text = paste0(
              paste0(
              "Window: ",
              x - window_size / 2,
              " - ",
              x + window_size / 2,
              "\nRegulation: ",
              ma,
              "\nChr Pos: ",
              x
            )
            )
          ), 
          size = 2, 
          height = 15
        )

        # Label significant regions
        sig_ch <- sort(table(moving_average$chNum), decreasing = TRUE)
        sig_ch <- names(ch)[ as.numeric(names(sig_ch))]
        if(length(sig_ch) <= 5) {
          # More than 5 just show 5
          sig_ch <- paste0("chr", sig_ch, collapse = ", ")
        } else {
          sig_ch <- sig_ch[1:5]
          sig_ch <- paste0("chr", sig_ch, collapse = ", ")                  
          sig_ch <- paste0(sig_ch, ", ...")
        }

        sig_ch <- paste(
          dim(moving_average)[1], 
          " enriched regions \n(",
          round(
            sum(ch_length_table$start_position) /
            window_size * steps * as.numeric(ch_region_p_val),
            2
          ),
          " expected)  detected on:\n ",
          sig_ch
        )
                 
        p <- p + ggplot2::annotate(
          geom = "text", 
          x = max(x$x) * 0.70,
          y = max(x$y) * 0.90,
          label = sig_ch
        )
      }
    } 		
  }
  plotly::ggplotly(p, tooltip = "text") 
}

#' PREDA PACKAGE ERRORR #########################
get_genome_plot <- function(
  genome_plot_data,
  regions_p_val_cutoff,
  statistic_cutoff
) {
  if(class(genome_plot_data) != "list") {
    plot.new()
    text(0.2, 1, "No significant regions found!")
  } else {
    ge_analysis_results <- genome_plot_data$mainResult
    if(is.null(ge_analysis_results)) {
      return(NULL)
    }	
	  genomic_regions_up <- PREDA::PREDAResults2GenomicRegions(
      ge_analysis_results,
      qval.threshold = regions_p_val_cutoff,
      smoothStatistic.tail = "upper",
      smoothStatistic.threshold = -1 * statistic_cutoff
    )
	  genomic_regions_down <- PREDA::PREDAResults2GenomicRegions(
      ge_analysis_results,
      qval.threshold = regions_p_val_cutoff,
      smoothStatistic.tail = "lower",
      smoothStatistic.threshold = -1 * statistic_cutoff
    )
	 
	  checkplot <- PREDA::genomePlot(
      ge_analysis_results,
      genomicRegions = c(genomic_regions_up, genomic_regions_down),
      grouping = c(1, 1),
      scale.positions = "Mb",
      region.colors = c("red", "blue")
    )
	  legend(
      x = genome_plot_data$legend.x,
      y = genome_plot_data$legend.y,
      legend = c("UP", "DOWN"),
      fill = c("red", "blue")
    )
	}
}

#' PREDA PACKAGE ERROR ###########################
get_genome_plot_data <- function(
  genome_plot_data_pre,
  all_gene_info,
  select_contrast,
  limma,
  regions_p_val_cutoff,
  statistic_cutoff
) {
	if(length(limma$comparisons) == 1) {
    top_1 <- limma$top_genes[[1]]  
	} else {
	  top <- limma$top_genes
	  ix <- match(select_contrast, names(top))
	  if(is.na(ix)) {
      return(NULL)
    }
	  top_1 <- top[[ix]] 
	}
	
  if(dim(top_1)[1] == 0) {
    return(NULL)
  }

  # Species in STRING-db do not have chr. location data
  if(sum(!is.na(all_gene_info$start_position)) < 5) { 
   return(NULL) 
  }
	colnames(top_1) <- c("Fold", "FDR")
	  
	x <- merge(
    top_1,
    all_gene_info,
    by.x = "row.names",
    by.y = "ensembl_gene_id"
  )
	 
	x <- x[order(x$chromosome_name, x$start_position), ]
	tem <- sort(table(x$chromosome_name), decreasing = T)
  # Chromosomes with less than 100 genes are excluded 
	chromosomes <- names(tem[tem > 100 ])
	if(length(chromosomes) > 50) {
    # At most 50 chromosomes
    chromosomes <- chromosomes[1:50]
  }

	chromosomes <- chromosomes[order(as.numeric(chromosomes))]
	chromosomes_numbers <- as.numeric(chromosomes)
	# Convert chr.x to numbers
	j <- max(chromosomes_numbers, na.rm = T) 
	for(i in 1:length(chromosomes)) {
	  if(is.na(chromosomes_numbers[i])) {
      chromosomes_numbers[i] <- j + 1
      j <- j + 1
    }
	}
	  
	x <- x[which(x$chromosome_name %in% chromosomes), ]
	x <- droplevels(x)
	 
	# Find the number coding for chromosome 
	get_chr_number <- function(
    chr_name,
    chromosomes_numbers,
    chromosomes
  ){
	  return(
      chromosomes_numbers[which(chromosomes == chr_name)]
    )
	}
  # Numeric coding
	x$chrNum <- 1
	x$chrNum <- unlist(lapply(
    x$chromosome_name,
    function(x) get_chr_number(
      chr_name = x,
      chromosomes_numbers = chromosomes_numbers,
      chromosomes = chromosomes
    )
  ))
	 
	x$Row.names <- as.character(x$Row.names)
	fold <- x$Fold
  names(fold) <- x$Row.names
	 
	
	x <- x[!duplicated(x$Row.names), ]

	ge_analysis_results <- genome_plot_data_pre 
	
	genomic_regions_up <- PREDA::PREDAResults2GenomicRegions(
    ge_analysis_results,
    qval.threshold = regions_p_val_cutoff,
    smoothStatistic.tail = "upper",
    smoothStatistic.threshold = statistic_cutoff
  )
	genomic_regions_down <- PREDA::PREDAResults2GenomicRegions(
    ge_analysis_results,
    qval.threshold = regions_p_val_cutoff,
    smoothStatistic.tail = "lower",
    smoothStatistic.threshold = -1 * statistic_cutoff
  )

	if(is.null(genomic_regions_up$Test) && 
     is.null(genomic_regions_down$Test)) {
    # No significant regions
    return(-1)
  }
	
	regions <- 0 
	if(!is.null(genomic_regions_up$Test)) { 
	  dataframe_up_regions<- PREDA::GenomicRegions2dataframe(
      genomic_regions_up[[1]]
    )
	  dataframe_up_regions$Regulation <- "Up"
	  regions <- dataframe_up_regions
	}

	if(!is.null(genomic_regions_down$Test)) { 
	  dataframe_down_regions <- PREDA::GenomicRegions2dataframe(
      genomic_regions_down[[1]]
    )
	  dataframe_down_regions$Regulation <- "Down"
	  if(class(regions) != "data.frame") {
      # if UP regions is NULL
      regions <- dataframe_down_regions 
    } else {
      regions <- rbind(regions, dataframe_down_regions)
    } 
	}

	if(class(regions) != "data.frame") {
    return(-1)
  }

	Regions <- regions[, c(4, 1:3)] 
	Regions <- Regions[which(Regions$end != Regions$start), ]
	Regions$size <- round((Regions$end - Regions$start) / 1000000, 3)
  # Convert from chr. number to names
	Regions$chr <- chromosomes[Regions$chr]

	# Find the gene indices in the up or down regulated regions
	regulated_genes <- function(
    i,
    Regions,
    x
  ) {
    # If the start position is within the region
	  ix <- which(
      Regions$chr == x$chromosome_name[i] &
		  Regions$start < x$start_position[i] &
      (Regions$end > x$start_position[i])
    )
	  if(length(ix) == 0 | length(ix) > 1) {
      return(NA)
    } else {
      return(ix)
    }
	}

	region_id <- unlist(lapply(
    1:dim(x)[1],
    function(w) regulated_genes(
      i = w,
      Regions = Regions,
      x = x
    )
  ))
	x1 <- x[which(!is.na(region_id)), ]
	region_id <- region_id[!is.na(region_id)]
	x1 <- cbind(
    region_id,
    Regions[region_id, ],
    x1[, c(
      "symbol", "Row.names", "Fold", "FDR", "band","start_position"
    )]
  )
	x1 <- x1[order(x1$region_id, x1$start_position), ]	
  colnames(x1)[8] <- "Ensembl"
  colnames(x1)[4] <- "Region.Start"
  colnames(x1)[5] <- "Region.End"
  colnames(x1)[12] <- "Gene.Start"

	# Number of genes
	tem <- table(x1$region_id)
	Regions$Ngenes <- 0
	Regions$Ngenes[as.integer(names(tem))] <- tem

	# Cytoband per region
	tem <- unique(x1[, c("region_id", "band")])
	tem$band <- gsub("\\..*", "", tem$band)
	tem <- unique(tem)
	Regions$band <- ""
	Regions$ID <- 1:dim(Regions)[1]
	for(i in 1:dim(Regions)[1]) {
    Regions$band[i] <- paste(
      tem[which(tem[, 1] == Regions$ID[i]), 2],
      collapse = ";"
    )
  }
	  
	# Genes
	Regions$Genes <- ""
	for(i in 1:dim(Regions)[1]) {
    Regions$Genes[i] <- paste(
      x1$symbol[which(x1[, 1] == Regions$ID[i])],
      collapse = " "
    )
  }
	return(list(
    mainResult = ge_analysis_results,
    Regions = Regions,
    Genes = x1,
    legend.x = max(x$start_position) * .6,
    legend.y = max(chromosomes_numbers) - 3
  ))
}

#' PREDA PACKAGE ERROR ########################
get_genome_plot_data_pre <- function(
  select_contrast,
  limma,
  all_gene_info
) {
	if(length(limma$comparisons) == 1) {
    top_1 <- limma$top_genes[[1]]  
	} else {
	  top <- limma$top_genes
	  ix <- match(select_contrast, names(top))
	  if(is.na(ix)) {
      return(NULL)
    }
	  top_1 <- top[[ix]] 
	}
	  
  if(dim(top_1)[1] == 0) {
    return(NULL)
  }
	colnames(top_1) <- c("Fold","FDR")
	
	
	x <- merge(
    top_1,
    all_gene_info,
    by.x = "row.names",
    by.y = "ensembl_gene_id"
  )
	
	
	# PREDA
	info_file <- system.file(
    "sampledata",
    "GeneExpression",
    "sampleinfoGE_PREDA.txt",
    package = "PREDAsampledata"
  )
	sample_info <-read.table(
    info_file,
    sep = "\t",
    header = TRUE
  )
	head(sample_info)

	data("ExpressionSetRCC", package = "PREDAsampledata")

	ge_statistics_for_preda <- PREDA::statisticsForPREDAfromEset(
    ExpressionSetRCC,
    statisticType = "tstatistic",
    referenceGroupLabel = "normal",
    classVector = sample_info[, "Class"]
  )
	PREDA::analysesNames(ge_statistics_for_preda)

  # Needs hgu133plus2.db package
	ge_genomic_annotations <- PREDA::eset2GenomicAnnotations(
    ExpressionSetRCC,
    retain.chrs = 1:22
  )
	ge_genomic_annotations_for_preda <- 
    PREDA::GenomicAnnotations2GenomicAnnotationsForPREDA(
      ge_genomic_annotations, reference_position_type = "median"
    )
	ge_data_for_preda <- PREDA::MergeStatisticAnnotations2DataForPREDA(
    ge_statistics_for_preda,
    ge_genomic_annotations_for_preda,
    sortAndCleanNA = TRUE
  )

	x <- x[order(x$chromosome_name, x$start_position), ]
	tem <- sort(table(x$chromosome_name), decreasing = T)
  # Chromosomes with less than 100 genes are excluded
	chromosomes <- names(tem[tem > 100 ])
	if(length(chromosomes) > 50) {
    # At most 50 chromosomes
    chromosomes <- chromosomes[1:50]
  }

	chromosomes <- chromosomes[order(as.numeric(chromosomes))]
	chromosomes_numbers <- as.numeric(chromosomes)
	# Convert chr.x to numbers
	j <- max(chromosomes_numbers, na.rm = T)
	for(i in 1:length(chromosomes)) {
	  if(is.na(chromosomes_numbers[i])) {
      chromosomes_numbers[i] <- j + 1
      j <- j + 1
    }
	}
	  
	x <- x[which(x$chromosome_name %in% chromosomes), ]
	x <- droplevels(x)
	 
	# Find the number coding for chromosome 
	get_chr_number <- function(
    chr_name,
    chromosomes,
    chromosomes_numbers
  ){
	  return(chromosomes_numbers[which(chromosomes == chr_name)])
	}
  # Numeric coding
	x$chrNum <- 1
	x$chrNum <- unlist(lapply(
    x$chromosome_name,
    function(x) get_chr_number(
      chr_name = x,
      chromosomes = chromosomes,
      chromosomes_numbers = chromosomes_numbers
    )
  ))
	 
	x$Row.names <- as.character(x$Row.names)
	fold <- x$Fold
  names(fold) <- x$Row.names
	
  x <- x[!duplicated(x$Row.names), ]
	
	my_data <- ge_data_for_preda
	  
	my_data@position <- as.integer(x$start_position + x$genomeSpan / 2)
	my_data@ids <- as.character(x$Row.names)
	my_data@chr <- as.integer(x$chrNum)
	my_data@start <- x$start_position
	my_data@end <- x$start_position + x$genomeSpan
	my_data@strand <- rep(1, dim(x)[1])
	my_data@chromosomesLabels <- chromosomes
	my_data@chromosomesNumbers <- as.integer(chromosomes_numbers)
	my_data@statistic <- as.matrix(x[, 2, drop = FALSE])
	my_data@analysesNames <- "Test"
	my_data@testedTail <- "both"
	my_data@optionalAnnotations <- as.matrix(x[, 1:2, drop = FALSE])
	my_data@optionalAnnotationsHeaders <- c("a", "b")

	set.seed(2)
	ge_analysis_results <- PREDA::PREDA_main(
    my_data,
    nperms = 1000
  )

	return(ge_analysis_results)
}