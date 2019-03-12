list.of.packages <- c("shiny","RSQLite","gplots","ggplot2",
    "e1071", # computing kurtosis
    "reshape2", # for melt correlation matrix in heatmap
    "DT", # for renderDataTable
    "plotly", # for interactive heatmap
    # bioconductor packages
    "limma", # Differential expression
    "DESeq2", # count data analysis
    "edgeR", # count data D.E.
    "gage", # pathway analysis
    "PGSEA", # pathway
    "fgsea", # fast GSEA
    "ReactomePA", # pathway analysis
    "pathview",
    "PREDA",  # showing expression on genome
    "PREDAsampledata","sfsmisc","lokern","multtest")
lapply(list.of.packages, require, character.only = TRUE)