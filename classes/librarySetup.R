# dplyr complains this required libraries: libudunits2-dev, libmariadb-client-lgpl-dev
# install.packages("plotly", repos="http://cran.rstudio.com/", dependencies=TRUE)

biocLibs = c( "limma", "DESeq2","edgeR","gage", "PGSEA", "fgsea", "ReactomePA", 
"pathview","PREDA","PREDAsampledata","sfsmisc","lokern","multtest","hgu133plus2.db")

list.of.packages <- c(
  "shiny", "shinyAce", "shinyBS", "plotly",
  "RSQLite", "gplots", 
  "ggplot2", "dplyr", #"tidyverse",
  "plotly",
  "e1071", "reshape2", "DT",
  "data.table", "Rcpp",
  #"RPostgreSQL"
)

list.of.bio.packages  <- c(
  "limma", "DESeq2", "edgeR", "gage", "PGSEA", "fgsea", "ReactomePA", "pathview", "PREDA",
  "PREDAsampledata", "sfsmisc", "lokern", "multtest", "hgu133plus2.db"
  # ,"org.Ag.eg.db","org.At.tair.db","org.Bt.eg.db","org.Ce.eg.db","org.Cf.eg.db",
  # "org.Dm.eg.db","org.Dr.eg.db","org.EcK12.eg.db","org.EcSakai.eg.db","org.Gg.eg.db",
  # "org.Hs.eg.db","org.Hs.ipi.db","org.Mm.eg.db","org.Mmu.eg.db","org.Pf.plasmo.db",
  # "org.Pt.eg.db","org.Rn.eg.db","org.Sc.sgd.db","org.Sco.eg.db","org.Ss.eg.db",
  # "org.Tgondii.eg.db","org.Xl.eg.db"
)

#Install Require packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="http://cran.rstudio.com/", dependencies=TRUE)

new.bio.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.bio.packages)){
  source("https://bioconductor.org/biocLite.R")
  biocLite(new.bio.packages, suppressUpdates = T)
}

#Load Packages
lapply(list.of.packages, require, character.only = TRUE)
lapply(list.of.bio.packages, require, character.only = TRUE)