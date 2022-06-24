
# Install packages for iDEP and associated apps
# revised 6/23/2022. 
# Install R packages, Bioconductor packages, then install ottoPlots and idepGolem. 
# idepGolem has dependencies. 

###############################################################################
# List of packages on CRAN
###############################################################################
list.of.packages <- c(
  "shiny", "shinyAce", "shinyBS", "plotly",
  "RSQLite", "gplots", 
  "ggplot2", "dplyr", #"tidyverse",
  "plotly",
  "e1071", "reshape2", "DT",
  "data.table", "Rcpp", "flashClust","statmod","biclust","igraph","Rtsne",
  "visNetwork", "feather","shinyjs","reactable", "remotes"
)

###############################################################################
# List of packages on Bioconductor
###############################################################################
list.of.bio.packages  <- c(
  "getDEE2", "limma", "DESeq2", "edgeR", "gage", "fgsea", "ReactomePA", "pathview", "PREDA",
  "impute", "runibic","QUBIC","rhdf5", "STRINGdb",
  "PREDAsampledata", "sfsmisc", "lokern", "multtest", "hgu133plus2.db",
  "BiocParallel", "ComplexHeatmap", "InteractiveComplexHeatmap", "SummarizedExperiment",
  "impute", "preprocessCore", # required by WGCNA
  "org.Ag.eg.db","org.At.tair.db","org.Bt.eg.db","org.Ce.eg.db","org.Cf.eg.db",
  "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db","org.EcSakai.eg.db","org.Gg.eg.db",
  "org.Hs.eg.db","org.Mm.eg.db","org.Mmu.eg.db","org.Pf.plasmo.db",
  "org.Pt.eg.db","org.Rn.eg.db","org.Sc.sgd.db","org.Ss.eg.db","org.Xl.eg.db"
)

###############################################################################
# Clean up existing packagges
###############################################################################
if(1) { # remove all old packages, to solve problem caused by Bioconductor upgrade
  # create a list of all installed packages
  ip <- as.data.frame(installed.packages())
  # head(ip)
  # if you use MRO, make sure that no packages in this library will be removed
  ip <- subset(ip, !grepl("MRO", ip$LibPath))
  # we don't want to remove base or recommended packages either\
  ip <- ip[!(ip[,"Priority"] %in% c("base", "recommended")),]
  # determine the library where the packages are installed
  path.lib <- unique(ip$LibPath)
  # create a vector with all the names of the packages you want to remove
  pkgs.to.remove <- ip[,1]
  # head(pkgs.to.remove)
  # remove the packages
  sapply(pkgs.to.remove, remove.packages, lib = path.lib)
}

###############################################################################
#  Install packages
###############################################################################
install.packages(c("remotes", "BiocManager"), ask = FALSE, dependencies=TRUE, quiet=TRUE)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
new.bio.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
notInstalledPackageCount <- length(new.packages) + length(new.bio.packages)

while(notInstalledPackageCount != 0){
  # CRAN
  if(length(new.packages)){ 
	install.packages(
	  new.packages, 
	  repos="http://cran.rstudio.com/", 
	  dependencies=TRUE, 
	  quiet=TRUE
	)
  }
  #Bioconductor
  if(length(new.bio.packages)){
    BiocManager::install(
	  new.bio.packages, 
	  ask = FALSE, 
	  dependencies=TRUE, 
	  quiet=TRUE
	)
  }
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  new.bio.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
  if( notInstalledPackageCount == length(new.packages) + length(new.bio.packages) ) {
    #no new package installed.
    break
  } else {
    notInstalledPackageCount = length(new.packages) + length(new.bio.packages)
  }
}

###############################################################################
#  Install special packages
###############################################################################
#PGSEA is deprecated since Bioconductor 3.12. So we have to install manually from source.
BiocManager::install(c("GO.db", "annaffy")) # required by PGSEA
install.packages("http://www.bioconductor.org/packages//2.11/data/annotation/src/contrib/KEGG.db_2.8.0.tar.gz", repos=NULL, type="source")
install.packages("https://bioconductor.org/packages/3.10/bioc/src/contrib/PGSEA_1.60.0.tar.gz", 
                 repos=NULL, type="source")
install.packages("WGCNA") # requires three Bioconductor packages: GO.db, impute, preprocessCore
list.of.bio.packages = c(list.of.bio.packages, "PGSEA") # add package for testing

# install from GitHub
remotes::install_github("espors/ottoPlots") # for download plots

if(0)
remotes::install_github(
  "espors/idepGolem@test_depolyment", 
  force = TRUE
  ask = FALSE, 
  dependencies = TRUE, 
  quiet = FALSE, 
  upgrade = "never"  # do not ask for upgrade for dependencies
  )
list.of.packages = c(list.of.packages, "ottoPlots", "idepGolem", "WGCNA") # add package for testing

###############################################################################
#  Test if packages are installed correctly
###############################################################################
suc = unlist ( lapply(list.of.packages, require, character.only = TRUE) )
if(sum(suc) < length(list.of.packages) )
  cat ("\n\nWarnning!!!!!! These R packages cannot be loaded:", list.of.packages[!suc] )

suc = unlist ( lapply(list.of.bio.packages, require, character.only = TRUE) )
if(sum(suc) < length(list.of.bio.packages) )
  cat ("\n\nWarnning!!!!!! These Bioconductor packages cannot be loaded:", list.of.bio.packages[!suc] )

