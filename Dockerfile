FROM rocker/shiny:latest
#FROM debian # for testing

MAINTAINER Kevin Son "eunwoo.son@sdstate.edu"
RUN apt-get update -qq && apt-get install -y \
  git-core \
  libcurl4-openssl-dev \
  libxml2-dev \
  libxml2  \
  libssl-dev \
  # https://stackoverflow.com/questions/42287164/install-udunits2-package-for-r3-3
  libudunits2-dev \ 
  libmariadb-client-lgpl-dev \
  wget \ 
  unzip

# COPY ./RSet /usr/local/src/myscripts
COPY ./classes /usr/local/src/myscripts
COPY ./shinyapps /srv/shiny-server

RUN mkdir -p /srv/data/geneInfo
RUN mkdir -p /srv/data/gmt
RUN mkdir -p /srv/data/motif
RUN mkdir -p /srv/data/pathwayDB
RUN mkdir -p /srv/data/data_go

# Install R libraries
RUN R -e 'install.packages(c("devtools"))'
RUN R -e 'install.packages(c("shiny", "shinyAce", "shinyBS", "plotly","RSQLite", "gplots", "ggplot2", "dplyr", "tidyverse","plotly","e1071", "reshape2", "DT","data.table", "Rcpp","flashClust","statmod","biclust","igraph","Rtsne"))' 
RUN R -e 'source("https://bioconductor.org/biocLite.R"); biocLite(c("limma", "DESeq2", "edgeR", "gage", "PGSEA", "fgsea", "ReactomePA", "pathview", "PREDA", "impute", "runibic","QUBIC","rhdf5", "PREDAsampledata", "sfsmisc", "lokern", "multtest", "hgu133plus2.db", "org.Ag.eg.db","org.At.tair.db","org.Bt.eg.db","org.Ce.eg.db","org.Cf.eg.db", "org.Dm.eg.db","org.Dr.eg.db","org.EcK12.eg.db","org.EcSakai.eg.db","org.Gg.eg.db","org.Hs.eg.db","org.Hs.ipi.db","org.Mm.eg.db","org.Mmu.eg.db","org.Pf.plasmo.db","org.Pt.eg.db","org.Rn.eg.db","org.Sc.sgd.db","org.Sco.eg.db","org.Ss.eg.db","org.Tgondii.eg.db","org.Xl.eg.db"), suppressUpdates = T)'
# WGCNA must be installed after AnnotationDbi, impute, GO.db, preprocessCore
RUN R -e 'install.packages(c("WGCNA"))' 

# Download Required Data
# geneInfo
# RUN wget -qO- -O tmp.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2Fupdate%2FgeneInfo.zip?alt=media&token=350524c5-a690-4253-8a12-a376ac69dc69'\
#   && unzip -f tmp.zip -d /srv/data && rm tmp.zip
# # gmt file 
# RUN wget -qO- -O tmp2.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2Fgmt.zip?alt=media&token=5ee100d1-e645-41ef-a591-7a9ba208ce3c' \
#   && unzip -f tmp2.zip -d /srv/data && rm tmp2.zip
# # motif
# RUN wget -qO- -O tmp3.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2Fmotif.zip?alt=media&token=dc1e5972-ffd9-43a1-bbcc-49b78da6f047' \
#   && unzip -f tmp3.zip -d /srv/data && rm tmp3.zip
# # pathwayDB
# RUN wget -O tmp4.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2FpathwayDB.zip?alt=media&token=e602f2f7-102a-4cc4-8412-be2b05997daa' \
#   && unzip -f tmp4.zip -d /srv/data && rm tmp4.zip
# # convertIDs
# RUN wget -qO- -O tmp5.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2Fupdate%2FconvertIDs.zip?alt=media&token=147fc26e-9199-45dc-8fce-9191d5d3a3a5' \
#   && unzip -f tmp5.zip -d /srv/data && rm tmp5.zip
# # data_go
# RUN wget -qO- -O tmp6.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2Fdata_go.zip?alt=media&token=96ddcf70-ead4-4386-b582-18afe0386b8d' \
#   && unzip -f tmp6.zip -d /srv/data && rm tmp6.zip

WORKDIR /usr/local/src/myscripts

EXPOSE 3838

# Install required R libraries
# CMD ["Rscript", "librarySetup.R"]
# CMD ["/usr/bin/shiny-server.sh"] #If you don't use docker-compose need to comment out