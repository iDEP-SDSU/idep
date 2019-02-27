# Required update
#sudo apt-get install libcurl4-openssl-dev libxml2-dev

# bioconductor packages
source("https://bioconductor.org/biocLite.R")
biocLite(c( "limma", "DESeq2","edgeR","gage", "PGSEA", "fgsea", "ReactomePA", "pathview", "PREDA", "PREDAsampledata","sfsmisc","lokern","multtest"))
# annotation packages needed by pathview
biocLite(c("org.Ag.eg.db","org.At.tair.db","org.Bt.eg.db","org.Ce.eg.db","org.Cf.eg.db",
  "org.Dm.eg.db","org.Dr.eg.db","org.EcK12.eg.db","org.EcSakai.eg.db","org.Gg.eg.db",
  "org.Hs.eg.db","org.Hs.ipi.db","org.Mm.eg.db","org.Mmu.eg.db","org.Pf.plasmo.db",
  "org.Pt.eg.db","org.Rn.eg.db","org.Sc.sgd.db","org.Sco.eg.db","org.Ss.eg.db",
  "org.Tgondii.eg.db","org.Xl.eg.db"))
