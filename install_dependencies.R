# Title     : Install R package dependencies
# Objective : prepare environment to run the R scripts
# Created by: todor
# Created on: 10.02.21

#strsplit(R.version.string," ")[[1]][3]

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
biocVersion=ifelse(args[1]=="4.1.0","3.14","3.16")
#usrLibPath=paste0(path.expand("~"),"/R/",R.Version()$platform,"-library/",R.Version()$major,".",gsub("\\..*$","",R.Version()$minor))

print(paste0("Using R version: ",args[1]," and biocVersion: ",biocVersion))

install_dependencies <- function(packages){
  for (package in packages){
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package)#,lib=usrLibPath)
    } else {
      print(paste("Package", package, " already installed."))
    }

  }
}


install_Bioc <- function(biocVersion){
    install.packages("BiocManager")#,lib=usrLibPath)
    BiocManager::install(version = biocVersion, ask=F)
}

install_Bioc_dependencies <- function(packages){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
    BiocManager::install(version = biocVersion, ask=F)#,lib=usrLibPath)
  }
  for (package in packages){
    if (!requireNamespace(package, quietly = TRUE)) {
      BiocManager::install(package)#,lib=usrLibPath)
    } else {
      print(paste("Bioconductor package", package, " already installed."))
    }

  }
}

install_dependencies(packages = c("devtools","Matrix","Rcpp","RColorBrewer", "ggplot2", "PoiClaClu", "pheatmap",  "tidyr", "gplots", "ggpubr",
                                  "R.utils", "forcats","ashr"))

install_Bioc(biocVersion)

install_Bioc_dependencies(packages = c("Organism.dplyr", "DESeq2", "GenomicRanges", "GenomicFeatures", "tximport",
                                       "BSgenome.Celegans.UCSC.ce11", "affy", "TxDb.Celegans.UCSC.ce11.refGene",
                                       "TxDb.Celegans.UCSC.ce11.ensGene", "rtracklayer", "apeglm"))

devtools::install_github('kevinblighe/EnhancedVolcano')
