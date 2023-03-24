if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("devtools")    # only if devtools not yet installed
BiocManager::install("pachterlab/sleuth")

BiocManager::install("COMBINE-lab/wasabi")


devtools::install_github("stephenturner/annotables")
