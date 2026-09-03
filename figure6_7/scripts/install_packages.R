#!/usr/bin/env Rscript

cran_pkgs <- c(
  "tidyverse", "irlba", "FNN", "cluster", "RobustRankAggreg",
  "pheatmap", "fgsea", "msigdbr", "gridExtra", "png", "svglite", "BiocManager"
)
bioc_pkgs <- c("DESeq2", "AnnotationDbi", "clusterProfiler", "org.Hs.eg.db")

install_if_missing <- function(pkgs, installer) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing)) installer(missing)
}

install_if_missing(setdiff(cran_pkgs, "BiocManager"), function(x) install.packages(x, repos = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
install_if_missing(bioc_pkgs, function(x) BiocManager::install(x, ask = FALSE, update = FALSE))

cat("Fig. 7/8 package installation check complete.\n")
