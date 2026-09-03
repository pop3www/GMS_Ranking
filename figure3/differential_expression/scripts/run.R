#!/usr/bin/env Rscript
# Figure 3 entry point. By default, this plots from a unified DE table.
# To build that table locally, run build_unified_de_results.R first, then rerun this script.

args <- commandArgs(trailingOnly = TRUE)
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(this_file)) dirname(normalizePath(this_file)) else getwd()
plot_script <- file.path(script_dir, "plot_figure3.R")
if (!file.exists(plot_script)) stop("Missing plot script: ", plot_script)

if ("--help" %in% args || "-h" %in% args) {
  cat("Figure 3 production entry point\n\n")
  cat("Typical usage:\n")
  cat("  Rscript figure3/differential_expression/scripts/build_unified_de_results.R \\\n")
  cat("    --contrast-ids DT_vs_D_4h,DT_vs_D_24h \\\n")
  cat("    --out figure3/differential_expression/source_data/unified_de_results.csv\n\n")
  cat("  Rscript figure3/differential_expression/scripts/run.R \\\n")
  cat("    --de-results figure3/differential_expression/source_data/unified_de_results.csv \\\n")
  cat("    --representative-contrast DT_vs_D_24h \\\n")
  cat("    --early-contrast DT_vs_D_4h\n\n")
  cat("The plotter removes RankCompV3 by default and expects the six manuscript methods: DESeq2, edgeR, limma/voom, RankProd, PenDA, CellComp.\n")
  quit(save = "no", status = 0)
}
status <- system2("Rscript", c(plot_script, args), stdout = "", stderr = "")
quit(save = "no", status = status)
