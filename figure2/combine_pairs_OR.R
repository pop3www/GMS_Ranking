#!/usr/bin/env Rscript

# combine_pairs_OR.R
#
# Combine per-comparison motif enrichment OR tables (Up vs Down) and plot a single summary figure.
#
# Inputs:
#   <outroot>/fimo_summary_4h_*/enrichment_odds_ratios_up_vs_down.csv
#
# Outputs written to <outroot>:
#   pairwise_OR_up_vs_down_combined.csv
#   pairwise_OR_up_vs_down_combined_all.png
#   pairwise_OR_up_vs_down_combined_all.svg
#
# Schema compatibility:
# - Newer 07_summarize_fimo.R (v26+) outputs columns:
#     odds_ratio, conf_low, conf_high, a, b, c, d
# - Older outputs may contain:
#     OR, lo, hi, n_up, n_down
# This script accepts both and standardizes internally.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: combine_pairs_OR.R <outroot>")
}

outroot <- args[[1L]]
if (!dir.exists(outroot)) {
  stop("[combine_pairs_OR] outroot does not exist: ", outroot)
}

# ---- helpers ----

standardize_or_cols <- function(dt) {
  # OR
  if (!"odds_ratio" %in% names(dt) && "OR" %in% names(dt)) dt[, odds_ratio := OR]
  # CI
  if (!"conf_low" %in% names(dt) && "lo" %in% names(dt)) dt[, conf_low := lo]
  if (!"conf_high" %in% names(dt) && "hi" %in% names(dt)) dt[, conf_high := hi]

  # Counts (for stability filter)
  if (!"n_up" %in% names(dt)) {
    if (all(c("a", "b") %in% names(dt))) {
      dt[, n_up := as.integer(a + b)]
    } else {
      dt[, n_up := NA_integer_]
    }
  }
  if (!"n_down" %in% names(dt)) {
    if (all(c("c", "d") %in% names(dt))) {
      dt[, n_down := as.integer(c + d)]
    } else {
      dt[, n_down := NA_integer_]
    }
  }
  if (!"min_cell" %in% names(dt)) {
    if (all(c("a", "b", "c", "d") %in% names(dt))) {
      dt[, min_cell := as.integer(pmin(a, b, c, d))]
    } else {
      dt[, min_cell := NA_integer_]
    }
  }

  # Required plotting cols
  req <- c("type", "band", "odds_ratio", "conf_low", "conf_high")
  missing_req <- setdiff(req, names(dt))
  if (length(missing_req)) {
    stop("[combine_pairs_OR] missing required columns: ", paste(missing_req, collapse = ", "))
  }

  dt[]
}

read_one <- function(summary_dir) {
  comp <- sub("^fimo_summary_4h_", "", basename(summary_dir))
  path <- file.path(summary_dir, "enrichment_odds_ratios_up_vs_down.csv")
  if (!file.exists(path)) return(NULL)

  dt <- fread(path)
  if (!nrow(dt)) return(NULL)

  dt <- standardize_or_cols(dt)
  dt[, comparison := comp]
  dt[]
}

# ---- collect ----

summ_dirs <- list.dirs(outroot, full.names = TRUE, recursive = FALSE)
summ_dirs <- summ_dirs[grepl("/fimo_summary_4h_", summ_dirs)]

if (!length(summ_dirs)) {
  stop("[combine_pairs_OR] No fimo_summary_4h_* directories under: ", outroot)
}

lst <- lapply(summ_dirs, read_one)
res <- rbindlist(lst, fill = TRUE)

if (!nrow(res)) {
  stop("[combine_pairs_OR] No enrichment_odds_ratios_up_vs_down.csv files found under: ", outroot)
}

# ---- stability filter ----

res[, stable := (!is.na(n_up) & !is.na(n_down) & !is.na(min_cell) & n_up >= 30 & n_down >= 30 & min_cell >= 5)]

# Write full combined table
out_csv <- file.path(outroot, "pairwise_OR_up_vs_down_combined.csv")
fwrite(res, out_csv)
message("[combine_pairs_OR] wrote: ", out_csv)

# Plot only stable + finite
plot_dt <- res[stable == TRUE]
plot_dt <- plot_dt[is.finite(odds_ratio) & odds_ratio > 0 & is.finite(conf_low) & is.finite(conf_high) & conf_low > 0 & conf_high > 0]

if (!nrow(plot_dt)) {
  message("[combine_pairs_OR] No stable rows to plot; skipping figure output.")
  quit(status = 0)
}

# Ordering helpers
plot_dt[, comparison := factor(comparison, levels = unique(comparison))]
if ("band" %in% names(plot_dt)) plot_dt[, band := factor(band, levels = c("head", "mid", "tail"))]
if ("type" %in% names(plot_dt)) plot_dt[, type := factor(type, levels = c("EBOX", "TETO"))]

p <- ggplot(plot_dt, aes(x = comparison, y = odds_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3) +
  geom_pointrange(aes(ymin = conf_low, ymax = conf_high), linewidth = 0.25) +
  scale_y_log10() +
  facet_grid(band ~ type, scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Odds ratio (Up vs Down)",
    title = "Motif enrichment OR by band (Up vs Down; stable rows only)"
  )

out_png <- file.path(outroot, "pairwise_OR_up_vs_down_combined_all.png")
out_svg <- file.path(outroot, "pairwise_OR_up_vs_down_combined_all.svg")

ggsave(out_png, p, width = 10, height = 6, dpi = 300)
ggsave(out_svg, p, width = 10, height = 6)
message("[combine_pairs_OR] wrote: ", out_png)
message("[combine_pairs_OR] wrote: ", out_svg)
