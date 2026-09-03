#!/usr/bin/env Rscript

# combine_pairs_rankjump_OR.R
#
# Combine per-contrast rank-jump OR tables (from 10_rankjump_motif_analysis.R)
# and make a single figure across all 4h contrasts.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript combine_pairs_rankjump_OR.R <outroot>")
}

outroot <- args[1]

dirs <- list.dirs(outroot, recursive = TRUE, full.names = TRUE)
dirs <- dirs[grepl("fimo_summary_4h_", basename(dirs))]

if (length(dirs) == 0) {
  stop("No fimo_summary_4h_* directories found under: ", outroot)
}

read_one <- function(d){
  f <- file.path(d, "rank_jump_or_by_band.csv")
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  dt[, dir := basename(d)]
  # Prefer explicit label field; otherwise derive from dirname
  if (!"label" %in% names(dt)) dt[, label := sub("^fimo_summary_4h_", "", basename(d))]
  if (all(is.na(dt$label)) || all(dt$label=="")) dt[, label := sub("^fimo_summary_4h_", "", basename(d))]
  dt
}

lst <- lapply(dirs, read_one)
lst <- lst[!vapply(lst, is.null, logical(1))]
if (length(lst) == 0) {
  stop("Found contrast directories, but none contained rank_jump_or_by_band.csv")
}

all <- rbindlist(lst, use.names=TRUE, fill=TRUE)

# Clean factors
all[, band := factor(band, levels=c("head","mid","tail"))]

# Stable flag: mimic per-contrast logic
all[, min_cell := pmin(a,b,c,d)]
all[, stable := (n_upjump >= 30 & n_downjump >= 30 & min_cell >= 5)]

# Contrast label order: keep in directory order if possible
lab_order <- unique(all$label)
all[, label := factor(label, levels=lab_order)]

p <- ggplot(all, aes(x=label, y=OR, alpha=stable)) +
  geom_hline(yintercept=1, linetype=2, linewidth=0.4) +
  geom_point(position=position_dodge(width=0.6)) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.18, position=position_dodge(width=0.6)) +
  facet_grid(type ~ band, scales="free_y") +
  scale_alpha_manual(values=c(`TRUE`=1.0, `FALSE`=0.25), guide="none") +
  scale_y_log10() +
  labs(
    title = "Motif enrichment (UpJump vs DownJump) across 4h contrasts",
    subtitle = "Gene-level presence (>=1 hit). Faded = unstable (low-count) estimates",
    x = "contrast",
    y = "OR (log scale)"
  ) +
  theme_minimal(base_size=12) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

out_png <- file.path(outroot, "pairwise_rankjump_OR_up_vs_downjump_combined_all.png")
ggsave(out_png, p, width=13.5, height=5.8, dpi=300)

fwrite(all, file.path(outroot, "pairwise_rankjump_OR_up_vs_downjump_combined_all.tsv"), sep="\t")

message("[combine_rankjump] wrote: ", out_png)
