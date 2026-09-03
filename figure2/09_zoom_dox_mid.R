#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(scales)
  library(stringr)
})

opt_list <- list(
  make_option("--tpm", type="character"),
  make_option("--id-col", type="character", dest="id_col", default="gene_id"),
  make_option("--ctrl-pattern", type="character", dest="ctrl_pattern"),
  make_option("--dox-pattern", type="character", dest="dox_pattern"),
  make_option("--bands", type="character"),
  make_option("--fimo-summary", type="character", dest="fimo_summary", default=NA),
  make_option("--out", type="character", default="out/fimo_summary")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

strip_quotes <- function(x) gsub('^"|"$', '', x)
normalize_rx <- function(rx){
  if (is.na(rx) || !nzchar(rx)) return(rx)
  gsub("\\\\", "\\", rx, fixed = TRUE)
}
to_mean <- function(dt, cols){
  m <- as.matrix(dt[, ..cols])
  storage.mode(m) <- "double"
  rowMeans(m, na.rm = TRUE)
}

rx_ctrl <- normalize_rx(opt$ctrl_pattern)
rx_dox  <- normalize_rx(opt$dox_pattern)

# Read TPM
tpm <- fread(opt$tpm, sep="\t", header=TRUE, check.names=FALSE, data.table=TRUE, quote="\"")
setnames(tpm, names(tpm), strip_quotes(names(tpm)))
stopifnot(opt$id_col %in% names(tpm))
setnames(tpm, opt$id_col, "gene")
tpm[, gene := strip_quotes(gene)]

num_cols <- setdiff(names(tpm), "gene")
tpm[, (num_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols=num_cols]

ix_ctrl <- grepl(rx_ctrl, num_cols, perl=TRUE)
ix_dox  <- grepl(rx_dox,  num_cols, perl=TRUE)

cat("[09] Columns matched: ctrl=", sum(ix_ctrl), ", dox=", sum(ix_dox), "\n", sep="")
stopifnot(sum(ix_ctrl)>0, sum(ix_dox)>0)

ctrl_tpm <- to_mean(tpm, num_cols[ix_ctrl])
dox_tpm  <- to_mean(tpm, num_cols[ix_dox])

dt <- data.table(
  gene = tpm$gene,
  ctrl = ctrl_tpm,
  dox  = dox_tpm,
  lfc  = log2((dox_tpm+1)/(ctrl_tpm+1))
)

bands <- fread(opt$bands)
dt <- merge(dt, bands[,.(gene, band)], by="gene", all.x=TRUE)
dt <- dt[band=="mid"]

p <- ggplot(dt, aes(x=ctrl, y=dox)) +
  geom_point(alpha=0.25, size=0.4) +
  geom_abline(slope=1, intercept=0, linetype=2, color="grey40") +
  scale_x_log10(labels=comma_format(accuracy=0.1)) +
  scale_y_log10(labels=comma_format(accuracy=0.1)) +
  coord_equal() +
  labs(title="DOX priming (MID band): DOX vs CTRL",
       x="CTRL TPM (log10)", y="DOX TPM (log10)") +
  theme_minimal(base_size=11)

ggsave(file.path(opt$out,"DOX_mid_uplift.png"), p, width=4.6, height=3.8, dpi=300)

# Optional overlay EBOX/TETO subsets
if (!is.na(opt$fimo_summary) && file.exists(opt$fimo_summary)) {
  hits <- fread(opt$fimo_summary)
  if (all(c("gene","type") %in% names(hits))) {
    hits <- unique(hits[, .(gene, type)])
    dtx <- merge(dt, hits, by="gene", all.x=FALSE)
    p2 <- ggplot(dtx, aes(x=ctrl, y=dox, color=type)) +
      geom_point(alpha=0.35, size=0.5) +
      geom_abline(slope=1, intercept=0, linetype=2, color="grey40") +
      scale_x_log10(labels=comma_format(accuracy=0.1)) +
      scale_y_log10(labels=comma_format(accuracy=0.1)) +
      coord_equal() +
      scale_color_brewer(palette="Set1") +
      labs(title="DOX priming (MID band) by motif",
           x="CTRL TPM (log10)", y="DOX TPM (log10)") +
      theme_minimal(base_size=11)
    ggsave(file.path(opt$out,"DOX_mid_uplift_EBOX_TETO.png"), p2, width=5.2, height=3.8, dpi=300)
  }
}
