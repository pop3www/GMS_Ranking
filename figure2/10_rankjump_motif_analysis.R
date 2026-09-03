#!/usr/bin/env Rscript

SCRIPT_VERSION <- "2026-01-09 v10 (robust optparse + env fallbacks; force numeric TPM; stable args)"

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(ggplot2)
})

message("[10] script version: ", SCRIPT_VERSION)

# ---------------- optparse ----------------
opt_list <- list(
  make_option(c("--tpm"), type = "character", dest = "tpm",
              help = "TPM table (tsv/txt) with header; one column is gene ID"),
  make_option(c("--id-col"), type = "character", dest = "id_col", default = "gene_id",
              help = "Gene ID column in TPM table [default %default]"),
  make_option(c("--groupA-pattern"), type = "character", dest = "groupA_pattern",
              help = "Regex for groupA sample columns (e.g. 4DT_RP[0-9]+)"),
  make_option(c("--groupB-pattern"), type = "character", dest = "groupB_pattern",
              help = "Regex for groupB sample columns (e.g. 4D_RP[0-9]+)"),
  make_option(c("--bands"), type = "character", dest = "bands",
              help = "rank_bands.csv from step 06"),
  make_option(c("--fimo-hits"), type = "character", dest = "fimo_hits",
              help = "fimo_hits_with_bands.csv from step 07"),
  make_option(c("--label"), type = "character", dest = "label", default = "",
              help = "Label for plots [optional]"),
  make_option(c("--top-frac"), type = "double", dest = "top_frac", default = 0.01,
              help = "Top fraction of genes to define UpJump/DownJump [default %default]"),
  make_option(c("--min-baseline-tpm"), type = "double", dest = "min_baseline_tpm", default = 1,
              help = "Require baseline (groupB) mean TPM >= this [default %default]"),
  make_option(c("--out"), type = "character", dest = "out",
              help = "Output directory")
)

parser <- OptionParser(option_list = opt_list)

# NOTE: Some clusters have older optparse/getopt behavior that can return
# dotted names, length-0 vectors, or nested lists. We defensively normalize.
opt0 <- parse_args(parser, positional_arguments = FALSE)
if (is.list(opt0) && !is.null(opt0$options) && is.list(opt0$options)) {
  # extremely defensive: handle parse_args(positional_arguments=TRUE) behavior
  opt0 <- opt0$options
}

get_opt1 <- function(key, default = NULL) {
  # 1) direct
  if (!is.null(opt0[[key]]) && length(opt0[[key]]) == 1) return(opt0[[key]])
  # 2) dotted variant (some optparse versions)
  key_dot <- gsub("_", ".", key, fixed = TRUE)
  if (!is.null(opt0[[key_dot]]) && length(opt0[[key_dot]]) == 1) return(opt0[[key_dot]])
  # 3) env fallback (optional)
  env_key <- toupper(key)
  env_key <- gsub("\\.", "_", env_key)
  if (nzchar(Sys.getenv(env_key, ""))) return(Sys.getenv(env_key))
  return(default)
}

# Normalize regex strings coming from bash (avoid doubled escaping).
norm_rx <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  x <- as.character(x)[1]
  x <- gsub("\\\\\\\\", "\\\\", x)  # \\\\d -> \\d
  x
}

opt <- list(
  tpm              = get_opt1("tpm", ""),
  id_col           = get_opt1("id_col", "gene_id"),
  groupA_pattern   = norm_rx(get_opt1("groupA_pattern", "")),
  groupB_pattern   = norm_rx(get_opt1("groupB_pattern", "")),
  bands            = get_opt1("bands", ""),
  fimo_hits        = get_opt1("fimo_hits", ""),
  label            = get_opt1("label", ""),
  top_frac         = suppressWarnings(as.numeric(get_opt1("top_frac", 0.01))),
  min_baseline_tpm = suppressWarnings(as.numeric(get_opt1("min_baseline_tpm", 1))),
  out              = get_opt1("out", "")
)

raw_args <- commandArgs(trailingOnly = TRUE)

is_blank <- function(x) is.null(x) || length(x) == 0 || !nzchar(as.character(x)[1])

required <- c("tpm", "groupA_pattern", "groupB_pattern", "bands", "fimo_hits", "out")
missing <- required[vapply(required, function(k) is_blank(opt[[k]]), logical(1))]
if (length(missing) > 0) {
  message("[10] Raw args: ", paste(raw_args, collapse = " "))
  stop("[10] Missing required args: ",
       paste(paste0("--", gsub("_", "-", missing)), collapse = ", "))
}

if (is.na(opt$top_frac) || !is.finite(opt$top_frac) || opt$top_frac <= 0 || opt$top_frac >= 1) {
  message("[10] WARNING: invalid --top-frac; using 0.01")
  opt$top_frac <- 0.01
}
if (is.na(opt$min_baseline_tpm) || !is.finite(opt$min_baseline_tpm) || opt$min_baseline_tpm < 0) {
  message("[10] WARNING: invalid --min-baseline-tpm; using 1")
  opt$min_baseline_tpm <- 1
}

if (!nzchar(opt$label)) opt$label <- "rankjump"

dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)

message("[10] label: ", opt$label)
message("[10] TPM:   ", opt$tpm)
message("[10] bands: ", opt$bands)
message("[10] hits:  ", opt$fimo_hits)
message("[10] out:   ", opt$out)
message("[10] Regex after normalization:\n  groupA: ", opt$groupA_pattern, "\n  groupB: ", opt$groupB_pattern)

# ---------------- load TPM ----------------
tpm <- fread(opt$tpm)
if (!(opt$id_col %in% names(tpm))) {
  stop("[10] id-col '", opt$id_col, "' not found in TPM. Available: ", paste(names(tpm), collapse = ", "))
}

cn <- names(tpm)
ixA <- which(grepl(opt$groupA_pattern, cn))
ixB <- which(grepl(opt$groupB_pattern, cn))

message("[10] Column matches")
message("  groupA (", length(ixA), "): ", paste(cn[ixA], collapse = ", "))
message("  groupB (", length(ixB), "): ", paste(cn[ixB], collapse = ", "))

if (length(ixA) < 1 || length(ixB) < 1) stop("[10] Need at least one sample in each group")

# Force numeric TPM columns (fread usually gets this right, but mixed types turn
# as.matrix() into character and break rowMeans).
num_cols <- setdiff(cn, opt$id_col)
for (cc in num_cols) {
  if (!is.numeric(tpm[[cc]])) {
    suppressWarnings(set(tpm, j = cc, value = as.numeric(tpm[[cc]])))
  }
}

colsA <- cn[ixA]
colsB <- cn[ixB]

A_mean <- rowMeans(as.matrix(tpm[, ..colsA]), na.rm = TRUE)
B_mean <- rowMeans(as.matrix(tpm[, ..colsB]), na.rm = TRUE)

dt <- data.table(gene = tpm[[opt$id_col]], A_mean = A_mean, B_mean = B_mean)
dt[, keep := is.finite(B_mean) & B_mean >= opt$min_baseline_tpm]
dt <- dt[keep == TRUE]

if (nrow(dt) < 10) {
  stop("[10] Too few genes after baseline filter (n=", nrow(dt), ")")
}

# Rank (1 = highest expression). Use ties.method='average' for stability.
dt[, rank_A := frank(-A_mean, ties.method = "average")]
dt[, rank_B := frank(-B_mean, ties.method = "average")]
dt[, rank_jump := rank_B - rank_A]  # + = moved up in A

n <- nrow(dt)
n_sel <- ceiling(opt$top_frac * n)
n_sel <- max(1, min(n_sel, n))
message("[10] Selecting top_frac=", opt$top_frac, " => n_sel=", n_sel, " of n=", n)

dt[, jump := "Other"]
dt[order(-rank_jump)][seq_len(n_sel), jump := "UpJump"]
dt[order(rank_jump)][seq_len(n_sel), jump := "DownJump"]

fwrite(dt, file.path(opt$out, "rankjump_gene_table.csv"))

# ---------------- load bands + fimo hits ----------------
bands <- fread(opt$bands)
if (!all(c("gene", "band") %in% names(bands))) {
  stop("[10] bands file must contain columns: gene, band")
}

fh <- fread(opt$fimo_hits)
if (nrow(fh) == 0) {
  message("[10] fimo_hits is empty; writing placeholder outputs and exiting")
  fwrite(data.table(), file.path(opt$out, "rankjump_motif_or.csv"))
  quit(status = 0)
}

req_fh <- c("motif", "gene", "band")
if (!all(req_fh %in% names(fh))) {
  stop("[10] fimo_hits must contain columns: ", paste(req_fh, collapse = ", "),
       ". Found: ", paste(names(fh), collapse = ", "))
}

# Keep only genes that survived TPM baseline filter
dt <- merge(dt, bands[, .(gene, band)], by = "gene", all.x = TRUE)
dt <- dt[!is.na(band)]

gene_motif <- unique(fh[, .(gene, motif)])
gene_motif[, has_motif := TRUE]

dt_filt <- merge(dt, gene_motif, by = "gene", all.x = TRUE)
dt_filt[is.na(has_motif), has_motif := FALSE]

# Count motif presence by band and jump class
counts <- dt_filt[, .(
  n = .N,
  n_has = sum(has_motif, na.rm = TRUE)
), by = .(motif, band, jump)]

wide <- dcast(counts, motif + band ~ jump, value.var = c("n", "n_has"), fill = 0)

calc_or <- function(a, b, c, d) {
  # Haldane-Anscombe correction for zeros
  a2 <- a + 0.5; b2 <- b + 0.5; c2 <- c + 0.5; d2 <- d + 0.5
  OR <- (a2 / b2) / (c2 / d2)
  se <- sqrt(1 / a2 + 1 / b2 + 1 / c2 + 1 / d2)
  lo <- exp(log(OR) - 1.96 * se)
  hi <- exp(log(OR) + 1.96 * se)
  list(OR = OR, lo = lo, hi = hi)
}

or_rows <- list()
for (i in seq_len(nrow(wide))) {
  r <- wide[i]
  # UpJump vs Other
  a <- r$n_has_UpJump; b <- r$n_UpJump - r$n_has_UpJump
  c <- r$n_has_Other;  d <- r$n_Other  - r$n_has_Other
  o1 <- calc_or(a, b, c, d)
  or_rows[[length(or_rows) + 1]] <- data.table(
    motif = r$motif, band = r$band, contrast = "UpJump_vs_Other",
    a = a, b = b, c = c, d = d,
    OR = o1$OR, lo = o1$lo, hi = o1$hi
  )
  # DownJump vs Other
  a <- r$n_has_DownJump; b <- r$n_DownJump - r$n_has_DownJump
  c <- r$n_has_Other;    d <- r$n_Other    - r$n_has_Other
  o2 <- calc_or(a, b, c, d)
  or_rows[[length(or_rows) + 1]] <- data.table(
    motif = r$motif, band = r$band, contrast = "DownJump_vs_Other",
    a = a, b = b, c = c, d = d,
    OR = o2$OR, lo = o2$lo, hi = o2$hi
  )
  # UpJump vs DownJump
  a <- r$n_has_UpJump;   b <- r$n_UpJump   - r$n_has_UpJump
  c <- r$n_has_DownJump; d <- r$n_DownJump - r$n_has_DownJump
  o3 <- calc_or(a, b, c, d)
  or_rows[[length(or_rows) + 1]] <- data.table(
    motif = r$motif, band = r$band, contrast = "UpJump_vs_DownJump",
    a = a, b = b, c = c, d = d,
    OR = o3$OR, lo = o3$lo, hi = o3$hi
  )
}

or_dt <- rbindlist(or_rows, fill = TRUE)
or_dt[, band := factor(band, levels = c("head", "mid", "tail"))]
fwrite(or_dt, file.path(opt$out, "rankjump_motif_or.csv"))

# ---------------- plot ----------------
if (nrow(or_dt) > 0) {
  p <- ggplot(or_dt, aes(band, OR)) +
    geom_hline(yintercept = 1, linetype = 2, linewidth = 0.4) +
    geom_point() +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
    facet_grid(contrast ~ motif, scales = "free_y") +
    scale_y_log10() +
    labs(
      title = paste0("Rank-jump motif enrichment (", opt$label, ")"),
      x = NULL, y = "OR (log scale)"
    ) +
    theme_minimal(base_size = 11)
  ggsave(file.path(opt$out, "Supp_rankjump_motif_OR.png"), p, width = 9.5, height = 6.5, dpi = 300)
}

message("[10] wrote: ", file.path(opt$out, "rankjump_gene_table.csv"), ", ", file.path(opt$out, "rankjump_motif_or.csv"))
