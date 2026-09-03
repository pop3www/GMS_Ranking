#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

SCRIPT_VERSION <- "2026-06-12 v39 (Fig2B/C naming; priming-dependent component; composite-ready PNG+PDF+SVG)"
message("[11] script version: ", SCRIPT_VERSION)

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(argv) {
  out <- list(); i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (startsWith(a, "--")) {
      a2 <- sub("^--", "", a)
      if (grepl("=", a2, fixed = TRUE)) {
        sp <- strsplit(a2, "=", fixed = TRUE)[[1]]
        key <- sp[[1]]; val <- paste(sp[-1], collapse = "=")
      } else {
        key <- a2
        if (i < length(argv) && !startsWith(argv[[i + 1L]], "--")) { val <- argv[[i + 1L]]; i <- i + 1L } else { val <- "TRUE" }
      }
      key <- tolower(gsub("_", "-", key))
      out[[key]] <- val
    }
    i <- i + 1L
  }
  out
}
amap <- parse_args(args)
get_arg <- function(keys, default = NA_character_, required = FALSE) {
  keys <- tolower(gsub("_", "-", keys))
  for (k in keys) if (!is.null(amap[[k]]) && nzchar(amap[[k]])) return(amap[[k]])
  if (required) stop("[11] Missing required arg: --", keys[[1]], call. = FALSE)
  default
}

extract_ensg <- function(x) {
  x <- as.character(x)
  has <- !is.na(x) & grepl("ENSG[0-9]+", x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  out[has] <- sub("^.*?(ENSG[0-9]+).*$", "\\1", x[has], perl = TRUE)
  out
}

save_plot_multi <- function(p, stem, out_dir, width = 6, height = 4) {
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(stem, ".svg")), p, width = width, height = height)
  invisible(TRUE)
}
copy_stem <- function(out_dir, from_stem, to_stem) {
  for (ext in c("png", "pdf", "svg")) {
    src <- file.path(out_dir, paste0(from_stem, ".", ext))
    dst <- file.path(out_dir, paste0(to_stem, ".", ext))
    if (file.exists(src)) file.copy(src, dst, overwrite = TRUE)
  }
}

rankA <- get_arg(c("ranka", "rank-a"), required = TRUE)
rankB <- get_arg(c("rankb", "rank-b"), required = TRUE)
labelA <- get_arg(c("labela", "label-a"), default = "DT_vs_D")
labelB <- get_arg(c("labelb", "label-b"), default = "Tam_vs_Ctrl")
band_source <- toupper(get_arg(c("band-source", "bandsource"), default = "B"))
motif_hits_path <- get_arg(c("motif-hits", "motifhits"), required = TRUE)
out_dir <- get_arg(c("out", "outdir", "out-dir"), required = TRUE)

for (p in c(rankA, rankB, motif_hits_path)) if (!file.exists(p)) stop("[11] input file not found: ", p, call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

A <- fread(rankA)
B <- fread(rankB)

id_col <- function(dt) {
  if ("ensg" %in% names(dt)) return("ensg")
  if ("gene" %in% names(dt)) return("gene")
  if ("gene_id" %in% names(dt)) return("gene_id")
  names(dt)[1]
}
lfc_col <- function(dt) {
  if ("lfc" %in% names(dt)) return("lfc")
  low <- tolower(names(dt))
  if ("log2fc" %in% low) return(names(dt)[low == "log2fc"][1])
  stop("[11] Could not find lfc/log2fc column", call. = FALSE)
}

idA <- id_col(A); idB <- id_col(B)
lfcA_col <- lfc_col(A); lfcB_col <- lfc_col(B)
A[, ensg := extract_ensg(get(idA))]
B[, ensg := extract_ensg(get(idB))]
A <- A[!is.na(ensg) & !is.na(get(lfcA_col)), .(ensg, lfcA = as.numeric(get(lfcA_col)))]
B <- B[!is.na(ensg) & !is.na(get(lfcB_col)), .(ensg, lfcB = as.numeric(get(lfcB_col)))]
A <- A[!duplicated(ensg)]
B <- B[!duplicated(ensg)]
AB <- merge(A, B, by = "ensg", all = FALSE)
if (!nrow(AB)) stop("[11] No overlap between rankA/rankB after ENSG parsing", call. = FALSE)

band_tbl <- if (band_source == "A") fread(rankA) else fread(rankB)
if (!"band" %in% names(band_tbl)) stop("[11] band-source table has no 'band' column", call. = FALSE)
id_band <- id_col(band_tbl)
band_tbl[, ensg := extract_ensg(get(id_band))]
band_tbl <- band_tbl[!is.na(ensg), .(ensg, band = as.character(band))]
band_tbl <- band_tbl[!duplicated(ensg)]
AB <- merge(AB, band_tbl, by = "ensg", all.x = TRUE)

# Motif classes from motif hits.
mh <- fread(motif_hits_path)
if (nrow(mh)) {
  mh_type_col <- if ("type" %in% names(mh)) "type" else names(mh)[2]
  mh_id_col <- if ("gene_core" %in% names(mh)) "gene_core" else if ("ensg" %in% names(mh)) "ensg" else if ("gene" %in% names(mh)) "gene" else names(mh)[1]
  if (mh_id_col %in% c("gene_core", "ensg")) {
    mh[, ensg := as.character(get(mh_id_col))]
  } else {
    mh[, ensg := extract_ensg(get(mh_id_col))]
  }
  mh[, type := toupper(as.character(get(mh_type_col)))]
  mh <- mh[!is.na(ensg) & type %in% c("EBOX", "TETO")]
  pres <- unique(mh[, .(ensg, type)])
  pres[, has := TRUE]
  wide <- dcast(pres, ensg ~ type, value.var = "has", fill = FALSE)
  AB <- merge(AB, wide, by = "ensg", all.x = TRUE)
}
if (!"EBOX" %in% names(AB)) AB[, EBOX := FALSE]
if (!"TETO" %in% names(AB)) AB[, TETO := FALSE]
AB[is.na(EBOX), EBOX := FALSE]
AB[is.na(TETO), TETO := FALSE]
AB[, motif_class := fifelse(EBOX & TETO, "Both", fifelse(EBOX, "EBOX", fifelse(TETO, "TETO", "None")))]

AB[, dd := lfcA - lfcB]
AB[, band := factor(tolower(trimws(band)), levels = c("head", "mid", "tail"))]
AB[, motif_class := factor(motif_class, levels = c("None", "EBOX", "TETO", "Both"))]

# ---- Fig2B: priming-dependent component ----
pB <- ggplot(AB[!is.na(band)], aes(x = band, y = dd)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45) +
  geom_violin(trim = TRUE, scale = "width", fill = "grey92", color = "grey55", linewidth = 0.35) +
  geom_boxplot(width = 0.16, outlier.size = 0.25, linewidth = 0.35) +
  facet_wrap(~motif_class, nrow = 1) +
  theme_bw(base_size = 15) +
  labs(
    title = "Priming-dependent component",
    subtitle = paste0("ΔΔlog2FC = log2FC(", labelA, ") − log2FC(", labelB, ")"),
    x = "Baseline rank band",
    y = "ΔΔlog2FC"
  ) +
  theme(plot.title = element_text(face = "bold", size = 18), plot.subtitle = element_text(size = 10.5, color = "grey25"), strip.text = element_text(face = "bold", size = 14), panel.grid.minor = element_blank(), plot.margin = margin(6, 8, 6, 12))

save_plot_multi(pB, "priming_dependent_component", out_dir, width = 10.8, height = 4.25)
copy_stem(out_dir, "priming_dependent_component", "Fig2B_priming_dependent_component")
# Backward-compatible source name, but not final manuscript name.
copy_stem(out_dir, "priming_dependent_component", "priming_boost_by_band")

# ---- Fig2C: gene-wise correspondence ----
spearman <- suppressWarnings(cor(AB$lfcA, AB$lfcB, method = "spearman", use = "pairwise.complete.obs"))
pearson <- suppressWarnings(cor(AB$lfcA, AB$lfcB, method = "pearson", use = "pairwise.complete.obs"))
ann <- sprintf("Spearman ρ = %.3f\nPearson r = %.3f", spearman, pearson)

# Use unprimed response on x and primed response on y.
pC <- ggplot(AB, aes(x = lfcB, y = lfcA)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey45") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey45") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", linewidth = 0.35, color = "grey60") +
  geom_point(alpha = 0.25, size = 0.45) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7, color = "#2C7FB8") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.08, vjust = -0.6, label = ann, size = 3.5) +
  theme_bw(base_size = 15) +
  labs(
    title = "Primed vs unprimed response",
    x = paste0("log2FC ", labelB, " (unprimed)"),
    y = paste0("log2FC ", labelA, " (primed)")
  ) +
  theme(plot.title = element_text(face = "bold", size = 18), panel.grid.minor = element_blank(), plot.margin = margin(6, 8, 6, 8))

save_plot_multi(pC, "priming_lfc_scatter", out_dir, width = 5.7, height = 4.25)
copy_stem(out_dir, "priming_lfc_scatter", "Fig2C_DTvsD_vs_TamvsCtrl_scatter")
# Backward-compatible old label.
copy_stem(out_dir, "priming_lfc_scatter", "Fig2B_priming_lfc_scatter")

stats <- data.table(n = nrow(AB), spearman_rho = spearman, pearson_r = pearson)
fwrite(stats, file.path(out_dir, "priming_component_summary.csv"))
fwrite(AB, file.path(out_dir, "priming_component_per_gene.tsv"), sep = "\t")
message(sprintf("[11] Done. n=%d | Spearman rho=%.3f | Pearson r=%.3f", nrow(AB), spearman, pearson))
