#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggseqlogo)
  library(gridExtra)
  library(scales)
})

SCRIPT_VERSION <- "2026-06-18 v42 (Fig2 D-F export; enlarged tetO weights row; All/Up/Down density; PNG+PDF+SVG)"
message("[08] script version: ", SCRIPT_VERSION)

# ---- manual args ----
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
get_arg <- function(keys, default = NA_character_, env = NULL, required = FALSE) {
  keys <- tolower(gsub("_", "-", keys))
  for (k in keys) if (!is.null(amap[[k]]) && nzchar(amap[[k]])) return(amap[[k]])
  if (!is.null(env)) {
    ev <- Sys.getenv(env, unset = "")
    if (nzchar(ev)) return(ev)
  }
  if (required) stop("[08] Missing required arg: --", keys[[1]], call. = FALSE)
  default
}

summary_dir <- get_arg(c("summary-dir", "summary"), required = TRUE)
motif_meme  <- get_arg(c("motif-meme"), env = "TETO_MOTIF", required = TRUE)
pos_weights <- get_arg(c("pos-weights"), env = "POS_WEIGHTS", required = TRUE)
outdir      <- get_arg(c("out", "outdir"), default = summary_dir)

if (!dir.exists(summary_dir)) stop("[08] summary-dir does not exist: ", summary_dir, call. = FALSE)
if (!file.exists(motif_meme)) stop("[08] motif-meme not found: ", motif_meme, call. = FALSE)
if (!file.exists(pos_weights)) stop("[08] pos-weights not found: ", pos_weights, call. = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- helpers ----
read_csv_dt <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  dt
}

save_plot_multi <- function(p, stem, w = 6, h = 4, dpi = 300) {
  ggsave(file.path(outdir, paste0(stem, ".png")), p, width = w, height = h, dpi = dpi)
  ggsave(file.path(outdir, paste0(stem, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(outdir, paste0(stem, ".svg")), p, width = w, height = h)
  invisible(TRUE)
}

copy_stem <- function(from_stem, to_stem) {
  for (ext in c("png", "pdf", "svg")) {
    src <- file.path(outdir, paste0(from_stem, ".", ext))
    dst <- file.path(outdir, paste0(to_stem, ".", ext))
    if (file.exists(src)) file.copy(src, dst, overwrite = TRUE)
  }
}

blank_plot <- function(title) {
  ggplot() + theme_void(base_size = 13) +
    annotate("text", x = 0, y = 0, label = title, size = 5) +
    xlim(-1, 1) + ylim(-1, 1)
}

std_group <- function(x) {
  y <- toupper(trimws(as.character(x)))
  out <- ifelse(y %in% c("ALL"), "All",
         ifelse(y %in% c("UP", "U"), "Up",
         ifelse(y %in% c("DOWN", "D"), "Down", as.character(x))))
  out
}

std_band <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y[!(y %in% c("head", "mid", "tail"))] <- NA_character_
  factor(y, levels = c("head", "mid", "tail"))
}

# ---- MEME parsing ----
read_meme_first_motif <- function(meme_path) {
  lines <- readLines(meme_path, warn = FALSE)
  idx <- grep("^letter-probability matrix", lines)
  if (!length(idx)) stop("[08] Failed to find letter-probability matrix in MEME file", call. = FALSE)
  i <- idx[[1]]
  header <- lines[[i]]
  w <- as.integer(sub(".*\\bw=\\s*([0-9]+).*", "\\1", header, perl = TRUE))
  if (is.na(w) || w <= 0) stop("[08] Failed to parse motif width (w)", call. = FALSE)
  motif_idx <- suppressWarnings(max(grep("^MOTIF\\s+", lines[seq_len(i)], perl = TRUE)))
  mname <- if (is.finite(motif_idx)) sub("^MOTIF\\s+", "", lines[[motif_idx]]) else "tetO motif"
  mat_lines <- character(0); j <- i + 1L
  while (j <= length(lines) && length(mat_lines) < w) {
    ln <- trimws(lines[[j]])
    if (nzchar(ln)) mat_lines <- c(mat_lines, ln)
    j <- j + 1L
  }
  if (length(mat_lines) != w) stop("[08] MEME ended before full PWM", call. = FALSE)
  mat <- do.call(rbind, lapply(mat_lines, function(ln) {
    vals <- suppressWarnings(as.numeric(strsplit(ln, "\\s+")[[1]]))
    vals <- vals[!is.na(vals)]
    if (length(vals) < 4) stop("[08] PWM row does not have 4 numeric values", call. = FALSE)
    vals[1:4]
  }))
  ppm <- t(mat)
  rownames(ppm) <- c("A", "C", "G", "T")
  colnames(ppm) <- seq_len(w)
  list(name = mname, w = w, ppm = ppm)
}

calc_ic <- function(ppm) {
  p <- pmax(ppm, 1e-9)
  H <- -colSums(p * log2(p))
  data.table(pos = seq_len(ncol(ppm)), IC = 2 - H)
}

read_pos_weights <- function(path) {
  dt <- fread(path)
  nms <- tolower(names(dt))
  if ("pos" %in% nms) dt[, pos := as.integer(get(names(dt)[which(nms == "pos")[1]]))]
  else if ("position" %in% nms) dt[, pos := as.integer(get(names(dt)[which(nms == "position")[1]]))]
  else stop("[08] weights TSV needs pos/position", call. = FALSE)
  if ("weight" %in% nms) dt[, weight := as.numeric(get(names(dt)[which(nms == "weight")[1]]))]
  else if ("w" %in% nms) dt[, weight := as.numeric(get(names(dt)[which(nms == "w")[1]]))]
  else stop("[08] weights TSV needs weight/w", call. = FALSE)
  if ("base" %in% nms) dt[, base := as.character(get(names(dt)[which(nms == "base")[1]]))] else dt[, base := NA_character_]
  dt <- dt[!is.na(pos) & !is.na(weight)]
  dt[order(pos)]
}

mot <- read_meme_first_motif(motif_meme)
wt <- read_pos_weights(pos_weights)
if (nrow(wt) != mot$w) message("[08] WARNING: weights rows (", nrow(wt), ") != motif width (", mot$w, ")")

# ---- motif/weight plots ----
logo_plot <- ggseqlogo(mot$ppm, method = "prob") +
  ggtitle("TRE3Gs-derived tetO motif") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 20), panel.grid = element_blank())
save_plot_multi(logo_plot, "Supp_tetO_logo", w = 14.8, h = 2.3)

wt_plot <- ggplot(wt, aes(x = pos, y = weight)) +
  geom_col(width = 0.75) +
  scale_x_continuous(breaks = seq_len(mot$w)) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), labels = c("0", "0.5", "1.0")) +
  coord_cartesian(ylim = c(0, max(1, max(wt$weight, na.rm = TRUE)))) +
  theme_bw(base_size = 17) +
  labs(title = "tetO positional weights", x = "Position", y = "Weight") +
  theme(
    plot.title = element_text(face = "bold", size = 19),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    panel.grid.minor = element_blank(),
    plot.margin = margin(4, 6, 4, 6)
  )
save_plot_multi(wt_plot, "Supp_tetO_weights", w = 14.8, h = 1.75)

ic_dt <- calc_ic(mot$ppm)
ic_plot <- ggplot(ic_dt, aes(x = pos, y = IC)) +
  geom_col(width = 0.75) +
  scale_x_continuous(breaks = seq_len(mot$w)) +
  theme_bw(base_size = 18) +
  labs(title = "tetO information content", x = "Position", y = "Information content (bits)") +
  theme(plot.title = element_text(face = "bold", size = 20), panel.grid.minor = element_blank())
save_plot_multi(ic_plot, "Supp_tetO_IC", w = 7.2, h = 2.8)

combo <- gridExtra::arrangeGrob(
  ggplotGrob(logo_plot + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(2, 4, 2, 4))),
  ggplotGrob(wt_plot + theme(plot.margin = margin(4, 6, 4, 6))),
  ncol = 1, heights = c(1.55, 1.35)
)
png(file.path(outdir, "Supp_tetO_logo_with_weights.png"), width = 14.8, height = 3.35, units = "in", res = 300)
grid::grid.draw(combo); dev.off()
pdf(file.path(outdir, "Supp_tetO_logo_with_weights.pdf"), width = 14.8, height = 3.35, useDingbats = FALSE)
grid::grid.draw(combo); dev.off()
svg(file.path(outdir, "Supp_tetO_logo_with_weights.svg"), width = 14.8, height = 3.35)
grid::grid.draw(combo); dev.off()
copy_stem("Supp_tetO_logo_with_weights", "Fig2D_tetO_logo_weights")

# ---- summary tables ----
hit_density <- read_csv_dt(file.path(summary_dir, "hit_density_by_band.csv"))
hit_gene_frac <- read_csv_dt(file.path(summary_dir, "hit_gene_fraction_by_band.csv"))
or_all <- read_csv_dt(file.path(summary_dir, "enrichment_odds_ratios.csv"))
or_ud <- read_csv_dt(file.path(summary_dir, "enrichment_odds_ratios_up_vs_down.csv"))
or_bb <- read_csv_dt(file.path(summary_dir, "enrichment_odds_ratios_by_band.csv"))

prep_band_dt <- function(dt) {
  if (is.null(dt)) return(NULL)
  if (!"band" %in% names(dt)) return(NULL)
  if (!"type" %in% names(dt)) dt[, type := "Motif"]
  if (!"updown_group" %in% names(dt)) dt[, updown_group := "All"]
  if (!"y" %in% names(dt)) {
    if (all(c("hits", "promoters") %in% names(dt))) dt[, y := ifelse(promoters > 0, hits / promoters * 1000, NA_real_)]
    else if (all(c("genes_hit", "genes_total") %in% names(dt))) dt[, y := ifelse(genes_total > 0, genes_hit / genes_total, NA_real_)]
    else return(NULL)
  }
  dt[, band := std_band(band)]
  dt[, type := factor(toupper(as.character(type)), levels = c("EBOX", "TETO"))]
  dt[, updown_group := factor(std_group(updown_group), levels = c("All", "Up", "Down"))]
  dt <- dt[!is.na(band) & !is.na(type) & !is.na(updown_group)]
  dt
}

hit_density <- prep_band_dt(hit_density)
hit_gene_frac <- prep_band_dt(hit_gene_frac)

# Fig2E / Supp hit density: All+Up+Down, log10 y as requested in legend.
if (!is.null(hit_density) && nrow(hit_density)) {
  plot_dt <- hit_density[!is.na(y) & y > 0]
  if (nrow(plot_dt)) {
    dodge <- position_dodge(width = 0.55)
    p <- ggplot(plot_dt, aes(x = band, y = y, color = updown_group, shape = updown_group)) +
      geom_point(position = dodge, size = 2.8, stroke = 1) +
      facet_wrap(~type, nrow = 1) +
      scale_y_log10(labels = label_number()) +
      theme_bw(base_size = 18) +
      labs(x = "Baseline rank band", y = "Hits per 1,000 promoters (log10)", color = "Group", shape = "Group") +
      theme(strip.text = element_text(face = "bold", size = 13), panel.grid.minor = element_blank())
    save_plot_multi(p, "Supp_hit_density_by_band", w = 7.9, h = 3.55)
    copy_stem("Supp_hit_density_by_band", "Fig2E_hit_density_by_band")
  } else {
    save_plot_multi(blank_plot("No positive hit-density rows"), "Supp_hit_density_by_band", w = 7.2, h = 3.8)
  }
} else {
  save_plot_multi(blank_plot("Missing hit density data"), "Supp_hit_density_by_band", w = 7.2, h = 3.8)
}

# Gene fraction is source QC, not final Fig2 panel, but keep/export.
if (!is.null(hit_gene_frac) && nrow(hit_gene_frac)) {
  plot_dt <- hit_gene_frac[!is.na(y)]
  dodge <- position_dodge(width = 0.55)
  p <- ggplot(plot_dt, aes(x = band, y = y, color = updown_group, shape = updown_group)) +
    geom_point(position = dodge, size = 2.8, stroke = 1) +
    facet_wrap(~type, nrow = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    theme_bw(base_size = 18) +
    labs(x = "Baseline rank band", y = "Genes with ≥1 hit", color = "Group", shape = "Group") +
    theme(strip.text = element_text(face = "bold", size = 13), panel.grid.minor = element_blank())
  save_plot_multi(p, "Supp_gene_fraction_by_band", w = 7.2, h = 3.8)
} else {
  save_plot_multi(blank_plot("Missing gene-fraction data"), "Supp_gene_fraction_by_band", w = 7.2, h = 3.8)
}

prep_or <- function(dt) {
  if (is.null(dt) || !nrow(dt)) return(NULL)
  if (!"odds_ratio" %in% names(dt) && "OR" %in% names(dt)) setnames(dt, "OR", "odds_ratio")
  if (!"conf_low" %in% names(dt) && "lo" %in% names(dt)) setnames(dt, "lo", "conf_low")
  if (!"conf_high" %in% names(dt) && "hi" %in% names(dt)) setnames(dt, "hi", "conf_high")
  req <- c("type", "band", "odds_ratio", "conf_low", "conf_high")
  if (length(setdiff(req, names(dt)))) return(NULL)
  dt[, band := std_band(band)]
  dt[, type := factor(toupper(as.character(type)), levels = c("EBOX", "TETO"))]
  dt <- dt[!is.na(band) & !is.na(type) & is.finite(odds_ratio) & odds_ratio > 0 & is.finite(conf_low) & conf_low > 0 & is.finite(conf_high) & conf_high > 0]
  dt
}

plot_or <- function(dt, stem, title = "Motif enrichment by baseline band") {
  dt <- prep_or(dt)
  if (is.null(dt) || !nrow(dt)) {
    save_plot_multi(blank_plot(paste("No data:", stem)), stem, w = 7.2, h = 3.8)
    return(invisible(FALSE))
  }
  p <- ggplot(dt, aes(x = band, y = odds_ratio)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
    geom_pointrange(aes(ymin = conf_low, ymax = conf_high), size = 0.55) +
    facet_wrap(~type, nrow = 1) +
    scale_y_log10() +
    theme_bw(base_size = 18) +
    labs(title = title, x = "Baseline rank band", y = "Odds ratio (band vs rest)") +
    theme(plot.title = element_text(face = "bold", size = 20), strip.text = element_text(face = "bold", size = 13), panel.grid.minor = element_blank())
  save_plot_multi(p, stem, w = 7.9, h = 3.55)
  invisible(TRUE)
}

plot_or(or_all, "Supp_enrichment_OR", "Motif enrichment (band vs rest)")
copy_stem("Supp_enrichment_OR", "Fig2F_enrichment_OR")
plot_or(or_ud, "Supp_enrichment_OR_up_vs_down", "Motif enrichment (Up vs Down)")
plot_or(or_bb, "Supp_enrichment_OR_by_band", "Motif enrichment (band vs rest)")

message("[08] wrote figures in: ", outdir)
