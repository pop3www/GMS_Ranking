#!/usr/bin/env Rscript
# Figure 1 production script - ported from 99_legacy/Main.R and synchronized to manuscript contract.
#
# Role: Figure 1 is descriptive. It establishes early near-parallel redistribution and
# a later head/shoulder deficit using transcriptome-wide log2(TPM+1) density summaries.
#
# Panels:
#   A - full matched-control delta-density profiles at 4 h and 24 h
#   B - head/shoulder zoom of the same delta-density profiles
#   C - condition-level raw density overlays, shoulder zoom + full distribution
#
# Important interpretation:
#   Ribbons are gene-resampling bootstrap bands around descriptive condition-level
#   transcriptome distributions. They are not biological-replicate confidence intervals.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(ggforce)
})

# -----------------------------
# Lightweight argument parser
# -----------------------------
script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", cmd, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = FALSE))
  normalizePath("figure1/scripts/run.R", mustWork = FALSE)
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    root = normalizePath(file.path(dirname(script_path()), "../.."), mustWork = FALSE),
    tpm = NA_character_,
    sample_sheet = NA_character_,
    contrasts = NA_character_,
    out_dir = NA_character_,
    source_dir = NA_character_,
    bootstrap_resamples = 5000L,
    seed = 123L,
    n_grid = 1024L,
    bandwidth_floor = 0.25,
    quick = FALSE
  )
  `%next%` <- function(dummy, i) if (i + 1L <= length(args)) args[[i + 1L]] else stop("Missing value for ", args[[i]])
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a == "--root") { out$root <- normalizePath(args %next% i, mustWork = FALSE); i <- i + 2L; next }
    if (a == "--tpm") { out$tpm <- args %next% i; i <- i + 2L; next }
    if (a == "--sample-sheet") { out$sample_sheet <- args %next% i; i <- i + 2L; next }
    if (a == "--contrasts") { out$contrasts <- args %next% i; i <- i + 2L; next }
    if (a == "--out-dir") { out$out_dir <- args %next% i; i <- i + 2L; next }
    if (a == "--source-dir") { out$source_dir <- args %next% i; i <- i + 2L; next }
    if (a == "--bootstrap-resamples") { out$bootstrap_resamples <- as.integer(args %next% i); i <- i + 2L; next }
    if (a == "--seed") { out$seed <- as.integer(args %next% i); i <- i + 2L; next }
    if (a == "--n-grid") { out$n_grid <- as.integer(args %next% i); i <- i + 2L; next }
    if (a == "--bandwidth-floor") { out$bandwidth_floor <- as.numeric(args %next% i); i <- i + 2L; next }
    if (a == "--quick") { out$quick <- TRUE; out$bootstrap_resamples <- min(out$bootstrap_resamples, 200L); i <- i + 1L; next }
    stop("Unknown argument: ", a)
  }
  out
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || is.na(a)) b else a
args <- parse_args()
root <- normalizePath(args$root, mustWork = FALSE)

fig_dir <- file.path(root, "figure1")
out_dir <- args$out_dir %||% file.path(fig_dir, "outputs")
src_dir <- args$source_dir %||% file.path(fig_dir, "source_data")
log_dir <- file.path(fig_dir, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

message("[Figure 1] repo root: ", root)
message("[Figure 1] bootstrap resamples: ", args$bootstrap_resamples)

# -----------------------------
# Inputs
# -----------------------------
find_first <- function(paths) {
  ok <- paths[file.exists(paths)]
  if (length(ok)) normalizePath(ok[[1]], mustWork = TRUE) else NA_character_
}

tpm_path <- if (!is.na(args$tpm)) args$tpm else find_first(c(
  file.path(root, "data/processed/TPMCountFile_rsemgenes.csv"),
  file.path(root, "data/processed/tpm_rsemgenes.tsv"),
  file.path(root, "data/processed/TPMCountFile_rsemgenes.tsv"),
  file.path(root, "data/raw/TPMCountFile_rsemgenes.csv"),
  file.path(root, "data/raw/TPMCountFile_rsemgenes.txt")
))
sample_sheet_path <- if (!is.na(args$sample_sheet)) args$sample_sheet else file.path(root, "config/sample_sheet.csv")
contrasts_path <- if (!is.na(args$contrasts)) args$contrasts else file.path(root, "config/contrasts.csv")

if (!file.exists(tpm_path)) stop("Missing TPM matrix: ", tpm_path)
if (!file.exists(sample_sheet_path)) stop("Missing sample sheet: ", sample_sheet_path)
if (!file.exists(contrasts_path)) stop("Missing contrasts file: ", contrasts_path)

read_expr <- function(path) {
  if (grepl("\\.(tsv|txt)$", tolower(path))) {
    x <- readr::read_tsv(path, show_col_types = FALSE)
  } else {
    x <- readr::read_csv(path, show_col_types = FALSE)
  }
  names(x)[1] <- "gene"
  # Defuse UTF-8 BOM if present.
  names(x)[1] <- "gene"
  x
}

sample_sheet <- readr::read_csv(sample_sheet_path, show_col_types = FALSE)
contrasts <- readr::read_csv(contrasts_path, show_col_types = FALSE)
expr_tpm <- read_expr(tpm_path)

canonical_condition <- function(group_label, time_h) {
  prefix <- as.character(time_h)
  dplyr::case_when(
    group_label == paste0(prefix, "_ctrl") ~ "Ctrl",
    group_label == paste0(prefix, "_Tam")  ~ "Tam",
    stringr::str_starts(group_label, prefix) ~ stringr::str_remove(group_label, paste0("^", prefix)),
    TRUE ~ as.character(group_label)
  )
}

sample_sheet <- sample_sheet %>%
  dplyr::mutate(
    time_h = as.integer(time_h),
    condition = canonical_condition(group_label, time_h)
  )

expected_conditions <- c("Ctrl", "Tam", "D", "DT", "D_L_CPT", "DT_L_CPT", "D_H_CPT", "DT_H_CPT")
unknown_conditions <- setdiff(unique(sample_sheet$condition), expected_conditions)
if (length(unknown_conditions)) {
  warning("Unexpected condition labels in sample sheet: ", paste(unknown_conditions, collapse = ", "))
}

sample_cols <- sample_sheet$sample_id
missing_samples <- setdiff(sample_cols, names(expr_tpm))
if (length(missing_samples)) stop("TPM matrix is missing sample columns: ", paste(missing_samples, collapse = ", "))

expr_tpm <- expr_tpm %>%
  dplyr::select(gene, tidyselect::all_of(sample_cols)) %>%
  dplyr::mutate(dplyr::across(-gene, ~ suppressWarnings(as.numeric(.x))))

# -----------------------------
# Condition-level mean TPM, then log2(mean TPM + 1)
# This matches the legacy Main.R figure path, which used *_mean columns before taking log2.
# -----------------------------
mean_tbl <- tibble::tibble(gene = expr_tpm$gene)
for (idx in seq_len(nrow(sample_sheet %>% dplyr::distinct(time_h, condition)))) {
  row <- sample_sheet %>% dplyr::distinct(time_h, condition) %>% dplyr::slice(idx)
  cols <- sample_sheet %>% dplyr::filter(time_h == row$time_h, condition == row$condition) %>% dplyr::pull(sample_id)
  cname <- paste0(row$time_h, "h__", row$condition)
  mean_tbl[[cname]] <- rowMeans(as.matrix(expr_tpm[, cols, drop = FALSE]), na.rm = TRUE)
}

long_mean <- mean_tbl %>%
  tidyr::pivot_longer(-gene, names_to = "profile", values_to = "TPM_mean") %>%
  tidyr::separate(profile, into = c("timepoint", "condition"), sep = "__", remove = FALSE) %>%
  dplyr::mutate(
    time_h = as.integer(stringr::str_remove(timepoint, "h$")),
    timepoint = factor(paste0(time_h, "h"), levels = c("4h", "24h")),
    condition = factor(condition, levels = expected_conditions),
    logTPM = log2(TPM_mean + 1)
  ) %>%
  dplyr::filter(is.finite(logTPM), !is.na(condition))

readr::write_tsv(long_mean %>% dplyr::select(gene, time_h, condition, TPM_mean, logTPM),
                 file.path(src_dir, "figure1_condition_mean_log2TPM.tsv"))

# -----------------------------
# Shared density grid/bandwidth
# -----------------------------
x_all <- long_mean$logTPM[is.finite(long_mean$logTPM)]
grid_lo <- floor(min(x_all))
grid_hi <- ceiling(max(x_all))
NGRID <- as.integer(args$n_grid)
grid_x <- seq(grid_lo, grid_hi, length.out = NGRID)

bw0 <- tryCatch(stats::bw.SJ(x_all), error = function(e) NA_real_)
if (!is.finite(bw0) || bw0 <= 0) bw0 <- stats::bw.nrd0(x_all)
bw0 <- max(bw0, as.numeric(args$bandwidth_floor))
message(sprintf("[Figure 1] density grid %.2f..%.2f n=%d bw=%.4f", grid_lo, grid_hi, NGRID, bw0))

density_y <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 5) return(rep(NA_real_, NGRID))
  stats::density(v, from = grid_lo, to = grid_hi, n = NGRID, bw = bw0)$y
}

# Raw condition-level densities.
dens_df <- long_mean %>%
  dplyr::group_by(time_h, timepoint, condition) %>%
  dplyr::summarise(y = list(density_y(logTPM)), .groups = "drop") %>%
  dplyr::mutate(i = list(seq_len(NGRID))) %>%
  tidyr::unnest(c(i, y)) %>%
  dplyr::mutate(x = grid_x[i]) %>%
  dplyr::arrange(time_h, condition, i)

readr::write_tsv(dens_df, file.path(src_dir, "figure1_raw_density_curves.tsv"))

# Matched-control delta density per config/contrasts.csv.
get_profile <- function(tp_h, cond) {
  out <- long_mean %>% dplyr::filter(time_h == tp_h, as.character(condition) == cond) %>% dplyr::arrange(gene)
  if (!nrow(out)) stop("No profile for time_h=", tp_h, " condition=", cond)
  out
}

calc_delta_one <- function(contrast_id, numerator, denominator, time_h) {
  a <- get_profile(time_h, numerator)
  b <- get_profile(time_h, denominator)
  common <- intersect(a$gene, b$gene)
  a <- a %>% dplyr::filter(gene %in% common) %>% dplyr::arrange(gene)
  b <- b %>% dplyr::filter(gene %in% common) %>% dplyr::arrange(gene)
  ya <- density_y(a$logTPM)
  yb <- density_y(b$logTPM)
  tibble::tibble(
    contrast_id = contrast_id,
    time_h = as.integer(time_h),
    timepoint = factor(paste0(time_h, "h"), levels = c("4h", "24h")),
    numerator = numerator,
    denominator = denominator,
    i = seq_len(NGRID),
    x = grid_x,
    density_numerator = ya,
    density_denominator = yb,
    delta = ya - yb
  )
}

delta_tc <- purrr::pmap_dfr(
  list(contrasts$contrast_id, contrasts$numerator, contrasts$denominator, contrasts$time_h),
  calc_delta_one
)

# Gene bootstrap bands for delta densities. Same resampled gene indices are applied to
# numerator and denominator to keep the contrast/baseline pair matched at the gene-resampling step.
set.seed(args$seed)
calc_ci_matched <- function(tp_h, numerator, denominator, contrast_id, B = args$bootstrap_resamples) {
  a <- get_profile(tp_h, numerator)
  b <- get_profile(tp_h, denominator)
  common <- intersect(a$gene, b$gene)
  a <- a %>% dplyr::filter(gene %in% common) %>% dplyr::arrange(gene)
  b <- b %>% dplyr::filter(gene %in% common) %>% dplyr::arrange(gene)
  xA <- a$logTPM; xB <- b$logTPM
  n <- length(xA)
  if (n < 5 || B <= 0) {
    return(tibble::tibble(contrast_id = contrast_id, i = seq_len(NGRID), x = grid_x,
                          ymin = NA_real_, ymax = NA_real_))
  }
  mat <- replicate(B, {
    idx <- sample.int(n, size = n, replace = TRUE)
    density_y(xA[idx]) - density_y(xB[idx])
  })
  tibble::tibble(
    contrast_id = contrast_id,
    time_h = as.integer(tp_h),
    timepoint = factor(paste0(tp_h, "h"), levels = c("4h", "24h")),
    numerator = numerator,
    denominator = denominator,
    i = seq_len(NGRID),
    x = grid_x,
    ymin = apply(mat, 1, stats::quantile, 0.025, na.rm = TRUE),
    ymax = apply(mat, 1, stats::quantile, 0.975, na.rm = TRUE),
    bootstrap_unit = "gene",
    bootstrap_pairing = "same resampled gene indices for numerator and denominator",
    n_bootstrap = as.integer(B)
  )
}

boot_delta_tc <- purrr::pmap_dfr(
  list(contrasts$time_h, contrasts$numerator, contrasts$denominator, contrasts$contrast_id),
  calc_ci_matched
)

readr::write_tsv(delta_tc, file.path(src_dir, "figure1_delta_density_curves.tsv"))
readr::write_tsv(boot_delta_tc, file.path(src_dir, "figure1_delta_density_gene_bootstrap.tsv"))

plot_delta <- delta_tc %>%
  dplyr::left_join(boot_delta_tc %>% dplyr::select(contrast_id, i, ymin, ymax), by = c("contrast_id", "i"))

# -----------------------------
# Plot styling
# -----------------------------
okabe_ito <- c(
  Ctrl = "#000000",
  Tam = "#E69F00",
  D = "#56B4E9",
  DT = "#009E73",
  D_L_CPT = "#F0E442",
  DT_L_CPT = "#0072B2",
  D_H_CPT = "#D55E00",
  DT_H_CPT = "#CC79A7"
)

# Show delta curves by numerator condition; source data carries denominator explicitly.
plot_delta <- plot_delta %>%
  dplyr::mutate(numerator = factor(numerator, levels = names(okabe_ito)),
                timepoint = factor(timepoint, levels = c("4h", "24h")))

dens_df <- dens_df %>%
  dplyr::mutate(condition = factor(condition, levels = names(okabe_ito)),
                timepoint = factor(timepoint, levels = c("4h", "24h")))

base_theme <- ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 16),
    axis.text = ggplot2::element_text(size = 13, colour = "black"),
    strip.text = ggplot2::element_text(size = 14, face = "bold"),
    legend.title = ggplot2::element_text(size = 12),
    legend.text = ggplot2::element_text(size = 11),
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(size = 14, face = "bold")
  )

p_a <- ggplot2::ggplot(plot_delta, ggplot2::aes(x = x, y = delta, colour = numerator)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.55, colour = "grey35") +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax, fill = numerator),
                       alpha = 0.08, colour = NA, show.legend = FALSE) +
  ggplot2::geom_line(linewidth = 0.75) +
  ggplot2::scale_colour_manual(values = okabe_ito, drop = FALSE, na.translate = FALSE) +
  ggplot2::scale_fill_manual(values = okabe_ito, drop = FALSE, na.translate = FALSE) +
  ggplot2::facet_grid(. ~ timepoint, scales = "fixed") +
  ggplot2::labs(x = "log2(TPM + 1)", y = "Δ density", colour = "Condition / numerator") +
  base_theme

x_zoom <- c(5, 15)
p_b <- ggplot2::ggplot(plot_delta, ggplot2::aes(x = x, y = delta, colour = numerator)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.55, colour = "grey35") +
  ggplot2::geom_ribbon(data = plot_delta %>% dplyr::filter(x >= x_zoom[1], x <= x_zoom[2]),
                       ggplot2::aes(ymin = ymin, ymax = ymax, fill = numerator),
                       alpha = 0.11, colour = NA, show.legend = FALSE) +
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::coord_cartesian(xlim = x_zoom) +
  ggplot2::scale_colour_manual(values = okabe_ito, drop = FALSE, na.translate = FALSE) +
  ggplot2::scale_fill_manual(values = okabe_ito, drop = FALSE, na.translate = FALSE) +
  ggplot2::facet_grid(. ~ timepoint, scales = "fixed") +
  ggplot2::labs(x = "log2(TPM + 1)", y = "Δ density", colour = "Condition / numerator") +
  base_theme

p_c <- ggplot2::ggplot(dens_df, ggplot2::aes(x = x, y = y, colour = condition, linetype = timepoint)) +
  ggplot2::geom_line(linewidth = 0.70, alpha = 0.95) +
  ggplot2::scale_colour_manual(values = okabe_ito, drop = FALSE, na.translate = FALSE) +
  ggplot2::scale_linetype_manual(values = c("4h" = "22", "24h" = "solid")) +
  ggforce::facet_zoom(xlim = x_zoom, ylim = c(0, 0.10), zoom.size = 0.70) +
  ggplot2::labs(x = "log2(TPM + 1)", y = "Density", colour = "Condition / numerator", linetype = "Timepoint") +
  base_theme

fig1 <- ggpubr::ggarrange(
  p_a, p_b, p_c,
  ncol = 1,
  heights = c(1.18, 1.18, 1.06),
  common.legend = TRUE,
  legend = "right",
  labels = c("A", "B", "C"),
  font.label = list(size = 18, face = "bold")
)

save_plot <- function(plot, stem, width = 8.2, height = 8.8) {
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in")
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 600, units = "in")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(file.path(out_dir, paste0(stem, ".svg")), plot, width = width, height = height, units = "in", device = svglite::svglite)
  } else {
    ggplot2::ggsave(file.path(out_dir, paste0(stem, ".svg")), plot, width = width, height = height, units = "in", device = grDevices::svg)
  }
}

save_plot(fig1, "figure1_main")
writeLines("ok", file.path(out_dir, "figure1.done"))

manifest_lines <- c(
  paste0("script=", script_path()),
  paste0("input_tpm=", normalizePath(tpm_path, mustWork = FALSE)),
  paste0("input_sample_sheet=", normalizePath(sample_sheet_path, mustWork = FALSE)),
  paste0("input_contrasts=", normalizePath(contrasts_path, mustWork = FALSE)),
  "expression_scale=log2(condition mean TPM + 1)",
  "delta_logic=matched contrast numerator density minus denominator density from config/contrasts.csv",
  paste0("density_bandwidth=", bw0),
  paste0("density_grid=", grid_lo, ",", grid_hi, ",", NGRID),
  paste0("bootstrap_resamples=", args$bootstrap_resamples),
  "bootstrap_interpretation=gene-resampling descriptive bands; not biological-replicate CIs",
  "panel_labels=uppercase A/B/C",
  "panel_C_linetype=4h dashed; 24h solid"
)
writeLines(manifest_lines, file.path(src_dir, "figure1_manifest.txt"))
writeLines(manifest_lines, file.path(log_dir, "run.log"))

message("[Figure 1] wrote outputs to: ", normalizePath(out_dir, mustWork = FALSE))
message("[Figure 1] wrote source data to: ", normalizePath(src_dir, mustWork = FALSE))
