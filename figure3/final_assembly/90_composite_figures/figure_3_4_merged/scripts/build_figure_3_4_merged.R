#!/usr/bin/env Rscript

# Build a journal-layout merged Figure 3 from existing Figure 3 and Figure 4
# source-data tables. This script is intentionally read-only with respect to
# the original figure directories. It writes only under
# 90_composite_figures/figure_3_4_merged/ by default.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggpubr)
})

parse_args <- function(args) {
  out <- list(
    root = ".",
    fig3_source = "figure3/differential_expression/source_data",
    fig4_source = "figure3/scalar_metrics/source_data",
    out_dir = "90_composite_figures/figure_3_4_merged/outputs",
    source_dir = "90_composite_figures/figure_3_4_merged/source_data",
    log_dir = "90_composite_figures/figure_3_4_merged/logs",
    representative_contrast = "DT_vs_D_24h",
    early_contrast = "DT_vs_D_4h",
    max_points_per_method = 2500,
    make_s5 = TRUE,
    width = 12,
    height = 12,
    dpi = 500
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    nm <- sub("^--", "", key)
    if (nm %in% c("make-s5")) nm <- "make_s5"
    if (nm %in% c("out-dir")) nm <- "out_dir"
    if (nm %in% c("source-dir")) nm <- "source_dir"
    if (nm %in% c("log-dir")) nm <- "log_dir"
    if (nm %in% c("fig3-source")) nm <- "fig3_source"
    if (nm %in% c("fig4-source")) nm <- "fig4_source"
    if (nm %in% c("representative-contrast")) nm <- "representative_contrast"
    if (nm %in% c("early-contrast")) nm <- "early_contrast"
    if (nm %in% c("max-points-per-method")) nm <- "max_points_per_method"
    if (!(nm %in% names(out))) stop("Unknown option: ", key)
    if (i == length(args)) stop("Missing value for ", key)
    val <- args[[i + 1]]
    if (nm %in% c("max_points_per_method", "width", "height", "dpi")) {
      val <- as.numeric(val)
    }
    if (nm %in% c("make_s5")) {
      val <- tolower(val) %in% c("true", "t", "1", "yes", "y")
    }
    out[[nm]] <- val
    i <- i + 2
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
root <- normalizePath(args$root, mustWork = TRUE)

path_in <- function(...) file.path(root, ...)
fig3_source <- path_in(args$fig3_source)
fig4_source <- path_in(args$fig4_source)
out_dir <- path_in(args$out_dir)
source_dir <- path_in(args$source_dir)
log_dir <- path_in(args$log_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

METHOD_ORDER <- c("DESeq2", "edgeR", "limma/voom", "RankProd", "PenDA", "RankCompV3")
DIR_COLS <- c("Down" = "#3B78B8", "NS" = "grey82", "Up" = "#D95F02")
TIME_COLS <- c("4 h" = "#3B78B8", "24 h" = "#D95F02")
HEAD_COLS <- c("Non-head" = "grey78", "Head" = "#7B3294")
BASE_THEME <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(5, 6, 5, 6)
  )

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required input not found: ", path,
         "\nRun the existing Figure 3/Figure 4 production scripts first. This composite builder does not regenerate their source data.")
  }
  path
}

clean_label <- function(x) {
  x %>%
    str_replace_all("_vs_", " vs ") %>%
    str_replace_all("_", " ") %>%
    str_replace(" 4h$", " (4 h)") %>%
    str_replace(" 24h$", " (24 h)")
}

get_time_label <- function(contrast_id, time_h = NULL) {
  if (!is.null(time_h)) {
    out <- ifelse(as.numeric(time_h) == 4, "4 h", ifelse(as.numeric(time_h) == 24, "24 h", paste0(time_h, " h")))
    return(out)
  }
  ifelse(str_detect(contrast_id, "_4h$"), "4 h",
         ifelse(str_detect(contrast_id, "_24h$"), "24 h", NA_character_))
}

pick_col <- function(df, candidates, required = TRUE, label = "column") {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    if (required) stop("Could not find ", label, ". Candidates: ", paste(candidates, collapse = ", "),
                       "; available: ", paste(names(df), collapse = ", "))
    return(NA_character_)
  }
  hit[[1]]
}

# -------------------------------------------------------------------------
# Read finalized Figure 3/4 source data.
# -------------------------------------------------------------------------
used_path <- require_file(file.path(fig3_source, "figure3_unified_de_results_used.tsv"))
status_path <- require_file(file.path(fig3_source, "unified_de_method_status.tsv"))
metrics_path <- require_file(file.path(fig4_source, "figure4_scalar_metrics.tsv"))

used <- read_tsv(used_path, show_col_types = FALSE) %>%
  mutate(method = as.character(method), contrast_id = as.character(contrast_id))
status <- read_tsv(status_path, show_col_types = FALSE)
metrics <- read_tsv(metrics_path, show_col_types = FALSE)

if (!all(METHOD_ORDER %in% unique(used$method))) {
  stop("Figure 3 used-results table does not contain the expected methods: ",
       paste(setdiff(METHOD_ORDER, unique(used$method)), collapse = ", "))
}
if (any(status$status != "ok")) {
  stop("Figure 3 method-status table contains non-ok status values. Inspect: ", status_path)
}
if ("CellComp" %in% unique(used$method)) {
  stop("CellComp appears in active used-results table. Current composite requires RankCompV3 instead.")
}
if (!("RankCompV3" %in% unique(used$method))) {
  stop("RankCompV3 is missing from active used-results table.")
}
if (!("score_for_agreement" %in% names(used))) {
  stop("score_for_agreement column missing from Figure 3 used-results. Run the score-normalized Figure 3 code first.")
}

# -------------------------------------------------------------------------
# Panel A: DE call counts by method, representative contrast.
# -------------------------------------------------------------------------
panelA_path <- file.path(fig3_source, "figure3_panelA_call_counts.tsv")
if (file.exists(panelA_path)) {
  panelA_counts <- read_tsv(panelA_path, show_col_types = FALSE)
} else {
  panelA_counts <- used %>%
    filter(contrast_id == args$representative_contrast, direction %in% c("Up", "Down")) %>%
    count(method, direction, name = "n") %>%
    group_by(method) %>%
    mutate(total = sum(n)) %>%
    ungroup()
}

panelA_counts <- panelA_counts %>%
  filter(method %in% METHOD_ORDER, direction %in% c("Up", "Down")) %>%
  mutate(
    method = factor(method, levels = rev(METHOD_ORDER)),
    direction = factor(direction, levels = c("Down", "Up"))
  )
panelA_totals <- panelA_counts %>%
  group_by(method) %>%
  summarise(total = sum(n), .groups = "drop")

pA <- ggplot(panelA_counts, aes(x = n, y = method, fill = direction)) +
  geom_col(width = 0.72) +
  geom_text(data = panelA_totals, aes(x = total, y = method, label = total),
            inherit.aes = FALSE, hjust = -0.15, size = 3.0) +
  scale_fill_manual(values = DIR_COLS[c("Down", "Up")], name = "Direction") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = "DE calls by method", x = "Number of genes", y = NULL) +
  BASE_THEME +
  theme(legend.position = "bottom")

# -------------------------------------------------------------------------
# Panel B: representative late contrast signed-score maps.
# -------------------------------------------------------------------------
set.seed(10)
panelB <- used %>%
  filter(contrast_id == args$representative_contrast, method %in% METHOD_ORDER) %>%
  mutate(
    method = factor(method, levels = METHOD_ORDER),
    direction = factor(direction, levels = c("Down", "NS", "Up"))
  )
if (nrow(panelB) == 0) stop("No rows for representative contrast: ", args$representative_contrast)

panelB_plot <- bind_rows(lapply(split(panelB, panelB$method), function(x) {
  if (nrow(x) > args$max_points_per_method) x[sample(seq_len(nrow(x)), args$max_points_per_method), , drop = FALSE] else x
}))

pB <- ggplot(panelB_plot, aes(x = baseline_expression, y = score_for_agreement, color = direction)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey50") +
  geom_point(size = 0.25, alpha = 0.35, na.rm = TRUE) +
  facet_wrap(~method, ncol = 3) +
  scale_color_manual(values = DIR_COLS, drop = FALSE, name = "Call") +
  labs(title = paste0("Representative late contrast: ", clean_label(args$representative_contrast)),
       x = "Matched-control baseline expression", y = "Rank-normalized signed method score") +
  BASE_THEME +
  theme(legend.position = "bottom")

# -------------------------------------------------------------------------
# Panel C: inter-method spread by early/late and head/non-head.
# Direct definition: IQR(score_for_agreement across methods) per gene.
# -------------------------------------------------------------------------
panelC <- used %>%
  filter(contrast_id %in% c(args$early_contrast, args$representative_contrast), method %in% METHOD_ORDER) %>%
  group_by(contrast_id, gene_id) %>%
  summarise(
    spread_iqr = IQR(score_for_agreement, na.rm = TRUE),
    baseline_expression = median(baseline_expression, na.rm = TRUE),
    n_methods = n_distinct(method),
    .groups = "drop"
  ) %>%
  filter(n_methods == length(METHOD_ORDER), is.finite(spread_iqr), is.finite(baseline_expression)) %>%
  group_by(contrast_id) %>%
  mutate(
    head_cutoff = quantile(baseline_expression, 0.90, na.rm = TRUE),
    rank_bin = ifelse(baseline_expression >= head_cutoff, "Head", "Non-head")
  ) %>%
  ungroup() %>%
  mutate(
    time_label = factor(get_time_label(contrast_id), levels = c("4 h", "24 h")),
    rank_bin = factor(rank_bin, levels = c("Non-head", "Head"))
  )

pC <- ggplot(panelC, aes(x = time_label, y = spread_iqr, fill = rank_bin)) +
  geom_boxplot(width = 0.65, outlier.alpha = 0.05, outlier.size = 0.25, linewidth = 0.3) +
  scale_fill_manual(values = HEAD_COLS, name = "Baseline rank") +
  labs(title = "Inter-method signed-score spread localizes to the expression head",
       x = NULL, y = "Inter-method score spread (IQR)") +
  BASE_THEME +
  theme(legend.position = "bottom")

# -------------------------------------------------------------------------
# Panels D/E: scalar summaries of directionality and proportionality.
# -------------------------------------------------------------------------
metrics2 <- metrics %>%
  mutate(
    time_label = if ("time_label" %in% names(.)) as.character(time_label) else get_time_label(contrast_id, time_h),
    time_label = factor(time_label, levels = c("4 h", "24 h")),
    contrast_label = if ("contrast_label" %in% names(.)) as.character(contrast_label) else clean_label(contrast_id),
    contrast_label = factor(contrast_label, levels = rev(unique(contrast_label)))
  )
if (!all(c("Up", "Down", "UpFrac", "slope_rlm") %in% names(metrics2))) {
  stop("Scalar metrics table is missing one of: Up, Down, UpFrac, slope_rlm")
}

pD <- ggplot(metrics2, aes(x = UpFrac, y = contrast_label, fill = time_label)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = TIME_COLS, name = "Time point") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Up fraction", x = "Up / (Up + Down)", y = NULL) +
  BASE_THEME +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 8))

pE <- ggplot(metrics2, aes(x = slope_rlm, y = contrast_label, fill = time_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey45") +
  geom_col(width = 0.72) +
  scale_fill_manual(values = TIME_COLS, name = "Time point") +
  labs(title = "Baseline-expression slope, beta", x = "Robust slope, beta", y = NULL) +
  BASE_THEME +
  theme(legend.position = "bottom", axis.text.y = element_blank(), axis.ticks.y = element_blank())

# -------------------------------------------------------------------------
# Assemble merged main Figure 3.
# -------------------------------------------------------------------------
row1 <- ggarrange(pA, pB, ncol = 2, widths = c(0.95, 2.05), labels = c("A", "B"),
                  font.label = list(size = 16, face = "bold"), align = "hv")
row2 <- ggarrange(pC, ncol = 1, labels = c("C"), font.label = list(size = 16, face = "bold"), align = "hv")
row3 <- ggarrange(pD, pE, ncol = 2, widths = c(1.15, 1), labels = c("D", "E"),
                  font.label = list(size = 16, face = "bold"), align = "hv")
fig_main <- ggarrange(row1, row2, row3, ncol = 1, heights = c(1.35, 0.95, 1.35), align = "v")

for (ext in c("pdf", "svg", "png")) {
  outfile <- file.path(out_dir, paste0("figure_3_4_merged.", ext))
  if (ext == "png") {
    ggsave(outfile, fig_main, width = args$width, height = args$height, dpi = args$dpi, bg = "white")
  } else if (ext == "svg") {
    ggsave(outfile, fig_main, width = args$width, height = args$height, device = svglite::svglite, bg = "white")
  } else {
    ggsave(outfile, fig_main, width = args$width, height = args$height, bg = "white")
  }
}

write_tsv(panelA_counts %>% mutate(method = as.character(method), direction = as.character(direction)),
          file.path(source_dir, "figure_3_4_panelA_call_counts.tsv"))
write_tsv(panelB %>% mutate(method = as.character(method), direction = as.character(direction)),
          file.path(source_dir, "figure_3_4_panelB_used_points_all_genes.tsv"))
write_tsv(panelC, file.path(source_dir, "figure_3_4_panelC_effect_spread.tsv"))
write_tsv(metrics2 %>% mutate(contrast_label = as.character(contrast_label), time_label = as.character(time_label)),
          file.path(source_dir, "figure_3_4_panelD_E_scalar_metrics.tsv"))

manifest <- c(
  paste0("Created: ", Sys.time()),
  "Composite figure: manuscript Figure 3 assembled from differential-expression and scalar-metric source data.",
  "This script is read-only with respect to figure3/differential_expression and figure3/scalar_metrics.",
  paste0("Representative contrast: ", args$representative_contrast),
  paste0("Early contrast for panel C: ", args$early_contrast),
  paste0("Differential-expression source: ", fig3_source),
  paste0("Scalar-metrics source: ", fig4_source),
  "Panel A: differential-expression call counts by method.",
  "Panel B: method-specific rank-normalized signed scores for the representative late contrast.",
  "Panel C: per-gene IQR of rank-normalized signed method scores, stratified by baseline head/non-head.",
  "Panel D: Up fraction among limma/voom DE genes.",
  "Panel E: baseline-expression slope, beta, using slope_rlm.",
  "Method set: DESeq2, edgeR, limma/voom, RankProd, PenDA, RankCompV3 (REO).",
  "RankCompV3 (REO) is the REO-family method used in the current comparison.",
  "Component analysis outputs are read without modification by this assembler."
)
writeLines(manifest, file.path(source_dir, "figure_3_4_merged_manifest.txt"))

# -------------------------------------------------------------------------
# Supplementary Figure S5: pairwise method-concordance panels.
# -------------------------------------------------------------------------
make_s5 <- isTRUE(args$make_s5)
if (make_s5) {
  s5_out <- path_in("90_composite_figures/figure_3_4_merged/outputs/supplementary_fig_s5_method_concordance")
  s5_src <- path_in("90_composite_figures/figure_3_4_merged/source_data/supplementary_fig_s5_method_concordance")
  dir.create(s5_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(s5_src, recursive = TRUE, showWarnings = FALSE)

  jacc_path <- file.path(fig3_source, "figure3_panelC_jaccard.tsv")
  spear_path <- file.path(fig3_source, "figure3_panelD_spearman.tsv")

  if (file.exists(jacc_path) && file.exists(spear_path)) {
    jacc0 <- read_tsv(jacc_path, show_col_types = FALSE)
    spear0 <- read_tsv(spear_path, show_col_types = FALSE)

    m1 <- pick_col(jacc0, c("method1", "method_1", "method_a", "Var1"), label = "Jaccard method1")
    m2 <- pick_col(jacc0, c("method2", "method_2", "method_b", "Var2"), label = "Jaccard method2")
    typ <- pick_col(jacc0, c("set_type", "kind", "call_set", "direction_set", "class"), label = "Jaccard set type")
    val <- pick_col(jacc0, c("jaccard", "Jaccard", "similarity", "value", "jaccard_value"), label = "Jaccard value")
    jacc <- tibble(
      method1 = as.character(jacc0[[m1]]),
      method2 = as.character(jacc0[[m2]]),
      set_type = as.character(jacc0[[typ]]),
      jaccard = as.numeric(jacc0[[val]])
    ) %>%
      mutate(
        method1 = factor(method1, levels = METHOD_ORDER),
        method2 = factor(method2, levels = rev(METHOD_ORDER)),
        label = case_when(
          is.na(jaccard) ~ "",
          jaccard > 0 & jaccard < 0.01 ~ "<0.01",
          TRUE ~ sprintf("%.2f", jaccard)
        )
      )

    s1 <- pick_col(spear0, c("method1", "method_1", "method_a", "Var1"), label = "Spearman method1")
    s2 <- pick_col(spear0, c("method2", "method_2", "method_b", "Var2"), label = "Spearman method2")
    sv <- pick_col(spear0, c("rho", "spearman", "cor", "value"), label = "Spearman value")
    spear <- tibble(
      method1 = as.character(spear0[[s1]]),
      method2 = as.character(spear0[[s2]]),
      rho = as.numeric(spear0[[sv]])
    ) %>%
      mutate(
        method1 = factor(method1, levels = METHOD_ORDER),
        method2 = factor(method2, levels = rev(METHOD_ORDER)),
        label = sprintf("%.2f", rho)
      )

    pS5A <- ggplot(jacc, aes(x = method1, y = method2, fill = jaccard)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = label), size = 2.7) +
      facet_wrap(~set_type) +
      scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1), na.value = "grey90") +
      labs(title = "Pairwise Jaccard overlap", x = NULL, y = NULL, fill = "Jaccard") +
      BASE_THEME +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

    pS5B <- ggplot(spear, aes(x = method1, y = method2, fill = rho)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = label, color = abs(rho) > 0.65), size = 3) +
      scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
      scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0,
                           limits = c(-1, 1), na.value = "grey90") +
      labs(title = "Pairwise Spearman concordance", x = NULL, y = NULL, fill = "Spearman \u03c1") +
      BASE_THEME +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

    fig_s5 <- ggarrange(pS5A, pS5B, ncol = 1, labels = c("A", "B"),
                        font.label = list(size = 16, face = "bold"), heights = c(1.15, 1))
    for (ext in c("pdf", "svg", "png")) {
      outfile <- file.path(s5_out, paste0("figure_s5_method_concordance.", ext))
      if (ext == "png") {
        ggsave(outfile, fig_s5, width = 10, height = 10, dpi = args$dpi, bg = "white")
      } else if (ext == "svg") {
        ggsave(outfile, fig_s5, width = 10, height = 10, device = svglite::svglite, bg = "white")
      } else {
        ggsave(outfile, fig_s5, width = 10, height = 10, bg = "white")
      }
    }
    write_tsv(jacc %>% mutate(method1 = as.character(method1), method2 = as.character(method2)),
              file.path(s5_src, "figure_s5_panelA_jaccard.tsv"))
    write_tsv(spear %>% mutate(method1 = as.character(method1), method2 = as.character(method2)),
              file.path(s5_src, "figure_s5_panelB_spearman.tsv"))
    writeLines(c(
      paste0("Created: ", Sys.time()),
      "Supplementary Figure S5 contains the pairwise Jaccard and Spearman concordance panels.",
      "Panel A: pairwise Jaccard overlap of thresholded call sets.",
      "Panel B: pairwise Spearman concordance of rank-normalized signed method scores.",
      "Small nonzero Jaccard values below 0.01 are displayed as <0.01."
    ), file.path(s5_src, "figure_s5_method_concordance_manifest.txt"))
  } else {
    warning("Skipping S5 because Jaccard or Spearman source table is missing.")
  }
}

cat("Wrote merged Figure 3 outputs to: ", out_dir, "\n", sep = "")
cat("Wrote merged Figure 3 source data to: ", source_dir, "\n", sep = "")
