#!/usr/bin/env Rscript
# Figure 3 production plotter
# Builds Figure 3 A-E from a unified DE result table.
# The DE table may be generated locally by build_unified_de_results.R or provided as a precomputed result table.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggpubr)
})

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      if (i == length(args) || startsWith(args[[i + 1]], "--")) {
        out[[key]] <- TRUE
        i <- i + 1
      } else {
        out[[key]] <- args[[i + 1]]
        i <- i + 2
      }
    } else {
      i <- i + 1
    }
  }
  out
}

repo_root <- function() {
  p <- normalizePath(getwd(), mustWork = FALSE)
  # If run from scripts directory, climb to repo root.
  if (basename(p) == "scripts") p <- normalizePath(file.path(p, "../../.."), mustWork = FALSE)
  p
}

coalesce_arg <- function(x, default) if (is.null(x) || isTRUE(x) || !nzchar(as.character(x))) default else as.character(x)

read_table_auto <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) readr::read_tsv(path, show_col_types = FALSE)
  else readr::read_csv(path, show_col_types = FALSE)
}

first_present <- function(nm, candidates) {
  hit <- candidates[candidates %in% nm]
  if (length(hit)) hit[[1]] else NA_character_
}

canonical_method <- function(x) {
  lx <- tolower(as.character(x))
  dplyr::case_when(
    stringr::str_detect(lx, "deseq2") ~ "DESeq2",
    stringr::str_detect(lx, "edger") ~ "edgeR",
    stringr::str_detect(lx, "limma|voom") ~ "limma/voom",
    stringr::str_detect(lx, "rankprod") ~ "RankProd",
    stringr::str_detect(lx, "penda") ~ "PenDA",
    stringr::str_detect(lx, "cellcomp") ~ "CellComp",
    stringr::str_detect(lx, "rankcomp") ~ "RankCompV3",
    TRUE ~ as.character(x)
  )
}

condition_from_sample_sheet <- function(sample_sheet) {
  sample_sheet %>%
    mutate(
      time_h = as.integer(time_h),
      condition = case_when(
        group_label %in% c("4_ctrl", "24_ctrl") ~ "Ctrl",
        group_label %in% c("4_Tam", "24_Tam") ~ "Tam",
        TRUE ~ str_remove(group_label, paste0("^", time_h))
      )
    )
}

samples_for <- function(sample_sheet, condition, time_h) {
  sample_sheet %>%
    filter(.data$condition == condition, .data$time_h == as.integer(time_h)) %>%
    arrange(.data$replicate) %>%
    pull(.data$sample_id)
}

read_expression_matrix <- function(path) {
  x <- read_table_auto(path)
  gcol <- first_present(names(x), c("gene_id", "gene", "name", "Gene", "Geneid"))
  if (is.na(gcol)) stop("Expression matrix needs a gene_id/gene/name column: ", path)
  x <- x %>% rename(gene_id = all_of(gcol))
  mat <- as.data.frame(x)
  rownames(mat) <- make.unique(as.character(mat$gene_id))
  mat$gene_id <- NULL
  mat[] <- lapply(mat, function(v) suppressWarnings(as.numeric(v)))
  as.matrix(mat)
}

median_ratio_log_expression <- function(counts) {
  counts <- as.matrix(counts)
  counts[counts < 0] <- NA_real_
  keep <- rowSums(is.finite(counts) & counts > 0) >= 2
  gm <- rep(NA_real_, nrow(counts))
  gm[keep] <- exp(rowMeans(log(pmax(counts[keep, , drop = FALSE], 1e-8)), na.rm = TRUE))
  ratios <- sweep(counts[is.finite(gm) & gm > 0, , drop = FALSE], 1, gm[is.finite(gm) & gm > 0], "/")
  sf <- apply(ratios, 2, function(z) median(z[is.finite(z) & z > 0], na.rm = TRUE))
  sf[!is.finite(sf) | sf <= 0] <- 1
  sf <- sf / exp(mean(log(sf)))
  log2(sweep(counts, 2, sf, "/") + 1)
}

baseline_expression_for_contrasts <- function(counts_path, sample_sheet_path, contrasts_path, contrast_ids) {
  if (is.null(counts_path) || isTRUE(counts_path) || !file.exists(counts_path)) {
    stop("Panel B/E baseline expression requires --counts when baseline_expression/Amean is absent in the DE table.")
  }
  counts <- read_expression_matrix(counts_path)
  sample_sheet <- condition_from_sample_sheet(readr::read_csv(sample_sheet_path, show_col_types = FALSE))
  contrasts <- readr::read_csv(contrasts_path, show_col_types = FALSE)
  log_expr <- median_ratio_log_expression(counts)
  purrr::map_dfr(contrast_ids, function(cid) {
    cc <- contrasts %>% filter(.data$contrast_id == cid)
    if (nrow(cc) != 1) stop("Cannot find contrast in contrasts.csv: ", cid)
    denom_samples <- samples_for(sample_sheet, cc$denominator[[1]], cc$time_h[[1]])
    denom_samples <- intersect(denom_samples, colnames(log_expr))
    if (!length(denom_samples)) stop("No denominator samples found for ", cid)
    tibble(
      contrast_id = cid,
      gene_id = rownames(log_expr),
      baseline_expression = rowMeans(log_expr[, denom_samples, drop = FALSE], na.rm = TRUE)
    )
  })
}

normalize_de_table <- function(path, fdr = 0.05, lfc_min = 0.58, expected_methods = NULL, exclude_methods = NULL,
                               counts_path = NULL, sample_sheet_path = NULL, contrasts_path = NULL,
                               baseline_contrast_ids = NULL) {
  x <- read_table_auto(path)
  nm <- names(x)
  gcol <- first_present(nm, c("gene_id", "gene", "Gene", "name", "gene_name", "target_id"))
  ccol <- first_present(nm, c("contrast_id", "contrast", "comparison"))
  mcol <- first_present(nm, c("method", "Method"))
  ecol <- first_present(nm, c("signed_effect", "method_log2FC", "common_effect", "harmonized_signed_effect", "log2FC", "logFC", "effect", "estimate"))
  scol <- first_present(nm, c("signed_score", "harmonized_signed_score", "score", "signed_statistic"))
  bcol <- first_present(nm, c("baseline_expression", "baseline_expr", "Amean", "AveExpr", "baseMean", "mean_expression"))
  pcol <- first_present(nm, c("pval", "pvalue", "P.Value", "PValue", "p_value"))
  qcol <- first_present(nm, c("padj", "p_adj", "FDR", "fdr", "adj.P.Val", "qvalue", "q_value", "pfp", "q_for_score"))
  dcol <- first_present(nm, c("direction", "Direction", "call", "Call"))
  notecol <- first_present(nm, c("method_note", ".method_note", "note"))
  if (any(is.na(c(gcol, ccol, mcol, ecol)))) {
    stop("Unified DE table must contain gene_id, contrast_id, method, and signed_effect/log2FC columns. File: ", path)
  }
  out <- x %>%
    transmute(
      gene_id = as.character(.data[[gcol]]),
      contrast_id = as.character(.data[[ccol]]),
      method_raw = as.character(.data[[mcol]]),
      method = canonical_method(.data[[mcol]]),
      signed_effect = suppressWarnings(as.numeric(.data[[ecol]])),
      signed_score_input = if (!is.na(scol)) suppressWarnings(as.numeric(.data[[scol]])) else NA_real_,
      baseline_expression = if (!is.na(bcol)) suppressWarnings(as.numeric(.data[[bcol]])) else NA_real_,
      pval = if (!is.na(pcol)) suppressWarnings(as.numeric(.data[[pcol]])) else NA_real_,
      padj = if (!is.na(qcol)) suppressWarnings(as.numeric(.data[[qcol]])) else NA_real_,
      direction_input = if (!is.na(dcol)) as.character(.data[[dcol]]) else NA_character_,
      method_note = if (!is.na(notecol)) as.character(.data[[notecol]]) else NA_character_
    )
  if (!is.null(exclude_methods) && length(exclude_methods)) {
    out <- out %>% filter(!.data$method %in% exclude_methods)
  }
  if (!is.null(expected_methods) && length(expected_methods)) {
    out <- out %>% filter(.data$method %in% expected_methods)
    out$method <- factor(out$method, levels = expected_methods)
  }
  need_base <- any(!is.finite(out$baseline_expression))
  if (need_base && !is.null(counts_path) && !is.null(baseline_contrast_ids)) {
    b <- baseline_expression_for_contrasts(counts_path, sample_sheet_path, contrasts_path, baseline_contrast_ids)
    out <- out %>% select(-baseline_expression) %>% left_join(b, by = c("contrast_id", "gene_id"))
  }
  out <- out %>%
    mutate(
      direction_input_norm = case_when(
        stringr::str_detect(tolower(.data$direction_input), "^up|upregulated|increase") ~ "Up",
        stringr::str_detect(tolower(.data$direction_input), "^down|downregulated|decrease") ~ "Down",
        stringr::str_detect(tolower(.data$direction_input), "^ns|not|none|false") ~ "NS",
        TRUE ~ NA_character_
      ),
      q_for_score = dplyr::coalesce(.data$padj, .data$pval),
      has_q = is.finite(.data$q_for_score),
      # Prefer harmonized FDR/effect calling whenever a method provides a p/q value.
      # Direction-only calls are used only for precomputed methods that provide no p/q field.
      # This prevents native/precomputed PenDA rows with non-significant directional annotations
      # from being counted automatically as DE genes.
      sig_from_q = .data$has_q & is.finite(.data$signed_effect) &
        (.data$q_for_score < fdr) & abs(.data$signed_effect) >= lfc_min,
      sig_from_direction_only = !.data$has_q & .data$direction_input_norm %in% c("Up", "Down"),
      sig = .data$sig_from_q | .data$sig_from_direction_only,
      direction = case_when(
        .data$sig_from_q & .data$signed_effect > 0 ~ "Up",
        .data$sig_from_q & .data$signed_effect < 0 ~ "Down",
        .data$sig_from_direction_only ~ .data$direction_input_norm,
        TRUE ~ "NS"
      ),
      neglog10 = -log10(pmax(.data$q_for_score, .Machine$double.xmin)),
      neglog10 = pmin(.data$neglog10, 50),
      signed_score = dplyr::coalesce(.data$signed_score_input, sign(.data$signed_effect) * .data$neglog10)
    )
  out <- out %>%
    group_by(.data$contrast_id, .data$method) %>%
    mutate(score_for_agreement = rank_normal_score(dplyr::coalesce(.data$signed_score, .data$signed_effect))) %>%
    ungroup()
  out
}

jaccard_pairs <- function(df, methods, kind = c("All", "Up", "Down")) {
  kind <- match.arg(kind)
  d <- df %>% filter(.data$sig)
  if (kind == "Up") d <- d %>% filter(.data$direction == "Up")
  if (kind == "Down") d <- d %>% filter(.data$direction == "Down")
  sets <- split(d$gene_id, d$method)
  sets <- lapply(sets, unique)
  grid <- tidyr::expand_grid(method1 = methods, method2 = methods)
  grid %>% rowwise() %>% mutate(
    inter = length(intersect(sets[[as.character(method1)]] %||% character(), sets[[as.character(method2)]] %||% character())),
    union = length(union(sets[[as.character(method1)]] %||% character(), sets[[as.character(method2)]] %||% character())),
    value = ifelse(.data$union > 0, .data$inter / .data$union, NA_real_),
    kind = kind
  ) %>% ungroup()
}

`%||%` <- function(a, b) if (is.null(a)) b else a

rank_normal_score <- function(x) {
  # Within-method rank-normalized signed score for cross-method display,
  # Spearman concordance, and inter-method spread. Ties get average ranks.
  out <- rep(NA_real_, length(x))
  idx <- is.finite(x)
  n <- sum(idx)
  if (n < 3) return(out)
  r <- rank(x[idx], ties.method = "average", na.last = "keep")
  p <- (r - 0.5) / n
  out[idx] <- stats::qnorm(p)
  out[idx] <- pmax(pmin(out[idx], 5), -5)
  out
}


save_triplet <- function(plot, stem, out_dir, width, height, dpi = 600) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", device = grDevices::cairo_pdf, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(file.path(out_dir, paste0(stem, ".svg")), plot, width = width, height = height, units = "in", device = svglite::svglite, bg = "white")
  } else {
    warning("svglite is not installed; SVG output skipped.")
  }
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
}

main <- function() {
  args <- parse_args()
  root <- repo_root()
  de_results <- coalesce_arg(args[["de-results"]], file.path(root, "figure3/differential_expression/source_data/unified_de_results.csv"))
  out_dir <- coalesce_arg(args[["output-dir"]], file.path(root, "figure3/differential_expression/outputs"))
  src_dir <- coalesce_arg(args[["source-dir"]], file.path(root, "figure3/differential_expression/source_data"))
  counts_path <- coalesce_arg(args[["counts"]], file.path(root, "data/processed/raw_counts_rsemgenes.tsv"))
  sample_sheet_path <- coalesce_arg(args[["sample-sheet"]], file.path(root, "config/sample_sheet.csv"))
  contrasts_path <- coalesce_arg(args[["contrasts"]], file.path(root, "config/contrasts.csv"))
  rep_contrast <- coalesce_arg(args[["representative-contrast"]], "DT_vs_D_24h")
  early_contrast <- coalesce_arg(args[["early-contrast"]], "DT_vs_D_4h")
  output_stem <- coalesce_arg(args[["output-stem"]], "figure3_main")
  fdr <- as.numeric(coalesce_arg(args[["fdr"]], "0.05"))
  lfc_min <- as.numeric(coalesce_arg(args[["lfc-min"]], "0.58"))
  max_points <- as.integer(coalesce_arg(args[["max-points-per-method"]], "8000"))
  expected_methods <- c("DESeq2", "edgeR", "limma/voom", "RankProd", "PenDA", "RankCompV3")
  exclude_methods <- unlist(strsplit(coalesce_arg(args[["exclude-methods"]], "CellComp"), ","))
  exclude_methods <- exclude_methods[nzchar(exclude_methods)]
  dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  baseline_ids <- unique(c(rep_contrast, early_contrast))
  de <- normalize_de_table(
    de_results, fdr = fdr, lfc_min = lfc_min,
    expected_methods = expected_methods, exclude_methods = exclude_methods,
    counts_path = counts_path, sample_sheet_path = sample_sheet_path,
    contrasts_path = contrasts_path, baseline_contrast_ids = baseline_ids
  )
  missing_methods <- setdiff(expected_methods, unique(as.character(de$method)))
  if (length(missing_methods)) {
    warning("Missing methods in DE table after filtering: ", paste(missing_methods, collapse = ", "))
  }
  readr::write_tsv(de, file.path(src_dir, "figure3_unified_de_results_used.tsv"))

  de_rep <- de %>% filter(.data$contrast_id == rep_contrast)
  if (!nrow(de_rep)) stop("No rows for representative contrast: ", rep_contrast)
  methods <- expected_methods[expected_methods %in% unique(as.character(de_rep$method))]
  de_rep$method <- factor(as.character(de_rep$method), levels = methods)

  counts_tbl <- de_rep %>%
    filter(.data$direction != "NS") %>%
    count(.data$method, .data$direction, name = "n") %>%
    complete(method = factor(methods, levels = methods), direction = c("Down", "Up"), fill = list(n = 0)) %>%
    group_by(.data$method) %>% mutate(total = sum(.data$n)) %>% ungroup()
  readr::write_tsv(counts_tbl, file.path(src_dir, "figure3_panelA_call_counts.tsv"))
  counts_totals_tbl <- counts_tbl %>%
    distinct(.data$method, .data$total)

  dir_pal <- c(Up = "#D55E00", Down = "#0072B2", NS = "grey75")
  p_a <- ggplot(counts_tbl, aes(x = reorder(method, total), y = n, fill = direction)) +
    geom_col(width = 0.70) +
    # v12: total call-count labels make very stringent methods such as RankProd interpretable.
    geom_text(
      data = counts_totals_tbl,
      aes(x = reorder(method, total), y = total, label = scales::comma(total)),
      inherit.aes = FALSE, hjust = -0.12, size = 3.2
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = dir_pal, breaks = c("Up", "Down")) +
    labs(title = "DE calls by method", x = NULL, y = "Number of DE genes", fill = "Direction") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12), axis.text = element_text(size = 10), legend.position = "right")

  # Randomly downsample panel-B points within each method without using dplyr::n()
  # inside slice_sample(n = ...). Newer dplyr requires n to be a constant,
  # so sample the full group order first and then keep at most max_points rows.
  set.seed(1)
  plot_b <- de_rep %>%
    group_by(.data$method) %>%
    slice_sample(prop = 1, replace = FALSE) %>%
    slice_head(n = max_points) %>%
    ungroup()
  readr::write_tsv(plot_b, file.path(src_dir, "figure3_panelB_plot_points.tsv"))
  p_b <- ggplot(plot_b, aes(x = baseline_expression, y = score_for_agreement, colour = direction)) +
    geom_point(alpha = 0.45, size = 0.45, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
    facet_wrap(~ method, ncol = 3) +
    scale_colour_manual(values = dir_pal, drop = FALSE) +
    labs(title = paste0("Representative late contrast: ", rep_contrast),
         x = "Matched-control baseline expression", y = "Rank-normalized signed method score", colour = "Call") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12), strip.text = element_text(face = "bold", size = 9),
          axis.text = element_text(size = 8), legend.position = "right")

  jac_tbl <- bind_rows(
    jaccard_pairs(de_rep, methods, "All"),
    jaccard_pairs(de_rep, methods, "Up"),
    jaccard_pairs(de_rep, methods, "Down")
  )
  readr::write_tsv(jac_tbl, file.path(src_dir, "figure3_panelC_jaccard.tsv"))
  p_c <- jac_tbl %>%
    mutate(kind = factor(kind, levels = c("All", "Up", "Down")),
           method1 = factor(method1, levels = methods), method2 = factor(method2, levels = rev(methods)),
           label_colour = ifelse(is.finite(value) & value > 0.55, "white", "black"),
           jaccard_label = dplyr::case_when(
             !is.finite(value) ~ "",
             value > 0 & value < 0.01 ~ "<0.01",
             TRUE ~ sprintf("%.2f", value)
           )) %>%
    ggplot(aes(x = method1, y = method2, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = jaccard_label, colour = label_colour), size = 3.0) +
    scale_colour_identity() +
    scale_fill_gradient(limits = c(0, 1), low = "#f7fbff", high = "#08306b", na.value = "grey90") +
    coord_equal() + facet_wrap(~ kind, nrow = 1) +
    labs(title = "Pairwise Jaccard overlap", x = NULL, y = NULL, fill = "Jaccard") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12), axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8), strip.text = element_text(face = "bold", size = 9))

  wide_fc <- de_rep %>% select(gene_id, method, score_for_agreement) %>%
    tidyr::pivot_wider(names_from = method, values_from = score_for_agreement)
  cor_tbl <- expand_grid(method1 = methods, method2 = methods) %>%
    rowwise() %>%
    mutate(value = suppressWarnings(cor(wide_fc[[as.character(method1)]], wide_fc[[as.character(method2)]], method = "spearman", use = "pairwise.complete.obs"))) %>%
    ungroup()
  readr::write_tsv(cor_tbl, file.path(src_dir, "figure3_panelD_spearman.tsv"))
  p_d <- cor_tbl %>%
    mutate(method1 = factor(method1, levels = methods), method2 = factor(method2, levels = rev(methods)),
           label_colour = ifelse(is.finite(value) & value > 0.55, "white", "black")) %>%
    ggplot(aes(x = method1, y = method2, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = ifelse(is.finite(value), sprintf("%.2f", value), ""), colour = label_colour), size = 3.2) +
    scale_colour_identity() +
    scale_fill_gradient2(limits = c(-1, 1), low = "#67001f", mid = "white", high = "#053061", midpoint = 0, na.value = "grey90") +
    coord_equal() +
    labs(title = "Pairwise Spearman concordance", x = NULL, y = NULL, fill = expression(rho)) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12), axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8))

  de_e <- de %>% filter(.data$contrast_id %in% c(early_contrast, rep_contrast), is.finite(.data$score_for_agreement), is.finite(.data$baseline_expression))
  head_bins <- de_e %>%
    distinct(.data$contrast_id, .data$gene_id, .data$baseline_expression) %>%
    group_by(.data$contrast_id) %>%
    mutate(head_cut = quantile(.data$baseline_expression, probs = 0.90, na.rm = TRUE),
           baseline_bin = ifelse(.data$baseline_expression >= .data$head_cut, "Head", "Non-head")) %>%
    ungroup() %>% select(contrast_id, gene_id, baseline_bin)
  robust_z <- function(x) {
    m <- median(x, na.rm = TRUE)
    s <- stats::mad(x, center = m, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(s) || s == 0) s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
    (x - m) / s
  }
  spread_tbl <- de_e %>%
    left_join(head_bins, by = c("contrast_id", "gene_id")) %>%
    group_by(.data$contrast_id, .data$gene_id, .data$baseline_bin) %>%
    summarise(n_methods = sum(is.finite(.data$score_for_agreement)),
              spread_iqr = IQR(.data$score_for_agreement, na.rm = TRUE), .groups = "drop") %>%
    filter(.data$n_methods >= 3) %>%
    mutate(phase = ifelse(.data$contrast_id == early_contrast, "Early", "Late"),
           phase = factor(.data$phase, levels = c("Early", "Late")),
           baseline_bin = factor(.data$baseline_bin, levels = c("Non-head", "Head")))
  readr::write_tsv(spread_tbl, file.path(src_dir, "figure3_panelE_effect_spread.tsv"))
  if (!nrow(spread_tbl)) stop("Panel E has no effect-spread data; check early/late contrasts and method coverage.")
  p_e <- ggplot(spread_tbl, aes(x = phase, y = spread_iqr, fill = baseline_bin)) +
    geom_boxplot(outlier.shape = NA, width = 0.65, position = position_dodge(width = 0.75), linewidth = 0.35) +
    coord_cartesian(ylim = quantile(spread_tbl$spread_iqr, probs = c(0.01, 0.98), na.rm = TRUE)) +
    scale_fill_manual(values = c("Non-head" = "grey80", "Head" = "#D55E00")) +
    labs(title = "Inter-method signed-score spread localizes to the expression head",
         x = NULL, y = "Inter-method score spread (IQR)", fill = "Baseline bin") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12), axis.text = element_text(size = 10), legend.position = "right")

  row1 <- ggpubr::ggarrange(p_a, p_b, ncol = 2, widths = c(0.85, 1.75), labels = c("A", "B"), font.label = list(size = 16, face = "bold"))
  row2 <- ggpubr::ggarrange(p_c, p_d, ncol = 2, widths = c(1.35, 0.85), labels = c("C", "D"), font.label = list(size = 16, face = "bold"))
  fig <- ggpubr::ggarrange(row1, row2, p_e, ncol = 1, heights = c(1.20, 1.05, 0.75), labels = c("", "", "E"), font.label = list(size = 16, face = "bold"))
  save_triplet(fig, output_stem, out_dir, width = 12.5, height = 11.2, dpi = 600)

  manifest <- c(
    paste0("created: ", as.character(Sys.time())),
    paste0("de_results: ", normalizePath(de_results, mustWork = FALSE)),
    paste0("representative_contrast: ", rep_contrast),
    paste0("early_contrast_for_panelE: ", early_contrast),
    paste0("methods: ", paste(methods, collapse = ", ")),
    paste0("excluded_methods: ", paste(exclude_methods, collapse = ", ")),
    paste0("fdr: ", fdr),
    paste0("lfc_min: ", lfc_min),
    "Panel D uses harmonized signed method scores; white text is used on dark positive cells.",
    "RankCompV3 is the REO-family method used in place of the old CellComp/REOA slot.",
    "Panel E uses per-gene IQR across rank-normalized signed method scores; head is top 10% of matched-control baseline expression.",
    "Calling rule: methods with p/q values use FDR plus effect-size threshold; direction-only calls are used only when a precomputed method lacks p/q values."
  )
  writeLines(manifest, file.path(src_dir, "figure3_manifest.txt"))
  message("Wrote Figure 3 outputs to: ", out_dir)
}

main()
