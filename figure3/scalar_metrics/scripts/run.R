#!/usr/bin/env Rscript
# Figure 4 production script: two compact scalar summaries across contrasts.
# A = Up fraction. B = global amplifier slope.
# Uses rounded RSEM expected counts with TMM/logCPM + limma/voom.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggpubr)
})

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(); i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      if (i == length(args) || startsWith(args[[i + 1]], "--")) { out[[key]] <- TRUE; i <- i + 1 }
      else { out[[key]] <- args[[i + 1]]; i <- i + 2 }
    } else i <- i + 1
  }
  out
}
coalesce_arg <- function(x, default) if (is.null(x) || isTRUE(x) || !nzchar(as.character(x))) default else as.character(x)
repo_root <- function() {
  p <- normalizePath(getwd(), mustWork = FALSE)
  if (basename(p) == "scripts") p <- normalizePath(file.path(p, "../../.."), mustWork = FALSE)
  p
}
read_table_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) readr::read_tsv(path, show_col_types = FALSE) else readr::read_csv(path, show_col_types = FALSE)
}
first_present <- function(nm, candidates) {
  hit <- candidates[candidates %in% nm]
  if (length(hit)) hit[[1]] else NA_character_
}
require_or_stop <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required R package is not installed: ", pkg, call. = FALSE)
}
condition_from_sample_sheet <- function(sample_sheet) {
  sample_sheet %>% mutate(
    time_h = as.integer(time_h),
    condition = dplyr::case_when(
      group_label %in% c("4_ctrl", "24_ctrl") ~ "Ctrl",
      group_label %in% c("4_Tam", "24_Tam") ~ "Tam",
      TRUE ~ stringr::str_remove(group_label, paste0("^", time_h))
    )
  )
}
samples_for <- function(sample_sheet, condition_value, time_value) {
  sample_sheet %>%
    dplyr::filter(
      .data$condition == .env$condition_value,
      .data$time_h == as.integer(.env$time_value)
    ) %>%
    dplyr::arrange(.data$replicate) %>%
    dplyr::pull(.data$sample_id)
}

read_counts_matrix <- function(path) {
  x <- read_table_auto(path)
  gcol <- first_present(names(x), c("gene_id", "gene", "name", "Gene", "Geneid"))
  if (is.na(gcol)) stop("Count matrix needs gene_id/gene/name column: ", path)
  x <- x %>% rename(gene_id = all_of(gcol))
  mat <- as.data.frame(x)
  rownames(mat) <- make.unique(as.character(mat$gene_id))
  mat$gene_id <- NULL
  mat[] <- lapply(mat, function(v) suppressWarnings(as.numeric(v)))
  as.matrix(mat)
}
safe_display_contrast <- function(cid) {
  cid %>%
    stringr::str_replace("_4h$", " (4 h)") %>%
    stringr::str_replace("_24h$", " (24 h)") %>%
    stringr::str_replace_all("_vs_", " vs ")
}
save_triplet <- function(plot, stem, out_dir, width, height, dpi = 600) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", device = grDevices::cairo_pdf)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(file.path(out_dir, paste0(stem, ".svg")), plot, width = width, height = height, units = "in", device = svglite::svglite)
  } else warning("svglite is not installed; SVG output skipped.")
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = dpi)
}

run_one_contrast <- function(counts_all, sample_sheet, contrast_row, alpha = 0.05, lfc_min = 0.58) {
  require_or_stop("edgeR"); require_or_stop("limma"); require_or_stop("MASS")
  cid <- contrast_row$contrast_id[[1]]
  
  selA <- samples_for(sample_sheet, contrast_row$denominator[[1]], contrast_row$time_h[[1]]) |> unique()
  selB <- samples_for(sample_sheet, contrast_row$numerator[[1]], contrast_row$time_h[[1]]) |> unique()
  
  selA <- intersect(selA, colnames(counts_all))
  selB <- intersect(selB, colnames(counts_all))
  
  overlap <- intersect(selA, selB)
  if (length(overlap) > 0) {
    stop(
      "Numerator and denominator groups overlap for contrast ", cid, ": ",
      paste(overlap, collapse = ", "),
      call. = FALSE
    )
  }
  
  dups <- c(selA, selB)[duplicated(c(selA, selB))]
  if (length(dups) > 0) {
    stop(
      "Duplicated selected samples for contrast ", cid, ": ",
      paste(unique(dups), collapse = ", "),
      call. = FALSE
    )
  }
  
  if (length(selA) < 2 || length(selB) < 2) {
    stop(
      "Need >=2 replicates per group for ", cid,
      ". denominator=", contrast_row$denominator[[1]],
      " n=", length(selA),
      "; numerator=", contrast_row$numerator[[1]],
      " n=", length(selB),
      call. = FALSE
    )
  }
  
  message(
    "Contrast ", cid,
    " | denominator ", contrast_row$denominator[[1]], ": ", paste(selA, collapse = ", "),
    " | numerator ", contrast_row$numerator[[1]], ": ", paste(selB, collapse = ", ")
  )
  
  X <- round(counts_all[, c(selA, selB), drop = FALSE])

    X[X < 0] <- 0
  group <- factor(c(rep("A", length(selA)), rep("B", length(selB))), levels = c("A", "B"))
  coldata <- data.frame(sample = colnames(X), group = group)
  rownames(coldata) <- coldata$sample

  y <- edgeR::DGEList(counts = X)
  keep <- rowSums(edgeR::cpm(y) > 1) >= 2
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  design <- model.matrix(~ group, data = coldata)
  v <- limma::voom(y, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")

  logCPM <- edgeR::cpm(y, prior.count = 1, log = TRUE)
  baseline <- rowMeans(logCPM[, group == "A", drop = FALSE], na.rm = TRUE)
  effect <- rowMeans(logCPM[, group == "B", drop = FALSE], na.rm = TRUE) - baseline
  de <- tibble(
    gene_id = rownames(tt),
    baseline_expression = baseline[rownames(tt)],
    signed_effect = effect[rownames(tt)],
    pval = tt$P.Value,
    padj = tt$adj.P.Val
  ) %>% mutate(
    sig = is.finite(.data$padj) & .data$padj < alpha & is.finite(.data$signed_effect) & abs(.data$signed_effect) >= lfc_min,
    direction = case_when(.data$sig & .data$signed_effect > 0 ~ "Up", .data$sig & .data$signed_effect < 0 ~ "Down", TRUE ~ "NS")
  )
  n_up <- sum(de$direction == "Up", na.rm = TRUE)
  n_down <- sum(de$direction == "Down", na.rm = TRUE)
  upfrac <- if ((n_up + n_down) > 0) n_up / (n_up + n_down) else NA_real_
  slope_df <- de %>% filter(is.finite(.data$baseline_expression), is.finite(.data$signed_effect))
  slope_rlm <- tryCatch(as.numeric(coef(MASS::rlm(signed_effect ~ baseline_expression, data = slope_df))[2]), error = function(e) NA_real_)
  slope_lm <- tryCatch(as.numeric(coef(lm(signed_effect ~ baseline_expression, data = slope_df))[2]), error = function(e) NA_real_)
  rho <- suppressWarnings(cor(slope_df$baseline_expression, slope_df$signed_effect, method = "spearman", use = "pairwise.complete.obs"))
  list(
    gene_table = de %>% mutate(contrast_id = cid, .before = 1),
    summary = tibble(
      contrast_id = cid,
      numerator = contrast_row$numerator[[1]], denominator = contrast_row$denominator[[1]],
      time_h = as.integer(contrast_row$time_h[[1]]), cpt_level = contrast_row$cpt_level[[1]],
      nA = length(selA), nB = length(selB), n_genes_tested = nrow(de),
      Up = n_up, Down = n_down, UpFrac = upfrac,
      slope_rlm = slope_rlm, slope_lm = slope_lm, spearman_baseline_effect = as.numeric(rho)
    )
  )
}

main <- function() {
  args <- parse_args()
  root <- repo_root()
  counts_path <- coalesce_arg(args[["counts"]], file.path(root, "data/processed/raw_counts_rsemgenes.tsv"))
  sample_sheet_path <- coalesce_arg(args[["sample-sheet"]], file.path(root, "config/sample_sheet.csv"))
  contrasts_path <- coalesce_arg(args[["contrasts"]], file.path(root, "config/contrasts.csv"))
  out_dir <- coalesce_arg(args[["output-dir"]], file.path(root, "figure3/scalar_metrics/outputs"))
  src_dir <- coalesce_arg(args[["source-dir"]], file.path(root, "figure3/scalar_metrics/source_data"))
  output_stem <- coalesce_arg(args[["output-stem"]], "figure4_main")
  alpha <- as.numeric(coalesce_arg(args[["fdr"]], "0.05"))
  lfc_min <- as.numeric(coalesce_arg(args[["lfc-min"]], "0.58"))
  min_total <- as.numeric(coalesce_arg(args[["min-total-count"]], "80"))
  dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  counts <- read_counts_matrix(counts_path)
  keep <- rowSums(counts, na.rm = TRUE) >= min_total
  counts <- counts[keep, , drop = FALSE]
  sample_sheet <- condition_from_sample_sheet(readr::read_csv(sample_sheet_path, show_col_types = FALSE)) %>%
    dplyr::distinct(.data$sample_id, .keep_all = TRUE)
  contrasts <- readr::read_csv(contrasts_path, show_col_types = FALSE)
  if (!"contrast_id" %in% names(contrasts)) stop("contrasts.csv must contain contrast_id")
  contrast_order <- contrasts$contrast_id

  runs <- purrr::map(seq_len(nrow(contrasts)), ~ run_one_contrast(counts, sample_sheet, contrasts[.x, ], alpha = alpha, lfc_min = lfc_min))
  gene_tbl <- bind_rows(purrr::map(runs, "gene_table"))
  summary_tbl <- bind_rows(purrr::map(runs, "summary")) %>%
    mutate(
      contrast_id = factor(.data$contrast_id, levels = contrast_order),
      contrast_label = safe_display_contrast(as.character(.data$contrast_id)),
      contrast_label = factor(.data$contrast_label, levels = safe_display_contrast(contrast_order)),
      time_label = factor(paste0(.data$time_h, " h"), levels = c("4 h", "24 h")),
      contrast_family = case_when(
        numerator == "DT" & denominator == "D" ~ "MYC on Dox",
        numerator == "Tam" & denominator == "Ctrl" ~ "Tam-only control",
        str_detect(numerator, "_CPT") & str_detect(numerator, "^DT") ~ "MYC on CPT+Dox",
        str_detect(numerator, "_CPT") & str_detect(numerator, "^D") ~ "CPT+Dox control",
        numerator == "D" & denominator == "Ctrl" ~ "Dox-only control",
        TRUE ~ "Other"
      )
    )
  readr::write_tsv(gene_tbl, file.path(src_dir, "figure4_limma_voom_gene_results.tsv"))
  readr::write_tsv(summary_tbl, file.path(src_dir, "figure4_scalar_metrics.tsv"))

  time_pal <- c("4 h" = "#0072B2", "24 h" = "#D55E00")
  p_a <- ggplot(summary_tbl, aes(x = contrast_label, y = UpFrac, fill = time_label)) +
    geom_col(width = 0.72) + coord_flip() +
    scale_fill_manual(values = time_pal, drop = FALSE) +
    labs(title = "Up fraction", x = NULL, y = "Up / (Up + Down)", fill = "Time point") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 14), axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 12), legend.position = "right")
  p_b <- ggplot(summary_tbl, aes(x = contrast_label, y = slope_rlm, fill = time_label)) +
    geom_col(width = 0.72) + coord_flip() +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
    scale_fill_manual(values = time_pal, drop = FALSE) +
    labs(title = "Global amplifier slope", x = NULL, y = "Robust slope of effect vs baseline", fill = "Time point") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 14), axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 12), legend.position = "right")
  fig <- ggpubr::ggarrange(p_a, p_b, ncol = 1, labels = c("A", "B"), heights = c(1, 1), font.label = list(size = 18, face = "bold"))
  save_triplet(fig, output_stem, out_dir, width = 8.2, height = 9.0, dpi = 600)
  manifest <- c(
    paste0("created: ", as.character(Sys.time())),
    paste0("counts: ", normalizePath(counts_path, mustWork = FALSE)),
    paste0("contrast_order: ", paste(contrast_order, collapse = ", ")),
    paste0("fdr: ", alpha),
    paste0("lfc_min: ", lfc_min),
    paste0("min_total_count_filter: ", min_total),
    "Input counts are rounded RSEM expected counts for limma/voom/edgeR count-scale analysis.",
    "Colors encode time point: 4 h = #0072B2, 24 h = #D55E00.",
    "Up fraction uses limma/voom DE calls. Slope is robust signed_effect ~ matched-control baseline_expression."
  )
  writeLines(manifest, file.path(src_dir, "figure4_manifest.txt"))
  message("Wrote Figure 4 outputs to: ", out_dir)
}

main()
