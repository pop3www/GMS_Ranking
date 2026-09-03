#!/usr/bin/env Rscript
# Supplementary Fig. S4: gene-length sensitivity of head-constrained distortion.
#
# Purpose:
#   Test whether Figure 1/3/4 expression-head distortion is reducible to gene length / TOP1-CPT effects.
#
# Outputs:
#   supplement/figureS4/outputs/Supplementary_Fig_S4_gene_length_sensitivity.{pdf,svg,png}
#   supplement/figureS4/source_data/Supplementary_Fig_S4_gene_length_sensitivity_source_data.csv
#   supplement/figureS4/source_data/Supplementary_Fig_S4_model_summary.csv
#
# Inputs:
#   - config/contrasts.csv
#   - config/sample_sheet.csv
#   - data/processed/raw_counts_rsemgenes.tsv, unless Figure 4 gene table already exists
#   - gene lengths via --gene-lengths or --gtf
#   - Figure 3 used-results table for panel D

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
  if (basename(p) == "scripts") p <- normalizePath(file.path(p, "../.."), mustWork = FALSE)
  p
}
read_table_auto <- function(path) {
  ext <- tolower(tools::file_ext(gsub("\\.gz$", "", path)))
  if (ext %in% c("tsv", "txt", "gtf", "gff", "gff3")) readr::read_tsv(path, show_col_types = FALSE, comment = if (ext %in% c("gtf", "gff", "gff3")) "#" else "")
  else readr::read_csv(path, show_col_types = FALSE)
}
first_present <- function(nm, candidates) {
  hit <- candidates[candidates %in% nm]
  if (length(hit)) hit[[1]] else NA_character_
}
require_or_stop <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required R package is not installed: ", pkg, call. = FALSE)
}
normalize_gene_key <- function(x) {
  x <- as.character(x)
  ens <- stringr::str_extract(x, "ENSG[0-9]+")
  out <- ifelse(!is.na(ens), ens, x)
  out <- stringr::str_replace(out, "\\..*$", "")
  out
}
infer_gene_name <- function(x) {
  x <- as.character(x)
  nm <- ifelse(stringr::str_detect(x, "^ENSG[0-9]+_"), stringr::str_replace(x, "^ENSG[0-9]+_", ""), NA_character_)
  nm
}
safe_display_contrast <- function(cid) {
  cid %>%
    stringr::str_replace("_4h$", " (4 h)") %>%
    stringr::str_replace("_24h$", " (24 h)") %>%
    stringr::str_replace_all("_vs_", " vs ")
}
save_triplet <- function(plot, stem, out_dir, width, height, dpi = 600) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", device = grDevices::cairo_pdf, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(file.path(out_dir, paste0(stem, ".svg")), plot, width = width, height = height, units = "in", device = svglite::svglite, bg = "white")
  } else warning("svglite is not installed; SVG output skipped.")
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
}

condition_from_sample_sheet <- function(sample_sheet) {
  sample_sheet %>% mutate(
    time_h = as.integer(.data$time_h),
    condition = dplyr::case_when(
      .data$group_label %in% c("4_ctrl", "24_ctrl") ~ "Ctrl",
      .data$group_label %in% c("4_Tam", "24_Tam") ~ "Tam",
      TRUE ~ stringr::str_remove(.data$group_label, paste0("^", .data$time_h))
    )
  )
}
samples_for <- function(sample_sheet, condition, time_h) {
  sample_sheet %>% filter(.data$condition == condition, .data$time_h == as.integer(time_h)) %>% arrange(.data$replicate) %>% pull(.data$sample_id)
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

parse_gtf_attr <- function(attr, key) {
  pat <- paste0(key, " \\\"([^\\\"]+)\\\"")
  out <- stringr::str_match(attr, pat)[, 2]
  missing <- is.na(out)
  if (any(missing)) {
    pat2 <- paste0(key, "=([^;]+)")
    out[missing] <- stringr::str_match(attr[missing], pat2)[, 2]
  }
  out
}
build_lengths_from_gtf <- function(gtf_path) {
  message("Reading GTF/GFF for transcription-unit gene spans: ", gtf_path)
  coln <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  gtf <- readr::read_tsv(gtf_path, col_names = coln, comment = "#", show_col_types = FALSE, progress = FALSE)
  if (!all(c("feature", "start", "end", "attribute") %in% names(gtf))) stop("Could not parse GTF/GFF: ", gtf_path)
  genes <- gtf %>%
    mutate(
      gene_id_raw = parse_gtf_attr(.data$attribute, "gene_id"),
      gene_name_gtf = parse_gtf_attr(.data$attribute, "gene_name")
    ) %>%
    filter(!is.na(.data$gene_id_raw))
  if (any(genes$feature == "gene", na.rm = TRUE)) {
    genes <- genes %>% filter(.data$feature == "gene")
  }
  genes %>%
    mutate(gene_key = normalize_gene_key(.data$gene_id_raw)) %>%
    group_by(.data$gene_key) %>%
    summarise(
      gene_name = dplyr::first(na.omit(.data$gene_name_gtf)),
      gene_start = min(as.numeric(.data$start), na.rm = TRUE),
      gene_end = max(as.numeric(.data$end), na.rm = TRUE),
      gene_length_bp = .data$gene_end - .data$gene_start + 1,
      length_source = ifelse(any(.data$feature == "gene"), "GTF gene span", "GTF feature span"),
      .groups = "drop"
    ) %>%
    mutate(gene_name = ifelse(is.na(.data$gene_name) | !nzchar(.data$gene_name), NA_character_, .data$gene_name))
}
read_gene_lengths <- function(path = NULL, gtf = NULL) {
  if (!is.null(path) && nzchar(path) && file.exists(path)) {
    x <- read_table_auto(path)
    id_col <- first_present(names(x), c("gene_key", "gene_id", "gene", "ensembl_gene_id", "Geneid", "Name", "name"))
    if (is.na(id_col)) stop("Gene-length file needs gene_id/gene_key/ensembl_gene_id column: ", path)
    nm_col <- first_present(names(x), c("gene_name", "symbol", "external_gene_name", "GeneName"))
    len_col <- first_present(names(x), c("gene_length_bp", "length_bp", "gene_length", "length", "Length", "gene_span_bp"))
    start_col <- first_present(names(x), c("gene_start", "start", "Start"))
    end_col <- first_present(names(x), c("gene_end", "end", "End"))
    if (!is.na(len_col)) {
      out <- x %>% transmute(
        gene_key = normalize_gene_key(.data[[id_col]]),
        gene_name = if (!is.na(nm_col)) as.character(.data[[nm_col]]) else NA_character_,
        gene_length_bp = as.numeric(.data[[len_col]]),
        length_source = paste0("table:", basename(path))
      )
    } else if (!is.na(start_col) && !is.na(end_col)) {
      out <- x %>% transmute(
        gene_key = normalize_gene_key(.data[[id_col]]),
        gene_name = if (!is.na(nm_col)) as.character(.data[[nm_col]]) else NA_character_,
        gene_length_bp = as.numeric(.data[[end_col]]) - as.numeric(.data[[start_col]]) + 1,
        length_source = paste0("table-span:", basename(path))
      )
    } else stop("Gene-length file needs a length column or start/end columns: ", path)
    out %>%
      filter(!is.na(.data$gene_key), is.finite(.data$gene_length_bp), .data$gene_length_bp > 0) %>%
      group_by(.data$gene_key) %>%
      summarise(
        gene_name = dplyr::first(na.omit(.data$gene_name)),
        gene_length_bp = max(.data$gene_length_bp, na.rm = TRUE),
        length_source = dplyr::first(.data$length_source),
        .groups = "drop"
      )
  } else if (!is.null(gtf) && nzchar(gtf) && file.exists(gtf)) {
    build_lengths_from_gtf(gtf)
  } else stop("No usable gene length input. Provide --gene-lengths table or --gtf GENCODE.gtf[.gz].")
}

# Limma/voom expression table generation. Used only if Figure 4 gene table is unavailable.
run_one_contrast_gene_table <- function(counts_all, sample_sheet, contrast_row, alpha = 0.05, lfc_min = 0.58) {
  require_or_stop("edgeR"); require_or_stop("limma")
  cid <- contrast_row$contrast_id[[1]]
  selA <- samples_for(sample_sheet, contrast_row$denominator[[1]], contrast_row$time_h[[1]])
  selB <- samples_for(sample_sheet, contrast_row$numerator[[1]], contrast_row$time_h[[1]])
  selA <- intersect(selA, colnames(counts_all)); selB <- intersect(selB, colnames(counts_all))
  if (length(selA) < 2 || length(selB) < 2) stop("Need >=2 replicates per group for ", cid)
  X <- round(counts_all[, c(selA, selB), drop = FALSE]); X[X < 0] <- 0
  group <- factor(c(rep("A", length(selA)), rep("B", length(selB))), levels = c("A", "B"))
  coldata <- data.frame(sample = colnames(X), group = group); rownames(coldata) <- coldata$sample
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
  tibble(
    contrast_id = cid,
    gene_id = rownames(tt),
    baseline_expression = baseline[rownames(tt)],
    signed_effect = effect[rownames(tt)],
    pval = tt$P.Value,
    padj = tt$adj.P.Val
  )
}
load_or_build_gene_results <- function(fig4_gene_path, counts_path, sample_sheet, contrasts, alpha, lfc_min, min_total) {
  if (!is.null(fig4_gene_path) && file.exists(fig4_gene_path)) {
    message("Using existing Figure 4 gene results: ", fig4_gene_path)
    x <- readr::read_tsv(fig4_gene_path, show_col_types = FALSE)
    if (!all(c("contrast_id", "gene_id", "baseline_expression", "signed_effect") %in% names(x))) {
      stop("Figure 4 gene result table lacks required columns: ", fig4_gene_path)
    }
    return(x)
  }
  message("Figure 4 gene results not found; rebuilding limma/voom gene table from counts.")
  counts <- read_counts_matrix(counts_path)
  keep <- rowSums(counts, na.rm = TRUE) >= min_total
  counts <- counts[keep, , drop = FALSE]
  purrr::map_dfr(seq_len(nrow(contrasts)), ~ run_one_contrast_gene_table(counts, sample_sheet, contrasts[.x, ], alpha = alpha, lfc_min = lfc_min))
}

robust_coef <- function(df, formula, term) {
  require_or_stop("MASS")
  tryCatch({
    fit <- MASS::rlm(formula, data = df, maxit = 100)
    cf <- stats::coef(fit)
    if (term %in% names(cf)) as.numeric(cf[[term]]) else NA_real_
  }, error = function(e) NA_real_)
}
lm_p_value <- function(df, formula, term) {
  tryCatch({
    fit <- lm(formula, data = df)
    co <- summary(fit)$coefficients
    if (term %in% rownames(co)) as.numeric(co[term, "Pr(>|t|)"]) else NA_real_
  }, error = function(e) NA_real_)
}
boot_ci <- function(df, fit_fun, B = 1000, seed = 1) {
  est <- fit_fun(df)
  if (!is.finite(est) || is.null(B) || B <= 0 || nrow(df) < 10) return(c(estimate = est, ci_low = NA_real_, ci_high = NA_real_, boot_p = NA_real_))
  set.seed(seed)
  vals <- rep(NA_real_, B)
  n <- nrow(df)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    vals[[b]] <- tryCatch(fit_fun(df[idx, , drop = FALSE]), error = function(e) NA_real_)
  }
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(c(estimate = est, ci_low = NA_real_, ci_high = NA_real_, boot_p = NA_real_))
  p <- 2 * min(mean(vals <= 0), mean(vals >= 0))
  c(estimate = est, ci_low = unname(stats::quantile(vals, 0.025, na.rm = TRUE)), ci_high = unname(stats::quantile(vals, 0.975, na.rm = TRUE)), boot_p = min(1, p))
}
assign_rank_bands <- function(df) {
  df %>%
    group_by(.data$contrast_id) %>%
    mutate(
      baseline_rank = rank(-.data$baseline_expression, ties.method = "average", na.last = "keep"),
      baseline_rank_pct = percent_rank(.data$baseline_expression),
      baseline_band = case_when(
        .data$baseline_rank_pct >= 0.90 ~ "head",
        .data$baseline_rank_pct <= 0.10 ~ "tail",
        .data$baseline_rank_pct >= 0.25 & .data$baseline_rank_pct <= 0.75 ~ "mid",
        TRUE ~ "other"
      ),
      baseline_bin_head = ifelse(.data$baseline_rank_pct >= 0.90, "Head", "Non-head")
    ) %>%
    ungroup()
}
length_tertiles <- function(log_len) {
  qs <- stats::quantile(log_len, c(1/3, 2/3), na.rm = TRUE, names = FALSE)
  cut(log_len, breaks = c(-Inf, qs[[1]], qs[[2]], Inf), labels = c("Short", "Medium", "Long"), include.lowest = TRUE)
}

main <- function() {
  args <- parse_args()
  root <- repo_root()
  contrasts_path <- coalesce_arg(args[["contrasts"]], file.path(root, "config/contrasts.csv"))
  sample_sheet_path <- coalesce_arg(args[["sample-sheet"]], file.path(root, "config/sample_sheet.csv"))
  counts_path <- coalesce_arg(args[["counts"]], file.path(root, "data/processed/raw_counts_rsemgenes.tsv"))
  fig4_gene_path <- coalesce_arg(args[["figure4-gene-results"]], file.path(root, "figure3/scalar_metrics/source_data/figure4_limma_voom_gene_results.tsv"))
  fig3_used_path <- coalesce_arg(args[["figure3-used"]], file.path(root, "figure3/differential_expression/source_data/figure3_unified_de_results_used.tsv"))
  out_dir <- coalesce_arg(args[["output-dir"]], file.path(root, "supplement/figureS4/outputs"))
  src_dir <- coalesce_arg(args[["source-dir"]], file.path(root, "supplement/figureS4/source_data"))
  output_stem <- coalesce_arg(args[["output-stem"]], "Supplementary_Fig_S4_gene_length_sensitivity")
  alpha <- as.numeric(coalesce_arg(args[["fdr"]], "0.05"))
  lfc_min <- as.numeric(coalesce_arg(args[["lfc-min"]], "0.58"))
  min_total <- as.numeric(coalesce_arg(args[["min-total-count"]], "80"))
  B <- as.integer(coalesce_arg(args[["bootstrap"]], "1000"))
  if (isTRUE(args[["quick"]])) B <- 100L
  seed <- as.integer(coalesce_arg(args[["seed"]], "20260616"))
  gene_lengths_path <- coalesce_arg(args[["gene-lengths"]], "")
  gtf_path <- coalesce_arg(args[["gtf"]], "")

  dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Auto-detect gene length files if not passed explicitly.
  if (!nzchar(gene_lengths_path)) {
    candidates <- c(
      file.path(root, "data/annotations/gene_lengths.tsv"),
      file.path(root, "data/annotations/gencode_gene_lengths.tsv"),
      file.path(root, "data/processed/gene_lengths.tsv"),
      file.path(root, "supplement/figureS4/source_data/gene_lengths.tsv")
    )
    hit <- candidates[file.exists(candidates)]
    if (length(hit)) gene_lengths_path <- hit[[1]]
  }
  if (!nzchar(gtf_path)) {
    gtf_candidates <- c(
      Sys.glob(file.path(root, "data/annotations/*.gtf")),
      Sys.glob(file.path(root, "data/annotations/*.gtf.gz")),
      Sys.glob(file.path(root, "data/annotations/*.gff3")),
      Sys.glob(file.path(root, "data/annotations/*.gff3.gz"))
    )
    if (!nzchar(gene_lengths_path) && length(gtf_candidates)) gtf_path <- gtf_candidates[[1]]
  }

  if (!nzchar(gene_lengths_path) && !nzchar(gtf_path)) {
    # Write a template and stop clearly.
    message("No gene-length table/GTF found. Writing template and stopping.")
    if (file.exists(fig4_gene_path)) {
      genes <- readr::read_tsv(fig4_gene_path, show_col_types = FALSE) %>% distinct(.data$gene_id)
    } else {
      counts_head <- read_table_auto(counts_path) %>% select(1)
      names(counts_head)[1] <- "gene_id"
      genes <- counts_head
    }
    template <- genes %>% mutate(
      gene_key = normalize_gene_key(.data$gene_id),
      gene_name = infer_gene_name(.data$gene_id),
      gene_length_bp = NA_real_,
      note = "Fill gene_length_bp using transcription-unit span from GENCODE if available."
    )
    tmpl <- file.path(src_dir, "Supplementary_Fig_S4_gene_length_template.tsv")
    readr::write_tsv(template, tmpl)
    stop("No gene-length input. Provide --gene-lengths <table> or --gtf <GENCODE.gtf.gz>. Template written: ", tmpl, call. = FALSE)
  }

  contrasts <- readr::read_csv(contrasts_path, show_col_types = FALSE) %>% mutate(time_h = as.integer(.data$time_h))
  contrast_order <- contrasts$contrast_id
  sample_sheet <- condition_from_sample_sheet(readr::read_csv(sample_sheet_path, show_col_types = FALSE))
  gene_lengths <- read_gene_lengths(if (nzchar(gene_lengths_path)) gene_lengths_path else NULL, if (nzchar(gtf_path)) gtf_path else NULL) %>%
    mutate(log10_gene_length = log10(.data$gene_length_bp))
  readr::write_tsv(gene_lengths, file.path(src_dir, "Supplementary_Fig_S4_gene_lengths_used.tsv"))

  gene_results <- load_or_build_gene_results(fig4_gene_path, counts_path, sample_sheet, contrasts, alpha, lfc_min, min_total) %>%
    mutate(gene_key = normalize_gene_key(.data$gene_id), gene_name_from_id = infer_gene_name(.data$gene_id)) %>%
    left_join(contrasts, by = "contrast_id") %>%
    left_join(gene_lengths, by = "gene_key") %>%
    mutate(
      gene_name = coalesce(.data$gene_name, .data$gene_name_from_id),
      log2FC = .data$signed_effect
    ) %>%
    filter(is.finite(.data$baseline_expression), is.finite(.data$log2FC), is.finite(.data$log10_gene_length))

  if (nrow(gene_results) == 0) stop("No genes remained after merging expression results with gene lengths.")
  matched_genes <- gene_results %>% distinct(.data$gene_key) %>% nrow()
  message("Merged gene-length data for ", matched_genes, " genes.")

  gene_data <- gene_results %>%
    assign_rank_bands() %>%
    mutate(length_tertile = length_tertiles(.data$log10_gene_length))

  # Panel A representative contrasts.
  panelA_contrasts <- c("DT_vs_D_4h", "DT_vs_D_24h", "DT_H_CPT_vs_D_H_CPT_4h", "DT_H_CPT_vs_D_H_CPT_24h")
  panelA_data <- gene_data %>%
    filter(.data$contrast_id %in% panelA_contrasts, .data$baseline_band %in% c("head", "mid", "tail")) %>%
    mutate(
      baseline_band = factor(.data$baseline_band, levels = c("head", "mid", "tail")),
      contrast_label = factor(safe_display_contrast(.data$contrast_id), levels = safe_display_contrast(panelA_contrasts))
    )
  panelA_spearman <- gene_data %>%
    filter(.data$contrast_id %in% panelA_contrasts) %>%
    group_by(.data$contrast_id, .data$time_h) %>%
    summarise(
      n_genes = n(),
      spearman_expr_length = suppressWarnings(cor(.data$baseline_expression, .data$log10_gene_length, method = "spearman", use = "pairwise.complete.obs")),
      .groups = "drop"
    )

  # Panel B slopes across all contrasts.
  fit_unadj <- function(df) robust_coef(df, log2FC ~ baseline_expression, "baseline_expression")
  fit_adj <- function(df) robust_coef(df, log2FC ~ baseline_expression + log10_gene_length, "baseline_expression")
  fit_gamma <- function(df) robust_coef(df, log2FC ~ baseline_expression + log10_gene_length, "log10_gene_length")

  panelB <- gene_data %>%
    group_by(.data$contrast_id, .data$time_h, .data$numerator, .data$denominator, .data$cpt_level) %>%
    group_modify(function(.x, .y) {
      ci_u <- boot_ci(.x, fit_unadj, B = B, seed = seed + as.integer(factor(.y$contrast_id)))
      ci_a <- boot_ci(.x, fit_adj, B = B, seed = seed + 1000L + as.integer(factor(.y$contrast_id)))
      ci_g <- boot_ci(.x, fit_gamma, B = B, seed = seed + 2000L + as.integer(factor(.y$contrast_id)))
      tibble(
        n_genes = nrow(.x),
        unadjusted_beta_expr = ci_u[["estimate"]],
        unadjusted_beta_CI_low = ci_u[["ci_low"]],
        unadjusted_beta_CI_high = ci_u[["ci_high"]],
        unadjusted_bootstrap_p = ci_u[["boot_p"]],
        length_adjusted_beta_expr = ci_a[["estimate"]],
        length_adjusted_beta_CI_low = ci_a[["ci_low"]],
        length_adjusted_beta_CI_high = ci_a[["ci_high"]],
        length_adjusted_bootstrap_p = ci_a[["boot_p"]],
        gamma_length = ci_g[["estimate"]],
        gamma_length_CI_low = ci_g[["ci_low"]],
        gamma_length_CI_high = ci_g[["ci_high"]],
        gamma_length_bootstrap_p = ci_g[["boot_p"]]
      )
    }) %>% ungroup() %>%
    mutate(
      contrast_id = factor(.data$contrast_id, levels = contrast_order),
      contrast_label = factor(safe_display_contrast(as.character(.data$contrast_id)), levels = safe_display_contrast(contrast_order)),
      time_label = factor(paste0(.data$time_h, " h"), levels = c("4 h", "24 h"))
    )

  # Panel C length-stratified slopes.
  panelC_contrasts <- c("DT_vs_D_4h", "DT_vs_D_24h", "DT_L_CPT_vs_D_L_CPT_24h", "DT_H_CPT_vs_D_H_CPT_24h", "D_H_CPT_vs_D_24h")
  panelC <- gene_data %>%
    filter(.data$contrast_id %in% panelC_contrasts, !is.na(.data$length_tertile)) %>%
    group_by(.data$contrast_id, .data$time_h, .data$length_tertile) %>%
    group_modify(function(.x, .y) {
      ci <- boot_ci(.x, fit_unadj, B = B, seed = seed + 3000L + as.integer(factor(paste(.y$contrast_id, .y$length_tertile))))
      tibble(n_genes = nrow(.x), beta_expr = ci[["estimate"]], beta_expr_CI_low = ci[["ci_low"]], beta_expr_CI_high = ci[["ci_high"]], p_value_or_bootstrap_p = ci[["boot_p"]])
    }) %>% ungroup() %>%
    mutate(
      contrast_id = factor(.data$contrast_id, levels = panelC_contrasts),
      contrast_label = factor(safe_display_contrast(as.character(.data$contrast_id)), levels = safe_display_contrast(panelC_contrasts)),
      length_tertile = factor(as.character(.data$length_tertile), levels = c("Short", "Medium", "Long")),
      time_label = factor(paste0(.data$time_h, " h"), levels = c("4 h", "24 h"))
    )

  # Panel D: Figure 3 signed-score spread, length-adjusted head effect.
  if (!file.exists(fig3_used_path)) stop("Panel D requires Figure 3 used-results table: ", fig3_used_path)
  fig3_used <- readr::read_tsv(fig3_used_path, show_col_types = FALSE)
  if (!all(c("gene_id", "contrast_id", "method", "score_for_agreement", "baseline_expression") %in% names(fig3_used))) {
    stop("Figure 3 used-results table needs gene_id, contrast_id, method, score_for_agreement, baseline_expression: ", fig3_used_path)
  }
  panelD_gene <- fig3_used %>%
    filter(.data$contrast_id %in% c("DT_vs_D_4h", "DT_vs_D_24h")) %>%
    mutate(gene_key = normalize_gene_key(.data$gene_id)) %>%
    group_by(.data$contrast_id, .data$gene_id, .data$gene_key) %>%
    summarise(
      baseline_expression = dplyr::first(.data$baseline_expression),
      spread_i = IQR(.data$score_for_agreement, na.rm = TRUE),
      n_methods = n_distinct(.data$method),
      .groups = "drop"
    ) %>%
    left_join(gene_lengths, by = "gene_key") %>%
    filter(.data$n_methods >= 6, is.finite(.data$spread_i), is.finite(.data$log10_gene_length), is.finite(.data$baseline_expression)) %>%
    group_by(.data$contrast_id) %>%
    mutate(
      baseline_rank_pct = percent_rank(.data$baseline_expression),
      baseline_bin = ifelse(.data$baseline_rank_pct >= 0.90, "Head", "Non-head"),
      head_indicator = ifelse(.data$baseline_bin == "Head", 1, 0),
      time_bin = ifelse(.data$contrast_id == "DT_vs_D_4h", "Early", "Late")
    ) %>% ungroup()

  fit_head <- function(df) robust_coef(df, spread_i ~ head_indicator + log10_gene_length, "head_indicator")
  fit_lenD <- function(df) robust_coef(df, spread_i ~ head_indicator + log10_gene_length, "log10_gene_length")
  panelD <- panelD_gene %>%
    group_by(.data$contrast_id, .data$time_bin) %>%
    group_modify(function(.x, .y) {
      ci_h <- boot_ci(.x, fit_head, B = B, seed = seed + 4000L + as.integer(factor(.y$contrast_id)))
      ci_l <- boot_ci(.x, fit_lenD, B = B, seed = seed + 5000L + as.integer(factor(.y$contrast_id)))
      tibble(
        n_genes = nrow(.x),
        length_adjusted_head_effect = ci_h[["estimate"]],
        head_effect_CI_low = ci_h[["ci_low"]],
        head_effect_CI_high = ci_h[["ci_high"]],
        head_effect_bootstrap_p = ci_h[["boot_p"]],
        gamma_length = ci_l[["estimate"]],
        gamma_length_CI_low = ci_l[["ci_low"]],
        gamma_length_CI_high = ci_l[["ci_high"]],
        gamma_length_bootstrap_p = ci_l[["boot_p"]]
      )
    }) %>% ungroup() %>%
    mutate(time_bin = factor(.data$time_bin, levels = c("Early", "Late")))

  # Model summary table.
  summary_A <- panelA_spearman %>% transmute(
    contrast = .data$contrast_id, time_point = .data$time_h, model = "A_spearman_baseline_expression_vs_log10_length",
    n_genes = .data$n_genes, beta_expr = NA_real_, beta_expr_CI_low = NA_real_, beta_expr_CI_high = NA_real_,
    gamma_length = NA_real_, gamma_length_CI_low = NA_real_, gamma_length_CI_high = NA_real_,
    spearman_expr_length = .data$spearman_expr_length, p_value_or_bootstrap_p = NA_real_,
    length_adjusted_head_effect = NA_real_, head_effect_CI_low = NA_real_, head_effect_CI_high = NA_real_
  )
  summary_B1 <- panelB %>% transmute(
    contrast = as.character(.data$contrast_id), time_point = .data$time_h, model = "B_unadjusted_amplifier_slope",
    n_genes = .data$n_genes, beta_expr = .data$unadjusted_beta_expr, beta_expr_CI_low = .data$unadjusted_beta_CI_low, beta_expr_CI_high = .data$unadjusted_beta_CI_high,
    gamma_length = NA_real_, gamma_length_CI_low = NA_real_, gamma_length_CI_high = NA_real_, spearman_expr_length = NA_real_, p_value_or_bootstrap_p = .data$unadjusted_bootstrap_p,
    length_adjusted_head_effect = NA_real_, head_effect_CI_low = NA_real_, head_effect_CI_high = NA_real_
  )
  summary_B2 <- panelB %>% transmute(
    contrast = as.character(.data$contrast_id), time_point = .data$time_h, model = "B_length_adjusted_amplifier_slope",
    n_genes = .data$n_genes, beta_expr = .data$length_adjusted_beta_expr, beta_expr_CI_low = .data$length_adjusted_beta_CI_low, beta_expr_CI_high = .data$length_adjusted_beta_CI_high,
    gamma_length = .data$gamma_length, gamma_length_CI_low = .data$gamma_length_CI_low, gamma_length_CI_high = .data$gamma_length_CI_high, spearman_expr_length = NA_real_, p_value_or_bootstrap_p = .data$length_adjusted_bootstrap_p,
    length_adjusted_head_effect = NA_real_, head_effect_CI_low = NA_real_, head_effect_CI_high = NA_real_
  )
  summary_C <- panelC %>% transmute(
    contrast = as.character(.data$contrast_id), time_point = .data$time_h, model = paste0("C_length_stratum_", as.character(.data$length_tertile)),
    n_genes = .data$n_genes, beta_expr = .data$beta_expr, beta_expr_CI_low = .data$beta_expr_CI_low, beta_expr_CI_high = .data$beta_expr_CI_high,
    gamma_length = NA_real_, gamma_length_CI_low = NA_real_, gamma_length_CI_high = NA_real_, spearman_expr_length = NA_real_, p_value_or_bootstrap_p = .data$p_value_or_bootstrap_p,
    length_adjusted_head_effect = NA_real_, head_effect_CI_low = NA_real_, head_effect_CI_high = NA_real_
  )
  summary_D <- panelD %>% transmute(
    contrast = .data$contrast_id, time_point = ifelse(.data$time_bin == "Early", 4, 24), model = "D_length_adjusted_head_effect_on_inter_method_spread",
    n_genes = .data$n_genes, beta_expr = NA_real_, beta_expr_CI_low = NA_real_, beta_expr_CI_high = NA_real_,
    gamma_length = .data$gamma_length, gamma_length_CI_low = .data$gamma_length_CI_low, gamma_length_CI_high = .data$gamma_length_CI_high, spearman_expr_length = NA_real_, p_value_or_bootstrap_p = .data$head_effect_bootstrap_p,
    length_adjusted_head_effect = .data$length_adjusted_head_effect, head_effect_CI_low = .data$head_effect_CI_low, head_effect_CI_high = .data$head_effect_CI_high
  )
  model_summary <- bind_rows(summary_A, summary_B1, summary_B2, summary_C, summary_D)

  # Source data table. Slopes are repeated per gene/contrast for convenience and provenance.
  source_data <- gene_data %>%
    left_join(panelB %>% select(contrast_id, unadjusted_beta_expr, length_adjusted_beta_expr, gamma_length, model_gene_count = n_genes), by = "contrast_id") %>%
    left_join(panelD_gene %>% select(contrast_id, gene_key, spread_i), by = c("contrast_id", "gene_key")) %>%
    left_join(panelD %>% select(contrast_id, length_adjusted_head_effect), by = "contrast_id") %>%
    transmute(
      gene_id = .data$gene_id,
      gene_name = .data$gene_name,
      gene_key = .data$gene_key,
      gene_length_bp = .data$gene_length_bp,
      log10_gene_length = .data$log10_gene_length,
      contrast = .data$contrast_id,
      time_point = .data$time_h,
      baseline_condition = .data$denominator,
      numerator_condition = .data$numerator,
      baseline_expr = .data$baseline_expression,
      baseline_rank = .data$baseline_rank,
      baseline_band = .data$baseline_band,
      length_tertile = as.character(.data$length_tertile),
      log2FC = .data$log2FC,
      unadjusted_beta_expr = .data$unadjusted_beta_expr,
      length_adjusted_beta_expr = .data$length_adjusted_beta_expr,
      gamma_length = .data$gamma_length,
      spread_i = .data$spread_i,
      length_adjusted_head_effect = .data$length_adjusted_head_effect,
      model_gene_count = .data$model_gene_count
    )

  readr::write_csv(source_data, file.path(src_dir, "Supplementary_Fig_S4_gene_length_sensitivity_source_data.csv"))
  readr::write_csv(model_summary, file.path(src_dir, "Supplementary_Fig_S4_model_summary.csv"))
  readr::write_tsv(panelA_spearman, file.path(src_dir, "Supplementary_Fig_S4_panelA_expr_length_spearman.tsv"))
  readr::write_tsv(panelB, file.path(src_dir, "Supplementary_Fig_S4_panelB_slopes.tsv"))
  readr::write_tsv(panelC, file.path(src_dir, "Supplementary_Fig_S4_panelC_length_stratified_slopes.tsv"))
  readr::write_tsv(panelD_gene, file.path(src_dir, "Supplementary_Fig_S4_panelD_gene_spread_with_length.tsv"))
  readr::write_tsv(panelD, file.path(src_dir, "Supplementary_Fig_S4_panelD_head_effect.tsv"))

  time_pal <- c("4 h" = "#0072B2", "24 h" = "#D55E00")
  length_pal <- c("Short" = "#56B4E9", "Medium" = "#999999", "Long" = "#CC79A7")
  band_pal <- c("tail" = "#BDBDBD", "mid" = "#737373", "head" = "#D55E00")

  pA <- ggplot(panelA_data, aes(x = baseline_band, y = log10_gene_length, fill = baseline_band)) +
    geom_violin(width = 0.85, scale = "width", color = NA, alpha = 0.55) +
    geom_boxplot(width = 0.16, outlier.size = 0.15, fill = "white", alpha = 0.75) +
    facet_wrap(~ contrast_label, ncol = 2) +
    scale_fill_manual(values = band_pal, drop = FALSE) +
    labs(title = "Gene length across expression-rank bands", x = "Baseline expression band", y = expression(log[10]*" gene length (bp)")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 11), strip.text = element_text(face = "bold", size = 8), axis.text.x = element_text(size = 8), axis.title = element_text(size = 9))

  pB_data <- panelB %>%
    select(contrast_id, contrast_label, time_label, unadjusted_beta_expr, length_adjusted_beta_expr, length_adjusted_beta_CI_low, length_adjusted_beta_CI_high) %>%
    mutate(contrast_label = factor(.data$contrast_label, levels = safe_display_contrast(contrast_order)))
  pB <- ggplot(pB_data, aes(x = contrast_label, y = length_adjusted_beta_expr, color = time_label)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
    geom_point(aes(y = unadjusted_beta_expr), color = "grey45", size = 1.6, alpha = 0.8, shape = 1) +
    geom_errorbar(aes(ymin = length_adjusted_beta_CI_low, ymax = length_adjusted_beta_CI_high), width = 0.18, linewidth = 0.35, na.rm = TRUE) +
    geom_point(size = 2.1) +
    coord_flip() +
    scale_color_manual(values = time_pal, drop = FALSE) +
    labs(title = "Length-adjusted amplifier slope", x = NULL, y = "Baseline-expression slope", color = "Time point") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11), axis.text.y = element_text(size = 7.2), axis.title = element_text(size = 9), legend.position = "right")

  pC <- ggplot(panelC, aes(x = contrast_label, y = beta_expr, color = length_tertile, group = length_tertile)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
    geom_errorbar(aes(ymin = beta_expr_CI_low, ymax = beta_expr_CI_high), position = position_dodge(width = 0.55), width = 0.18, linewidth = 0.35, na.rm = TRUE) +
    geom_point(position = position_dodge(width = 0.55), size = 2.0) +
    coord_flip() +
    scale_color_manual(values = length_pal, drop = FALSE) +
    labs(title = "Amplifier slope within length strata", x = NULL, y = "Baseline-expression slope", color = "Length tertile") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11), axis.text.y = element_text(size = 7.5), axis.title = element_text(size = 9), legend.position = "right")

  pD <- ggplot(panelD, aes(x = time_bin, y = length_adjusted_head_effect, color = time_bin)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
    geom_errorbar(aes(ymin = head_effect_CI_low, ymax = head_effect_CI_high), width = 0.12, linewidth = 0.4, na.rm = TRUE) +
    geom_point(size = 2.6) +
    scale_color_manual(values = c("Early" = time_pal[["4 h"]], "Late" = time_pal[["24 h"]]), guide = "none") +
    labs(title = "Length-adjusted head effect on method spread", x = NULL, y = "Head effect on spread") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11), axis.title = element_text(size = 9), axis.text = element_text(size = 8))

  top <- ggpubr::ggarrange(pA, pB, ncol = 2, widths = c(1.08, 1.0), labels = c("A", "B"), font.label = list(size = 16, face = "bold"))
  bottom <- ggpubr::ggarrange(pC, pD, ncol = 2, widths = c(1.25, 0.75), labels = c("C", "D"), font.label = list(size = 16, face = "bold"))
  fig <- ggpubr::ggarrange(top, bottom, ncol = 1, heights = c(1.05, 0.95))
  save_triplet(fig, output_stem, out_dir, width = 12.5, height = 10.2, dpi = 600)

  manifest <- c(
    paste0("created: ", Sys.time()),
    paste0("contrasts: ", normalizePath(contrasts_path, mustWork = FALSE)),
    paste0("sample_sheet: ", normalizePath(sample_sheet_path, mustWork = FALSE)),
    paste0("figure4_gene_results: ", normalizePath(fig4_gene_path, mustWork = FALSE)),
    paste0("figure3_used_results: ", normalizePath(fig3_used_path, mustWork = FALSE)),
    paste0("gene_lengths: ", if (nzchar(gene_lengths_path)) normalizePath(gene_lengths_path, mustWork = FALSE) else ""),
    paste0("gtf: ", if (nzchar(gtf_path)) normalizePath(gtf_path, mustWork = FALSE) else ""),
    paste0("n_gene_keys_with_lengths_used: ", matched_genes),
    paste0("bootstrap_gene_resamples: ", B),
    paste0("fdr: ", alpha),
    paste0("lfc_min: ", lfc_min),
    "Gene length is transcription-unit/gene span if derived from GTF gene features; if a length table was supplied, see length_source in Supplementary_Fig_S4_gene_lengths_used.tsv.",
    "Panel B: gray open circles show unadjusted robust slopes; colored points/error bars show length-adjusted robust slopes.",
    "Panel D: head effect estimated by robust regression spread_i ~ head_indicator + log10_gene_length using Figure 3 score_for_agreement."
  )
  writeLines(manifest, file.path(src_dir, "Supplementary_Fig_S4_manifest.txt"))
  message("Wrote Supplementary Fig. S4 outputs to: ", out_dir)
}

main()
