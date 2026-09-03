#!/usr/bin/env Rscript
# Optional local builder for the Figure 3 unified DE table.

# This script is intentionally explicit about dependencies and fallbacks.
# Preferred production mode: run all six methods locally, then plot with plot_figure3.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
})

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1
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
require_or_stop <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package is not installed: ", pkg, call. = FALSE)
  }
}
first_present <- function(nm, candidates) {
  hit <- candidates[candidates %in% nm]
  if (length(hit)) hit[[1]] else NA_character_
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
  condition_value <- as.character(condition_value)
  time_value <- as.integer(time_value)
  sample_sheet %>%
    dplyr::filter(
      .data$condition == .env$condition_value,
      .data$time_h == .env$time_value
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
safe_neglog10 <- function(p) {
  p <- suppressWarnings(as.numeric(p)); p <- pmax(p, .Machine$double.xmin)
  pmin(-log10(p), 50)
}
median_ratio_log_expression <- function(counts) {
  counts <- as.matrix(counts)
  keep <- rowSums(is.finite(counts) & counts > 0) >= 2
  gm <- rep(NA_real_, nrow(counts))
  gm[keep] <- exp(rowMeans(log(pmax(counts[keep, , drop = FALSE], 1e-8)), na.rm = TRUE))
  ratios <- sweep(counts[is.finite(gm) & gm > 0, , drop = FALSE], 1, gm[is.finite(gm) & gm > 0], "/")
  sf <- apply(ratios, 2, function(z) median(z[is.finite(z) & z > 0], na.rm = TRUE))
  sf[!is.finite(sf) | sf <= 0] <- 1
  sf <- sf / exp(mean(log(sf)))
  log2(sweep(counts, 2, sf, "/") + 1)
}
permutation_fallback <- function(log_expr, group, note = "permutation fallback", nperm = 1999, seed = 42) {
  set.seed(seed)
  group <- factor(group, levels = c("A", "B"))
  A <- which(group == "A"); B <- which(group == "B")
  obs <- rowMeans(log_expr[, B, drop = FALSE], na.rm = TRUE) - rowMeans(log_expr[, A, drop = FALSE], na.rm = TRUE)
  n <- length(group)
  ge <- rep(0L, nrow(log_expr))
  for (i in seq_len(nperm)) {
    gp <- sample(group, n, replace = FALSE)
    Ap <- which(gp == "A"); Bp <- which(gp == "B")
    st <- rowMeans(log_expr[, Bp, drop = FALSE], na.rm = TRUE) - rowMeans(log_expr[, Ap, drop = FALSE], na.rm = TRUE)
    ge <- ge + as.integer(abs(st) >= abs(obs))
  }
  p <- (ge + 1) / (nperm + 1)
  tibble(gene_id = rownames(log_expr), pval = p, padj = p.adjust(p, "BH"), method_note = note)
}
read_precomputed_method <- function(dir, contrast_id, method) {
  if (is.null(dir) || isTRUE(dir) || !dir.exists(dir)) return(NULL)
  cand <- list.files(dir, pattern = paste0(contrast_id, ".*", method, ".*\\.(csv|tsv|txt)$"), full.names = TRUE, ignore.case = TRUE)
  if (!length(cand)) return(NULL)
  x <- read_table_auto(cand[[1]])
  gcol <- first_present(names(x), c("gene_id", "gene", "Gene", "name"))
  pcol <- first_present(names(x), c("pval", "pvalue", "P.Value", "PValue", "p_value"))
  qcol <- first_present(names(x), c("padj", "p_adj", "FDR", "fdr", "adj.P.Val", "qvalue", "q_value", "pfp"))
  ecol <- first_present(names(x), c("method_log2FC", "log2FC", "logFC", "signed_effect", "effect", "estimate"))
  scol <- first_present(names(x), c("signed_score", "harmonized_signed_score", "score", "signed_statistic"))
  dcol <- first_present(names(x), c("direction", "Direction", "call", "Call", "status", "regulation"))
  if (is.na(gcol)) stop("Precomputed ", method, " file needs gene column: ", cand[[1]])
  tibble(
    gene_id = as.character(x[[gcol]]),
    pval = if (!is.na(pcol)) suppressWarnings(as.numeric(x[[pcol]])) else NA_real_,
    padj = if (!is.na(qcol)) suppressWarnings(as.numeric(x[[qcol]])) else NA_real_,
    method_log2FC = if (!is.na(ecol)) suppressWarnings(as.numeric(x[[ecol]])) else NA_real_,
    precomputed_signed_score = if (!is.na(scol)) suppressWarnings(as.numeric(x[[scol]])) else NA_real_,
    direction = if (!is.na(dcol)) as.character(x[[dcol]]) else NA_character_,
    method_note = paste0(method, " precomputed: ", basename(cand[[1]]))
  )
}
collect_rows <- function(base_tbl, res_tbl, method, note_default) {
  out <- base_tbl %>% left_join(res_tbl, by = "gene_id")
  if (!"pval" %in% names(out)) out$pval <- NA_real_
  if (!"padj" %in% names(out)) out$padj <- NA_real_
  if (!"method_log2FC" %in% names(out)) out$method_log2FC <- NA_real_
  if (!"precomputed_signed_score" %in% names(out)) out$precomputed_signed_score <- NA_real_
  if (!"direction" %in% names(out)) out$direction <- NA_character_
  if (!"method_note" %in% names(out)) out$method_note <- note_default
  out %>% mutate(
    pval = suppressWarnings(as.numeric(.data$pval)),
    padj = suppressWarnings(as.numeric(.data$padj)),
    method_log2FC = suppressWarnings(as.numeric(.data$method_log2FC)),
    precomputed_signed_score = suppressWarnings(as.numeric(.data$precomputed_signed_score)),
    signed_effect = dplyr::coalesce(.data$method_log2FC, .data$common_effect),
    q_for_score = dplyr::coalesce(.data$padj, .data$pval),
    signed_score = dplyr::coalesce(.data$precomputed_signed_score, sign(.data$signed_effect) * safe_neglog10(.data$q_for_score)),
    method = method,
    method_note = ifelse(is.na(.data$method_note), note_default, .data$method_note)
  )
}
run_contrast <- function(counts_all, sample_sheet, contrasts, contrast_id, allow_fallbacks = FALSE, precomputed_dir = NULL, nperm = 1999) {
  cid <- stringr::str_trim(gsub("\ufeff", "", as.character(contrast_id)))
  contrasts <- contrasts %>%
    dplyr::mutate(
      contrast_id = stringr::str_trim(gsub("\ufeff", "", as.character(.data$contrast_id)))
    )
  cc <- contrasts %>% dplyr::filter(.data$contrast_id == .env$cid)
  if (nrow(cc) != 1) {
    message("Available contrast_ids:")
    message(paste(contrasts$contrast_id, collapse = "\n"))
    stop("Unknown or non-unique contrast_id: ", cid, " (matched n=", nrow(cc), ")", call. = FALSE)
  }

  selA <- unique(samples_for(sample_sheet, cc$denominator[[1]], cc$time_h[[1]]))
  selB <- unique(samples_for(sample_sheet, cc$numerator[[1]], cc$time_h[[1]]))
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
      ". denominator=", cc$denominator[[1]], " n=", length(selA),
      "; numerator=", cc$numerator[[1]], " n=", length(selB),
      call. = FALSE
    )
  }
  message(
    "Contrast ", cid,
    " | denominator ", cc$denominator[[1]], ": ", paste(selA, collapse = ", "),
    " | numerator ", cc$numerator[[1]], ": ", paste(selB, collapse = ", ")
  )

  X_expected <- counts_all[, c(selA, selB), drop = FALSE]
  X <- round(X_expected)
  X[X < 0] <- 0
  group <- factor(c(rep("A", length(selA)), rep("B", length(selB))), levels = c("A", "B"))
  coldata <- data.frame(sample = colnames(X), group = group)
  rownames(coldata) <- coldata$sample

  # Common harmonized effect for all method rows.
  require_or_stop("edgeR")
  y <- edgeR::DGEList(counts = X)
  y <- edgeR::calcNormFactors(y)
  logCPM <- edgeR::cpm(y, prior.count = 1, log = TRUE)
  baseline_expression <- rowMeans(logCPM[, group == "A", drop = FALSE], na.rm = TRUE)
  signed_effect <- rowMeans(logCPM[, group == "B", drop = FALSE], na.rm = TRUE) - baseline_expression
  base_tbl <- tibble(gene_id = rownames(X), contrast_id = cid,
                     baseline_expression = baseline_expression[rownames(X)],
                     common_effect = signed_effect[rownames(X)])

  status <- tibble(contrast_id = character(), method = character(), status = character(), note = character())
  add_status <- function(method, stat, note) {
    status <<- bind_rows(status, tibble(contrast_id = cid, method = method, status = stat, note = note))
  }
  rows <- list()

  # DESeq2
  if (requireNamespace("DESeq2", quietly = TRUE)) {
    res <- tryCatch({
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = X, colData = coldata, design = ~ group)
      dds <- dds[rowSums(DESeq2::counts(dds) >= 10) >= 2, ]
      dds <- DESeq2::DESeq(dds, quiet = TRUE)
      rr <- DESeq2::results(dds, contrast = c("group", "B", "A"))
      tibble(gene_id = rownames(rr), pval = rr$pvalue, padj = rr$padj,
             method_log2FC = rr$log2FoldChange,
             method_note = "DESeq2 native; rounded RSEM expected counts")
    }, error = function(e) { add_status("DESeq2", "error", conditionMessage(e)); NULL })
    if (!is.null(res)) { rows[["DESeq2"]] <- collect_rows(base_tbl, res, "DESeq2", "DESeq2 native"); add_status("DESeq2", "ok", "native") }
  } else add_status("DESeq2", "missing", "Package DESeq2 not installed")

  # edgeR
  if (requireNamespace("edgeR", quietly = TRUE)) {
    res <- tryCatch({
      y2 <- edgeR::DGEList(counts = X, group = group)
      keep <- rowSums(edgeR::cpm(y2) > 1) >= 2
      y2 <- y2[keep, , keep.lib.sizes = FALSE]
      y2 <- edgeR::calcNormFactors(y2)
      design <- model.matrix(~ group, data = coldata)
      y2 <- edgeR::estimateDisp(y2, design)
      fit <- edgeR::glmQLFit(y2, design)
      ql <- edgeR::glmQLFTest(fit, coef = 2)
      tt <- edgeR::topTags(ql, n = Inf)$table
      tibble(gene_id = rownames(tt), pval = tt$PValue, padj = tt$FDR,
             method_log2FC = tt$logFC,
             method_note = "edgeR QLF native; rounded RSEM expected counts")
    }, error = function(e) { add_status("edgeR", "error", conditionMessage(e)); NULL })
    if (!is.null(res)) { rows[["edgeR"]] <- collect_rows(base_tbl, res, "edgeR", "edgeR QLF native"); add_status("edgeR", "ok", "native") }
  }

  # limma/voom
  if (requireNamespace("limma", quietly = TRUE) && requireNamespace("edgeR", quietly = TRUE)) {
    res <- tryCatch({
      y3 <- edgeR::DGEList(counts = X)
      y3 <- edgeR::calcNormFactors(y3)
      design <- model.matrix(~ group, data = coldata)
      v <- limma::voom(y3, design, plot = FALSE)
      fit <- limma::lmFit(v, design)
      fit <- limma::eBayes(fit)
      tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
      tibble(gene_id = rownames(tt), pval = tt$P.Value, padj = tt$adj.P.Val,
             method_log2FC = tt$logFC,
             method_note = "limma/voom native; rounded RSEM expected counts")
    }, error = function(e) { add_status("limma/voom", "error", conditionMessage(e)); NULL })
    if (!is.null(res)) { rows[["limma/voom"]] <- collect_rows(base_tbl, res, "limma/voom", "limma/voom native"); add_status("limma/voom", "ok", "native") }
  } else add_status("limma/voom", "missing", "Package limma or edgeR not installed")

  # RankProd on logCPM.
  if (requireNamespace("RankProd", quietly = TRUE)) {
    res <- tryCatch({
      log_rank <- logCPM
      row_sd <- apply(log_rank, 1, stats::sd, na.rm = TRUE)
      keep <- is.finite(row_sd) & row_sd > 1e-4 & rowSums(log_rank > 0.5, na.rm = TRUE) >= 2
      set.seed(123)
      cl <- ifelse(group == "B", 1, 0)
      RP <- suppressWarnings(RankProd::RankProducts(log_rank[keep, , drop = FALSE], cl, logged = TRUE, na.rm = TRUE, rand = 123))
      tg <- suppressWarnings(RankProd::topGene(RP, cutoff = 1, method = "pfp"))
      grab <- function(Ta) {
        if (is.null(Ta) || !NROW(Ta)) return(NULL)
        df <- as.data.frame(Ta)
        gcol <- intersect(names(df), c("GeneName", "Gene", "gene", "Symbol"))
        gene <- if (length(gcol)) as.character(df[[gcol[[1]]]]) else {
          idx <- suppressWarnings(as.integer(rownames(df)))
          if (all(!is.na(idx))) rownames(log_rank[keep, , drop = FALSE])[idx] else rownames(df)
        }
        pcol <- intersect(names(df), c("pfp", "PFP", "pfp.up", "pfp.down"))
        pfp <- if (length(pcol)) suppressWarnings(as.numeric(df[[pcol[[1]]]])) else NA_real_
        tibble(gene_id = gene, padj = pfp)
      }
      U <- bind_rows(grab(tg$Table1), grab(tg$Table2))
      if (!nrow(U)) stop("RankProd returned empty tables")
      U %>% group_by(.data$gene_id) %>% summarise(padj = min(.data$padj, na.rm = TRUE), .groups = "drop") %>%
        mutate(pval = NA_real_, method_log2FC = NA_real_, method_note = "RankProd native PFP")
    }, error = function(e) { add_status("RankProd", "error", conditionMessage(e)); NULL })
    if (!is.null(res)) { rows[["RankProd"]] <- collect_rows(base_tbl, res, "RankProd", "RankProd native PFP"); add_status("RankProd", "ok", "native") }
  } else add_status("RankProd", "missing", "Package RankProd not installed")

  # PenDA: precomputed > native is difficult/version-specific > fallback if explicitly allowed.
  pre_penda <- read_precomputed_method(precomputed_dir, cid, "PenDA")
  if (!is.null(pre_penda)) {
    rows[["PenDA"]] <- collect_rows(base_tbl, pre_penda, "PenDA", "PenDA precomputed")
    add_status("PenDA", "ok", "precomputed")
  } else if (requireNamespace("penda", quietly = TRUE)) {
    # Native PenDA APIs differ across versions. To avoid silent incorrect calls, require precomputed native output.
    add_status("PenDA", "available_not_run", "Package penda installed, but this builder requires precomputed PenDA output unless --allow-fallbacks is used.")
  } else if (allow_fallbacks) {
    A <- logCPM[, group == "A", drop = FALSE]
    B <- logCPM[, group == "B", drop = FALSE]
    qa <- t(apply(A, 1, quantile, probs = c(0.25, 0.75), na.rm = TRUE))
    indB <- (B < qa[, 1]) | (B > qa[, 2])
    k <- rowSums(indB, na.rm = TRUE)
    n_cases <- ncol(B)
    ctrl_out <- rowMeans((A < qa[, 1]) | (A > qa[, 2]), na.rm = TRUE)
    p0 <- pmin(pmax(ctrl_out, 0.01), 0.20)
    pval <- stats::pbinom(q = pmax(k, 0) - 1, size = n_cases, prob = p0, lower.tail = FALSE)
    res <- tibble(gene_id = rownames(A), pval = pval, padj = p.adjust(pval, "BH"), method_log2FC = NA_real_, method_note = "PenDA IQR/binomial fallback; not native PenDA")
    rows[["PenDA"]] <- collect_rows(base_tbl, res, "PenDA", "PenDA fallback")
    add_status("PenDA", "fallback", "IQR/binomial fallback used because native/precomputed PenDA output was unavailable")
  } else add_status("PenDA", "missing", "No precomputed PenDA output and --allow-fallbacks not set")

  # RankCompV3 (REO-family): precomputed preferred; permutation fallback only if explicitly allowed for layout QA.
  # Older files may have been written into a CellComp compatibility slot; accept them but label the method RankCompV3.
  pre_rank <- read_precomputed_method(precomputed_dir, cid, "RankCompV3")
  rank_note <- "RankCompV3 precomputed"
  if (is.null(pre_rank)) {
    pre_rank <- read_precomputed_method(precomputed_dir, cid, "CellComp")
    rank_note <- "RankCompV3/REO imported from its archived precomputed result table"
  }
  if (!is.null(pre_rank)) {
    rows[["RankCompV3"]] <- collect_rows(base_tbl, pre_rank, "RankCompV3", rank_note)
    add_status("RankCompV3", "ok", rank_note)
  } else if (allow_fallbacks) {
    res <- permutation_fallback(logCPM, group, note = "RankCompV3/REO permutation fallback; not native RankCompV3", nperm = nperm) %>% mutate(method_log2FC = NA_real_)
    rows[["RankCompV3"]] <- collect_rows(base_tbl, res, "RankCompV3", "RankCompV3 fallback")
    add_status("RankCompV3", "fallback", "Permutation fallback used because native/precomputed RankCompV3 output was unavailable")
  } else add_status("RankCompV3", "missing", "No precomputed RankCompV3 output and --allow-fallbacks not set")

  expected <- c("DESeq2", "edgeR", "limma/voom", "RankProd", "PenDA", "RankCompV3")
  have <- names(rows)
  missing <- setdiff(expected, have)
  if (length(missing) && !allow_fallbacks) {
    message("Missing methods for ", cid, ": ", paste(missing, collapse = ", "))
  }
  list(results = bind_rows(rows), status = status)
}

main <- function() {
  args <- parse_args()
  root <- repo_root()
  counts_path <- coalesce_arg(args[["counts"]], file.path(root, "data/processed/raw_counts_rsemgenes.tsv"))
  sample_sheet_path <- coalesce_arg(args[["sample-sheet"]], file.path(root, "config/sample_sheet.csv"))
  contrasts_path <- coalesce_arg(args[["contrasts"]], file.path(root, "config/contrasts.csv"))
  out_path <- coalesce_arg(args[["out"]], file.path(root, "figure3/differential_expression/source_data/unified_de_results.csv"))
  status_path <- coalesce_arg(args[["status-out"]], file.path(root, "figure3/differential_expression/source_data/unified_de_method_status.tsv"))
  precomputed_dir <- coalesce_arg(args[["precomputed-dir"]], file.path(root, "figure3/differential_expression/source_data/precomputed_method_outputs"))
  contrast_ids <- unlist(strsplit(coalesce_arg(args[["contrast-ids"]], "DT_vs_D_4h,DT_vs_D_24h"), ","))
  contrast_ids <- contrast_ids[nzchar(contrast_ids)]
  allow_fallbacks <- isTRUE(args[["allow-fallbacks"]]) || tolower(coalesce_arg(args[["allow-fallbacks"]], "false")) %in% c("true", "1", "yes")
  nperm <- as.integer(coalesce_arg(args[["fallback-permutations"]], "1999"))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  counts <- read_counts_matrix(counts_path)
  # Manuscript low-count filter already recorded in config; apply it here too.
  keep <- rowSums(counts, na.rm = TRUE) >= as.numeric(coalesce_arg(args[["min-total-count"]], "80"))
  counts <- counts[keep, , drop = FALSE]
  sample_sheet <- condition_from_sample_sheet(readr::read_csv(sample_sheet_path, show_col_types = FALSE)) %>%
    dplyr::distinct(.data$sample_id, .keep_all = TRUE)
  contrasts <- readr::read_csv(contrasts_path, show_col_types = FALSE) %>%
    dplyr::mutate(
      contrast_id = stringr::str_trim(gsub("\ufeff", "", as.character(.data$contrast_id)))
    )
  contrast_ids <- stringr::str_trim(gsub("\ufeff", "", as.character(contrast_ids)))
  contrast_ids <- contrast_ids[nzchar(contrast_ids)]
  message("Loaded contrast IDs:")
  message(paste(contrasts$contrast_id, collapse = "\n"))

  runs <- purrr::map(contrast_ids, ~ run_contrast(counts, sample_sheet, contrasts, .x, allow_fallbacks = allow_fallbacks, precomputed_dir = precomputed_dir, nperm = nperm))
  results <- bind_rows(purrr::map(runs, "results"))
  status <- bind_rows(purrr::map(runs, "status"))
  readr::write_csv(results, out_path)
  readr::write_tsv(status, status_path)
  message("Wrote unified DE table: ", out_path)
  message("Wrote method-status table: ", status_path)
  if (any(status$status %in% c("missing", "available_not_run", "error"))) {
    warning("Some methods were missing/not run. Inspect method-status table before using this for the manuscript.")
  }
}

main()
