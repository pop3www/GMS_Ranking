#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales) })
SCRIPT_VERSION <- "2026-06-12 v39 (rank bands from expression matrix; Fig2 uses normalized RSEM expected counts)"
message("[06] script version: ", SCRIPT_VERSION)
args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(argv) { out <- list(); i <- 1L; while (i <= length(argv)) { a <- argv[[i]]; if (startsWith(a, "--")) { a2 <- sub("^--", "", a); if (grepl("=", a2, fixed = TRUE)) { sp <- strsplit(a2, "=", fixed = TRUE)[[1]]; key <- sp[[1]]; val <- paste(sp[-1], collapse = "=") } else { key <- a2; if (i < length(argv) && !startsWith(argv[[i + 1L]], "--")) { val <- argv[[i + 1L]]; i <- i + 1L } else val <- "TRUE" }; key <- tolower(gsub("_", "-", key)); out[[key]] <- val }; i <- i + 1L }; out }
opt <- parse_args(args)
get_arg <- function(keys, default = NA_character_, env = NULL, required = FALSE) { keys <- tolower(gsub("_", "-", keys)); for (k in keys) if (!is.null(opt[[k]]) && nzchar(opt[[k]])) return(opt[[k]]); if (!is.null(env)) { ev <- Sys.getenv(env, unset = ""); if (nzchar(ev)) return(ev) }; if (required) stop("[06] Missing required arg: --", keys[[1]], call. = FALSE); default }
get_num <- function(keys, default, env = NULL) { raw <- get_arg(keys, default = as.character(default), env = env); val <- suppressWarnings(as.numeric(raw)); if (!is.finite(val)) default else val }
strip_quotes <- function(x) gsub('^"|"$', '', x)
normalize_rx <- function(rx) { if (is.na(rx) || !nzchar(rx)) return(rx); gsub("\\\\", "\\", rx, fixed = TRUE) }
extract_ensg <- function(x) { x <- as.character(x); out <- rep(NA_character_, length(x)); has <- !is.na(x) & grepl("ENSG[0-9]+", x, perl = TRUE); out[has] <- sub("^.*?(ENSG[0-9]+).*$", "\\1", x[has], perl = TRUE); out }
expr_path <- get_arg(c("expr", "tpm", "counts"), required = TRUE)
id_col <- get_arg("id-col", default = "gene_id")
baseline_pattern <- get_arg("baseline-pattern", required = TRUE)
groupA_pattern <- get_arg("groupA-pattern", required = TRUE)
groupB_pattern <- get_arg("groupB-pattern", required = TRUE)
head_frac <- get_num("head-frac", default = 0.10, env = "HEAD_FRAC")
tail_frac <- get_num("tail-frac", default = 0.10, env = "TAIL_FRAC")
lfc_thresh <- get_num("lfc-thresh", default = 0.32, env = "LFC_THRESH")
fdr_thresh <- get_num("fdr-thresh", default = 0.05, env = "FDR_THRESH")
input_scale <- tolower(get_arg(c("input-scale", "expr-scale", "scale"), default = Sys.getenv("EXPR_SCALE", unset = "normcounts")))
de_method <- tolower(get_arg("de-method", default = Sys.getenv("DE_METHOD", unset = "none")))
out_dir <- get_arg("out", default = "out/fimo_summary")
if (!file.exists(expr_path)) stop("[06] expression matrix not found: ", expr_path, call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("[06] Reading expression matrix: ", expr_path); message("[06] input scale: ", input_scale)
expr <- fread(expr_path, sep = "\t", header = TRUE, check.names = FALSE, data.table = TRUE, quote = "\"")
setnames(expr, names(expr), strip_quotes(names(expr)))
if (!(id_col %in% names(expr))) stop("[06] id-col not found. --id-col=", id_col, "; available: ", paste(head(names(expr), 8), collapse = ", "), call. = FALSE)
setnames(expr, id_col, "gene"); expr[, gene := strip_quotes(as.character(gene))]; expr[, ensg := extract_ensg(gene)]
num_cols <- setdiff(names(expr), c("gene", "ensg")); expr[, (num_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols = num_cols]
rx_base <- normalize_rx(baseline_pattern); rx_A <- normalize_rx(groupA_pattern); rx_B <- normalize_rx(groupB_pattern)
message("[06] Regex after normalization:"); message("  baseline: ", rx_base); message("  groupA:   ", rx_A); message("  groupB:   ", rx_B)
ix <- function(rx) if (is.na(rx) || !nzchar(rx)) rep(FALSE, length(num_cols)) else grepl(rx, num_cols, perl = TRUE)
ix_base <- ix(rx_base); ix_A <- ix(rx_A); ix_B <- ix(rx_B)
cat("[06] Column matches \n", "  baseline (", sum(ix_base), "): ", paste(num_cols[ix_base], collapse = ", "), "\n", "  groupA   (", sum(ix_A), "): ", paste(num_cols[ix_A], collapse = ", "), "\n", "  groupB   (", sum(ix_B), "): ", paste(num_cols[ix_B], collapse = ", "), "\n", sep = "")
if (sum(ix_base) == 0 || sum(ix_A) == 0 || sum(ix_B) == 0) stop("[06] No columns matched baseline/groupA/groupB regex patterns.", call. = FALSE)
cn_base <- num_cols[ix_base]; cn_A <- num_cols[ix_A]; cn_B <- num_cols[ix_B]
to_mean <- function(dt, cols) { m <- as.matrix(dt[, ..cols]); storage.mode(m) <- "double"; rowMeans(m, na.rm = TRUE) }
expr$baseline_expr <- to_mean(expr, cn_base); expr$A_expr <- to_mean(expr, cn_A); expr$B_expr <- to_mean(expr, cn_B)
expr$baseline_tpm <- expr$baseline_expr; expr$A_tpm <- expr$A_expr; expr$B_tpm <- expr$B_expr
lfc_mean <- log2(expr$A_expr + 1) - log2(expr$B_expr + 1); p_value <- rep(NA_real_, nrow(expr)); fdr <- rep(NA_real_, nrow(expr)); lfc_use <- lfc_mean
if (de_method == "limma") { if (!requireNamespace("limma", quietly = TRUE)) stop("[06] --de-method limma requested, but package 'limma' is not installed.", call. = FALSE); message("[06] Running limma on log2(expression+1) for DE stats ..."); cols_de <- c(cn_A, cn_B); mat <- as.matrix(expr[, ..cols_de]); storage.mode(mat) <- "double"; mat <- log2(mat + 1); rownames(mat) <- make.unique(as.character(expr$gene)); group <- factor(c(rep("A", length(cn_A)), rep("B", length(cn_B))), levels = c("B", "A")); design <- model.matrix(~0 + group); colnames(design) <- levels(group); fit <- limma::lmFit(mat, design); cont <- limma::makeContrasts(A_vs_B = A - B, levels = design); fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont)); tt <- limma::topTable(fit2, coef = "A_vs_B", number = Inf, sort.by = "none"); lfc_use <- tt$logFC; p_value <- tt$P.Value; fdr <- tt$adj.P.Val; message("[06] limma done. Using limma logFC and BH-FDR for Up/Down calls.") }
expr$lfc <- lfc_use; expr$p_value <- p_value; expr$fdr <- fdr
ord <- order(-expr$baseline_expr, expr$gene); expr <- expr[ord]; expr$rank <- seq_len(nrow(expr)); N <- nrow(expr)
head_n <- max(1L, as.integer(ceiling(head_frac * N))); tail_n <- max(1L, as.integer(ceiling(tail_frac * N)))
if (head_n + tail_n >= N) { head_n <- max(1L, floor(N / 3)); tail_n <- max(1L, floor(N / 3)); warning("[06] head/tail fractions too large; adjusted to avoid overlap.") }
band <- rep("mid", N); band[seq_len(head_n)] <- "head"; band[seq.int(N - tail_n + 1L, N)] <- "tail"; expr$band <- factor(band, levels = c("head", "mid", "tail"))
if (de_method == "limma") expr$updown <- ifelse(!is.na(expr$fdr) & expr$fdr <= fdr_thresh & expr$lfc >= lfc_thresh, "Up", ifelse(!is.na(expr$fdr) & expr$fdr <= fdr_thresh & expr$lfc <= -lfc_thresh, "Down", "Neutral")) else expr$updown <- ifelse(expr$lfc >= lfc_thresh, "Up", ifelse(expr$lfc <= -lfc_thresh, "Down", "Neutral"))
bands <- expr[, .(gene, ensg, baseline_expr, baseline_tpm, A_expr, B_expr, lfc, rank, band, updown, p_value, fdr, expr_scale = input_scale)]
fwrite(bands, file.path(out_dir, "rank_bands.csv"))
cnt <- bands[, .N, keyby = .(band, updown)]; fwrite(cnt, file.path(out_dir, "rank_band_counts.tsv"), sep = "\t")
p_counts <- ggplot(cnt, aes(x = band, y = N, fill = updown)) + geom_col(position = "stack") + scale_y_continuous(labels = comma) + scale_fill_brewer(palette = "Set2") + theme_minimal(base_size = 11) + labs(title = "Band counts by Up/Down", x = NULL, y = "Genes")
ggsave(file.path(out_dir, "rank_bands_qc.png"), p_counts, width = 6, height = 4, dpi = 300)
message("[06] Wrote: ", file.path(out_dir, "rank_bands.csv"))
