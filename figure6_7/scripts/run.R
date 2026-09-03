#!/usr/bin/env Rscript

## ============================================================
## Figures 7-8: program-level reweighting bundle
## Patch level: pubready v10 text-alignment + annotation polish
##
## Production workflow for Figures 6-7 and Figure S3.
##
## Key decisions in this bundle:
## - Uses validated RSEM expected counts by default for DESeq2::vst
##   (or an explicitly supplied precomputed VST matrix from core)
## - Refuses TPM/CPM-like matrices and summarized means in the count path
## - Uses sample_sheet.csv as metadata source of truth, including the compatibility
##   schema supported in this repository: sample_id/time_h/dox/myc/cpt_level
## - Keeps parser-derived metadata only as a QC cross-check
## - Applies the manuscript count filter: genes with >= 80 total counts
## - Uses the manuscript module workflow: PC loading + NRI features,
##   PAM clustering, silhouette selection over k = 3:8, and top genes by
##   the magnitude of the surrogate temporal shift
## - Anchors program-level baseline ranking to shared 4 h CTRL using
##   size-factor-normalized count-like expression when counts are supplied
## - Uses mean VST differences for Fig. 8 fgsea ranks, matching Methods
## - Uses the manuscript distortion-associated PC definition, oriented
##   so positive loadings track the later / more distorted state
## - Produces the manuscript-aligned 4-panel Figure 7 by default, plus a
##   simplified A/B/D alternate composite for readability review
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(irlba)
  library(FNN)
  library(cluster)
  library(RobustRankAggreg)
  library(pheatmap)
  library(fgsea)
  library(msigdbr)
  library(grid)
  library(gridExtra)
  library(png)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

## ----------------------------
## Helpers
## ----------------------------
stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

required_packages <- c(
  "tidyverse","DESeq2","irlba","FNN","cluster",
  "RobustRankAggreg","pheatmap","fgsea","msigdbr","gridExtra","png",
  "AnnotationDbi","clusterProfiler","org.Hs.eg.db"
)

check_required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing)) {
    stopf("Missing required packages: %s\nRun install_packages.R first.",
          paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

norm_key <- function(s) gsub("[^A-Za-z0-9]+", "_", toupper(s))

normalize_token <- function(s) {
  s <- tolower(s)
  s <- gsub("[^a-z0-9]+", "_", s)
  s <- gsub("_+", "_", s)
  s <- gsub("^_|_$", "", s)
  s
}

canonicalize_time <- function(x) {
  x0 <- normalize_token(x)
  case_when(
    x0 %in% c("4", "4h", "4_hr", "4_hrs", "4_hour", "4_hours") ~ "4 hours",
    x0 %in% c("24", "24h", "24_hr", "24_hrs", "24_hour", "24_hours") ~ "24 hours",
    grepl("^4($|_)", x0) ~ "4 hours",
    grepl("^24($|_)", x0) ~ "24 hours",
    TRUE ~ as.character(x)
  )
}

canonicalize_induction <- function(x) {
  x0 <- normalize_token(x)
  case_when(
    x0 %in% c("ctrl", "control") ~ "CTRL",
    x0 %in% c("dox", "d", "doxy") ~ "DOX",
    x0 %in% c("tam", "tamox") ~ "TAM",
    x0 %in% c("dt", "dox_tam", "dox_tamox", "dox_tamoxifen", "dox_tamoxifen_treatment",
              "dox_tam_combo", "dox_tam_combination", "dox_tam_combined", "dox_tamox_combined",
              "dox_tamox_combo", "dox_plus_tam", "dox_tam_plus", "dox_tam_treatment",
              "dox_tamox_treatment", "dox_tamoxifen_combo", "dox_tamoxifen_combination") ~ "DOX+TAM",
    grepl("^dt($|_)", x0) ~ "DOX+TAM",
    grepl("^dox.*tam|^tam.*dox", x0) ~ "DOX+TAM",
    TRUE ~ toupper(as.character(x))
  )
}

canonicalize_stress <- function(x) {
  x0 <- normalize_token(x)
  case_when(
    x0 %in% c("nocpt", "no_cpt", "none", "unstressed", "ctrl", "na", "") ~ "NoCPT",
    x0 %in% c("cpt") ~ "CPT",
    x0 %in% c("cpt_low", "low_cpt", "lcpt", "l_cpt") ~ "CPT_low",
    x0 %in% c("cpt_high", "high_cpt", "hcpt", "h_cpt") ~ "CPT_high",
    grepl("hcpt|high_cpt|cpt_high", x0) ~ "CPT_high",
    grepl("lcpt|low_cpt|cpt_low", x0) ~ "CPT_low",
    grepl("cpt", x0) ~ "CPT",
    TRUE ~ as.character(x)
  )
}

short_time <- function(t) {
  ifelse(t == "4 hours", "4h", ifelse(t == "24 hours", "24h", t))
}

short_stress <- function(s) {
  case_when(
    is.na(s) ~ NA_character_,
    s == "NoCPT" ~ "NoCPT",
    s == "CPT_low" ~ "LCPT",
    s == "CPT_high" ~ "HCPT",
    s == "CPT" ~ "CPT",
    TRUE ~ as.character(s)
  )
}

resolve_time_levels <- function(x) {
  lv <- unique(as.character(x))
  lv <- lv[!is.na(lv)]
  ord <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", lv)))
  if (all(!is.na(ord))) lv[order(ord)] else sort(lv)
}

baseline_ranks01 <- function(vec) {
  r <- rank(-vec, ties.method = "average")
  r / max(r)
}

flatten_fgsea_for_csv <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(df)
  out <- as_tibble(df)
  if ("leadingEdge" %in% names(out)) {
    out <- out |> mutate(
      leadingEdge = vapply(
        leadingEdge,
        function(x) paste(as.character(x), collapse = ";"),
        FUN.VALUE = character(1)
      )
    )
  }
  out
}

collapse_named_stat <- function(values, genes, mode = c("maxabs", "mean")) {
  mode <- match.arg(mode)
  tb <- tibble(gene = as.character(genes), value = as.numeric(values)) |>
    filter(!is.na(gene), gene != "", is.finite(value))
  if (!nrow(tb)) return(numeric(0))
  if (mode == "maxabs") {
    out <- tb |>
      group_by(gene) |>
      slice_max(order_by = abs(value), n = 1, with_ties = FALSE) |>
      ungroup()
  } else {
    out <- tb |>
      group_by(gene) |>
      summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  }
  ranks <- out$value
  names(ranks) <- out$gene
  ranks
}

make_progress <- function(total) {
  i <- 0L
  t0 <- Sys.time()
  pb <- utils::txtProgressBar(0, total, style = 3)
  step <- function(label) {
    i <<- i + 1L
    utils::setTxtProgressBar(pb, i)
    cat(sprintf("\n[%s] Step %02d/%02d | %s\n",
                format(Sys.time(), "%H:%M"), i, total, label))
    flush.console()
  }
  done <- function() {
    utils::setTxtProgressBar(pb, total)
    close(pb)
    cat(sprintf("[Done] in %.1fs\n",
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  list(step = step, done = done)
}

png_to_svg <- function(path_png) {
  if (is.na(path_png) || !nzchar(path_png)) return(path_png)
  if (grepl("\\.[Pp][Nn][Gg]$", path_png)) {
    sub("\\.[Pp][Nn][Gg]$", ".svg", path_png)
  } else {
    paste0(path_png, ".svg")
  }
}

png_to_pdf <- function(path_png) {
  if (is.na(path_png) || !nzchar(path_png)) return(path_png)
  if (grepl("\\.[Pp][Nn][Gg]$", path_png)) {
    sub("\\.[Pp][Nn][Gg]$", ".pdf", path_png)
  } else {
    paste0(path_png, ".pdf")
  }
}

draw_grid_object_with_padding <- function(g, padding = NULL) {
  if (is.null(padding)) padding <- if (!is.null(CONFIG$COMPOSITE_CANVAS_PADDING)) CONFIG$COMPOSITE_CANVAS_PADDING else 0.012
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    x = 0.5, y = 0.5,
    width = grid::unit(1 - 2 * padding, "npc"),
    height = grid::unit(1 - 2 * padding, "npc"),
    clip = "off"
  ))
  grid::grid.draw(g)
  grid::popViewport()
}

open_svg_device <- function(filename, width, height) {
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file = filename, width = width, height = height)
  } else {
    grDevices::svg(filename = filename, width = width, height = height)
  }
}

save_gtable_png_svg <- function(gtable, out_png, width, height, res = NULL) {
  if (is.null(res)) res <- CONFIG$PLOT_DPI
  dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
  grDevices::png(filename = out_png, width = width, height = height, units = "in", res = res)
  draw_grid_object_with_padding(gtable, padding = 0.010)
  grDevices::dev.off()
  saveRDS(gtable, sub("[.]png$", ".grob.rds", out_png))
  if (isTRUE(CONFIG$WRITE_SVG)) {
    out_svg <- png_to_svg(out_png)
    open_svg_device(out_svg, width = width, height = height)
    draw_grid_object_with_padding(gtable, padding = 0.012)
    grDevices::dev.off()
  }
  if (isTRUE(CONFIG$WRITE_PDF)) {
    out_pdf <- png_to_pdf(out_png)
    tryCatch({
      grDevices::pdf(file = out_pdf, width = width, height = height, useDingbats = FALSE)
      draw_grid_object_with_padding(gtable, padding = 0.010)
      grDevices::dev.off()
    }, error = function(e) {
      message(sprintf("[Plot] PDF export failed for %s: %s", basename(out_png), conditionMessage(e)))
      try(grDevices::dev.off(), silent = TRUE)
    })
  }
  invisible(TRUE)
}

save_grob_png_svg <- function(grob, out_png, width, height, res = NULL) {
  if (is.null(res)) res <- CONFIG$PLOT_DPI
  dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
  grDevices::png(filename = out_png, width = width, height = height, units = "in", res = res)
  draw_grid_object_with_padding(grob, padding = 0.010)
  grDevices::dev.off()
  saveRDS(grob, sub("[.]png$", ".grob.rds", out_png))
  if (isTRUE(CONFIG$WRITE_SVG)) {
    out_svg <- png_to_svg(out_png)
    open_svg_device(out_svg, width = width, height = height)
    draw_grid_object_with_padding(grob, padding = 0.012)
    grDevices::dev.off()
  }
  if (isTRUE(CONFIG$WRITE_PDF)) {
    out_pdf <- png_to_pdf(out_png)
    tryCatch({
      grDevices::pdf(file = out_pdf, width = width, height = height, useDingbats = FALSE)
      draw_grid_object_with_padding(grob, padding = 0.010)
      grDevices::dev.off()
    }, error = function(e) {
      message(sprintf("[Plot] PDF export failed for %s: %s", basename(out_png), conditionMessage(e)))
      try(grDevices::dev.off(), silent = TRUE)
    })
  }
  invisible(TRUE)
}

save_ggplot_png_svg <- function(p, out_png, width, height, dpi = NULL) {
  if (is.null(dpi)) dpi <- CONFIG$PLOT_DPI
  dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(filename = out_png, plot = p, width = width, height = height, dpi = dpi, limitsize = FALSE)
  saveRDS(p, sub("[.]png$", ".plot.rds", out_png))
  if (isTRUE(CONFIG$WRITE_SVG)) {
    out_svg <- png_to_svg(out_png)
    if (requireNamespace("svglite", quietly = TRUE)) {
      ggplot2::ggsave(filename = out_svg, plot = p, width = width, height = height, device = svglite::svglite, limitsize = FALSE)
    } else {
      ggplot2::ggsave(filename = out_svg, plot = p, width = width, height = height, limitsize = FALSE)
    }
  }
  if (isTRUE(CONFIG$WRITE_PDF)) {
    out_pdf <- png_to_pdf(out_png)
    tryCatch(
      ggplot2::ggsave(filename = out_pdf, plot = p, width = width, height = height, device = grDevices::pdf, limitsize = FALSE),
      error = function(e) message(sprintf("[Plot] PDF export failed for %s: %s", basename(out_png), conditionMessage(e)))
    )
  }
  invisible(TRUE)
}

## ----------------------------
## Config
## ----------------------------

## ----------------------------
## Repository / path resolution
## ----------------------------
get_script_path <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- ca[grep("^--file=", ca)]
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

is_abs_path <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\])", path)
}

SCRIPT_PATH <- get_script_path()
SCRIPT_DIR <- dirname(SCRIPT_PATH)
MODULE_DIR <- if (basename(SCRIPT_DIR) == "scripts") dirname(SCRIPT_DIR) else SCRIPT_DIR
REPO_ROOT <- dirname(MODULE_DIR)

first_existing_path <- function(paths) {
  for (p in unique(paths)) {
    if (!is.na(p) && nzchar(p) && file.exists(p)) {
      return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }
  ""
}

resolve_input_path <- function(path, fallback_paths = character()) {
  auto_tokens <- c("", "AUTO", "AUTO_COUNTS", "AUTO_SAMPLE_SHEET", "AUTO_VST")
  if (is.null(path) || path %in% auto_tokens) {
    resolved <- first_existing_path(fallback_paths)
    return(if (nzchar(resolved)) resolved else fallback_paths[1])
  }
  if (is_abs_path(path)) return(normalizePath(path, winslash = "/", mustWork = FALSE))
  cands <- c(
    path,
    file.path(getwd(), path),
    file.path(MODULE_DIR, path),
    file.path(REPO_ROOT, path)
  )
  resolved <- first_existing_path(cands)
  if (nzchar(resolved)) return(resolved)
  normalizePath(file.path(MODULE_DIR, path), winslash = "/", mustWork = FALSE)
}

resolve_output_dir <- function(path) {
  if (is.null(path) || !nzchar(path)) path <- "outputs"
  if (is_abs_path(path)) return(normalizePath(path, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(MODULE_DIR, path), winslash = "/", mustWork = FALSE)
}

finalize_config_paths <- function() {
  CONFIG$COUNTS_CSV <<- resolve_input_path(
    CONFIG$COUNTS_CSV,
    c(
      file.path(REPO_ROOT, "data", "processed", "raw_counts_rsemgenes.tsv"),
      file.path(REPO_ROOT, "data", "processed", "RawCountFile_rsemgenes.tsv"),
      file.path(REPO_ROOT, "core", "outputs", "counts_replicate_level.csv"),
      file.path(REPO_ROOT, "core", "outputs", "counts_filtered.tsv.gz"),
      file.path(REPO_ROOT, "core", "outputs", "counts_filtered.csv"),
      file.path(REPO_ROOT, "data", "processed", "counts_replicate_level.csv"),
      file.path(REPO_ROOT, "data", "raw", "RawCountFile_rsemgenes.txt"),
      file.path(REPO_ROOT, "data", "raw", "rsem", "RawCountFile_rsemgenes.txt"),
      file.path(MODULE_DIR, "source_data", "raw_counts_rsemgenes.tsv"),
      file.path(MODULE_DIR, "source_data", "counts_replicate_level.csv")
    )
  )
  CONFIG$VST_CSV <<- resolve_input_path(
    CONFIG$VST_CSV,
    c(
      file.path(REPO_ROOT, "core", "outputs", "vst_matrix.tsv.gz"),
      file.path(REPO_ROOT, "core", "outputs", "vst_matrix.csv"),
      file.path(MODULE_DIR, "source_data", "vst_matrix.tsv.gz"),
      file.path(MODULE_DIR, "source_data", "vst_matrix.csv")
    )
  )
  CONFIG$SAMPLE_SHEET_CSV <<- resolve_input_path(
    CONFIG$SAMPLE_SHEET_CSV,
    c(
      file.path(REPO_ROOT, "config", "sample_sheet.csv"),
      file.path(MODULE_DIR, "source_data", "sample_sheet.csv")
    )
  )
  CONFIG$OUT_DIR <<- resolve_output_dir(CONFIG$OUT_DIR)
  invisible(CONFIG)
}

CONFIG <- list(
  COUNTS_CSV = "AUTO_COUNTS",
  VST_CSV = "AUTO_VST",
  SAMPLE_SHEET_CSV = "AUTO_SAMPLE_SHEET",
  GENE_COL = "gene_id",
  GENE_ID_TYPE = "auto",
  COUNT_INPUT_MODE = "auto",  # auto | integer_counts | rsem_expected_counts
  OUT_DIR = "outputs",
  GENES_PER_MODULE_HEATMAP = 100,
  MAX_GENES_PER_SET_UNION = 6,
  FIG7B_TOP_TERMS = 1,
  FIG7_COMPOSITE_WIDTH = 25.5,
  FIG7_COMPOSITE_HEIGHT = 29.6,
  FIG8_COMPOSITE_WIDTH = 21.2,
  FIG8_COMPOSITE_HEIGHT = 16.6,
  PLOT_DPI = 450,
  COMPOSITE_PANEL_LETTER_SIZE = 30,
  SUPP3_COMPOSITE_WIDTH = 18.6,
  SUPP3_COMPOSITE_HEIGHT = 8.8,
  MIN_TOTAL_COUNTS = 80,
  NCOMP_PCA = 10,
  K_MODULES = 7,
  MODULE_SEED = 1771,
  MODULE_K_MIN = 3,
  MODULE_K_MAX = 8,
  REQUIRE_SEVEN_MODULES = TRUE,
  MODULE_K_RULE = "fixed_7_with_silhouette_audit",
  MODULE_CLUSTER_ENGINE = "clara",  # clara | exact_pam | auto
  MODULE_EXACT_PAM_MAX_N = 5000,
  MODULE_SILHOUETTE_SAMPLE_N = 4000,
  MODULE_CLARA_SAMPLES = 20,
  MODULE_CLARA_SAMPSIZE = 5000,
  HEAD_FRAC = 0.10,
  FGSEA_NPERM = 10000,
  FGSEA_MIN_SIZE = 10,
  FGSEA_MAX_SIZE = 500,
  RUN_FIG7 = TRUE,
  FIG7_INCLUDE_PANEL_C = TRUE,
  RUN_FIG8 = TRUE,
  RUN_MITO_CONTROL = TRUE,
  META_ONLY = FALSE,
  WRITE_SVG = TRUE,
  WRITE_PDF = TRUE,
  COMPOSITE_CANVAS_PADDING = 0.030,
  AXIS_MODE = "distortion_associated_time_correlated_pc",
  BASELINE_RANK_SOURCE = "size_factor_normalized_counts_when_available",
  REUSE_MODULE_TABLES = FALSE,
  REUSE_FGSEA_TABLES = FALSE
)

parse_cli_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) return(invisible(NULL))
  for (a in args) {
    if (!grepl("^--", a)) next
    kv <- sub("^--", "", a)
    key <- sub("=.*$", "", kv)
    val <- sub("^[^=]*=", "", kv)
    if (identical(key, val)) {
      if (key %in% c("no-fig7", "no_fig7")) CONFIG$RUN_FIG7 <<- FALSE
      if (key %in% c("no-fig8", "no_fig8")) CONFIG$RUN_FIG8 <<- FALSE
      if (key %in% c("no-mito", "no_mito")) CONFIG$RUN_MITO_CONTROL <<- FALSE
      if (key %in% c("meta-only", "meta_only", "meta")) CONFIG$META_ONLY <<- TRUE
      if (key %in% c("no-svg", "no_svg")) CONFIG$WRITE_SVG <<- FALSE
      if (key %in% c("svg", "write-svg", "write_svg")) CONFIG$WRITE_SVG <<- TRUE
      if (key %in% c("no-pdf", "no_pdf")) CONFIG$WRITE_PDF <<- FALSE
      if (key %in% c("pdf", "write-pdf", "write_pdf")) CONFIG$WRITE_PDF <<- TRUE
      next
    }
    key_norm <- toupper(gsub("-", "_", key))
    if (key_norm %in% names(CONFIG)) {
      if (is.logical(CONFIG[[key_norm]])) {
        CONFIG[[key_norm]] <<- tolower(val) %in% c("1", "true", "t", "yes", "y")
      } else if (is.numeric(CONFIG[[key_norm]])) {
        CONFIG[[key_norm]] <<- as.numeric(val)
      } else {
        CONFIG[[key_norm]] <<- val
      }
    } else {
      if (key %in% c("counts", "counts_csv")) CONFIG$COUNTS_CSV <<- val
      if (key %in% c("vst", "vst_csv", "vst_matrix")) CONFIG$VST_CSV <<- val
      if (key %in% c("sample_sheet", "sample_sheet_csv")) CONFIG$SAMPLE_SHEET_CSV <<- val
      if (key %in% c("outdir", "out_dir")) CONFIG$OUT_DIR <<- val
      if (key %in% c("gene_col", "gene")) CONFIG$GENE_COL <<- val
      if (key %in% c("gene_id_type")) CONFIG$GENE_ID_TYPE <<- val
      if (key %in% c("nperm", "fgsea_nperm")) CONFIG$FGSEA_NPERM <<- as.numeric(val)
      if (key %in% c("head_frac")) CONFIG$HEAD_FRAC <<- as.numeric(val)
      if (key %in% c("k_modules")) CONFIG$K_MODULES <<- as.numeric(val)
    }
  }
  invisible(NULL)
}

## ----------------------------
## Gene-set definitions
## ----------------------------
TARGET_GSETS <- c(
  "HALLMARK_ADIPOGENESIS","HALLMARK_ANDROGEN_RESPONSE","HALLMARK_ANGIOGENESIS",
  "HALLMARK_APICAL_JUNCTION","HALLMARK_APICAL_SURFACE","HALLMARK_APOPTOSIS",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS","HALLMARK_COAGULATION","HALLMARK_COMPLEMENT",
  "HALLMARK_DNA_REPAIR","HALLMARK_E2F_TARGETS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY","HALLMARK_ESTROGEN_RESPONSE_LATE",
  "HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_G2M_CHECKPOINT","HALLMARK_GLYCOLYSIS",
  "HALLMARK_HEDGEHOG_SIGNALING","HALLMARK_HEME_METABOLISM","HALLMARK_HYPOXIA",
  "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_KRAS_SIGNALING_DN","HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_MITOTIC_SPINDLE","HALLMARK_MTORC1_SIGNALING","HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_NOTCH_SIGNALING","HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_P53_PATHWAY",
  "HALLMARK_PANCREAS_BETA_CELLS","HALLMARK_PEROXISOME","HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_PROTEIN_SECRETION","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_TGF_BETA_SIGNALING","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE","HALLMARK_UV_RESPONSE_DN","HALLMARK_UV_RESPONSE_UP",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING","HALLMARK_XENOBIOTIC_METABOLISM",
  "PID_MYC_ACTIV_PATHWAY","GOMF_HEAT_SHOCK_PROTEIN_BINDING","GOMF_RNA_ENDONUCLEASE_ACTIVITY",
  "REACTOME_CHROMATIN_MODIFYING_ENZYMES","REACTOME_MRNA_SPLICING","REACTOME_GLUCOSE_METABOLISM",
  "GOBP_MITOCHONDRIAL_GENE_EXPRESSION","GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY",
  "WP_ADIPOGENESIS","WP_CELL_CYCLE","WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA",
  "WP_AMINO_ACID_METABOLISM","WP_APOPTOSIS_MODULATION_AND_SIGNALING","WP_APOPTOSIS",
  "WP_DNA_REPLICATION","WP_PROTEASOME_DEGRADATION","WP_RAS_SIGNALING","WP_PYRIMIDINE_METABOLISM",
  "GOBP_MITOCHONDRIAL_TRANSLATION",
  "REACTOME_MITOCHONDRIAL_TRANSLATION",
  "REACTOME_SIGNALING_BY_ATM",
  "REACTOME_SIGNALING_BY_ATR",
  "REACTOME_CELL_CYCLE_CHECKPOINTS",
  "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS"
)

FIG8A_MYC_SETS <- c(
  "HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2","HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT","HALLMARK_GLYCOLYSIS"
)

FIG8B_DDR_SETS <- c(
  "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS",
  "REACTOME_SIGNALING_BY_ATM",
  "REACTOME_SIGNALING_BY_ATR",
  "REACTOME_CELL_CYCLE_CHECKPOINTS",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "WP_DNA_REPLICATION"
)

FIG8D_MITO_SETS <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_MITOCHONDRIAL_TRANSLATION",
  "GOBP_MITOCHONDRIAL_TRANSLATION",
  "GOBP_MITOCHONDRIAL_GENE_EXPRESSION",
  "GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY",
  "WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA"
)

## Main-text Fig. 7C is a readability-limited union heatmap.  The full
## source table is still written, but the rendered panel is restricted to
## a balanced set of growth/metabolic/checkpoint/stress/mitochondrial
## signatures so it does not become an unreadable ontology wall.
FIG7C_UNION_SETS <- c(
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS",
  "WP_DNA_REPLICATION"
)


## Publication Fig. 7B uses a compact module-by-pathway matrix, not the
## full ontology wall. These selected terms preserve the manuscript's
## growth/biosynthesis/mitochondrial/checkpoint/stress axis while keeping
## labels readable in the main-text figure.
FIG7B_PATHWAY_SETS <- c(
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_MITOCHONDRIAL_TRANSLATION",
  "HALLMARK_DNA_REPAIR",
  "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
)

FIG8C_ES_SETS <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_DNA_REPAIR")

FIG8_CONTRASTS <- tibble::tribble(
  ~Contrast,      ~Numerator,  ~Denominator,
  "D_vs_Ctrl",    "DOX",       "CTRL",
  "DT_vs_D",      "DOX+TAM",   "DOX",
  "Tam_vs_Ctrl",  "TAM",       "CTRL"
)

## ----------------------------
## Sample parsing / metadata
## ----------------------------
derive_sample_meta <- function(samples) {
  tibble(Sample = samples) |>
    mutate(
      tok0 = normalize_token(Sample),
      tok = gsub("_mean$", "", tok0),
      tok = gsub("_+$", "", tok),
      Time = case_when(
        grepl("^24", tok) ~ "24 hours",
        grepl("^4", tok) ~ "4 hours",
        grepl("(^|_)24($|_)", tok) ~ "24 hours",
        grepl("(^|_)4($|_)", tok) ~ "4 hours",
        TRUE ~ NA_character_
      ),
      rest = case_when(
        Time == "24 hours" ~ sub("^24", "", tok),
        Time == "4 hours" ~ sub("^4", "", tok),
        TRUE ~ tok
      ),
      rest = sub("^_+", "", rest),
      Stress = case_when(
        grepl("(^|_)h(_)?cpt|hcpt|cpt_high|high_cpt|highcpt|cpt_h", rest) ~ "CPT_high",
        grepl("(^|_)l(_)?cpt|lcpt|cpt_low|low_cpt|lowcpt|cpt_l", rest) ~ "CPT_low",
        grepl("cpt", rest) ~ "CPT",
        TRUE ~ "NoCPT"
      ),
      has_ctrl = grepl("(^|_)ctrl($|_)|(^|_)control($|_)", rest),
      has_dt = grepl("(^|_)dt($|_)", rest),
      has_dox = grepl("(^|_)dox($|_)|(^|_)doxy($|_)|(^|_)d($|_)", rest),
      has_tam = grepl("(^|_)tam($|_)|(^|_)tamox($|_)", rest),
      Induction = case_when(
        has_dt ~ "DOX+TAM",
        has_dox & has_tam ~ "DOX+TAM",
        has_dox ~ "DOX",
        has_tam ~ "TAM",
        has_ctrl ~ "CTRL",
        TRUE ~ "OTHER"
      ),
      Condition = if_else(Stress == "NoCPT", Induction, paste0(Induction, "+", Stress))
    ) |>
    select(Sample, Time, Induction, Stress, Condition)
}

read_delim_auto <- function(path) {
  if (grepl("\\.tsv(\\.gz)?$|\\.txt(\\.gz)?$", path, ignore.case = TRUE)) {
    suppressMessages(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE))
  } else {
    suppressMessages(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  }
}

normalize_sample_sheet_schema <- function(ss) {
  nms <- names(ss)
  legacy_req <- c("sample_id", "time_h", "dox", "myc", "cpt_level")
  if (!all(c("Sample", "Time", "Induction", "Stress") %in% nms) && all(legacy_req %in% nms)) {
    message("[Sample sheet] Detected compatibility schema; deriving Sample/Time/Induction/Stress from sample_id/time_h/dox/myc/cpt_level.")
    ss <- ss |>
      mutate(
        Sample = as.character(sample_id),
        Time = paste0(as.character(time_h), " hours"),
        Induction = case_when(
          dox != "None_Dox" & myc == "MYC" ~ "DOX+TAM",
          dox != "None_Dox" & myc != "MYC" ~ "DOX",
          dox == "None_Dox" & myc == "MYC" ~ "TAM",
          dox == "None_Dox" & myc != "MYC" ~ "CTRL",
          TRUE ~ NA_character_
        ),
        Stress = case_when(
          cpt_level %in% c("None_CPT", "none", "None", "NoCPT", "") ~ "NoCPT",
          cpt_level %in% c("LCPT", "low", "CPT_low", "L_CPT") ~ "CPT_low",
          cpt_level %in% c("HCPT", "high", "CPT_high", "H_CPT") ~ "CPT_high",
          grepl("lcpt|low", tolower(cpt_level)) ~ "CPT_low",
          grepl("hcpt|high", tolower(cpt_level)) ~ "CPT_high",
          grepl("cpt", tolower(cpt_level)) ~ "CPT",
          TRUE ~ as.character(cpt_level)
        )
      )
  }
  ss
}

read_sample_sheet <- function(path) {
  if (!file.exists(path)) stopf("Sample sheet not found: %s", path)
  ss <- read_delim_auto(path)
  ss <- normalize_sample_sheet_schema(ss)
  req <- c("Sample", "Time", "Induction", "Stress")
  miss <- setdiff(req, names(ss))
  if (length(miss)) stopf("sample_sheet.csv is missing required canonical columns: %s. Accepted compatibility columns are sample_id,time_h,dox,myc,cpt_level.", paste(miss, collapse = ", "))
  ss <- ss |>
    mutate(
      Sample = as.character(Sample),
      Time = canonicalize_time(Time),
      Induction = canonicalize_induction(Induction),
      Stress = canonicalize_stress(Stress)
    )
  if (!"Condition" %in% names(ss)) {
    ss <- ss |> mutate(Condition = if_else(Stress == "NoCPT", Induction, paste0(Induction, "+", Stress)))
  } else {
    ss <- ss |> mutate(Condition = as.character(Condition))
  }
  if (anyDuplicated(ss$Sample)) {
    dup <- unique(ss$Sample[duplicated(ss$Sample)])
    stopf("Sample sheet contains duplicated Sample values: %s", paste(dup, collapse = ", "))
  }
  if (any(is.na(ss$Time) | ss$Time == "")) stop("sample_sheet.csv has missing Time values.", call. = FALSE)
  if (any(is.na(ss$Induction) | ss$Induction == "")) stop("sample_sheet.csv has missing Induction values.", call. = FALSE)
  if (any(is.na(ss$Stress) | ss$Stress == "")) stop("sample_sheet.csv has missing Stress values.", call. = FALSE)
  allowed_induction <- c("CTRL", "DOX", "TAM", "DOX+TAM")
  bad_induction <- setdiff(unique(ss$Induction), allowed_induction)
  if (length(bad_induction)) {
    stopf("Unsupported Induction values in sample_sheet.csv: %s", paste(bad_induction, collapse = ", "))
  }
  allowed_stress <- c("NoCPT", "CPT", "CPT_low", "CPT_high")
  bad_stress <- setdiff(unique(ss$Stress), allowed_stress)
  if (length(bad_stress)) {
    stopf("Unsupported Stress values in sample_sheet.csv: %s", paste(bad_stress, collapse = ", "))
  }
  rownames(ss) <- ss$Sample
  as.data.frame(ss, stringsAsFactors = FALSE)
}

write_sample_parse_qc <- function(sample_sheet, out_dir) {
  parsed <- derive_sample_meta(sample_sheet$Sample)
  ss <- tibble::as_tibble(sample_sheet) |>
    select(Sample, Time, Induction, Stress, Condition)
  qc <- ss |>
    left_join(parsed, by = "Sample", suffix = c("_sheet", "_parsed")) |>
    mutate(
      Time_match = Time_sheet == Time_parsed,
      Induction_match = Induction_sheet == Induction_parsed,
      Stress_match = Stress_sheet == Stress_parsed,
      Condition_match = Condition_sheet == Condition_parsed
    )
  readr::write_csv(qc, file.path(out_dir, "sample_meta_qc_compare.csv"))
  mism <- qc |>
    filter(!(Time_match & Induction_match & Stress_match & Condition_match))
  readr::write_csv(mism, file.path(out_dir, "sample_meta_qc_mismatches.csv"))
  invisible(qc)
}

## ----------------------------
## Input / VST / gene ID handling
## ----------------------------
infer_gene_id_type <- function(gene_ids, config_value = "auto") {
  if (!identical(config_value, "auto")) return(tolower(config_value))
  gene_ids <- as.character(gene_ids)
  gene_ids0 <- sub("\\..*$", "", gene_ids)
  frac_ensg_symbol <- mean(grepl("^ENSG[0-9]+(\\.[0-9]+)?_.+", gene_ids))
  frac_ensg <- mean(grepl("^ENSG[0-9]+$", gene_ids0))
  if (is.finite(frac_ensg_symbol) && frac_ensg_symbol > 0.5) {
    "ensembl_symbol"
  } else if (is.finite(frac_ensg) && frac_ensg > 0.5) {
    "ensembl"
  } else {
    "symbol"
  }
}

read_counts <- function(path, sample_sheet, gene_col = "name", min_total_counts = 80,
                        out_dir = NULL, count_input_mode = "auto") {
  if (!file.exists(path)) stopf("Counts CSV not found: %s", path)
  raw <- as.data.frame(read_delim_auto(path), check.names = FALSE)
  if (!gene_col %in% names(raw)) {
    cand <- intersect(c("name", "gene", "gene_id", "GeneID", "version", "Gene"), names(raw))
    if (!length(cand)) stop("Counts file needs a gene identifier column.", call. = FALSE)
    message(sprintf("[Counts] Using fallback gene column '%s'.", cand[1]))
    gene_col <- cand[1]
  } else {
    message(sprintf("[Counts] Using gene column '%s'.", gene_col))
  }
  sample_cols <- sample_sheet$Sample
  miss <- setdiff(sample_cols, names(raw))
  if (length(miss)) {
    stopf("Counts CSV is missing sample columns listed in sample_sheet.csv: %s",
          paste(miss, collapse = ", "))
  }
  gene_ids <- as.character(raw[[gene_col]])
  df <- raw[, sample_cols, drop = FALSE]
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(x)))
  if (anyNA(df)) stop("Counts CSV contains non-numeric values in sample columns.", call. = FALSE)
  m <- as.matrix(df)
  if (any(!is.finite(m))) stop("Counts CSV contains non-finite values in sample columns.", call. = FALSE)
  if (any(m < 0)) stop("Counts CSV contains negative values; raw counts are required.", call. = FALSE)
  frac_nonint <- mean(abs(m - round(m)) > 1e-6)
  lib_sums <- colSums(m)
  looks_tpm <- all(is.finite(lib_sums)) && stats::median(abs(lib_sums - 1e6)) < 5e3 && frac_nonint > 0.01
  filename_looks_tpm <- grepl("TPM|tpm|cpm|CPM", basename(path))
  if (isTRUE(looks_tpm) || isTRUE(filename_looks_tpm)) {
    stopf("Input matrix appears to be TPM/CPM-normalized or is named as a TPM/CPM file (path=%s; median column sum=%.2f; fraction non-integer=%.4f). Use data/processed/raw_counts_rsemgenes.tsv or provide a precomputed VST matrix with --vst_csv.",
          path, stats::median(lib_sums), frac_nonint)
  }
  count_input_mode <- tolower(as.character(count_input_mode))
  if (!count_input_mode %in% c("auto", "integer_counts", "rsem_expected_counts")) {
    stopf("Unsupported COUNT_INPUT_MODE='%s'. Use auto, integer_counts, or rsem_expected_counts.", count_input_mode)
  }
  resolved_count_mode <- count_input_mode
  rounded_for_deseq2 <- FALSE
  if (count_input_mode == "auto") {
    resolved_count_mode <- if (is.finite(frac_nonint) && frac_nonint > 1e-4) "rsem_expected_counts" else "integer_counts"
  }
  if (resolved_count_mode == "integer_counts" && is.finite(frac_nonint) && frac_nonint > 1e-4) {
    stopf("Counts matrix contains non-integer values (fraction non-integer = %.4f) but COUNT_INPUT_MODE=integer_counts. If this is the validated RSEM expected-count matrix, rerun with --count_input_mode=rsem_expected_counts or leave COUNT_INPUT_MODE=auto.",
          frac_nonint)
  }
  if (resolved_count_mode == "rsem_expected_counts") {
    rounded_for_deseq2 <- is.finite(frac_nonint) && frac_nonint > 1e-4
    message(sprintf("[Counts] Treating input as RSEM expected counts (fraction non-integer = %.4f). Values will be rounded explicitly for DESeq2 VST after the total-count filter.", frac_nonint))
  }
  if (any(grepl("_mean$", colnames(m), ignore.case = TRUE))) {
    stop("Input sample columns include *_mean names; summarized condition means are not valid VST input.", call. = FALSE)
  }
  input_qc <- tibble(
    path = path,
    resolved_count_mode = resolved_count_mode,
    rounded_for_deseq2 = rounded_for_deseq2,
    fraction_noninteger = frac_nonint,
    column_sum_min = min(lib_sums),
    column_sum_median = stats::median(lib_sums),
    column_sum_max = max(lib_sums),
    looks_tpm_or_cpm_by_column_sum = looks_tpm,
    filename_looks_tpm_or_cpm = filename_looks_tpm
  )
  if (!is.null(out_dir)) {
    readr::write_csv(input_qc, file.path(out_dir, "count_input_handling.csv"))
  }
  total_counts <- rowSums(m)
  keep <- total_counts >= min_total_counts
  if (!any(keep)) {
    stopf("No genes remained after total-count filter >= %s.", min_total_counts)
  }
  filter_summary <- tibble(
    filter = sprintf("total_counts >= %s", min_total_counts),
    genes_before = nrow(m),
    genes_after = sum(keep),
    genes_removed = sum(!keep)
  )
  if (!is.null(out_dir)) {
    readr::write_csv(filter_summary, file.path(out_dir, "gene_count_filter_summary.csv"))
  }
  m_expected <- m[keep, , drop = FALSE]
  gene_ids <- gene_ids[keep]
  rownames(m_expected) <- make.unique(gene_ids)

  m_for_vst <- m_expected
  if (resolved_count_mode == "rsem_expected_counts") {
    m_for_vst <- round(m_for_vst)
  }
  storage.mode(m_for_vst) <- "integer"
  rownames(m_for_vst) <- rownames(m_expected)
  message(sprintf("[Counts] %d genes x %d samples after total-count filter >= %s", nrow(m_for_vst), ncol(m_for_vst), min_total_counts))
  list(counts = m_for_vst, counts_expected = m_expected, gene_ids = gene_ids, gene_col = gene_col,
       filter_summary = filter_summary, input_qc = input_qc,
       resolved_count_mode = resolved_count_mode,
       rounded_for_deseq2 = rounded_for_deseq2)
}

vst_from_counts <- function(counts, meta, counts_for_normalized_baseline = NULL, out_dir = NULL) {
  cd <- as.data.frame(meta)
  rownames(cd) <- cd$Sample
  if (!identical(colnames(counts), rownames(cd))) {
    stop("Counts columns and sample-sheet rows are not aligned.", call. = FALSE)
  }
  if (!is.null(counts_for_normalized_baseline) &&
      !identical(colnames(counts_for_normalized_baseline), rownames(cd))) {
    stop("Expected-count columns and sample-sheet rows are not aligned.", call. = FALSE)
  }

  ## Estimate size factors explicitly before VST. This keeps program-level
  ## baseline ranks on the same size-factor-normalized count scale used for
  ## the distributional analysis, while VST remains the scale for heatmaps,
  ## PCA, and fgsea mean-difference ranks.
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = cd, design = ~ 1)
  dds <- DESeq2::estimateSizeFactors(dds)
  sf <- DESeq2::sizeFactors(dds)
  if (is.null(sf) || any(!is.finite(sf)) || any(sf <= 0)) {
    stop("DESeq2 size-factor estimation failed or returned non-positive values.", call. = FALSE)
  }
  if (!is.null(out_dir)) {
    readr::write_csv(
      tibble(Sample = names(sf), size_factor = as.numeric(sf)),
      file.path(out_dir, "deseq2_size_factors.csv")
    )
  }

  vs <- DESeq2::vst(dds, blind = TRUE)
  if (is.null(counts_for_normalized_baseline)) {
    counts_for_normalized_baseline <- counts
  }
  counts_for_normalized_baseline <- counts_for_normalized_baseline[, rownames(cd), drop = FALSE]
  norm_counts <- sweep(as.matrix(counts_for_normalized_baseline), 2, sf[colnames(counts_for_normalized_baseline)], "/")
  list(
    mat = SummarizedExperiment::assay(vs),
    coldata = cd,
    size_factors = sf,
    normalized_counts = norm_counts
  )
}

read_precomputed_vst <- function(path, meta, gene_col = "gene_id") {
  if (!file.exists(path)) stopf("Precomputed VST matrix not found: %s", path)
  raw <- as.data.frame(read_delim_auto(path), check.names = FALSE)
  if (!gene_col %in% names(raw)) {
    cand <- intersect(c("gene_id", "name", "gene", "GeneID", "version", "Gene"), names(raw))
    if (!length(cand)) stop("VST matrix needs a gene identifier column.", call. = FALSE)
    message(sprintf("[VST] Using fallback gene column '%s'.", cand[1]))
    gene_col <- cand[1]
  }
  sample_cols <- meta$Sample
  miss <- setdiff(sample_cols, names(raw))
  if (length(miss)) {
    stopf("VST matrix is missing sample columns listed in sample_sheet.csv: %s", paste(miss, collapse = ", "))
  }
  gene_ids <- as.character(raw[[gene_col]])
  df <- raw[, sample_cols, drop = FALSE]
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(x)))
  if (anyNA(df)) stop("VST matrix contains non-numeric values in sample columns.", call. = FALSE)
  m <- as.matrix(df)
  if (any(!is.finite(m))) stop("VST matrix contains non-finite values in sample columns.", call. = FALSE)
  rownames(m) <- make.unique(gene_ids)
  cd <- as.data.frame(meta)
  rownames(cd) <- cd$Sample
  list(mat = m, coldata = cd, gene_ids = gene_ids, gene_col = gene_col)
}

gene_ids_to_symbols <- function(gene_ids, gene_id_type = "auto") {
  id_type <- infer_gene_id_type(gene_ids, gene_id_type)
  ids0 <- as.character(gene_ids)
  if (id_type == "ensembl_symbol") {
    symbols <- sub("^ENSG[0-9]+(\\.[0-9]+)?_", "", ids0)
  } else if (id_type == "ensembl") {
    ens <- sub("\\..*$", "", ids0)
    sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
    symbols <- unname(sym)
  } else {
    symbols <- ids0
  }
  as.character(symbols)
}

collapse_matrix_to_symbols <- function(mat, gene_ids, gene_id_type = "auto", tie_breaker = NULL) {
  if (length(gene_ids) != nrow(mat)) stop("gene_ids length does not match matrix rows.", call. = FALSE)
  symbols <- gene_ids_to_symbols(gene_ids, gene_id_type)
  keep <- !is.na(symbols) & symbols != ""
  if (!any(keep)) stop("No genes remained after mapping IDs to symbols.", call. = FALSE)
  mat_keep <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (is.null(tie_breaker)) tie_breaker <- rowMeans(mat_keep, na.rm = TRUE) else tie_breaker <- as.numeric(tie_breaker[keep])
  tie_breaker[!is.finite(tie_breaker)] <- -Inf
  idx_tbl <- tibble(row_idx = seq_len(nrow(mat_keep)), symbol = symbols, tie_breaker = tie_breaker) |>
    arrange(symbol, desc(tie_breaker)) |>
    distinct(symbol, .keep_all = TRUE)
  mat2 <- mat_keep[idx_tbl$row_idx, , drop = FALSE]
  rownames(mat2) <- idx_tbl$symbol
  mat2
}

collapse_vst_to_symbols <- function(vst_mat, gene_ids, gene_id_type = "auto") {
  collapse_matrix_to_symbols(vst_mat, gene_ids, gene_id_type = gene_id_type, tie_breaker = rowMeans(vst_mat, na.rm = TRUE))
}

## ----------------------------
## Core statistics
## ----------------------------
safe_pr_t <- function(Xt, ncomp) {
  n <- nrow(Xt)
  p <- ncol(Xt)
  ncomp <- max(1, min(ncomp, min(n, p) - 1))
  pr <- NULL
  if (ncomp / min(n, p) <= 0.3) {
    pr <- tryCatch(
      irlba::prcomp_irlba(Xt, n = ncomp, center = TRUE, scale. = FALSE),
      error = function(e) NULL
    )
  }
  if (is.null(pr)) {
    pr_full <- stats::prcomp(Xt, center = TRUE, scale. = FALSE)
    pr <- list(
      rotation = pr_full$rotation[, seq_len(ncomp), drop = FALSE],
      x = pr_full$x[, seq_len(ncomp), drop = FALSE]
    )
  }
  pr
}

time_aligned_axis <- function(vst_mat, coldata, time_var = "Time", ncomp = 10) {
  valid <- !is.na(coldata[[time_var]])
  if (sum(valid) < 4) stop("Too few samples with non-missing Time for distortion-associated PC.", call. = FALSE)
  X <- vst_mat[, valid, drop = FALSE]
  cd <- coldata[valid, , drop = FALSE]
  tt <- factor(cd[[time_var]], levels = resolve_time_levels(cd[[time_var]]))
  if (length(levels(tt)) < 2) stop("Need at least two time levels to compute the distortion-associated PC.", call. = FALSE)
  y <- as.numeric(tt) - 1
  Xc <- X - rowMeans(X)
  pr <- safe_pr_t(t(Xc), ncomp = ncomp)
  cors <- vapply(seq_len(ncol(pr$x)), function(k) {
    suppressWarnings(cor(pr$x[, k], y, method = "spearman", use = "complete.obs"))
  }, numeric(1))
  k <- which.max(abs(cors))
  sgn <- ifelse(is.na(cors[k]) || cors[k] >= 0, 1, -1)
  tibble(
    gene = rownames(X),
    loading = sgn * as.numeric(pr$rotation[, k]),
    comp = k,
    comp_time_cor = sgn * cors[k],
    axis_mode = "distortion_associated_time_correlated_pc_oriented",
    axis_polarity = "positive_loading=later_more_distorted"
  )
}

mean_by_time <- function(vst_mat, coldata, time_var = "Time") {
  valid <- !is.na(coldata[[time_var]])
  cd <- coldata[valid, , drop = FALSE]
  if (nrow(cd) < 2) stop("Too few samples with valid Time labels.", call. = FALSE)
  tt <- factor(cd[[time_var]], levels = resolve_time_levels(cd[[time_var]]))
  lv <- levels(tt)
  if (length(lv) < 2) stop("Need at least two time levels.", call. = FALSE)
  X <- vst_mat[, rownames(cd), drop = FALSE]
  m1 <- rowMeans(X[, tt == lv[1], drop = FALSE], na.rm = TRUE)
  m2 <- rowMeans(X[, tt == lv[2], drop = FALSE], na.rm = TRUE)
  list(levels = lv[1:2], mean1 = m1, mean2 = m2)
}

nri_table <- function(vst_mat, coldata, time_var = "Time", k = NULL) {
  ms <- mean_by_time(vst_mat, coldata, time_var)
  r1 <- baseline_ranks01(ms$mean1)
  r2 <- baseline_ranks01(ms$mean2)
  z1 <- (ms$mean1 - median(ms$mean1)) / (mad(ms$mean1) + 1e-8)
  z2 <- (ms$mean2 - median(ms$mean2)) / (mad(ms$mean2) + 1e-8)
  F1 <- cbind(r1, z1)
  F2 <- cbind(r2, z2)
  n <- nrow(vst_mat)
  if (is.null(k)) k <- min(20L, max(5L, floor(0.01 * n)))
  nn1 <- FNN::get.knn(F1, k = k)$nn.index
  nn2 <- FNN::get.knn(F2, k = k)$nn.index
  nri <- vapply(seq_len(n), function(i) {
    a <- nn1[i, ]
    b <- nn2[i, ]
    1 - length(intersect(a, b)) / length(unique(c(a, b)))
  }, numeric(1))
  tibble(
    gene = rownames(vst_mat),
    NRI = as.numeric(nri),
    rank_early = as.numeric(r1),
    rank_late = as.numeric(r2),
    delta_rank = as.numeric(r2 - r1)
  )
}

condition_profile_matrix <- function(vst_mat, cd) {
  cd2 <- as.data.frame(cd)
  rownames(cd2) <- rownames(cd)
  ord <- order_samples(cd2)
  cd2 <- cd2[ord, , drop = FALSE]
  cd2$Profile <- paste(short_time(cd2$Time), short_stress(cd2$Stress), cd2$Induction, sep = "__")
  prof_levels <- unique(cd2$Profile)
  mats <- lapply(prof_levels, function(p) {
    samples <- rownames(cd2)[cd2$Profile == p]
    rowMeans(vst_mat[, samples, drop = FALSE], na.rm = TRUE)
  })
  out <- do.call(cbind, mats)
  colnames(out) <- prof_levels
  rownames(out) <- rownames(vst_mat)
  out
}

compute_baseline <- function(expr_mat, cd) {
  idx <- which(cd$Time == "4 hours" & cd$Induction == "CTRL" & cd$Stress == "NoCPT")
  source_rule <- "4h_CTRL_NoCPT"
  if (!length(idx)) {
    idx <- which(cd$Time == "4 hours" & cd$Induction == "CTRL")
    source_rule <- "4h_CTRL_any_CPT"
  }
  if (!length(idx)) {
    idx <- which(cd$Time == "4 hours")
    source_rule <- "4h_any_state"
  }
  if (!length(idx)) stop("No 4 h baseline samples found.", call. = FALSE)
  baseline <- rowMeans(expr_mat[, idx, drop = FALSE], na.rm = TRUE)
  attr(baseline, "baseline_samples") <- colnames(expr_mat)[idx]
  attr(baseline, "baseline_rule") <- source_rule
  baseline
}

write_baseline_source_audit <- function(baseline_vec, source_matrix_name, out_dir) {
  readr::write_csv(
    tibble(
      baseline_source_matrix = source_matrix_name,
      baseline_rule = attr(baseline_vec, "baseline_rule"),
      baseline_samples = paste(attr(baseline_vec, "baseline_samples"), collapse = ";"),
      n_genes = length(baseline_vec)
    ),
    file.path(out_dir, "baseline_rank_source.csv")
  )
}

band_from_baseline <- function(baseline_vec) {
  n <- length(baseline_vec)
  rk <- rank(-baseline_vec, ties.method = "average")
  pct <- (rk - 0.5) / n
  out <- dplyr::case_when(
    pct <= 0.10 ~ "head",
    pct >= 0.90 ~ "tail",
    pct >= 0.25 & pct <= 0.75 ~ "mid",
    TRUE ~ NA_character_
  )
  factor(out, levels = c("head", "mid", "tail"))
}

head_genes_from_baseline <- function(baseline_vec, head_frac = 0.10) {
  n <- length(baseline_vec)
  head_n <- max(1, floor(n * head_frac))
  names(sort(baseline_vec, decreasing = TRUE))[seq_len(head_n)]
}


safe_pam <- function(X, k) {
  ## cluster::pam arguments vary a little across cluster versions.  Try the
  ## faster modern variant first and fall back to the basic call if unavailable.
  tryCatch(
    cluster::pam(X, k = k, metric = "euclidean", stand = FALSE, pamonce = 6),
    error = function(e) cluster::pam(X, k = k, metric = "euclidean", stand = FALSE)
  )
}

approx_membership_to_medoids <- function(X, medoids, clustering) {
  if (is.null(medoids)) return(rep(NA_real_, nrow(X)))
  medoids <- as.matrix(medoids)
  if (!nrow(medoids) || !ncol(medoids)) return(rep(NA_real_, nrow(X)))
  dmat <- vapply(seq_len(nrow(medoids)), function(j) {
    sqrt(rowSums((sweep(X, 2, medoids[j, ], FUN = "-"))^2))
  }, numeric(nrow(X)))
  if (is.null(dim(dmat))) dmat <- matrix(dmat, ncol = 1)
  k <- ncol(dmat)
  if (k < 2) return(rep(NA_real_, nrow(X)))
  a <- dmat[cbind(seq_len(nrow(X)), pmax(1, pmin(k, as.integer(clustering))))]
  tmp <- dmat
  tmp[cbind(seq_len(nrow(X)), pmax(1, pmin(k, as.integer(clustering))))] <- Inf
  b <- apply(tmp, 1, min, na.rm = TRUE)
  out <- 1 - a / (b + 1e-8)
  out[!is.finite(out)] <- NA_real_
  pmax(-1, pmin(1, out))
}

module_silhouette_audit <- function(X, k_range, sample_n = 4000L, seed = 1771L) {
  n <- nrow(X)
  sample_n <- as.integer(min(n, max(100L, sample_n)))
  set.seed(seed)
  idx <- if (sample_n < n) sort(sample(seq_len(n), sample_n)) else seq_len(n)
  X_sub <- X[idx, , drop = FALSE]
  message(sprintf("[Modules] Silhouette audit uses exact PAM on %s/%s genes.", nrow(X_sub), n))
  scores <- vapply(k_range, function(k) {
    if (k >= nrow(X_sub)) return(NA_real_)
    fit <- safe_pam(X_sub, k = k)
    as.numeric(fit$silinfo$avg.width)
  }, numeric(1))
  list(scores = scores, sample_n = length(idx), sample_total = n)
}

run_module_clustering <- function(X, k, engine = "auto", seed = 1771L,
                                  exact_pam_max_n = 5000L,
                                  clara_samples = 20L,
                                  clara_sampsize = 5000L) {
  n <- nrow(X)
  engine <- tolower(engine)
  if (engine == "auto") engine <- if (n <= exact_pam_max_n) "exact_pam" else "clara"

  if (engine == "exact_pam") {
    message(sprintf("[Modules] Running exact PAM for final k=%s on %s genes.", k, n))
    fit <- safe_pam(X, k = k)
    sil_width <- rep(NA_real_, n)
    if (!is.null(fit$silinfo$widths) && nrow(fit$silinfo$widths) == n) {
      sil_width <- as.numeric(fit$silinfo$widths[, "sil_width"])
    }
    return(list(
      clustering = as.integer(fit$clustering),
      membership = sil_width,
      engine_used = "exact_pam",
      final_avg_width = as.numeric(fit$silinfo$avg.width),
      clara_sampsize = NA_integer_,
      clara_samples = NA_integer_
    ))
  }

  if (engine != "clara") stopf("Unknown MODULE_CLUSTER_ENGINE: %s", engine)
  clara_sampsize <- as.integer(min(n, max(40L + 2L * k, clara_sampsize)))
  clara_samples <- as.integer(max(5L, clara_samples))
  message(sprintf(
    "[Modules] Running low-memory CLARA/PAM approximation for final k=%s on %s genes (samples=%s, sampsize=%s).",
    k, n, clara_samples, clara_sampsize
  ))
  set.seed(seed)
  fit <- tryCatch(
    cluster::clara(
      X, k = k, metric = "euclidean", stand = FALSE,
      samples = clara_samples, sampsize = clara_sampsize,
      rngR = TRUE, pamLike = TRUE, medoids.x = TRUE, keep.data = FALSE
    ),
    error = function(e) {
      message(sprintf("[Modules] CLARA call with modern arguments failed (%s); retrying with basic cluster::clara arguments.", e$message))
      cluster::clara(
        X, k = k, metric = "euclidean", stand = FALSE,
        samples = clara_samples, sampsize = clara_sampsize,
        medoids.x = TRUE, keep.data = FALSE
      )
    }
  )
  med <- tryCatch(as.matrix(fit$medoids), error = function(e) NULL)
  membership <- approx_membership_to_medoids(X, med, fit$clustering)
  list(
    clustering = as.integer(fit$clustering),
    membership = membership,
    engine_used = "clara_pam_approximation",
    final_avg_width = suppressWarnings(as.numeric(fit$silinfo$avg.width)),
    clara_sampsize = clara_sampsize,
    clara_samples = clara_samples
  )
}

build_modules <- function(vst_mat, coldata, load_tbl, nri_tbl, baseline_vec,
                          head_frac = 0.10, k_min = 3, k_max = 8,
                          final_k = 7L, require_seven = TRUE, seed = 1771L,
                          cluster_engine = "clara",
                          exact_pam_max_n = 5000L,
                          silhouette_sample_n = 4000L,
                          clara_samples = 20L,
                          clara_sampsize = 5000L) {
  ms <- mean_by_time(vst_mat, coldata)
  lfc_sur <- ms$mean2 - ms$mean1

  features <- load_tbl |>
    select(gene, loading) |>
    left_join(nri_tbl |> select(gene, NRI, rank_early, rank_late, delta_rank), by = "gene") |>
    mutate(lfc_sur = as.numeric(lfc_sur[gene])) |>
    filter(!is.na(gene))

  feature_cols <- c("loading", "NRI", "rank_early", "rank_late", "delta_rank")
  X <- as.matrix(features[, feature_cols, drop = FALSE])
  rownames(X) <- features$gene
  X[!is.finite(X)] <- 0
  X <- scale(X)
  X[!is.finite(X)] <- 0

  k_range <- seq.int(as.integer(k_min), as.integer(k_max))
  k_range <- k_range[k_range >= 2 & k_range < nrow(X)]
  if (!length(k_range)) stop("No valid module k values for PAM clustering.", call. = FALSE)

  final_k <- as.integer(final_k)
  if (is.na(final_k) || final_k < 2L) stop("final_k must be an integer >= 2.", call. = FALSE)
  if (!(final_k %in% k_range)) {
    stopf("Requested final module k=%s is outside the available PAM k-range [%s].", final_k, paste(k_range, collapse = ","))
  }
  if (isTRUE(require_seven) && final_k != 7L) {
    stopf("Final module k=%s, but frozen manuscript expects seven modules. Use --k_modules=7 or --require_seven_modules=false.", final_k)
  }

  message(sprintf("[Modules] Feature matrix: %s genes x %s module features.", nrow(X), ncol(X)))
  audit <- module_silhouette_audit(X, k_range, sample_n = silhouette_sample_n, seed = seed)
  scores <- audit$scores
  silhouette_best_k <- as.integer(k_range[which.max(scores)])

  clustering <- run_module_clustering(
    X, k = final_k, engine = cluster_engine, seed = seed,
    exact_pam_max_n = exact_pam_max_n,
    clara_samples = clara_samples,
    clara_sampsize = clara_sampsize
  )
  cl <- as.integer(unlist(clustering$clustering, use.names = FALSE))
  sil_width <- as.numeric(unlist(clustering$membership, use.names = FALSE))
  if (length(cl) != nrow(features)) {
    stopf("Module clustering returned %s labels for %s genes.", length(cl), nrow(features))
  }
  if (length(sil_width) != nrow(features)) {
    sil_width <- rep(NA_real_, nrow(features))
  }
  chosen_k <- final_k

  if (silhouette_best_k != chosen_k) {
    message(sprintf(
      "[Modules] PAM silhouette audit optimum is k=%s, but using manuscript/display resolution k=%s. See module_k_selection.csv.",
      silhouette_best_k, chosen_k
    ))
  }

  bands <- as.character(band_from_baseline(baseline_vec))
  names(bands) <- names(baseline_vec)

  assign_raw <- tibble(
    gene = features$gene,
    cluster_raw = cl,
    module_membership = sil_width,
    lfc_sur = as.numeric(features$lfc_sur),
    pc_loading = as.numeric(features$loading),
    NRI = as.numeric(features$NRI),
    rank_early = as.numeric(features$rank_early),
    rank_late = as.numeric(features$rank_late),
    delta_rank = as.numeric(features$delta_rank)
  )

  mdi_raw <- assign_raw |>
    mutate(band = bands[gene]) |>
    filter(!is.na(band)) |>
    group_by(cluster_raw, band) |>
    summarise(med = median(lfc_sur, na.rm = TRUE), .groups = "drop") |>
    tidyr::complete(cluster_raw, band = c("head", "mid", "tail"), fill = list(med = NA_real_)) |>
    pivot_wider(names_from = band, values_from = med) |>
    mutate(MDI = head - tail)

  module_map <- mdi_raw |>
    arrange(desc(MDI), cluster_raw) |>
    mutate(cluster = paste0("M", row_number())) |>
    select(cluster_raw, cluster)

  assignments <- assign_raw |>
    left_join(module_map, by = "cluster_raw") |>
    select(gene, cluster, module_membership, lfc_sur, pc_loading, NRI, rank_early, rank_late, delta_rank)

  mdi <- mdi_raw |>
    left_join(module_map, by = "cluster_raw") |>
    select(cluster, head, mid, tail, MDI) |>
    arrange(match(cluster, paste0("M", seq_len(nrow(module_map)))))

  k_selection <- tibble(k = k_range, silhouette_avg_width = scores) |>
    mutate(
      silhouette_best = k == silhouette_best_k,
      used_for_figure = k == chosen_k,
      selected = used_for_figure,
      module_k_rule = "fixed_7_with_silhouette_audit_low_memory",
      silhouette_audit_sample_n = audit$sample_n,
      silhouette_audit_total_genes = audit$sample_total,
      final_clustering_engine = clustering$engine_used,
      final_clara_samples = clustering$clara_samples,
      final_clara_sampsize = clustering$clara_sampsize,
      final_avg_width = clustering$final_avg_width,
      note = if_else(
        k == chosen_k & silhouette_best_k != chosen_k,
        paste0("Used for frozen manuscript figure; silhouette audit optimum was k=", silhouette_best_k),
        ""
      )
    )

  list(
    assignments = assignments,
    mdi = mdi,
    k_selection = k_selection,
    selected_k = chosen_k,
    silhouette_best_k = silhouette_best_k,
    clustering_engine = clustering$engine_used
  )
}

## ----------------------------
## Sample annotations / heatmaps
## ----------------------------
sample_annotation <- function(cd) {
  ## Short, publication-facing annotation labels keep legends from being
  ## clipped in Fig. 7A/C and make the heatmaps self-decoding.
  ann <- as.data.frame(cd[, c("Time", "Induction", "Stress"), drop = FALSE])
  ann <- ann |>
    mutate(
      Time = dplyr::case_when(
        Time == "4 hours" ~ "4 h",
        Time == "24 hours" ~ "24 h",
        TRUE ~ as.character(Time)
      ),
      State = dplyr::case_when(
        Induction == "CTRL" ~ "Ctrl",
        Induction == "DOX" ~ "D",
        Induction == "TAM" ~ "Tam",
        Induction == "DOX+TAM" ~ "DT",
        TRUE ~ as.character(Induction)
      ),
      CPT = dplyr::case_when(
        Stress == "NoCPT" ~ "No CPT",
        Stress == "CPT_low" ~ "Low CPT",
        Stress == "CPT_high" ~ "High CPT",
        Stress == "CPT" ~ "CPT",
        TRUE ~ as.character(Stress)
      )
    ) |>
    select(Time, State, CPT)
  ann$Time <- factor(ann$Time, levels = c("4 h", "24 h"))
  ann$State <- factor(ann$State, levels = c("Ctrl", "D", "Tam", "DT"))
  ann$CPT <- factor(ann$CPT, levels = c("No CPT", "Low CPT", "High CPT", "CPT"))
  rownames(ann) <- rownames(cd)
  ann
}

pheatmap_annotation_colors <- function(...) {
  ann_list <- list(...)
  ann_list <- ann_list[vapply(ann_list, function(x) is.data.frame(x) && ncol(x) > 0, logical(1))]
  if (!length(ann_list)) return(NULL)
  col_names <- unique(unlist(lapply(ann_list, names), use.names = FALSE))
  out <- list()
  for (nm in col_names) {
    vals <- unique(unlist(lapply(ann_list, function(df) {
      if (nm %in% names(df)) as.character(df[[nm]]) else character()
    }), use.names = FALSE))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (!length(vals)) next
    if (nm == "Time") vals <- intersect(c("4 h", "24 h"), vals)
    if (nm == "State") vals <- intersect(c("Ctrl", "D", "Tam", "DT"), vals)
    if (nm == "CPT") vals <- intersect(c("No CPT", "Low CPT", "High CPT", "CPT"), vals)
    if (nm == "Band") vals <- intersect(c("head", "mid", "tail"), vals)
    if (nm == "Module") vals <- intersect(paste0("M", seq_len(30)), vals)
    preset <- switch(
      nm,
      Time = c("4 h" = "#80B1D3", "24 h" = "#FDB462"),
      State = c("Ctrl" = "#66C2A5", "D" = "#FC8D62", "Tam" = "#8DA0CB", "DT" = "#E78AC3"),
      CPT = c("No CPT" = "#BDBDBD", "Low CPT" = "#FDBF6F", "High CPT" = "#E31A1C", "CPT" = "#FDBF6F"),
      Band = c("head" = "#D01C8B", "mid" = "#F1B6DA", "tail" = "#4DAC26"),
      NULL
    )
    if (nm == "Module") preset <- setNames(grDevices::hcl.colors(max(3, length(vals)), palette = "Dark 3")[seq_along(vals)], vals)
    if (nm == "TopSet") preset <- setNames(grDevices::hcl.colors(max(3, length(vals)), palette = "Set 3")[seq_along(vals)], vals)
    if (is.null(preset)) {
      out[[nm]] <- setNames(grDevices::hcl.colors(max(3, length(vals)), palette = "Dark 3")[seq_along(vals)], vals)
    } else {
      kept <- preset[names(preset) %in% vals]
      missing <- setdiff(vals, names(kept))
      if (length(missing)) {
        kept <- c(kept, setNames(grDevices::hcl.colors(max(3, length(missing)), palette = "Dark 3")[seq_along(missing)], missing))
      }
      out[[nm]] <- kept[vals]
    }
  }
  out
}

write_annotation_key <- function(annotation_colors, out_csv) {
  if (is.null(annotation_colors) || !length(annotation_colors)) return(invisible(NULL))
  key <- purrr::imap_dfr(annotation_colors, function(cols, nm) {
    tibble::tibble(annotation = nm, level = names(cols), color = unname(cols))
  })
  readr::write_csv(key, out_csv)
  invisible(key)
}

order_samples <- function(cd) {
  tt <- factor(cd$Time, levels = resolve_time_levels(cd$Time))
  ss <- factor(cd$Stress, levels = c("NoCPT", "CPT", "CPT_low", "CPT_high"))
  ii <- factor(cd$Induction, levels = c("CTRL", "DOX", "TAM", "DOX+TAM", "OTHER"))
  order(tt, ss, ii, rownames(cd))
}

sample_correlation_heatmap <- function(vst_mat, cd, out_file) {
  cm <- stats::cor(vst_mat, method = "spearman", use = "pairwise.complete.obs")
  ann <- sample_annotation(cd)
  ph <- pheatmap::pheatmap(
    cm,
    annotation_col = ann,
    clustering_method = "average",
    main = "Sample correlation (VST)",
    silent = TRUE
  )
  g <- if (is.list(ph) && !is.null(ph$gtable)) ph$gtable else ph
  save_gtable_png_svg(g, out_file, width = 8.5, height = 7.0, res = 300)
}

## ----------------------------
## MSigDB helpers / enrichment
## ----------------------------
msigdbr_get <- function(collection, subcollection = NULL) {
  ## msigdbr changed argument names across releases.  Try the current
  ## collection/subcollection API first and fall back to category/subcategory
  ## for older environments.
  if (is.null(subcollection)) {
    out <- try(msigdbr::msigdbr(species = "Homo sapiens", collection = collection), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    out <- try(msigdbr::msigdbr(species = "Homo sapiens", category = collection), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
  } else {
    out <- try(msigdbr::msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    out <- try(msigdbr::msigdbr(species = "Homo sapiens", category = collection, subcategory = subcollection), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
  }
  stopf("Unable to query msigdbr collection=%s subcollection=%s. Check the installed msigdbr version.",
        collection, ifelse(is.null(subcollection), "<none>", subcollection))
}

build_msig_term2gene <- function(target_sets = NULL) {
  ## Robust across msigdbr releases and partial local installations.
  ## We try several collection/subcollection spellings and keep whichever
  ## calls succeed. Missing collections are reported via MSIGDB_unmatched_sets.txt.
  calls <- list(
    c("H", NA),
    c("C2", "CP:REACTOME"),
    c("C2", "CP:PID"),
    c("C2", "CP:WIKIPATHWAYS"),
    c("C5", "GO:BP"),
    c("C5", "BP"),
    c("C5", "GO:MF"),
    c("C5", "MF")
  )
  parts <- list()
  for (cc in calls) {
    coll <- cc[[1]]
    subc <- cc[[2]]
    tmp <- tryCatch(
      if (is.na(subc)) msigdbr_get(coll) else msigdbr_get(coll, subc),
      error = function(e) NULL
    )
    if (!is.null(tmp) && nrow(tmp)) {
      parts[[length(parts) + 1L]] <- tmp |> select(gs_name, gene_symbol) |> distinct()
    }
  }
  if (!length(parts)) stop("No MSigDB gene sets were available from msigdbr.", call. = FALSE)
  all <- bind_rows(parts) |> distinct(gs_name, gene_symbol)
  if (!is.null(target_sets) && length(target_sets)) {
    keep <- norm_key(target_sets)
    all <- all |>
      mutate(gs_norm = norm_key(gs_name)) |>
      filter(gs_norm %in% keep) |>
      select(gs_name, gene_symbol)
  }
  all
}

report_unmatched_sets <- function(requested_sets, t2g, out_file) {
  req <- unique(norm_key(requested_sets))
  got <- unique(norm_key(unique(t2g$gs_name)))
  miss_norm <- setdiff(req, got)
  if (length(miss_norm)) {
    miss <- requested_sets[norm_key(requested_sets) %in% miss_norm] |> unique()
    writeLines(miss, con = out_file)
  } else {
    writeLines(character(0), con = out_file)
  }
  invisible(miss_norm)
}

t2g_to_genesets <- function(t2g) {
  t2g <- t2g |>
    distinct(gs_name, gene_symbol) |>
    filter(!is.na(gs_name), !is.na(gene_symbol))
  split(t2g$gene_symbol, t2g$gs_name) |>
    lapply(unique)
}

label_modules_with_msigdb <- function(assign_df, universe_genes, t2g, out_csv) {
  gene2cl <- assign_df |>
    as_tibble() |>
    transmute(
      SYMBOL = as.character(unlist(gene, use.names = FALSE)),
      Cluster = as.character(unlist(cluster, use.names = FALSE))
    ) |>
    distinct(SYMBOL, Cluster) |>
    filter(!is.na(SYMBOL), !is.na(Cluster), nzchar(SYMBOL), nzchar(Cluster))
  clv <- sort(unique(gene2cl$Cluster))
  tab <- purrr::map_dfr(clv, function(cl) {
    genes <- unique(gene2cl$SYMBOL[gene2cl$Cluster == cl])
    if (!length(genes)) {
      return(tibble(Cluster = cl, TopTerm = NA_character_, FDR = NA_real_))
    }
    enr <- try(
      clusterProfiler::enricher(
        gene = genes,
        TERM2GENE = t2g,
        universe = universe_genes,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      ),
      silent = TRUE
    )
    if (inherits(enr, "try-error") || is.null(enr) || !nrow(as.data.frame(enr))) {
      return(tibble(Cluster = cl, TopTerm = NA_character_, FDR = NA_real_))
    }
    as.data.frame(enr) |>
      arrange(p.adjust, desc(Count)) |>
      slice_head(n = 1) |>
      transmute(Cluster = cl, TopTerm = Description, FDR = p.adjust)
  })
  readr::write_csv(tab, out_csv)
  tab
}


short_pathway_label <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^REACTOME_", "", x)
  x <- gsub("^GOBP_", "", x)
  x <- gsub("^GOMF_", "", x)
  x <- gsub("^WP_", "", x)
  x <- gsub("^PID_", "", x)
  x <- gsub("_", " ", x)
  x <- stringr::str_to_title(tolower(x))
  x <- gsub("Myc", "MYC", x)
  x <- gsub("Dna", "DNA", x)
  x <- gsub("Rna", "RNA", x)
  x <- gsub("G2m", "G2M", x)
  x <- gsub("E2f", "E2F", x)
  x <- gsub("Atm", "ATM", x)
  x <- gsub("Atr", "ATR", x)
  x <- gsub("Cpt", "CPT", x)
  x <- gsub("Oxphos", "OXPHOS", x)
  stringr::str_trunc(x, width = 42)
}


pathway_source_prefix <- function(x) {
  case_when(
    grepl("^HALLMARK_", x) ~ "Hallmark: ",
    grepl("^REACTOME_", x) ~ "Reactome: ",
    grepl("^GOBP_", x) ~ "GO BP: ",
    grepl("^GOMF_", x) ~ "GO MF: ",
    grepl("^WP_", x) ~ "WikiPathways: ",
    grepl("^PID_", x) ~ "PID: ",
    TRUE ~ ""
  )
}

unique_pathway_display_labels <- function(pathways, wrap_width = 30) {
  ## short_pathway_label() deliberately removes source prefixes to improve
  ## readability.  That can create duplicated labels for mitochondrial sets,
  ## e.g. REACTOME_MITOCHONDRIAL_TRANSLATION and
  ## GOBP_MITOCHONDRIAL_TRANSLATION.  ggplot factors cannot have duplicated
  ## levels, so add source prefixes only for ambiguous labels and use
  ## make.unique() as a final safety net.
  pathways <- unique(as.character(pathways))
  base <- short_pathway_label(pathways)
  ambiguous <- base %in% base[duplicated(base)]
  labels <- base
  labels[ambiguous] <- paste0(pathway_source_prefix(pathways[ambiguous]), base[ambiguous])
  labels <- stringr::str_wrap(labels, width = wrap_width)
  if (any(duplicated(labels))) {
    labels <- make.unique(labels, sep = " ")
  }
  names(labels) <- pathways
  labels
}

figure7c_set_label <- function(x) {
  ## Compact, manuscript-facing labels for the Fig. 7C row annotation legend.
  ## The source table preserves full MSigDB names; this display label keeps the
  ## heatmap self-decoding without turning the legend into an ontology wall.
  key <- norm_key(x)
  dplyr::case_when(
    key == norm_key("HALLMARK_MYC_TARGETS_V1") ~ "MYC targets",
    key == norm_key("HALLMARK_MTORC1_SIGNALING") ~ "mTORC1",
    key == norm_key("HALLMARK_G2M_CHECKPOINT") ~ "G2M checkpoint",
    key == norm_key("HALLMARK_P53_PATHWAY") ~ "p53 pathway",
    key == norm_key("HALLMARK_GLYCOLYSIS") ~ "Glycolysis",
    key == norm_key("HALLMARK_APOPTOSIS") ~ "Apoptosis",
    key == norm_key("HALLMARK_OXIDATIVE_PHOSPHORYLATION") ~ "OXPHOS",
    key == norm_key("HALLMARK_E2F_TARGETS") ~ "E2F targets",
    key == norm_key("WP_DNA_REPLICATION") ~ "DNA replication",
    key == norm_key("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY") ~ "ROS pathway",
    key == norm_key("HALLMARK_DNA_REPAIR") ~ "DNA repair",
    key == norm_key("REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS") ~ "ATR/replication stress",
    TRUE ~ stringr::str_trunc(short_pathway_label(x), width = 28)
  )
}

modules_terms_heatmap <- function(assign_df, universe_genes, t2g, out_png,
                                  top_terms_per_module = CONFIG$FIG7B_TOP_TERMS,
                                  selected_terms = FIG7B_PATHWAY_SETS) {
  ## Publication-facing panel B: readable module-by-pathway matrix.
  ## The full enrichment table is written to CSV, but the rendered figure is
  ## restricted to a curated program axis so it matches the manuscript while
  ## remaining legible.
  gene2cl <- assign_df |>
    as_tibble() |>
    transmute(
      SYMBOL = as.character(unlist(gene, use.names = FALSE)),
      Cluster = as.character(unlist(cluster, use.names = FALSE))
    ) |>
    distinct(SYMBOL, Cluster) |>
    filter(!is.na(SYMBOL), !is.na(Cluster), nzchar(SYMBOL), nzchar(Cluster))

  clv <- sort(unique(gene2cl$Cluster))
  if (all(grepl("^M[0-9]+$", clv))) clv <- clv[order(as.integer(sub("^M", "", clv)))]

  term_tabs <- lapply(clv, function(cl) {
    genes <- unique(gene2cl$SYMBOL[gene2cl$Cluster == cl])
    enr <- try(
      clusterProfiler::enricher(
        gene = genes,
        TERM2GENE = t2g,
        universe = universe_genes,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1
      ),
      silent = TRUE
    )
    if (inherits(enr, "try-error") || is.null(enr) || !nrow(as.data.frame(enr))) {
      return(tibble(Cluster = cl, Description = character(), score = numeric(), FDR = numeric(), Count = integer()))
    }
    as_tibble(as.data.frame(enr)) |>
      mutate(
        Description = as.character(Description),
        score = -log10(p.adjust + 1e-300),
        Cluster = cl,
        FDR = p.adjust
      ) |>
      arrange(desc(score), desc(Count), Description) |>
      select(Cluster, Description, score, FDR, Count)
  })

  full_df <- bind_rows(term_tabs)
  readr::write_csv(full_df, sub("\\.png$", "_full_enrichment.csv", out_png))

  if (!nrow(full_df)) {
    warning("No module enrichment results available for Fig. 7B; writing placeholder.", call. = FALSE)
    full_df <- tibble(Cluster = clv, Description = NA_character_, score = 0, FDR = NA_real_, Count = NA_integer_)
  }

  selected_norm <- norm_key(selected_terms)
  selected_map <- t2g |>
    distinct(gs_name) |>
    mutate(gs_norm = norm_key(gs_name)) |>
    filter(gs_norm %in% selected_norm) |>
    mutate(order = match(gs_norm, selected_norm)) |>
    arrange(order) |>
    distinct(gs_norm, .keep_all = TRUE)

  if (!nrow(selected_map)) {
    selected_map <- full_df |>
      distinct(Description) |>
      filter(!is.na(Description)) |>
      slice_head(n = 10) |>
      transmute(gs_name = Description, gs_norm = norm_key(Description), order = row_number())
  }

  selected_terms_present <- selected_map$gs_name
  grid_df <- tidyr::expand_grid(
    Cluster = clv,
    Description = selected_terms_present
  ) |>
    left_join(full_df, by = c("Cluster", "Description")) |>
    mutate(
      score = if_else(is.na(score), 0, score),
      FDR = if_else(is.na(FDR), NA_real_, FDR),
      Count = if_else(is.na(Count), 0L, as.integer(Count)),
      TermLabel = stringr::str_wrap(short_pathway_label(Description), width = 16)
    )

  dominant <- full_df |>
    group_by(Cluster) |>
    arrange(desc(score), FDR, .by_group = TRUE) |>
    slice_head(n = 1) |>
    ungroup() |>
    transmute(
      Cluster,
      Dominant = short_pathway_label(Description),
      ModuleLabel = paste0(Cluster, " - ", stringr::str_trunc(Dominant, 28))
    )
  if (!nrow(dominant)) dominant <- tibble(Cluster = clv, ModuleLabel = clv)

  grid_df <- grid_df |>
    left_join(dominant, by = "Cluster") |>
    mutate(
      ModuleLabel = if_else(is.na(ModuleLabel), as.character(Cluster), ModuleLabel),
      ModuleLabel = factor(ModuleLabel, levels = rev(dominant$ModuleLabel[match(clv, dominant$Cluster)])),
      TermLabel = factor(TermLabel, levels = unique(stringr::str_wrap(short_pathway_label(selected_terms_present), width = 16)))
    )

  readr::write_csv(grid_df, sub("\\.png$", ".csv", out_png))

  max_score <- suppressWarnings(max(grid_df$score, na.rm = TRUE))
  if (!is.finite(max_score) || max_score <= 0) max_score <- 1
  cap_score <- min(max_score, 30)
  grid_df <- grid_df |> mutate(score_plot = pmin(score, cap_score))

  plt <- ggplot(grid_df, aes(x = TermLabel, y = ModuleLabel, fill = score_plot)) +
    geom_tile(color = "white", linewidth = 0.55) +
    geom_text(aes(label = if_else(score >= 1.3, sprintf("%.1f", pmin(score, 99)), "")),
              size = 3.2, fontface = "bold", color = "grey10") +
    scale_fill_gradient(name = "-log10(FDR)", low = "grey96", high = "#2166ac", limits = c(0, cap_score)) +
    labs(
      x = NULL,
      y = NULL,
      title = "Module-by-pathway annotation",
      subtitle = "Curated program axis; full enrichment table exported"
    ) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12.5, angle = 38, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 13.5, face = "bold"),
      plot.title = element_text(face = "bold", size = 20),
      plot.subtitle = element_text(size = 12.5, color = "grey30"),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10.5),
      plot.margin = margin(6, 8, 6, 6)
    )
  save_ggplot_png_svg(plt, out_png, width = 13.4, height = 7.2)
  invisible(TRUE)
}


## ----------------------------
## Figure 7
## ----------------------------
gene_cluster_heatmap <- function(vst_mat, cd, assign_df_labeled, module_labels,
                                 baseline_vec, genes_per_module = 100,
                                 out_file = "Heatmap_Modules_Genes.png") {
  ordc <- order_samples(cd)
  mat <- vst_mat[, ordc, drop = FALSE]
  cd2 <- cd[ordc, , drop = FALSE]

  ## Defensive coercion: CLARA/PAM output and downstream joins can leave
  ## cluster-like columns as list-ish objects on some tibble/dplyr versions.
  ## Publication plotting should never fail because a grouping label is stored
  ## as a list column, so collapse to atomic vectors here and avoid bare count().
  df0 <- assign_df_labeled |>
    as_tibble() |>
    transmute(
      gene = as.character(unlist(gene, use.names = FALSE)),
      cluster = as.character(unlist(cluster, use.names = FALSE)),
      module_membership = suppressWarnings(as.numeric(unlist(module_membership, use.names = FALSE))),
      lfc_sur = suppressWarnings(as.numeric(unlist(lfc_sur, use.names = FALSE))),
      pc_loading = if ("pc_loading" %in% names(assign_df_labeled)) suppressWarnings(as.numeric(unlist(pc_loading, use.names = FALSE))) else NA_real_,
      NRI = if ("NRI" %in% names(assign_df_labeled)) suppressWarnings(as.numeric(unlist(NRI, use.names = FALSE))) else NA_real_
    ) |>
    filter(!is.na(gene), !is.na(cluster), nzchar(gene), nzchar(cluster), gene %in% rownames(mat)) |>
    mutate(
      baseline_value = suppressWarnings(as.numeric(baseline_vec[gene])),
      abs_lfc_sur = abs(lfc_sur)
    )

  if (!nrow(df0)) {
    warning("No module genes overlap the VST matrix; skipping Fig. 7A heatmap.", call. = FALSE)
    return(invisible(NULL))
  }

  module_levels <- unique(df0$cluster)
  if (all(grepl("^M[0-9]+$", module_levels))) {
    module_levels <- module_levels[order(as.integer(sub("^M", "", module_levels)))]
  } else {
    module_levels <- sort(module_levels)
  }

  pick_df <- df0 |>
    mutate(cluster = factor(cluster, levels = module_levels)) |>
    group_by(cluster) |>
    arrange(desc(abs_lfc_sur), desc(baseline_value), .by_group = TRUE) |>
    slice_head(n = genes_per_module) |>
    ungroup() |>
    arrange(cluster, desc(baseline_value), gene)

  readr::write_csv(
    pick_df |>
      select(gene, cluster, module_membership, lfc_sur, pc_loading, NRI, baseline_value, abs_lfc_sur),
    file.path(CONFIG$OUT_DIR, "Fig7A_module_heatmap_display_genes.csv")
  )

  gene_order <- pick_df$gene[pick_df$gene %in% rownames(mat)]
  gene_order <- gene_order[!duplicated(gene_order)]
  if (!length(gene_order)) return(invisible(NULL))

  sub <- mat[gene_order, , drop = FALSE]
  z <- t(scale(t(sub)))
  z[!is.finite(z)] <- 0

  baseline_band_all <- as.character(band_from_baseline(baseline_vec))
  names(baseline_band_all) <- names(baseline_vec)

  ann_row <- pick_df |>
    filter(gene %in% rownames(z)) |>
    distinct(gene, cluster, baseline_value) |>
    mutate(
      Module = as.character(cluster),
      Band = baseline_band_all[gene]
    ) |>
    select(gene, Module, Band) |>
    column_to_rownames("gene")
  ann_row <- ann_row[rownames(z), c("Module", "Band"), drop = FALSE]
  ann_col <- sample_annotation(cd2)
  ann_colors <- pheatmap_annotation_colors(ann_row, ann_col)
  write_annotation_key(ann_colors, file.path(CONFIG$OUT_DIR, "Fig7A_annotation_key.csv"))

  module_sizes <- as.integer(table(factor(pick_df$cluster, levels = module_levels)))
  module_sizes <- module_sizes[module_sizes > 0]
  gaps_row <- cumsum(module_sizes)
  gaps_row <- gaps_row[gaps_row < nrow(z)]

  ph <- pheatmap::pheatmap(
    z,
    show_rownames = FALSE,
    show_colnames = FALSE,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    gaps_row = gaps_row,
    annotation_row = ann_row,
    annotation_col = ann_col,
    main = sprintf("Seven modules: top %d genes/module", genes_per_module),
    fontsize = 13,
    fontsize_col = 10,
    fontsize_row = 5,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    annotation_colors = ann_colors,
    legend = TRUE,
    annotation_legend = TRUE,
    border_color = NA,
    silent = TRUE
  )
  g <- if (is.list(ph) && !is.null(ph$gtable)) ph$gtable else ph
  save_gtable_png_svg(g, file.path(CONFIG$OUT_DIR, out_file), width = 22.8, height = 14.4)
}
gene_set_union_heatmap <- function(vst_mat, cd, t2g, baseline_vec,
                                   max_per_set = 12,
                                   out_file = "Heatmap_GeneSets_Union.png",
                                   focus_sets = FIG7C_UNION_SETS) {
  if (!is.null(focus_sets) && length(focus_sets)) {
    focus_norm <- norm_key(focus_sets)
    t2g <- t2g |> mutate(gs_norm = norm_key(gs_name)) |> filter(gs_norm %in% focus_norm) |> select(-gs_norm)
  }
  if (!nrow(t2g)) {
    warning("No gene sets available for Fig. 7C union heatmap after filtering.")
    return(invisible(NULL))
  }

  t2g2 <- distinct(t2g, gs_name, gene_symbol) |>
    mutate(bl = baseline_vec[gene_symbol]) |>
    filter(!is.na(bl)) |>
    group_by(gs_name) |>
    arrange(desc(bl), .by_group = TRUE) |>
    slice_head(n = max_per_set) |>
    ungroup() |>
    select(-bl)
  genes <- intersect(unique(t2g2$gene_symbol), rownames(vst_mat))
  if (!length(genes)) return(invisible(NULL))

  ## Source table keeps the exact compact gene-set union used for the panel.
  readr::write_csv(
    t2g2 |> filter(gene_symbol %in% genes) |> arrange(gs_name, gene_symbol),
    file.path(CONFIG$OUT_DIR, "Fig7C_gene_set_union_membership.csv")
  )

  ordc <- order_samples(cd)
  mat <- vst_mat[genes, ordc, drop = FALSE]
  cd2 <- cd[ordc, , drop = FALSE]
  mat <- mat[order(-baseline_vec[rownames(mat)]), , drop = FALSE]
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0

  baseline_band_all <- as.character(band_from_baseline(baseline_vec))
  names(baseline_band_all) <- names(baseline_vec)

  g2s <- t2g2 |>
    group_by(gene_symbol) |>
    slice_head(n = 1) |>
    ungroup() |>
    transmute(
      gene = gene_symbol,
      TopSet = figure7c_set_label(gs_name),
      Band = baseline_band_all[gene_symbol]
    ) |>
    column_to_rownames("gene")
  ann_row <- g2s[rownames(z), c("TopSet", "Band"), drop = FALSE]
  ann_col <- sample_annotation(cd2)
  ann_colors <- pheatmap_annotation_colors(ann_row, ann_col)
  write_annotation_key(ann_colors, file.path(CONFIG$OUT_DIR, "Fig7C_annotation_key.csv"))

  ph <- pheatmap::pheatmap(
    z,
    show_rownames = FALSE,
    show_colnames = FALSE,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    annotation_row = ann_row,
    annotation_col = ann_col,
    main = "Presorted gene-set union: shared 4 h Ctrl order",
    fontsize = 13,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    annotation_colors = ann_colors,
    annotation_legend = TRUE,
    legend = TRUE,
    border_color = NA,
    silent = TRUE
  )
  g <- if (is.list(ph) && !is.null(ph$gtable)) ph$gtable else ph
  save_gtable_png_svg(g, file.path(CONFIG$OUT_DIR, out_file), width = 23.8, height = 14.2)
}

gsea_time_aligned_axis <- function(load_tbl, out_csv, out_png) {
  ranks <- collapse_named_stat(load_tbl$loading, load_tbl$gene, mode = "maxabs")
  hallmark <- msigdbr_get("H") |>
    select(gs_name, gene_symbol) |>
    distinct()
  genesets <- split(hallmark$gene_symbol, hallmark$gs_name)
  fg <- fgsea::fgsea(
    genesets,
    ranks,
    nperm = CONFIG$FGSEA_NPERM,
    minSize = CONFIG$FGSEA_MIN_SIZE,
    maxSize = CONFIG$FGSEA_MAX_SIZE
  )
  readr::write_csv(flatten_fgsea_for_csv(fg), out_csv)
  fg_tbl <- as_tibble(fg)
  top_pos <- fg_tbl |> filter(NES >= 0) |> arrange(padj, desc(abs(NES))) |> slice_head(n = 8)
  top_neg <- fg_tbl |> filter(NES < 0) |> arrange(padj, desc(abs(NES))) |> slice_head(n = 8)
  top <- bind_rows(top_neg, top_pos) |> distinct(pathway, .keep_all = TRUE)
  if (!nrow(top)) top <- fg_tbl |> arrange(padj) |> slice_head(n = min(16, nrow(fg_tbl)))
  if (nrow(top)) {
    top <- top |>
      mutate(
        Direction = if_else(NES >= 0, "Positive PC pole", "Negative PC pole"),
        pathway_label = short_pathway_label(pathway)
      )
    plt <- ggplot(top, aes(reorder(pathway_label, NES), NES, fill = Direction)) +
      geom_col(width = 0.78) +
      geom_hline(yintercept = 0, linewidth = 0.4, color = "grey35") +
      coord_flip() +
      labs(
        x = NULL,
        y = "NES on distortion-associated PC",
        title = "Hallmark polarity along distortion-associated PC",
        subtitle = "Signed Hallmark NES on the distortion-associated PC",
        fill = NULL
      ) +
      theme_minimal(base_size = 20) +
      theme(
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 17),
        plot.title = element_text(face = "bold", size = 22),
        plot.subtitle = element_text(size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 14)
      )
    save_ggplot_png_svg(plt, out_png, width = 13.2, height = 8.8)
  }
}

## ----------------------------
## Figure 8
## ----------------------------
get_samples_for_group <- function(cd, time, stress, induction) {
  rownames(cd)[cd$Time == time & cd$Stress == stress & cd$Induction == induction]
}


make_rank_vector <- function(vst_mat, cd, time, stress, numerator, denominator) {
  s_num <- get_samples_for_group(cd, time, stress, numerator)
  s_den <- get_samples_for_group(cd, time, stress, denominator)
  if (!length(s_num) || !length(s_den)) return(NULL)
  v_num <- rowMeans(vst_mat[, s_num, drop = FALSE], na.rm = TRUE)
  v_den <- rowMeans(vst_mat[, s_den, drop = FALSE], na.rm = TRUE)
  ranks <- v_num - v_den
  ranks <- ranks[is.finite(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  list(
    ranks = ranks,
    rank_stat = "vst_mean_diff",
    samples_num = s_num,
    samples_den = s_den
  )
}


leading_edge_head_fraction <- function(leading_edge, head_genes) {
  if (is.null(leading_edge) || !length(leading_edge)) return(NA_real_)
  mean(leading_edge %in% head_genes)
}

run_fgsea_contrasts <- function(vst_mat, cd, genesets, contrast_defs,
                                nperm = 10000, minSize = 10, maxSize = 500,
                                head_genes = NULL) {
  times <- resolve_time_levels(na.omit(unique(cd$Time)))
  stresses <- na.omit(unique(cd$Stress))
  stresses <- stresses[order(match(stresses, c("NoCPT", "CPT", "CPT_low", "CPT_high")), na.last = TRUE)]
  used <- list()
  res_all <- list()
  idx <- 0L
  for (i in seq_len(nrow(contrast_defs))) {
    con <- contrast_defs$Contrast[i]
    num <- contrast_defs$Numerator[i]
    den <- contrast_defs$Denominator[i]
    for (t in times) {
      for (s in stresses) {
        rk <- make_rank_vector(vst_mat, cd, t, s, num, den)
        if (is.null(rk)) next
        idx <- idx + 1L
        used[[idx]] <- tibble(
          Contrast = con,
          Time = t,
          Stress = s,
          Numerator = num,
          Denominator = den,
          N_num = length(rk$samples_num),
          N_den = length(rk$samples_den),
          RankStatistic = rk$rank_stat,
          Samples_num = paste(rk$samples_num, collapse = ";"),
          Samples_den = paste(rk$samples_den, collapse = ";")
        )
        fg <- fgsea::fgsea(
          genesets,
          rk$ranks,
          nperm = nperm,
          minSize = minSize,
          maxSize = maxSize
        )
        fg <- as_tibble(fg) |>
          mutate(
            Contrast = con,
            Time = t,
            Stress = s,
            Numerator = num,
            Denominator = den,
            N_num = length(rk$samples_num),
            N_den = length(rk$samples_den),
            RankStatistic = rk$rank_stat,
            sig = -log10(padj + 1e-300),
            LeadingEdgeSize = lengths(leadingEdge)
          )
        if (!is.null(head_genes)) {
          fg <- fg |>
            mutate(LeadingEdgeHeadFrac = vapply(
              leadingEdge,
              function(x) leading_edge_head_fraction(x, head_genes),
              FUN.VALUE = numeric(1)
            ))
        } else {
          fg <- fg |> mutate(LeadingEdgeHeadFrac = NA_real_)
        }
        res_all[[length(res_all) + 1L]] <- fg
      }
    }
  }
  list(
    results = bind_rows(res_all),
    contrasts_used = bind_rows(used)
  )
}

prep_panel_df <- function(fg_df) {
  str_levels <- expand.grid(
    TimeShort = c("4h", "24h"),
    StressShort = c("NoCPT", "CPT", "LCPT", "HCPT"),
    stringsAsFactors = FALSE
  ) |>
    mutate(Stratum = paste(TimeShort, StressShort, sep = "_")) |>
    pull(Stratum)
  fg_df |>
    mutate(
      TimeShort = short_time(Time),
      StressShort = short_stress(Stress),
      Stratum = paste(TimeShort, StressShort, sep = "_"),
      Stratum = factor(Stratum, levels = str_levels),
      Contrast = factor(Contrast, levels = FIG8_CONTRASTS$Contrast),
      ContrastLabel = dplyr::recode(as.character(Contrast),
                                    D_vs_Ctrl = "D vs Ctrl",
                                    DT_vs_D = "DT vs D",
                                    Tam_vs_Ctrl = "Tam vs Ctrl",
                                    .default = as.character(Contrast)),
      ContrastLabel = factor(ContrastLabel, levels = c("D vs Ctrl", "DT vs D", "Tam vs Ctrl"))
    )
}

plot_dotmap <- function(df, pathways, title, out_png,
                        facet_by_contrast = FALSE,
                        color_var = "NES",
                        size_var = "sig",
                        color_name = "NES",
                        size_name = "-log10(FDR)",
                        wrap_width = 30) {
  pathways <- unique(as.character(pathways))
  pathways <- pathways[pathways %in% unique(as.character(df$pathway))]
  if (!length(pathways)) return(invisible(NULL))

  label_map <- unique_pathway_display_labels(pathways, wrap_width = wrap_width)
  df2 <- df |>
    filter(pathway %in% pathways) |>
    mutate(
      pathway = as.character(pathway),
      pathway_label = unname(label_map[pathway]),
      pathway_label = factor(pathway_label, levels = rev(unname(label_map[pathways])))
    )
  if (!nrow(df2)) return(invisible(NULL))

  p <- ggplot(df2, aes(x = Stratum, y = pathway_label)) +
    geom_point(aes(size = .data[[size_var]], color = .data[[color_var]]), alpha = 0.94) +
    scale_size_continuous(name = size_name, range = c(3.0, 8.4)) +
    scale_color_continuous(name = color_name) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 17) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 12.5),
      axis.text.y = element_text(size = 12.8, lineheight = 0.9),
      strip.text = element_text(size = 15, face = "bold"),
      plot.title = element_text(face = "bold", size = 18, margin = margin(b = 4)),
      legend.title = element_text(size = 11.5),
      legend.text = element_text(size = 10.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.35, color = "grey88"),
      plot.margin = margin(8, 14, 8, 10)
    )
  if (facet_by_contrast) {
    ## Free-x facets remove empty CPT/time columns from contrasts that are
    ## genuinely absent (e.g. D_vs_Ctrl controls may only exist for NoCPT),
    ## while keeping the observed strata explicit.
    p <- p + facet_grid(. ~ ContrastLabel, scales = "free_x", space = "free_x") +
      scale_x_discrete(drop = TRUE)
  }
  h <- max(ifelse(facet_by_contrast, 5.0, 4.2), 0.50 * length(pathways) + 2.6)
  w <- ifelse(facet_by_contrast, 15.8, 12.4)
  save_ggplot_png_svg(p, out_png, width = w, height = h)
  invisible(TRUE)
}


make_mito_contrast_layout <- function(df) {
  df |>
    filter(Contrast %in% c("D_vs_Ctrl", "DT_vs_D")) |>
    mutate(
      TimeLabel = if_else(Time == "4 hours", "4 h", if_else(Time == "24 hours", "24 h", as.character(Time))),
      CPTLabel = case_when(
        Stress == "NoCPT" ~ "no CPT",
        Stress == "CPT_low" ~ "LCPT",
        Stress == "CPT_high" ~ "HCPT",
        Stress == "CPT" ~ "CPT",
        TRUE ~ as.character(Stress)
      ),
      ComparisonBlock = case_when(
        Contrast == "D_vs_Ctrl" ~ "Dox baseline",
        Contrast == "DT_vs_D" ~ "MYC-ER activation\nDox held constant",
        TRUE ~ as.character(Contrast)
      ),
      ExactContrast = case_when(
        Contrast == "D_vs_Ctrl" ~ "D vs Ctrl",
        Contrast == "DT_vs_D" & Stress == "NoCPT" ~ "DT vs D",
        Contrast == "DT_vs_D" & Stress == "CPT_low" ~ "DT_L_CPT vs D_L_CPT",
        Contrast == "DT_vs_D" & Stress == "CPT_high" ~ "DT_H_CPT vs D_H_CPT",
        Contrast == "DT_vs_D" & Stress == "CPT" ~ "DT_CPT vs D_CPT",
        TRUE ~ as.character(Contrast)
      ),
      ColumnLabel = paste(TimeLabel, CPTLabel, sep = "\n"),
      KeepForMitoControl = case_when(
        Contrast == "D_vs_Ctrl" & Stress == "NoCPT" ~ TRUE,
        Contrast == "DT_vs_D" & Stress %in% c("NoCPT", "CPT_low", "CPT_high", "CPT") ~ TRUE,
        TRUE ~ FALSE
      )
    ) |>
    filter(KeepForMitoControl) |>
    mutate(
      ComparisonBlock = factor(ComparisonBlock, levels = c("Dox baseline", "MYC-ER activation\nDox held constant")),
      ColumnLabel = factor(ColumnLabel, levels = c("4 h\nno CPT", "24 h\nno CPT", "4 h\nLCPT", "24 h\nLCPT", "4 h\nHCPT", "24 h\nHCPT", "4 h\nCPT", "24 h\nCPT"))
    )
}

plot_mito_matched_contrast_dotmap <- function(df, pathways, title, out_png,
                                              color_var = "NES",
                                              size_var = "sig",
                                              color_name = "NES",
                                              size_name = "-log10(FDR)",
                                              wrap_width = 31) {
  pathways <- unique(as.character(pathways))
  pathways <- pathways[pathways %in% unique(as.character(df$pathway))]
  if (!length(pathways)) return(invisible(NULL))

  label_map <- unique_pathway_display_labels(pathways, wrap_width = wrap_width)
  df2 <- make_mito_contrast_layout(df) |>
    filter(pathway %in% pathways) |>
    mutate(
      pathway = as.character(pathway),
      pathway_label = unname(label_map[pathway]),
      pathway_label = factor(pathway_label, levels = rev(unname(label_map[pathways])))
    )
  if (!nrow(df2)) return(invisible(NULL))

  layout_csv <- sub("\\.png$", "_displayed_contrast_layout.csv", out_png)
  readr::write_csv(
    df2 |>
      distinct(ComparisonBlock, ColumnLabel, Contrast, Time, Stress, ExactContrast) |>
      arrange(ComparisonBlock, ColumnLabel, ExactContrast),
    layout_csv
  )

  p <- ggplot(df2, aes(x = ColumnLabel, y = pathway_label)) +
    geom_point(aes(size = .data[[size_var]], color = .data[[color_var]]), alpha = 0.94) +
    facet_grid(. ~ ComparisonBlock, scales = "free_x", space = "free_x") +
    scale_x_discrete(drop = TRUE, guide = ggplot2::guide_axis(n.dodge = 2)) +
    scale_size_continuous(name = size_name, range = c(3.0, 8.6)) +
    scale_color_continuous(name = color_name) +
    labs(
      x = NULL,
      y = NULL,
      title = title
    ) +
    theme_minimal(base_size = 17) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 10.8, lineheight = 0.88),
      axis.text.y = element_text(size = 12.8, lineheight = 0.9),
      strip.text.x = element_text(size = 15.0, face = "bold", lineheight = 0.95, margin = margin(t = 6, b = 8)),
      plot.title = element_text(face = "bold", size = 18, margin = margin(b = 8)),
      plot.caption = element_text(size = 10.2, color = "grey28", hjust = 0, margin = margin(t = 8)),
      legend.title = element_text(size = 11.5),
      legend.text = element_text(size = 10.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.35, color = "grey88"),
      panel.spacing.x = unit(1.2, "lines"),
      plot.margin = margin(26, 66, 58, 72)
    )
  h <- max(5.2, 0.52 * length(pathways) + 2.9)
  save_ggplot_png_svg(p, out_png, width = 20.8, height = h + 1.15)
  invisible(TRUE)
}


plot_es_small_multiples <- function(vst_mat, cd, genesets, contrast_defs,
                                    contrast_name, stress = "NoCPT",
                                    sets = FIG8C_ES_SETS,
                                    out_png = "Fig8C_ES_curves.png") {
  con_row <- contrast_defs[contrast_defs$Contrast == contrast_name, , drop = FALSE]
  if (!nrow(con_row)) return(invisible(NULL))
  num <- con_row$Numerator[[1]]
  den <- con_row$Denominator[[1]]
  times <- intersect(c("4 hours", "24 hours"), unique(cd$Time))
  if (length(times) < 2) return(invisible(NULL))
  plots <- list()
  for (gs in sets) {
    if (!gs %in% names(genesets)) next
    for (t in times) {
      rk <- make_rank_vector(vst_mat, cd, t, stress, num, den)
      if (is.null(rk)) {
        plots[[length(plots) + 1L]] <- ggplot() + theme_void() +
          ggtitle(sprintf("Missing: %s | %s | %s", gs, short_time(t), stress))
        next
      }
      p <- fgsea::plotEnrichment(genesets[[gs]], rk$ranks) +
        labs(title = sprintf("%s | %s | %s", short_pathway_label(gs), short_time(t), short_stress(stress)),
             x = "Rank", y = "ES") +
        theme_minimal(base_size = 16) +
        theme(
          plot.title = element_text(size = 13.5, face = "bold", margin = margin(b = 3)),
          axis.text = element_text(size = 11.5),
          axis.title = element_text(size = 12.5),
          panel.grid = element_line(linewidth = 0.25, color = "grey90"),
          plot.margin = margin(6, 8, 6, 6)
        )
      plots[[length(plots) + 1L]] <- p
    }
  }
  if (!length(plots)) return(invisible(NULL))
  g <- gridExtra::arrangeGrob(grobs = plots, ncol = length(times))
  save_grob_png_svg(g, file.path(CONFIG$OUT_DIR, out_png), width = 13.8, height = 8.4)
  invisible(TRUE)
}

## ----------------------------
## Figure composition
## ----------------------------
save_grob_png <- function(g, file, w = 16, h = 12, res = NULL) {
  save_grob_png_svg(g, file, width = w, height = h, res = res)
}

# Composite SVG files made from imported PNG panels are fragile in Illustrator
# and can render as clipped/duplicated raster fragments.  For final composites,
# write PNG plus PDF only; panel-level PDF files remain available for vector
# editing and manuscript layout.
save_composite_png_pdf <- function(g, file, w = 16, h = 12, res = NULL) {
  ## Final composite writer with explicit canvas padding. The padding prevents
  ## panel labels, facet headers, captions, and pheatmap legends from being cut
  ## at the artboard edge in PNG/PDF/SVG export.
  if (is.null(res)) res <- CONFIG$PLOT_DPI
  pad <- if (!is.null(CONFIG$COMPOSITE_CANVAS_PADDING)) CONFIG$COMPOSITE_CANVAS_PADDING else 0.018
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  grDevices::png(filename = file, width = w, height = h, units = "in", res = res)
  draw_grid_object_with_padding(g, padding = pad)
  grDevices::dev.off()
  saveRDS(g, sub("[.]png$", ".grob.rds", file))
  if (isTRUE(CONFIG$WRITE_SVG)) {
    out_svg <- png_to_svg(file)
    tryCatch({
      open_svg_device(out_svg, width = w, height = h)
      draw_grid_object_with_padding(g, padding = pad)
      grDevices::dev.off()
    }, error = function(e) {
      message(sprintf("[Plot] Composite SVG export failed for %s: %s", basename(file), conditionMessage(e)))
      try(grDevices::dev.off(), silent = TRUE)
    })
  }
  if (isTRUE(CONFIG$WRITE_PDF)) {
    out_pdf <- png_to_pdf(file)
    tryCatch({
      grDevices::pdf(file = out_pdf, width = w, height = h, useDingbats = FALSE)
      draw_grid_object_with_padding(g, padding = pad)
      grDevices::dev.off()
    }, error = function(e) {
      message(sprintf("[Plot] Composite PDF export failed for %s: %s", basename(file), conditionMessage(e)))
      try(grDevices::dev.off(), silent = TRUE)
    })
  }
  invisible(TRUE)
}

xml_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

panel_href_for_svg <- function(outdir, panel_png) {
  panel_svg <- png_to_svg(file.path(outdir, panel_png))
  panel_png_path <- file.path(outdir, panel_png)
  if (file.exists(panel_svg)) return(basename(panel_svg))
  if (file.exists(panel_png_path)) return(basename(panel_png_path))
  NA_character_
}

write_linked_svg_composite <- function(out_svg, outdir, panel_files, labels,
                                       layout_matrix, widths, heights,
                                       width_in, height_in,
                                       label_size = 34,
                                       padding = 0.012) {
  ## SVG-safe composite writer. This creates an SVG assembly that links the
  ## panel-level SVG files (or PNG fallback) as whole objects. It avoids the
  ## broken R rasterGrob-on-SVG pathway while still providing an SVG artboard
  ## that opens cleanly in Illustrator/Inkscape. Keep the panel SVG/PNG files in
  ## the same directory as the composite SVG.
  ## Deprecated fallback retained for compatibility. Composite SVGs are now
  ## written as true vector grid output from the saved panel grobs/RDS sidecars.
  ## Do not overwrite them with linked/raster assemblies.
  return(invisible(FALSE))
  if (!isTRUE(CONFIG$WRITE_SVG)) return(invisible(FALSE))
  dir.create(dirname(out_svg), showWarnings = FALSE, recursive = TRUE)
  U <- 100
  W <- width_in * U
  H <- height_in * U
  nr <- nrow(layout_matrix)
  nc <- ncol(layout_matrix)
  widths <- rep(widths, length.out = nc)
  heights <- rep(heights, length.out = nr)
  col_w <- W * widths / sum(widths)
  row_h <- H * heights / sum(heights)
  x0 <- c(0, cumsum(col_w))[seq_len(nc)]
  y0 <- c(0, cumsum(row_h))[seq_len(nr)]
  ids <- sort(unique(as.vector(layout_matrix)))
  ids <- ids[!is.na(ids)]
  lines <- c(
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="%.3fin" height="%.3fin" viewBox="0 0 %.3f %.3f">', width_in, height_in, W, H),
    '<rect x="0" y="0" width="100%" height="100%" fill="white"/>'
  )
  manifest <- list()
  for (id in ids) {
    pos <- which(layout_matrix == id, arr.ind = TRUE)
    r1 <- min(pos[, "row"]); r2 <- max(pos[, "row"])
    c1 <- min(pos[, "col"]); c2 <- max(pos[, "col"])
    x <- x0[c1]
    y <- y0[r1]
    w <- sum(col_w[c1:c2])
    h <- sum(row_h[r1:r2])
    nm <- names(panel_files)[id]
    if (is.na(nm) || !nzchar(nm)) nm <- as.character(id)
    href <- panel_href_for_svg(outdir, panel_files[[id]])
    lab <- labels[[id]]
    if (is.na(href)) {
      lines <- c(lines, sprintf('<rect x="%.3f" y="%.3f" width="%.3f" height="%.3f" fill="#f6f6f6" stroke="#999"/>', x, y, w, h),
                 sprintf('<text x="%.3f" y="%.3f" font-family="Arial, Helvetica, sans-serif" font-size="24" fill="#555">Missing: %s</text>', x + 20, y + 40, xml_escape(panel_files[[id]])))
    } else {
      lines <- c(lines,
                 sprintf('<image x="%.3f" y="%.3f" width="%.3f" height="%.3f" href="%s" xlink:href="%s" preserveAspectRatio="xMidYMid meet"/>', x, y, w, h, xml_escape(href), xml_escape(href)),
                 sprintf('<text x="%.3f" y="%.3f" font-family="Arial, Helvetica, sans-serif" font-size="%.1f" font-weight="700" fill="#000000">%s</text>', x + W * padding, y + H * padding + label_size * 0.35, label_size, xml_escape(lab)))
      manifest[[length(manifest) + 1L]] <- data.frame(panel = nm, label = lab, source = href, x = x, y = y, width = w, height = h)
    }
  }
  lines <- c(lines, '</svg>')
  writeLines(lines, out_svg, useBytes = TRUE)
  if (length(manifest)) {
    readr::write_csv(bind_rows(manifest), sub('[.]svg$', '_svg_manifest.csv', out_svg))
  }
  invisible(TRUE)
}


panel_letter_grob <- function(g, label, fontsize = CONFIG$COMPOSITE_PANEL_LETTER_SIZE) {
  ## Put panel letters in a dedicated header row instead of overlaying them on
  ## the panel. Overlay labels were clipping/obscuring facet headers such as
  ## "Dox baseline" in Fig. 8D and can also create edge-cut artifacts in SVG.
  label_g <- grid::textGrob(
    label,
    x = grid::unit(0, "npc"),
    y = grid::unit(0.45, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontsize = fontsize, fontface = "bold")
  )
  gridExtra::arrangeGrob(
    label_g,
    g,
    ncol = 1,
    heights = grid::unit.c(grid::unit(max(0.24, fontsize / 72 * 1.05), "in"), grid::unit(1, "null"))
  )
}

read_panel_png <- function(fp) {
  ## Historical name kept for compatibility. Prefer vector/editable panel
  ## sidecars generated during panel plotting; fall back to PNG only if needed.
  grob_rds <- sub("[.]png$", ".grob.rds", fp)
  plot_rds <- sub("[.]png$", ".plot.rds", fp)
  if (file.exists(grob_rds)) return(readRDS(grob_rds))
  if (file.exists(plot_rds)) return(ggplot2::ggplotGrob(readRDS(plot_rds)))
  if (!file.exists(fp)) return(textGrob(sprintf("Missing: %s", basename(fp))))
  warning(sprintf("Composing %s from PNG fallback because no panel RDS sidecar was found; rerun panel generation for fully vector SVG.", basename(fp)), call. = FALSE)
  rasterGrob(png::readPNG(fp), interpolate = TRUE, width = unit(1, "npc"), height = unit(1, "npc"))
}

compose_figure7 <- function(outdir, include_panel_c = TRUE) {
  files4 <- c(
    A = "Heatmap_Modules_Genes.png",
    B = "MSIGDB_modulesXterms_heatmap.png",
    C = "Heatmap_GeneSets_Union.png",
    D = "GSEA_time_aligned_axis_Hallmark_barplot.png"
  )
  panel4 <- lapply(names(files4), function(nm) panel_letter_grob(read_panel_png(file.path(outdir, files4[[nm]])), nm, fontsize = 28))
  names(panel4) <- names(files4)

  ## Publication layout for the R2 manuscript: the two heatmaps are given
  ## full-width lanes so their annotation keys remain attached and legible;
  ## B and D are paired as the interpretive middle row.
  g4 <- arrangeGrob(
    grobs = list(panel4$A, panel4$B, panel4$C, panel4$D),
    layout_matrix = rbind(c(1, 1), c(2, 4), c(3, 3)),
    widths = c(1.0, 1.0),
    heights = c(1.15, 0.80, 1.40)
  )
  if (isTRUE(include_panel_c)) {
    save_composite_png_pdf(
      g4,
      file.path(outdir, "Figure7_Composite.png"),
      w = CONFIG$FIG7_COMPOSITE_WIDTH,
      h = CONFIG$FIG7_COMPOSITE_HEIGHT,
      res = 300
    )
    write_linked_svg_composite(
      file.path(outdir, "Figure7_Composite.svg"), outdir,
      panel_files = files4,
      labels = c(A = "A", B = "B", C = "C", D = "D"),
      layout_matrix = rbind(c(1, 1), c(2, 4), c(3, 3)),
      widths = c(1.0, 1.0),
      heights = c(1.15, 0.80, 1.40),
      width_in = CONFIG$FIG7_COMPOSITE_WIDTH,
      height_in = CONFIG$FIG7_COMPOSITE_HEIGHT,
      label_size = 34
    )
  } else {
    save_composite_png_pdf(
      g4,
      file.path(outdir, "Figure7_Composite_withPanelC.png"),
      w = CONFIG$FIG7_COMPOSITE_WIDTH,
      h = CONFIG$FIG7_COMPOSITE_HEIGHT,
      res = 300
    )
    write_linked_svg_composite(
      file.path(outdir, "Figure7_Composite_withPanelC.svg"), outdir,
      panel_files = files4,
      labels = c(A = "A", B = "B", C = "C", D = "D"),
      layout_matrix = rbind(c(1, 1), c(2, 4), c(3, 3)),
      widths = c(1.0, 1.0),
      heights = c(1.15, 0.80, 1.40),
      width_in = CONFIG$FIG7_COMPOSITE_WIDTH,
      height_in = CONFIG$FIG7_COMPOSITE_HEIGHT,
      label_size = 34
    )
  }

  panel3 <- panel4[c("A", "B", "D")]
  g3 <- arrangeGrob(
    panel3$A,
    arrangeGrob(panel3$B, panel3$D, ncol = 2, widths = c(1.02, 1.05)),
    heights = c(1.16, 1.0),
    ncol = 1
  )
  save_composite_png_pdf(g3, file.path(outdir, "Figure7_ABD_Readability_Composite.png"), w = 19.2, h = 13.4)
  write_linked_svg_composite(
    file.path(outdir, "Figure7_ABD_Readability_Composite.svg"), outdir,
    panel_files = files3,
    labels = c(A = "A", B = "B", D = "D"),
    layout_matrix = rbind(c(1, 1), c(2, 3)),
    widths = c(1.02, 1.05),
    heights = c(1.16, 1.0),
    width_in = 19.2,
    height_in = 13.4,
    label_size = 32
  )
  if (!isTRUE(include_panel_c)) {
    file.copy(file.path(outdir, "Figure7_ABD_Readability_Composite.png"), file.path(outdir, "Figure7_Composite.png"), overwrite = TRUE)
    file.copy(file.path(outdir, "Figure7_ABD_Readability_Composite.svg"), file.path(outdir, "Figure7_Composite.svg"), overwrite = TRUE)
  }
}

compose_figure8 <- function(outdir) {
  files <- c(
    A = "Fig8A_MYC_targets_dotmap.png",
    B = "Fig8B_DDR_dotmap.png",
    C = "Fig8C_ES_curves.png",
    D = "Fig8D_MitoControl_dotmap.png"
  )
  panel <- lapply(names(files), function(nm) panel_letter_grob(read_panel_png(file.path(outdir, files[[nm]])), nm, fontsize = 30))
  names(panel) <- names(files)
  ## Figure 8 is a quantitative-panel figure, so reduce dead space and let
  ## panel C/D occupy a slightly larger row than A/B.
  g <- arrangeGrob(
    grobs = list(panel$A, panel$B, panel$C, panel$D),
    layout_matrix = rbind(c(1, 2), c(3, 4)),
    widths = c(1.05, 1.03),
    heights = c(0.82, 1.30)
  )
  save_composite_png_pdf(g, file.path(outdir, "Figure8_Composite.png"), w = CONFIG$FIG8_COMPOSITE_WIDTH, h = CONFIG$FIG8_COMPOSITE_HEIGHT)
  write_linked_svg_composite(
    file.path(outdir, "Figure8_Composite.svg"), outdir,
    panel_files = files,
    labels = c(A = "A", B = "B", C = "C", D = "D"),
    layout_matrix = rbind(c(1, 2), c(3, 4)),
    widths = c(1.05, 1.03),
    heights = c(0.82, 1.30),
    width_in = CONFIG$FIG8_COMPOSITE_WIDTH,
    height_in = CONFIG$FIG8_COMPOSITE_HEIGHT,
    label_size = 34
  )
}

compose_supplementary_figure3 <- function(outdir) {
  fp <- file.path(outdir, "SuppFig3_MitoControl_LeadingEdgeHeadFrac.png")
  if (!file.exists(fp)) return(invisible(NULL))
  ## The supplementary panel already carries the correct title and legends.
  ## Do not wrap it in a second title, which caused clipping/duplicate headers.
  targets <- c(
    png = file.path(outdir, "SupplementaryFigure3_Composite.png"),
    pdf = file.path(outdir, "SupplementaryFigure3_Composite.pdf")
  )
  file.copy(fp, targets[["png"]], overwrite = TRUE)
  src_svg <- png_to_svg(fp)
  if (file.exists(src_svg)) file.copy(src_svg, file.path(outdir, "SupplementaryFigure3_Composite.svg"), overwrite = TRUE)
  src_pdf <- png_to_pdf(fp)
  if (file.exists(src_pdf)) file.copy(src_pdf, targets[["pdf"]], overwrite = TRUE)
  invisible(TRUE)
}


write_run_manifest <- function(config, outdir, counts_mat, meta_df) {
  lines <- c(
    "# Fig7/Fig8 program-level reweighting run manifest",
    sprintf("timestamp: %s", format(Sys.time(), "%F %T")),
    sprintf("working_directory: %s", normalizePath(getwd())),
    sprintf("script_path: %s", SCRIPT_PATH),
    sprintf("module_dir: %s", MODULE_DIR),
    sprintf("repo_root: %s", REPO_ROOT),
    sprintf("counts_csv: %s", config$COUNTS_CSV),
    sprintf("vst_csv: %s", config$VST_CSV),
    sprintf("sample_sheet_csv: %s", config$SAMPLE_SHEET_CSV),
    sprintf("out_dir: %s", normalizePath(outdir)),
    sprintf("gene_col: %s", config$GENE_COL),
    sprintf("gene_id_type: %s", config$GENE_ID_TYPE),
    sprintf("count_input_mode: %s", config$COUNT_INPUT_MODE),
    sprintf("axis_mode: %s", config$AXIS_MODE),
    sprintf("baseline_rank_source: %s", config$BASELINE_RANK_SOURCE),
    sprintf("module_k_range: %s-%s", config$MODULE_K_MIN, config$MODULE_K_MAX),
    sprintf("module_k_final: %s", config$K_MODULES),
    sprintf("module_k_rule: %s", config$MODULE_K_RULE),
    sprintf("require_seven_modules: %s", config$REQUIRE_SEVEN_MODULES),
    sprintf("min_total_counts: %s", config$MIN_TOTAL_COUNTS),
    sprintf("head_frac: %s", config$HEAD_FRAC),
    sprintf("fgsea_nperm: %s", config$FGSEA_NPERM),
    sprintf("plot_dpi: %s", config$PLOT_DPI),
    sprintf("write_svg: %s", config$WRITE_SVG),
    sprintf("write_pdf: %s", config$WRITE_PDF),
    sprintf("counts_dim: %d x %d", nrow(counts_mat), ncol(counts_mat)),
    sprintf("sample_count: %d", nrow(meta_df)),
    "",
    "# sessionInfo()",
    capture.output(sessionInfo())
  )
  writeLines(lines, con = file.path(outdir, "run_manifest.txt"))
}

## ----------------------------
## Main
## ----------------------------
run_pipeline <- function() {
  parse_cli_args()
  finalize_config_paths()
  check_required_packages(required_packages)
  dir.create(CONFIG$OUT_DIR, showWarnings = FALSE, recursive = TRUE)

  cat(sprintf("[ %s ] Starting %s in %s\n",
              format(Sys.time(), "%F %T"), "run.R", normalizePath(getwd())))
  cat(sprintf("[Config] counts=%s | vst=%s | sample_sheet=%s | outdir=%s | count_mode=%s | nperm=%s | fig7=%s | fig7C=%s | fig8=%s | mito=%s | dpi=%s | pdf=%s\n",
              CONFIG$COUNTS_CSV, CONFIG$VST_CSV, CONFIG$SAMPLE_SHEET_CSV, CONFIG$OUT_DIR,
              CONFIG$COUNT_INPUT_MODE, CONFIG$FGSEA_NPERM, CONFIG$RUN_FIG7, CONFIG$FIG7_INCLUDE_PANEL_C, CONFIG$RUN_FIG8, CONFIG$RUN_MITO_CONTROL, CONFIG$PLOT_DPI, CONFIG$WRITE_PDF))
  cat(sprintf("[Config] module_engine=%s | module_k=%s | silhouette_sample_n=%s | clara_samples=%s | clara_sampsize=%s\n",
              CONFIG$MODULE_CLUSTER_ENGINE, CONFIG$K_MODULES, CONFIG$MODULE_SILHOUETTE_SAMPLE_N, CONFIG$MODULE_CLARA_SAMPLES, CONFIG$MODULE_CLARA_SAMPSIZE))

  pr <- make_progress(15)

  pr$step("Read sample sheet")
  meta <- read_sample_sheet(CONFIG$SAMPLE_SHEET_CSV)
  readr::write_csv(as_tibble(meta), file.path(CONFIG$OUT_DIR, "sample_meta.csv"))
  meta_counts <- dplyr::count(as_tibble(meta), Time, Stress, Induction, sort = TRUE)
  readr::write_csv(meta_counts, file.path(CONFIG$OUT_DIR, "sample_meta_counts.csv"))

  pr$step("Write parser QC against sample sheet")
  write_sample_parse_qc(meta, CONFIG$OUT_DIR)

  use_precomputed_vst <- file.exists(CONFIG$VST_CSV) && !file.exists(CONFIG$COUNTS_CSV)

  counts <- matrix(nrow = 0, ncol = 0)
  counts_info <- NULL
  vst_all <- NULL
  if (isTRUE(use_precomputed_vst)) {
    pr$step("Read precomputed VST matrix")
    message(sprintf("[VST] Using precomputed VST matrix because count matrix is absent: %s", CONFIG$VST_CSV))
    vst_all <- read_precomputed_vst(CONFIG$VST_CSV, meta, gene_col = CONFIG$GENE_COL)
    counts_info <- list(gene_ids = vst_all$gene_ids, gene_col = vst_all$gene_col)
    counts <- matrix(NA_integer_, nrow = nrow(vst_all$mat), ncol = ncol(vst_all$mat), dimnames = dimnames(vst_all$mat))
    readr::write_csv(
      tibble(filter = "precomputed_vst_input", genes_before = nrow(vst_all$mat), genes_after = nrow(vst_all$mat), genes_removed = NA_integer_),
      file.path(CONFIG$OUT_DIR, "gene_count_filter_summary.csv")
    )
  } else {
    pr$step("Read counts")
    counts_info <- read_counts(CONFIG$COUNTS_CSV, meta, gene_col = CONFIG$GENE_COL,
                               min_total_counts = CONFIG$MIN_TOTAL_COUNTS,
                               out_dir = CONFIG$OUT_DIR,
                               count_input_mode = CONFIG$COUNT_INPUT_MODE)
    counts <- counts_info$counts
  }

  if (isTRUE(CONFIG$META_ONLY)) {
    message("[Meta] META_ONLY=TRUE: stopping after metadata validation and input checks.")
    write_run_manifest(CONFIG, CONFIG$OUT_DIR, counts, meta)
    pr$done()
    return(invisible(TRUE))
  }

  if (is.null(vst_all)) {
    pr$step("VST")
    vst_all <- vst_from_counts(counts, meta, counts_for_normalized_baseline = counts_info$counts_expected, out_dir = CONFIG$OUT_DIR)
  } else {
    pr$step("VST")
    message("[VST] Skipping DESeq2::vst because a precomputed VST matrix was supplied.")
  }
  vst_raw <- vst_all$mat
  cd <- vst_all$coldata

  pr$step("Collapse VST matrix to gene symbols")
  vst_mat <- collapse_vst_to_symbols(vst_raw, counts_info$gene_ids, gene_id_type = CONFIG$GENE_ID_TYPE)

  baseline_source_matrix <- "VST_fallback"
  baseline_mat <- vst_mat
  if (!is.null(vst_all$normalized_counts)) {
    norm_counts_mat <- collapse_matrix_to_symbols(
      vst_all$normalized_counts,
      counts_info$gene_ids,
      gene_id_type = CONFIG$GENE_ID_TYPE,
      tie_breaker = rowMeans(vst_all$normalized_counts, na.rm = TRUE)
    )
    common_genes <- intersect(rownames(vst_mat), rownames(norm_counts_mat))
    if (!length(common_genes)) {
      stop("No overlap between VST genes and normalized-count genes after symbol collapse.", call. = FALSE)
    }
    vst_mat <- vst_mat[common_genes, , drop = FALSE]
    norm_counts_mat <- norm_counts_mat[common_genes, , drop = FALSE]
    baseline_mat <- norm_counts_mat
    baseline_source_matrix <- "size_factor_normalized_counts"
  } else {
    message("[Baseline] WARNING: no count matrix was supplied; using VST matrix for program-level baseline ranks.")
  }

  readr::write_csv(
    tibble(gene = rownames(vst_mat)),
    file.path(CONFIG$OUT_DIR, "analysis_gene_symbols.csv")
  )

  pr$step("Distortion-associated PC / NRI / module clustering")
  message("[Step 6] Computing shared 4 h CTRL baseline ranks...")
  baseline_vec <- compute_baseline(baseline_mat, cd)
  write_baseline_source_audit(baseline_vec, baseline_source_matrix, CONFIG$OUT_DIR)
  readr::write_csv(
    tibble(
      gene = names(baseline_vec),
      baseline_value = as.numeric(baseline_vec),
      baseline_band = as.character(band_from_baseline(baseline_vec)),
      baseline_source_matrix = baseline_source_matrix
    ),
    file.path(CONFIG$OUT_DIR, "program_baseline_ranks_and_bands.csv")
  )
  head_genes <- head_genes_from_baseline(baseline_vec, head_frac = CONFIG$HEAD_FRAC)

  reuse_files <- c(
    load = file.path(CONFIG$OUT_DIR, "time_aligned_axis_loadings.csv"),
    nri = file.path(CONFIG$OUT_DIR, "neighborhood_rank_instability.csv"),
    mdi = file.path(CONFIG$OUT_DIR, "module_MDI.csv"),
    assignments = file.path(CONFIG$OUT_DIR, "module_assignments.csv"),
    k_selection = file.path(CONFIG$OUT_DIR, "module_k_selection.csv")
  )
  can_reuse_modules <- isTRUE(CONFIG$REUSE_MODULE_TABLES) && all(file.exists(reuse_files))

  if (can_reuse_modules) {
    message("[Step 6] Reusing existing module/axis tables from output directory.")
    load_tbl <- readr::read_csv(reuse_files[["load"]], show_col_types = FALSE)
    nri_tbl <- readr::read_csv(reuse_files[["nri"]], show_col_types = FALSE)
    mod <- list(
      mdi = readr::read_csv(reuse_files[["mdi"]], show_col_types = FALSE),
      assignments = readr::read_csv(reuse_files[["assignments"]], show_col_types = FALSE) |>
        mutate(
          gene = as.character(gene),
          cluster = as.character(cluster),
          module_membership = suppressWarnings(as.numeric(module_membership)),
          lfc_sur = suppressWarnings(as.numeric(lfc_sur))
        ),
      k_selection = readr::read_csv(reuse_files[["k_selection"]], show_col_types = FALSE),
      selected_k = CONFIG$K_MODULES,
      clustering_engine = "reused_from_output_tables"
    )
  } else {
    if (isTRUE(CONFIG$REUSE_MODULE_TABLES)) {
      message("[Step 6] Reuse requested but one or more module tables are missing; recomputing modules.")
    }
    message("[Step 6] Computing distortion-associated PC loadings...")
    load_tbl <- time_aligned_axis(vst_mat, cd, ncomp = CONFIG$NCOMP_PCA)
    message("[Step 6] Computing neighborhood-rank instability...")
    nri_tbl <- nri_table(vst_mat, cd)
    message("[Step 6] Building seven-module solution with low-memory clustering...")
    mod <- build_modules(
      vst_mat,
      cd,
      load_tbl = load_tbl,
      nri_tbl = nri_tbl,
      baseline_vec = baseline_vec,
      head_frac = CONFIG$HEAD_FRAC,
      k_min = CONFIG$MODULE_K_MIN,
      k_max = CONFIG$MODULE_K_MAX,
      final_k = CONFIG$K_MODULES,
      require_seven = CONFIG$REQUIRE_SEVEN_MODULES,
      seed = CONFIG$MODULE_SEED,
      cluster_engine = CONFIG$MODULE_CLUSTER_ENGINE,
      exact_pam_max_n = CONFIG$MODULE_EXACT_PAM_MAX_N,
      silhouette_sample_n = CONFIG$MODULE_SILHOUETTE_SAMPLE_N,
      clara_samples = CONFIG$MODULE_CLARA_SAMPLES,
      clara_sampsize = CONFIG$MODULE_CLARA_SAMPSIZE
    )
    readr::write_csv(load_tbl, file.path(CONFIG$OUT_DIR, "time_aligned_axis_loadings.csv"))
    readr::write_csv(nri_tbl, file.path(CONFIG$OUT_DIR, "neighborhood_rank_instability.csv"))
    readr::write_csv(mod$mdi, file.path(CONFIG$OUT_DIR, "module_MDI.csv"))
    readr::write_csv(mod$assignments, file.path(CONFIG$OUT_DIR, "module_assignments.csv"))
    readr::write_csv(mod$k_selection, file.path(CONFIG$OUT_DIR, "module_k_selection.csv"))
  }


  pr$step("Sample correlation")
  sample_correlation_heatmap(vst_mat, cd, file.path(CONFIG$OUT_DIR, "Heatmap_SampleCorr_ALL.png"))

  pr$step("MSigDB term2gene")
  t2g_all <- build_msig_term2gene(TARGET_GSETS)
  report_unmatched_sets(TARGET_GSETS, t2g_all, file.path(CONFIG$OUT_DIR, "MSIGDB_unmatched_sets.txt"))
  genesets_all <- t2g_to_genesets(t2g_all)

  pr$step("Module labels")
  module_labels <- label_modules_with_msigdb(
    mod$assignments,
    universe_genes = rownames(vst_mat),
    t2g = t2g_all,
    out_csv = file.path(CONFIG$OUT_DIR, "module_labels.csv")
  )
  mod_assign_labeled <- mod$assignments |>
    left_join(module_labels, by = c("cluster" = "Cluster"))

  pr$step("Figure 7 panels")
  if (isTRUE(CONFIG$RUN_FIG7)) {
    gene_cluster_heatmap(
      vst_mat,
      cd,
      mod_assign_labeled,
      module_labels,
      baseline_vec,
      genes_per_module = CONFIG$GENES_PER_MODULE_HEATMAP,
      out_file = "Heatmap_Modules_Genes.png"
    )

    modules_terms_heatmap(
      mod_assign_labeled,
      universe_genes = rownames(vst_mat),
      t2g = t2g_all,
      out_png = file.path(CONFIG$OUT_DIR, "MSIGDB_modulesXterms_heatmap.png"),
      top_terms_per_module = CONFIG$FIG7B_TOP_TERMS
    )

    t2g_for_heat <- t2g_all |>
      mutate(bl = baseline_vec[gene_symbol]) |>
      filter(!is.na(bl)) |>
      arrange(desc(bl)) |>
      select(-bl)

    gene_set_union_heatmap(
      vst_mat,
      cd,
      t2g_for_heat,
      baseline_vec,
      max_per_set = CONFIG$MAX_GENES_PER_SET_UNION,
      out_file = "Heatmap_GeneSets_Union.png"
    )

    gsea_time_aligned_axis(
      load_tbl,
      out_csv = file.path(CONFIG$OUT_DIR, "GSEA_time_aligned_axis_Hallmark.csv"),
      out_png = file.path(CONFIG$OUT_DIR, "GSEA_time_aligned_axis_Hallmark_barplot.png")
    )

    compose_figure7(CONFIG$OUT_DIR, include_panel_c = CONFIG$FIG7_INCLUDE_PANEL_C)
  } else {
    message("[Fig7] Skipped (RUN_FIG7=FALSE)")
  }

  pr$step("Figure 8 fgsea")
  fg8 <- NULL
  fg8_results_file <- file.path(CONFIG$OUT_DIR, "Fig8_fgsea_results.csv")
  fg8_used_file <- file.path(CONFIG$OUT_DIR, "Fig8_contrasts_used.csv")
  can_reuse_fgsea <- isTRUE(CONFIG$REUSE_FGSEA_TABLES) && file.exists(fg8_results_file) && file.exists(fg8_used_file)
  if (isTRUE(CONFIG$RUN_FIG8)) {
    if (can_reuse_fgsea) {
      message("[Fig8] Reusing existing fgsea result tables from output directory.")
      fg8_res <- readr::read_csv(fg8_results_file, show_col_types = FALSE)
      if (!"sig" %in% names(fg8_res)) fg8_res$sig <- -log10(suppressWarnings(as.numeric(fg8_res$padj)) + 1e-300)
      if (!"LeadingEdgeSize" %in% names(fg8_res)) fg8_res$LeadingEdgeSize <- NA_real_
      if (!"LeadingEdgeHeadFrac" %in% names(fg8_res)) fg8_res$LeadingEdgeHeadFrac <- NA_real_
      fg8_res <- fg8_res |>
        mutate(
          NES = suppressWarnings(as.numeric(NES)),
          padj = suppressWarnings(as.numeric(padj)),
          sig = suppressWarnings(as.numeric(sig)),
          LeadingEdgeSize = suppressWarnings(as.numeric(LeadingEdgeSize)),
          LeadingEdgeHeadFrac = suppressWarnings(as.numeric(LeadingEdgeHeadFrac))
        )
      fg8 <- list(
        results = fg8_res,
        contrasts_used = readr::read_csv(fg8_used_file, show_col_types = FALSE)
      )
    } else {
      if (isTRUE(CONFIG$REUSE_FGSEA_TABLES)) {
        message("[Fig8] Reuse requested but Fig8_fgsea_results.csv/Fig8_contrasts_used.csv are missing; recomputing fgsea.")
      }
      fg8 <- run_fgsea_contrasts(
        vst_mat, cd,
        genesets = genesets_all,
        contrast_defs = FIG8_CONTRASTS,
        nperm = CONFIG$FGSEA_NPERM,
        minSize = CONFIG$FGSEA_MIN_SIZE,
        maxSize = CONFIG$FGSEA_MAX_SIZE,
        head_genes = head_genes
      )
      readr::write_csv(fg8$contrasts_used, fg8_used_file)
      readr::write_csv(flatten_fgsea_for_csv(fg8$results), fg8_results_file)
    }
    if (!is.data.frame(fg8$results) || nrow(fg8$results) == 0) {
      message("[Fig8] WARNING: fgsea produced 0 results. Check sample_meta_counts.csv and sample_meta_qc_mismatches.csv.")
    }
  } else {
    message("[Fig8] Skipped (RUN_FIG8=FALSE)")
  }

  pr$step("Figure 8 panels A/B/C")
  if (isTRUE(CONFIG$RUN_FIG8) &&
      !is.null(fg8) &&
      is.data.frame(fg8$results) &&
      nrow(fg8$results) > 0 &&
      all(c("Time", "Stress", "Contrast") %in% names(fg8$results))) {
    dfp <- prep_panel_df(fg8$results)

    df_myc <- dfp |> filter(Contrast == "DT_vs_D")
    plot_dotmap(
      df_myc,
      intersect(FIG8A_MYC_SETS, unique(df_myc$pathway)),
      title = "MYC/growth programs (DT vs D)",
      out_png = file.path(CONFIG$OUT_DIR, "Fig8A_MYC_targets_dotmap.png"),
      facet_by_contrast = FALSE,
      color_name = "NES",
      size_name = "-log10(FDR)"
    )

    df_ddr <- dfp |> filter(Contrast == "DT_vs_D")
    plot_dotmap(
      df_ddr,
      intersect(FIG8B_DDR_SETS, unique(df_ddr$pathway)),
      title = "DDR/checkpoint programs (DT vs D)",
      out_png = file.path(CONFIG$OUT_DIR, "Fig8B_DDR_dotmap.png"),
      facet_by_contrast = FALSE,
      color_name = "NES",
      size_name = "-log10(FDR)"
    )

    plot_es_small_multiples(
      vst_mat, cd, genesets_all,
      contrast_defs = FIG8_CONTRASTS,
      contrast_name = "DT_vs_D",
      stress = "NoCPT",
      sets = intersect(FIG8C_ES_SETS, names(genesets_all)),
      out_png = "Fig8C_ES_curves.png"
    )
  } else if (isTRUE(CONFIG$RUN_FIG8)) {
    message("[Fig8] Skipping panels A/B/C (no fgsea results).")
  }

  pr$step("Figure 8 panel D + leading-edge control")
  if (isTRUE(CONFIG$RUN_FIG8) &&
      !is.null(fg8) &&
      is.data.frame(fg8$results) &&
      nrow(fg8$results) > 0 &&
      all(c("Time", "Stress", "Contrast") %in% names(fg8$results)) &&
      isTRUE(CONFIG$RUN_MITO_CONTROL)) {
    dfp <- prep_panel_df(fg8$results)
    mito_sets_present <- intersect(FIG8D_MITO_SETS, unique(dfp$pathway))
    df_mito <- dfp |> filter(Contrast %in% c("D_vs_Ctrl", "DT_vs_D"))
    readr::write_csv(
      make_mito_contrast_layout(df_mito) |>
        distinct(ComparisonBlock, ColumnLabel, Contrast, Time, Stress, ExactContrast) |>
        arrange(ComparisonBlock, ColumnLabel, ExactContrast),
      file.path(CONFIG$OUT_DIR, "Fig8D_mito_observed_strata.csv")
    )
    plot_mito_matched_contrast_dotmap(
      df_mito,
      mito_sets_present,
      title = "Mitochondrial control: enrichment/significance",
      out_png = file.path(CONFIG$OUT_DIR, "Fig8D_MitoControl_dotmap.png"),
      color_name = "NES",
      size_name = "-log10(FDR)"
    )
    if ("LeadingEdgeHeadFrac" %in% names(df_mito)) {
      plot_mito_matched_contrast_dotmap(
        df_mito,
        mito_sets_present,
        title = "Mitochondrial leading-edge localization across matched contrasts",
        out_png = file.path(CONFIG$OUT_DIR, "SuppFig3_MitoControl_LeadingEdgeHeadFrac.png"),
        color_var = "LeadingEdgeHeadFrac",
        size_var = "LeadingEdgeSize",
        color_name = "Fraction of leading edge in expression head",
        size_name = "Leading-edge genes, n"
      )
    }
  } else if (isTRUE(CONFIG$RUN_FIG8) && isTRUE(CONFIG$RUN_MITO_CONTROL)) {
    message("[Fig8] Skipping mito control panel (no fgsea results).")
  }

  pr$step("Figure 8 composite")
  if (isTRUE(CONFIG$RUN_FIG8)) {
    compose_figure8(CONFIG$OUT_DIR)
    compose_supplementary_figure3(CONFIG$OUT_DIR)
  }

  pr$step("Write run manifest")
  write_run_manifest(CONFIG, CONFIG$OUT_DIR, counts, meta)
  writeLines("ok", con = file.path(CONFIG$OUT_DIR, "figure7_8.done"))

  pr$done()
  invisible(TRUE)
}

run_pipeline()
