#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

SCRIPT_VERSION <- "2026-06-12 v39 (integrated Fig2 contract + normalized-count Figure 2 inputs; All/Up/Down densities)"
message("[07] script version: ", SCRIPT_VERSION)

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Usage: 07_summarize_fimo.R \\\n  --bands <rank_bands.csv> \\\n  --out <output_dir> \\\n  [--id-col gene] \\\n  [--fimo-ebox <fimo_dir_or_tsv>] [--fimo-teto <fimo_dir_or_tsv>] \\\n  [--pval-max 1e-4] [--qval-max 1] [--tss-window 1000]\n\n",
    "Env fallbacks:\n",
    "  FIMO_EBOX_DIR, FIMO_TETO_DIR, PVAL_MAX, QVAL_MAX, TSS_WINDOW\n",
    sep = ""
  )
}

# Minimal robust arg parser: supports --key value and --key=value.
parse_args <- function(argv) {
  out <- list()
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (startsWith(a, "--")) {
      a2 <- sub("^--", "", a)
      if (grepl("=", a2, fixed = TRUE)) {
        sp <- strsplit(a2, "=", fixed = TRUE)[[1]]
        key <- sp[[1]]
        val <- paste(sp[-1], collapse = "=")
      } else {
        key <- a2
        if (i < length(argv) && !startsWith(argv[[i + 1L]], "--")) {
          val <- argv[[i + 1L]]
          i <- i + 1L
        } else {
          val <- "TRUE"
        }
      }
      key <- tolower(gsub("_", "-", key))
      out[[key]] <- val
    }
    i <- i + 1L
  }
  out
}

amap <- parse_args(args)
get_arg <- function(keys, default = NA_character_, required = FALSE, env = NULL) {
  keys <- tolower(gsub("_", "-", keys))
  for (k in keys) {
    if (!is.null(amap[[k]]) && length(amap[[k]]) == 1L && nzchar(amap[[k]])) return(amap[[k]])
  }
  if (!is.null(env)) {
    ev <- Sys.getenv(env, unset = "")
    if (nzchar(ev)) return(ev)
  }
  if (required) stop("[07] Missing required arg: --", keys[[1]], call. = FALSE)
  default
}
get_num <- function(keys, default = NA_real_, required = FALSE, env = NULL) {
  raw <- get_arg(keys, default = NA_character_, required = required, env = env)
  if (is.na(raw) || !nzchar(raw)) return(default)
  val <- suppressWarnings(as.numeric(raw))
  if (!is.finite(val)) return(default)
  val
}

if (length(args) == 0) {
  usage()
  quit(save = "no", status = 1)
}

bands_file <- get_arg(c("bands", "rank-bands", "rank_bands"), required = TRUE)
out_dir <- get_arg(c("out", "outdir", "out-dir"), required = TRUE)
id_col <- get_arg(c("id-col", "id_col"), default = "gene")

fimo_ebox <- get_arg(c("fimo-ebox", "fimo_ebox"), default = NA_character_, env = "FIMO_EBOX_DIR")
fimo_teto <- get_arg(c("fimo-teto", "fimo_teto"), default = NA_character_, env = "FIMO_TETO_DIR")

pval_max <- get_num(c("pval-max", "pval_max"), default = 1e-4, env = "PVAL_MAX")
qval_max <- get_num(c("qval-max", "qval_max"), default = 1, env = "QVAL_MAX")
tss_window <- get_num(c("tss-window", "tss_window"), default = 1000, env = "TSS_WINDOW")

if (!file.exists(bands_file)) stop("[07] bands file not found: ", bands_file, call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

has_ebox <- !(is.na(fimo_ebox) || !nzchar(fimo_ebox))
has_teto <- !(is.na(fimo_teto) || !nzchar(fimo_teto))
if (!has_ebox && !has_teto) stop("[07] Need at least one FIMO input", call. = FALSE)

# ---- gene helpers ----
extract_ensg <- function(x) {
  x <- as.character(x)
  has <- !is.na(x) & grepl("ENSG[0-9]+", x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  out[has] <- sub("^.*?(ENSG[0-9]+).*$", "\\1", x[has], perl = TRUE)
  out
}

normalize_gene_label <- function(x) {
  x <- as.character(x)
  x <- gsub("^>", "", x)
  x <- trimws(x)
  x <- sub("\\([+-]\\)$", "", x)                # strand suffix
  x <- sub("^(ENSG[0-9]+)\\.[0-9]+", "\\1", x) # ENSG version
  x
}

std_band <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[!(x %in% c("head", "mid", "tail"))] <- NA_character_
  x
}

std_ud <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  y <- toupper(x)
  out <- ifelse(y %in% c("UP", "U"), "Up",
         ifelse(y %in% c("DOWN", "D"), "Down",
         ifelse(y %in% c("NEUTRAL", "UNCHANGED", "NS"), "Neutral", x)))
  out[is.na(out) | !nzchar(out)] <- "Neutral"
  out
}

# ---- rank bands ----
bands <- fread(bands_file)
if (!(id_col %in% names(bands))) {
  alt <- intersect(c("gene", "gene_id", "id", "ensg"), names(bands))
  if (length(alt)) {
    message("[07] WARNING: id-col '", id_col, "' not found; using '", alt[[1]], "'")
    id_col <- alt[[1]]
  } else {
    stop("[07] id-col not found in bands file: ", id_col, call. = FALSE)
  }
}
if (!"band" %in% names(bands)) stop("[07] bands file missing 'band' column", call. = FALSE)

bands[, gene_label := normalize_gene_label(get(id_col))]
bands[, gene_core := extract_ensg(gene_label)]
bands[, band := std_band(band)]

ud_col <- NA_character_
for (cc in c("de_class", "UpDown", "updown", "up_down", "regulation", "call", "direction", "de_call")) {
  if (cc %in% names(bands)) { ud_col <- cc; break }
}
if (!is.na(ud_col)) {
  bands[, updown_group := std_ud(get(ud_col))]
} else {
  bands[, updown_group := "Neutral"]
}

bands_u <- unique(bands[!is.na(gene_core) & !is.na(band), .(gene_core, gene = gene_label, band, updown_group)])
if (nrow(bands_u) == 0) stop("[07] No valid gene_core/band rows in rank_bands", call. = FALSE)

# denominator tables for All/Up/Down, by band
prom_all <- bands_u[, .(promoters = uniqueN(gene_core)), by = .(band)]
prom_all[, updown_group := "All"]
prom_ud <- bands_u[updown_group %in% c("Up", "Down"), .(promoters = uniqueN(gene_core)), by = .(band, updown_group)]
prom_den <- rbindlist(list(prom_all, prom_ud), fill = TRUE)

# ---- FIMO readers ----
resolve_fimo <- function(path) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  if (dir.exists(path)) {
    cand <- file.path(path, c("fimo.tsv", "fimo.txt", "fimo.tsv.gz", "fimo.txt.gz"))
    cand <- cand[file.exists(cand)]
    if (!length(cand)) stop("[07] No fimo.tsv/fimo.txt found in directory: ", path, call. = FALSE)
    return(cand[[1]])
  }
  if (!file.exists(path)) stop("[07] FIMO path not found: ", path, call. = FALSE)
  path
}

read_fimo_any <- function(path) {
  f <- resolve_fimo(path)
  dt <- suppressWarnings(fread(f))
  # Drop MEME comment footer if fread somehow preserved it; robust to empty.
  if (!nrow(dt)) return(dt)
  if (!"sequence_name" %in% names(dt)) {
    alt <- intersect(c("sequence_name", "sequence", "seqname"), names(dt))
    if (!length(alt)) stop("[07] FIMO file missing sequence_name column: ", f, call. = FALSE)
    setnames(dt, alt[[1]], "sequence_name")
  }
  pcol <- intersect(c("p-value", "p_value", "p.value", "pvalue", "pval"), names(dt))
  qcol <- intersect(c("q-value", "q_value", "q.value", "qvalue", "qval"), names(dt))
  if (!length(pcol)) stop("[07] FIMO file missing p-value column: ", f, call. = FALSE)
  setnames(dt, pcol[[1]], "pval")
  if (length(qcol)) setnames(dt, qcol[[1]], "qval") else dt[, qval := NA_real_]
  if (!"start" %in% names(dt) || !"stop" %in% names(dt)) stop("[07] FIMO file missing start/stop columns: ", f, call. = FALSE)
  dt
}

make_hits <- function(dt, type_label) {
  if (is.null(dt) || !nrow(dt)) return(data.table())
  dt <- copy(dt)
  dt[, gene_raw := as.character(sequence_name)]
  dt[, gene := normalize_gene_label(gene_raw)]
  dt[, gene_core := extract_ensg(gene)]
  dt <- dt[!is.na(gene_core)]
  if (!nrow(dt)) return(data.table())
  dt[, pval := suppressWarnings(as.numeric(pval))]
  dt[, qval := suppressWarnings(as.numeric(qval))]
  dt <- dt[!is.na(pval) & pval <= pval_max]
  dt <- dt[is.na(qval) | qval <= qval_max]
  if (!nrow(dt)) return(data.table())
  dt[, start_num := suppressWarnings(as.numeric(start))]
  dt[, stop_num := suppressWarnings(as.numeric(stop))]
  dt[, start_center := (start_num + stop_num) / 2]
  if (is.finite(tss_window)) {
    tss_pos <- floor(max(dt$stop_num, na.rm = TRUE) / 2 + 1)
    dt[, dist_to_tss := start_center - tss_pos]
    dt <- dt[is.na(dist_to_tss) | abs(dist_to_tss) <= tss_window]
  } else {
    dt[, dist_to_tss := NA_real_]
  }
  if (!nrow(dt)) return(data.table())
  dt[, type := type_label]
  keep <- intersect(c("gene_core", "gene", "gene_raw", "type", "start", "stop", "strand", "score", "pval", "qval", "dist_to_tss"), names(dt))
  dt[, ..keep]
}

hits_list <- list()
if (has_ebox) hits_list[[length(hits_list) + 1L]] <- make_hits(read_fimo_any(fimo_ebox), "EBOX")
if (has_teto) hits_list[[length(hits_list) + 1L]] <- make_hits(read_fimo_any(fimo_teto), "TETO")
hits <- rbindlist(hits_list, fill = TRUE)
if (!nrow(hits)) stop("[07] No FIMO hits after filtering; check thresholds and inputs", call. = FALSE)

# Join to bands; keep unmatched rows in full hit table but summarize only matched.
hits <- merge(hits, bands_u[, .(gene_core, band, updown_group, gene_rank = gene)], by = "gene_core", all.x = TRUE)
hits[!is.na(gene_rank) & nzchar(gene_rank), gene := gene_rank]
hits[, gene_rank := NULL]

ng_bands <- uniqueN(bands_u$gene_core)
ng_hits <- uniqueN(hits$gene_core)
ng_overlap <- uniqueN(intersect(bands_u$gene_core, hits$gene_core))
message(sprintf("[07] gene_core overlap: %d (bands=%d, hits=%d)", ng_overlap, ng_bands, ng_hits))
n_miss <- sum(is.na(hits$band))
if (n_miss > 0) message(sprintf("[07] WARNING: %d hit-rows did not match rank_bands (band=NA).", n_miss))

fwrite(hits, file.path(out_dir, "fimo_hits_with_bands.csv"))
hits_in <- hits[!is.na(band) & band %in% c("head", "mid", "tail")]

# Helper for complete grid and row-order.
types <- sort(unique(hits$type))
if (!length(types)) types <- c("EBOX", "TETO")
bands_levels <- c("head", "mid", "tail")
groups_levels <- c("All", "Up", "Down")

# ---- hit density: hits per 1,000 promoters ----
make_density <- function() {
  out <- list()
  # All rows
  h_all <- hits_in[, .(hits = .N), by = .(type, band)]
  h_all[, updown_group := "All"]
  h_all <- merge(h_all, prom_all, by = c("band", "updown_group"), all.x = TRUE)
  out[[1]] <- h_all
  # Up/Down rows
  h_ud <- hits_in[updown_group %in% c("Up", "Down"), .(hits = .N), by = .(type, band, updown_group)]
  h_ud <- merge(h_ud, prom_ud, by = c("band", "updown_group"), all.x = TRUE)
  out[[2]] <- h_ud
  dt <- rbindlist(out, fill = TRUE)
  grid <- CJ(type = types, band = bands_levels, updown_group = groups_levels, unique = TRUE)
  dt <- merge(grid, dt, by = c("type", "band", "updown_group"), all.x = TRUE)
  dt[is.na(hits), hits := 0L]
  # denominators for groups without Up/Down genes remain NA/0; keep y=NA
  dt[is.na(promoters), promoters := 0L]
  dt[, y := ifelse(promoters > 0, hits / promoters * 1000, NA_real_)]
  dt[]
}

density <- make_density()
fwrite(density, file.path(out_dir, "hit_density_by_band.csv"))

# ---- gene fraction: genes with >=1 hit / genes in stratum ----
make_fraction <- function() {
  pres <- unique(hits_in[, .(gene_core, type, band, updown_group)])
  # All rows: presence irrespective of updown
  p_all <- unique(pres[, .(gene_core, type, band)])[, .(genes_hit = uniqueN(gene_core)), by = .(type, band)]
  p_all[, updown_group := "All"]
  p_all <- merge(p_all, prom_all, by = c("band", "updown_group"), all.x = TRUE)
  setnames(p_all, "promoters", "genes_total")
  # Up/Down rows
  p_ud <- pres[updown_group %in% c("Up", "Down"), .(genes_hit = uniqueN(gene_core)), by = .(type, band, updown_group)]
  p_ud <- merge(p_ud, prom_ud, by = c("band", "updown_group"), all.x = TRUE)
  setnames(p_ud, "promoters", "genes_total")
  dt <- rbindlist(list(p_all, p_ud), fill = TRUE)
  grid <- CJ(type = types, band = bands_levels, updown_group = groups_levels, unique = TRUE)
  dt <- merge(grid, dt, by = c("type", "band", "updown_group"), all.x = TRUE)
  dt[is.na(genes_hit), genes_hit := 0L]
  dt[is.na(genes_total), genes_total := 0L]
  dt[, y := ifelse(genes_total > 0, genes_hit / genes_total, NA_real_)]
  dt[]
}

frac <- make_fraction()
fwrite(frac, file.path(out_dir, "hit_gene_fraction_by_band.csv"))

# ---- OR band vs rest, gene-level presence ----
calc_or <- function(a, b, c, d) {
  # a: in-band hit, b: in-band no hit, c: out-band hit, d: out-band no hit
  if (any(c(a, b, c, d) == 0)) {
    a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5
  }
  or <- (a / b) / (c / d)
  se <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
  list(or = or, lo = exp(log(or) - 1.96 * se), hi = exp(log(or) + 1.96 * se))
}

all_genes <- unique(bands_u$gene_core)
presence <- unique(hits_in[, .(gene_core, type)])
or_list <- list()
for (tt in types) {
  hit_genes <- unique(presence[type == tt, gene_core])
  for (bb in bands_levels) {
    in_band <- unique(bands_u[band == bb, gene_core])
    out_band <- setdiff(all_genes, in_band)
    a <- length(intersect(in_band, hit_genes))
    b <- length(in_band) - a
    c <- length(intersect(out_band, hit_genes))
    d <- length(out_band) - c
    est <- calc_or(a, b, c, d)
    or_list[[length(or_list) + 1L]] <- data.table(type = tt, band = bb, a = a, b = b, c = c, d = d, odds_ratio = est$or, conf_low = est$lo, conf_high = est$hi, y = est$or)
  }
}
or_dt <- rbindlist(or_list, fill = TRUE)
fwrite(or_dt, file.path(out_dir, "enrichment_odds_ratios.csv"))
# Alias for plotting/export compatibility.
fwrite(or_dt, file.path(out_dir, "enrichment_odds_ratios_by_band.csv"))

# ---- Up vs Down OR within band ----
or_ud_list <- list()
for (tt in types) {
  hit_genes <- unique(presence[type == tt, gene_core])
  for (bb in bands_levels) {
    up_genes <- unique(bands_u[band == bb & updown_group == "Up", gene_core])
    dn_genes <- unique(bands_u[band == bb & updown_group == "Down", gene_core])
    if (!length(up_genes) || !length(dn_genes)) {
      or_ud_list[[length(or_ud_list) + 1L]] <- data.table(type = tt, band = bb, a = length(intersect(up_genes, hit_genes)), b = length(up_genes) - length(intersect(up_genes, hit_genes)), c = length(intersect(dn_genes, hit_genes)), d = length(dn_genes) - length(intersect(dn_genes, hit_genes)), odds_ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_, y = NA_real_)
      next
    }
    a <- length(intersect(up_genes, hit_genes))
    b <- length(up_genes) - a
    c <- length(intersect(dn_genes, hit_genes))
    d <- length(dn_genes) - c
    est <- calc_or(a, b, c, d)
    pv <- suppressWarnings(fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value)
    or_ud_list[[length(or_ud_list) + 1L]] <- data.table(type = tt, band = bb, a = a, b = b, c = c, d = d, odds_ratio = est$or, conf_low = est$lo, conf_high = est$hi, p_value = pv, y = est$or)
  }
}
or_ud <- rbindlist(or_ud_list, fill = TRUE)
fwrite(or_ud, file.path(out_dir, "enrichment_odds_ratios_up_vs_down.csv"))

message("[07] Wrote: ", paste(file.path(out_dir, c("fimo_hits_with_bands.csv", "hit_density_by_band.csv", "hit_gene_fraction_by_band.csv", "enrichment_odds_ratios.csv", "enrichment_odds_ratios_up_vs_down.csv")), collapse = ", "))
