#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(grid)
})

SCRIPT_VERSION <- "2026-06-18 v42 (publication compositor; aspect-preserving panels; larger final canvas; PNG+PDF+SVG)"
message("[assemble] script version: ", SCRIPT_VERSION)

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
opt <- parse_args(args)
get_arg <- function(keys, default = NA_character_, required = FALSE) {
  keys <- tolower(gsub("_", "-", keys))
  for (k in keys) if (!is.null(opt[[k]]) && nzchar(opt[[k]])) return(opt[[k]])
  if (required) stop("[assemble] Missing required arg: --", keys[[1]], call. = FALSE)
  default
}
get_num <- function(keys, default) {
  raw <- get_arg(keys, default = as.character(default))
  val <- suppressWarnings(as.numeric(raw))
  if (!is.finite(val)) default else val
}
get_logical <- function(keys, default = TRUE) {
  raw <- tolower(get_arg(keys, default = if (default) "true" else "false"))
  raw %in% c("true", "t", "1", "yes", "y")
}

figure_dir <- get_arg(c("figure-dir", "figdir", "dir"), default = ".", required = TRUE)
width_in  <- get_num(c("width", "width-in"), 18.0)
height_in <- get_num(c("height", "height-in"), 16.5)
dpi       <- get_num("dpi", 300)
label_size <- get_num("label-size", 24)
preserve_aspect <- get_logical("preserve-aspect", TRUE)

if (!requireNamespace("png", quietly = TRUE)) {
  stop("[assemble] Package 'png' is required. Please install.packages('png').", call. = FALSE)
}
if (!dir.exists(figure_dir)) stop("[assemble] figure-dir not found: ", figure_dir, call. = FALSE)

panel_files <- c(
  A = file.path(figure_dir, "Fig2A_delta_density_controls.png"),
  B = file.path(figure_dir, "Fig2B_priming_dependent_component.png"),
  C = file.path(figure_dir, "Fig2C_DTvsD_vs_TamvsCtrl_scatter.png"),
  D = file.path(figure_dir, "Fig2D_tetO_logo_weights.png"),
  E = file.path(figure_dir, "Fig2E_hit_density_by_band.png"),
  F = file.path(figure_dir, "Fig2F_enrichment_OR.png")
)
missing <- panel_files[!file.exists(panel_files)]
if (length(missing)) stop("[assemble] missing panel PNG(s): ", paste(missing, collapse = ", "), call. = FALSE)

imgs <- lapply(panel_files, png::readPNG, native = FALSE)

# Four-row manuscript layout with less distortion and more breathing room.
margin_x <- 0.025
margin_y <- 0.016
gap_x <- 0.026
gap_y <- 0.018

hA  <- 0.245
hBC <- 0.255
hD  <- 0.220
hEF <- 0.200

yEF <- margin_y
yD  <- yEF + hEF + gap_y
yBC <- yD + hD + gap_y
yA  <- yBC + hBC + gap_y

xA <- margin_x; wA <- 1 - 2 * margin_x

xB <- margin_x
wB <- 0.620
xC <- xB + wB + gap_x
wC <- 1 - margin_x - xC

xD <- margin_x; wD <- 1 - 2 * margin_x

wE <- (1 - 2 * margin_x - gap_x) / 2
wF <- wE
xE <- margin_x
xF <- xE + wE + gap_x

placements <- list(
  A = c(x = xA, y = yA,  w = wA, h = hA),
  B = c(x = xB, y = yBC, w = wB, h = hBC),
  C = c(x = xC, y = yBC, w = wC, h = hBC),
  D = c(x = xD, y = yD,  w = wD, h = hD),
  E = c(x = xE, y = yEF, w = wE, h = hEF),
  F = c(x = xF, y = yEF, w = wF, h = hEF)
)

fit_rect <- function(img, pl) {
  if (!preserve_aspect) return(pl)
  ih <- dim(img)[1]
  iw <- dim(img)[2]
  img_aspect <- ih / iw
  slot_aspect <- (pl[["h"]] * height_in) / (pl[["w"]] * width_in)
  if (img_aspect > slot_aspect) {
    # Image is relatively taller than slot; use full slot height.
    h <- pl[["h"]]
    w <- h * height_in / (img_aspect * width_in)
  } else {
    # Image is relatively wider than slot; use full slot width.
    w <- pl[["w"]]
    h <- w * width_in * img_aspect / height_in
  }
  c(x = pl[["x"]] + (pl[["w"]] - w) / 2,
    y = pl[["y"]] + (pl[["h"]] - h) / 2,
    w = w,
    h = h)
}

draw_panel <- function(img, pl) {
  rr <- fit_rect(img, pl)
  grid.raster(img,
    x = unit(rr[["x"]], "npc"),
    y = unit(rr[["y"]], "npc"),
    width = unit(rr[["w"]], "npc"),
    height = unit(rr[["h"]], "npc"),
    just = c("left", "bottom"), interpolate = TRUE)
}

draw_label <- function(label, pl) {
  grid.text(label,
    x = unit(max(0.004, pl[["x"]] - 0.016), "npc"),
    y = unit(min(0.995, pl[["y"]] + pl[["h"]] + 0.002), "npc"),
    just = c("left", "bottom"),
    gp = gpar(fontface = "bold", fontsize = label_size))
}

draw_all <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  for (nm in names(imgs)) draw_panel(imgs[[nm]], placements[[nm]])
  for (nm in names(imgs)) draw_label(nm, placements[[nm]])
}

out_png <- file.path(figure_dir, "Figure2_combined.png")
out_pdf <- file.path(figure_dir, "Figure2_combined.pdf")
out_svg <- file.path(figure_dir, "Figure2_combined.svg")

png(filename = out_png, width = width_in, height = height_in, units = "in", res = dpi, bg = "white")
draw_all(); dev.off()

pdf(file = out_pdf, width = width_in, height = height_in, bg = "white", useDingbats = FALSE)
draw_all(); dev.off()

svg(filename = out_svg, width = width_in, height = height_in, bg = "white")
draw_all(); dev.off()

message("[assemble] wrote: ", out_png)
message("[assemble] wrote: ", out_pdf)
message("[assemble] wrote: ", out_svg)
