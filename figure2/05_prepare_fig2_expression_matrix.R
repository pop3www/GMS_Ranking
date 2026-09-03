#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table) })
SCRIPT_VERSION <- "2026-06-12 v39 (prepare size-factor-normalized RSEM expected counts for Figure 2)"
message("[05] script version: ", SCRIPT_VERSION)
args <- commandArgs(trailingOnly=TRUE)
parse_args <- function(argv){ out<-list(); i<-1L; while(i<=length(argv)){ a<-argv[[i]]; if(startsWith(a,"--")){ a2<-sub("^--","",a); if(grepl("=",a2,fixed=TRUE)){sp<-strsplit(a2,"=",fixed=TRUE)[[1]]; key<-sp[[1]]; val<-paste(sp[-1],collapse="=")} else {key<-a2; if(i<length(argv)&&!startsWith(argv[[i+1L]],"--")){val<-argv[[i+1L]]; i<-i+1L}else val<-"TRUE"}; key<-tolower(gsub("_","-",key)); out[[key]]<-val}; i<-i+1L}; out }
opt <- parse_args(args)
get_arg <- function(keys, default=NA_character_, required=FALSE){ keys<-tolower(gsub("_","-",keys)); for(k in keys) if(!is.null(opt[[k]])&&nzchar(opt[[k]])) return(opt[[k]]); if(required) stop("[05] Missing required arg: --",keys[[1]],call.=FALSE); default }
get_num <- function(keys, default){ v<-suppressWarnings(as.numeric(get_arg(keys, as.character(default)))); if(!is.finite(v)) default else v }
counts_path <- get_arg(c("counts","raw-counts","raw"), required=TRUE)
id_col <- get_arg(c("id-col","id"), default="gene_id")
out_dir <- get_arg("out", required=TRUE)
min_total <- get_num(c("min-total","min-count-total"), 80)
if(!file.exists(counts_path)) stop("[05] counts file not found: ",counts_path,call.=FALSE)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dt <- fread(counts_path, header=TRUE, check.names=FALSE, data.table=TRUE, quote="\"")
if(!(id_col %in% names(dt))){ message("[05] id-col '",id_col,"' not found; using first column: ",names(dt)[1]); id_col <- names(dt)[1] }
setnames(dt, id_col, "gene_id"); dt[, gene_id := as.character(gene_id)]
sample_cols <- setdiff(names(dt), "gene_id")
dt[, (sample_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols=sample_cols]
mat <- as.matrix(dt[, ..sample_cols]); storage.mode(mat) <- "double"; mat[!is.finite(mat)] <- NA_real_
if(any(mat < 0, na.rm=TRUE)) stop("[05] negative count values detected", call.=FALSE)
row_total <- rowSums(mat, na.rm=TRUE); keep <- is.finite(row_total) & row_total >= min_total
message("[05] genes before filter: ",nrow(dt),"; after total expected-count >= ",min_total,": ",sum(keep))
dt_f <- dt[keep]; mat_f <- mat[keep,,drop=FALSE]
geom_pos <- function(x){ x<-x[is.finite(x)&x>0]; if(!length(x)) return(NA_real_); exp(mean(log(x))) }
gm <- apply(mat_f, 1, geom_pos); valid <- is.finite(gm) & gm > 0
sf <- rep(NA_real_, ncol(mat_f))
for(j in seq_len(ncol(mat_f))){ ratios <- mat_f[valid,j] / gm[valid]; ratios <- ratios[is.finite(ratios)&ratios>0]; sf[j] <- if(length(ratios)) median(ratios,na.rm=TRUE) else NA_real_ }
if(any(!is.finite(sf)|sf<=0)){ lib <- colSums(mat_f,na.rm=TRUE); sf2 <- lib/exp(mean(log(lib[lib>0]))); sf[!is.finite(sf)|sf<=0] <- sf2[!is.finite(sf)|sf<=0] }
sf <- sf / exp(mean(log(sf[sf>0])))
norm <- sweep(mat_f, 2, sf, FUN="/"); colnames(norm) <- sample_cols
out_dt <- cbind(data.table(gene_id=dt_f$gene_id), as.data.table(norm))
out_expr <- file.path(out_dir,"Fig2_expression_normcounts.tsv")
out_sf <- file.path(out_dir,"Fig2_expression_size_factors.tsv")
out_manifest <- file.path(out_dir,"Fig2_expression_normcounts_manifest.txt")
fwrite(out_dt, out_expr, sep="\t")
fwrite(data.table(sample=sample_cols, size_factor=sf), out_sf, sep="\t")
cat("Figure 2 normalized expression matrix\nScript version: ",SCRIPT_VERSION,"\nInput counts: ",counts_path,"\nNormalization: median-ratio size factors (positive-count geometric means)\nFilter: total expected counts >= ",min_total,"\nGenes before filter: ",nrow(dt),"\nGenes after filter: ",nrow(out_dt),"\nOutput matrix: ",out_expr,"\nSize factors: ",out_sf,"\n",sep="",file=out_manifest)
message("[05] wrote: ",out_expr); message("[05] wrote: ",out_sf)
