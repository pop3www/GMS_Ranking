#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

SCRIPT_VERSION <- "2026-06-18 v42 (Fig2A enlarged-font normalized-count controls; no head inset; shared y-axis; PNG+PDF+SVG)"
message("[12] script version: ", SCRIPT_VERSION)

args <- commandArgs(trailingOnly=TRUE)
parse_args <- function(argv){
  out<-list(); i<-1L
  while(i<=length(argv)){
    a<-argv[[i]]
    if(startsWith(a,"--")){
      a2<-sub("^--","",a)
      if(grepl("=",a2,fixed=TRUE)){
        sp<-strsplit(a2,"=",fixed=TRUE)[[1]]; key<-sp[[1]]; val<-paste(sp[-1],collapse="=")
      } else {
        key<-a2
        if(i<length(argv)&&!startsWith(argv[[i+1L]],"--")){ val<-argv[[i+1L]]; i<-i+1L } else val<-"TRUE"
      }
      key<-tolower(gsub("_","-",key)); out[[key]]<-val
    }
    i<-i+1L
  }
  out
}
opt <- parse_args(args)
get_arg <- function(keys, default=NA_character_, required=FALSE){
  keys<-tolower(gsub("_","-",keys))
  for(k in keys) if(!is.null(opt[[k]])&&nzchar(opt[[k]])) return(opt[[k]])
  if(required) stop("[12] Missing required arg: --",keys[[1]],call.=FALSE)
  default
}
get_num <- function(keys, default){ v<-suppressWarnings(as.numeric(get_arg(keys, as.character(default)))); if(!is.finite(v)) default else v }
get_int <- function(keys, default) as.integer(round(get_num(keys, default)))

expr_path <- get_arg(c("expr","tpm"), required=TRUE)
id_col <- get_arg("id-col", default="gene_id")
out_dir <- get_arg("out", required=TRUE)
B <- get_int(c("bootstrap","B"), default=as.integer(Sys.getenv("FIG2A_BOOTSTRAP", unset="5000")))
seed <- get_int("seed", 1)
n_grid <- get_int("n-grid", 512)
fixed_ymax <- get_num("fixed-ymax", NA_real_)
min_ymax <- get_num("min-ymax", .002)
y_pad <- get_num("y-pad", 1.08)
width_in <- get_num("width", 16.2)
height_in <- get_num("height", 4.05)

if(!file.exists(expr_path)) stop("[12] expression file not found: ",expr_path,call.=FALSE)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
expr <- fread(expr_path, sep="\t", header=TRUE, check.names=FALSE, quote="\"")
if(!(id_col %in% names(expr))) stop("[12] id-col not found: ",id_col,call.=FALSE)
setnames(expr,id_col,"gene")
num_cols <- setdiff(names(expr),"gene")
expr[, (num_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols=num_cols]

match_cols <- function(rx){
  cc <- num_cols[grepl(rx, num_cols, perl=TRUE)]
  message("[12] regex ",rx," matched ",length(cc),": ",paste(cc,collapse=", "))
  if(!length(cc)) stop("[12] no columns matched regex: ",rx,call.=FALSE)
  cc
}
mean_expr <- function(cols){ m<-as.matrix(expr[, ..cols]); storage.mode(m)<-"double"; rowMeans(m, na.rm=TRUE) }
dens_on_grid <- function(x, grid){
  x<-x[is.finite(x)]
  if(length(unique(x))<2L) return(rep(0,length(grid)))
  d<-density(x, from=min(grid), to=max(grid), n=length(grid), na.rm=TRUE)
  approx(d$x,d$y,xout=grid,rule=2)$y
}

contrasts <- list(
  list(time="4 h", contrast="Dox-only", group_rx='(^|_)4D_RP[0-9]+($|_)', ref_rx='(^|_)4_ctrl_RP[0-9]+($|_)'),
  list(time="24 h", contrast="Dox-only", group_rx='(^|_)24D_RP[0-9]+($|_)', ref_rx='(^|_)24_ctrl_RP[0-9]+($|_)'),
  list(time="4 h", contrast="Tam-only", group_rx='(^|_)4_Tam_RP[0-9]+($|_)', ref_rx='(^|_)4_ctrl_RP[0-9]+($|_)'),
  list(time="24 h", contrast="Tam-only", group_rx='(^|_)24_Tam_RP[0-9]+($|_)', ref_rx='(^|_)24_ctrl_RP[0-9]+($|_)')
)

vecs<-list(); all_x<-c()
for(co in contrasts){
  grp<-mean_expr(match_cols(co$group_rx)); ref<-mean_expr(match_cols(co$ref_rx))
  xg<-log2(grp+1); xr<-log2(ref+1)
  nm<-paste(co$time,co$contrast,sep="__")
  vecs[[nm]]<-list(co=co,group=xg,ref=xr)
  all_x<-c(all_x,xg,xr)
}
all_x<-all_x[is.finite(all_x)]
grid<-seq(max(0,floor(min(all_x))), ceiling(max(all_x)), length.out=n_grid)

set.seed(seed)
res<-list(); boot_res<-list()
for(nm in names(vecs)){
  vv<-vecs[[nm]]; co<-vv$co
  delta<-dens_on_grid(vv$group,grid)-dens_on_grid(vv$ref,grid)
  res[[length(res)+1L]]<-data.table(x=grid,delta_density=delta,time=co$time,contrast=co$contrast)
  if(B>0){
    n<-min(length(vv$group),length(vv$ref))
    mat_boot<-matrix(NA_real_, nrow=length(grid), ncol=B)
    for(b in seq_len(B)){
      idx<-sample.int(n,size=n,replace=TRUE)
      mat_boot[,b]<-dens_on_grid(vv$group[idx],grid)-dens_on_grid(vv$ref[idx],grid)
    }
    boot_res[[length(boot_res)+1L]]<-data.table(
      x=grid,
      lo=apply(mat_boot,1,quantile,probs=.025,na.rm=TRUE),
      hi=apply(mat_boot,1,quantile,probs=.975,na.rm=TRUE),
      time=co$time,contrast=co$contrast)
  }
}
plot_dt<-rbindlist(res)
if(length(boot_res)) plot_dt<-merge(plot_dt,rbindlist(boot_res),by=c("x","time","contrast"),all.x=TRUE) else plot_dt[,`:=`(lo=NA_real_,hi=NA_real_)]
plot_dt[, time:=factor(time, levels=c("4 h","24 h"))]
plot_dt[, contrast:=factor(contrast, levels=c("Dox-only","Tam-only"))]
vals<-c(plot_dt$delta_density,plot_dt$lo,plot_dt$hi)
ymax<-if(is.finite(fixed_ymax)) fixed_ymax else max(min_ymax,max(abs(vals[is.finite(vals)]),na.rm=TRUE)*y_pad)
ylim<-c(-ymax,ymax)

p <- ggplot(plot_dt, aes(x=x, y=delta_density)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=.55, color="grey35") +
  geom_ribbon(aes(ymin=lo, ymax=hi), fill="#9ECAE1", alpha=.30, na.rm=TRUE) +
  geom_line(linewidth=.90, color="#2C7FB8", na.rm=TRUE) +
  facet_grid(contrast ~ time) +
  coord_cartesian(ylim=ylim) +
  theme_bw(base_size=18) +
  labs(
    title="Dox-only and Tam-only Δ-density controls",
    x="log2(size-factor-normalized RSEM expected count + 1)",
    y="Δ density (condition − Ctrl)"
  ) +
  theme(
    plot.title=element_text(face="bold",size=22,hjust=.5, margin=margin(b=4)),
    strip.text=element_text(face="bold",size=17),
    axis.title=element_text(size=17),
    axis.text=element_text(size=14),
    panel.grid.minor=element_blank(),
    plot.caption=element_text(size=12,color="grey25",hjust=.5, margin=margin(t=4)),
    plot.margin=margin(8,10,8,10)
  )

out_png<-file.path(out_dir,"Fig2A_delta_density_controls.png")
out_pdf<-file.path(out_dir,"Fig2A_delta_density_controls.pdf")
out_svg<-file.path(out_dir,"Fig2A_delta_density_controls.svg")
ggsave(out_png,p,width=width_in,height=height_in,dpi=300)
ggsave(out_pdf,p,width=width_in,height=height_in,device=cairo_pdf)
ggsave(out_svg,p,width=width_in,height=height_in)
fwrite(plot_dt,file.path(out_dir,"Fig2A_delta_density_controls_source.csv"))
fwrite(data.table(input=expr_path,scale="size-factor-normalized RSEM expected counts",bootstrap_B=B),file.path(out_dir,"Fig2A_delta_density_controls_scale.tsv"),sep="\t")
message("[12] wrote: ",out_png)
message("[12] wrote: ",out_pdf)
message("[12] wrote: ",out_svg)
