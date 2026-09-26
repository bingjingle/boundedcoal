# Invoke with Rscript; archived input paths may be overridden with an environment variable.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_file <- gsub("~+~", " ", sub("^--file=", "", script_arg), fixed = TRUE)
script_dir <- dirname(normalizePath(script_file))
data_dir <- Sys.getenv("BOUNDEDCOAL_CELL_LINEAGE_DATA_DIR", script_dir)
output_dir <- Sys.getenv("BOUNDEDCOAL_CELL_LINEAGE_OUTPUT_DIR", file.path(script_dir, "../../outputs/cell_lineage"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library("phylodyn")
#library("INLA")
#install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
set.seed(123)

load(file.path(data_dir, "UGPMA3257new.rda"))
points_inhomo=bt_adj
ntip=3257
bins=100
out<-BNPR(list(coal_times=points_inhomo,samp_times=c(0),n_sampled=c(ntip)), prec_alpha = .1,
          prec_beta = .1)
plot_BNPR(out)
x=seq(0,max(points_inhomo),length.out=101)[-1]
low=out$summary$quant0.025
med=out$summary$quant0.5
high=out$summary$quant0.975

plot(x,med,type='l',ylim=c(0,1000000))
lines(x,high)
lines(x,low)
points(points_inhomo, rep(0.1,ntip-1))
