# Invoke with Rscript; locate inputs relative to this script.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_file <- gsub("~+~", " ", sub("^--file=", "", script_arg), fixed = TRUE)
script_dir <- dirname(normalizePath(script_file))
synthetic_root <- normalizePath(file.path(script_dir, "../.."))
data_dir <- Sys.getenv("BOUNDEDCOAL_SYNTHETIC_DATA_DIR", file.path(synthetic_root, "data"))
output_dir <- Sys.getenv("BOUNDEDCOAL_SYNTHETIC_OUTPUT_DIR", file.path(synthetic_root, "../../outputs/synthetic"))
output_file <- function(name) {
  path <- file.path(output_dir, name)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

pdf(output_file("standard/INLA_diagnostics.pdf"))
library("phylodyn")
#library("INLA")
#install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
set.seed(123)

load(file.path(data_dir, "syn1_std_ntip100.rda"))
l2_dist1=rep(0,30)
width1=rep(0,30)
coverage1=rep(0,30)
time_1=rep(0,30)
inten2<-function(t){
  return (1+t-t)}
for (jjj in 1:30){
  points_inhomo=coal_time[jjj,]
  ntip=100
  bins=100
  out<-BNPR(list(coal_times=points_inhomo,samp_times=c(0),n_sampled=c(ntip)), prec_alpha = .1,
            prec_beta = .1)
  plot_BNPR(out)
  x=seq(0,max(points_inhomo),length.out=101)[-1]
  y=sapply(x,inten2)
  lines(x,y)
  low=out$summary$quant0.025
  med=out$summary$quant0.5
  high=out$summary$quant0.975
  width=sum(high-low)/bins
  intensity=sapply(x,inten2)
  plot(x, intensity, type = "l", ylim = c(0, 3), col = "red")
  lines(x,low)
  lines(x,med)
  lines(x,high)
  points(points_inhomo, rep(0.1,ntip-1))
  coverage=sum((intensity>=low)*(intensity<=high))/bins
  l2dist_med=sum((intensity-med)**2)
  coverage1[jjj]=coverage
  l2_dist1[jjj]=l2dist_med
  width1[jjj]=width
  print(c(jjj,l2dist_med,coverage,width,max(points_inhomo)))
  printed_output <- capture.output(out$result)
  out_str <- paste(printed_output, collapse = " ")
  # Use a regex to capture the number after "Total = "
  total_time <- sub(".*Total = ([0-9\\.]+).*", "\\1", out_str)
  time_1[jjj]=as.numeric(total_time)
}
round(quantile(l2_dist1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(quantile(100*coverage1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),0)
round(quantile(width1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_1),2)
round(sd(time_1),2)

##################################################################################
load(file.path(data_dir, "syn1_std_ntip50.rda"))
l2_dist1=rep(0,30)
width1=rep(0,30)
coverage1=rep(0,30)
time_1=rep(0,30)

for (jjj in 1:30){
  points_inhomo=coal_time[jjj,]
  ntip=50
  bins=100
  out<-BNPR(list(coal_times=points_inhomo,samp_times=c(0),n_sampled=c(ntip)), prec_alpha = .1,
            prec_beta = .1)
  plot_BNPR(out)
  x=seq(0,max(points_inhomo),length.out=101)[-1]
  y=sapply(x,inten2)
  lines(x,y)
  low=out$summary$quant0.025
  med=out$summary$quant0.5
  high=out$summary$quant0.975
  width=sum(high-low)/bins
  intensity=sapply(x,inten2)
  plot(x, intensity, type = "l", ylim = c(0, 3), col = "red")
  lines(x,low)
  lines(x,med)
  lines(x,high)
  points(points_inhomo, rep(0.1,ntip-1))
  coverage=sum((intensity>=low)*(intensity<=high))/bins
  l2dist_med=sum((intensity-med)**2)
  coverage1[jjj]=coverage
  l2_dist1[jjj]=l2dist_med
  width1[jjj]=width
  print(c(jjj,l2dist_med,coverage,width,max(points_inhomo)))
  printed_output <- capture.output(out$result)
  out_str <- paste(printed_output, collapse = " ")
  # Use a regex to capture the number after "Total = "
  total_time <- sub(".*Total = ([0-9\\.]+).*", "\\1", out_str)
  time_1[jjj]=as.numeric(total_time)
}
round(quantile(l2_dist1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(quantile(100*coverage1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),0)
round(quantile(width1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_1),2)
round(sd(time_1),2)

##################################################################################
load(file.path(data_dir, "syn2_std_ntip50.rda"))
l2_dist1=rep(0,30)
width1=rep(0,30)
coverage1=rep(0,30)
time_1=rep(0,30)
inten2 <- function(t) {
  3 * exp(-t)
}
for (jjj in 1:30){
  points_inhomo=coal_time[jjj,]
  ntip=50
  bins=100
  out<-BNPR(list(coal_times=points_inhomo,samp_times=c(0),n_sampled=c(ntip)), prec_alpha = .1,
            prec_beta = .1)
  plot_BNPR(out)
  x=seq(0,max(points_inhomo),length.out=101)[-1]
  y=sapply(x,inten2)
  lines(x,y)
  low=out$summary$quant0.025
  med=out$summary$quant0.5
  high=out$summary$quant0.975
  width=sum(high-low)/bins
  intensity=sapply(x,inten2)
  plot(x, intensity, type = "l", ylim = c(0, 3), col = "red")
  lines(x,low)
  lines(x,med)
  lines(x,high)
  points(points_inhomo, rep(0.1,ntip-1))
  coverage=sum((intensity>=low)*(intensity<=high))/bins
  l2dist_med=sum((intensity-med)**2)
  coverage1[jjj]=coverage
  l2_dist1[jjj]=l2dist_med
  width1[jjj]=width
  print(c(jjj,l2dist_med,coverage,width,max(points_inhomo)))
  printed_output <- capture.output(out$result)
  out_str <- paste(printed_output, collapse = " ")
  # Use a regex to capture the number after "Total = "
  total_time <- sub(".*Total = ([0-9\\.]+).*", "\\1", out_str)
  time_1[jjj]=as.numeric(total_time)
}

round(quantile(l2_dist1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(quantile(100*coverage1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),0)
round(quantile(width1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_1),2)
round(sd(time_1),2)


##################################################################################
load(file.path(data_dir, "syn2_std_ntip100.rda"))
l2_dist1=rep(0,30)
width1=rep(0,30)
coverage1=rep(0,30)
time_1=rep(0,30)

for (jjj in 1:30){
  points_inhomo=coal_time[jjj,]
  ntip=100
  bins=100
  out<-BNPR(list(coal_times=points_inhomo,samp_times=c(0),n_sampled=c(ntip)), prec_alpha = .1,
            prec_beta = .1)
  plot_BNPR(out)
  x=seq(0,max(points_inhomo),length.out=101)[-1]
  y=sapply(x,inten2)
  lines(x,y)
  low=out$summary$quant0.025
  med=out$summary$quant0.5
  high=out$summary$quant0.975
  width=sum(high-low)/bins
  intensity=sapply(x,inten2)
  plot(x, intensity, type = "l", ylim = c(0, 3), col = "red")
  lines(x,low)
  lines(x,med)
  lines(x,high)
  points(points_inhomo, rep(0.1,ntip-1))
  coverage=sum((intensity>=low)*(intensity<=high))/bins
  l2dist_med=sum((intensity-med)**2)
  coverage1[jjj]=coverage
  l2_dist1[jjj]=l2dist_med
  width1[jjj]=width
  print(c(jjj,l2dist_med,coverage,width,max(points_inhomo)))
  printed_output <- capture.output(out$result)
  out_str <- paste(printed_output, collapse = " ")
  # Use a regex to capture the number after "Total = "
  total_time <- sub(".*Total = ([0-9\\.]+).*", "\\1", out_str)
  time_1[jjj]=as.numeric(total_time)
}

round(quantile(l2_dist1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(quantile(100*coverage1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),0)
round(quantile(width1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_1),2)
round(sd(time_1),2)

##################################################################################
load(file.path(data_dir, "syn3_std_ntip100.rda"))
l2_dist1=rep(0,30)
width1=rep(0,30)
coverage1=rep(0,30)
time_1=rep(0,30)
inten2 <- function(t) {
  25 * exp(-5*t)
}
jjj=1
for (jjj in 1:30){
  points_inhomo=coal_time[jjj,]
  ntip=100
  bins=100
  out<-BNPR(list(coal_times=points_inhomo,samp_times=c(0),n_sampled=c(ntip)), prec_alpha = .1,
            prec_beta = .1)
  plot_BNPR(out)
  x=seq(0,max(points_inhomo),length.out=101)[-1]
  y=sapply(x,inten2)
  lines(x,y)
  low=out$summary$quant0.025
  med=out$summary$quant0.5
  high=out$summary$quant0.975
  width=sum(high-low)/bins
  intensity=sapply(x,inten2)
  plot(x, intensity, type = "l", ylim = c(0, 100), col = "red")
  lines(x,low)
  lines(x,med)
  lines(x,high)
  points(points_inhomo, rep(0.1,ntip-1))
  coverage=sum((intensity>=low)*(intensity<=high))/bins
  l2dist_med=sum((intensity-med)**2)
  coverage1[jjj]=coverage
  l2_dist1[jjj]=l2dist_med
  width1[jjj]=width
  print(c(jjj,l2dist_med,coverage,width,max(points_inhomo)))
  printed_output <- capture.output(out$result)
  out_str <- paste(printed_output, collapse = " ")
  # Use a regex to capture the number after "Total = "
  total_time <- sub(".*Total = ([0-9\\.]+).*", "\\1", out_str)
  time_1[jjj]=as.numeric(total_time)
}

round(quantile(l2_dist1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(quantile(100*coverage1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),0)
round(quantile(width1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_1),2)
round(sd(time_1),2)

##################################################################################
load(file.path(data_dir, "syn3_std_ntip50.rda"))
l2_dist1=rep(0,30)
width1=rep(0,30)
coverage1=rep(0,30)
time_1=rep(0,30)
for (jjj in 1:30){
  points_inhomo=coal_time[jjj,]
  ntip=50
  bins=100
  out<-BNPR(list(coal_times=points_inhomo,samp_times=c(0),n_sampled=c(ntip)), prec_alpha = .1,
            prec_beta = .1)
  plot_BNPR(out)
  x=seq(0,max(points_inhomo),length.out=101)[-1]
  y=sapply(x,inten2)
  lines(x,y)
  low=out$summary$quant0.025
  med=out$summary$quant0.5
  high=out$summary$quant0.975
  width=sum(high-low)/bins
  intensity=sapply(x,inten2)
  plot(x, intensity, type = "l", ylim = c(0, 3), col = "red")
  lines(x,low)
  lines(x,med)
  lines(x,high)
  points(points_inhomo, rep(0.1,ntip-1))
  coverage=sum((intensity>=low)*(intensity<=high))/bins
  l2dist_med=sum((intensity-med)**2)
  coverage1[jjj]=coverage
  l2_dist1[jjj]=l2dist_med
  width1[jjj]=width
  print(c(jjj,l2dist_med,coverage,width,max(points_inhomo)))
  printed_output <- capture.output(out$result)
  out_str <- paste(printed_output, collapse = " ")
  # Use a regex to capture the number after "Total = "
  total_time <- sub(".*Total = ([0-9\\.]+).*", "\\1", out_str)
  time_1[jjj]=as.numeric(total_time)
}

round(quantile(l2_dist1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(quantile(100*coverage1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),0)
round(quantile(width1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_1),2)
round(sd(time_1),2)



dev.off()
