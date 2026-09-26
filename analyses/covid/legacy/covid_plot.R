##real data analysis

library("ape")
tree<-read.nexus(Sys.getenv("BC_COVID_TREE", "analyses/covid/data/median_ccd0.tree"))
factor<-.001/.0012
data<-summarize_phylo(tree)
data$coal_times<-data$coal_times*factor
tree2 <- tree
tree2$edge.length <- tree$edge.length * factor   # e.g., k = 0.5 or k = 365, etc.
##tree2 or data, already have the scaled time



##Skyline plot without bound
ski<-skyline(tree2)
t0<-as.Date("2020-06-08") ##this is the date of sampling
plot(ski$time,c(ski$population.size),type="S",xlim=c(0.6,0),xaxt="n",ylab="Ne",xlab="Time")
x_ticks<-pretty(c(0.6,0))
tick_dates<-t0-round((x_ticks*365))
axis(1,at=x_ticks,labels=paste0(format(tick_dates,"%b"),"'",format(tick_dates,"%y")))

#skyline with bound
bound<-max(data$coal_times)
library("phylodyn")
res1<-phylodyn:::bounded_skyline_ascent(data, bound =bound)
points(res1$grid,res1$Ne,type="S",ylab="Ne(t)",xlab="Time",ylim=c(0,20),col="red")

##report from Washington state of public heatlh
#Page 7
#https://doh.wa.gov/sites/default/files/2023-05/421038-2020Covid19AnnualReport.pdf

df <- data.frame(
  month = factor(c("Jan","Feb","Mar","Apr","May","Jun"),
                 levels = c("Jan","Feb","Mar","Apr","May","Jun"),
                 ordered = TRUE),
  cases = c(20, 44, 6947, 8698, 7315, 13662)
)

barplot(df$cases, names.arg=df$month, col="steelblue", border=NA,
        xlab="Month of 2020", ylab="Cases", cex.names=0.8)


##Discretized
res1_mcmc<-mcmc_sampling(data,alg="bound_ESS",nsamp=30000,nburnin=3000,ngrid=100,bound=bound)
res2_mcmc<-mcmc_sampling(data,alg="ESS",nsamp=30000,nburnin=3000,ngrid=100)

data1<-data.frame(x=res1_mcmc$x,bc_med=res1_mcmc$med,bc_low_eq=res1_mcmc$low,bc_hi_eq=res1_mcmc$hi,sc_med=res2_mcmc$med,sc_low_eq=res2_mcmc$low,sc_hi_eq=res2_mcmc$hi)
write.csv(data1, Sys.getenv("BC_COVID_COORDS", "analyses/reproduction/coords/panel_covid_DIS.csv"), row.names = FALSE)

par(mfrow=c(1,1))
plot(ski$time,c(ski$population.size),type="S",xlim=c(0.6,0),xaxt="n",ylab="Ne",xlab="Time",ylim=c(0,10),col="white")
x_ticks<-pretty(c(0.6,0))
tick_dates<-t0-round((x_ticks*365))
axis(1,at=x_ticks,labels=paste0(format(tick_dates,"%b"),"'",format(tick_dates,"%y")))


polygon(c(res1_mcmc$x, rev(res1_mcmc$x)), c(res1_mcmc$low, rev(res1_mcmc$hi)),
        col = adjustcolor("blue", alpha.f = 0.25),
        border = NA)                                 # shaded 95% region
lines(res1_mcmc$x, res1_mcmc$med, col="blue", lwd=2)                  # median/mean curve






grid = BNPR_out$grid
if (is.null(xlim)) {
  xlim = c(max(grid), min(grid))
}
mask = BNPR_out$x >= min(xlim) & BNPR_out$x <= max(xlim)
t = BNPR_out$x[mask]
y = BNPR_out$effpop[mask] * yscale
yhi = BNPR_out$effpop975[mask] * yscale
ylo = BNPR_out$effpop025[mask] * yscale
if (newplot) {
  if (is.null(ylim)) {
    ymax = max(yhi)
    ymin = min(ylo)
  }
  else {
    ymin = min(ylim)
    ymax = max(ylim)
  }
  if (heatmaps) {
    yspan = ymax/ymin
    yextra = yspan^(1/10)
    ylim = c(ymin/(yextra^1.35), ymax)
  }
  else {
    ylim = c(ymin, ymax)
  }
  if (is.null(axlabs)) {
    graphics::plot(1, 1, type = "n", log = log, xlab = xlab, 
                   ylab = ylab, main = main, xlim = xlim, ylim = ylim, 
                   ...)
  }
  else {
    graphics::plot(1, 1, type = "n", log = log, xlab = "", 
                   ylab = ylab, main = main, xlim = xlim, ylim = ylim, 
                   xaxt = "n", ...)
    graphics::axis(1, at = axlabs$x, labels = axlabs$labs, 
                   las = 2)
    graphics::mtext(text = xlab, side = 1, line = xmarline)
  }
}
if (credible_region) {
  shade_band(x = t, ylo = ylo, yhi = yhi, col = "lightgray")
}
if (!is.null(traj)) {
  graphics::lines(t, traj(t), lwd = traj_lwd, lty = traj_lty, 
                  col = traj_col)
}
if (newplot) {
  if (heatmaps) {
    samps = rep(BNPR_out$samp_times, BNPR_out$n_sampled)
    samps = samps[samps <= max(xlim) & samps >= min(xlim)]
    coals = BNPR_out$coal_times
    coals = coals[coals <= max(xlim) & coals >= min(xlim)]
    breaks = seq(min(xlim), max(xlim), length.out = nbreaks)
    h_samp = graphics::hist(samps, breaks = breaks, plot = FALSE)
    h_coal = graphics::hist(coals, breaks = breaks, plot = FALSE)
    hist2heat(h_samp, y = ymin/yextra^0.5, wd = heatmap_width)
    hist2heat(h_coal, y = ymin/yextra, wd = heatmap_width)
    if (heatmap_labels) {
      if (heatmap_labels_side == "left") {
        lab_x = max(xlim)
        lab_adj = 0
      }
      else if (heatmap_labels_side == "right") {
        lab_x = min(xlim)
        lab_adj = 1
      }
      else {
        warning("heatmap_labels_side not \"left\" or \"right\", defaulting to right")
        lab_x = min(xlim)
        lab_adj = 1
      }
      graphics::text(x = lab_x, y = ymin/(yextra^0.2), 
                     labels = "Sampling events", adj = c(lab_adj, 
                                                         0), cex = heatmap_labels_cex)
      graphics::text(x = lab_x, y = ymin/(yextra^1.25), 
                     labels = "Coalescent events", adj = c(lab_adj, 
                                                           1), cex = heatmap_labels_cex)
    }
  }
}
graphics::lines(t, y, lwd = lwd, col = col, lty = lty)
}

    