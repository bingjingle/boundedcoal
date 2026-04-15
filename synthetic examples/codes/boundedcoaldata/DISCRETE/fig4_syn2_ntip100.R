
library(expm)
library(ape)
remotes::install_github("JuliaPalacios/phylodyn")

data_ori=c(
  0.0006457, 0.00094795, 0.00097818, 0.00115723, 0.00136032,
  0.0014971, 0.00188022, 0.00229853, 0.00275317, 0.00414087,
  0.00510088, 0.00546055, 0.00724715, 0.0075598, 0.00843643,
  0.01121559, 0.01140645, 0.01150385, 0.01385105, 0.01563069,
  0.01782013, 0.01795062, 0.01827905, 0.01898047, 0.01989265,
  0.01997209, 0.02033843, 0.02223489, 0.02405914, 0.02455844,
  0.02490612, 0.02514106, 0.02523667, 0.0254744, 0.02547566,
  0.02562561, 0.02600492, 0.02673282, 0.02833923, 0.02908865,
  0.03005297, 0.03110919, 0.03234614, 0.03744378, 0.03752536,
  0.03869813, 0.04311578, 0.04386097, 0.04909827, 0.05170783,
  0.05189045, 0.05292001, 0.05726541, 0.05837406, 0.06083225,
  0.06454715, 0.06758013, 0.07431144, 0.07451043, 0.07548091,
  0.07577013, 0.07728782, 0.08168786, 0.08383543, 0.08460874,
  0.0851407, 0.08686686, 0.08687484, 0.08774483, 0.08805298,
  0.09739123, 0.10208261, 0.10602849, 0.11395216, 0.13507485,
  0.14032089, 0.18690144, 0.19247492, 0.19934897, 0.22431027,
  0.22508206, 0.2307226, 0.23352532, 0.23890136, 0.24700717,
  0.25998017, 0.27098572, 0.28242977, 0.28818101, 0.2978475,
  0.33223878, 0.40785599, 0.42174163, 0.48693148, 0.49802362,
  0.52905337, 0.58797086, 0.62521005, 0.66886555
)
data <- list(coal_times=data_ori,lineages= 100:2,intercoal_times= c(
  data_ori[1],
  data_ori[-1] - data_ori[-length(data_ori)]  ), samp_times=0,n_sampled=100)
bound=0.7
inten2 <- function(t) {
  3 * exp(-t)
}


x = seq(0, bound+1e-4, length.out = 101)
intensity=sapply(x,inten2)




res1<-phylodyn:::bounded_skyline_ascent(data, bound =bound)
tree<-phylodyn:::generate_newick(data)$newick

pdf("plot20.pdf")
skylineplot(tree)
points(x, intensity, type = "l",
       col = "red", lwd = 3,xlab = "Time",
       ylab = expression(N[e](t)))
lines(res1$grid, res1$Ne, col="blue", lwd=2, type="s")

legend("topleft",
       legend = c("Truth","SC-MLE",  "BC-MLE"),
       col    = c("red","black", "blue"),
       lwd    = c(1.5, 3, 2, 2),
       lty    = c(1, 1, 1, 1))
points(data_ori,
       rep(3e-3, length(data_ori)),
       pch = 124,      # vertical line marker '|'
       col = "black",
       cex = 0.6)      # controls size (≈ markersize)
dev.off()


res1_mcmc<-phylodyn:::mcmc_sampling(data,alg="bound_ESS",nsamp= 200000,nburnin=100000,ngrid=100,f_init=log(bnpr$effpop)[-1],bound=bound)
res2_mcmc<-phylodyn:::mcmc_sampling(data,alg="ESS",nsamp=200000,ngrid=100,nburnin=100000)


pdf("plot21.pdf")
# plot truth
plot(x, intensity, type = "l",
     ylim = c(0, 5),
     col = "black", lwd = 3,xlab = "Time",
     ylab = expression(N[e](t)))
# magenta color with alpha (for lines)
lines(res1_mcmc$med_fun,pch="",col="magenta",lwd=2.5)
lines(res1_mcmc$low_fun,pch="",col="magenta",lwd=2,lty=3)
lines(res1_mcmc$hi_fun,pch="",col="magenta",lwd=2,lty=3)

lines(res2_mcmc$med_fun,pch="",col="cyan",lwd=2.5)
lines(res2_mcmc$low_fun,pch="",col="cyan",lwd=2,lty=3)
lines(res2_mcmc$hi_fun,pch="",col="cyan",lwd=2,lty=3)


legend("topright",
       legend = c("Truth", "BC-DISCRETIZED", "SC-DISCRETIZED"),
       col    = c("black", "magenta", "cyan"),
       lwd    = c(3, 2.5, 2.5),
       lty    = c(1, 1, 1))
points(data_ori,
       rep(0, length(data_ori)),
       pch = 124,      # vertical line marker '|'
       col = "black",
       cex = 0.6)      # controls size (≈ markersize)

dev.off()


pdf("plot22.pdf")
plot(res1_mcmc$samp[,50],type = "l")
#plot(res1_mcmc$pos_summ$loglik)
dev.off()





