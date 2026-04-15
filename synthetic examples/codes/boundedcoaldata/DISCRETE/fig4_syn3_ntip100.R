library(expm)
library(ape)
remotes::install_github("JuliaPalacios/phylodyn")

########################################################Ne_3
data_ori=c(
  0.00394216, 0.00476293, 0.00683948, 0.0200167, 0.02287877,
  0.02494215, 0.03330005, 0.04279202, 0.04464086, 0.05133851,
  0.06159588, 0.06361727, 0.06367681, 0.07908421, 0.08059573,
  0.09299994, 0.09443599, 0.09549568, 0.09588555, 0.10405885,
  0.11171124, 0.11418793, 0.12270596, 0.12950722, 0.13497745,
  0.14037071, 0.1523515, 0.15279899, 0.15396341, 0.16147663,
  0.16673145, 0.17032533, 0.17050465, 0.17246076, 0.17604802,
  0.17896829, 0.17919375, 0.1830164, 0.18728525, 0.18885518,
  0.19792073, 0.20166918, 0.20608097, 0.20841984, 0.21629147,
  0.21767764, 0.22033493, 0.22826348, 0.22856353, 0.24337846,
  0.24869935, 0.25186855, 0.25229314, 0.25921076, 0.26590992,
  0.27000209, 0.2832244, 0.28461947, 0.29129408, 0.30368242,
  0.3070937, 0.30828093, 0.31094398, 0.31160093, 0.31274821,
  0.33824175, 0.34456313, 0.34837455, 0.34847975, 0.35691756,
  0.36103201, 0.36145358, 0.36420368, 0.37016777, 0.38164857,
  0.3841142, 0.39028157, 0.40819351, 0.40954445, 0.41293516,
  0.41318948, 0.41327992, 0.41510723, 0.44342549, 0.44921248,
  0.47056455, 0.47564584, 0.48771615, 0.50843228, 0.54641561,
  0.55410468, 0.57840887, 0.61957205, 0.64408916, 0.65992476,
  0.66517345, 0.68430376, 0.70009952, 0.70644973
)
data <- list(coal_times=data_ori,lineages= 100:2,intercoal_times= c(
  data_ori[1],
  data_ori[-1] - data_ori[-length(data_ori)] ), samp_times=0,n_sampled=100)
bound=0.71

inten2 <- function(t) {
  25 * exp(-5*t)
}

x = seq(0, bound+1e-4, length.out = 101)
intensity=sapply(x,inten2)

res1<-phylodyn:::bounded_skyline_ascent(data, bound =bound)
tree<-phylodyn:::generate_newick(data)$newick
pdf("plot30.pdf")
skylineplot(tree)
points(x, intensity, type = "l",
       col = "red", lwd = 3,xlab = "Time",
       ylab = expression(N[e](t)))
lines(res1$grid, res1$Ne, col="blue", lwd=2, type="s")

legend("topleft",
       legend = c("Truth","SC-MLE", "BC-MLE"),
       col = c("red","black", "blue"),
       lwd = c(1.5, 3, 2, 2),
       lty = c(1, 1, 1, 1))
points(data_ori,
       rep(5e-3, length(data_ori)),
       pch = 124, # vertical line marker '|'
       col = "black",
       cex = 0.6) # controls size (≈ markersize)
dev.off()



res1_mcmc<-phylodyn:::mcmc_sampling(data,alg="bound_ESS",nsamp=200000,nburnin=100000,ngrid=100,bound=bound)
res2_mcmc<-phylodyn:::mcmc_sampling(data,alg="ESS",nsamp=200000,ngrid=100,nburnin=100000)


pdf("plot31.pdf")
# plot truth
plot(x, intensity, type = "l",
     ylim = c(0, 70),
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


pdf("plot32.pdf")
plot(res1_mcmc$samp[,50],type = "l")
#plot(res1_mcmc$pos_summ$loglik)
dev.off()




