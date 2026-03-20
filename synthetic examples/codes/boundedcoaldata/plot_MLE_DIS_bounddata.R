install.packages("devtools")
library(devtools)
install.packages("remotes")
install.packages("expm")
library(expm)
library(ape)
remotes::install_github("JuliaPalacios/phylodyn")

##################################################Ne_1
data_ori=c(
  0.0001266, 0.00032191, 0.00044205, 0.00088579, 0.00104198,
  0.00113909, 0.0012439, 0.00130914, 0.00146451, 0.00152415,
  0.00153834, 0.00155985, 0.00166426, 0.00202518, 0.00212264,
  0.0021716, 0.00264311, 0.00275249, 0.00347892, 0.00352986,
  0.00357284, 0.00369775, 0.00407869, 0.00460924, 0.00517073,
  0.00582611, 0.00601573, 0.00619605, 0.00683839, 0.0073865,
  0.00755511, 0.00771401, 0.00815352, 0.00818267, 0.0085899,
  0.00865567, 0.00889056, 0.00890924, 0.00997735, 0.01011257,
  0.01026728, 0.0104445, 0.01112186, 0.01130533, 0.01135701,
  0.0120875, 0.01263979, 0.01282922, 0.01304582, 0.01318416,
  0.01320701, 0.01428112, 0.01625349, 0.01693358, 0.01964695,
  0.0207158, 0.02517708, 0.02550328, 0.0275826, 0.02949231,
  0.03126972, 0.03651252, 0.03702314, 0.0382631, 0.04201485,
  0.04752628, 0.04792048, 0.04853976, 0.05543565, 0.05551773,
  0.05566007, 0.05586548, 0.05881217, 0.06054298, 0.06377044,
  0.06630486, 0.07089239, 0.07276823, 0.07817595, 0.08104758,
  0.08890323, 0.10267324, 0.10405011, 0.1133765, 0.12528137,
  0.1276583, 0.14280372, 0.14369328, 0.1956738, 0.20957733,
  0.22259239, 0.23262752, 0.23351797, 0.23974517, 0.41920302,
  0.60871934, 0.69770135, 0.70704948, 0.91230628
)
data <- list(coal_times=data_ori,lineages= 100:2,intercoal_times= c(
  data_ori[1],
  data_ori[-1] - data_ori[-length(data_ori)]  ), samp_times=0,n_sampled=100)
bound=1


inten2<-function(t){
  return (1+t-t)}

x = seq(0, bound+1e-4, length.out = 101)
intensity=sapply(x,inten2)
res1<-phylodyn:::bounded_skyline_ascent(data, bound =bound)
tree<-phylodyn:::generate_newick(data)$newick
skylineplot(tree)
pdf("plot10.pdf")
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
       rep(0.025, length(data_ori)),
       pch = 124,      # vertical line marker '|'
       col = "black",
       cex = 0.6)      # controls size (≈ markersize)
dev.off()


bnpr<-phylodyn:::BNPR(tree)
phylodyn:::plot_BNPR(bnpr)
abline(h=1, col="red", lwd=2)

res1_mcmc<-phylodyn:::mcmc_sampling(data,alg="bound_ESS",nsamp=50000,nburnin=500,ngrid=100,f_init=log(bnpr$effpop)[-1],bound=1)
res2_mcmc<-phylodyn:::mcmc_sampling(data,alg="ESS",nsamp=50000,ngrid=100,nburnin=500)

pdf("plot11.pdf")
# plot truth
plot(x, intensity, type = "l",
     ylim = c(0, 10),
     col = "black", lwd = 3,xlab = "Time",
     ylab = expression(N[e](t)))
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


#######################################Ne_2
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
skylineplot(tree)
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


bnpr<-phylodyn:::BNPR(tree)
phylodyn:::plot_BNPR(bnpr)
res1_mcmc<-phylodyn:::mcmc_sampling(data,alg="bound_ESS",nsamp=50000,nburnin=500,ngrid=100,f_init=log(bnpr$effpop)[-1],bound=1)
res2_mcmc<-phylodyn:::mcmc_sampling(data,alg="ESS",nsamp=50000,ngrid=100,nburnin=500)

pdf("plot21.pdf")
# plot truth
plot(x, intensity, type = "l",
     ylim = c(0, 18),
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
  data_ori[-1] - data_ori[-length(data_ori)]  ), samp_times=0,n_sampled=100)
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
       legend = c("Truth","SC-MLE",  "BC-MLE"),
       col    = c("red","black", "blue"),
       lwd    = c(1.5, 3, 2, 2),
       lty    = c(1, 1, 1, 1))
points(data_ori,
       rep(5e-3, length(data_ori)),
       pch = 124,      # vertical line marker '|'
       col = "black",
       cex = 0.6)      # controls size (≈ markersize)

dev.off()


bnpr<-phylodyn:::BNPR(tree)
phylodyn:::plot_BNPR(bnpr)
points(seq(0, max(res1$grid), by = .05),
       25 * exp(-5*seq(0, max(res1$grid), by = .05)),
       type = "l", col = "red"
res1_mcmc<-phylodyn:::mcmc_sampling(data,alg="bound_ESS",nsamp=50000,nburnin=500,ngrid=100,f_init=log(bnpr$effpop)[-1],bound=1)
res2_mcmc<-phylodyn:::mcmc_sampling(data,alg="ESS",nsamp=50000,ngrid=100,nburnin=500)

pdf("plot31.pdf")
# plot truth
plot(x, intensity, type = "l",
     ylim = c(0, 120),
     col = "black", lwd = 3,xlab = "Time",
     ylab = expression(N[e](t)))
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


