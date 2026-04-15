
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

res1_mcmc<-phylodyn:::mcmc_sampling(data,alg="bound_ESS",nsamp= 200000,nburnin=100000,ngrid=100,f_init=log(bnpr$effpop)[-1],bound=bound)
res2_mcmc<-phylodyn:::mcmc_sampling(data,alg="ESS",nsamp=200000,ngrid=100,nburnin=100000)


pdf("plot11.pdf")
# plot truth
plot(x, intensity, type = "l",
     ylim = c(0, 8),
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


pdf("plot12.pdf")
plot(res1_mcmc$samp[,50],type = "l")
#plot(res1_mcmc$pos_summ$loglik)
dev.off()





