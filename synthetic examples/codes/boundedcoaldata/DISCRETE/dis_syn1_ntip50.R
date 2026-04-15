library(expm)
library(ape)
remotes::install_github("JuliaPalacios/phylodyn")

########################################################Ne_3
args = commandArgs(trailingOnly=TRUE)

set.seed(123)
nam <- paste0("/home/groups/juliapr/Bingjing/code/boundedcoal/syndatafinal/syn1_bounded_ntip50.rda")#args[1] range from 1 to 300
load(nam)

idx <- as.integer(commandArgs(trailingOnly = TRUE)[1])   # 1 to 30
data_ori <- coal_time[idx,]
ntip=50
Ngrid=100
data <- list(coal_times=data_ori,lineages= ntip:2,intercoal_times= c(
  data_ori[1],
  data_ori[-1] - data_ori[-length(data_ori)] ), samp_times=0,n_sampled=ntip)


bound=1

inten2<-function(t){
  return (1+t-t)}

x4 <- seq(0, bound, length.out = Ngrid + 1)[-1]
truth=sapply(x4,inten2)

start_time <- Sys.time()
res1_mcmc<-phylodyn:::mcmc_sampling(data,alg="bound_ESS",nsamp=200000,nburnin=100000,ngrid=100,bound=bound)
timerun <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
timerun_10000_bounded <- timerun / 200000 * 10000

y_low<-res1_mcmc$low_fun(x4)
y_med<-res1_mcmc$med_fun(x4)
y_hi<-res1_mcmc$hi_fun(x4)

l2_dist1_bound <- sum((y_med - truth)^2)
coverage1_bound <- sum((truth >= y_low) & (truth <= y_hi)) / length(x4)
width1_bound <- sum(y_hi - y_low) / Ngrid



start_time <- Sys.time()
res2_mcmc<-phylodyn:::mcmc_sampling(data,alg="ESS",nsamp=200000,ngrid=100,nburnin=100000)
timerun <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
timerun_10000_std <- timerun / 200000 * 10000


y_low<-res2_mcmc$low_fun(x4)
y_med<-res2_mcmc$med_fun(x4)
y_hi<-res2_mcmc$hi_fun(x4)

l2_dist1_std <- sum((y_med - truth)^2)
coverage1_std <- sum((truth >= y_low) & (truth <= y_hi)) / length(x4)
width1_std <- sum(y_hi - y_low) / Ngrid



outfile <- paste0(  "/scratch/groups/juliapr/output_Bingjing/boundedcoaldata/dis/syn1ntip50/",  "mcmc_results_", idx, ".rda")
save(  timerun_10000_bounded, res1_mcmc, l2_dist1_bound, coverage1_bound, width1_bound,  timerun_10000_std, res2_mcmc, l2_dist1_std, coverage1_std, width1_std,  file = outfile)















