# Set BC_SIMULATION_OUTPUT to choose where generated benchmark data are saved.
output_dir <- Sys.getenv("BC_SIMULATION_OUTPUT", "generated")
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

####partial examples (rest are repeated, skipped)
########################Carson####################################################################
########################syn1######################################################################
########################ntip=100##################################################################
########################tau=0.5###################################################################
library(BoundedCoalescent)
library(ape)
Nsim_ori=3000
Nsim=3000
ntip=100
tau=0.5
start_time <- Sys.time()
coal_time1=replicate(ntip-1, numeric(Nsim) )
while(Nsim){
  t=bounded_sample_phylo(c(tau),c(ntip),1,0)
  coal_time1[Nsim,]=cumsum(coalescent.intervals(t$phylo)$interval.length)
  Nsim=Nsim-1
}
end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken <- as.numeric(time_taken, units = "secs")
time_taken/Nsim_ori
nam <- paste("carsonsyn1data_ntip", ntip, "_samplesize", Nsim_ori,"_tau", tau, sep = "")
save(coal_time1, file=file.path(output_dir, paste(nam,'.rda',sep='')))

##################Thinning#########################################################################
##################syn1#############################################################################
########################ntip=100###################################################################
########################tau=0.5####################################################################
ntip=100
samp_times = c(0)
n_sampled = c(ntip)
tau=0.5
Nsim_ori=3000
Nsim=3000
coal_time2=replicate(ntip-1, numeric(Nsim) )

r_func <- function(k, j) {
  if(j == 1) return(1)
  prod <- 1
  for(m in 1:(j-1)) {
    prod <- prod * ((k - m)/(k + m))
  }
  return((-1)^(j-1)*(2*j-1) * prod)
}

inten2 <- function(t) {
  1+t-t
}

inten2_inv <- function(t) {
  1+t-t
}
inten2_inv_cum<- function(t) {
  t
}

inten2_inv_cum_inv<- function(t) {
  t
}

result_list <- lapply(seq_len(ntip), function(k) {
  sapply(seq_len(k), function(i) r_func(k, i))
})

com_vec <- choose(1:ntip, 2)
start_time <- Sys.time()
while(Nsim){
  simulation4<-rep(0,ntip-1)
  t_upper<-0
  for (j in ntip:2){
    while(1){
      ori<-rexp(n=1, rate=choose(j,2))
      t_upper<-inten2_inv_cum_inv(inten2_inv_cum(tau)+ori-log(exp(inten2_inv_cum(tau)-inten2_inv_cum(t_upper))-1+exp(ori)))
      kk=j
      con=inten2_inv_cum(t_upper)-inten2_inv_cum(tau) #gamma(t)-gamma(\tau)
      denom=sum(result_list[[kk]]*exp(com_vec[1:(kk)]*con))
      numer=sum(result_list[[kk-1]]*exp(com_vec[1:(kk-1)]*con))

      if(stats::runif(1) <= numer/denom*(1-exp(con))) {
        simulation4[ntip+1-j]<-t_upper
        break
      }
    }
  }

  coal_time2[Nsim,]=simulation4
  Nsim=Nsim-1
}
end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken <- as.numeric(time_taken, units = "secs")
time_taken/Nsim_ori
nam <- paste("thinsyn1data_ntip",ntip, "_samplesize", Nsim_ori,"_tau", tau, sep = "")
save(coal_time2, file=file.path(output_dir, paste(nam,'.rda',sep='')))


##################Naive rejection#########################################################################
##################syn1#############################################################################
########################ntip=100###################################################################
########################tau=0.5####################################################################
ntip=100
tau=0.5
Nsim_ori=3000
Nsim=3000
coal_time3=replicate(ntip-1, numeric(Nsim) )
start_time <- Sys.time()
while(Nsim){
  tree_s<-rcoal(ntip) #from coalescent with Ne=1, for more general Ne trajectories, use https://github.com/JuliaPalacios/phylodyn/blob/master/vignettes/Simulation.Rmd
  times_s<-coalescent.intervals(tree_s)
  simulation<-cumsum(times_s$interval.length)
  if(simulation[ntip-1]<tau){
    coal_time3[Nsim,]=simulation
    Nsim=Nsim-1
  }
}
end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken <- as.numeric(time_taken, units = "secs")
time_taken/Nsim_ori
nam <- paste("naivesyn1data_ntip",ntip, "_samplesize", Nsim_ori,"_tau", tau, sep = "")
save(coal_time3, file=file.path(output_dir, paste(nam,'.rda',sep='')))


##################Thinning#########################################################################
##################syn3#############################################################################
########################ntip=100###################################################################
########################tau=0.71####################################################################
ntip=100
samp_times = c(0)
n_sampled = c(ntip)
tau=0.71
Nsim_ori=3000
Nsim=3000
coal_time4=replicate(ntip-1, numeric(Nsim) )


r_func <- function(k, j) {
  if(j == 1) return(1)
  prod <- 1
  for(m in 1:(j-1)) {
    prod <- prod * ((2*m + 1)/(2*m - 1)) * ((k - m)/(k + m))
  }
  return((-1)^(j-1) * prod)
}

inten2 <- function(t) {
  25 * exp(-5 * t)
}
inten2_inv<- function(t) {
  exp(5 * t)/25
}

inten2_inv_cum<- function(t) {
  (exp(5*t)-1)/125
}


inten2_inv_cum_inv<- function(t) {
  log(125*t+1)/5
}

result_list <- lapply(seq_len(ntip), function(k) {
  sapply(seq_len(k), function(i) r_func(k, i))
})

com_vec <- choose(1:ntip, 2)
start_time <- Sys.time()
while(Nsim){
  simulation4<-rep(0,ntip-1)
  t_upper<-0
  for (j in ntip:2){
    while(1){
      ori<-rexp(n=1, rate=choose(j,2))
      t_upper<-inten2_inv_cum_inv(inten2_inv_cum(tau)+ori-log(exp(inten2_inv_cum(tau)-inten2_inv_cum(t_upper))-1+exp(ori)))
      kk=j
      con=inten2_inv_cum(t_upper)-inten2_inv_cum(tau) #gamma(t)-gamma(\tau)
      denom=sum(result_list[[kk]]*exp(com_vec[1:(kk)]*con))
      numer=sum(result_list[[kk-1]]*exp(com_vec[1:(kk-1)]*con))

      if(stats::runif(1) <= numer/denom*(1-exp(con))) {
        simulation4[ntip+1-j]<-t_upper
        break
      }
    }
  }
  coal_time4[Nsim,]=simulation4
  Nsim=Nsim-1
}
end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken <- as.numeric(time_taken, units = "secs")
time_taken/Nsim_ori

nam <- paste("thinsyn3data_ntip",ntip, "_samplesize", Nsim_ori,"_tau", tau, sep = "")
save(coal_time4, file=file.path(output_dir, paste(nam,'.rda',sep='')))

##################Naive rejection#########################################################################
##################syn3#############################################################################
########################ntip=50###################################################################
########################tau=0.71####################################################################
ntip=100
n_sampled = c(ntip)
tau=0.71
Nsim_ori=3000
Nsim=3000
coal_time5=replicate(ntip-1, numeric(Nsim) )
#require(combinat)
com=rep(0,ntip-1)
for (k in 1:(ntip-1)){
  com[k]=choose(ntip+1-k,2)}

    #dim(combn(ntip+1-k,2))[2]}
com[ntip-1]=1
start_time <- Sys.time()
while(Nsim){
  simulation3<-rep(0,(ntip-1))
  simulation4<-rep(0,ntip)
  for (j in 1:(ntip-1)){
    simulation3[j]<-rexp(n=1, rate=com[j])
    simulation4[j+1]<-log(125*simulation3[j]+exp(5*simulation4[j]))/5

  }
  if(simulation4[ntip]<tau){
    coal_time5[Nsim,]=simulation4[2:ntip]
    Nsim=Nsim-1
  }
}
end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken <- as.numeric(time_taken, units = "secs")
time_taken/Nsim_ori
nam <- paste("naivesyn3data_ntip",ntip, "_samplesize", Nsim_ori,"_tau", tau, sep = "")
save(coal_time5, file=file.path(output_dir, paste(nam,'.rda',sep='')))

###############################comparison of histogram and mean&std###################################
#################################################histogram comparison
#histogram of the median tip
jj=50
par(mfrow = c(2,1))
#coal_time2 is syn1_tau0.5_ntip100_thinning
hist(coal_time2[,jj],
     breaks = 20,
     col = rgb(1, 0, 0, 0.5), # semi-transparent red
     main = expression("Thinning: a histogram of " ~ T[51]^B),
     xlab = "Time",
     freq = FALSE)
#hist(coal_time1[,jj],
#     breaks = 20,
#     col = rgb(0, 0, 1, 0.5), # semi-transparent blue
#     main = expression("Carson: a histogram of " ~ T[51]^B),
#     xlab = "Time",
#     freq = FALSE           # match freq=FALSE or =TRUE from first
#)
#coal_time3 is syn1_tau0.5_ntip100_naiverejection
hist(coal_time3[,jj],
     breaks = 20,
     col = rgb(0, 1, 0, 0.5), # semi-transparent blue
     main = expression("Naive rejection: a histogram of " ~ T[51]^B),
     xlab = "Time",
     freq = FALSE           # match freq=FALSE or =TRUE from first
)
#mean(coal_time1[,jj])
mean(coal_time2[,jj])
mean(coal_time3[,jj])
#sd(coal_time1[,jj])
sd(coal_time2[,jj])
sd(coal_time3[,jj])




