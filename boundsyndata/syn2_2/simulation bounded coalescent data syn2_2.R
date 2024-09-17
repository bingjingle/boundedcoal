library("phylodyn")

#####################################################################################################
#Ne(t)=25exp(-5t)+1
#Ne(t)=scale*exp(-rate*t)
bottleneck_traj<-function(t){
  result=exp_traj(t,scale = 25, rate = 5)+1
  return(result)
}
#################################transformation#####################################################################
ntip=20
samp_times = c(0)
n_sampled = c(ntip)
tau=1
Nsim=10
coal_time=replicate(ntip-1, numeric(Nsim) ) 
while(Nsim){
  simulation3<-coalsim(samp_times = samp_times, n_sampled = n_sampled, traj = bottleneck_traj)
  if(simulation3$coal_times[ntip-1]<tau){
    coal_time[Nsim,]=simulation3$coal_times
    Nsim=Nsim-1
  }
  else{print('fail')}
}

getwd()
setwd("C:/Users/bingj/OneDrive - Stanford/bounded coalescent/syndata")
nam <- paste("syn2data_",2, sep = "")
save(coal_time, file=paste(nam,'.rda',sep=''))

