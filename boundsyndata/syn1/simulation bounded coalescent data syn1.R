library("phylodyn")

###########################syn1#####################################################
constant<-function(x){
  return (rep(1,length(x)))
}
ntip=20

tau=1

Nsim=10
coal_time=replicate(ntip-1, numeric(Nsim) ) 
while(Nsim){
  simulation1<-coalsim(0,ntip,constant)
  if(simulation1$coal_times[ntip-1]<tau){
    coal_time[Nsim,]=simulation1$coal_times
    Nsim=Nsim-1
  }
  else{print('fail')}
}

setwd("C:/Users/bingj/OneDrive - Stanford/bounded coalescent/syndata")
nam <- paste("syn1data_",1, sep = "")
save(coal_time, file=paste(nam,'.rda',sep=''))
