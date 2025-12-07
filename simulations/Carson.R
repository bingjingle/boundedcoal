library(BoundedCoalescent)
library(ape)
Nsim=3000
ntip=100
tau=0.5
start_time <- Sys.time()
coal_time2=replicate(ntip-1, numeric(Nsim) ) 
while(Nsim){
  
  t=bounded_sample_phylo(c(tau),c(ntip),1,0)
  coal_time2[Nsim,]=cumsum(coalescent.intervals(t$phylo)$interval.length)
  Nsim=Nsim-1
}
end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken/Nsim

setwd("C:/Users/bingj/OneDrive - Stanford/bounded coalescent/simulationcomparisonnew")
nam <- paste("carsonsyn1data_ntip", ntip, "_samplesize", Nsim, sep = "")
save(coal_time2, file=paste(nam,'.rda',sep=''))

