#Synethic datasets generation

library("phylodyn")
library("ape")
setwd("C:/Users/bingj/OneDrive - Stanford/bounded coalescent/syndatafinal")
set.seed(123)

##################################################################lambda_1 N_e(t)=1##########################################################################
#############30 synthetic datasets generating from standard coalescent model

# ntip=100
ntip=100
Nsim=30
coal_time=replicate(ntip-1, numeric(Nsim) ) 
while(Nsim){
  #simulation1<-coalsim(0,ntip,constant)
  tree_s<-rcoal(ntip) #from coalescent with Ne=1, for more general Ne trajectories, use https://github.com/JuliaPalacios/phylodyn/blob/master/vignettes/Simulation.Rmd
  times_s<-coalescent.intervals(tree_s)
  simulation<-cumsum(times_s$interval.length)
  coal_time[Nsim,]=simulation
  Nsim=Nsim-1
  
}
nam <- paste("syn1_std_ntip",ntip, sep = "")
save(coal_time, file=paste(nam,'.rda',sep=''))

#ntip=50
ntip=50
Nsim=30
coal_time=replicate(ntip-1, numeric(Nsim) ) 
while(Nsim){
  #simulation1<-coalsim(0,ntip,constant)
  tree_s<-rcoal(ntip) #from coalescent with Ne=1, for more general Ne trajectories, use https://github.com/JuliaPalacios/phylodyn/blob/master/vignettes/Simulation.Rmd
  times_s<-coalescent.intervals(tree_s)
  simulation<-cumsum(times_s$interval.length)
  coal_time[Nsim,]=simulation
  Nsim=Nsim-1
  
}
nam <- paste("syn1_std_ntip",ntip, sep = "")
save(coal_time, file=paste(nam,'.rda',sep=''))

#################30 synthetic datasets generating from bounded coalescent model, \tau=1
#ntip=100
ntip=100
tau=1
Nsim=30
coal_time=replicate(ntip-1, numeric(Nsim) ) 
while(Nsim){
  tree_s<-rcoal(ntip) #from coalescent with Ne=1, for more general Ne trajectories, use https://github.com/JuliaPalacios/phylodyn/blob/master/vignettes/Simulation.Rmd
  times_s<-coalescent.intervals(tree_s)
  simulation<-cumsum(times_s$interval.length)
  if(simulation[ntip-1]<tau){
    coal_time[Nsim,]=simulation
    Nsim=Nsim-1
  }
  else{print('fail')}
}
nam <- paste("syn1_bounded_ntip",ntip, sep = "")
save(coal_time, file=paste(nam,'.rda',sep=''))


#ntip=50
ntip=50
tau=1
Nsim=30
coal_time=replicate(ntip-1, numeric(Nsim) ) 
while(Nsim){
  tree_s<-rcoal(ntip) #from coalescent with Ne=1, for more general Ne trajectories, use https://github.com/JuliaPalacios/phylodyn/blob/master/vignettes/Simulation.Rmd
  times_s<-coalescent.intervals(tree_s)
  simulation<-cumsum(times_s$interval.length)
  if(simulation[ntip-1]<tau){
    coal_time[Nsim,]=simulation
    Nsim=Nsim-1
  }
  else{print('fail')}
}
nam <- paste("syn1_bounded_ntip",ntip, sep = "")
save(coal_time, file=paste(nam,'.rda',sep=''))


##################################################################lambda_2: N_e(t) = 3exp{−t} ##########################################################################
#############30 synthetic datasets generating from standard coalescent model
#ntip=100
ntip=100
Nsim=30
coal_time=replicate(ntip-1, numeric(Nsim) ) 

com=rep(0,ntip-1)
for (k in 1:(ntip-1)){
  com[k]=choose(ntip+1-k,2)}


while(Nsim){
  simulation3<-rep(0,(ntip-1))
  simulation4<-rep(0,ntip)
  for (j in 1:(ntip-1)){
    simulation3[j]<-rexp(n=1, rate=com[j])
    simulation4[j+1]<-log(3*simulation3[j]+exp(simulation4[j]))
  }
  coal_time[Nsim,]=simulation4[2:ntip]
  Nsim=Nsim-1
}
nam <- paste("syn2_std_ntip",ntip, sep = "")
save(coal_time, file=paste(nam,'.rda',sep=''))



#ntip=50
ntip=50
Nsim=30
coal_time=replicate(ntip-1, numeric(Nsim) ) 

com=rep(0,ntip-1)
for (k in 1:(ntip-1)){
  com[k]=choose(ntip+1-k,2)}


while(Nsim){
  simulation3<-rep(0,(ntip-1))
  simulation4<-rep(0,ntip)
  for (j in 1:(ntip-1)){
    simulation3[j]<-rexp(n=1, rate=com[j])
    simulation4[j+1]<-log(3*simulation3[j]+exp(simulation4[j]))
  }
  coal_time[Nsim,]=simulation4[2:ntip]
  Nsim=Nsim-1
}
nam <- paste("syn2_std_ntip",ntip, sep = "")
save(coal_time, file=paste(nam,'.rda',sep=''))


#################30 synthetic datasets generating from bounded coalescent model, \tau=0.7
#ntip=100
ntip=100
tau=.7
Nsim=30
coal_time=replicate(ntip-1, numeric(Nsim) ) 

com=rep(0,ntip-1)
for (k in 1:(ntip-1)){
  com[k]=choose(ntip+1-k,2)}

while(Nsim){
  simulation3<-rep(0,(ntip-1))
  simulation4<-rep(0,ntip)
  for (j in 1:(ntip-1)){
    simulation3[j]<-rexp(n=1, rate=com[j])
    simulation4[j+1]<-log(3*simulation3[j]+exp(simulation4[j]))
    
  }
  if(simulation4[ntip]<tau){
    coal_time[Nsim,]=simulation4[2:ntip]
    Nsim=Nsim-1
  }
  else{print('fail')}
}

nam <- paste("syn2_bounded_ntip",ntip, sep = "")
save(coal_time, file=paste(nam,'.rda',sep=''))


#ntip=50
ntip=50
tau=.7
Nsim=30
coal_time=replicate(ntip-1, numeric(Nsim) ) 

com=rep(0,ntip-1)
for (k in 1:(ntip-1)){
  com[k]=choose(ntip+1-k,2)}

while(Nsim){
  simulation3<-rep(0,(ntip-1))
  simulation4<-rep(0,ntip)
  for (j in 1:(ntip-1)){
    simulation3[j]<-rexp(n=1, rate=com[j])
    simulation4[j+1]<-log(3*simulation3[j]+exp(simulation4[j]))
    
  }
  if(simulation4[ntip]<tau){
    coal_time[Nsim,]=simulation4[2:ntip]
    Nsim=Nsim-1
  }
  else{print('fail')}
}

nam <- paste("syn2_bounded_ntip",ntip, sep = "")
save(coal_time, file=paste(nam,'.rda',sep=''))



##################################################################lambda_3: N_e(t) = exp{−5t} ##########################################################################
#############30 synthetic datasets generating from standard coalescent model
#ntip=100
ntip=100
Nsim=30
coal_time=replicate(ntip-1, numeric(Nsim) ) 

com=rep(0,ntip-1)
for (k in 1:(ntip-1)){
  com[k]=choose(ntip+1-k,2)}

while(Nsim){
  simulation3<-rep(0,(ntip-1))
  simulation4<-rep(0,ntip)
  for (j in 1:(ntip-1)){
    simulation3[j]<-rexp(n=1, rate=com[j])
    simulation4[j+1]<-log(125*simulation3[j]+exp(5*simulation4[j]))/5
  }
  coal_time[Nsim,]=simulation4[2:ntip]
  Nsim=Nsim-1
}
nam <- paste("syn3_std_ntip",ntip, sep = "")
save(coal_time, file=paste(nam,'.rda',sep=''))


#ntip=50
ntip=50
Nsim=30
coal_time=replicate(ntip-1, numeric(Nsim) ) 

com=rep(0,ntip-1)
for (k in 1:(ntip-1)){
  com[k]=choose(ntip+1-k,2)}

while(Nsim){
  simulation3<-rep(0,(ntip-1))
  simulation4<-rep(0,ntip)
  for (j in 1:(ntip-1)){
    simulation3[j]<-rexp(n=1, rate=com[j])
    simulation4[j+1]<-log(125*simulation3[j]+exp(5*simulation4[j]))/5
  }
  coal_time[Nsim,]=simulation4[2:ntip]
  Nsim=Nsim-1
}
nam <- paste("syn3_std_ntip",ntip, sep = "")
save(coal_time, file=paste(nam,'.rda',sep=''))


#################30 synthetic datasets generating from bounded coalescent model, \tau=0.71
#ntip=100
ntip=100
tau=0.71
Nsim=30
coal_time=replicate(ntip-1, numeric(Nsim) ) 

com=rep(0,ntip-1)
for (k in 1:(ntip-1)){
  com[k]=choose(ntip+1-k,2)}

while(Nsim){
  simulation3<-rep(0,(ntip-1))
  simulation4<-rep(0,ntip)
  for (j in 1:(ntip-1)){
    simulation3[j]<-rexp(n=1, rate=com[j])
    simulation4[j+1]<-log(125*simulation3[j]+exp(5*simulation4[j]))/5
    
  }
  if(simulation4[ntip]<tau){
    coal_time[Nsim,]=simulation4[2:ntip]
    Nsim=Nsim-1
  }
  else{print('fail')}
}

nam <- paste("syn3_bounded_ntip",ntip, sep = "")
save(coal_time, file=paste(nam,'.rda',sep=''))

#ntip=50
ntip=50
tau=0.71
Nsim=30
coal_time=replicate(ntip-1, numeric(Nsim) ) 

com=rep(0,ntip-1)
for (k in 1:(ntip-1)){
  com[k]=choose(ntip+1-k,2)}

while(Nsim){
  simulation3<-rep(0,(ntip-1))
  simulation4<-rep(0,ntip)
  for (j in 1:(ntip-1)){
    simulation3[j]<-rexp(n=1, rate=com[j])
    simulation4[j+1]<-log(125*simulation3[j]+exp(5*simulation4[j]))/5
    
  }
  if(simulation4[ntip]<tau){
    coal_time[Nsim,]=simulation4[2:ntip]
    Nsim=Nsim-1
  }
  else{print('fail')}
}
coal_time[30,]

nam <- paste("syn3_bounded_ntip",ntip, sep = "")
save(coal_time, file=paste(nam,'.rda',sep=''))

