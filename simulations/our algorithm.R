##################Thinning#########################################################################
##################syn3#############################################################################
########################ntip=100###################################################################
########################tau=0.71####################################################################
####################user defined#####################################
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


ntip=100
tau=0.71
Nsim=3000
####################################################################

sim_bc<-function(inten2,inten2_inv,inten2_inv_cum,inten2_inv_cum_inv,ntip,tau,Nsim){

  n_sampled = c(ntip)
  coal_time1=replicate(ntip-1, numeric(Nsim) ) 

  r_func <- function(k,j) {
    if(j == 1) return(1)
    prod <- 1
    for(m in 1:(j-1)) {
      prod <- prod * ((2*m + 1)/(2*m - 1)) * ((k - m)/(k + m))
    }
    return((-1)^(j-1) * prod)
  }

  result_list <- lapply(seq_len(ntip), function(k) {
    sapply(seq_len(k), function(i) r_func(k, i))
  })

  com_vec <- choose(1:ntip, 2)
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
    coal_time1[Nsim,]=simulation4
    Nsim=Nsim-1
  }
  return(coal_time1)
}

sim_bc(inten2,inten2_inv,inten2_inv_cum,inten2_inv_cum_inv,ntip,tau,Nsim)
