args = commandArgs(trailingOnly=TRUE)

library(parallel)
library(DEoptim)
library(pracma)
library(TruncatedNormal)
library(matrixcalc) # is.positive.definite
#library(matlib)
library(Matrix)
library(base)
library(mvtnorm)
library(truncnorm)

points_inhomo=c(0.002404551, 0.003284109, 0.014545722, 0.015916560, 0.016345294, 0.024471511, 0.028874182, 0.031693528, 0.043305040, 0.045027608, 
                0.048955524, 0.049016560, 0.051282521, 0.051827613, 0.058178058, 0.059350846, 0.061520275, 0.063898512, 0.064510277, 0.073763681, 
                0.074190927, 0.075921131, 0.084351262, 0.123951058, 0.127376496, 0.130710786, 0.131076997, 0.134619519, 0.135516560, 0.137263625, 
                0.146624658, 0.176979844, 0.177407090, 0.179539485, 0.181650978, 0.189057353, 0.190807508, 0.194124995, 0.210141724, 0.210507935, 
                0.213908250, 0.224499077, 0.228386124, 0.232965850, 0.235036115, 0.238789020, 0.243985010, 0.248129450, 0.257160114, 0.262849516, 
                0.272764811, 0.273868512, 0.277445977, 0.279938943, 0.282193146, 0.284454867, 0.284515902, 0.305034539, 0.310960039, 0.327718416, 
                0.328885495, 0.348102452, 0.350663863, 0.352375918, 0.357017910, 0.358052490, 0.362508677, 0.363135938, 0.372395330, 0.377740297, 
                0.378707698, 0.395722036, 0.401328539, 0.403371599, 0.427609114, 0.450158221, 0.464668061, 0.488975108, 0.494816686, 0.499820467, 
                0.506044204, 0.516651251, 0.522448434, 0.537757906, 0.557675146, 0.578180530, 0.581644027, 0.588355957, 0.603176719, 0.610313729, 
                0.639500622, 0.652845651, 0.702422905, 0.760440888, 0.799546757, 0.910557807, 0.959638933, 1.058403095, 1.244499064)

days=c(0,points_inhomo)


N=length(points_inhomo)
com=rep(0,N)
for (k in 1:(N-1)){
  print(k)
  print(dim(combn(N+2-k,2))[2])
  com[k]=dim(combn(N+2-k,2))[2]}
com[N]=1


T=max(points_inhomo)


expo_quad_kernel<-function(theta00,theta11,xn,xm){ # 1,0.1
  return(theta00*exp(-theta11/2*sum((xn - xm)**2)))
}

expo_quad_kernel2<-function(theta00,theta11,xn,t1,t2){ # 1,0.1
  return(sqrt(pi/2/theta11)*theta00*(erf(sqrt(theta11/2)*(t2-xn))-erf(sqrt(theta11/2)*(t1-xn))))
}

expo_quad_kernel3<-function(theta00,theta11,t1,t2,t3,t4){ # 1,0.1
  return(sqrt(pi/2/theta11)*theta00*(    (t2-t3)*erf(  sqrt(theta11/2) *(t2-t3)    )- (t2-t4)*erf(sqrt(theta11/2) *(t2-t4)    )
                                         -(t1-t3)*erf(  sqrt(theta11/2) *(t1-t3)        )+  (t1-t4)*erf(sqrt(theta11/2) *(t1-t4)    ))+theta00/theta11*( exp(  -theta11/2*(t2-t3)**2    ) -exp(  -theta11/2*(t1-t3)**2    )  - exp(  -theta11/2*(t2-t4)**2    ) 
                                                                                                                                                         +exp(  -theta11/2*(t1-t4)**2    ))
  )
  
}

ccc=100
mm=as.numeric(args[1])
delta_m=T/(mm-1)
t=seq(0,T,length.out=mm)
negloglikelihood31<-function(samps,par){
  theta00=par[1]
  theta11=par[2]
  m=length(samps)
  t=seq(0,T,length.out=m)
  delta_m=T/(m-1)
  c_vec=rep(delta_m,m)
  c_vec[1]=delta_m/2
  c_vec[m]=delta_m/2
  #comp1=-sum(c_vec*samps)
  comp1=c(0,cumsum(c_vec*samps))
  comp2=c(0,cumsum(c_vec))
  samp_event=rep(0,N)
  fin<-function(x){
    index_L=floor(x/delta_m)+1
    index_R=ceiling(x/delta_m)+1
    if (x-t[index_L]<=delta_m/2){
      fin=index_L
    }
    else{
      fin=index_R
    }
    return(samps[fin])
  }
  
  samp_event=sapply(points_inhomo,fin)
  
  
  samp_event_idx=rep(0,N)
  
  fin2<-function(x){
    index_L=floor(x/delta_m)+1
    index_R=ceiling(x/delta_m)+1
    if (x-t[index_L]<=delta_m/2){
      fin=index_L
    }
    else{
      fin=index_R
    }
    return(fin)
  }
  
  
  
  
  samp_event_idx=sapply(points_inhomo,fin2)
  comp3=comp1[samp_event_idx]+samp_event*(points_inhomo-comp2[samp_event_idx])#double check later
  comp4=diff(c(0, comp3))
  
  comp_total=sum(log(samp_event))-sum(com*comp4)#loglikelihood
  
  
  mmm=m+N
  cov=matrix(0,mmm,mmm)
  sampinteg=c(samps,comp4)
  for (i in 1:mmm)
  {
    for (j in i:mmm){
      if (j<=m & i<=m){
        cov[i,j]=expo_quad_kernel(theta00,theta11,t[i],t[j])
      }
      if(j>m&i<=m){
        cov[i,j]=expo_quad_kernel2(theta00,theta11,t[i],days[j-m],days[j-m+1])
      }
      if(i>m){
        cov[i,j]=expo_quad_kernel3(theta00,theta11,days[i-m],days[i-m+1],days[j-m],days[j-m+1])  
      }
      if(j!=i){
        cov[j,i]=cov[i,j]
      }
    }
  }
  cov_noise=cov+diag(mmm)*1e-5
  finalvalue=dmvnorm(sampinteg, mean=rep(0,mmm),sigma = cov_noise,log=TRUE)-log(pmvnorm(lower=rep(0,mmm),
                                                                                        upper=rep(Inf,mmm),
                                                                                        mean=rep(0,mmm),
                                                                                        sigma = cov_noise)[1])
  
  if(is.positive.definite(cov_noise)==FALSE){finalvalue=-1e10}
  #print('comp_total')
  #print(-comp_total)
  #print('finalvalue')
  #print(-finalvalue/ccc)
  return(-comp_total-finalvalue/ccc)
}

hypeasy<-function(par){
  return(negloglikelihood31(par[-c(1,2)],par[c(1,2)]))
}


outnew <- DEoptim(fn=hypeasy, lower = rep(0,2+mm), upper = c(1e15,1000,rep(200,mm)),control=list(trace=TRUE,itermax=10))#,fnMap=Mapfun)

b=c(outnew$optim$bestval,outnew$optim$bestmem)
save(b, file =paste0( "/home/groups/juliapr/Bingjing/code/boundedcoal/syndata/hyperpara/stdsyn2data_1/ccc100_mm",args[1],".rda"))