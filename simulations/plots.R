##Plots for Bounded coalescent paper
##Experiment with constant pop. size
##Validation that the new function with large bound is equivalent to standard coalescent

constant<-function(x){
  return (rep(1,length(x)))
}

iters<-1000
nsamp<-10
bound<-0.5
result<-matrix(0,nrow=iters,ncol=nsamp-1)
for (j in 1:iters){
  cc<-coalsim_bounded(c(0),c(nsamp),constant,bound=bound)
  result[j,]<-cc$coal_times
}


result2<-matrix(0,nrow=iters,ncol=nsamp-1)
for (j in 1:iters){
  cc<-phylodyn:::coalsim_tt(c(0),c(nsamp),constant)
  result2[j,]<-cc$coal_times
}

par(mfrow=c(1,2))
idx <- seq(10,2)
labs <- as.expression(lapply(idx, function(i) bquote(t[.(i)])))
boxplot(result,2,names=labs,ylim=c(0,10),main="Bounded Coalescent")
boxplot(result2,2,names=labs,ylim=c(0,10),main="Coalescent")



library(ggplot2)


overlay_boxplots <- function(X, Y,
                             name_x = "BC", name_y = "C",
                             blue = "#1f77b4", red = "#d62728") {
  stopifnot(is.matrix(X), is.matrix(Y), ncol(X) == ncol(Y))
  n <- ncol(X)
  
  dfX <- data.frame(
    value = as.vector(X),
    col   = factor(rep(seq_len(n), each = nrow(X)), levels = seq_len(n))
  )
  dfY <- data.frame(
    value = as.vector(Y),
    col   = factor(rep(seq_len(n), each = nrow(Y)), levels = seq_len(n))
  )
  
  ggplot() +
    # map colour to a constant label so ggplot creates a legend
    geom_boxplot(data = dfX, aes(x = col, y = value, colour = name_x),
                 fill = blue, alpha = 0.25, width = 0.65, outlier.alpha = 0.25) +
    geom_boxplot(data = dfY, aes(x = col, y = value, colour = name_y),
                 fill = NA, linewidth = 1.2, width = 0.45,
                 outlier.colour = red, outlier.alpha = 0.5) +
    scale_x_discrete(labels = parse(text = paste0("t[", seq(n+1,2), "]"))) +
    scale_colour_manual(values = c(setNames(blue, name_x),
                                   setNames(red,  name_y)),
                        name = NULL) +
    theme_classic() +
    theme( axis.text.x = element_text(size = 14))+
    theme(
      panel.border = element_rect(colour = "black", fill = NA),
      legend.position   = c(0.1, 0.8),
      legend.background = element_rect(fill = "white", colour = "black")
    ) +
    labs(title ="Ne=1, Bound=0.5", x = "Coalescent Time", y = NULL)
}


plot1<-overlay_boxplots(result,result2)




#Experiment with exponential growth
exp_traj = function(t, scale=1000, rate=1)
{
  return(scale * exp(-t*rate))
}

iters<-1000
nsamp<-10
bound<-0.5
result<-matrix(0,nrow=iters,ncol=nsamp-1)
for (j in 1:iters){
  cc<-coalsim_bounded(c(0),c(nsamp),exp_traj,bound=bound,scale=25,rate=5)
  result[j,]<-cc$coal_times
}


result2<-matrix(0,nrow=iters,ncol=nsamp-1)
for (j in 1:iters){
  cc<-phylodyn:::coalsim_tt(c(0),c(nsamp),exp_traj,scale=25,rate=5)
  result2[j,]<-cc$coal_times
}

par(mfrow=c(1,2))
idx <- seq(10,2)
labs <- as.expression(lapply(idx, function(i) bquote(t[.(i)])))
boxplot(result,2,names=labs,ylim=c(0,2),main="Bounded Coalescent")
boxplot(result2,2,names=labs,ylim=c(0,2),main="Coalescent")



library(ggplot2)


overlay_boxplots <- function(X, Y,
                             name_x = "BC", name_y = "C",
                             blue = "#1f77b4", red = "#d62728") {
  stopifnot(is.matrix(X), is.matrix(Y), ncol(X) == ncol(Y))
  n <- ncol(X)
  
  dfX <- data.frame(
    value = as.vector(X),
    col   = factor(rep(seq_len(n), each = nrow(X)), levels = seq_len(n))
  )
  dfY <- data.frame(
    value = as.vector(Y),
    col   = factor(rep(seq_len(n), each = nrow(Y)), levels = seq_len(n))
  )
  
  ggplot() +
    # map colour to a constant label so ggplot creates a legend
    geom_boxplot(data = dfX, aes(x = col, y = value, colour = name_x),
                 fill = blue, alpha = 0.25, width = 0.65, outlier.alpha = 0.25) +
    geom_boxplot(data = dfY, aes(x = col, y = value, colour = name_y),
                 fill = NA, linewidth = 1.2, width = 0.45,
                 outlier.colour = red, outlier.alpha = 0.5) +
    scale_x_discrete(labels = parse(text = paste0("t[", seq(n+1,2), "]"))) +
    scale_colour_manual(values = c(setNames(blue, name_x),
                                   setNames(red,  name_y)),
                        name = NULL) +
    theme_classic() +
    theme( axis.text.x = element_text(size = 14))+
    theme(
      panel.border = element_rect(colour = "black", fill = NA),
      legend.position   = c(0.1, 0.8),
      legend.background = element_rect(fill = "white", colour = "black")
    ) +
    labs(title = "Ne=25 exp(-5t), Bound=0.5", x = "Coalescent Time", y = NULL)
}


plot2<-overlay_boxplots(result,result2)

# Comparison with rejection sampling t_{3} Ne=1 and Ne=exponential


constant<-function(x){
  return (rep(1,length(x)))
}

iters<-1000
nsamp<-10
bound<-0.5
result<-matrix(0,nrow=iters,ncol=nsamp-1)
for (j in 1:iters){
  cc<-coalsim_bounded(c(0),c(nsamp),constant,bound=bound)
  result[j,]<-cc$coal_times
}

iters<-1000
nsamp<-10
bound<-0.5

#rejection sampling
library("ape")
result3<-result
i<-1
iter<-0
while(i<=iters){
  iter<-iter+1
  tree_s<-phylodyn:::coalsim_tt(c(0),nsamp,constant)
  simulation<-tree_s$coal_times
  if(simulation[length(simulation)]<bound){
    result3[i,]<-simulation
    i<-i+1
    print(i)
  }
}
print(i/iter) ##Acceptance probability
mean(result3[,8])
hist(result3[,8])

overlay_hist <- function(x, y,
                         name_x = "X", name_y = "Y",
                         blue = "#1f77b4", red = "#d62728",
                         bins = 30) {
  df <- rbind(
    data.frame(value = x, group = name_x),
    data.frame(value = y, group = name_y)
  )
  
  ggplot(df, aes(x = value, fill = group, colour = group)) +
    geom_histogram(position = "identity", bins = bins, alpha = 0.35) +
    scale_fill_manual(values = c(setNames(blue, name_x), setNames(red, name_y)),
                      name = NULL) +
    scale_colour_manual(values = c(setNames(blue, name_x), setNames(red, name_y)),
                        name = NULL) +
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA),
      legend.position = c(0.2, 0.8),  # inside panel (0-1, 0-1)
      legend.background = element_rect(fill = "white", colour = "black")
    ) +
    labs(title = "Validation: Algorithm1 vs Rejection, Ne=1, Bound=0.5 ",
         x = expression(T[3]), y = "Count")
}
plot3<-overlay_hist(result[,8],result3[,8],name_x="Algorithm 1",name_y="Rejection")

#ks.test(result[,8],result3[,8])
library("covalchemy")
calculate_tv_distance_empirical(result[,8], result3[,8])

constant<-function(x){
  return (rep(1,length(x)))
}

distance<-rep(0,100)
result_t1<-matrix(0,nrow=1000,ncol=9)
result_t2<-matrix(0,nrow=1000,ncol=9)
for (i in 1:100){
  
  iters<-1000
  nsamp<-10
  bound<-0.5
  for (j in 1:iters){
    cc<-coalsim_bounded(c(0),c(nsamp),constant,bound=bound)
    cc2<-coalsim_bounded(c(0),c(nsamp),constant,bound=bound)
    result_t1[j,]<-cc$coal_times
    result_t2[j,]<-cc2$coal_times
  }
  distance[i]<-calculate_tv_distance_empirical(result_t1[,8], result_t2[,8])

}
mean(distance)

calculate_tv_distance_empirical(resulta[,8], resultb[,8])

mean(result[,8])
sd(result[,8])
mean(result3[,8])
sd(result3[,8])

# Comparison with rejection sampling t_{3} Ne=exponential


iters<-1000
nsamp<-10
bound<-0.9
result<-matrix(0,nrow=iters,ncol=nsamp-1)
for (j in 1:iters){
  cc<-coalsim_bounded(c(0),c(nsamp),exp_traj,bound=bound,scale=25,rate=5)
  result[j,]<-cc$coal_times
}

mean(result[,8])
sd(result[,8])

library("ape")
result4<-matrix(0,nrow=iters,ncol=nsamp-1)
i<-1
iter<-0
iters<-1000
bound<-0.9
while(i<=iters){
  iter<-iter+1
  tree_s<-phylodyn:::coalsim_tt(c(0),nsamp,exp_traj,scale=25,rate=5)
  simulation<-tree_s$coal_times
  if(simulation[length(simulation)]<bound){
    result4[i,]<-simulation
    i<-i+1
    print(i)
  }
}
print(i/iter) ##Acceptance probability
mean(result4[,8])
hist(result4[,8])

calculate_tv_distance_empirical(result[,8], result4[,8])


distance<-rep(0,100)
result_t1<-matrix(0,nrow=1000,ncol=9)
result_t2<-matrix(0,nrow=1000,ncol=9)
for (i in 1:100){
  
  iters<-1000
  nsamp<-10
  bound<-0.9
  for (j in 1:iters){
    cc<-coalsim_bounded(c(0),c(nsamp),exp_traj,bound=bound,scale=25,rate=5)
    cc2<-coalsim_bounded(c(0),c(nsamp),exp_traj,bound=bound,scale=25,rate=5)
    result_t1[j,]<-cc$coal_times
    result_t2[j,]<-cc2$coal_times
  }
  distance[i]<-calculate_tv_distance_empirical(result_t1[,8], result_t2[,8])
  
}
mean(distance)



overlay_hist <- function(x, y,
                         name_x = "X", name_y = "Y",
                         blue = "#1f77b4", red = "#d62728",
                         bins = 30) {
  df <- rbind(
    data.frame(value = x, group = name_x),
    data.frame(value = y, group = name_y)
  )
  
  ggplot(df, aes(x = value, fill = group, colour = group)) +
    geom_histogram(position = "identity", bins = bins, alpha = 0.35) +
    scale_fill_manual(values = c(setNames(blue, name_x), setNames(red, name_y)),
                      name = NULL) +
    scale_colour_manual(values = c(setNames(blue, name_x), setNames(red, name_y)),
                        name = NULL) +
    theme_classic() +
    theme( axis.text.x = element_text(size = 14))+
    theme(
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA)
    ) +
    labs(title = "Algorithm 1 vs Rejection, Ne=25 exp(-5t), Bound=0.9 ",
         x = expression(T[3]), y = "Count")
}
plot4<-overlay_hist(result[,8],result4[,8],name_x="BC",name_y="Rejection")
ks.test(result[,8], result4[,8])

mean(result[,8])
sd(result[,8])
mean(result4[,8])
sd(result4[,8])

library(cowplot)
plot_grid(plot1, plot2, plot3, plot4, ncol = 2)

plot_grid(plot1,plot2,ncol=2)

plot_grid(plot1, plot2, plot3, plot4,
          ncol = 2,
          labels = c("(A)", "(B)", "(C)", "(D)"),
          label_size = 12)

