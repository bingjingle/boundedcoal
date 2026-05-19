ntip <- 100
bound=1.5
Ngrid=1000
noise_var <- 1e-3#jitter
com_vec <- choose(seq_len(ntip), 2)

a_coeffs_kmax <- function(kmax) {
  if (kmax == 1 || kmax == 2) {
    return(c(1))
  }
  
  a_prev <- c(1)  # k = 2
  
  for (k in 3:kmax) {
    Mk <- choose(k - 1, 2)
    Mk_prev <- choose(k - 2, 2)
    
    a_cur <- numeric(Mk + 1)
    
    for (i in 0:Mk) {
      # R index is i + 1
      
      if (i == 0) {
        a_cur[i + 1] <- 1
        next
      }
      
      a_im1_k <- a_cur[i]  # a_{i-1,k}
      
      if (i <= Mk_prev) {
        a_i_km1 <- a_prev[i + 1]  # a_{i,k-1}
      } else {
        a_i_km1 <- 0
      }
      
      numerator <- 
        (choose(k - 1, 2) - i + 1) * a_im1_k +
        choose(k, 2) * a_i_km1
      
      denominator <- choose(k, 2) - i
      
      a_cur[i + 1] <- numerator / denominator
    }
    
    a_prev <- a_cur
  }
  
  return(a_prev)
}
a <- a_coeffs_kmax(ntip)

rhs_value <- function(x) {
  poly <- 0
  for (i in length(a):1) {
    poly <- poly * x + a[i]
  }
  return(poly)
}

load('/home/users/bingjing/test100_samplesize2_tau1.5.rda')
data=coal_time4[1,]

df <- list(coal_times=data,lineages= ntip:2,intercoal_times= c(
  data[1],
  data[-1] - data[-length(data)]  ), samp_times=0,n_sampled=ntip)



inten2 <- function(t) {
  3 * exp(-t)
}


grid<- seq(0, bound, length.out = Ngrid + 1)
t_vec <-df$coal_times
idx <- findInterval(t_vec, grid, rightmost.closed = TRUE)

t_vec_new<-c(0,t_vec)
I <- numeric(length(t_vec_new) - 1)




integral_between <- function(a, b, grid) { 
  points <- sort(unique(c(a, grid[grid > a & grid < b], b))) 
  
  left <- points[-length(points)] 
  right <- points[-1] 
  mid <- (left + right) / 2 
  
  grid_idx <- findInterval(mid, grid, rightmost.closed = TRUE) 
  
  return(list(
    left = left,
    right = right,
    grid_idx = grid_idx
  ))
}


interval_info <- vector("list", length(t_vec_new) - 1)

for (i in 1:(length(t_vec_new) - 1)) {
  interval_info[[i]] <- integral_between(t_vec_new[i], t_vec_new[i + 1], grid)
}




initC=rev(com_vec)
initC=initC[-length(initC)]

log_lik <- function(f) {
  #f is a vector of size ntip-1
  sum_f<- sum(f[idx])
  
  
  for (i in 1:(length(t_vec_new) - 1)) {
    idx_i <- interval_info[[i]]$grid_idx
    len_i <- interval_info[[i]]$right - interval_info[[i]]$left
    
    I[i] <- sum(exp(-f[idx_i]) * len_i)
  }
  
  llnocoal  <-  initC * I 
  
  Lambda <- sum(exp(-f))*bound/Ngrid
  
  x<-exp(-Lambda)
  val=rhs_value( x) 
  logboundprob=((ntip-1)*log(1-x)+log(val))
  #logboundprob=((ntip-1)*log(1-x)+log(val))
  
  ll <- -sum_f-sum(llnocoal)- logboundprob
  bound_prob <- exp(logboundprob)
  #print(c(2,Lambda,sum(llnocoal),bound_prob,ll))
  print(c(Lambda,bound_prob,logboundprob,ll))
  if (!is.finite(logboundprob)) {
    print(f)
  }
  return(ll)
}


library(MASS) # For mvrnorm

elliptical_slice <- function(initial_theta, prior, lnpdf, pdf_params = list(),
                             cur_lnpdf = NULL, angle_range = NULL) {
  
  D <- length(initial_theta)
  
  if (is.null(cur_lnpdf)) {
    cur_lnpdf <- do.call(lnpdf, c(list(initial_theta), pdf_params))
  }
  
  # Set up the ellipse and the slice threshold
  if (is.null(dim(prior))) { # prior = prior sample
    nu <- prior
  } else { # prior = cholesky decomp
    if (nrow(prior) != D || ncol(prior) != D) {
      stop("Prior must be given by a D-element sample or DxD chol(Sigma)")
    }
    nu <- as.numeric(prior %*% rnorm(D))
  }
  
  hh <- log(runif(1)) + cur_lnpdf
  
  # Set up a bracket of angles and pick a first proposal.
  # "phi = (theta'-theta)" is a change in angle.
  if (is.null(angle_range) || angle_range == 0) {
    # Bracket whole ellipse with both edges at first proposed point
    phi <- runif(1) * 2 * pi
    phi_min <- phi - 2 * pi
    phi_max <- phi
  } else {
    # Randomly center bracket on current point
    phi_min <- -angle_range * runif(1)
    phi_max <- phi_min + angle_range
    phi <- runif(1) * (phi_max - phi_min) + phi_min
  }
  
  # Slice sampling loop
  repeat {
    # Compute xx for proposed angle difference and check if it's on the slice
    xx_prop <- initial_theta * cos(phi) + nu * sin(phi)
    cur_lnpdf <- do.call(lnpdf, c(list(xx_prop), pdf_params))
    
    if (cur_lnpdf > hh) {
      # New point is on slice, ** EXIT LOOP **
      break
    }
    
    # Shrink slice to rejected point
    if (phi > 0) {
      phi_max <- phi
    } else if (phi < 0) {
      phi_min <- phi
    } else {
      stop("BUG DETECTED: Shrunk to current position and still not acceptable.")
    }
    
    # Propose new angle difference
    phi <- runif(1) * (phi_max - phi_min) + phi_min
  }
  
  return(list(xx_prop, cur_lnpdf))
}

expo_quad_kernel <- function(xn, xm) {
  return(min(xn, xm))
}

chol_sample <- function(mean, cov_chol) {
  mean + as.numeric(cov_chol %*% rnorm(length(mean)))
}

loc=grid[-length(grid)]

L <- length(loc)

g_mk <- 3 + rep(0, L)
# 0--(Ngrid-1): function values at grids
# Ngrid--(Ngrid+K-1): K observations
# last: integral term


theta <-1
rem <- 1
alpha <- 0.1
beta <- 0.1

nsim1 <- 10000
nsim2 <- 50000

g_mk_list <- list()
theta_list <- c()









Q_matrix <- function(input, s_noise = 0, signal = 1)
{
  n2 <- length(input)
  diff1 <- diff(input)
  diff1[diff1==0] <- s_noise #correction for dividing over 0
  diff <- (1/(signal*diff1))
  
  Q<-spam::spam(0,n2,n2)  
  if (n2>2)
  {
    Q[cbind(seq(1,n2),seq(1,n2))] <- c(diff[1], diff[1:(n2-2)] + diff[2:(n2-1)],
                                       diff[n2-1]) + (1/signal)*rep(s_noise, n2)
  }
  else
  {
    Q[cbind(seq(1,n2),seq(1,n2))] <- c(diff[1],diff[n2-1])+(1/signal)*rep(s_noise,n2)
  }
  Q[cbind(seq(1,n2-1),seq(2,n2))] <- -diff[1:(n2-1)]
  Q[cbind(seq(2,n2),seq(1,n2-1))] <- -diff[1:(n2-1)]
  
  return(Q)
}

cov_K_mod_inv =Q_matrix(loc)



cov_K_mod_inv_add <- cov_K_mod_inv + diag(L) * noise_var

cov_K_mod <- solve(cov_K_mod_inv_add)
cov_K_mod_chol <- t(chol(cov_K_mod))


I0 <- diag(nrow(cov_K_mod))
max(abs(cov_K_mod%*%cov_K_mod_inv - I0))
max(abs(cov_K_mod%*%cov_K_mod_inv_add - I0))


for (ite in 1:(nsim1 + nsim2)) {
  if (ite %% 100 == 0) {
    print(ite)
  }
  cov_K_chol_final <- cov_K_mod_chol / sqrt(theta)
  
  prior <- chol_sample(
    mean = rep(0, L),
    cov_chol = cov_K_chol_final
  )
  
  out <- elliptical_slice(
    initial_theta = as.numeric(g_mk),
    prior = prior,
    lnpdf = log_lik,
    pdf_params = list(),
    cur_lnpdf = NULL,
    angle_range = NULL
  )
  
  g_mk <- as.numeric(out[[1]])
  
  curloglike <- out[[2]]
  
  g_mk_list[[ite]] <- g_mk
  
  alpha_pos <- alpha + L / 2
  
  beta_pos <- beta + 0.5 * as.numeric(
    t(g_mk) %*% cov_K_mod_inv_add %*% g_mk
  )
  
  theta <- rgamma(
    n = 1,
    shape = alpha_pos,
    rate = beta_pos
  )
  
  theta_list[ite] <- theta
}



exp_g_mk_list <- lapply(g_mk_list, exp)

g_post <- do.call(rbind, exp_g_mk_list[(nsim1 + 1):length(g_mk_list)])






low  <- apply(g_post, 2, quantile, probs = 0.025)
high <- apply(g_post, 2, quantile, probs = 0.975)
med  <- apply(g_post, 2, quantile, probs = 0.5)

length(loc)
groundtruth=sapply(c(loc,bound),inten2)

par(mfrow = c(1, 3))
plot(c(loc,bound), groundtruth, type = "l",
     ylim = c(0, 10),
     col = "black", lwd = 3,
     xlab = "Time",
     ylab = expression(N[e](t)),
     main = paste0("Ne2,ntips=100,bound=1.5,jitter=", noise_var,"Ngrid=",Ngrid))

# magenta color with alpha (for lines)
lines(c(loc,bound),c(med,med[length(med)]),type = "s",pch="",col="magenta",lwd=2.5)
lines(c(loc,bound),c(high,high[length(high)]),type = "s",pch="",col="magenta",lwd=2,lty=3)
lines(c(loc,bound),c(low,low[length(low)]),type = "s",pch="",col="magenta",lwd=2,lty=3)
lines(c(loc,bound),c(med2,med2[length(med2)]),type = "s",pch="",col="cyan",lwd=2.5)
lines(c(loc,bound),c(high2,high2[length(high2)]),type = "s",pch="",col="cyan",lwd=2,lty=3)
lines(c(loc,bound),c(low2,low2[length(low2)]),type = "s",pch="",col="cyan",lwd=2,lty=3)

legend("topright",
       legend = c("Truth", "BC-DISCRETIZED", "SC-DISCRETIZED"),
       col    = c("black", "magenta", "cyan"),
       lwd    = c(3, 2.5, 2.5),
       lty    = c(1, 1, 1))
points(coal_time4[1,],
       rep(0, length(coal_time4[1,])),
       pch = 124,      # vertical line marker '|'
       col = "black",
       cex = 0.6)      # controls size (≈ markersize)
g_mat <- do.call(rbind, g_mk_list[(nsim1 + 1):length(g_mk_list)])
plot(g_mat[, 500], type = "l",main = "midpoint traceplot under BC")
g_mat <- do.call(rbind, g_mk_list2[(nsim1 + 1):length(g_mk_list2)])
plot(g_mat[, 500], type = "l",main = "midpoint traceplot under SC")

