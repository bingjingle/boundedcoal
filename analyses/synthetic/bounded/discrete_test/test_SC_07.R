ntip <- 100
Ngrid=1000

noise_var <- 0.001
bound=0.7
com_vec <- choose(seq_len(ntip), 2)



data_ori=c(
  0.0006457, 0.00094795, 0.00097818, 0.00115723, 0.00136032,
  0.0014971, 0.00188022, 0.00229853, 0.00275317, 0.00414087,
  0.00510088, 0.00546055, 0.00724715, 0.0075598, 0.00843643,
  0.01121559, 0.01140645, 0.01150385, 0.01385105, 0.01563069,
  0.01782013, 0.01795062, 0.01827905, 0.01898047, 0.01989265,
  0.01997209, 0.02033843, 0.02223489, 0.02405914, 0.02455844,
  0.02490612, 0.02514106, 0.02523667, 0.0254744, 0.02547566,
  0.02562561, 0.02600492, 0.02673282, 0.02833923, 0.02908865,
  0.03005297, 0.03110919, 0.03234614, 0.03744378, 0.03752536,
  0.03869813, 0.04311578, 0.04386097, 0.04909827, 0.05170783,
  0.05189045, 0.05292001, 0.05726541, 0.05837406, 0.06083225,
  0.06454715, 0.06758013, 0.07431144, 0.07451043, 0.07548091,
  0.07577013, 0.07728782, 0.08168786, 0.08383543, 0.08460874,
  0.0851407, 0.08686686, 0.08687484, 0.08774483, 0.08805298,
  0.09739123, 0.10208261, 0.10602849, 0.11395216, 0.13507485,
  0.14032089, 0.18690144, 0.19247492, 0.19934897, 0.22431027,
  0.22508206, 0.2307226, 0.23352532, 0.23890136, 0.24700717,
  0.25998017, 0.27098572, 0.28242977, 0.28818101, 0.2978475,
  0.33223878, 0.40785599, 0.42174163, 0.48693148, 0.49802362,
  0.52905337, 0.58797086, 0.62521005, 0.66886555
)
data <- list(coal_times=data_ori,lineages= ntip:2,intercoal_times= c(
  data_ori[1],
  data_ori[-1] - data_ori[-length(data_ori)]  ), samp_times=0,n_sampled=ntip)


inten2 <- function(t) {
  3 * exp(-t)
}


grid<- seq(0, bound, length.out = Ngrid + 1)
t_vec <-data_ori
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
  
  
  
  ll <- -sum_f-sum(llnocoal)
  
  #print(c(Lambda,bound_prob,logboundprob,ll))
  
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

g_mk_list2 <- list()
theta_list2 <- c()









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
  
  g_mk_list2[[ite]] <- g_mk
  
  alpha_pos <- alpha + L / 2
  
  beta_pos <- beta + 0.5 * as.numeric(
    t(g_mk) %*% cov_K_mod_inv_add %*% g_mk
  )
  
  theta <- rgamma(
    n = 1,
    shape = alpha_pos,
    rate = beta_pos
  )
  
  theta_list2[ite] <- theta
}



exp_g_mk_list2 <- lapply(g_mk_list2, exp)
g_post2 <- do.call(rbind, exp_g_mk_list2[(nsim1 + 1):length(g_mk_list2)])
low2  <- apply(g_post2, 2, quantile, probs = 0.025)
high2 <- apply(g_post2, 2, quantile, probs = 0.975)
med2 <- apply(g_post2, 2, quantile, probs = 0.5)
